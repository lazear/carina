use carina::database::{IndexedDatabase, PeptideIx, Theoretical};
use carina::ion_series::Kind;
use carina::mass::{Tolerance, PROTON};
use carina::peptide::TargetDecoy;
use carina::spectrum::{read_ms2, ProcessedSpectrum, SpectrumProcessor};
use clap::{Arg, Command};
use indicatif::ParallelProgressIterator;
use log::info;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use std::time::{self, Instant};

#[derive(Copy, Clone)]
pub struct Score {
    peptide: PeptideIx,
    matched_b: u32,
    matched_y: u32,
    summed_b: f32,
    summed_y: f32,
    q_value: f32,
    hyperscore: f32,
}

#[derive(Serialize)]
pub struct Percolator<'db> {
    peptide: String,
    proteins: &'db str,
    specid: usize,
    scannr: u32,
    label: i32,
    expmass: f32,
    calcmass: f32,
    charge: u8,
    rt: f32,
    delta_mass: f32,
    hyperscore: f32,
    deltascore: f32,
    matched_peaks: u32,
    percent_matched_peaks: f32,
    matched_intensity: f32,
    percent_matched_intensity: f32,
    total_candidates: usize,
    spectrum_q_value: f32,
    q_value: f32,
}

impl Score {
    /// Calculate the X!Tandem hyperscore
    /// * `fact_table` is a precomputed vector of factorials
    fn hyperscore(&self, fact_table: &[f32]) -> f32 {
        let i = (self.summed_b + 1.0) * (self.summed_y + 1.0);
        let m = fact_table[(self.matched_b as usize).min(fact_table.len() - 2)]
            * fact_table[(self.matched_y as usize).min(fact_table.len() - 2)];

        let score = i.ln() + m.ln();
        if score.is_finite() {
            score
        } else {
            f32::MAX
        }
    }

    pub fn new(peptide: &Theoretical) -> Self {
        Score {
            peptide: peptide.peptide_index,
            matched_b: 0,
            matched_y: 0,
            summed_b: 0.0,
            summed_y: 0.0,
            q_value: 1.0,
            hyperscore: 0.0,
        }
    }
}

pub struct Scorer<'db> {
    db: &'db IndexedDatabase,
    search: &'db Search,
    factorial: [f32; 32],
}

impl<'db> Scorer<'db> {
    pub fn new(db: &'db IndexedDatabase, search: &'db Search) -> Self {
        let mut factorial = [1.0f32; 32];
        for i in 1..32 {
            factorial[i] = factorial[i - 1] * i as f32;
        }

        debug_assert!(factorial[3] == 6.0);
        debug_assert!(factorial[4] == 24.0);

        Scorer {
            db,
            search,
            factorial,
        }
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score<'s>(&self, query: &ProcessedSpectrum) -> Vec<Percolator<'db>> {
        let mut scores: HashMap<PeptideIx, Score> = HashMap::new();

        // Create a new `IndexedQuery`
        let candidates = self
            .db
            .query(query, self.search.precursor_tol, self.search.fragment_tol);

        let mut total_intensity = 0.0;
        for (fragment_mz, intensity) in query.peaks.iter() {
            total_intensity += intensity;
            for frag in candidates.page_search(*fragment_mz) {
                let mut sc = scores
                    .entry(frag.peptide_index)
                    .or_insert_with(|| Score::new(frag));

                match frag.kind {
                    Kind::B => {
                        sc.matched_b += 1;
                        sc.summed_b += intensity
                    }
                    Kind::Y => {
                        sc.matched_y += 1;
                        sc.summed_y += intensity
                    }
                }
            }
        }

        if scores.is_empty() {
            return Vec::new();
        }

        // Now that we have processed all candidates, calculate the hyperscore
        let mut scores = scores
            .into_values()
            .map(|mut sc| {
                sc.hyperscore = sc.hyperscore(&self.factorial);
                sc
            })
            .collect::<Vec<_>>();

        scores.sort_by(|b, a| a.hyperscore.total_cmp(&b.hyperscore));

        let mut decoy = 1;
        let mut target = 0;
        for score in scores.iter_mut() {
            match self.db[score.peptide] {
                TargetDecoy::Target(_) => target += 1,
                TargetDecoy::Decoy(_) => decoy += 1,
            }
            score.q_value = decoy as f32 / target as f32;
        }
        let mut q_min = 1.0f32;
        for score in scores.iter_mut().rev() {
            q_min = q_min.min(score.q_value);
            score.q_value = q_min;
        }

        let mut reporting = Vec::new();

        for idx in 0..self.search.report_psms.min(scores.len()) {
            let better = scores[idx];
            let next = scores
                .get(idx + 1)
                .map(|score| score.hyperscore)
                .unwrap_or_default();

            let peptide = self.db[better.peptide].peptide();
            reporting.push(Percolator {
                peptide: peptide.to_string(),
                proteins: &peptide.protein,
                specid: 0,
                scannr: query.scan,
                label: self.db[better.peptide].label(),
                expmass: query.monoisotopic_mass + PROTON,
                calcmass: peptide.monoisotopic + PROTON,
                charge: query.charge,
                rt: query.rt,
                delta_mass: (query.monoisotopic_mass - peptide.monoisotopic),
                hyperscore: better.hyperscore,
                deltascore: better.hyperscore - next,
                matched_peaks: better.matched_b + better.matched_y,
                percent_matched_peaks: (better.matched_b + better.matched_y) as f32
                    / query.peaks.len() as f32,
                matched_intensity: better.summed_b + better.summed_y,
                percent_matched_intensity: (better.summed_b + better.summed_y) / total_intensity,
                total_candidates: scores.len(),
                spectrum_q_value: better.q_value,
                q_value: 1.0,
            })
        }
        reporting
    }

    /// Assign q_values in place to a set of PSMs
    ///
    /// # Invariants
    /// * `scores` must be sorted in descending order (e.g. best PSM is first)
    pub fn assign_q_values(&self, scores: &mut [Percolator]) -> usize {
        // FDR Calculation:
        // * Sort by score, descending
        // * Estimate FDR
        // * Calculate q-value

        let mut decoy = 1;
        let mut target = 0;

        for score in scores.iter_mut() {
            match score.label == -1 {
                true => decoy += 1,
                false => target += 1,
            }
            score.q_value = decoy as f32 / target as f32;
        }

        // Reverse slice, and calculate the cumulative minimum
        let mut q_min = 1.0f32;
        let mut passing = 0;
        for score in scores.iter_mut().rev() {
            q_min = q_min.min(score.q_value);
            score.q_value = q_min;
            if q_min <= 0.01 {
                passing += 1;
            }
        }
        passing
    }
}

#[derive(Serialize)]
pub struct Search {
    database: carina::database::Parameters,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    max_fragment_charge: u8,
    min_peaks: usize,
    max_peaks: usize,
    report_psms: usize,
    ms2_paths: Vec<String>,
    pin_paths: Vec<String>,
    search_time: f32,
}

#[derive(Deserialize)]
struct Input {
    database: carina::database::Builder,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    report_psms: Option<usize>,
    min_peaks: Option<usize>,
    max_peaks: Option<usize>,
    max_fragment_charge: Option<u8>,
    ms2_paths: Vec<String>,
}

impl Search {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let mut file = std::fs::File::open(path)?;
        let request: Input = serde_json::from_reader(&mut file)?;
        let database = request.database.make_parameters();
        Ok(Search {
            database,
            precursor_tol: request.precursor_tol,
            fragment_tol: request.fragment_tol,
            report_psms: request.report_psms.unwrap_or(1),
            max_peaks: request.max_peaks.unwrap_or(150),
            min_peaks: request.min_peaks.unwrap_or(15),
            max_fragment_charge: request.max_fragment_charge.unwrap_or(3),
            pin_paths: Vec::new(),
            ms2_paths: request.ms2_paths,
            search_time: 0.0,
        })
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let env = env_logger::Env::default().filter_or("CARINA_LOG", "info");
    env_logger::init_from_env(env);

    let start = time::Instant::now();

    let matches = Command::new("carina")
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .arg(Arg::new("parameters").required(true))
        .get_matches();

    let path = matches
        .get_one::<String>("parameters")
        .expect("required parameters");
    let mut search = Search::load(path)?;

    let db = search.database.clone().build()?;

    let buckets = db.buckets();
    let mut avg_delta = 0.0;

    for i in 1..buckets.len() {
        let delta = buckets[i] - buckets[i - 1];
        avg_delta += delta;
    }
    dbg!(avg_delta / buckets.len() as f32);

    info!(
        "generated {} fragments in {}ms",
        db.size(),
        (Instant::now() - start).as_millis()
    );

    let scorer = Scorer::new(&db, &search);
    let sp = SpectrumProcessor::new(
        search.max_peaks,
        search.max_fragment_charge,
        search.database.fragment_max_mz,
    );

    let mut pin_paths = Vec::with_capacity(search.ms2_paths.len());
    for ms2_path in &search.ms2_paths {
        let start = Instant::now();
        let mut scores = read_ms2(ms2_path)?
            .into_par_iter()
            .filter(|spec| spec.peaks.len() >= search.min_peaks)
            .flat_map(|spectra| scorer.score(&sp.process(spectra)))
            .collect::<Vec<_>>();
        let duration = Instant::now() - start;

        let pin_path = format!("{}.carina.pin", ms2_path);
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&pin_path)?;

        scores.sort_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));
        let passing_psms = scorer.assign_q_values(&mut scores);
        // let passing_psms = 0;
        let total_psms = scores.len();

        for (idx, mut score) in scores.into_iter().enumerate() {
            score.specid = idx;
            writer.serialize(score)?;
        }

        pin_paths.push(pin_path);

        info!(
            "{:?}: assigned {} PSMs ({} with 1% FDR) in {} ms ({} PSMs/sec)",
            ms2_path,
            total_psms,
            passing_psms,
            duration.as_millis(),
            total_psms as f32 / duration.as_secs_f32()
        );
    }

    search.search_time = (Instant::now() - start).as_secs_f32();
    search.pin_paths = pin_paths;
    let results = serde_json::to_string_pretty(&search)?;

    eprintln!("{}", &results);
    std::fs::write("results.json", results)?;

    Ok(())
}
