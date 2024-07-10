use crate::database::{IndexedDatabase, PeptideIx};
use crate::heap::bounded_min_heapify;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::spectrum::{Precursor, ProcessedSpectrum};
use serde::Serialize;
use std::ops::AddAssign;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default, Debug)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    longest_b: usize,
    longest_y: usize,
    hyperscore: f64,
    ppm_difference: f32,
    precursor_charge: u8,
    isotope_error: i8,
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PreScore {
    matched: u16,
    peptide: PeptideIx,
    precursor_charge: u8,
    isotope_error: i8,
}

/// Store preliminary scores & stats for first pass search for a query spectrum
#[derive(Clone, Default)]
struct InitialHits {
    matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    scored_candidates: usize,
    preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;

        self.preliminary.extend(rhs.preliminary);
    }
}

#[derive(Serialize, Clone, Debug)]
/// Features of a candidate peptide spectrum match
pub struct Feature {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    // psm_id help to match with matched fragments table.
    pub psm_id: usize,
    pub peptide_len: usize,
    /// Spectrum id
    pub spec_id: String,
    /// Precursor id
    pub prec_id: usize,
    /// File identifier
    pub file_id: usize,
    /// PSM rank
    pub rank: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass
    pub expmass: f32,
    /// Calculated mass
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Globally aligned retention time
    pub aligned_rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt_model: f32,
    /// Ion mobility
    pub ims: f32,
    /// Predicted ion mobility, if enabled
    pub predicted_ims: f32,
    /// Difference between predicted & observed ion mobility
    pub delta_ims_model: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// C13 isotope error
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_next: f64,
    /// Difference between hyperscore of this candidate, and the best candidate
    pub delta_best: f64,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Longest b-ion series
    pub longest_b: u32,
    /// Longest y-ion series
    pub longest_y: u32,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,

    pub ms2_intensity: f32,

    pub fragments: Option<Fragments>,
}

/// Matching Fragment details
#[derive(Serialize, Default, Clone, Debug)]
pub struct Fragments {
    #[serde(skip_serializing)]
    pub charges: Vec<i32>,
    pub kinds: Vec<Kind>,
    pub fragment_ordinals: Vec<i32>,
    pub intensities: Vec<f32>,
    pub mz_calculated: Vec<f32>,
    pub mz_experimental: Vec<f32>,
}

static PSM_COUNTER: AtomicUsize = AtomicUsize::new(1);

fn increment_psm_counter() -> usize {
    PSM_COUNTER.fetch_add(1, Ordering::Relaxed)
}

/// Stirling's approximation for log factorial
fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

impl Score {
    /// Calculate the X!Tandem hyperscore
    /// * `fact_table` is a precomputed vector of factorials
    fn hyperscore(&self) -> f64 {
        let i = (self.summed_b + 1.0) as f64 * (self.summed_y + 1.0) as f64;
        let score = i.ln() + lnfact(self.matched_b) + lnfact(self.matched_y);
        if score.is_finite() {
            score
        } else {
            255.0
        }
    }
}

pub struct Scorer<'db> {
    pub db: &'db IndexedDatabase,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    /// What is the minimum number of matched b and y ion peaks to report PSMs for?
    pub min_matched_peaks: u16,
    /// Precursor isotope error lower bounds (e.g. -1)
    pub min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub max_fragment_charge: Option<u8>,
    pub min_fragment_mass: f32,
    pub max_fragment_mass: f32,
    pub chimera: bool,
    pub report_psms: usize,

    // Rather than use a fixed precursor tolerance, dynamically alter
    // the precursor tolerance window based on MS2 isolation window and charge
    pub wide_window: bool,
    pub annotate_matches: bool,
}

#[inline(always)]
/// Calculate upper bound (excluded) of the charge state range to use for
/// searching fragment ions (1..N)
/// If user has configured max_fragment_charge, potentially override precursor
/// charge
fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

impl<'db> Scorer<'db> {
    pub fn score(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera_fast(query),
            false => self.score_standard(query),
        }
    }

    /// Perform a k-select and truncation of an [`InitialHits`] list.
    ///
    /// Determine how many candidates to actually calculate hyperscore for.
    /// Hyperscore is relatively computationally expensive, so we don't want
    /// to calculate it for every possible candidate (100s - 10,000s depending on search)
    /// when we are only going to report a few PSMs. But we also want to calculate
    /// it for enough candidates that we don't accidentally miss the best hit!
    ///
    /// Given that hyperscore is dominated by the number of matched peaks, it seems
    /// reasonable to assume that the highest hyperscore will belong to one of the
    /// top 50 candidates sorted by # of matched peaks.
    fn trim_hits(&self, hits: &mut InitialHits) {
        let k = 50.clamp(
            (self.report_psms * 2).min(hits.preliminary.len()),
            hits.preliminary.len(),
        );
        bounded_min_heapify(&mut hits.preliminary, k);
        hits.preliminary.truncate(k);
        // Q: could default peptide idx be done here?
    }

    /// Preliminary Score, return # of matched peaks per candidate
    /// Returned hits are guaranteed to be the top-K hits (see above comment)
    /// from among all potential candidates, but the returned vector is not
    /// in sorted order.
    fn matched_peaks_with_isotope(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
    ) -> InitialHits {
        let candidates = self.db.query(
            precursor_mass - isotope_error as f32 * NEUTRON,
            precursor_tol,
            self.fragment_tol,
        );

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        // Allocate space for all potential candidates - many potential candidates
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut hits = InitialHits {
            matched_peaks: 0,
            scored_candidates: 0,
            preliminary: vec![PreScore::default(); potential],
        };

        for peak in query.peaks.iter() {
            for charge in 1..max_fragment_charge {
                let mass = peak.mass * charge as f32;
                for frag in candidates.page_search(mass) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    let sc = &mut hits.preliminary[idx];
                    if sc.matched == 0 {
                        hits.scored_candidates += 1;
                        sc.precursor_charge = precursor_charge;
                        sc.peptide = frag.peptide_index;
                        sc.isotope_error = isotope_error;
                    }
                    sc.matched += 1;
                    hits.matched_peaks += 1;
                }
            }
        }
        if hits.matched_peaks == 0 {
            return hits;
        }

        self.trim_hits(&mut hits);
        hits
    }

    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        query,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                        isotope,
                    );
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        } else {
            self.matched_peaks_with_isotope(
                query,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                0,
            )
        }
    }

    fn initial_hits(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> InitialHits {
        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        // Search in wide-window/DIA mode
        if self.wide_window {
            let mut hits = (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    let precursor_tol = precursor
                        .isolation_window
                        .unwrap_or(Tolerance::Da(-2.4, 2.4))
                        * precursor_charge as f32;
                    hits +=
                        self.matched_peaks(query, precursor_mass, precursor_charge, precursor_tol);
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        } else if let Some(charge) = precursor.charge {
            // Charge state is already annotated for this precusor, only search once
            let precursor_mass = mz * charge as f32;
            self.matched_peaks(query, precursor_mass, charge, self.precursor_tol)
        } else {
            // Not all selected ion precursors have charge states annotated -
            // assume it could be z=2, z=3, z=4 and search all three
            let mut hits = (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks(
                        query,
                        precursor_mass,
                        precursor_charge,
                        self.precursor_tol,
                    );
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        }
    }

    fn notched_initial_hits(&self, query: &ProcessedSpectrum) -> (InitialHits, Vec<usize>) {
        let mut cumsum = 0;
        let mut hits = InitialHits::default();
        // Match lengths is pre cumulative sum of the number of hits for each precursor
        let mut cum_match_lengths = Vec::with_capacity(query.precursors.len());

        for precursor in &query.precursors {
            let precursor_hits = self.initial_hits(query, precursor);
            cumsum += precursor_hits.preliminary.len();
            hits += precursor_hits;
            cum_match_lengths.push(cumsum);
        }

        (hits, cum_match_lengths)
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        if query.precursors.is_empty() {
            panic!("missing precursor for {}", query.id);
        }
        let (hits, match_lens) = self.notched_initial_hits(query);

        let mut features = Vec::with_capacity(self.report_psms);
        self.build_features(query, &query.precursors, &hits, self.report_psms, &mut features, &match_lens);
        features
    }

    /// Given a set of [`InitialHits`] against a query spectrum, prepare N=`report_psms`
    /// best PSMs ([`Feature`])
    fn build_features(
        &self,
        query: &ProcessedSpectrum,
        precursor: &[Precursor],
        hits: &InitialHits,
        report_psms: usize,
        features: &mut Vec<Feature>,
        index_rle: &[usize],
    ) {
        let mut score_vector = Vec::with_capacity(hits.preliminary.len());
        for (idx, pre) in (0u8..).zip(hits.preliminary.iter()) {
            if pre.peptide != PeptideIx::default() {
                let score = self.score_candidate(query, pre);
                if (score.0.matched_b + score.0.matched_y) >= self.min_matched_peaks {
                    score_vector.push((score, idx));
                }
            }
        };

        // Hyperscore is our primary score function for PSMs
        score_vector.sort_by(|(a, _i), (b, _ii)| b.0.hyperscore.total_cmp(&a.0.hyperscore));

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        for idx in 0..report_psms.min(score_vector.len()) {
            let score = score_vector[idx].0.0;
            let score_index = score_vector[idx].1 as usize;
            let mut precursor_index: usize = 0;
            while score_index >= index_rle[precursor_index] {
                precursor_index += 1;
            }

            let fragments: Option<Fragments> = score_vector[idx].0.1.take();
            let psm_id = increment_psm_counter();

            let peptide = &self.db[score.peptide];

            let next = score_vector
                .get(idx + 1)
                .map(|score| score.0.0.hyperscore)
                .unwrap_or_default();

            let best = score_vector
                .get(0)
                .map(|score| score.0.0.hyperscore)
                .expect("we know that index 0 is valid");

            // Poisson distribution probability mass function
            let k = score.matched_b + score.matched_y;
            let mut poisson = lambda.powi(k as i32) * f64::exp(-lambda) / lnfact(k).exp();

            if poisson.is_infinite() {
                // Approximately the smallest positive non-zero value representable by f64
                poisson = 1E-325;
            }

            // Sage operates on masses without protons; [M] instead of [MH+]
            let mz = precursor[precursor_index].mz - PROTON;
            let precursor_mass = mz * score.precursor_charge as f32;
            let isotope_error = score.isotope_error as f32 * NEUTRON;
            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error).abs() * 2E6
                / (precursor_mass - isotope_error + peptide.monoisotopic);

            // let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            features.push(Feature {
                // Identifiers
                psm_id,
                peptide_idx: score.peptide,
                spec_id: query.id.clone(),
                prec_id: precursor_index,
                file_id: query.file_id,
                rank: idx as u32 + 1,
                label: peptide.label(),
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: score.precursor_charge,
                rt: query.scan_start_time,
                ims: query
                    .precursors
                    .get(0)
                    .unwrap()
                    .inverse_ion_mobility
                    .unwrap_or(0.0),
                delta_mass,
                isotope_error,
                average_ppm: score.ppm_difference,
                hyperscore: score.hyperscore,
                delta_next: score.hyperscore - next,
                delta_best: best - score.hyperscore,
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (score.summed_b + score.summed_y)
                    / query.total_ion_current,
                poisson: poisson.log10(),
                longest_b: score.longest_b as u32,
                longest_y: score.longest_y as u32,
                longest_y_pct: score.longest_y as f32 / (peptide.sequence.len() as f32),
                peptide_len: peptide.sequence.len(),
                scored_candidates: hits.scored_candidates as u32,
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                spectrum_q: 1.0,
                protein_q: 1.0,
                peptide_q: 1.0,
                predicted_rt: 0.0,
                predicted_ims: 0.0,
                aligned_rt: query.scan_start_time,
                delta_rt_model: 0.999,
                delta_ims_model: 0.999,
                ms2_intensity: score.summed_b + score.summed_y,

                //Fragments
                fragments,
            })
        }
    }

    /// Remove peaks matching a PSM from a query spectrum
    fn remove_matched_peaks(&self, query: &mut ProcessedSpectrum, psm: &Feature) {
        let peptide = &self.db[psm.peptide_idx];
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind));

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, psm.charge);

        // Remove MS2 peaks matched by previous match
        let mut to_remove = Vec::new();
        for frag in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                if let Some(peak) = crate::spectrum::select_most_intense_peak(
                    &query.peaks,
                    frag.monoisotopic_mass / charge as f32,
                    self.fragment_tol,
                    None,
                ) {
                    to_remove.push(*peak);
                }
            }
        }

        query.peaks = query
            .peaks
            .drain(..)
            .filter(|peak| !to_remove.contains(peak))
            .collect();
        query.total_ion_current = query.peaks.iter().map(|peak| peak.intensity).sum::<f32>();
    }

    /// Return multiple PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed, etc
    pub fn score_chimera_fast(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        if query.precursors.is_empty() {
            panic!("missing precursor for {}", query.id);
        }

        let mut query = query.clone();
        let (hits, match_lengths) = self.notched_initial_hits(&query);

        let mut candidates: Vec<Feature> = Vec::with_capacity(self.report_psms);

        let mut prev = 0;
        while candidates.len() < self.report_psms {
            self.build_features(&query, &query.precursors, &hits, 1, &mut candidates, &match_lengths);
            if candidates.len() > prev {
                if let Some(feat) = candidates.get_mut(prev) {
                    self.remove_matched_peaks(&mut query, feat);
                    feat.rank = prev as u32 + 1;
                }
                prev = candidates.len()
            } else {
                break;
            }
        }
        candidates
    }

    /// Calculate full hyperscore for a given PSM
    fn score_candidate(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
    ) -> (Score, Option<Fragments>) {
        let mut score = Score {
            peptide: pre_score.peptide,
            precursor_charge: pre_score.precursor_charge,
            isotope_error: pre_score.isotope_error,
            ..Default::default()
        };
        let peptide = &self.db[score.peptide];
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, score.precursor_charge);

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind).enumerate());

        let mut b_run = Run::default();
        let mut y_run = Run::default();

        let mut fragments_details = Fragments::default();

        for (idx, frag) in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                let mz = frag.monoisotopic_mass / charge as f32;
                if let Some(peak) = crate::spectrum::select_most_intense_peak(
                    &query.peaks,
                    mz,
                    self.fragment_tol,
                    None,
                ) {
                    score.ppm_difference +=
                        peak.intensity * (mz - peak.mass).abs() * 2E6 / (mz + peak.mass);

                    let exp_mz = peak.mass + PROTON;
                    let calc_mz = mz + PROTON;

                    match frag.kind {
                        Kind::A | Kind::B | Kind::C => {
                            score.matched_b += 1;
                            score.summed_b += peak.intensity;
                            b_run.matched(idx);
                        }
                        Kind::X | Kind::Y | Kind::Z => {
                            score.matched_y += 1;
                            score.summed_y += peak.intensity;
                            y_run.matched(idx);
                        }
                    }

                    if self.annotate_matches {
                        let idx = match frag.kind {
                            Kind::A | Kind::B | Kind::C => idx as i32 + 1,
                            Kind::X | Kind::Y | Kind::Z => {
                                peptide.sequence.len().saturating_sub(1) as i32 - idx as i32
                            }
                        };
                        fragments_details.kinds.push(frag.kind);
                        fragments_details.charges.push(charge as i32);
                        fragments_details.mz_experimental.push(exp_mz);
                        fragments_details.mz_calculated.push(calc_mz);
                        fragments_details.fragment_ordinals.push(idx);
                        fragments_details.intensities.push(peak.intensity);
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore();
        score.longest_b = b_run.longest;
        score.longest_y = y_run.longest;
        score.ppm_difference /= score.summed_b + score.summed_y;

        if self.annotate_matches {
            (score, Some(fragments_details))
        } else {
            // drop(fragments_details);
            (score, None)
        }
    }
}

/// Maintain information about the longest continous ion ladder for a series
#[derive(Default)]
struct Run {
    start: usize,
    length: usize,
    last: usize,
    pub longest: usize,
}

impl Run {
    pub fn matched(&mut self, index: usize) {
        if self.last == index {
            return;
        } else if self.start + self.length == index {
            self.length += 1;
            self.longest = self.longest.max(self.length);
        } else {
            self.start = index;
            self.length = 1;
            self.longest = self.longest.max(self.length);
        }
        self.last = index;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn longest_series() {
        let mut run = Run::default();

        run.matched(1);
        run.matched(2);
        run.matched(3);
        run.matched(3);
        run.matched(3);

        assert_eq!(run.length, 3);
        assert_eq!(run.longest, 3);

        run.matched(5);
        run.matched(5);
        assert_eq!(run.length, 1);
        assert_eq!(run.longest, 3);
        run.matched(6);
        assert_eq!(run.length, 2);
    }

    #[test]
    fn test_max_fragment_charge() {
        assert_eq!(max_fragment_charge(None, 1), 2);
        assert_eq!(max_fragment_charge(None, 2), 2);
        assert_eq!(max_fragment_charge(None, 3), 3);
        assert_eq!(max_fragment_charge(None, 4), 4);
        assert_eq!(max_fragment_charge(Some(1), 2), 2);
        assert_eq!(max_fragment_charge(Some(1), 3), 2);
        assert_eq!(max_fragment_charge(Some(2), 4), 3);
        assert_eq!(max_fragment_charge(Some(4), 1), 2);
    }
}
