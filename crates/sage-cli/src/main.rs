use anyhow::Context;
use clap::{value_parser, Arg, Command, ValueHint};
use input::{Input, Search};
use log::info;
use rayon::prelude::*;
use sage_cloudpath::CloudPath;
use sage_core::database::{IndexedDatabase, Parameters, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::mass::Tolerance;
use sage_core::peptide::Peptide;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{ProcessedSpectrum, SpectrumProcessor};
use sage_core::tmt::TmtQuant;
use std::collections::HashSet;
use std::time::Instant;

mod input;
mod output;
mod telemetry;

struct Runner {
    database: IndexedDatabase,
    parameters: input::Search,
    start: Instant,
}

#[derive(Default)]
struct SageResults {
    ms1: Vec<ProcessedSpectrum>,
    features: Vec<Feature>,
    quant: Vec<TmtQuant>,
}

impl FromParallelIterator<SageResults> for SageResults {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = SageResults>,
    {
        par_iter
            .into_par_iter()
            .reduce(SageResults::default, |mut acc, x| {
                acc.features.extend(x.features);
                acc.quant.extend(x.quant);
                acc.ms1.extend(x.ms1);
                acc
            })
    }
}

impl FromIterator<SageResults> for SageResults {
    fn from_iter<I>(par_iter: I) -> Self
    where
        I: IntoIterator<Item = SageResults>,
    {
        par_iter
            .into_iter()
            .fold(SageResults::default(), |mut acc, x| {
                acc.features.extend(x.features);
                acc.quant.extend(x.quant);
                acc.ms1.extend(x.ms1);
                acc
            })
    }
}

impl Runner {
    pub fn new(parameters: Search, parallel: usize) -> anyhow::Result<Self> {
        let start = Instant::now();
        let fasta = sage_cloudpath::util::read_fasta(
            &parameters.database.fasta,
            &parameters.database.decoy_tag,
            parameters.database.generate_decoys,
        )
        .with_context(|| {
            format!(
                "Failed to build database from `{}`",
                parameters.database.fasta
            )
        })?;

        let database = match parameters.database.prefilter {
            false => parameters.database.clone().build(fasta),
            true => {
                let mini_runner = Self {
                    database: IndexedDatabase::default(),
                    parameters: parameters.clone(),
                    start,
                };
                let peptides = mini_runner.prefilter_peptides(parallel, fasta);
                parameters.database.clone().build_from_peptides(peptides)
            }
        };
        info!(
            "generated {} fragments, {} peptides in {}ms",
            database.fragments.len(),
            database.peptides.len(),
            (Instant::now() - start).as_millis()
        );
        Ok(Self {
            database,
            parameters,
            start,
        })
    }

    pub fn prefilter_peptides(self, parallel: usize, fasta: Fasta) -> Vec<Peptide> {
        let spectra: Option<Vec<ProcessedSpectrum>> =
            match parallel >= self.parameters.mzml_paths.len() {
                true => Some(self.read_processed_spectra(&self.parameters.mzml_paths, 0, 0)),
                false => None,
            };
        let mut all_peptides: Vec<Peptide> = fasta
            .iter_chunks(self.parameters.database.prefilter_chunk_size)
            .enumerate()
            .flat_map(|(chunk_id, fasta_chunk)| {
                let start = Instant::now();
                info!("pre-filtering fasta chunk {}", chunk_id,);
                let db = &self.parameters.database.clone().build(fasta_chunk);
                info!(
                    "generated {} fragments, {} peptides in {}ms",
                    db.fragments.len(),
                    db.peptides.len(),
                    (Instant::now() - start).as_millis()
                );
                let scorer = Scorer {
                    db,
                    precursor_tol: self.parameters.precursor_tol,
                    fragment_tol: self.parameters.fragment_tol,
                    min_matched_peaks: self.parameters.min_matched_peaks,
                    min_isotope_err: self.parameters.isotope_errors.0,
                    max_isotope_err: self.parameters.isotope_errors.1,
                    min_precursor_charge: self.parameters.precursor_charge.0,
                    max_precursor_charge: self.parameters.precursor_charge.1,
                    override_precursor_charge: self.parameters.override_precursor_charge,
                    max_fragment_charge: self.parameters.max_fragment_charge,
                    chimera: self.parameters.chimera,
                    report_psms: self.parameters.report_psms + 1,
                    wide_window: self.parameters.wide_window,
                    annotate_matches: self.parameters.annotate_matches,
                };
                let peptide_idxs: HashSet<PeptideIx> = match &spectra {
                    Some(spectra) => self.peptide_filter_processed_spectra(&scorer, spectra),
                    None => self
                        .parameters
                        .mzml_paths
                        .chunks(parallel)
                        .enumerate()
                        .flat_map(|(chunk_idx, chunk)| {
                            let spectra_chunk =
                                self.read_processed_spectra(chunk, chunk_idx, parallel);
                            self.peptide_filter_processed_spectra(&scorer, &spectra_chunk)
                        })
                        .collect(),
                }
                .into_iter()
                .collect();
                let peptides: Vec<Peptide> = peptide_idxs
                    .into_iter()
                    .map(|idx| db[idx].clone())
                    .collect();
                info!(
                    "found {} pre-filtered peptides for fasta chunk {}",
                    peptides.len(),
                    chunk_id,
                );
                peptides
            })
            .collect();
        Parameters::reorder_peptides(&mut all_peptides);
        all_peptides
    }

    fn peptide_filter_processed_spectra(
        &self,
        scorer: &Scorer,
        spectra: &Vec<ProcessedSpectrum>,
    ) -> Vec<PeptideIx> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        let peptide_idxs: Vec<_> = spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .map(|x| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev % 10_000 == 0 {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                x
            })
            .flat_map(|spec| {
                scorer.quick_score(spec, self.parameters.database.prefilter_low_memory)
            })
            .collect();

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!("- search:  {:8} ms ({} spectra/s)", duration, rate);
        peptide_idxs
    }

    fn spectrum_fdr(&self, features: &mut [Feature]) -> usize {
        if sage_core::ml::linear_discriminant::score_psms(features, self.parameters.precursor_tol)
            .is_none()
        {
            log::warn!("linear model fitting failed, falling back to heuristic discriminant score");
            features.par_iter_mut().for_each(|feat| {
                feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
            });
        }
        features.par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));
        sage_core::ml::qvalue::spectrum_q_value(features)
    }

    // Create a path for `file_name` in the specified output directory, if it exists,
    // otherwise, write to current directory
    fn make_path<S: AsRef<str>>(&self, file_name: S) -> CloudPath {
        let mut path = self.parameters.output_directory.clone();
        path.push(file_name);
        path
    }

    fn search_processed_spectra(
        &self,
        scorer: &Scorer,
        spectra: &Vec<ProcessedSpectrum>,
    ) -> Vec<Feature> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let counter = AtomicUsize::new(0);
        let start = Instant::now();

        let features: Vec<_> = spectra
            .par_iter()
            .filter(|spec| spec.peaks.len() >= self.parameters.min_peaks && spec.level == 2)
            .map(|x| {
                let prev = counter.fetch_add(1, Ordering::Relaxed);
                if prev > 0 && prev % 10_000 == 0 {
                    let duration = Instant::now().duration_since(start).as_millis() as usize;

                    let rate = prev * 1000 / (duration + 1);
                    log::trace!("- searched {} spectra ({} spectra/s)", prev, rate);
                }
                x
            })
            .flat_map(|spec| scorer.score(spec))
            .collect();

        let duration = Instant::now().duration_since(start).as_millis() as usize;
        let prev = counter.load(Ordering::Relaxed);
        let rate = prev * 1000 / (duration + 1);
        log::info!("- search:  {:8} ms ({} spectra/s)", duration, rate);
        features
    }

    fn complete_features(
        &self,
        spectra: Vec<ProcessedSpectrum>,
        features: Vec<Feature>,
    ) -> SageResults {
        let quant = self
            .parameters
            .quant
            .tmt
            .as_ref()
            .map(|isobaric| {
                let level = self.parameters.quant.tmt_settings.level;
                if level != 2 && level != 3 {
                    log::warn!("TMT quant level set at {}, is this correct?", level);
                }
                sage_core::tmt::quantify(&spectra, isobaric, Tolerance::Ppm(-20.0, 20.0), level)
            })
            .unwrap_or_default();
        let ms1 = spectra.into_iter().filter(|s| s.level == 1).collect();

        SageResults {
            features,
            quant,
            ms1,
        }
    }

    fn process_chunk(
        &self,
        scorer: &Scorer,
        chunk: &[String],
        chunk_idx: usize,
        batch_size: usize,
    ) -> SageResults {
        let spectra = self.read_processed_spectra(chunk, chunk_idx, batch_size);
        let features = self.search_processed_spectra(scorer, &spectra);
        self.complete_features(spectra, features)
    }

    fn read_processed_spectra(
        &self,
        chunk: &[String],
        chunk_idx: usize,
        batch_size: usize,
    ) -> Vec<ProcessedSpectrum> {
        // Read all of the spectra at once - this can help prevent memory over-consumption issues
        info!(
            "processing files {} .. {} ",
            batch_size * chunk_idx,
            batch_size * chunk_idx + chunk.len()
        );
        let start = Instant::now();

        let sn = self
            .parameters
            .quant
            .tmt_settings
            .sn
            .then_some(self.parameters.quant.tmt_settings.level);

        let min_deisotope_mz = match &self.parameters.quant.tmt {
            Some(i) => match self.parameters.quant.tmt_settings.level {
                2 => i.reporter_masses().last().map(|x| x * (1.0 + 20E-6)),
                _ => None,
            },
            None => None,
        };

        let sp = SpectrumProcessor::new(
            self.parameters.max_peaks,
            self.parameters.deisotope,
            min_deisotope_mz.unwrap_or(0.0),
        );

        let bruker_extensions = [".d", ".tdf", ".tdf_bin", "ms2", "raw"];
        let spectra = chunk
            .par_iter()
            .enumerate()
            .flat_map(|(idx, path)| {
                let file_id = chunk_idx * batch_size + idx;

                let path_lower = path.to_lowercase();
                let res = if path_lower.ends_with(".mgf.gz") || path_lower.ends_with(".mgf") {
                    sage_cloudpath::util::read_mgf(path, file_id)
                } else if bruker_extensions.iter().any(|ext| {
                    if path_lower.ends_with(std::path::MAIN_SEPARATOR) {
                        path_lower
                            .strip_suffix(std::path::MAIN_SEPARATOR)
                            .unwrap()
                            .ends_with(ext)
                    } else {
                        path_lower.ends_with(ext)
                    }
                }) {
                    sage_cloudpath::util::read_tdf(
                        path,
                        file_id,
                        self.parameters.bruker_spectrum_processor,
                    )
                } else {
                    sage_cloudpath::util::read_mzml(path, file_id, sn)
                };

                match res {
                    Ok(s) => {
                        log::trace!("- {}: read {} spectra", path, s.len());
                        Ok(s)
                    }
                    Err(e) => {
                        log::error!("- {}: {}", path, e);
                        Err(e)
                    }
                }
            })
            .flat_map_iter(|spectra| spectra.into_iter().map(|s| sp.process(s)))
            .collect::<Vec<_>>();

        let io_time = Instant::now() - start;
        info!("- file IO: {:8} ms", io_time.as_millis());

        spectra
    }

    pub fn batch_files(&self, scorer: &Scorer, batch_size: usize) -> SageResults {
        self.parameters
            .mzml_paths
            .chunks(batch_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| self.process_chunk(scorer, chunk, chunk_idx, batch_size))
            .collect::<SageResults>()
    }

    pub fn run(mut self, parallel: usize, parquet: bool) -> anyhow::Result<telemetry::Telemetry> {
        let scorer = Scorer {
            db: &self.database,
            precursor_tol: self.parameters.precursor_tol,
            fragment_tol: self.parameters.fragment_tol,
            min_matched_peaks: self.parameters.min_matched_peaks,
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            min_precursor_charge: self.parameters.precursor_charge.0,
            max_precursor_charge: self.parameters.precursor_charge.1,
            override_precursor_charge: self.parameters.override_precursor_charge,
            max_fragment_charge: self.parameters.max_fragment_charge,
            chimera: self.parameters.chimera,
            report_psms: self.parameters.report_psms,
            wide_window: self.parameters.wide_window,
            annotate_matches: self.parameters.annotate_matches,
            score_type: self.parameters.score_type,
        };

        //Collect all results into a single container
        let mut outputs = self.batch_files(&scorer, parallel);

        let alignments = if self.parameters.predict_rt {
            // Poisson probability is usually the best single feature for refining FDR.
            // Take our set of 1% FDR filtered PSMs, and use them to train a linear
            // regression model for predicting retention time
            outputs
                .features
                .par_sort_unstable_by(|a, b| a.poisson.total_cmp(&b.poisson));
            sage_core::ml::qvalue::spectrum_q_value(&mut outputs.features);

            let alignments = sage_core::ml::retention_alignment::global_alignment(
                &mut outputs.features,
                self.parameters.mzml_paths.len(),
            );
            let _ = sage_core::ml::retention_model::predict(&self.database, &mut outputs.features);
            let _ = sage_core::ml::mobility_model::predict(&self.database, &mut outputs.features);
            Some(alignments)
        } else {
            None
        };

        let q_spectrum = self.spectrum_fdr(&mut outputs.features);
        let q_peptide = sage_core::fdr::picked_peptide(&self.database, &mut outputs.features);
        let q_protein = sage_core::fdr::picked_protein(&self.database, &mut outputs.features);

        let filenames = self
            .parameters
            .mzml_paths
            .iter()
            .map(|s| {
                s.parse::<CloudPath>()
                    .ok()
                    .and_then(|c| c.filename().map(|s| s.to_string()))
                    .unwrap_or_else(|| s.clone())
            })
            .collect::<Vec<_>>();

        let areas = alignments.and_then(|alignments| {
            if self.parameters.quant.lfq {
                let mut areas = sage_core::lfq::build_feature_map(
                    self.parameters.quant.lfq_settings,
                    self.parameters.precursor_charge,
                    &outputs.features,
                )
                .quantify(&self.database, &outputs.ms1, &alignments);

                let q_precursor = sage_core::fdr::picked_precursor(&mut areas);

                log::info!("discovered {} target MS1 peaks at 5% FDR", q_precursor);
                Some(areas)
            } else {
                None
            }
        });

        log::info!(
            "discovered {} target peptide-spectrum matches at 1% FDR",
            q_spectrum
        );
        log::info!("discovered {} target peptides at 1% FDR", q_peptide);
        log::info!("discovered {} target proteins at 1% FDR", q_protein);
        log::trace!("writing outputs");

        // Write either a single parquet file, or multiple tsv files
        if parquet {
            log::warn!("parquet output format is currently unstable! There may be failures or schema changes!");

            let bytes = sage_cloudpath::parquet::serialize_features(
                &outputs.features,
                &outputs.quant,
                &filenames,
                &self.database,
            )?;

            let path = self.make_path("results.sage.parquet");
            path.write_bytes_sync(bytes)?;
            self.parameters.output_paths.push(path.to_string());

            if self.parameters.annotate_matches {
                let bytes =
                    sage_cloudpath::parquet::serialize_matched_fragments(&outputs.features)?;
                let path = self.make_path("matched_fragments.sage.parquet");
                path.write_bytes_sync(bytes)?;
                self.parameters.output_paths.push(path.to_string());
            }

            if let Some(areas) = &areas {
                let bytes =
                    sage_cloudpath::parquet::serialize_lfq(areas, &filenames, &self.database)?;

                let path = self.make_path("lfq.parquet");
                path.write_bytes_sync(bytes)?;
                self.parameters.output_paths.push(path.to_string());
            }
        } else {
            self.parameters
                .output_paths
                .push(self.write_features(&outputs.features, &filenames)?);

            if self.parameters.annotate_matches {
                self.parameters
                    .output_paths
                    .push(self.write_fragments(&outputs.features)?);
            }

            if !outputs.quant.is_empty() {
                self.parameters
                    .output_paths
                    .push(self.write_tmt(&outputs.quant, &filenames)?);
            }
            if let Some(areas) = areas {
                self.parameters
                    .output_paths
                    .push(self.write_lfq(areas, &filenames)?);
            }
        }

        // Write percolator input file if requested
        if self.parameters.write_pin {
            self.parameters
                .output_paths
                .push(self.write_pin(&outputs.features, &filenames)?);
        }

        let path = self.make_path("results.json");
        self.parameters.output_paths.push(path.to_string());
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);

        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        path.write_bytes_sync(bytes)?;

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);
        info!("cite: \"Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale\" https://doi.org/10.1021/acs.jproteome.3c00486");

        let telemetry = telemetry::Telemetry::new(
            self.parameters,
            self.database.peptides.len(),
            self.database.fragments.len(),
            parquet,
            run_time,
        );

        Ok(telemetry)
    }
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .arg(
            Arg::new("parameters")
                .required(true)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to configuration parameters (JSON file)")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("mzml_paths")
                .num_args(1..)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Paths to mzML files to process. Overrides mzML files listed in the \
                     configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to FASTA database. Overrides the FASTA file \
                     specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output_directory")
                .short('o')
                .long("output_directory")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path where search and quant results will be written. \
                     Overrides the directory specified in the configuration file.",
                )
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("batch-size")
                .long("batch-size")
                .value_parser(value_parser!(u16).range(1..))
                .help("Number of files to load and search in parallel (default = # of CPUs/2)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("parquet")
                .long("parquet")
                .action(clap::ArgAction::SetTrue)
                .help("Write search output in parquet format instead of tsv"),
        )
        .arg(
            Arg::new("annotate-matches")
                .long("annotate-matches")
                .action(clap::ArgAction::SetTrue)
                .help("Write matched fragments output file."),
        )
        .arg(
            Arg::new("write-pin")
                .long("write-pin")
                .action(clap::ArgAction::SetTrue)
                .help("Write percolator-compatible `.pin` output files"),
        )
        .arg(
            Arg::new("disable-telemetry")
                .long("disable-telemetry-i-dont-want-to-improve-sage")
                .action(clap::ArgAction::SetFalse)
                .help("Disable sending telemetry data"),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    let parallel = matches
        .get_one::<u16>("batch-size")
        .copied()
        .unwrap_or_else(|| num_cpus::get() as u16 / 2) as usize;

    let parquet = matches.get_one::<bool>("parquet").copied().unwrap_or(false);
    let send_telemetry = matches
        .get_one::<bool>("disable-telemetry")
        .copied()
        .unwrap_or(true);

    let input = Input::from_arguments(matches)?;

    let runner = input
        .build()
        .and_then(|parameters| Runner::new(parameters, parallel))?;

    let tel = runner.run(parallel, parquet)?;

    if send_telemetry {
        tel.send();
    }

    Ok(())
}
