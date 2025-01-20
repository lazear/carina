# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.15.0-alpha] (unreleased)
### Added
- Initial support for LFQ on data with ion mobility.
- Speedup on the generation of databases when large number of peptides are redundant.
- Initial support for searching diaPASEF data
- `override_precursor_charge` setting that forces multiple charge states to be searched
### Breaking Changes
- `precursor_ppm` field reports the non-absoluted average mass error, rather than the absoluted average mass error.
- Don't deisotope reporter ion regions if MS2-based TMT/iTRAQ is used
- Removed `fragment_min_mz` and `fragment_max_mz` parameters. These were decreasing the accuracy of preliminary scoring estimation when attempting to annotate multiply-charged, high-m/z ions.

## [v0.14.7]
### Added
- Added columns missing from parquet output: `semi_enzymatic` and `missed_cleavages`
### Changed
- Fixed ion mobility parsing from some mzMLs
- MGF paths were being lowercased prior to parsing

## [v0.14.6]
### Added
- Support for MGF files
- Support for writing ion mobility measurements to output files: `ion_mobility`, `predicted_mobility`, `delta_mobility` added to primary tsv and parquet reports. Ion mobility is predicted in a similar manner to RT, using a linear model trained on the data from the search.

## [v0.14.5]
### Added
- Support for semi-enzymatic digests (`database.enzyme.semi_enzymatic` parameter)
- Ability to directly export matched fragment ions (e.g. for spectral library or rescoring) with the `--annotate-matches` CLI option. This is compatible with the `--parquet` CLI option as well. Annotations will be written to `matched_fragments.sage.tsv` or `matched_fragments.sage.parquet`
- Sage sends basic telemetry data (version of Sage, run time, OS, # of CPU cores, # of peptides in database, whether LFQ is used) to a remote server. No information about your actual data is sent - e.g. identifications, quantities, organism, or modifications are NOT tracked or reported.  This data will be used to help focus efforts on improving Sage and figuring which features are most used. Please take a look at `crates/sage-cli/src/telemetry.rs` to see exactly what is sent! You can disable sending telemetry data  by using the `--disable-telemetry-i-dont-want-to-improve-sage` CLI flag.
### Changed
- Modified visibility on some crate internals to support the [sagepy project](https://github.com/theGreatHerrLebert/sagepy)
- Added `psm_id` field to various output files to match the new `--annotate-matches` option.
### Removed
- Removed the `ms1_intensity` field from CSV output, since it is essentially useless


## [v0.14.4]
### Added
- **Unstable feature**: Preliminary support for reading Bruker .d folders (ddaPASEF; no MS1/LFQ support yet)
### Changed
- Retention times are converted to minutes
### Fixed
- Fixed bug where charge state 1 would never be searched

## [v0.14.3]
### Fixed
- Hotfix for bug in parquet LFQ writer

## [v0.14.2]
### Added
- `quant.lfq_settings.combine_charge_state` boolean option. By default this is set to `true`, and LFQ is performed on the peptide-level, where all charge states are treated as the same precursor. Setting this to `false` performs LFQ on the peptide-charge-level, where each charge state will be treated separately.
### Changed
- Percolator output format now contains the integer-valued charge state encoded in the `z=other` column, if the charge state is outside the range 2-6 (e.g. a value of 7 will appear in the `z=other` column, rather than it being one-hot encoded)
- LFQ uses the the charge state range from the `precursor_charge` configuration option for tracing MS1 peaks

## [v0.14.1]
### Added
- Added additional output showing search progress if `SAGE_LOG=trace` environment variable is set
- Added additional warnings about precursor tolerances
- Added configuration option `precursor_charge` to make it explicit what charge states are being searched in the case where the mzML does not contain charge state information, or where `wide_window` is turned on.
### Changed
- Added a warning message if variable modifications are specified as single values (e.g. `15.9949`) instead of lists of values (e.g. `[15.9949]`). By v0.15 this will become a hard error and will not parse, to simply some of the internal logic.

## [v0.14.0]
### Added
- Support for parquet file format output. Search results and reporter ion quantification will be written to one file (`results.sage.parquet`) and label-free quant will be written to another (`lfq.parquet`). Parquet files tend to be significantly smaller than TSV files, faster to parse, and are compatible with a variety of distributed SQL engines.
### Changed
- Implement heapselect algorithm for faster sorting of candidate matches (#80). This is a backwards-incompatible change with respect to output - small changes in PSM ranks will be present between v0.13.4 and v0.14.0

## [v0.13.4]
### Fixed
- Bug in mzML parser, where some older specification-compliant mzMLs would not parse. If your mzMLs previously parsed, then there will be no change in behavior. Added a test case

## [v0.13.3]
### Fixed
- Bug in `database.enzyme.restrict` parameter, where `null` values were being overriden with "P" (causing Trypsin/P to behave like Trypsin)

## [v0.13.2]
### Changed
- Subtle change to TMT integration tolerance, and selection of which ion to quantify (most intense). As a result, TMT integration should be more in agreement (if not 100% so) with ProteomeDiscover/FragPipe/etc
- Remove `delta_mass` (precursor ppm) LDA feature - instead, build a delta mass (or ppm) profile using KDE/posterior error calculation code, and use the P(decoy) as a feature for LDA.

## [v0.13.1]
### Changed
- Internal performance and stability improvements for RT prediction & LDA

## [v0.13.0]
### Added
- Better error reporting thanks to @Elendol
- Added support for multiple variable mods for the same amino acid
- Added support for N/C-terminal modifications specific to an individual amino acid

New syntax:
```json
"variable_mods": {
    "M": [15.9949],
    "^Q": -17.026549,
    "^E": -18.010565,
    "[": 42.010565
}
```

Either a single floating point number (-18.0) or a list of floating point numbers ([-18.0, -15.2]) can be supplied as modifications. Support for single values may eventually be phased out to simplify the parser.

### Changed
- Changed "_fdr" columns to "_q" (e.g. "spectrum_q") in "results.sage.tsv" file
- Changed internal data representation of `Peptide` struct to allow for sharing of sequences (using `Arc`) among modified peptides
- Fragment index creation should now be faster

## [0.12.0]
### Added
- Add `wide_window` option to configuration file. This option turns off `precursor_tol`, instead using the isolation window written in the mzML file.
### Changed
- Changed internal calculation of precursor tolerances when searching with `isotope_errors`. The new version should be more accurate. This change also enables a significant boost to search speed for open searches.

## [0.11.2]
### Added
- Add rank & charge features to LDA
### Changed
- One-hot encode charge state information for percolator `.pin` files
- Change PSMId -> SpecId for Mokapot compatibility with `.pin` files

## [0.11.1]
### Added
- Support for additional fragment ion types, via the "database.ion_kinds" configuration option. Valid values are "a", "b", "c", "x", "y", "z"
### Changed
- Sort protein names alphanumerically for each peptide entry. This should enhance stability across runs, and fixes a bug with picked-protein group FDR
- Fix another bug where picked-FDR approaches assume internal decoy generation

### Changed
- Modify order of operations during deisotoping. Deisotoped peaks can contribute intensity to only 1 parent peak now, rather than potentially multiple parent peaks

## [0.11.0]
### Added
- Support for percolator output files (`--write-pin` CLI flag)
- Support for modifying file batch size (`--batch-size N` CLI flag)
- Add `delta_best` feature, which reports the delta hyperscore from the best match to current ranked PSM
- Add Sage version to `results.json` files

### Changed
- Breaking changes to `quant` section of the configuration file format
- Rename `delta_hyperscore` to `delta_next`
- Altered internal scoring algorithm. Rather than consider all MS2 peaks within a m/z tolerance window to be matches to a theoretical spectrum, consider only the closest peak. This should increase the accuracy of # of matched peaks, and subsequent scores
- Overhaul of chimeric scoring, `report_psms` can now be used to search for multiple chimeric spectra
- Completely overhauled the LFQ algorithm: added match-between runs, peak scoring using normalized spectral angle relative to theoretical isotopic envelope, target decoy scoring of MS1 integration
- Fixed bug in picked-peptide FDR that could lead to liberal FDR
- Fixed bug in picked-protein FDR that could lead to conservative FDR
- Fixed bug where using variable protein terminal (e.g. protein N-terminal acetylation) modifications could cause some determinism. This also improves the accuracy of peptide => protein assignment. Unfortunately this fix has performance implications, causing creation of the fragment index to take up to ~2x as long.

### Removed
- Remove `no-parallel` CLI flag, and `parallel` configuration file entry

## [0.10.0]
### Added
- Retention times are now globally aligned across files
- RT prediction is then performed on all files at once (on aligned RTs), rather than one file at a time - previously, there were many instances where some files in a search could not have RTs predicted, decreasing the effectiveness of delta_rt as a feature for LDA.

### Changed
- Peptide sequences within a protein are now deduplicated - previously, repeated peptides would be called multiple times for the same protein (e.g. num_proteins > 1 even if the peptide was unique)

## [0.9.4]
### Changed
- Fix issues with RT prediction (and occasionally LDA) that arise from 0's being present on the diagonals of the covariance matrix (small amount of regularization added)

## [0.9.3]
### Added
- Allow users to set minimum number of matched b+y ions for reporting PSMs (`min_matched_peaks`)

### Changed
- Internal code for calculating factorials

## [0.9.2]
### Added
- Added option for TMT signal/noise quantification, if noise values are present in mzML

## [0.9.1]
### Changed
- FASTA file path, JSON configuration file can now be specified as "s3://" paths, allowing Sage to run completely disk-free

## [0.9.0]
### Added
- Support for non-specific digests, N-terminal enzymatic digestion

## [0.8.1]
### Added
- `quant.tmt_level` configuration option to enable MS2 (or MSn) isobaric quantification

## [0.8.0]
### Added
- Support for protein N-terminal ('['), C-terminal (']') as well as peptide C-terminal ('$') modifications
- Support for k-combinations of variable modifications. This can be specified with the `database.max_variable_mods` parameter

## [0.7.1] - 2022-11-04
### Changed
- Fix bug with in silico digest: Logic around overwriting decoys with target sequences was incorrect peptides shared between targets/decoys were being annotated as decoy peptides but assigned to non-decoy proteins. We now make sure that they are assigned to non-decoy proteins and also annotated as target sequences.

## [0.7.0] - 2022-11-03
### Added
- Add support for user-specified enzymes to JSON file.  `database.enzyme.sites` and `database.enzyme.restrict` are limited to valid amino acids
- Sage can now search MS2 spectra without annotated precursor charge states. Default behavior is to search with z=2, z=3, z=4, and then merge the PSMs for scoring

### Changed
- Configuration file schema changed. `peptide_min_len`, `peptide_max_len`, `missed_cleavages` are now specified under `database.enzyme` in the JSON file
- Internal behavior of Sage was changed to enable deterministic searching
- Docker file changed from Alpine to Debian


## [0.6.0] - 2022-11-01
### Added
- Changelog
- `rank` column added to output file
- `database.generate_decoys` parameter, which turns off internal decoy generation. This enables the use of FASTA databases for SearchGUI/PeptideShaker

### Changed
- Base ProForma v2 notation is used for peptide modifications, i.e. "\[+304.2071\]-PEPTIDEM\[+15.9949\]AAC\[+57.0214\]H"
- `scannr` column now contains the full nativeID/spectrum title from the mzML file, i.e. "controllerType=0 controllerNumber=1 scan=30069"
- `discriminant_score` column renamed to `sage_discriminant_score` for PeptideShaker recognition
- `database.decoy_prefix` JSON option changed to `database.decoy_tag`. This allows decoy tagging to occur anywhere within the accession: "sp|P01234_REVERSED|HUMAN"
- Output file renamed:  `results.pin` to `results.sage.tsv`
- Output file renamed: `quant.csv` to `quant.tsv`
- Rename `pin_paths` to `output_paths` in results.json file


## [0.5.1] - 2022-10-31
### Added
- Support for selenocysteine and pyrrolysine amino acids

## [0.5.0] - 2022-10-28
### Added
- Ability to directly read/write files from AWS S3

### Changed
- Processing files in parallel processes them in batches of `num_cpus / 2` to avoid memory issues
- Fixed bug where `protein_fdr` was erroneously assigned to `peptide_fdr` output field
- Additional parallelization for assignment of PEP, FDR, writing output files

## [0.4.0] - 2022-10-18
### Added
- Label free quantification can be enabled by turning on `quant.lfq` JSON parameter 
- Commmand line arguments can be used to override configuration file

## [0.3.1] - 2022-10-06
### Added
- Workflow contributions from [@wfondrie](https://github.com/wfondrie).

### Changed
- Don't parse empty MS2 spectra

## [0.3.0] - 2015-09-15
### Added
- Retention time prediction
- Ability to filter low-number b/y-ions for faster preliminary scoring (`database.min_ion_index` option)
- Ability to toggle retention time prediction (`predict_rt`)
