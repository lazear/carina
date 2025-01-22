//! TMT quantification
#![allow(clippy::excessive_precision)]
#![allow(unused_imports)]
use crate::database::binary_search_slice;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, H2O, NH3, PROTON};
use crate::peptide::Peptide;
use crate::scoring::{Feature, Scorer};
use crate::spectrum::{self, Peak, Precursor, ProcessedSpectrum};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum Isobaric {
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
    User(Vec<f32>),
}

impl Isobaric {
    /// Return the monoisotopic mass of reporter ions
    pub fn reporter_masses(&self) -> &[f32] {
        match self {
            Isobaric::Tmt6 => &TMT6PLEX,
            Isobaric::Tmt10 => &TMT11PLEX[0..10],
            Isobaric::Tmt11 => &TMT11PLEX,
            Isobaric::Tmt16 => &TMT18PLEX[0..16],
            Isobaric::Tmt18 => &TMT18PLEX,
            Isobaric::User(labels) => labels,
        }
    }

    /// Return the monoisotopic mass of tag
    pub fn modification_mass(&self) -> Option<f32> {
        match self {
            Isobaric::Tmt6 | Isobaric::Tmt10 | Isobaric::Tmt11 => Some(229.162932),
            Isobaric::Tmt16 => Some(304.2071),
            Isobaric::Tmt18 => Some(304.2135),
            Isobaric::User(_) => None,
        }
    }

    /// Return a column name for each tag
    pub fn headers(&self) -> Vec<String> {
        match self {
            Isobaric::User(v) => v
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("user_{}", idx + 1))
                .collect(),
            _ => self
                .reporter_masses()
                .iter()
                .enumerate()
                .map(|(idx, _)| format!("tmt_{}", idx + 1))
                .collect(),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Purity {
    pub ratio: f32,
    pub correct_precursors: usize,
    pub incorrect_precursors: usize,
}

/// Calculate SPS purity stats - which SPS precursor ions actually correspond
/// to theoretical b/y ions, and percentage of total SPS precursor MS2 intensity
/// that is explained by theoretical b/y ions
// fn purity_of_match(
//     precursors: &[Precursor],
//     ms2: &ProcessedSpectrum,
//     theoretical_peaks: &[Peak],
//     max_charge: u8,
//     fragment_tolerance: Tolerance,
// ) -> Purity {
//     let mut explained_intensity = 0.0;
//     let mut interference = 0.0;
//     let mut correct_precursors = 0;
//     let mut incorrect_precursors = 0;

//     for precursor in precursors {
//         // Spurious SPS precursor introduced by MSConvert
//         // https://github.com/ProteoWizard/pwiz/issues/2202
//         if precursor.scan.unwrap_or(0) != ms2.scan {
//             continue;
//         }

//         // This is really a m/z window, we haven't performed charge state deconvolution!
//         // Select a window of MS2 peaks that have been sampled for MS3
//         let isolation_tolerance = precursor
//             .isolation_window
//             .unwrap_or_else(|| Tolerance::Da(-0.5, 0.5));
//         let isolation_window = isolation_tolerance.bounds(precursor.mz - PROTON);

//         let (idx_lo, idx_hi) = binary_search_slice(
//             &ms2.peaks,
//             |peak, key| peak.mass.total_cmp(key),
//             isolation_window.0,
//             isolation_window.1,
//         );

//         // Create slice surrounding SPS peaks selected
//         let window = &ms2.peaks[idx_lo..idx_hi];
//         interference += window
//             .iter()
//             .filter(|peak| peak.mass >= isolation_window.0 && peak.mass <= isolation_window.1)
//             .map(|peak| peak.intensity)
//             .sum::<f32>();

//         // Pick the most intense peak within the isolation window, and decide whether it is
//         // assignable to the best candidate peptide
//         let assigned_to_candidate = if let Some(best_peak) =
//             spectrum::select_closest_peak(window, precursor.mz - PROTON, isolation_tolerance)
//         {
//             let mut correct = false;
//             'outer: for charge in 1..max_charge {
//                 for loss in [0.0, NH3, H2O] {
//                     let mass = best_peak.mass * charge as f32 + loss;
//                     if spectrum::select_closest_peak(theoretical_peaks, mass, fragment_tolerance)
//                         .is_some()
//                     {
//                         correct = true;
//                         interference -= best_peak.intensity;
//                         explained_intensity += best_peak.intensity;

//                         break 'outer;
//                     }
//                 }
//             }
//             correct
//         } else {
//             false
//         };

//         if assigned_to_candidate {
//             correct_precursors += 1;
//         } else {
//             incorrect_precursors += 1;
//         }
//     }

//     Purity {
//         ratio: explained_intensity / (explained_intensity + interference),
//         correct_precursors,
//         incorrect_precursors,
//     }
// }

// fn mk_theoretical(peptide: &Peptide) -> Vec<Peak> {
//     let mut theoretical_peaks = IonSeries::new(peptide, Kind::B)
//         .map(|ion| Peak {
//             mass: ion.monoisotopic_mass,
//             intensity: 0.0,
//         })
//         .collect::<Vec<_>>();

//     if let Some(crate::mass::Residue::Mod('K', _)) = peptide.sequence.last() {
//         theoretical_peaks.extend(IonSeries::new(peptide, Kind::Y).map(|ion| Peak {
//             mass: ion.monoisotopic_mass,
//             intensity: 0.0,
//         }));
//     }

//     theoretical_peaks.sort_unstable_by(|a, b| a.mass.total_cmp(&b.mass));
//     theoretical_peaks
// }

#[derive(Debug)]
pub struct Quant<'ms3> {
    /// Top hit for this MS3 spectrum
    pub hit: Feature,
    /// Top chimeric/co-fragmenting hit for this spectrum
    pub chimera: Option<Feature>,
    /// SPS precursor purity for the top hit
    pub hit_purity: Purity,
    /// SPS precursor purity for the chimeric hit
    pub chimera_purity: Option<Purity>,
    /// Quanitified TMT reporter ion intensities
    pub intensities: Vec<Option<&'ms3 Peak>>,
    /// MS3 spectrum
    pub spectrum: &'ms3 ProcessedSpectrum<Peak>,
}

/// Return a vector containing the peaks closest to the m/zs defined in
/// `labels`, within a given tolerance window.
/// This function is MS-level agnostic, so it can be used for either MS2 or MS3
/// quant.
pub fn find_reporter_ions<'a>(
    peaks: &'a [Peak],
    labels: &[f32],
    label_tolerance: Tolerance,
) -> Vec<Option<&'a Peak>> {
    labels
        .iter()
        .map(|&label| {
            spectrum::select_most_intense_peak(peaks, label, label_tolerance, Some(-PROTON))
        })
        .collect()
}

const TMT6PLEX: [f32; 6] = [
    126.127726, 127.124761, 128.134436, 129.131471, 130.141145, 131.138180,
];

const TMT11PLEX: [f32; 11] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144499,
];

const TMT18PLEX: [f32; 18] = [
    126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
    130.141145, 131.138180, 131.144500, 132.141535, 132.147855, 133.144890, 133.151210, 134.148245,
    134.154565, 135.15160,
];

/// Search MS/MS and quantify isobaric tag intensities from an SPS-MS3 spectrum
///
/// * `scorer`: used for searching/scoring precursor MS2 spectrum
/// * `spectra`: a slice (generally entire mzML) of spectra, that can be searched
///     for precursor spectra
/// * `ms3`: The MS3 spectrum to search and quantify
/// * `isobaric_labels`: specify label m/zs to be used
/// * `isobaric_tolerance`: specify label tolerance
// pub fn quantify_sps<'a, 'b>(
//     scorer: &'a Scorer<'a>,
//     spectra: &[ProcessedSpectrum],
//     ms3: &'b ProcessedSpectrum,
//     isobaric_labels: &Isobaric,
//     isobaric_tolerance: Tolerance,
// ) -> Option<Quant<'b>> {
//     let first_precursor = ms3
//         .precursors
//         .first()
//         .expect("MS3 scan without at least one precursor!");

//     let ms2 = spectrum::find_spectrum_by_id(
//         spectra,
//         first_precursor
//             .scan
//             .expect("MS3 scan without a MS2 precursor scan ID"),
//     )
//     .expect("Couldn't locate parent MS2 scan!");

//     let ms1_charge = ms2
//         .precursors
//         .get(0)
//         .and_then(|p| p.charge)
//         .unwrap_or(2)
//         .saturating_sub(1);

//     let scores = scorer.score_chimera(ms2);
//     let hit = scores.first()?.clone();
//     let peptide = &scorer.db[hit.peptide_idx];
//     let hit_purity = purity_of_match(
//         &ms3.precursors,
//         ms2,
//         &mk_theoretical(peptide),
//         ms1_charge,
//         scorer.fragment_tol,
//     );

//     let chimera = scores.get(1).cloned();
//     let chimera_purity = chimera.as_ref().map(|score| {
//         purity_of_match(
//             &ms3.precursors,
//             ms2,
//             &mk_theoretical(&scorer.db[score.peptide_idx]),
//             ms1_charge,
//             scorer.fragment_tol,
//         )
//     });

//     Some(Quant {
//         hit,
//         hit_purity,
//         chimera,
//         chimera_purity,
//         intensities: find_reporter_ions(
//             &ms3.peaks,
//             isobaric_labels.reporter_masses(),
//             isobaric_tolerance,
//         ),
//         spectrum: ms3,
//     })
// }

#[derive(Clone)]
pub struct TmtQuant {
    pub spec_id: String,
    pub file_id: usize,
    pub ion_injection_time: f32,
    pub peaks: Vec<f32>,
}

/// Quantify isobaric tags from an MS2 or MS3 spectrum
///
/// * `spectra`: a slice (generally entire mzML) of spectra, that can be searched
///     for precursor spectra
/// * `isobaric_labels`: specify label m/zs to be used
/// * `isobaric_tolerance`: specify label tolerance
/// * `level`: MSn level to extract isobaric peaks from
pub fn quantify(
    spectra: &[ProcessedSpectrum<Peak>],
    isobaric_labels: &Isobaric,
    isobaric_tolerance: Tolerance,
    level: u8,
) -> Vec<TmtQuant> {
    spectra
        .par_iter()
        .filter(|spectrum| spectrum.level == level)
        .filter_map(|spectrum| {
            let spec_id = match level {
                1 => return None,
                2 => spectrum.id.clone(),
                _ => spectrum
                    .precursors
                    .first()
                    .and_then(|precursor| precursor.spectrum_ref.clone())
                    .unwrap_or_default(),
            };

            let peaks = find_reporter_ions(
                &spectrum.peaks,
                isobaric_labels.reporter_masses(),
                isobaric_tolerance,
            )
            .into_iter()
            .map(|peak| peak.map(|p| p.intensity).unwrap_or_default())
            .collect();

            Some(TmtQuant {
                spec_id,
                file_id: spectrum.file_id,
                ion_injection_time: spectrum.ion_injection_time,
                peaks,
            })
        })
        .collect()
}
