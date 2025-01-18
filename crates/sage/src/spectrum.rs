use crate::database::binary_search_slice;
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::spectrum;

/// A charge-less peak at monoisotopic mass
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct Peak {
    pub intensity: f32,
    pub mass: f32,
}

impl Eq for Peak {}

impl PartialOrd for Peak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Peak {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.intensity
            .total_cmp(&other.intensity)
            .then_with(|| self.mass.total_cmp(&other.mass))
    }
}

/// A charge-less peak at monoisotopic mass with ion mobility
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct IMPeak {
    pub intensity: f32,
    pub mass: f32,
    pub mobility: f32, // I would use f16 but its not stable
}

impl Eq for IMPeak {}

impl PartialOrd for IMPeak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IMPeak {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.intensity
            .total_cmp(&other.intensity)
            .then_with(|| self.mass.total_cmp(&other.mass))
            .then_with(|| self.mobility.total_cmp(&other.mobility))
    }
}

/// A de-isotoped peak, that might have some charge state information
#[derive(PartialEq, PartialOrd, Debug, Copy, Clone)]
pub struct Deisotoped {
    pub mz: f32,
    // Cumulative intensity of all isotopic peaks in the envelope higher than this one
    pub intensity: f32,
    // Assigned charge
    pub charge: Option<u8>,
    // If `Some(idx)`, idx is the index of the parent isotopic envelope
    pub envelope: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct SpectrumProcessor {
    pub take_top_n: usize,
    pub min_deisotope_mz: f32,
    pub deisotope: bool,
}

#[derive(Default, Debug, Clone)]
pub struct Precursor {
    pub mz: f32,
    pub intensity: Option<f32>,
    pub charge: Option<u8>,
    // pub scan: Option<usize>,
    pub spectrum_ref: Option<String>,
    pub isolation_window: Option<Tolerance>,
    pub inverse_ion_mobility: Option<f32>,
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum<T> {
    /// MSn level
    pub level: u8,
    /// Scan ID
    pub id: String,
    /// File ID
    pub file_id: usize,
    /// Retention time in minutes
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Selected ions for precursors, if `level > 1`
    pub precursors: Vec<Precursor>,
    /// MS peaks, sorted by mass in ascending order
    pub peaks: Vec<T>,
    /// Total ion current
    pub total_ion_current: f32,
}

#[derive(Default, Debug, Clone)]
/// An unprocessed mass spectrum, as returned by a parser
/// *CRITICAL*: Users must set all fields manually, including `file_id`
pub struct RawSpectrum {
    pub file_id: usize,
    /// MSn level
    pub ms_level: u8,
    /// Spectrum identifier
    pub id: String,
    /// Vector of precursors associated with this spectrum
    pub precursors: Vec<Precursor>,
    /// Profile or Centroided data
    pub representation: Representation,
    /// Scan start time in minutes
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Total ion current
    pub total_ion_current: f32,
    /// M/z array
    pub mz: Vec<f32>,
    /// Intensity array
    pub intensity: Vec<f32>,
    /// Mobility array
    pub mobility: Option<Vec<f32>>,
}

impl RawSpectrum {
    /// Return a [`RawSpectrum`] with default values, but with the `file_id` field
    /// properly set
    pub fn default_with_file_id(file_id: usize) -> Self {
        Self {
            file_id,
            ..Default::default()
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Representation {
    Profile,
    Centroid,
}

impl Default for Representation {
    fn default() -> Self {
        Self::Profile
    }
}

pub enum MS1Spectra {
    NoMobility(Vec<ProcessedSpectrum<Peak>>),
    WithMobility(Vec<ProcessedSpectrum<IMPeak>>),
    Empty,
}

impl Default for MS1Spectra {
    fn default() -> Self {
        Self::Empty
    }
}

/// Binary search followed by linear search to select the most intense peak within `tolerance` window
/// * `offset` - this parameter allows for a static adjustment to the lower and upper bounds of the search window.
///     Sage subtracts a proton (and assumes z=1) for all experimental peaks, and stores all fragments as monoisotopic
///     masses. This simplifies downstream calculations at multiple charge states, but it also subtly changes tolerance
///     bounds. For most applications this is completely OK to ignore - however, for exact similarity of TMT reporter ion
///     measurements with ProteomeDiscoverer, FragPipe, etc, we need to account for this minor difference (which has an impact
///     perhaps 0.01% of the time)
pub fn select_most_intense_peak(
    peaks: &[Peak],
    center: f32,
    tolerance: Tolerance,
    offset: Option<f32>,
) -> Option<&Peak> {
    let (lo, hi) = tolerance.bounds(center);
    let (lo, hi) = (
        lo + offset.unwrap_or_default(),
        hi + offset.unwrap_or_default(),
    );

    let (i, j) = binary_search_slice(peaks, |peak, query| peak.mass.total_cmp(query), lo, hi);

    let mut best_peak = None;
    let mut max_int = 0.0;
    for peak in peaks[i..j]
        .iter()
        .filter(|peak| peak.mass >= lo && peak.mass <= hi)
    {
        if peak.intensity >= max_int {
            max_int = peak.intensity;
            best_peak = Some(peak);
        }
    }
    best_peak
}

// pub fn find_spectrum_by_id(
//     spectra: &[ProcessedSpectrum],
//     scan_id: usize,
// ) -> Option<&ProcessedSpectrum> {
//     // First try indexing by scan
//     if let Some(first) = spectra.get(scan_id.saturating_sub(1)) {
//         if first.scan == scan_id {
//             return Some(first);
//         }
//     }
//     // Fall back to binary search
//     let idx = spectra
//         .binary_search_by(|spec| spec.scan.cmp(&scan_id))
//         .ok()?;
//     spectra.get(idx)
// }

/// Deisotope a set of peaks by attempting to find C13 peaks under a given `ppm` tolerance
pub fn deisotope(
    mz: &[f32],
    int: &[f32],
    max_charge: u8,
    ppm: f32,
    min_mz: f32,
) -> Vec<Deisotoped> {
    let mut peaks = mz
        .iter()
        .zip(int.iter())
        .map(|(mz, int)| Deisotoped {
            mz: *mz,
            intensity: *int,
            envelope: None,
            charge: None,
        })
        .collect::<Vec<_>>();

    // Is the peak at index `i` an isotopic peak?
    for i in (0..mz.len()).rev() {
        // Two pointer approach, j is fast pointer
        let mut j = i.saturating_sub(1);
        while mz[i] - mz[j] <= NEUTRON + Tolerance::ppm_to_delta_mass(mz[i], ppm) && mz[j] >= min_mz
        {
            let delta = mz[i] - mz[j];
            let tol = Tolerance::ppm_to_delta_mass(mz[i], ppm);
            for charge in 1..=max_charge {
                let iso = NEUTRON / charge as f32;
                if (delta - iso).abs() <= tol && int[i] < int[j] {
                    // Make sure this peak isn't already part of an isotopic envelope
                    if let Some(existing) = peaks[i].charge {
                        if existing != charge {
                            continue;
                        }
                    }
                    peaks[j].intensity += peaks[i].intensity;
                    peaks[j].charge = Some(charge);
                    peaks[i].charge = Some(charge);
                    peaks[i].envelope = Some(j);
                }
            }
            j = j.saturating_sub(1);
            if j == 0 {
                break;
            }
        }
    }
    peaks
}

/// Path compression of isotopic envelope links
pub fn path_compression(peaks: &mut [Deisotoped]) {
    for idx in 0..peaks.len() {
        if let Some(parent) = peaks[idx].envelope {
            if let Some(upper) = peaks[parent].envelope {
                peaks[idx].envelope = Some(upper);
            }
            peaks[idx].intensity = 0.0;
        }
    }
}

impl<T> ProcessedSpectrum<T> {
    pub fn extract_ms1_precursor(&self) -> Option<(f32, u8)> {
        let precursor = self.precursors.first()?;
        let charge = precursor.charge?;
        let mass = (precursor.mz - PROTON) * charge as f32;
        Some((mass, charge))
    }
    pub fn in_isolation_window(&self, mz: f32) -> Option<bool> {
        let precursor = self.precursors.first()?;
        let (lo, hi) = precursor.isolation_window?.bounds(precursor.mz - PROTON);
        Some(mz >= lo && mz <= hi)
    }
}

impl SpectrumProcessor {
    /// Create a new [`SpectrumProcessor`]
    ///
    /// # Arguments
    /// * `take_top_n`: Keep only the top N most intense peaks from the spectrum
    /// * `min_fragment_mz`: Keep only fragments >= this m/z
    /// * `max_fragment_mz`: Keep only fragments <= this m/z
    /// * `deisotope`: Perform deisotoping & charge state deconvolution
    pub fn new(take_top_n: usize, deisotope: bool, min_deisotope_mz: f32) -> Self {
        Self {
            take_top_n,
            min_deisotope_mz,
            deisotope,
        }
    }

    fn process_ms2(&self, should_deisotope: bool, spectrum: &RawSpectrum) -> Vec<Peak> {
        if spectrum.representation != Representation::Centroid {
            // Panic, because there's really nothing we can do with profile data
            panic!(
                "Scan {} contains profile data! Please convert to centroid",
                spectrum.id
            );
        }

        // If there is no precursor charge from the mzML file, then deisotope fragments up to z=3
        let charge = spectrum
            .precursors
            .first()
            .and_then(|p| p.charge)
            .unwrap_or(3);

        if should_deisotope {
            let mut peaks = deisotope(
                &spectrum.mz,
                &spectrum.intensity,
                charge,
                10.0,
                self.min_deisotope_mz,
            );
            peaks.sort_unstable_by(|a, b| {
                b.intensity
                    .total_cmp(&a.intensity)
                    .then_with(|| a.mz.total_cmp(&b.mz))
            });

            peaks
                .into_iter()
                .filter(|peak| peak.envelope.is_none())
                .map(|peak| {
                    // Convert from MH* to M
                    let mass = (peak.mz - PROTON) * peak.charge.unwrap_or(1) as f32;
                    Peak {
                        mass,
                        intensity: peak.intensity,
                    }
                })
                .take(self.take_top_n)
                .collect::<Vec<Peak>>()
        } else {
            let mut peaks = spectrum
                .mz
                .iter()
                .zip(spectrum.intensity.iter())
                .map(|(mz, &intensity)| {
                    let mass = (mz - PROTON) * 1.0;
                    Peak { mass, intensity }
                })
                .collect::<Vec<_>>();
            crate::heap::bounded_min_heapify(&mut peaks, self.take_top_n);
            peaks.truncate(self.take_top_n);
            peaks
        }
    }

    pub fn process(&self, spectrum: RawSpectrum) -> ProcessedSpectrum<Peak> {
        let mut peaks = match spectrum.ms_level {
            2 => self.process_ms2(self.deisotope, &spectrum),
            _ => spectrum
                .mz
                .iter()
                .zip(spectrum.intensity.iter())
                .map(|(&mass, &intensity)| {
                    let mass = (mass - PROTON) * 1.0;
                    Peak { mass, intensity }
                })
                .collect::<Vec<_>>(),
        };

        peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
        let total_ion_current = peaks.iter().map(|peak| peak.intensity).sum::<f32>();

        ProcessedSpectrum {
            level: spectrum.ms_level,
            id: spectrum.id,
            file_id: spectrum.file_id,
            scan_start_time: spectrum.scan_start_time,
            ion_injection_time: spectrum.ion_injection_time,
            precursors: spectrum.precursors,
            peaks,
            total_ion_current,
        }
    }

    pub fn process_with_mobility(&self, spectrum: RawSpectrum) -> ProcessedSpectrum<IMPeak> {
        assert!(
            spectrum.ms_level == 1,
            "Logic error, mobility processing should only be used for MS1"
        );
        let mut peaks = spectrum
            .mz
            .iter()
            .zip(
                spectrum
                    .intensity
                    .iter()
                    .zip(spectrum.mobility.unwrap().iter()),
            )
            .map(|(&mass, (&intensity, &mobility))| {
                let mass = (mass - PROTON) * 1.0;
                IMPeak {
                    mass,
                    intensity,
                    mobility,
                }
            })
            .collect::<Vec<_>>();

        peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
        let total_ion_current = peaks.iter().map(|peak| peak.intensity).sum::<f32>();

        ProcessedSpectrum {
            level: spectrum.ms_level,
            id: spectrum.id,
            file_id: spectrum.file_id,
            scan_start_time: spectrum.scan_start_time,
            ion_injection_time: spectrum.ion_injection_time,
            precursors: spectrum.precursors,
            peaks,
            total_ion_current,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_deisotope() {
        let mut mz = [
            800.9,
            800.9 + NEUTRON * 1.0,
            800.9 + NEUTRON * 2.0,
            803.4080,
            804.4108,
            805.4106,
            806.4116,
            810.0,
            812.0,
            812.0 + NEUTRON / 2.0,
        ];
        let mut int = [2., 1.5, 1., 4., 3., 2., 1., 1., 9.0, 4.5];
        let mut peaks = deisotope(&mut mz, &mut int, 2, 5.0, 800.91);

        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 2.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 1.0,
                    intensity: 2.5,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 2.0,
                    intensity: 1.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 6.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 3.0,
                    charge: Some(1),
                    envelope: Some(4),
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 1.0,
                    charge: Some(1),
                    envelope: Some(5),
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 4.5,
                    charge: Some(2),
                    envelope: Some(8),
                }
            ]
        );

        path_compression(&mut peaks);
        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 2.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 1.0,
                    intensity: 2.5,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 2.0,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 0.0,
                    charge: Some(2),
                    envelope: Some(8),
                }
            ]
        );
    }
}
