use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use timsrust::converters::ConvertableDomain;
pub use timsrust::readers::SpectrumReaderConfig as BrukerSpectrumProcessor;

pub struct TdfReader;

use std::cmp::Ordering;

fn squash_frame(mz_array: &[f32], intensity_array: &[f32], tol_ppm: f32) -> (Vec<f32>, Vec<f32>) {
    // Make sure the mz array is sorted
    assert!(mz_array.windows(2).all(|x| x[0] <= x[1]));

    let arr_len = mz_array.len();
    let mut touched = vec![false; arr_len];
    let mut global_num_touched = 0;

    let mut order: Vec<usize> = (0..arr_len).collect();
    order.sort_unstable_by(|&a, &b| {
        intensity_array[b]
            .partial_cmp(&intensity_array[a])
            .unwrap_or(Ordering::Equal)
    });

    let mut agg_mz = vec![0.0; arr_len];
    let mut agg_intensity = vec![0.0; arr_len];

    let utol = tol_ppm / 1e6;

    for &idx in &order {
        if touched[idx] {
            continue;
        }

        let mz = mz_array[idx];
        let da_tol = mz * utol;
        let left_e = mz - da_tol;
        let right_e = mz + da_tol;

        let ss_start = mz_array.partition_point(|&x| x < left_e);
        let ss_end = mz_array.partition_point(|&x| x <= right_e);

        let slice_width = ss_end - ss_start;
        let local_num_touched = touched[ss_start..ss_end].iter().filter(|&&x| x).count();
        let local_num_untouched = slice_width - local_num_touched;

        if local_num_touched == slice_width {
            continue;
        }

        let mut curr_intensity = 0.0;
        let mut curr_weighted_mz = 0.0;

        for i in ss_start..ss_end {
            if !touched[i] && intensity_array[i] > 0.0 {
                curr_intensity += intensity_array[i];
                curr_weighted_mz += mz_array[i] * intensity_array[i];
            }
        }

        if curr_intensity > 0.0 {
            curr_weighted_mz /= curr_intensity;

            agg_intensity[idx] = curr_intensity;
            agg_mz[idx] = curr_weighted_mz;

            touched[ss_start..ss_end].iter_mut().for_each(|x| *x = true);
            global_num_touched += local_num_untouched;
        }

        if global_num_touched == arr_len {
            break;
        }
    }

    // Drop the zeros and sort
    let mut result: Vec<(f32, f32)> = agg_mz
        .into_iter()
        .zip(agg_intensity.into_iter())
        .filter(|&(mz, intensity)| mz > 0.0 && intensity > 0.0)
        .collect();

    result.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    result.into_iter().unzip()
}

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<str>,
        file_id: usize,
        bruker_spectrum_processor: BrukerSpectrumProcessor,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let tdf_path = std::path::Path::new(path_name.as_ref()).join("analysis.tdf");
        let spectrum_reader = timsrust::readers::SpectrumReader::build()
            .with_path(path_name.as_ref())
            .with_config(bruker_spectrum_processor)
            .finalize()?;
        let frame_reader = timsrust::readers::FrameReader::new(path_name.as_ref())?;
        let metadata = timsrust::readers::MetadataReader::new(tdf_path)?;
        let mz_converter = metadata.mz_converter;

        // Get ms2 spectra
        let spectra: Vec<RawSpectrum> = (0..spectrum_reader.len())
            .into_par_iter()
            .filter_map(|index| match spectrum_reader.get(index) {
                Ok(dda_spectrum) => match dda_spectrum.precursor {
                    Some(dda_precursor) => {
                        let mut precursor = Self::parse_precursor(dda_precursor);
                        precursor.isolation_window = Option::from(Tolerance::Da(
                            -dda_spectrum.isolation_width as f32 / 2.0,
                            dda_spectrum.isolation_width as f32 / 2.0,
                        ));
                        let spectrum: RawSpectrum = RawSpectrum {
                            file_id,
                            precursors: vec![precursor],
                            representation: Representation::Centroid,
                            scan_start_time: dda_precursor.rt as f32 / 60.0,
                            // I am pretty sure this is wrong ...
                            // The retention time is most certainly not what we want for the
                            // injection time.
                            ion_injection_time: dda_precursor.rt as f32,
                            total_ion_current: 0.0,
                            mz: dda_spectrum.mz_values.iter().map(|&x| x as f32).collect(),
                            ms_level: 2,
                            id: dda_spectrum.index.to_string(),
                            intensity: dda_spectrum.intensities.iter().map(|&x| x as f32).collect(),
                        };
                        Some(spectrum)
                    }
                    None => {
                        log::warn!("No precursor found for spectrum {:?}", dda_spectrum.index);
                        None
                    }
                },
                // Q: should we raise/propagate/log an error here?
                Err(x) => {
                    log::error!("error parsing spectrum: {:?}", x);
                    None
                }
            })
            .collect();

        // Get MS1 spectra
        let ms1_spectra = frame_reader
            .get_all_ms1()
            .into_iter()
            .filter_map(|frame| match frame {
                Ok(frame) => {
                    let mz: Vec<f32> = frame
                        .tof_indices
                        .iter()
                        .map(|&x| mz_converter.convert(x as f64) as f32)
                        .collect();
                    let intensity: Vec<f32> = frame.intensities.iter().map(|&x| x as f32).collect();

                    // Sort the mzs and intensities by mz
                    let mut indices: Vec<usize> = (0..mz.len()).collect();
                    indices.sort_by(|&i, &j| {
                        mz[i]
                            .partial_cmp(&mz[j])
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let sorted_mz: Vec<f32> = indices.iter().map(|&i| mz[i].clone()).collect();
                    let sorted_inten: Vec<f32> =
                        indices.iter().map(|&i| intensity[i].clone()).collect();

                    // Squash the mobility dimension
                    let (mz, intensity) = squash_frame(&sorted_mz, &sorted_inten, 20.0);
                    let scan_start_time = frame.rt as f32 / 60.0;
                    let ion_injection_time = 100.0;
                    let total_ion_current = 0.0;
                    let id = frame.index.to_string();

                    let spec = RawSpectrum {
                        file_id,
                        precursors: vec![],
                        representation: Representation::Centroid,
                        scan_start_time,
                        ion_injection_time,
                        mz,
                        ms_level: 1,
                        id,
                        intensity,
                        total_ion_current,
                    };
                    Some(spec)
                }
                Err(x) => {
                    log::error!("error parsing spectrum: {:?}", x);
                    None
                }
            });

        // Merge the two
        let spectra: Vec<RawSpectrum> = ms1_spectra.chain(spectra.into_iter()).collect();

        Ok(spectra)
    }

    fn parse_precursor(dda_precursor: timsrust::Precursor) -> Precursor {
        let mut precursor: Precursor = Precursor::default();
        precursor.mz = dda_precursor.mz as f32;
        precursor.charge = match dda_precursor.charge {
            Some(x) => Some(x as u8),
            None => None,
        };
        precursor.intensity = match dda_precursor.intensity {
            Some(x) => Some(x as f32),
            None => None,
        };
        precursor.spectrum_ref = Option::from(dda_precursor.frame_index.to_string());
        precursor.inverse_ion_mobility = Option::from(dda_precursor.im as f32);
        precursor
    }
}
