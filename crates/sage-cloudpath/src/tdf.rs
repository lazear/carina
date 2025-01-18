use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use std::cmp::Ordering;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::SpectrumReader;
pub use timsrust::readers::SpectrumReaderConfig as BrukerSpectrumProcessor;

pub struct TdfReader;

impl TdfReader {
    pub fn parse(
        &self,
        path_name: impl AsRef<str>,
        file_id: usize,
        bruker_spectrum_processor: BrukerSpectrumProcessor,
        requires_ms1: bool,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let spectrum_reader = timsrust::readers::SpectrumReader::build()
            .with_path(path_name.as_ref())
            .with_config(bruker_spectrum_processor)
            .finalize()?;
        let mut spectra = self.read_msn_spectra(file_id, &spectrum_reader)?;
        if requires_ms1 {
            let ms1s = self.read_ms1_spectra(&path_name, file_id, &spectrum_reader)?;
            spectra.extend(ms1s);
        }

        Ok(spectra)
    }

    fn read_ms1_spectra(
        &self,
        path_name: impl AsRef<str>,
        file_id: usize,
        spectrum_reader: &SpectrumReader,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
        let start = std::time::Instant::now();
        let frame_reader = timsrust::readers::FrameReader::new(path_name.as_ref())?;
        let tdf_path = std::path::Path::new(path_name.as_ref()).join("analysis.tdf");
        let metadata = timsrust::readers::MetadataReader::new(tdf_path)?;
        let mz_converter = metadata.mz_converter;
        let ims_converter = metadata.im_converter;

        let ms1_spectra: Vec<RawSpectrum> = frame_reader
            .parallel_filter(|f| f.ms_level == timsrust::MSLevel::MS1)
            .filter_map(|frame| match frame {
                Ok(frame) => {
                    let mz: Vec<f32> = frame
                        .tof_indices
                        .iter()
                        .map(|&x| mz_converter.convert(x as f64) as f32)
                        .collect();
                    let intensity: Vec<f32> = frame.intensities.iter().map(|&x| x as f32).collect();
                    let mut imss: Vec<f32> = vec![0.0; mz.len()];
                    // TODO: This is getting pretty big ... I should refactor this block.
                    frame
                        .scan_offsets
                        .windows(2)
                        .enumerate()
                        .map(|(i, w)| {
                            let num = w[1] - w[0];
                            if num == 0 {
                                return None;
                            }
                            let lo = w[0];
                            let hi = w[1];

                            let im = ims_converter.convert(i as f64) as f32;
                            Some((im, lo, hi))
                        })
                        .for_each(|x| {
                            if let Some((im, lo, hi)) = x {
                                for i in lo..hi {
                                    imss[i] = im;
                                }
                            }
                        });

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
                    let sorted_imss: Vec<f32> = indices.iter().map(|&i| imss[i].clone()).collect();

                    // Squash the mobility dimension
                    let tol_ppm = 15.0;
                    let im_tol_pct = 2.0;
                    let (mz, (intensity, mobility)): (Vec<f32>, (Vec<f32>, Vec<f32>)) =
                        dumbcentroid_frame(
                            &sorted_mz,
                            &sorted_inten,
                            &sorted_imss,
                            tol_ppm,
                            im_tol_pct,
                        );

                    let scan_start_time = frame.rt as f32 / 60.0;
                    let ion_injection_time = 100.0; // This is made up, in theory we can read
                                                    // if from the tdf file
                    let total_ion_current = sorted_inten.iter().sum::<f32>();
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
                        mobility: Some(mobility),
                    };
                    Some(spec)
                }
                Err(x) => {
                    log::error!("error parsing spectrum: {:?}", x);
                    None
                }
            })
            .collect();
        log::info!(
            "read {} ms1 spectra in {:#?}",
            ms1_spectra.len(),
            start.elapsed()
        );
        Ok(ms1_spectra)
    }

    fn read_msn_spectra(
        &self,
        file_id: usize,
        spectrum_reader: &SpectrumReader,
    ) -> Result<Vec<RawSpectrum>, timsrust::TimsRustError> {
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
                            ion_injection_time: dda_precursor.rt as f32,
                            total_ion_current: 0.0,
                            mz: dda_spectrum.mz_values.iter().map(|&x| x as f32).collect(),
                            ms_level: 2,
                            id: dda_spectrum.index.to_string(),
                            intensity: dda_spectrum.intensities.iter().map(|&x| x as f32).collect(),
                            mobility: None,
                        };
                        Some(spectrum)
                    }
                    None => None,
                },
                Err(_) => None,
            })
            .collect();
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

fn dumbcentroid_frame(
    mz_array: &[f32],
    intensity_array: &[f32],
    ims_array: &[f32],
    mz_tol_ppm: f32,
    im_tol_pct: f32,
) -> (Vec<f32>, (Vec<f32>, Vec<f32>)) {
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

    #[derive(Debug, Clone, Copy, Default)]
    struct ImsPeak {
        mz: f32,
        intensity: f32,
        im: f32,
    }
    let mut agg_buff = Vec::with_capacity(10_000.min(arr_len));
    const TOUCH_BUFF_SIZE: usize = 1000;
    let mut touch_buff = [false; TOUCH_BUFF_SIZE];

    let utol = mz_tol_ppm / 1e6;
    let im_tol = im_tol_pct / 100.0;

    for &idx in &order {
        if touched[idx] {
            continue;
        }

        let mz = mz_array[idx];
        let im = ims_array[idx];
        let da_tol = mz * utol;
        let left_e = mz - da_tol;
        let right_e = mz + da_tol;

        let mut ss_start = mz_array.partition_point(|&x| x < left_e);
        let mut ss_end = mz_array.partition_point(|&x| x <= right_e);

        let mut slice_width = ss_end - ss_start;
        if slice_width >= 1000 {
            // It is EXCEEDINGLY UNLIKELY that more than 1000 points
            // will be aggregated or in the range of a single mz peak.
            // Here we just handle those edge cases by making sure the 'center'
            // will be aggregated.

            let new_ss_start = idx.saturating_sub(TOUCH_BUFF_SIZE / 2);
            let new_ss_end = new_ss_start + TOUCH_BUFF_SIZE - 1;
            // TODO: make a better warning message here.
            log::warn!(
                "More than {} points are in the mz range of a point, limiting span",
                TOUCH_BUFF_SIZE
            );
            ss_start = new_ss_start;
            ss_end = new_ss_end;
            slice_width = ss_end - ss_start;
        }
        let local_num_touched = touched[ss_start..ss_end].iter().filter(|&&x| x).count();
        let local_num_untouched = slice_width - local_num_touched;

        if local_num_touched == slice_width {
            continue;
        }

        let abs_im_tol = im * im_tol;
        let left_im = im - abs_im_tol;
        let right_im = im + abs_im_tol;

        let mut curr_intensity = 0.0;

        let mut num_touchable = 0;
        for i in ss_start..ss_end {
            let im_i = ims_array[i];
            if !touched[i] && intensity_array[i] > 0.0 && im_i >= left_im && im_i <= right_im {
                curr_intensity += intensity_array[i];
                num_touchable += 1;
                touch_buff[i - ss_start] = true;
            }
        }

        if curr_intensity > 0.0 {
            agg_buff.push(ImsPeak {
                mz,
                intensity: curr_intensity,
                im,
            });
            touched[ss_start..ss_end]
                .iter_mut()
                .zip(touch_buff.iter_mut().take(slice_width))
                .for_each(|(t, tb)| {
                    *t = true;
                    *tb = false;
                });
            global_num_touched += num_touchable;
            const MAX_PEAKS: usize = 10000;
            if agg_buff.len() > MAX_PEAKS {
                let curr_loc_int = intensity_array[idx];
                if curr_loc_int > 200.0 {
                    log::debug!(
                        "Reached limit of the agg buffer at index {}/{} curr int={}",
                        idx,
                        arr_len,
                        curr_loc_int
                    );
                }
                break;
            }
        }

        if global_num_touched == arr_len {
            break;
        }
    }

    assert!(touch_buff.iter().all(|x| !x), "{:?}", touch_buff);

    // Drop the zeros and sort
    // I could in theory truncate instead of filtering.
    let mut result: Vec<(f32, (f32, f32))> = agg_buff
        .iter()
        .filter(|&x| x.mz > 0.0 && x.intensity > 0.0)
        .map(|x| (x.mz, (x.intensity, x.im as f32)))
        .collect();

    result.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    // println!("Centroiding: Start len: {}; end len: {};", arr_len, result.len());
    // Ultra data is usually start: 40k end 10k,
    // HT2 data is usually start 400k end 40k

    result.into_iter().unzip()
}
