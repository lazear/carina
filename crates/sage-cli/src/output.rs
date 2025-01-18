use rayon::prelude::*;
use sage_core::spectrum::MS1Spectra;
use sage_core::{scoring::Feature, tmt::TmtQuant};

#[derive(Default)]
pub struct SageResults {
    pub ms1: MS1Spectra,
    pub features: Vec<Feature>,
    pub quant: Vec<TmtQuant>,
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
                match (acc.ms1, x.ms1) {
                    (MS1Spectra::NoMobility(mut a), MS1Spectra::NoMobility(b)) => {
                        a.extend(b);
                        acc.ms1 = MS1Spectra::NoMobility(a);
                    }
                    (MS1Spectra::WithMobility(mut a), MS1Spectra::WithMobility(b)) => {
                        a.extend(b);
                        acc.ms1 = MS1Spectra::WithMobility(a);
                    }
                    (MS1Spectra::Empty, MS1Spectra::Empty) => {
                        acc.ms1 = MS1Spectra::Empty;
                    }
                    (MS1Spectra::Empty, MS1Spectra::WithMobility(a))
                    | (MS1Spectra::WithMobility(a), MS1Spectra::Empty) => {
                        acc.ms1 = MS1Spectra::WithMobility(a);
                    }
                    (MS1Spectra::Empty, MS1Spectra::NoMobility(a))
                    | (MS1Spectra::NoMobility(a), MS1Spectra::Empty) => {
                        acc.ms1 = MS1Spectra::NoMobility(a);
                    }
                    _ => {
                        unreachable!()
                    }
                };
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
                match (acc.ms1, x.ms1) {
                    (MS1Spectra::NoMobility(mut a), MS1Spectra::NoMobility(b)) => {
                        a.extend(b);
                        acc.ms1 = MS1Spectra::NoMobility(a);
                    }
                    (MS1Spectra::WithMobility(mut a), MS1Spectra::WithMobility(b)) => {
                        a.extend(b);
                        acc.ms1 = MS1Spectra::WithMobility(a);
                    }
                    (MS1Spectra::Empty, MS1Spectra::Empty) => {
                        acc.ms1 = MS1Spectra::Empty;
                    }
                    (MS1Spectra::Empty, MS1Spectra::WithMobility(a))
                    | (MS1Spectra::WithMobility(a), MS1Spectra::Empty) => {
                        acc.ms1 = MS1Spectra::WithMobility(a);
                    }
                    (MS1Spectra::Empty, MS1Spectra::NoMobility(a))
                    | (MS1Spectra::NoMobility(a), MS1Spectra::Empty) => {
                        acc.ms1 = MS1Spectra::NoMobility(a);
                    }
                    _ => {
                        unreachable!()
                    }
                };
                acc
            })
    }
}
