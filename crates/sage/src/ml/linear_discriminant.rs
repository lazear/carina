//! Linear Discriminant Analysis for FDR refinement
//!
//! "What I cannot create, I do not understand" - Richard Feynman
//!
//! One of the major reasons for the creation of Sage is to develop a search
//! engine from first principles - And when I mean first principles, I mean
//! first principles - we are going to implement a basic linear algebra system
//! (complete with Gauss-Jordan elimination and eigenvector calculation) from scratch
//! to enable LDA.

use super::gauss::Gauss;
use super::matrix::Matrix;
use rayon::prelude::*;

use crate::mass::Tolerance;
use crate::scoring::Feature;

// Declare, so that we have compile time checking of matrix dimensions
const FEATURES: usize = 20;
const FEATURE_NAMES: [&str; FEATURES] = [
    "rank",
    "charge",
    "ln1p(hyperscore)",
    "ln1p(delta_next)",
    "ln1p(delta_best)",
    "delta_mass_model",
    "isotope_error",
    "average_ppm",
    "ln1p(-poisson)",
    "ln1p(matched_intensity_pct)",
    "ln1p(matched_peaks)",
    "ln1p(longest_b)",
    "ln1p(longest_y)",
    "longest_y_pct",
    "ln1p(peptide_len)",
    "missed_cleavages",
    "rt",
    "ims",
    "sqrt(delta_rt_model)",
    "sqrt(delta_ims_model)",
];

struct Features<'a>(&'a [f64]);

impl std::fmt::Debug for Features<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_map()
            .entries(FEATURE_NAMES.iter().zip(self.0))
            .finish()
    }
}

pub struct LinearDiscriminantAnalysis {
    eigenvector: Vec<f64>,
}

impl LinearDiscriminantAnalysis {
    pub fn train(features: &Matrix, decoy: &[bool]) -> Option<LinearDiscriminantAnalysis> {
        assert_eq!(features.rows, decoy.len());

        // Calculate class means, and overall mean
        let x_bar = features.mean();
        let mut scatter_within = Matrix::zeros(features.cols, features.cols);
        let mut scatter_between = Matrix::zeros(features.cols, features.cols);

        let mut class_means = Vec::new();

        for class in [true, false] {
            let count = decoy.iter().filter(|&label| *label == class).count();

            let class_data = (0..features.rows)
                .zip(decoy)
                .filter(|&(_, label)| *label == class)
                .flat_map(|(row, _)| features.row(row))
                .collect::<Vec<_>>();

            let mut class_data = Matrix::new(class_data, count, features.cols);
            let class_mean = class_data.mean();

            for row in 0..class_data.rows {
                for col in 0..class_data.cols {
                    class_data[(row, col)] -= class_mean[col];
                }
            }

            let cov = class_data.transpose().dot(&class_data) / class_data.rows as f64;
            scatter_within += cov;

            let diff = Matrix::col_vector(
                class_mean
                    .iter()
                    .zip(x_bar.iter())
                    .map(|(x, y)| x - y)
                    .collect::<Vec<_>>(),
            );

            scatter_between += diff.dot(&diff.transpose());
            class_means.extend(class_mean);
        }

        // Use overall mean as the initial vector for power method... seems
        // unlikely to be the actual best eigenvector!
        let mut evec =
            Gauss::solve(scatter_within, scatter_between).map(|mat| mat.power_method(&x_bar))?;

        // In some cases, power method can return eigenvector with signs flipped -
        // Make it so that Target class scores are higher than Decoy, so that
        // we can make assumptions about this for ranking
        let class_means = Matrix::new(class_means, 2, features.cols);
        let coef = class_means.dotv(&evec);
        if coef[1] < coef[0] {
            evec.iter_mut().for_each(|c| *c *= -1.0);
        }

        log::trace!("- linear model fit with {:?}", Features(&evec));

        Some(LinearDiscriminantAnalysis { eigenvector: evec })
    }

    pub fn score(&self, features: &Matrix) -> Vec<f64> {
        features.dotv(&self.eigenvector)
    }
}

pub fn score_psms(scores: &mut [Feature], precursor_tol: Tolerance) -> Option<()> {
    log::trace!("fitting linear discriminant model...");
    let decoys = scores
        .par_iter()
        .map(|sc| sc.label == -1)
        .collect::<Vec<_>>();

    let mass_error = match precursor_tol {
        Tolerance::Ppm(_, _) => |feat: &Feature| feat.delta_mass as f64,
        Tolerance::Pct(_, _) => unreachable!("Pct tolerance should never be used on mz"),
        Tolerance::Da(_, _) => |feat: &Feature| (feat.expmass - feat.calcmass) as f64,
    };

    let (bw_adjust, bin_size) = match precursor_tol {
        Tolerance::Ppm(lo, hi) => (2.0f64, (hi - lo).max(100.0)),
        Tolerance::Pct(_, _) => unreachable!("Pct tolerance should never be used on mz"),
        Tolerance::Da(lo, hi) => (0.1f64, (hi - lo).max(1000.0)),
    };

    let delta_mass = scores.par_iter().map(mass_error).collect::<Vec<_>>();

    let mass_model = super::kde::Builder::default()
        .monotonic(false)
        .bw_adjust(move |x| x * bw_adjust)
        .bins(bin_size.ceil().abs() as usize)
        .build(&delta_mass, &decoys);

    let features = scores
        .into_par_iter()
        .flat_map_iter(|perc| {
            let poisson = match (-perc.poisson).ln_1p() {
                x if x.is_finite() => x,
                _ => 3.5,
            };

            // Transform features - LDA requires that each feature is normally
            // distributed. This is not true for all of our inputs, so we log
            // transform many of them to get them closer to a gaussian distr.
            let x: [f64; FEATURES] = [
                (perc.rank as f64),
                (perc.charge as f64),
                (perc.hyperscore).ln_1p(),
                (perc.delta_next).ln_1p(),
                (perc.delta_best).ln_1p(),
                mass_model.posterior_error(mass_error(perc)),
                (perc.isotope_error as f64),
                (perc.average_ppm as f64),
                (poisson),
                (perc.matched_intensity_pct as f64).ln_1p(),
                (perc.matched_peaks as f64),
                (perc.longest_b as f64).ln_1p(),
                (perc.longest_y as f64).ln_1p(),
                (perc.longest_y as f64 / perc.peptide_len as f64),
                (perc.peptide_len as f64).ln_1p(),
                (perc.missed_cleavages as f64),
                (perc.aligned_rt as f64),
                (perc.ims as f64),
                (perc.delta_rt_model as f64).clamp(0.001, 0.999).sqrt(),
                (perc.delta_ims_model as f64).clamp(0.001, 0.999).sqrt(),
            ];
            x
        })
        .collect::<Vec<_>>();

    let features = Matrix::new(features, scores.len(), FEATURES);
    let lda = LinearDiscriminantAnalysis::train(&features, &decoys)?;
    if !lda.eigenvector.iter().all(|f| f.is_finite()) {
        log::error!(
            "linear model eigenvector includes NaN: this likely indicates a bug, please report!"
        );
        for row in 0..features.rows {
            if features.row(row).any(|f| !f.is_finite()) {
                let row = features.row(row).collect::<Vec<_>>();
                log::error!("example feature vector with NaN: {:?}", row);
                break;
            }
        }
        return None;
    }
    let discriminants = lda.score(&features);

    log::trace!("- fitting non-parametric model for posterior error probabilities");
    let kde = super::kde::Builder::default().build(&discriminants, &decoys);

    scores
        .par_iter_mut()
        .zip(&discriminants)
        .for_each(|(perc, score)| {
            perc.discriminant_score = *score as f32;
            perc.posterior_error = kde.posterior_error(*score).log10() as f32;
            if perc.posterior_error.is_infinite() {
                // This is approximately the log10 of the smallest positive
                // non-zero f64
                perc.posterior_error = -324.0;
            }
        });

    Some(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::ml::*;

    #[test]
    fn linear_discriminant() {
        let a = Matrix::new([1., 2., 3., 4.], 2, 2);
        let eigenvector = [0.4159736, 0.90937671];
        assert!(all_close(
            &a.power_method(&[0.54, 0.34]),
            &eigenvector,
            1E-5
        ));

        #[rustfmt::skip]
        let feats = Matrix::new(
            [
                5., 4., 3., 2., 
                4., 5., 4., 3., 
                6., 3., 4., 5., 
                1., 0., 2., 9., 
                5., 4., 4., 3., 
                2., 1., 1., 9.5, 
                1., 0., 2., 8., 
                3., 2., -2., 10.,
            ],
            8,
            4,
        );

        let lda = LinearDiscriminantAnalysis::train(
            &feats,
            &[false, false, false, true, false, true, true, true],
        )
        .expect("error training LDA");

        let mut scores = lda.score(&feats);
        let norm = norm(&scores);
        scores = scores.into_iter().map(|s| s / norm).collect();

        let expected = [
            0.49706043,
            0.48920177,
            0.48920177,
            -0.07209359,
            0.51204672,
            -0.02849527,
            -0.04924864,
            -0.06055943,
        ];

        assert!(
            all_close(&scores, &expected, 1E-8),
            "{:?} {:?}",
            scores,
            expected
        );
    }
}
