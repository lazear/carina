//! False discovery rate control using double-competition (picked-peptide &
//! picked-protein) approaches
//!
//! Lin et al., https://pubmed.ncbi.nlm.nih.gov/36166314/
//! Savitski et al., https://pubmed.ncbi.nlm.nih.gov/25987413/

use crate::database::{IndexedDatabase, PeptideIx};
use crate::lfq::PrecursorId;
use crate::ml::kde::Estimator;
use crate::scoring::Feature;
use fnv::FnvHashMap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::BuildHasher;

#[derive(Copy, Clone, Debug)]
pub struct Competition<Ix> {
    pub forward: f32,
    pub foward_ix: Option<Ix>,
    pub reverse: f32,
    pub reverse_ix: Option<Ix>,
}

struct Row<Ix> {
    ix: Ix,
    decoy: bool,
    score: f32,
    q: f32,
}

impl<Ix: Default + Send> Default for Competition<Ix> {
    fn default() -> Self {
        Self {
            forward: f32::MIN,
            reverse: f32::MIN,
            foward_ix: None,
            reverse_ix: None,
        }
    }
}

impl<Ix: Default + Send> Competition<Ix> {
    fn score(&self) -> f32 {
        self.forward.max(self.reverse)
    }

    fn is_decoy(&self) -> bool {
        self.reverse >= self.forward
    }

    fn fit_kde<K, B>(scores: &HashMap<K, Self, B>) -> Estimator {
        let (scores, decoys): (Vec<f64>, Vec<bool>) = scores
            .values()
            .map(|score| (score.score() as f64, score.is_decoy()))
            .unzip();
        crate::ml::kde::Builder::default().build(&scores, &decoys)
    }

    fn assign_q_value<K, B>(
        scores: HashMap<K, Self, B>,
        threshold: f32,
    ) -> (HashMap<Ix, f32, B>, usize)
    where
        K: Eq + std::hash::Hash + Send,
        Ix: Eq + std::hash::Hash,
        B: BuildHasher + Default + Send,
    {
        let estimator = Self::fit_kde(&scores);
        let mut scores = scores
            .into_par_iter()
            .flat_map(|(_, comp)| {
                [
                    (comp.foward_ix, false, comp.forward),
                    (comp.reverse_ix, true, comp.reverse),
                ]
            })
            .filter_map(|(ix, decoy, score)| {
                ix.map(|ix| Row {
                    ix,
                    decoy,
                    score,
                    q: 1.0,
                })
            })
            .collect::<Vec<Row<Ix>>>();

        scores.par_sort_by(|a, b| b.score.total_cmp(&a.score));

        let mut decoy = 1.0;
        let mut target = 0.0;
        for score in scores.iter_mut() {
            let pep = estimator.posterior_error(score.score as f64) as f32;

            // Cumulative sum of PEP ~ # of decoys
            decoy += pep;
            if !score.decoy {
                target += 1.0;
            }
            score.q = decoy / target;
        }
        // Q-value is the minimum q-value at any given score threshold
        // `q = q[::-1].cummin()[::-1] in python`
        let mut q_min = 1.0f32;
        let mut passing = 0;
        for score in scores.iter_mut().rev() {
            q_min = q_min.min(score.q);
            score.q = q_min;
            if q_min <= threshold && !score.decoy {
                passing += 1;
            }
        }

        (
            scores
                .into_par_iter()
                .map(|score| (score.ix, score.q))
                .collect(),
            passing,
        )
    }
}

pub fn picked_peptide(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    let mut map: FnvHashMap<String, Competition<PeptideIx>> = FnvHashMap::default();
    for feat in features.iter() {
        let peptide = &db[feat.peptide_idx];
        // Only reverse the peptide sequence if we generated decoys ourselves
        let key = match db.generate_decoys && peptide.decoy {
            true => peptide.reverse().to_string(),
            false => peptide.to_string(),
        };

        let entry = map.entry(key).or_default();
        match peptide.decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some(feat.peptide_idx);
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some(feat.peptide_idx);
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features.par_iter_mut().for_each(|feat| {
        feat.peptide_q = scores[&feat.peptide_idx];
    });

    passing
}

pub fn picked_protein(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    let mut map: FnvHashMap<_, Competition<String>> = FnvHashMap::default();
    for feat in features.iter() {
        let decoy = db[feat.peptide_idx].decoy;
        let entry = map.entry(&db[feat.peptide_idx].proteins).or_default();
        let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        match decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some(proteins);
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some(proteins);
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features.par_iter_mut().for_each(|feat| {
        let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        feat.protein_q = scores[&proteins];
    });

    passing
}

pub fn picked_precursor(
    peaks: &mut FnvHashMap<(PrecursorId, bool), (crate::lfq::Peak, Vec<f64>)>,
) -> usize {
    // let mut map: FnvHashMap<PeptideIx, Competition<(PeptideIx, bool)>> = FnvHashMap::default();
    // for (key, (peak, _)) in peaks.iter() {
    //     let entry = map.entry(key.0).or_default();
    //     match key.1 {
    //         true => {
    //             entry.reverse = entry.reverse.max(peak.score as f32);
    //             entry.reverse_ix = Some(*key);
    //         }
    //         false => {
    //             entry.forward = entry.forward.max(peak.score as f32);
    //             entry.foward_ix = Some(*key);
    //         }
    //     }
    // }
    let mut scores = peaks
        .par_iter()
        .map(|(&(ix, decoy), (peak, _))| Row {
            ix,
            decoy,
            score: peak.score as f32,
            q: 1.0,
        })
        .collect::<Vec<_>>();

    scores.par_sort_by(|a, b| b.score.total_cmp(&a.score));

    let mut decoy = 1.0;
    let mut target = 0.0;
    for score in scores.iter_mut() {
        match score.decoy {
            true => decoy += 1.0,
            false => target += 1.0,
        };
        score.q = decoy / target;
    }
    // Q-value is the minimum q-value at any given score threshold
    // `q = q[::-1].cummin()[::-1] in python`
    let mut q_min = 1.0f32;
    let mut passing = 0;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.q);
        score.q = q_min;
        if q_min <= 0.05 && !score.decoy {
            passing += 1;
        }
    }

    let scores = scores
        .into_par_iter()
        .map(|score| ((score.ix, score.decoy), score.q))
        .collect::<FnvHashMap<_, _>>();

    peaks.par_iter_mut().for_each(|(ix, (peak, _))| {
        peak.q_value = scores[ix];
    });
    passing
}
