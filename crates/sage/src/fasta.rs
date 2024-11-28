use crate::enzyme::{Digest, EnzymeParameters};
use rayon::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
pub struct Fasta {
    pub targets: Vec<(Arc<String>, String)>,
    decoy_tag: String,
    // Should we ignore decoys in the fasta database
    // and generate them internally?
    generate_decoys: bool,
}

impl Fasta {
    // Parse a string into a fasta database
    pub fn parse<S: Into<String>>(contents: String, decoy_tag: S, generate_decoys: bool) -> Fasta {
        let decoy_tag = decoy_tag.into();

        let mut targets = Vec::new();
        let mut last_id = "";
        let mut s = String::new();

        for line in contents.as_str().lines() {
            if line.is_empty() {
                continue;
            }
            let line = line.trim();
            if let Some(id) = line.strip_prefix('>') {
                if !s.is_empty() {
                    let acc: Arc<String> =
                        Arc::new(last_id.split_ascii_whitespace().next().unwrap().to_string());
                    let seq = std::mem::take(&mut s);
                    if !acc.contains(&decoy_tag) || !generate_decoys {
                        targets.push((acc, seq));
                    }
                }
                last_id = id;
            } else {
                s.push_str(line);
            }
        }

        if !s.is_empty() {
            let acc: Arc<String> =
                Arc::new(last_id.split_ascii_whitespace().next().unwrap().to_string());
            if !acc.contains(&decoy_tag) || !generate_decoys {
                targets.push((acc, s));
            }
        }

        Fasta {
            targets,
            decoy_tag,
            generate_decoys,
        }
    }

    pub fn digest(&self, enzyme: &EnzymeParameters) -> Vec<Digest> {
        self.targets
            .par_iter()
            .flat_map_iter(|(protein, sequence)| {
                enzyme
                    .digest(sequence, protein.clone())
                    .into_iter()
                    .filter_map(|mut digest| {
                        if protein.contains(&self.decoy_tag) {
                            if !self.generate_decoys {
                                digest.decoy = true;
                                Some(digest)
                            } else {
                                None
                            }
                        } else {
                            Some(digest)
                        }
                    })
            })
            .collect()
    }

    pub fn iter_chunks(&self, chunk_size: usize) -> impl Iterator<Item = Self> + '_ {
        self.targets
            .chunks(chunk_size)
            .map(move |target_chunk| Self {
                targets: target_chunk.to_vec(),
                decoy_tag: self.decoy_tag.clone(),
                generate_decoys: self.generate_decoys,
            })
    }
}
