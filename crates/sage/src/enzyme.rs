use fnv::FnvHashSet;
use regex::Regex;
use std::sync::Arc;

use crate::mass::VALID_AA;

#[derive(Clone, PartialOrd, Ord, Debug, Default)]
/// An enzymatic digest
///
/// # Important invariant about [`Digest`]:
/// * two digests are equal if and only if their sequences and position are equal
///   i.e., decoy status is ignored for equality and hashing
pub struct Digest {
    /// Is this a decoy peptide?
    pub decoy: bool,
    /// Semi-enzymatic?
    pub semi_enzymatic: bool,
    /// Cleaved peptide sequence
    pub sequence: String,
    /// Protein accession
    pub protein: Arc<str>,
    /// Missed cleavages
    pub missed_cleavages: u8,
    /// Is this an N-terminal peptide of the protein?
    pub position: Position,
}

pub struct DigestGroup {
    pub reference: Digest,
    pub proteins: Vec<Arc<str>>,
}

pub fn group_digests(mut digests: Vec<Digest>) -> Vec<DigestGroup> {
    let mut groups = Vec::new();
    digests.sort_unstable_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then(a.decoy.cmp(&b.decoy))
            .then(a.sequence.cmp(&b.sequence))
    });
    let mut curr_group = DigestGroup {
        reference: digests[0].clone(),
        proteins: Vec::new(),
    };
    for digest in digests {
        if digest.decoy == curr_group.reference.decoy
            && digest.position == curr_group.reference.position
            && digest.sequence == curr_group.reference.sequence
        {
            curr_group.proteins.push(digest.protein);
        } else {
            curr_group.proteins.sort_unstable();
            groups.push(curr_group);
            curr_group = DigestGroup {
                reference: digest.clone(),
                proteins: vec![digest.protein],
            };
        }
    }
    groups.push(curr_group);
    groups
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum Position {
    Nterm,
    Cterm,
    Full,
    Internal,
}

impl Default for Position {
    fn default() -> Self {
        Self::Internal
    }
}

impl Digest {
    /// Generate an internal decoy sequence by reversing the sequence
    /// inside the first and last amino acids
    pub fn reverse(&self) -> Self {
        if self.decoy {
            return self.clone();
        }

        let mut sequence = self.sequence.chars().rev().collect::<Vec<char>>();
        let n = sequence.len().saturating_sub(1);
        sequence.swap(0, n);

        Digest {
            decoy: true,
            semi_enzymatic: self.semi_enzymatic,
            protein: self.protein.clone(),
            sequence: sequence.into_iter().collect(),
            missed_cleavages: self.missed_cleavages,
            position: self.position,
        }
    }
}

impl PartialEq for Digest {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.position == other.position
    }
}

impl Eq for Digest {
    fn assert_receiver_is_total_eq(&self) {}
}

impl std::hash::Hash for Digest {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
        self.position.hash(state);
    }
}

pub struct EnzymeParameters {
    /// Number of missed cleavages to produce
    pub missed_cleavages: u8,
    /// Inclusive
    pub min_len: usize,
    /// Inclusive
    pub max_len: usize,
    pub enyzme: Option<Enzyme>,
}

#[derive(Clone)]
pub struct Enzyme {
    // Skip cleaving if the site is followed matching this AA
    pub skip_suffix: Option<char>,
    // Regex for matching cleavage sites
    regex: Regex,
    // Cleave at c-terminal?
    pub c_terminal: bool,
    // Semi-enzymatic cleavage?
    pub semi_enzymatic: bool,
}

#[derive(Clone)]
pub struct DigestSite {
    // Range defining cleavage position
    pub site: std::ops::Range<usize>,
    // Number of missed cleavages
    pub missed_cleavages: u8,

    pub semi_enzymatic: bool,
}

impl Enzyme {
    pub fn new(
        cleave: &str,
        skip_suffix: Option<char>,
        c_terminal: bool,
        semi_enzymatic: bool,
    ) -> Option<Self> {
        assert!(
            cleave.chars().all(|x| VALID_AA.contains(&(x as u8))) || cleave == "$",
            "Enzyme cleavage sequence contains non-amino acid characters: {}",
            cleave
        );
        assert!(
            skip_suffix
                .map(|x| VALID_AA.contains(&(x as u8)))
                .unwrap_or(true),
            "Enzyme cleavage restriction is non-amino acid character: {}",
            skip_suffix.unwrap(),
        );

        // At this point, cleave can be three things: empty, "$", or a string of valid AA's
        match cleave {
            "" => None,
            "$" => Some(Enzyme {
                regex: Regex::new("$").unwrap(),
                skip_suffix: None,
                // Allowing this to be set to false could cause unexpected behavior
                c_terminal: true,
                // Do not allow strange behavior
                semi_enzymatic: false,
            }),
            _ => Some(Enzyme {
                regex: Regex::new(&format!("[{}]", cleave.replace('?', ""))).unwrap(),
                skip_suffix,
                c_terminal,
                semi_enzymatic,
            }),
        }
    }

    pub fn cleavage_sites(&self, sequence: &str) -> Vec<DigestSite> {
        let mut sites = Vec::new();
        let mut left = 0;
        for mat in self.regex.find_iter(sequence) {
            let right = match self.c_terminal {
                true => mat.end(),
                false => mat.start(),
            };
            if let Some(skip) = self.skip_suffix {
                if right < sequence.len() && sequence[right..].starts_with(skip) {
                    continue;
                }
            }
            sites.push(DigestSite {
                site: left..right,
                missed_cleavages: 0,
                semi_enzymatic: false,
            });
            left = right;
        }
        sites.push(DigestSite {
            site: left..sequence.len(),
            missed_cleavages: 0,
            semi_enzymatic: false,
        });
        sites
    }
}

impl EnzymeParameters {
    pub fn cleavage_sites(&self, sequence: &str) -> Vec<DigestSite> {
        match &self.enyzme {
            Some(enzyme) => enzyme.cleavage_sites(sequence),
            None => {
                // Perform a non-specific digest
                let mut v = Vec::new();
                for len in self.min_len..=self.max_len {
                    for i in 0..=sequence.len().saturating_sub(len) {
                        v.push(DigestSite {
                            site: i..i + len,
                            missed_cleavages: 0,
                            semi_enzymatic: false,
                        });
                    }
                }
                v
            }
        }
    }

    fn missed_cleavage_sites(
        &self,
        sites: &mut Vec<DigestSite>,
        missed_cleavages: u8,
    ) -> Vec<DigestSite> {
        let mut missed_cleavage_sites = Vec::new();
        for cleavage in 1..=(1 + missed_cleavages) {
            // Generate missed cleavages
            for win in sites.windows(cleavage as usize) {
                let start = win[0].site.start;
                let end = win[cleavage as usize - 1].site.end;
                missed_cleavage_sites.push(DigestSite {
                    site: start..end,
                    missed_cleavages: cleavage - 1,
                    semi_enzymatic: false,
                });
            }
        }
        sites.append(&mut missed_cleavage_sites);
        sites.to_vec()
    }

    fn is_semi_enzymatic(&self) -> bool {
        match &self.enyzme {
            Some(enzyme) => enzyme.semi_enzymatic,
            None => false,
        }
    }

    fn semi_enzymatic_sites(&self, sites: &mut Vec<DigestSite>) -> Vec<DigestSite> {
        let mut semi_enzymatic_sites = Vec::new();
        for site in sites.iter_mut() {
            let start = site.site.start;
            let end = site.site.end;
            for cut_site in start..end {
                // Missed cleavages should actually be split across the sites,
                // but we rely on the fact that we generate missed cleavages in ascending
                // order, and eliminate previously seen digests
                semi_enzymatic_sites.push(DigestSite {
                    site: start..cut_site,
                    missed_cleavages: site.missed_cleavages,
                    semi_enzymatic: true,
                });
                semi_enzymatic_sites.push(DigestSite {
                    site: cut_site..end,
                    missed_cleavages: site.missed_cleavages,
                    semi_enzymatic: true,
                });
            }
        }
        sites.append(&mut semi_enzymatic_sites);
        sites.to_vec()
    }

    pub fn digest(&self, sequence: &str, protein: Arc<str>) -> Vec<Digest> {
        let n = sequence.len();
        let mut digests = Vec::new();
        let mut sites = self.cleavage_sites(sequence);
        // Allowing missed_cleavages with non-specific digest causes OOB panics
        // in the below indexing code
        let missed_cleavages = match self.enyzme {
            None => 0,
            _ => self.missed_cleavages,
        };

        let mut sites = match missed_cleavages {
            0 => sites,
            _ => self.missed_cleavage_sites(&mut sites, missed_cleavages),
        };

        let mut sites = match self.is_semi_enzymatic() {
            false => sites,
            true => self.semi_enzymatic_sites(&mut sites),
        };

        // Keep a set of peptides that have been digested from this sequence
        // - handles cases where the same peptide occurs multiple times in a protein
        let mut seen = FnvHashSet::default();

        for site in sites.iter_mut() {
            let start = site.site.start;
            let end = site.site.end;

            let sequence = match sequence.get(start..end) {
                Some(sequence) => sequence,
                None => continue,
            };

            let len = sequence.len();

            let position = match (start == 0, end == n) {
                (true, true) => Position::Full,
                (true, false) => Position::Nterm,
                (false, true) => Position::Cterm,
                (false, false) => Position::Internal,
            };

            if len >= self.min_len && len <= self.max_len && len > 0 && seen.insert(sequence) {
                digests.push(Digest {
                    sequence: sequence.into(),
                    missed_cleavages: site.missed_cleavages,
                    decoy: false,
                    semi_enzymatic: site.semi_enzymatic,
                    position,
                    protein: protein.clone(),
                });
            }
        }
        digests
    }
}

#[cfg(test)]
mod test {
    use quickcheck_macros::quickcheck;
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn hash_digest() {
        let mut digests = vec![
            Digest {
                decoy: false,
                semi_enzymatic: false,
                sequence: "MADEEK".into(),
                missed_cleavages: 0,
                position: Position::Nterm,
                protein: Arc::from(String::default()),
            },
            Digest {
                decoy: false,
                semi_enzymatic: false,
                sequence: "MADEEK".into(),
                missed_cleavages: 0,
                position: Position::Nterm,
                protein: Arc::from(String::default()),
            },
        ];

        // Make sure hashing a digest works
        let set = digests.drain(..).collect::<HashSet<_>>();
        assert_eq!(set.len(), 1);

        let mut digests = vec![
            Digest {
                decoy: false,
                semi_enzymatic: false,
                sequence: "MADEEK".into(),
                missed_cleavages: 0,
                position: Position::Nterm,
                protein: Arc::from(String::default()),
            },
            Digest {
                decoy: false,
                semi_enzymatic: false,
                sequence: "MADEEK".into(),
                missed_cleavages: 0,
                position: Position::Internal,
                protein: Arc::from(String::default()),
            },
        ];

        // // Make sure hashing a digest works
        let set = digests.drain(..).collect::<HashSet<_>>();
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn trypsin() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            ("MADEEK".into(), Position::Nterm),
            ("LPPGWEK".into(), Position::Internal),
            ("MSR".into(), Position::Internal),
            ("SSGR".into(), Position::Internal),
            ("VYYFNHITNASQWERPSGN".into(), Position::Cterm),
        ];

        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", Some('P'), true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| (d.sequence, d.position))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn trypsin_missed_cleavage() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "R",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWERPSGN",
            "MADEEKLPPGWEK",
            "LPPGWEKR",
            "RMSR",
            "MSRSSGR",
            "SSGRVYYFNHITNASQWERPSGN",
        ];

        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 1,
            enyzme: Enzyme::new("KR", Some('P'), true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn trypsin_missed_cleavage_2() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "R",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWERPSGN",
            "MADEEKLPPGWEK",
            "LPPGWEKR",
            "RMSR",
            "MSRSSGR",
            "SSGRVYYFNHITNASQWERPSGN",
            "MADEEKLPPGWEKR",
            "LPPGWEKRMSR",
            "RMSRSSGR",
            "MSRSSGRVYYFNHITNASQWERPSGN",
        ];

        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: 50,
            missed_cleavages: 2,
            enyzme: Enzyme::new("KR", Some('P'), true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_trypsin_pro() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",
            "MSR",
            "SSGR",
            "VYYFNHITNASQWER",
            "PSGN",
        ];

        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", None, true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_asp_n() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
        let expected = vec!["MA", "DEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW"];

        let tryp = EnzymeParameters {
            min_len: 1,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("D", None, false, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_chymotrypsin_pro() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
        let expected = vec![
            "MADEEKL",
            "PPGW",
            "EKRMSRSSGRVY",
            "Y",
            "F",
            "NHITNASQW",
            "ERPSGNW",
        ];

        let tryp = EnzymeParameters {
            min_len: 1,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("FYWL", None, true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nonspecific_digest_5() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";

        let expected = sequence
            .as_bytes()
            .windows(5)
            .flat_map(std::str::from_utf8)
            .collect::<Vec<_>>();

        let tryp = EnzymeParameters {
            min_len: 5,
            max_len: 5,
            missed_cleavages: 0,
            enyzme: None,
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn nonspecific_digest_5_7() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";

        let expected = (5..=7)
            .flat_map(|window| {
                sequence
                    .as_bytes()
                    .windows(window)
                    .flat_map(std::str::from_utf8)
            })
            .collect::<Vec<_>>();

        let tryp = EnzymeParameters {
            min_len: 5,
            max_len: 7,
            missed_cleavages: 0,
            enyzme: Enzyme::new("", None, true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn no_digest() {
        let sequence = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNW";
        let expected = vec![sequence];

        let tryp = EnzymeParameters {
            min_len: 0,
            max_len: usize::MAX,
            missed_cleavages: 0,
            enyzme: Enzyme::new("$", None, true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn ensure_unique() {
        let sequence = "KVEGAQNQGKKVEGAQNQGK";
        let expected = vec!["VEGAQNQGK"];

        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: usize::MAX,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", None, true, false),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn mini_semi_trypsin() {
        let sequence = "MADEEK";
        let expected = vec![
            "MADEEK", "ADEEK", "MA", "DEEK", "MAD", "EEK", "MADE", "EK", "MADEE",
        ];

        let tryp = EnzymeParameters {
            min_len: 2,
            max_len: 50,
            missed_cleavages: 0,
            enyzme: Enzyme::new("KR", None, true, true),
        };

        assert_eq!(
            expected,
            tryp.digest(sequence, Arc::default())
                .into_iter()
                .map(|d| d.sequence)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn semi_trypsin_trypsin_missed_cleavage() {
        let sequence = "MADEEKLPPGWEK";
        let expected = vec![
            "MADEEK",
            "LPPGWEK",       // normal KR
            "MADEEKLPPGWEK", // one missed cleavage
            "ADEEK",
            "DEEK",
            "MAD",
            "EEK",
            "MADE",
            "MADEE", // normal half-tryptics
            "PPGWEK",
            "PGWEK",
            "LPP",
            "GWEK",
            "LPPG",
            "WEK",
            "LPPGW",
            "LPPGWE",
            "ADEEKLPPGWEK",
            "DEEKLPPGWEK",
            "EEKLPPGWEK", // one missed cleavage half-tryptics
            "EKLPPGWEK",
            "KLPPGWEK",
            "MADEEKL",
            "MADEEKLP",
            "MADEEKLPP",
            "MADEEKLPPG",
            "MADEEKLPPGW",
            "MADEEKLPPGWE",
        ];

        let tryp = EnzymeParameters {
            min_len: 3,
            max_len: 50,
            missed_cleavages: 1,
            enyzme: Enzyme::new("KR", None, true, true),
        };

        for (digest, expected) in tryp
            .digest(sequence, Arc::default())
            .into_iter()
            .zip(expected)
        {
            assert_eq!(digest.sequence, expected);
            // reverse and skip the first (C-terminal) AA, counting interior missed cleavages
            let missed_cleavages = digest
                .sequence
                .as_bytes()
                .iter()
                .rev()
                .skip(1)
                .map(|s| (*s == b'K' || *s == b'R') as u8)
                .sum::<u8>();
            assert_eq!(
                missed_cleavages, digest.missed_cleavages,
                "{}",
                digest.sequence
            );

            if digest.sequence.starts_with("MAD") && digest.sequence != sequence {
                assert_eq!(digest.position, Position::Nterm);
            }
        }
    }

    /// Helper struct for generation of random sequences of valid amino acids
    #[derive(Clone, Debug)]
    struct RandomSequence {
        sequence: String,
    }

    impl quickcheck::Arbitrary for RandomSequence {
        fn arbitrary(g: &mut quickcheck::Gen) -> Self {
            let bytes = (0..g.size())
                .filter_map(|_| g.choose(&VALID_AA))
                .copied()
                .collect();
            Self {
                sequence: String::from_utf8(bytes).unwrap(),
            }
        }
    }

    #[quickcheck]
    /// Check that our strict ordering of missed cleavage generation is not
    /// broken for arbitrary peptide sequences
    fn quickcheck_semi_missed_cleavages(RandomSequence { sequence }: RandomSequence) {
        let tryp = EnzymeParameters {
            min_len: 3,
            max_len: 50,
            missed_cleavages: 2,
            enyzme: Enzyme::new("KR", None, true, true),
        };

        for digest in tryp.digest(&sequence, Arc::default()) {
            // reverse and skip the first (C-terminal) AA, counting interior missed cleavages
            let missed_cleavages = digest
                .sequence
                .as_bytes()
                .iter()
                .rev()
                .skip(1)
                .map(|s| (*s == b'K' || *s == b'R') as u8)
                .sum::<u8>();
            assert_eq!(
                missed_cleavages, digest.missed_cleavages,
                "{}",
                digest.sequence
            );

            assert!(digest.missed_cleavages <= 2);
        }
    }
}
