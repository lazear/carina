use std::{
    collections::HashMap,
    fmt::{Display, Write},
    str::FromStr,
};

use serde::Serialize;

use crate::mass::VALID_AA;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModificationSpecificity {
    PeptideN(Option<u8>),
    PeptideC(Option<u8>),
    ProteinN(Option<u8>),
    ProteinC(Option<u8>),
    Residue(u8),
}

impl Display for ModificationSpecificity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = match self {
            ModificationSpecificity::PeptideN(r) => {
                f.write_char('^')?;
                *r
            }
            ModificationSpecificity::PeptideC(r) => {
                f.write_char('$')?;
                *r
            }
            ModificationSpecificity::ProteinN(r) => {
                f.write_char('[')?;
                *r
            }
            ModificationSpecificity::ProteinC(r) => {
                f.write_char(']')?;
                *r
            }
            ModificationSpecificity::Residue(r) => Some(*r),
        };

        if let Some(r) = r {
            f.write_char(r as char)?;
        }

        Ok(())
    }
}

impl Serialize for ModificationSpecificity {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InvalidModification {
    Empty,
    InvalidResidue(char),
    TooLong(String),
}

impl FromStr for ModificationSpecificity {
    type Err = InvalidModification;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() > 2 {
            return Err(InvalidModification::TooLong(s.into()));
        }
        if let Some(rest) = s.strip_prefix('^') {
            return Ok(ModificationSpecificity::PeptideN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('$') {
            return Ok(ModificationSpecificity::PeptideC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('[') {
            return Ok(ModificationSpecificity::ProteinN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix(']') {
            return Ok(ModificationSpecificity::ProteinC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        match s.chars().next() {
            Some(c) => {
                if VALID_AA.contains(&(c as u8)) {
                    Ok(ModificationSpecificity::Residue(c as u8))
                } else {
                    Err(InvalidModification::InvalidResidue(c))
                }
            }
            None => Err(InvalidModification::Empty),
        }
    }
}

pub fn validate_mods(input: Option<HashMap<String, f32>>) -> HashMap<ModificationSpecificity, f32> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, mass) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, mass);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!("Invalid modification string: unrecognized residue ({})", c)
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

pub fn validate_var_mods(
    input: Option<HashMap<String, Vec<f32>>>,
) -> HashMap<ModificationSpecificity, Vec<f32>> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, mass) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, mass);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Skipping invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!(
                        "Skipping invalid modification string: unrecognized residue ({})",
                        c
                    )
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Skipping invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parse_modifications() {
        use InvalidModification::*;
        use ModificationSpecificity::*;
        assert_eq!("[".parse::<ModificationSpecificity>(), Ok(ProteinN(None)));
        assert_eq!(
            "[M".parse::<ModificationSpecificity>(),
            Ok(ProteinN(Some(b'M')))
        );
        assert_eq!(
            "]M".parse::<ModificationSpecificity>(),
            Ok(ProteinC(Some(b'M')))
        );
        assert_eq!("M".parse::<ModificationSpecificity>(), Ok(Residue(b'M')));
        assert_eq!(
            "Z".parse::<ModificationSpecificity>(),
            Err(InvalidResidue('Z'))
        );
    }
}
