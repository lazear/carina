use std::fmt::Write;

use serde::{Deserialize, Serialize};

pub const H2O: f32 = 18.010565;
pub const PROTON: f32 = 1.007_276_4;
pub const NH3: f32 = 17.026548;

#[derive(Copy, Clone, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Tolerance {
    Ppm(f32),
    Th(f32),
}

impl Tolerance {
    /// Compute the (`lower`, `upper`) window (in Da) for for a monoisotopic
    /// mass and a given tolerance
    pub fn bounds(&self, center: f32) -> (f32, f32) {
        match self {
            Tolerance::Ppm(ppm) => {
                let delta = center * ppm / 1_000_000.0;
                (center - delta, center + delta)
            }
            Tolerance::Th(th) => (center - th, center + th),
        }
    }
}

pub trait Mass {
    fn monoisotopic(&self) -> f32;
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize)]
pub enum Residue {
    // Standard amino acid residue
    Just(char),
    // Amino acid residue with a mass modification
    Mod(char, f32),
}

impl Mass for Residue {
    fn monoisotopic(&self) -> f32 {
        match self {
            Residue::Just(c) => c.monoisotopic(),
            Residue::Mod(c, m) => c.monoisotopic() + m,
        }
    }
}

pub const VALID_AA: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

impl Mass for char {
    fn monoisotopic(&self) -> f32 {
        match self {
            'A' => 71.037_12,
            'R' => 156.101_1,
            'N' => 114.042_93,
            'D' => 115.026_94,
            'C' => 103.009_186,
            'E' => 129.042_59,
            'Q' => 128.058_58,
            'G' => 57.021_465,
            'H' => 137.058_91,
            'I' => 113.084_06,
            'L' => 113.084_06,
            'K' => 128.094_96,
            'M' => 131.040_48,
            'F' => 147.068_42,
            'P' => 97.052_765,
            'S' => 87.032_03,
            'T' => 101.047_676,
            'W' => 186.079_32,
            'Y' => 163.063_32,
            'V' => 99.068_41,
            _ => unreachable!("BUG: invalid amino acid"),
        }
    }
}

impl std::fmt::Display for Residue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Residue::Just(c) => f.write_char(*c),
            Residue::Mod(c, m) => write!(f, "{}({})", c, m),
        }
    }
}

#[cfg(test)]
mod test {
    use super::{Mass, Tolerance, VALID_AA};

    #[test]
    fn smoke() {
        for ch in VALID_AA {
            assert!(ch.monoisotopic() > 0.0);
        }
    }

    #[test]
    fn tolerances() {
        assert_eq!(Tolerance::Ppm(10.0).bounds(1000.0), (999.99, 1000.01));
        assert_eq!(Tolerance::Ppm(10.0).bounds(487.0), (486.99513, 487.00487));

        assert_eq!(Tolerance::Ppm(50.0).bounds(1000.0), (999.95, 1000.05));
    }
}
