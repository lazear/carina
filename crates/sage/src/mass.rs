use std::{iter::Sum, ops::Mul};

use serde::{Deserialize, Serialize};

pub const H2O: f32 = 18.010565;
pub const PROTON: f32 = 1.0072764;
pub const NEUTRON: f32 = 1.00335;
pub const NH3: f32 = 17.026548;

#[derive(Copy, Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd)]
#[serde(rename_all = "lowercase")]
pub enum Tolerance {
    Ppm(f32, f32),
    Pct(f32, f32),
    Da(f32, f32),
}

impl Tolerance {
    /// Compute the (`lower`, `upper`) window (in Da) for for a monoisotopic
    /// mass and a given tolerance
    pub fn bounds(&self, center: f32) -> (f32, f32) {
        match self {
            Tolerance::Ppm(lo, hi) => {
                let delta_lo = center * lo / 1_000_000.0;
                let delta_hi = center * hi / 1_000_000.0;
                (center + delta_lo, center + delta_hi)
            }
            Tolerance::Pct(lo, hi) => {
                let delta_lo = center * lo / 100.0;
                let delta_hi = center * hi / 100.0;
                (center + delta_lo, center + delta_hi)
            }
            Tolerance::Da(lo, hi) => (center + lo, center + hi),
        }
    }

    pub fn contains(&self, center: f32, rhs: f32) -> bool {
        let (lo, hi) = self.bounds(center);
        rhs >= lo && rhs <= hi
    }

    pub fn ppm_to_delta_mass(center: f32, ppm: f32) -> f32 {
        ppm * center / 1_000_000.0
    }
}

impl Mul<f32> for Tolerance {
    type Output = Tolerance;

    fn mul(self, rhs: f32) -> Self::Output {
        match self {
            Tolerance::Ppm(lo, hi) => Tolerance::Ppm(lo * rhs, hi * rhs),
            Tolerance::Pct(lo, hi) => Tolerance::Pct(lo * rhs, hi * rhs),
            Tolerance::Da(lo, hi) => Tolerance::Da(lo * rhs, hi * rhs),
        }
    }
}

pub const VALID_AA: [u8; 22] = [
    b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S',
    b'T', b'V', b'W', b'Y', b'U', b'O',
];

pub const MONOISOTOPIC_MASSES: [f32; 26] = [
    71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891, 113.08406, 0.0,
    128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276, 128.05858, 156.1011, 87.03203,
    101.04768, 150.95363, 99.06841, 186.07932, 0.0, 163.06332, 0.0,
];

pub const fn monoisotopic(aa: u8) -> f32 {
    if aa.is_ascii_uppercase() {
        MONOISOTOPIC_MASSES[(aa - b'A') as usize]
    } else {
        0.0
    }
}

pub const fn composition(aa: u8) -> Composition {
    match aa {
        b'A' => Composition::new(3, 2, 0),
        b'R' => Composition::new(6, 2, 0),
        b'N' => Composition::new(4, 3, 0),
        b'D' => Composition::new(4, 4, 0),
        b'C' => Composition::new(3, 2, 1),
        b'E' => Composition::new(5, 4, 0),
        b'Q' => Composition::new(5, 3, 0),
        b'G' => Composition::new(2, 2, 0),
        b'H' => Composition::new(6, 2, 0),
        b'I' => Composition::new(6, 2, 0),
        b'L' => Composition::new(6, 2, 0),
        b'K' => Composition::new(6, 2, 0),
        b'M' => Composition::new(5, 2, 1),
        b'F' => Composition::new(9, 2, 0),
        b'P' => Composition::new(5, 2, 0),
        b'S' => Composition::new(3, 3, 0),
        b'T' => Composition::new(4, 3, 0),
        b'W' => Composition::new(11, 2, 0),
        b'Y' => Composition::new(9, 3, 0),
        b'V' => Composition::new(5, 2, 0),
        b'U' => Composition::new(3, 2, 0),
        b'O' => Composition::new(12, 3, 0),
        _ => Composition::new(0, 0, 0),
    }
}

#[derive(Clone, Debug)]
pub struct Composition {
    pub carbon: u16,
    pub sulfur: u16,
}

impl Composition {
    pub const fn new(carbon: u16, _oxygen: u8, sulfur: u16) -> Self {
        Self { carbon, sulfur }
    }
}

impl Sum for Composition {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut comp = Composition::new(0, 0, 0);
        for i in iter {
            comp.carbon += i.carbon;
            // comp.oxygen += i.oxygen;
            comp.sulfur += i.sulfur;
        }
        comp
    }
}

#[cfg(test)]
mod test {
    use crate::mass::monoisotopic;

    use super::{Tolerance, VALID_AA};

    #[test]
    fn smoke() {
        for ch in VALID_AA {
            assert!(monoisotopic(ch) > 0.0);
        }
    }

    #[test]
    fn tolerances() {
        assert_eq!(
            Tolerance::Ppm(-10.0, 20.0).bounds(1000.0),
            (999.99, 1000.02)
        );
        assert_eq!(
            Tolerance::Ppm(-10.0, 10.0).bounds(487.0),
            (486.99513, 487.00487)
        );
        assert_eq!(
            Tolerance::Ppm(-50.0, 50.0).bounds(1000.0),
            (999.95, 1000.05)
        );
    }
}
