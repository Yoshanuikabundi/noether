use crate::units::f64;
use crate::result::*;

use crate::boundaries::BoundaryConditions;

pub mod simple;

#[derive(Debug, Clone, Copy)]
/// A cutoff for computing a non-bonded potential
pub enum Cutoff {
    At(f64::Length),
    None,
}

impl From<Cutoff> for Result<f64::Length> {
    fn from(cutoff: Cutoff) -> Result<f64::Length> {
        match cutoff {
            Cutoff::At(l) => Ok(l),
            Cutoff::None => Err(CutoffRequired),
        }
    }
}

impl From<Cutoff> for Option<f64::Length> {
    fn from(cutoff: Cutoff) -> Option<f64::Length> {
        match cutoff {
            Cutoff::At(l) => Some(l),
            Cutoff::None => None,
        }
    }
}

/// A pair of atom indices and the squared distance between them
pub type AtomPair = ([usize; 2], f64::Area);

/// Trait for producing pairlists to compute potentials from
pub trait Pairlist
    where
        Self: std::marker::Sized
{
    /// Update the pairlist based on the positions of atoms
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        cutoff: Cutoff,
        boundaries: &impl BoundaryConditions
    ) -> Result<()>;

    /// Iterate over the atom pairs in the pairlist
    fn iter(&self) -> std::slice::Iter<AtomPair>;
}


