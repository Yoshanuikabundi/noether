use crate::units::f64;
use crate::result::*;

use crate::boundaries::BoundaryConditions;

pub mod simple;

pub type Cutoff = Option<f64::Length>;

/// A pair of atom indices and the squared distance between them
pub type AtomPair = ([usize; 2], f64::Area);

/// Trait for producing pairlists to compute potentials from
pub trait Pairlist<B: BoundaryConditions>
{
    /// Update the pairlist based on the positions of atoms
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &B
    ) -> Result<()>;

    /// Iterate over the atom pairs in the pairlist
    fn iter(&self) -> std::slice::Iter<AtomPair>;

    fn new(cutoff: Cutoff) -> Result<Self> where Self: Sized;

    fn cutoff(&self) -> Cutoff;
}


