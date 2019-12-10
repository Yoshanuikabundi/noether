use crate::units::f64;
use crate::result::*;

use crate::boundaries::BoundaryConditions;

pub mod simple;
pub mod verlet;

pub type Cutoff = Option<f64::Length>;

/// A pair of atom indices and the squared distance between them
pub type AtomPair = ([usize; 2], f64::Area);

/// Trait for producing neighbourlists to compute potentials from
pub trait Neighbourlist<B: BoundaryConditions>
{
    /// Regenerate the neighbourlist based on the positions of atoms
    ///
    /// This is an expensive computation guaranteed to produce a
    /// correct neighbourlist (within whatever error bounds the
    /// implementation provides. Check the implementaton docs
    /// for details.
    fn regenerate(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &B
    ) -> Result<()>;

    /// Update the neighbourlist based on the positions of atoms
    ///
    /// This is a cheap computation that usually makes weaker
    /// guarantees than regenerate. Check the implementation
    /// docs for details.
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &B
    ) -> Result<()>;

    /// Iterate over the atom pairs in the neighbourlist
    fn iter(&self) -> std::slice::Iter<AtomPair>;

    /// Get the parameters of the neighbourlist
    fn neighbourlist_params(&self) -> NeighbourlistParams;
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum NeighbourlistParams {
    NonbondedCutoff(f64::Length)
}
