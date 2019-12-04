use crate::pairlist::{Pairlist, Cutoff};
use crate::units::f64;
use crate::result::*;

pub mod lj;


/// Trait for types that represent particular potentials
/// in a topology
trait Potential<P>
    where
        P: Pairlist
{
    /// Compute the potential, with aid of a pairlist if required
    fn compute_potential(
        &self,
        pairlist: &P
    ) -> f64::Energy;

    /// Get the number of atoms accounted for by the potential.
    ///
    /// Should always be the same number (in the same order) as
    /// the topology
    fn num_atoms(&self) -> usize;

    /// The cutoff required for the potential. None if no pairs
    /// should be computed.
    fn cutoff(&self) -> Cutoff {
        Cutoff::None
    }
}

/// Stores the parameters of all atoms as vectors
///
/// It's important that we can stream parameters into
/// the CPU, and ideally fit the entire topology in cache.
/// We could simply have a vector of Atoms, and have each
/// `Atom` store all its parameters, but then the Atom
/// struct would have to be $64 \times n_\mathrm{params}$
/// bits large, and when we computed a particular potential
/// memory access would be strided as we skipped over the
/// parameters we don't care about.
///
/// We can solve both problems by storing parameters in
/// their own vectors, and pointing to them with indices.
/// Each potential has its own vectors. The `Topology`
/// `atom_names` field forms a master list whose indices
/// are universal atom indices.
///
/// We're using vectors here because we don't know the sizes at
/// compile time; Topology shouldn't ever have to reallocate its
/// fields.
pub struct Topology<P: Pairlist> {
    /// Every atom has a unique atom_name, index not in Atom
    atom_names: Vec<String>,
    /// All the potentials in the topology
    potentials: Vec<Box<dyn Potential<P>>>,
    /// The cutoff for the overall topology
    cutoff: Cutoff
}

#[allow(clippy::len_without_is_empty)]
impl<P: Pairlist> Topology<P> {
    fn new(
        atom_names: Vec<String>,
        potentials: Vec<Box<dyn Potential<P>>>
    ) -> Result<Self>{
        let mut cutoff = Cutoff::None;
        for potential in potentials.iter() {
            // Check that each potential has the right number of Atoms
            if potential.num_atoms() != atom_names.len() {
                return Err(IllegalTopology);
            }
            // Compute the overall topology cutoff
            cutoff = match (potential.cutoff(), cutoff) {
                (Cutoff::None, Cutoff::None) => Cutoff::None,
                (Cutoff::At(l), Cutoff::None) => Cutoff::At(l),
                (Cutoff::None, Cutoff::At(l)) => Cutoff::At(l),
                (Cutoff::At(a), Cutoff::At(b)) if a == b => Cutoff::At(a),
                (Cutoff::At(_), Cutoff::At(_)) => return Err(InconsistentCutoffs)
            }
        }

        Ok(Topology {
            atom_names,
            potentials,
            cutoff
        })

    }

    /// Create a topology for a homogenous Lennard Jones fluid
    ///
    /// # Arguments
    ///
    /// * `name: String` - Name to use for the atoms
    /// * `num_atoms: usize` - Number of atoms in the topology
    /// * `sigma: Length` - Value of sigma to use for all atoms
    /// * `epsilon: Energy` - Value of epsilon to use for all atoms
    /// * `cutoff: Length` - Distance beyond which potential is zero
    pub fn lj_fluid(
        name: String,
        num_atoms: usize,
        sigma: f64::Length,
        epsilon: f64::Energy,
        cutoff: f64::Length
    ) -> Result<Self> {
        let lj_potential = lj::LjPotential::lj_fluid(
            num_atoms,
            sigma,
            epsilon,
            cutoff
        );

        Topology::new(
            vec![name; num_atoms],
            vec![Box::new(lj_potential?)]
        )
    }

    /// Get the name of atom `index`
    pub fn get_atom_names(&self, index: usize) -> Option<&String> {
        self.atom_names.get(index)
    }

    /// Iterate over the atom names
    pub fn iter_atom_names(&self) -> std::slice::Iter<String> {
        self.atom_names.iter()
    }

    /// Return the number of atoms in the topology
    pub fn len(&self) -> usize {
        self.atom_names.len()
    }

    pub fn compute_potential(&self, pairlist: &P) -> f64::Energy {
        self
            .potentials
            .iter()
            .fold(
                0.0 * f64::KJPERMOL,
                |acc, pot| {
                    acc + pot.compute_potential(pairlist)
                }
            )
    }

    pub fn cutoff(&self) -> Cutoff {
        self.cutoff
    }
}
