use crate::boundaries::BoundaryConditions;
use crate::neighbourlist::{Neighbourlist, NeighbourlistParams};
use crate::units::f64;
use crate::result::*;

pub mod lj;


/// Trait for types that represent particular potentials
/// in a topology
trait Potential<B: BoundaryConditions>
{

    /// Get the number of atoms accounted for by the potential.
    ///
    /// Should always be the same number (in the same order) as
    /// the topology
    fn num_atoms(&self) -> usize;

    /// Get the neighbourlist parameters of the potential
    fn neighbourlist_params(&self) -> NeighbourlistParams;

    /// Compute the potential for a single pair of atom indices and squared distance
    fn potential(&self, i: usize, j: usize, r_squared: f64::Area) -> f64::Energy;

    /// Compute the potential for the entire neighbourlist
    fn compute_potential(&self, neighbourlist: &dyn Neighbourlist<B>) -> Result<f64::Energy> {
        if neighbourlist.neighbourlist_params() != self.neighbourlist_params() {
            return Err(NeighbourlistNotCompatible)
        }

        Ok(neighbourlist
            .iter()
            .fold(
                0.0 * f64::KJPERMOL,
                |energy, ([i, j], r_squared)| {
                    energy
                    + self.potential(*i, *j, *r_squared)
                }
            )
        )
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
pub struct Topology<B: BoundaryConditions> {
    /// Every atom has a unique atom_name, index not in Atom
    atom_names: Vec<String>,
    /// All the potentials in the topology
    potentials: Vec<Box<dyn Potential<B>>>
}

#[allow(clippy::len_without_is_empty)]
impl<B: BoundaryConditions> Topology<B> {
    fn new(
        atom_names: &[String],
        potentials: Vec<Box<dyn Potential<B>>>
    ) -> Result<Self>{
        for potential in potentials.iter() {
            // Check that each potential has the right number of Atoms
            if potential.num_atoms() != atom_names.len() {
                return Err(IllegalTopology);
            }
        }

        Ok(Topology {
            atom_names: atom_names.to_vec(),
            potentials
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
        )?;

        Ok(Topology::new(
            &vec![name; num_atoms],
            vec![Box::new(lj_potential)],
        )?)
    }

    /// Get the name of atom `index`
    pub fn get_atom_name(&self, index: usize) -> Option<&String> {
        self.atom_names.get(index)
    }

    /// Iterate over the atom names
    pub fn iter_atom_names(&self) -> std::slice::Iter<String> {
        self.atom_names.iter()
    }

    /// Iterate over the potential's neighbourlist parameters
    pub fn iter_neighbourlists<'a>(&'a self) -> Box<dyn Iterator<Item=NeighbourlistParams> + 'a> {
        Box::new(self
            .potentials
            .iter()
            .map(|pot| pot.neighbourlist_params())
        )
    }

    /// Return the number of atoms in the topology
    pub fn len(&self) -> usize {
        self.atom_names.len()
    }

    /// Compute the total potential of the topology from the neighbourlists
    ///
    /// # Arguments
    ///
    /// * `neighbourlists`: Vector of references to the neighbourlists, in the same order as the potentials.
    ///
    /// # Panics
    ///
    /// Panics if the vector of neighbourlists doesn't match the potentials
    pub fn compute_potential(&self, neighbourlists: &Vec<&dyn Neighbourlist<B>>) -> f64::Energy {
        let panicstr = "Neighbourlists must match potentials";

        if neighbourlists.len() != self.potentials.len() {
            panic!(panicstr)
        }

        self
            .potentials
            .iter()
            .zip(neighbourlists)
            .fold(
                0.0 * f64::KJPERMOL,
                |acc, (pot, neighbourlist)| {
                    acc + pot.compute_potential(*neighbourlist).expect(panicstr)
                }
            )
    }
}
