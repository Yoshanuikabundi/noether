use crate::boundaries::BoundaryConditions;
use crate::pairlist::Pairlist;
use crate::units::f64;
use crate::result::*;

pub mod lj;


/// Trait for types that represent particular potentials
/// in a topology
trait Potential
{
    /// Compute the potential
    ///
    /// Potentials must figure out their own access to
    /// the simulation's state
    fn compute_potential(&self) -> f64::Energy;

    /// Get the number of atoms accounted for by the potential.
    ///
    /// Should always be the same number (in the same order) as
    /// the topology
    fn num_atoms(&self) -> usize;
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
pub struct Topology<'a> {
    /// Every atom has a unique atom_name, index not in Atom
    atom_names: Vec<String>,
    /// All the potentials in the topology
    potentials: Vec<Box<dyn Potential + 'a>>
}

#[allow(clippy::len_without_is_empty)]
impl<'a> Topology<'a> {
    fn new(
        atom_names: &[String],
        potentials: Vec<Box<dyn Potential + 'a>>
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
    pub fn lj_fluid<P: 'a + Pairlist<B>, B: 'a + BoundaryConditions>(
        name: String,
        num_atoms: usize,
        sigma: f64::Length,
        epsilon: f64::Energy,
        cutoff: f64::Length
    ) -> Result<(Self, P)> {
        let (lj_potential, pairlist) = lj::LjPotential::lj_fluid(
            num_atoms,
            sigma,
            epsilon,
            cutoff
        )?;

        Ok((
            Topology::new(
                &vec![name; num_atoms],
                vec![Box::new(lj_potential)],
            )?,
            pairlist
        ))
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

    pub fn compute_potential(&self) -> f64::Energy {
        self
            .potentials
            .iter()
            .fold(
                0.0 * f64::KJPERMOL,
                |acc, pot| {
                    acc + pot.compute_potential()
                }
            )
    }
}
