use crate::units::f64;

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
/// their own vectors, and pointing to them in the `Atom`
/// struct. This way, when we compute a potential the
/// memory streams through that vector. This also lets
/// us minimise the size of the `Atom` struct in two ways:
///
/// 1. Combine parameters that are used together into their
/// own structs, even if this means duplicating values in
/// memory
/// 2. Use indexes smaller than `usize`. This means we need
/// to group different atoms with the same parameters
/// together. u16 should be a good starting point - it gives us
/// four-fold smaller Atoms while still allowing ~66k different
/// parameters. Since systems are mostly copies of identical
/// atoms (eg, water), this should be plenty.
///
/// We're using vectors here because we don't know the sizes at
/// compile time; Topology shouldn't ever have to reallocate its
/// fields.
pub struct Topology {
    /// The atoms in the topology
    ///
    /// Contains indices for other parameters
    atoms: Vec<Atom>,
    /// Every atom has a unique atom_name, index not in Atom
    atom_names: Vec<String>,
    /// lj_params is an atomtype in GROMACS
    ///
    /// Indexed by u16
    lj_params: Vec<LjParams>,
}

#[allow(clippy::len_without_is_empty)]
impl Topology {
    /// Create a topology for a homogenous Lennard Jones fluid
    ///
    /// # Arguments
    ///
    /// * `name: String` - Name to use for the atoms
    /// * `num_atoms: usize` - Number of atoms in the topology
    /// * `sigma: Length` - Value of sigma to use for all atoms
    /// * `epsilon: Energy` - Value of epsilon to use for all atoms
    pub fn lj_fluid(
        name: String,
        num_atoms: usize,
        sigma: f64::Length,
        epsilon: f64::Energy
    ) -> Topology {
        let lj_param = LjParams { sigma, epsilon };
        let atom = Atom { lj_params: 0 };
        Topology {
            atoms: vec![atom; num_atoms],
            atom_names: vec![name; num_atoms],
            lj_params: vec![lj_param],
        }
    }

    /// Get the LJ parameters of atom `index`
    pub fn get_lj_params(&self, index: usize) -> Option<&LjParams> {
        let index = self.atoms.get(index)?.lj_params as usize;
        self.lj_params.get(index)
    }

    /// Get the name of atom `index`
    pub fn get_atom_names(&self, index: usize) -> Option<&String> {
        self.atom_names.get(index)
    }

    /// Iterate over the Lennard-Jones parameters of the atoms
    pub fn iter_lj_params(&self) -> TopolIterator<LjParams> {
        TopolIterator {
            iterator: self.atoms.iter(),
            vector: &self.lj_params
        }
    }

    /// Iterate over the atom names
    pub fn iter_atom_names(&self) -> std::slice::Iter<String> {
        self.atom_names.iter()
    }

    /// Return the number of atoms in the topology
    pub fn len(&self) -> usize {
        self.atoms.len()
    }
}

/// An atom.
#[derive(Clone, Debug)]
struct Atom {
    /// 66k unique Lennard-Jones parameters
    lj_params: u16
}

/// A pair of LJ parameters.
#[derive(Clone, Debug)]
pub struct LjParams {
    pub sigma: f64::Length,
    pub epsilon: f64::Energy
}

pub struct TopolIterator<'a, T> {
    iterator: std::slice::Iter<'a, Atom>,
    vector: &'a [T]
}

impl<'a> Iterator for TopolIterator<'a, LjParams> {
    type Item = &'a LjParams;

    fn next(&mut self) -> Option<Self::Item>  {
        let atom = self.iterator.next()?;
        self.vector.get(atom.lj_params as usize)
    }
}
