use crate::units::f64;
use crate::result::*;
use uom::typenum::consts::*;

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
    /// Precomputed mixed LJ parameters
    ///
    /// Indexed by u16
    lj_params_pairs: Vec<Vec<LjParamsPair>>,
    /// Cutoff really is a feature of the topology
    pub lj_cutoff: f64::Length,
}

#[allow(clippy::len_without_is_empty)]
impl Topology {
    fn new(
        atoms: Vec<Atom>,
        atom_names: Vec<String>,
        lj_params: Vec<LjParams>,
        lj_cutoff: f64::Length
    ) -> Result<Topology>{
        // Check that every atom has a name
        if atoms.len() != atom_names.len() {
            return Err(IllegalTopology);
        }

        // Check that every atom has LJ parameters
        for atom in atoms.iter() {
            if lj_params.get(atom.lj_params as usize).is_none() {
                return Err(IllegalTopology);
            }
        }

        // Pre-compute all LJ pairs
        // TODO: Replace Vec with some ndarray-style type
        let mut lj_params_pairs = Vec::with_capacity(lj_params.len());
        for a in lj_params.iter() {
            let mut lj_params_pairs_inner = Vec::with_capacity(lj_params.len());
            for b in lj_params.iter() {
                lj_params_pairs_inner.push(LjParamsPair::new(a, b, lj_cutoff));
            }
            lj_params_pairs.push(lj_params_pairs_inner);
        }

        Ok(Topology {
            atoms,
            atom_names,
            lj_params,
            lj_params_pairs,
            lj_cutoff,
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
    ) -> Result<Topology> {
        let lj_param = LjParams { sigma, epsilon };
        let atom = Atom { lj_params: 0 };
        Topology::new(
            vec![atom; num_atoms],
            vec![name; num_atoms],
            vec![lj_param],
            cutoff
        )
    }

    /// Get the LJ parameters of atom `index`
    pub fn get_lj_params(&self, index: usize) -> Option<&LjParams> {
        let index = self.atoms.get(index)?.lj_params as usize;
        self.lj_params.get(index)
    }

    /// Get the LJ pair parameters of atom `i` mixed with `j`
    pub fn get_lj_pair(&self, i: usize, j: usize) -> Option<&LjParamsPair> {
        let i = self.atoms.get(i)?.lj_params as usize;
        let j = self.atoms.get(j)?.lj_params as usize;
        self.lj_params_pairs.get(i)?.get(j)
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

    /// Compute the Lennard-Jones energy of a single frame
    ///
    /// $$ V_\mathrm{frame} = \sum_\mathrm{atoms} 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] - V_\mathrm{cutoff} $$
    ///
    /// # Arguments
    ///
    /// * `top` - Topology of the system
    /// * `pairlist` - Pairlist to compute for, includes atom indices and distances
    /// * `cutoff_squared` - the LJ cutoff squared for shift
    ///
    /// # Panics
    ///
    /// Panics if an index in the pairlist isn't present in the topology
    pub fn lj_potential(&self, pairlist: crate::boundaries::Pairlist) -> f64::Energy {
        pairlist
            .iter()
            .fold(
                0.0 * f64::KJPERMOL,
                |energy, ([i, j], r_squared)| {
                    energy
                    + self
                        .get_lj_pair(*i, *j)
                        .unwrap()
                        .lj_potential(*r_squared)
                }
            )
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

/// LJ parameters for a pair of atoms
#[derive(Clone, Debug)]
pub struct LjParamsPair {
    pub sigma_squared: f64::Area,
    pub epsilon: f64::Energy,
    pub shift: f64::Energy
}

impl LjParamsPair {
    /// Mix LjParams with geometric average for epsilon and sigma
    fn new(a: &LjParams, b: &LjParams, cutoff: f64::Length) -> LjParamsPair {
        let LjParams {
            sigma: sigma_a,
            epsilon: epsilon_a
        } = a;

        let LjParams {
            sigma: sigma_b,
            epsilon: epsilon_b
        } = b;

        let sigma_squared = *sigma_a * *sigma_b;
        let epsilon = (*epsilon_a * *epsilon_b).sqrt();
        let shift = LjParamsPair {
            sigma_squared,
            epsilon,
            shift: 0.0 * f64::KJPERMOL
        }.lj_potential(cutoff * cutoff);

        LjParamsPair {
            sigma_squared,
            epsilon,
            shift
        }
    }

    /// Compute the LJ potential of a pair of atoms
    ///
    /// $$ V = 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] - V_\mathrm{cutoff} $$
    ///
    /// # Arguments
    ///
    /// * `r_squared: Area`: Square of the distance between particles ($r^2$)
    ///
    /// # Other values
    ///
    /// * `sigma_squared: Length`: Finite distance squared at which potential is zero ($\sigma^2$)
    /// * `epsilon: Energy`: Depth of potential well ($\epsilon$)
    /// * `shift: Energy`: Potential value at cutoff ($V_\mathrm{cutoff}$)
    #[inline]
    pub fn lj_potential(&self, r_squared: f64::Area) -> f64::Energy {
        let six = (self.sigma_squared / r_squared).powi(P3::new());
        let twelve = six.powi(P2::new());

        4.0 * self.epsilon * (twelve - six) - self.shift
    }
}

/// Iterator struct for iterating over topology parameters
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

