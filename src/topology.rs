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
pub struct Topology {
    atoms: Vec<Atom>,
    atom_names: Vec<String>,
    /// lj_params is an atomtype in GROMACS
    lj_params: Vec<LjParams>,
}

impl Topology {
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

    pub fn iter_lj(&self) -> TopolIterator<LjParams> {
        TopolIterator::new(self, &self.lj_params)
    }
}

/// An atom.
#[derive(Clone, Debug, Copy)]
struct Atom {
    /// 66k unique Lennard-Jones parameters
    lj_params: u16
}

/// A pair of LJ parameters.
pub struct LjParams {
    pub sigma: f64::Length,
    pub epsilon: f64::Energy
}

pub struct TopolIterator<'a, T> {
    topol: &'a Topology,
    vector: &'a Vec<T>,
    index: usize
}

impl<'a, T> TopolIterator<'a, T> {
    #[allow(clippy::ptr_arg)]
    fn new(topol: &'a Topology, vector: &'a Vec<T>) -> Self {
        TopolIterator {
            topol,
            vector,
            index: 0
        }
    }
}

impl<T> Iterator for TopolIterator<'_, T>
    where T: Copy
{
    type Item = T;

    fn next(&mut self) -> Option<T>  {
        let atom = self.topol.atoms[self.index];
        self.index += 1;
        self.vector.get(atom.lj_params as usize).copied()
    }

}
