use crate::boundaries::BoundaryConditions;
use crate::pairlist::{Pairlist};
use crate::units::f64;
use crate::result::*;
use uom::typenum::consts::*;
use super::Potential;

/// u16 means ~66k unique LJ parameters
type Atom = u16;

pub struct LjPotential<'a, P: 'a + Pairlist<B>, B: BoundaryConditions> {
    /// LJ index for each atom.
    atoms: Vec<Atom>,
    /// lj_params is an atomtype in GROMACS
    _lj_params: Vec<LjParams>,
    /// Precomputed mixed LJ parameters
    lj_params_pairs: Vec<Vec<LjParamsPair>>,
    /// Cutoff really is a feature of the topology
    lj_cutoff: f64::Length,
    /// The pairlist
    pairlist: &'a P,
    _marker: std::marker::PhantomData<B>
}

impl<'a, P: 'a +  Pairlist<B>, B: BoundaryConditions> LjPotential<'a, P, B> {
    fn new(
        atoms: Vec<Atom>,
        lj_params: Vec<LjParams>,
        lj_cutoff: f64::Length
    ) -> Result<(Self, P)> {
        let pairlist = P::new(Some(lj_cutoff))?;
        let lj_pot = Self::with_pairlist(
            atoms,
            lj_params,
            &pairlist
        )?;
        Ok((lj_pot, pairlist))
    }

    fn with_pairlist(
        atoms: Vec<Atom>,
        lj_params: Vec<LjParams>,
        pairlist: &'a P
    ) -> Result<Self> {
        // Check that every atom has LJ parameters
        for atom in atoms.iter() {
            if lj_params.get(*atom as usize).is_none() {
                return Err(IllegalTopology);
            }
        }

        let lj_cutoff = match pairlist.cutoff() {
            Some(f) => f,
            None => return Err(CutoffRequired)
        };

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

        Ok(Self {
            atoms,
            _lj_params: lj_params,
            lj_params_pairs,
            lj_cutoff,
            pairlist,
            _marker: std::marker::PhantomData
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
        num_atoms: usize,
        sigma: f64::Length,
        epsilon: f64::Energy,
        cutoff: f64::Length
    ) -> Result<(Self, P)> {
        let lj_param = LjParams { sigma, epsilon };
        let atom = 0;
        Self::new(
            vec![atom; num_atoms],
            vec![lj_param],
            cutoff
        )
    }

    /// Get the LJ pair parameters of atom `i` mixed with `j`
    #[inline]
    fn get_lj_pair(&self, i: usize, j: usize) -> Option<&LjParamsPair> {
        let i = *self.atoms.get(i)? as usize;
        let j = *self.atoms.get(j)? as usize;
        self.lj_params_pairs.get(i)?.get(j)
    }
}

impl<P: Pairlist<B>, B: BoundaryConditions> Potential for LjPotential<'_, P, B> {
    /// Compute the Lennard-Jones energy of a single frame
    ///
    /// $$ V_\mathrm{frame} = \sum_\mathrm{atoms} 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] - V_\mathrm{cutoff} $$
    ///
    /// # Arguments
    ///
    /// * `pairlist` - Pairlist to compute for, includes atom indices and distances
    ///
    /// # Panics
    ///
    /// Panics if an index in the pairlist isn't present in the topology
    fn compute_potential(&self) -> f64::Energy {
        self.pairlist
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

    /// Return the number of atoms in the topology
    fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
}

/// A pair of LJ parameters.
#[derive(Clone, Debug)]
struct LjParams {
    sigma: f64::Length,
    epsilon: f64::Energy
}

/// LJ parameters for a pair of atoms
#[derive(Clone, Debug)]
struct LjParamsPair {
    sigma_squared: f64::Area,
    epsilon: f64::Energy,
    shift: f64::Energy
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


