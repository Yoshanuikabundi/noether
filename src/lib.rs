#[macro_use]
extern crate uom;

#[cfg(test)]
#[macro_use]
extern crate more_asserts;

use std::collections::HashMap;

/// Keep our units straight (Based on `uom` crate)
pub mod units;
use units::f64;

/// Read from and write to disk
pub mod io;

/// Store a simulation's boundary conditions and perform relevant calculations (eg, distances)
pub mod boundaries;

/// Store and compute the simulations pairlist
pub mod pairlist;
use pairlist::{Pairlist, PairlistParams};
use pairlist::simple::SimplePairlist;

/// Store the topology of a simulation
///
/// This module has to be super-optimized and follows an ECS-style layout
pub mod topology;
use crate::topology::*;

/// Error handling for this crate
pub mod result;
use result::*;
pub use result::FriendlyResult;

use boundaries::BoundaryConditions;

/// Current version of noether
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Compute the potential energy of uncorrelated frames of the system described by `top` and `boundaries`
///
/// $$ V_\mathrm{frame} = \sum_\mathrm{atoms} V_\mathrm{atom}(r)
///
/// Returns a list of values $V_\mathrm{frame}$ for each frame.
///
/// The pairlist is regenerated every frame. This is very expensive, but allows uncorrelated frames to be computed.
///
/// # Arguments
///
/// * `frames`: List of frames, which are lists of atomic position 3-vectors; shape $(n_\mathrm{frames}, n_\mathrm{atoms}, 3)$
/// * `top`: Topology of the system
/// * `boundaries`: Boundary conditions; determines how distance is calculated
pub fn pot_from_positions<B: BoundaryConditions>(
    frames: &[Vec<[f64::Length; 3]>],
    top: &Topology<B>,
    boundaries: &B
) -> Result<Vec<f64::Energy>> {
    let mut energies = Vec::with_capacity(frames.len());

    // Construct the pairlists for the topology and store them together
    let mut pairlists_dict: HashMap<String, SimplePairlist<B>> = HashMap::new();
    for PairlistParams::NonbondedCutoff(cutoff) in top.iter_pairlists() {
        pairlists_dict
            .entry(format!("{:?}", cutoff))
            .or_insert(SimplePairlist::new(Some(cutoff))?);
    }

    for frame in frames {
        if frame.len() != top.len() {
            return Err(PositionTopologyMismatch);
        }

        // Update all the pairlists
        pairlists_dict
            .iter_mut()
            .try_for_each(|(_, pairlist)| pairlist.regenerate(frame, boundaries))?
        ;

        // Construct a list of pairlists that matches the topology
        let mut pairlists: Vec<&dyn Pairlist<B>> = vec![];
        for PairlistParams::NonbondedCutoff(cutoff) in top.iter_pairlists() {
            let pairlist = pairlists_dict.get(&format!("{:?}", cutoff)).unwrap();
            boundaries.pairlist_checks(pairlist)?;
            pairlists.push(pairlist)
        }

        let energy = top.compute_potential(&pairlists);

        energies.push(energy);
    }

    Ok(energies)
}

#[cfg(test)]
mod tests {
    use crate::units::f64::KJPERMOL;
    use crate::units::f64::NM;
    use crate::io;
    use crate::pot_from_positions;
    use crate::boundaries::*;
    use crate::topology::Topology;

    #[test]
    fn lj_potential_2atoms() {
        let positions = io::read_positions(
            "test_targets/2_atoms/2_atoms_frommax.trr"
        ).unwrap();

        let topol = Topology::lj_fluid(
            "Me".to_string(),
            2,
            0.3405 * NM,
            1.0000 * KJPERMOL,
            1.2 * NM
        ).unwrap();

        let energies = pot_from_positions(
            &positions,
            &topol,
            &NoBounds
        ).unwrap();

        let energies_ref = io::read_xvg(
            "test_targets/2_atoms/2_atoms_frommax.xvg",
            1
        ).unwrap();
        for (a, b) in energies.iter().zip(energies_ref.iter()) {
            let b: f64 = b.parse().unwrap();
            assert_lt!(*a - b * KJPERMOL, 1e-5 * KJPERMOL);
        }
    }

    #[test]
    fn lj_potential_100atoms() {
        let positions = io::read_positions(
            "test_targets/100_atoms/100_atoms.trr"
        ).unwrap();

        let topol = Topology::lj_fluid(
            "Me".to_string(),
            100,
            0.3405 * NM,
            1.0000 * KJPERMOL,
            1.2 * NM
        ).unwrap();


        let energies = pot_from_positions(
            &positions,
            &topol,
            &Pbc::cubic(5.0*NM)
        ).unwrap();

        let energies_ref = io::read_xvg(
            "test_targets/100_atoms/100_atoms.xvg",
            1
        ).unwrap();
        for (a, b) in energies.iter().zip(energies_ref.iter()) {
            let b: f64 = b.parse().unwrap();
            assert_lt!(*a - b * KJPERMOL, 2e-4 * KJPERMOL);
        }
    }
}
