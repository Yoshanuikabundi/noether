#[macro_use]
extern crate uom;

#[cfg(test)]
#[macro_use]
extern crate more_asserts;

/// Keep our units straight (Based on `uom` crate)
pub mod units;
use units::f64;

/// Read from and write to disk
pub mod io;

/// Store a simulation's boundary conditions and perform relevant calculations (eg, distances)
pub mod boundaries;

/// Store the topology of a simulation
///
/// This module has to be super-optimized and follows an ECS-style layout
pub mod topology;
use crate::topology::*;

/// Error handling for this crate
pub mod result;
use result::*;
pub use result::FriendlyResult;

use boundaries::{BoundaryConditions, Pairlist};

/// Current version of noether
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

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
pub fn frame_lj(
    top: &Topology,
    pairlist: Pairlist
) -> f64::Energy {
    pairlist
        .iter()
        .fold(
            0.0 * f64::KJPERMOL,
            |energy, ([i, j], r_squared)| {
                energy
                + top
                    .get_lj_pair(*i, *j)
                    .unwrap()
                    .lj_potential(*r_squared)
            }
        )

}

/// Compute the Lennard-Jones energy of uncorrelated frames of a fluid of particles
///
/// $$ V_\mathrm{frame} = \sum_\mathrm{atoms} 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] - V_\mathrm{cutoff} $$
///
/// Returns a list of values $V_\mathrm{frame}$ for each frame.
///
/// The pairlist is regenerated every frame. This is very expensive, but allows uncorrelated frames to be computed.
///
/// # Arguments
///
/// * `frames`: List of frames, which are lists of atomic position 3-vectors; shape $(n_\mathrm{frames}, n_\mathrm{atoms}, 3)$
/// * `top`: Topology of the system
/// * `cutoff`: The cutoff beyond which the potential is zero. Potential is shifted by a constant $V_\mathrm{cutoff}$ to be zero at this distance.
/// * `boundaries`: Boundary conditions; determines how distance is calculated
pub fn lj_from_positions(
    frames: &[Vec<[f64::Length; 3]>],
    top: &Topology,
    boundaries: &impl BoundaryConditions
) -> Result<Vec<f64::Energy>> {
    let mut energies = Vec::new();

    for frame in frames {
        if frame.len() != top.len() {
            return Err(ValueError(concat!(
                "lj_from_positions requires topology and ",
                "frames to have equal numbers of atoms"
            )));
        }

        let pairlist = boundaries.construct_pairlist(frame, top)?;

        let energy = frame_lj(
            top,
            pairlist
        );

        energies.push(energy);
    }

    Ok(energies)
}

#[cfg(test)]
mod tests {
    use crate::units::f64::KJPERMOL;
    use crate::units::f64::NM;
    use crate::io;
    use crate::lj_from_positions;
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

        let energies = lj_from_positions(
            &positions,
            &topol,
            &NoBounds,
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

        let energies = lj_from_positions(
            &positions,
            &topol,
            &Pbc::cubic(5.0*NM),
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
