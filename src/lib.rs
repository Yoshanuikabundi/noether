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

/// Error handling for this crate
pub mod result;
use result::*;
pub use result::FriendlyResult;

use boundaries::BoundaryConditions;

use uom::typenum::consts::*;

/// Current version of noether
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Compute the LJ potential of a pair of atoms
///
/// $$ V = 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] $$
///
/// # Arguments
///
/// * `r_squared: Area`: Square of the distance between particles ($r^2$)
/// * `sigma_squared: Length`: Finite distance squared at which potential is zero ($\sigma^2$)
/// * `epsilon: Energy`: Depth of potential well ($\epsilon$)
#[inline]
pub fn lj_potential(
    r_squared: f64::Area,
    sigma_squared: f64::Area,
    epsilon: f64::Energy
) -> f64::Energy {
    let six = (sigma_squared / r_squared).powi(P3::new());
    let twelve = six.powi(P2::new());

    4.0 * epsilon * (twelve - six)
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
/// * `sigmas`: Sigma ($\sigma$) values for each atom in `frames`; shape $(n_\mathrm{atoms})$
/// * `epsilons`: Epsilon ($\epsilon$) values for each atom in `frames`; shape $(n_\mathrm{atoms})$
/// * `cutoff`: The cutoff beyond which the potential is zero. Potential is shifted by a constant $V_\mathrm{cutoff}$ to be zero at this distance.
/// * `boundaries`: Boundary conditions; determines how distance is calculated
pub fn lj_from_positions(
    frames: &[Vec<[f64::Length; 3]>],
    sigmas: &[f64::Length],
    epsilons: &[f64::Energy],
    cutoff: f64::Length,
    boundaries: &impl BoundaryConditions
) -> Result<Vec<f64::Energy>> {
    if sigmas.len() != epsilons.len() {
        return Err(ValueError("lj_from_positions requires equal numbers of sigmas and epsilons"));
    }

    let mut energies = Vec::new();
    let cutoff_squared = cutoff.powi(P2::new());

    for frame in frames {
        let mut energy = 0.0 * f64::KJPERMOL;

        if frame.len() != sigmas.len() {
            return Err(ValueError("lj_from_positions requires frames with equal numbers of atoms to sigmas and epsilons"));
        }

        for ([i, j], r_squared) in boundaries.construct_pairlist(frame, cutoff)? {
            let sigma_squared = sigmas[i] * sigmas[j];
            let epsilon = (epsilons[i] * epsilons[j]).sqrt();

            energy += lj_potential(r_squared, sigma_squared, epsilon) - lj_potential(cutoff_squared, sigma_squared, epsilon);
        }

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

    #[test]
    fn lj_potential_2atoms() {
        let positions = io::read_positions(
            "test_targets/2_atoms/2_atoms_frommax.trr"
        ).unwrap();

        let energies = lj_from_positions(
            &positions,
            &[0.3405 * NM; 2],
            &[1.0000 * KJPERMOL; 2],
            1.2 * NM,
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

        let energies = lj_from_positions(
            &positions,
            &[0.3405 * NM; 100],
            &[1.0000 * KJPERMOL; 100],
            1.2 * NM,
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
