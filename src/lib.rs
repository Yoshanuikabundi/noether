#[macro_use]
extern crate uom;

#[cfg(test)]
#[macro_use]
extern crate more_asserts;

pub mod units;
pub mod io;
pub mod boundaries;

use boundaries::BoundaryConditions;

use uom::typenum::consts::*;

fn lj_potential(
    r_squared: units::f64::Area,
    sigma: units::f64::Length,
    epsilon: units::f64::Energy
) -> units::f64::Energy {
    let sigma_over_r_squared = sigma.powi(P2::new()) / r_squared;
    let six = sigma_over_r_squared.powi(P3::new());
    let twelve = six.powi(P2::new());

    4.0 * epsilon * (twelve - six)
}

/// Construct a list of pairs of indices for atoms within cutoff of each other
fn construct_pairlist(
    positions: &[[units::f64::Length; 3]],
    cutoff: units::f64::Length,
    boundaries: &impl BoundaryConditions
) -> Vec<([usize; 2], units::f64::Area)> {
    let cutoff_squared = cutoff * cutoff;
    let mut out = Vec::new();

    for (i, a) in positions.iter().enumerate() {
        for (j, b) in positions[i+1..].iter().enumerate() {
            let r_squared = boundaries.dist2(*a, *b);

            if r_squared < cutoff_squared {
                out.push(([i, j], r_squared))
            };
        }
    }
    out
}

/// Compute the LJ energy of many frames of a homogenous fluid of particles
pub fn lj_from_positions(
    frames: &[units::f64::Positions],
    sigmas: &[units::f64::Length],
    epsilons: &[units::f64::Energy],
    cutoff: units::f64::Length,
    boundaries: &impl BoundaryConditions
) -> Vec<units::f64::Energy> {
    if sigmas.len() != epsilons.len() {
        panic!("lj_from_positions requires equal numbers of sigmas and epsilons")
    }

    let mut energies = Vec::new();
    let cutoff_squared = cutoff.powi(P2::new());

    for frame in frames {
        let mut energy = 0.0 * units::f64::KJPERMOL;

        if frame.len() != sigmas.len() {
            panic!("lj_from_positions requires frames with equal numbers of atoms to sigmas and epsilons")
        }

        for ([i, j], r_squared) in construct_pairlist(frame, cutoff, boundaries) {
            let sigma = (sigmas[i] * sigmas[j]).sqrt();
            let epsilon = (epsilons[i] * epsilons[j]).sqrt();

            energy += lj_potential(r_squared, sigma, epsilon) - lj_potential(cutoff_squared, sigma, epsilon);
        }

        energies.push(energy);
    }

    energies
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
        );

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
        );

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
