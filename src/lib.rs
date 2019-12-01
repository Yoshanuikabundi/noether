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
    epsilon: units::f64::Energy,
    shift: units::f64::Energy,
) -> units::f64::Energy {
    let sigma_over_r_squared = sigma.powi(P2::new()) / r_squared;
    let six = sigma_over_r_squared.powi(P3::new());
    let twelve = six.powi(P2::new());

    4.0 * epsilon * (twelve - six) - shift
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
    sigma: units::f64::Length,
    epsilon: units::f64::Energy,
    cutoff: units::f64::Length,
    boundaries: &impl BoundaryConditions
) -> Vec<units::f64::Energy> {
    let mut energies = Vec::new();
    let shift = lj_potential(cutoff.powi(P2::new()), sigma, epsilon, 0.0 * units::f64::KJPERMOL);

    for frame in frames {
        let mut energy = 0.0 * units::f64::KJPERMOL;

        for ([_i, _j], r_squared) in construct_pairlist(frame, cutoff, boundaries) {
            energy += lj_potential(r_squared, sigma, epsilon, shift);
        }

        energies.push(energy);
    }

    energies
}

// #[cfg(test)]
// mod tests {
//     use crate::units::f64::KJPERMOL;
//     use crate::units::f64::NANOMETER;
//     use super::io;

//     #[test]
//     fn lj_potential() {
//         let positions = io::read_positions("test_targets/equil.trr").unwrap();
//         let energies = super::lj_from_positions(positions, 0.006_221_27 * KJPERMOL, 9.695_76_E-6 * NANOMETER);

//         let energies_ref = io::read_xvg("test_targets/energy.xvg", 1);

//         for (a, b) in energies.iter().zip(energies_ref) {
//             assert_lt!((*a - *b).abs(), 1.0E-4 * KJPERMOL);
//         }
//     }
// }
