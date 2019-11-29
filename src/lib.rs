#[macro_use]
extern crate uom;

#[cfg(test)]
#[macro_use]
extern crate more_asserts;

pub mod units;
pub mod io;

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

/// Compute the LJ energy of many frames of a homogenous fluid of particles
pub fn lj_from_positions(
    positions: Vec<units::f64::Positions>,
    sigma: units::f64::Length,
    epsilon: units::f64::Energy,
    cutoff: units::f64::Length,
) -> Vec<units::f64::Energy> {
    let mut energies = Vec::new();

    let r_max_squared = cutoff * cutoff;
    let shift = lj_potential(cutoff.powi(P2::new()), sigma, epsilon, 0.0 * units::f64::KJPERMOL);

    for frame in positions {
        let mut energy = 0.0 * units::f64::KJPERMOL;
        for (i, [x1, y1, z1]) in frame.iter().enumerate() {
            for (j, [x2, y2, z2]) in frame.iter().enumerate() {
                if i == j {continue};

                let r_squared = (*x2-*x1).powi(P2::new()) + (*y2-*y1).powi(P2::new()) + (*z2-*z1).powi(P2::new());

                if r_squared < r_max_squared {
                    energy += lj_potential(r_squared, sigma, epsilon, shift);
                };
            }
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
