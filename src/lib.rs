#[macro_use]
extern crate uom;

#[cfg(test)]
#[macro_use]
extern crate more_asserts;

pub mod units;
pub mod io;

/// Compute the LJ energy of many frames of a homogenous fluid of particles
pub fn lj_from_positions(
    positions: Vec<units::f64::Positions>,
    sigma: units::f64::Length,
    epsilon: units::f64::Energy,
    cutoff: units::f64::Length,
) -> Vec<units::f64::Energy> {
    let mut energies = Vec::new();

    let two = uom::typenum::P2::new();
    let three = uom::typenum::P3::new();

    let r_max_squared = cutoff * cutoff;

    for frame in positions {
        let mut energy = 0.0 * units::f64::KJPERMOL;
        for (i, [x1, y1, z1]) in frame.iter().enumerate() {
            for (j, [x2, y2, z2]) in frame.iter().enumerate() {
                if i == j {continue};

                let r_squared = (*x2-*x1).powi(two) + (*y2-*y1).powi(two) + (*z2-*z1).powi(two);

                if r_squared > r_max_squared {continue};

                let sigma_over_r_squared = sigma.powi(two) / r_squared;
                let six = sigma_over_r_squared.powi(three);
                let twelve = six.powi(two);
                energy += 4.0 * epsilon * (twelve - six);
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
