extern crate noether;

use noether::units::f64::KJPERMOL;
use noether::units::f64::NANOMETER;
use noether::io;
use noether::lj_from_positions;

fn main() {
    let positions = io::read_positions("test_targets/equil.trr").unwrap();

    let energies = lj_from_positions(
        positions,
        0.3405 * NANOMETER,
        1.0000 * KJPERMOL,
        1.2 * NANOMETER
    );

    let energies_ref = io::read_xvg("test_targets/equil.xvg", 1).unwrap();
    for (a, b) in energies.iter().zip(energies_ref.iter()) {
        let b: f64 = b.parse().unwrap();
        println!("diff = {:?}", *a - (b * KJPERMOL));
    }

    // for a in energies.iter() {
    //     println!("{:?}", a);
    // }

}
