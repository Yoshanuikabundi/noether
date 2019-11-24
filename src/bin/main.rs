extern crate noether;

use noether::units::f64::KJPERMOL;
use noether::units::f64::NANOMETER;
use noether::io;
use noether::lj_from_positions;

fn main() {
    let positions = io::read_positions("test_targets/equil.trr").unwrap();

    let energies = lj_from_positions(
        positions,
        0.997_967_163 * KJPERMOL,
        0.340_500_014 * NANOMETER,
        1.2 * NANOMETER
    );

    // let _energies_ref = io::read_xvg("test_targets/energy.xvg", 1);

    // for (a, b) in energies.iter().zip(energies_ref.iter()) {
    //     println!("{:?} == {:?} ?", a, b);
    // }

    for a in energies.iter() {
        println!("{:?}", a);
    }

}
