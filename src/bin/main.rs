extern crate noether;

use noether::units::f64::KJPERMOL;
use noether::units::f64::NM;
use noether::io;
use noether::lj_from_positions;
use noether::boundaries::*;

fn main() {
    let positions = io::read_positions("test_targets/2_atoms/2_atoms_frommax.trr").unwrap();

    lj_from_positions(
        &positions,
        &[0.3405 * NM; 1],
        &[1.0000 * KJPERMOL; 2],
        1.2 * NM,
        &NoBounds,
    ).err().unwrap().explain();

    let energies = lj_from_positions(
        &positions,
        &[0.3405 * NM; 2],
        &[1.0000 * KJPERMOL; 2],
        1.2 * NM,
        &NoBounds,
    ).unwrap();

    let energies_ref = io::read_xvg("test_targets/2_atoms/2_atoms_frommax.xvg", 1).unwrap();
    for (a, b) in energies.iter().zip(energies_ref.iter()) {
        let b: f64 = b.parse().unwrap();
        println!("noether = {:?}, gromacs = {:?}", *a, b * KJPERMOL);
    }


    let positions = io::read_positions("test_targets/100_atoms/100_atoms.trr").unwrap();

    let energies = lj_from_positions(
        &positions,
        &[0.3405 * NM; 100],
        &[1.0000 * KJPERMOL; 100],
        1.2 * NM,
        &Pbc::cubic(5.0*NM),
    ).unwrap();

    let energies_ref = io::read_xvg("test_targets/100_atoms/100_atoms.xvg", 1).unwrap();
    for (a, b) in energies.iter().zip(energies_ref.iter()) {
        let b: f64 = b.parse().unwrap();
        println!("noether = {:?}, gromacs = {:?}", *a, b * KJPERMOL);
    }


}
