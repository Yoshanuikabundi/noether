extern crate noether;

use noether::units::f64::KJPERMOL;
use noether::units::f64::NM;
use noether::io;
use noether::lj_from_positions;
use noether::boundaries::*;
use noether::FriendlyResult;
use noether::topology::Topology;

fn main() {
    let topol_2atoms = Topology::lj_fluid(
        "Me".to_string(),
        2,
        0.3405 * NM,
        1.0000 * KJPERMOL,
        1.2 * NM,
    ).unwrap();

    let topol_100atoms = Topology::lj_fluid(
        "Me".to_string(),
        100,
        0.3405 * NM,
        1.0000 * KJPERMOL,
        1.2 * NM,
    ).unwrap();

    let positions = io::read_positions("test_targets/2_atoms/2_atoms_frommax.trr").unwrap();

    let energies = lj_from_positions(
        &positions,
        &topol_2atoms,
        &NoBounds,
    ).unwrap_nicely();

    let energies_ref = io::read_xvg("test_targets/2_atoms/2_atoms_frommax.xvg", 1).unwrap();
    for (a, b) in energies.iter().zip(energies_ref.iter()) {
        let b: f64 = b.parse().unwrap();
        println!("noether = {:?}, gromacs = {:?}", *a, b * KJPERMOL);
    }

    let positions = io::read_positions("test_targets/100_atoms/100_atoms.trr").unwrap();

    let energies = lj_from_positions(
        &positions,
        &topol_100atoms,
        &Pbc::cubic(5.0*NM),
    ).unwrap_nicely();

    let energies_ref = io::read_xvg("test_targets/100_atoms/100_atoms.xvg", 1).unwrap();
    for (a, b) in energies.iter().zip(energies_ref.iter()) {
        let b: f64 = b.parse().unwrap();
        println!("noether = {:?}, gromacs = {:?}", *a, b * KJPERMOL);
    }

    lj_from_positions(
        &positions,
        &topol_2atoms,
        &NoBounds,
    ).unwrap_nicely();
}
