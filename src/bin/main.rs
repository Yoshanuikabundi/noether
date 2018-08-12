extern crate noether;

use noether::topology::Top;
use noether::state::State;
use noether::units::f32consts::*;
use noether::geom::*;

extern crate rand;
extern crate chemfiles;
use chemfiles::{Trajectory, Frame, Atom, UnitCell};

use rand::Rng;

fn main() {
    println!("OK let's go");

    const SYS_SIZE: usize = 10_000;

    let top = Top::gen_lj_fluid(SYS_SIZE, 12.0*DA, 0.3*KJPM, 0.30*NM);

    let mut rng = rand::thread_rng();

    let l = 30.0_f32;
    let boxvecs = (
        PosVec::from(  l, 0.0, 0.0),
        PosVec::from(0.0,   l, 0.0),
        PosVec::from(0.0, 0.0,   l)
    );

    let state = State::new(
        &top,
        (0..SYS_SIZE)
            .map(|_| PosVec::from(
                rng.gen_range(0.0, l * 0.5),
                rng.gen_range(0.0, l * 0.5),
                rng.gen_range(0.0, l * 0.5)
            )).collect(),
        (0..SYS_SIZE)
            .map(|_| VelocVec::from(
                0.0,
                0.0,
                0.0
            )).collect(),
        boxvecs,
        "trajout.pdb".to_string()
    );

    // println!("State generated with {} LJ spheres, calculating energy...", SYS_SIZE);

    // let energy = state.calc_energy();

    let nsteps = 10_000_000;
    let temp = 300.0 * K;

    write_state(String::from("start.pdb"), &state, l);

    println!("Monte carlo time! Let's do {} steps at {}.", nsteps, temp);

    let final_state = state.sample(nsteps, temp);

    let energy = final_state.calc_energy();

    let out = String::from("finish.pdb");
    println!("Found a state with energy {}! Writing to {}. Bye bye!", energy, out);

    write_state(out, &final_state, l);


}

fn write_state(filename:String, state:&State, box_l:f32) {
    let mut frame = Frame::new().unwrap();
    for posvec in state.positions.iter() {
        let (x, y, z) = (
            posvec.x.value_unsafe as f64 * 10.0,
            posvec.y.value_unsafe as f64 * 10.0,
            posvec.z.value_unsafe as f64 * 10.0
        );
        frame.add_atom(&Atom::new("LJ").unwrap(), [x, y, z], None).unwrap();
    }
    let unit_cell = UnitCell::new([box_l as f64 * 10.0; 3]).unwrap();
    frame.set_cell(&unit_cell).unwrap();

    let mut trajectory = Trajectory::open(filename, 'w').unwrap();
    trajectory.write(&frame).unwrap();
}
