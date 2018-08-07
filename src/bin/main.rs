extern crate noether;

use noether::topology::Top;
use noether::units::f32;
use noether::units::f64;

fn main() {
    println!("OK let's go");

    println!("{:?}", f64::NM * f64::NM);

    let top = Top::gen_lj_fluid(1000000usize, 12.0*f32::DA, 1.0*f32::KJPM, 0.10*f32::NM);
}
