#[macro_use]
extern crate uom;

extern crate typenum;

pub mod geom;
pub mod units;

mod potentials {
    mod bonded {
    }

    mod nonbonded {
    }
}

mod samplers {
    mod mc {
        // Monte Carlo sampler
    }

    mod ld {
        // Langevin Dynamics sampler
    }

    mod point_energy{
    }
}

pub mod topology {
    use units::{
        f32,
        // f64,
        // vec3
    };
    use units::charge::{
        elem_charge,
    };

    #[derive(Debug)]
    pub struct Top {
        atoms: Vec<Atom>,
    }

    impl Top {
        pub fn gen_lj_fluid(num:usize, mass: f32::Mass, epsilon: f32::Energy, sigma: f32::Length) -> Top {
            let atom = Atom {
                mass,
                epsilon,
                sigma,
                charge: f32::Charge::new::<elem_charge>(0.0),
            };
            let atoms = vec![atom.clone(); num];
            Top { atoms }
        }
    }

    #[derive(Debug, Clone)]
    struct Atom {
        mass: f32::Mass,
        charge: f32::Charge,
        epsilon: f32::Energy,
        sigma: f32::Length,
    }
}




#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
