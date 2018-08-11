#[macro_use]
extern crate dimensioned as dim;

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
    use units::*;

    #[derive(Debug)]
    pub struct Top {
        atoms: Vec<Atom>,
    }

    impl Top {
        pub fn gen_lj_fluid(
            num:usize,
            mass: Dalton<f32>,
            epsilon: KilojoulePerMole<f32>,
            sigma: Nanometer<f32>
        ) -> Top {
            let atom = Atom {
                mass,
                epsilon,
                sigma,
                charge: 0.0 * f32consts::E,
            };
            let atoms = vec![atom.clone(); num];
            Top { atoms }
        }
    }

    #[derive(Debug, Clone)]
    struct Atom {
        mass: Dalton<f32>,
        charge: ElemCharge<f32>,
        epsilon: KilojoulePerMole<f32>,
        sigma: Nanometer<f32>,
    }
}




#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
