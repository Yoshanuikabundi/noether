#[macro_use]
extern crate dimensioned as dim;

extern crate rand;
extern crate itertools;
extern crate chemfiles;
extern crate rayon;

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

pub mod state {
    use units::*;
    use units::f32consts::*;

    use geom::{
        PosVec,
        VelocVec
    };
    use topology::Top;
    use rand;
    use rand::Rng;
    use itertools::Itertools;
    use rayon::prelude::*;
    use chemfiles;
    use chemfiles::{Trajectory, Frame, Atom, UnitCell};

    pub struct State<'a> {
        pub topology: &'a Top,
        pub positions: Vec<PosVec>,
        pub velocities: Vec<VelocVec>,
        pairlist: Vec<(usize, usize)>,
        boxvecs: (PosVec, PosVec, PosVec),
        trajout: String
    }

    impl<'a> State<'a> {
        pub fn new(
            topology:&Top,
            positions: Vec<PosVec>,
            velocities: Vec<VelocVec>,
            boxvecs: (PosVec, PosVec, PosVec),
            filename: String
        ) -> State {

            {
                Trajectory::open(&filename, 'w').unwrap();
            }

            let mut state = State {
                topology,
                positions,
                velocities,
                boxvecs,
                trajout: filename,
                pairlist: vec![]
            };

            println!("Generating pairlist without cutoff");

            state.gen_pairs(0.0 * NM);

            state
        }

        pub fn dist2(&self, first:&PosVec, second:&PosVec) -> Nanometer2<f32> {
            let (vec_a, vec_b, vec_c) = self.boxvecs.clone();
            let mut offset: PosVec;
            let mut r2s: Vec<Nanometer2<f32>> = vec![];

            for a in [-1.0, 0.0, 1.0].iter() {
                for b in [-1.0, 0.0, 1.0].iter() {
                    for c in [-1.0, 0.0, 1.0].iter() {
                        offset = vec_a.clone() * a.clone()
                            + vec_b.clone() * b.clone()
                            + vec_c.clone() * c.clone();
                        r2s.push((first - second + offset).norm2());
                    }
                }
            }
            r2s.into_iter()
                .min_by(|a, b| a.partial_cmp(b).expect("Tried to compare a NaN"))
                .unwrap()
        }

        /// Generate a verlet pairlist
        pub fn gen_pairs(&mut self, cutoff:Nanometer<f32>) {
            let cutoff2 = cutoff * cutoff;

            println!("Collecting combinations");
            let pair_vec:Vec<((_, _),(_, _))> = self.positions.iter()
                .enumerate()
                .tuple_combinations()
                .collect();
            println!("Generating pairlist");
            self.pairlist = pair_vec.par_iter()
                .filter(|((_, ri), (_, rj))| cutoff == 0.0 * NM || self.dist2(&ri, &rj) <= cutoff2)
                .map(|((i, _), (j, _))| (i.clone(), j.clone()))
                .collect();
            println!("New pairlist has {} entries", self.pairlist.len());
        }

        pub fn calc_energy(&self) -> KilojoulePerMole<f32> {
            self.topology.calc_energy(
                &self.positions,
                &self.pairlist,
                |ri, rj| self.dist2(&ri, &rj)
            )
        }

        pub fn write_traj(&self) -> chemfiles::Result<()> {
            let mut frame = Frame::new()?;

            for posvec in self.positions.iter() {
                let (x, y, z) = (
                    posvec.x.value_unsafe as f64 * 10.0,
                    posvec.y.value_unsafe as f64 * 10.0,
                    posvec.z.value_unsafe as f64 * 10.0
                );
                let atom = Atom::new("LJ")?;
                frame.add_atom(&atom, [x, y, z], None)?;
            }

            // TODO: Set trajectory unit cell rationally
            let box_l = self.boxvecs.0.x.value_unsafe;

            let unit_cell = UnitCell::new([box_l as f64 * 10.0; 3])?;
            frame.set_cell(&unit_cell)?;

            let mut trajout = Trajectory::open(&self.trajout, 'a')?;
            trajout.write(&frame)?;
            Ok(())
        }

        pub fn sample(mut self, nsteps: usize, temp: Kelvin<f32>) -> Self {
            let mut rng = rand::thread_rng();

            let move_std_dev = 0.001f32;
            let distrib = rand::distributions::Normal::new(0.0, move_std_dev as f64);

            let steps_between_pairlist_updates = 1000;

            let pairlist_cutoff = self.topology.lj_cutoff + (steps_between_pairlist_updates/5) as f32 *  move_std_dev * NM;

            println!("Generating pairlist with cutoff {}", pairlist_cutoff);
            self.gen_pairs(pairlist_cutoff);
            println!("Calculating initial energy");
            let mut prev_energy = self.calc_energy();
            let mut new_energy;
            let mut accepts_since_pairlist_regen = 0;

            for n in 0..nsteps {
                if accepts_since_pairlist_regen % steps_between_pairlist_updates == 0 && accepts_since_pairlist_regen != 0 {
                    println!("Regenerating pairlist at step {}", n);
                    self.gen_pairs(pairlist_cutoff);
                    new_energy = self.calc_energy();
                    if new_energy != prev_energy {
                        println!("New pairlist changed energies: {:?}, {:?}", prev_energy, new_energy);
                    }
                    prev_energy = new_energy;
                    accepts_since_pairlist_regen = 0;
                }

                if n % 5000 == 0 {
                    print!("Step {}, energy is {}, ", n, prev_energy);
                    match self.write_traj() {
                        Ok(()) => println!("frame written to file {}", &self.trajout),
                        Err(e) => println!("frame could not be written to file: {}", e)
                    }
                }

                let mut attempt_pos = self.positions.clone();
                for mut pos in attempt_pos.iter_mut() {
                    *pos += PosVec::from(
                        rng.sample(distrib) as f32,
                        rng.sample(distrib) as f32,
                        rng.sample(distrib) as f32
                    );
                }

                let attempt_energy = self.topology.calc_energy(
                    &attempt_pos,
                    &self.pairlist,
                    |ri, rj| self.dist2(ri, rj)
                );
                let energy_diff = prev_energy - attempt_energy;
                let accept_prob = (energy_diff / (KB * temp)).exp();

                if accept_prob >= rng.gen() {
                    // println!("Accepted move with delta {}, P {:.1}.", -energy_diff, accept_prob);
                    self.positions = attempt_pos;
                    prev_energy = self.calc_energy();
                    accepts_since_pairlist_regen += 1;
                } else {
                    // println!("Rejected move with delta {}, P {:.1}.", -energy_diff, accept_prob);
                }
            }
            self
        }
    }
}

pub mod topology {
    use units::*;
    use units::f32consts::*;
    use geom::PosVec;
    use rayon::prelude::*;
    use std;

    #[derive(Debug)]
    pub struct Top {
        pub atoms: Vec<Atom>,
        pub lj_cutoff: Nanometer<f32>
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
                charge: 0.0 * E
            };
            let atoms = vec![atom.clone(); num];
            Top {
                atoms,
                lj_cutoff: 1.0 * NM
            }
        }

        pub fn calc_energy<F>(&self, positions: &Vec<PosVec>, pairlist: &Vec<(usize, usize)>, dist2: F) -> KilojoulePerMole<f32>
            where
                F: Fn(&PosVec, &PosVec) -> Nanometer2<f32> + std::marker::Sync
        {
            let atoms = &self.atoms;
            let lj_cutoff_squared = self.lj_cutoff * self.lj_cutoff;

            pairlist
                .par_iter()
                .map(|(i, j)| {
                    let r2 = dist2(&positions[i.clone()], &positions[j.clone()]);
                    (i, j, r2)
                }).filter(|(_, _, r2)| r2 <= &lj_cutoff_squared)
                .map(|(i, j, r2)| {
                    // TODO: Allow other LJ combination rules than averaging
                    let eps = (atoms[i.clone()].epsilon + atoms[j.clone()].epsilon) / 2.0;
                    let sig = (atoms[i.clone()].sigma + atoms[j.clone()].sigma) / 2.0;

                    let sig12 = sig.value_unsafe.powi(12);
                    let r12 = r2.value_unsafe.powi(6);
                    let sig6 = sig.value_unsafe.powi(6);
                    let r6 = r2.value_unsafe.powi(3);

                    4.0*ONE * eps * ((sig12 / r12) - (sig6 / r6))
                }).reduce(
                    || 0.0 * KJPM,
                    |acc, e| acc + e
                )

        }
    }

    #[derive(Debug, Clone)]
    pub struct Atom {
        pub mass: Dalton<f32>,
        pub charge: ElemCharge<f32>,
        pub epsilon: KilojoulePerMole<f32>,
        pub sigma: Nanometer<f32>,
    }
}




#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
