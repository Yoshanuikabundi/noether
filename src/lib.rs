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

        /// Get the minimum image convention distance between two vectors,
        /// and the corresponding distance vector
        pub fn dist2(&self, first:&PosVec, second:&PosVec) -> (PosVec, Nanometer2<f32>) {
            let (vec_a, vec_b, vec_c) = self.boxvecs.clone();
            let mut r2s: Vec<(PosVec, Nanometer2<f32>)> = vec![];

            for a in [-1.0, 0.0, 1.0].iter() {
                for b in [-1.0, 0.0, 1.0].iter() {
                    for c in [-1.0, 0.0, 1.0].iter() {
                        let offset = vec_a.clone() * a.clone()
                            + vec_b.clone() * b.clone()
                            + vec_c.clone() * c.clone();
                        let diff = first - second + offset;
                        let diff2 = diff.clone().norm2();
                        r2s.push((diff, diff2));
                    }
                }
            }
            r2s.into_iter()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).expect("Tried to compare a NaN"))
                .unwrap()
        }

        /// Generate a verlet pairlist
        pub fn gen_pairs(&mut self, cutoff:Nanometer<f32>) {
            let cutoff2 = cutoff * cutoff;

            let pair_vec:Vec<((_, _),(_, _))> = self.positions.iter()
                .enumerate()
                .tuple_combinations()
                .collect();
            self.pairlist = pair_vec.par_iter()
                .filter(|((_, ri), (_, rj))| cutoff == 0.0 * NM || self.dist2(&ri, &rj).1 <= cutoff2)
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

            // TODO: Stop assuming cubic box
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

        pub fn simulate(mut self, nsteps: usize, timestep: Picosecond<f32>) -> Self {
            // let mut rng = rand::thread_rng();

            let steps_between_pairlist_updates = 10;

            let pairlist_cutoff = self.topology.lj_cutoff + 1.0 * NM;

            let buffer_tolerance = 0.005 * KJPM;

            println!("Generating pairlist with cutoff {}", pairlist_cutoff);
            self.gen_pairs(pairlist_cutoff);
            println!("Calculating initial energy");

            let dt = timestep;
            let n_atoms = self.velocities.len();

            assert_eq!(n_atoms, self.positions.len());
            assert_eq!(n_atoms, self.topology.atoms.len());

            for n in 0..nsteps {
                if n % steps_between_pairlist_updates == 0 && n != 0 {
                    println!("Regenerating pairlist at step {}", n);
                    let prev_energy = self.calc_energy();
                    self.gen_pairs(pairlist_cutoff);
                    let new_energy = self.calc_energy();
                    if new_energy != prev_energy {
                        println!("New pairlist changed energies: {:?}, {:?}", prev_energy, new_energy);
                    }
                    if (new_energy - prev_energy).value_unsafe.abs() > buffer_tolerance.value_unsafe {
                        panic!("Energy difference greater than buffer tolerance!");
                    }
                }

                if n % 10 == 0 {
                    let temp:Kelvin<f32> = self.velocities.iter()
                        .zip(&self.topology.atoms)
                        .map(|(v, atom)| {
                            let mass = atom.mass;
                            0.5 * mass * v.norm2()
                        }).fold(
                            0.0 * KJPM,
                            |acc, k| acc + k
                        ) / n_atoms as f32 / KB;

                    print!(
                        "Step {}, potential energy is {}, temperature is {}, ",
                        n,
                        self.calc_energy(),
                        temp
                    );

                    match self.write_traj() {
                        Ok(()) => println!("frame written to file {}", &self.trajout),
                        Err(e) => println!("frame could not be written to file: {}", e)
                    }
                }

                self.positions = self.velocities.iter()
                    .zip(self.positions)
                    .map(|(v, r)| {
                        r + v.clone() * dt/2.0
                    }).collect();

                let forces = self.topology.calc_forces(
                    &self.positions,
                    &self.pairlist,
                    |ri, rj| self.dist2(ri, rj)
                );

                self.velocities = forces.into_iter()
                    .zip(&self.velocities)
                    .zip(&self.topology.atoms)
                    .map(|((f, v), atom)| {
                        let mass = atom.mass;
                        v.clone() + f * dt / mass
                    }).collect();

                self.positions = self.velocities.iter()
                    .zip(self.positions)
                    .map(|(v, r)| {
                        r + v.clone() * dt/2.0
                    }).collect();


                // TODO: Stop assuming rectangular box
                let box_x = self.boxvecs.0.x;
                let box_y = self.boxvecs.1.y;
                let box_z = self.boxvecs.2.z;

                self.positions = self.positions.into_iter()
                    .map(|mut pos| {
                        while pos.x < 0.0 * NM {pos.x += box_x};
                        while pos.y < 0.0 * NM {pos.y += box_y};
                        while pos.z < 0.0 * NM {pos.z += box_z};
                        pos.x %= box_x;
                        pos.y %= box_y;
                        pos.z %= box_z;
                        pos
                    }).collect();

            }
            self
        }
    }
}

pub mod topology {
    use units::*;
    use units::f32consts::*;
    use geom::{
        PosVec,
        ForceVec
    };
    use rayon::prelude::*;
    use std;
    use dim::Sqrt;

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
                F: Fn(&PosVec, &PosVec) -> (PosVec, Nanometer2<f32>) + std::marker::Sync
        {
            let atoms = &self.atoms;
            let lj_cutoff_squared = self.lj_cutoff * self.lj_cutoff;

            pairlist
                .par_iter()
                .map(|(i, j)| {
                    let (_, r2) = dist2(&positions[*i], &positions[*j]);
                    (i, j, r2)
                }).filter(|(_, _, r2)| r2 <= &lj_cutoff_squared)
                .map(|(i, j, r2)| {
                    // TODO: Allow other LJ combination rules than averaging
                    let eps = (atoms[*i].epsilon + atoms[*j].epsilon) / 2.0;
                    let sig = (atoms[*i].sigma + atoms[*j].sigma) / 2.0;

                    let sig6 = sig.value_unsafe.powi(6);
                    let r6 = r2.value_unsafe.powi(3);
                    let sig12 = sig6 * sig6;
                    let r12 = r6 * r6;

                    4.0*ONE * eps * ((sig12 / r12) - (sig6 / r6))
                }).reduce(
                    || 0.0 * KJPM,
                    |acc, e| acc + e
                )

        }

        pub fn calc_forces<F>(&self, positions: &Vec<PosVec>, pairlist: &Vec<(usize, usize)>, dist2: F) -> Vec<ForceVec>
            where
                F: Fn(&PosVec, &PosVec) -> (PosVec, Nanometer2<f32>) + std::marker::Sync
        {
            let atoms = &self.atoms;
            let lj_cutoff_squared = self.lj_cutoff * self.lj_cutoff;
            let mut forces:Vec<ForceVec> = vec![ForceVec::zero(); atoms.len()];


            for (i, j) in pairlist.iter() {
                let (r, r2) = dist2(&positions[*i], &positions[*j]);
                if r2 <= lj_cutoff_squared {
                    // TODO: Allow other LJ combination rules than averaging
                    let eps = (atoms[*i].epsilon + atoms[*j].epsilon) / 2.0;
                    let sig = (atoms[*i].sigma + atoms[*j].sigma) / 2.0;

                    let sig6 = sig.value_unsafe.powi(6);
                    let r6 = r2.value_unsafe.powi(3);
                    let sig12 = sig6 * sig6;
                    let r12 = r6 * r6;

                    let f = r.normalize_into()/NM * 48.0*ONE * (eps / r2.sqrt()) * ((sig12 / r12) - (sig6 / r6));

                    // println!("{:?}", f);
                    forces[*i] += f.clone();
                    forces[*j] -= f;
                }
            }

            forces
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
