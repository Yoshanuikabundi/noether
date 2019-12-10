use super::*;

/// SimplePairlist must recompute the entire pairlist whenever
/// the positions change and stores all pairs within a cutoff
/// in its neighbourlist
pub struct SimplePairlist<B: BoundaryConditions> {
    pairs: Vec<AtomPair>,
    cutoff: f64::Length,
    _marker: std::marker::PhantomData<B>
}

impl<B: BoundaryConditions> SimplePairlist<B> {
    pub fn new(cutoff: f64::Length) -> Self where Self: Sized {
        Self {
            pairs: vec![],
            cutoff,
            _marker: std::marker::PhantomData
        }
    }
}

impl<B: BoundaryConditions> Neighbourlist<B> for SimplePairlist<B> {
    /// Update the pairlist based on the positions of atoms
    ///
    /// This completely rebuilds the pairlist from scratch.
    /// It is an $\mathcal{O}(n^2)$ operation. SimplePairlist
    /// should generally not be used in simulations, but may
    /// be useful when evaluating energies of fixed states.
    fn regenerate(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &B
    ) -> Result<()> {
        self.pairs.clear();

        let cutoff_squared = self.cutoff * self.cutoff;

        for (i, a) in positions.iter().enumerate() {
            for (j, b) in (i+1..).zip(positions[i+1..].iter()) {
                let r_squared = boundaries.dist2(*a, *b);

                if r_squared < cutoff_squared {
                    self.pairs.push(([i, j], r_squared))
                };
            }
        }

        Ok(())
    }

    /// Update the pairlist based on the positions of atoms
    ///
    /// This updates the positions of atoms based on the
    /// existing pairlist. It assumes no atoms have crossed
    /// the cutoff boundary and is probably not useful at all.
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &B
    ) -> Result<()> {
        for ([i, j], r_squared) in self.pairs.iter_mut() {
            let a = positions[*i];
            let b = positions[*j];
            *r_squared = boundaries.dist2(a, b);
        }

        Ok(())
    }

    /// Iterate over the atom pairs in the pairlist
    fn iter(&self) -> std::slice::Iter<AtomPair> {
        self.pairs.iter()
    }

    fn neighbourlist_params(&self) -> NeighbourlistParams {
        NeighbourlistParams::NonbondedCutoff(self.cutoff)
    }
}


#[cfg(test)]
mod tests {
    use crate::neighbourlist::AtomPair;
    use crate::neighbourlist::simple::SimplePairlist;
    use crate::neighbourlist::Neighbourlist;
    use crate::units::f64;
    use crate::result::Result;

    fn assert_pairs_eq(pairs_a: &[AtomPair], pairs_b: &[AtomPair]) -> Result<()> {
        match pairs_a
            .iter()
            .zip(pairs_b)
            .try_for_each(|(([i_a, j_a], r_a), ([i_b, j_b], r_b))| {
                assert_eq!(i_a, i_b);
                assert_eq!(j_a, j_b);
                assert_eq!(r_a, r_b);
                Ok(())
            })
        {
            Ok(_) => Ok(()),
            Err(e) => Err(e)
        }
    }

    #[test]
    fn check_update_regenerate() {
        let mut pairlist = SimplePairlist::new(1.2 * f64::NM);
        let boundaries = crate::boundaries::NoBounds;
        pairlist.regenerate(
            &vec![
                [0.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
                [1.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
                [0.0 * f64::NM, 1.0 * f64::NM, 0.0 * f64::NM],
                [0.0 * f64::NM, 0.0 * f64::NM, 1.0 * f64::NM],
                [2.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
            ],
            &boundaries
        ).unwrap();

        assert_pairs_eq(&pairlist.pairs, &[
            ([0, 1], 1.00 * f64::NM * f64::NM),
            ([0, 2], 1.00 * f64::NM * f64::NM),
            ([0, 3], 1.00 * f64::NM * f64::NM),
            ([1, 4], 1.00 * f64::NM * f64::NM),
        ]).unwrap();

        pairlist.update(
            &vec![
                [0.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
                [1.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
                [0.0 * f64::NM, 1.5 * f64::NM, 0.0 * f64::NM],
                [0.0 * f64::NM, 0.0 * f64::NM, 1.0 * f64::NM],
                [1.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
            ],
            &boundaries
        ).unwrap();

        assert_pairs_eq(&pairlist.pairs, &[
            ([0, 1], 1.00 * f64::NM * f64::NM),
            ([0, 2], 2.25 * f64::NM * f64::NM),
            ([0, 3], 1.00 * f64::NM * f64::NM),
            ([1, 4], 0.00 * f64::NM * f64::NM),
        ]).unwrap();

        pairlist.regenerate(
            &vec![
                [0.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
                [1.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
                [0.0 * f64::NM, 1.5 * f64::NM, 0.0 * f64::NM],
                [0.0 * f64::NM, 0.0 * f64::NM, 1.0 * f64::NM],
                [1.0 * f64::NM, 0.0 * f64::NM, 0.0 * f64::NM],
            ],
            &boundaries
        ).unwrap();

        assert_pairs_eq(&pairlist.pairs, &[
            ([0, 1], 1.00 * f64::NM * f64::NM),
            ([0, 3], 1.00 * f64::NM * f64::NM),
            ([0, 4], 1.00 * f64::NM * f64::NM),
            ([1, 4], 0.00 * f64::NM * f64::NM),
        ]).unwrap();
    }
}
