use super::*;

/// SimplePairlist must recompute the entire pairlist whenever
/// the positions change
pub struct SimplePairlist<B: BoundaryConditions> {
    pairs: Vec<AtomPair>,
    cutoff: f64::Length,
    _marker: std::marker::PhantomData<B>
}

impl<B: BoundaryConditions> Pairlist<B> for SimplePairlist<B> {
    /// Update the pairlist based on the positions of atoms
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &B
    ) -> Result<()> {
        self.pairs.clear();

        let cutoff_squared = self.cutoff * self.cutoff;

        for (i, a) in positions.iter().enumerate() {
            for (j, b) in positions[i+1..].iter().enumerate() {
                let r_squared = boundaries.dist2(*a, *b);

                if r_squared < cutoff_squared {
                    self.pairs.push(([i, j], r_squared))
                };
            }
        }

        Ok(())
    }

    /// Iterate over the atom pairs in the pairlist
    fn iter(&self) -> std::slice::Iter<AtomPair> {
        self.pairs.iter()
    }

    fn new(cutoff: Cutoff) -> Result<Self> where Self: Sized {
        if let Some(cutoff) = cutoff {
            Ok(Self {
                pairs: vec![],
                cutoff,
                _marker: std::marker::PhantomData
            })
        } else {
            Err(CutoffRequired)
        }
    }

    fn cutoff(&self) -> Cutoff {
        Some(self.cutoff)
    }
}
