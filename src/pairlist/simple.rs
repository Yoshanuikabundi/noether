use super::*;

/// SimplePairlist must recompute the entire pairlist whenever
/// the positions change
#[derive(Default)]
pub struct SimplePairlist (Vec<AtomPair>);

impl SimplePairlist {
    pub fn new() -> Self {
        Self (Vec::new())
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self (Vec::with_capacity(capacity))
    }
}

impl Pairlist for SimplePairlist {
    /// Update the pairlist based on the positions of atoms
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        cutoff: Cutoff,
        boundaries: &impl BoundaryConditions
    ) -> Result<()> {
        let Self (pairs) = self;
        pairs.clear();

        let cutoff_squared = match cutoff {
            Cutoff::None => return Ok(()),
            Cutoff::At(cutoff) => cutoff * cutoff
        };

        for (i, a) in positions.iter().enumerate() {
            for (j, b) in positions[i+1..].iter().enumerate() {
                let r_squared = boundaries.dist2(*a, *b);

                if r_squared < cutoff_squared {
                    pairs.push(([i, j], r_squared))
                };
            }
        }

        Ok(())
    }

    /// Iterate over the atom pairs in the pairlist
    fn iter(&self) -> std::slice::Iter<AtomPair> {
        let Self (pairs) = self;
        pairs.iter()
    }
}
