use crate::boundaries::NoBounds;
use super::*;

/// Verlet pairlist
#[derive(Debug)]
pub struct VerletPairlist<B: BoundaryConditions> {
    pairs: Vec<AtomPair>,
    cutoff: f64::Length,
    buffer: f64::Length,
    _marker: std::marker::PhantomData<B>,
    cell_list: CellList<B>

}

impl VerletPairlist<NoBounds> {
    pub fn with_buffer(
        cutoff: f64::Length,
        buffer: f64::Length,
        n_cells: usize
    ) -> Self
    {
        Self {
            pairs: Vec::new(),
            cutoff,
            buffer,
            _marker: std::marker::PhantomData,
            cell_list: CellList::new(
                NoBounds,
                [n_cells; 3],
                [cutoff + buffer; 3],
                &[]
            )

        }
    }
}

impl Neighbourlist<NoBounds> for VerletPairlist<NoBounds> {
    /// Regenerate the neighbourlist based on the positions of atoms
    ///
    /// Regenerates the cell list, and then rewrites the pairlist.
    /// This is a terrible API and not an efficient thing to do
    fn regenerate(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &NoBounds
    ) -> Result<()> {
        self.cell_list.regenerate(positions);
        self.pairs = self.cell_list.iter_pairs().map(|(i, j)| {
            let a = positions[i];
            let b = positions[j];
            let r_squared = boundaries.dist2(a, b);
            ([i, j], r_squared)
        }).collect();
        Ok(())
    }

    /// Update the neighbourlist based on the positions of atoms
    ///
    /// Takes the existing pairs and regenerates the positions.
    /// This is also not very useful.
    fn update(
        &mut self,
        positions: &[[f64::Length; 3]],
        boundaries: &NoBounds,
    ) -> Result<()> {
        for ([i, j], r_squared) in self.pairs.iter_mut() {
            *r_squared = boundaries.dist2(positions[*i], positions[*j]);
        }
        Ok(())
    }

    /// Iterate over the atom pairs in the neighbourlist
    fn iter(&self) -> std::slice::Iter<AtomPair> {
        unimplemented!()
    }

    /// Get the parameters of the neighbourlist
    fn neighbourlist_params(&self) -> NeighbourlistParams {
        unimplemented!()
    }
}

/// Stores atom indices according to their positions
/// in the unit cell
#[derive(Debug)]
struct CellList<B: BoundaryConditions> {
    /// Number of cells per unit cell in the x, y, z directions
    n_cells: [usize; 3],
    /// The length of each cell in the x, y, z directions
    cell_length: [f64::Length; 3],
    /// The boundary conditions
    boundaries: B,
    /// The indices of atoms in each cell
    ///
    /// Cells are stored as (0,0,0), (0,0,1), (0,0,2), (0,1,0) etc...
    /// Indices of the start of each cell are stored in the head field
    indices: Vec<Option<usize>>,
    /// The indices of the first atoms of each cell
    head: Vec<Option<usize>>,
}

impl CellList<NoBounds> {
    fn new(
        boundaries: NoBounds,
        n_cells: [usize; 3],
        cell_length: [f64::Length; 3],
        positions: &[[f64::Length; 3]]
    ) -> Self {
        let [a, b, c] = n_cells;

        let mut new = Self {
            n_cells,
            cell_length,
            boundaries,
            indices: vec![None; positions.len()],
            head: vec![None; a * b * c],
        };

        new.regenerate(positions);

        new
    }

    fn regenerate(&mut self, positions: &[[f64::Length; 3]]) {
        for (i, [x, y, z]) in positions.iter().enumerate() {
            let cell_index = self.position_to_cell(&[*x, *y, *z]);
            self.add_atom(i, cell_index);
        }
    }

    #[allow(clippy::many_single_char_names)]
    fn position_to_cell(&self, position: &[f64::Length; 3]) -> usize {
        let [x, y, z] = position;
        let [i, j, k] = self.n_cells;
        let [l_x, l_y, l_z] = self.cell_length;

        // Cells on the edge of the box are open to the rest of infinity
        let x_idx = std::cmp::min(f64::from(*x / l_x) as usize, i - 1);
        let y_idx = std::cmp::min(f64::from(*y / l_y) as usize, j - 1);
        let z_idx = std::cmp::min(f64::from(*z / l_z) as usize, k - 1);

        z_idx + k * y_idx + j * k * x_idx
    }

    fn iter_neighbouring_cells(&self, cell: usize, neighbours: usize) -> impl Iterator<Item=usize> {
        let [i, j, k] = self.n_cells;
        let x_idx = cell / (j * k);
        let y_idx = (cell - x_idx) / k;
        let z_idx = cell - x_idx - y_idx;

        let neighbours = neighbours as isize;

        (-neighbours..neighbours).flat_map(move |dx| {
            (-neighbours..neighbours).flat_map(move |dy| {
                (-neighbours..neighbours).map(move |dz| {
                    let z_idx_d: usize = match z_idx as isize + dz {
                        z if z >= k as isize => k - 1,
                        z if z < 0 => 0,
                        z => z as usize
                    };
                    let y_idx_d: usize = match y_idx as isize + dy {
                        y if y >= j as isize => j - 1,
                        y if y < 0 => 0,
                        y => y as usize
                    };
                    let x_idx_d: usize = match x_idx as isize + dx {
                        x if x >= i as isize => i - 1,
                        x if x < 0 => 0,
                        x => x as usize
                    };
                    z_idx_d + k * y_idx_d + j * k * x_idx_d
                })
            })
        })
    }

    fn iter_pairs(&self) -> impl Iterator<Item=(usize, usize)> + '_ {
        self.iter_cells().flat_map(move |cell_index| {
            self.iter_cell_atoms(cell_index).flat_map(move |atom_a_idx| {
                self.iter_neighbouring_cells(cell_index, 1).flat_map(move |neighbour_cell_index| {
                    self.iter_cell_atoms(neighbour_cell_index).map(move |atom_b_idx| {
                        (atom_a_idx, atom_b_idx)
                    })
                })
            })
        })
    }
}

impl<B: BoundaryConditions> CellList<B> {

    /// Add `atom_index` to the end of cell `cell_index`
    ///
    /// # Panics
    ///
    /// Panics if `atom_index` or `cell_index` are invalid
    fn add_atom(&mut self, atom_index: usize, cell_index: usize) {
        match self.head[cell_index] {
            None => self.head[cell_index] = Some(atom_index),
            Some(mut i) => {
                i = loop {
                    i = match self.indices[i] {
                        Some(i) => i,
                        None => break i
                    }
                };
                self.indices[i] = Some(atom_index);
                self.indices[atom_index] = None;
            }
        }
    }

    /// Remove `atom_index`, currently in `cell_index`, from the cell list
    ///
    /// # Panics
    ///
    /// Panics if `atom_index` or `cell_index` are invalid
    fn remove_atom(&mut self, atom_index: usize, cell_index: usize) {
        let next_i = self.indices[atom_index];
        match self.head[cell_index] {
            None => panic!("Could not find atom in CellList at old_cell"),
            Some(i) if i == atom_index => self.head[cell_index] = next_i,
            Some(mut prev_i) => {
                while let Some(i) = self.indices[prev_i] {
                    if i == atom_index {
                        self.indices[prev_i] = next_i;
                        self.indices[i] = None;
                        return;
                    }
                    prev_i = i;
                }
                panic!("Could not find atom in CellList at old_cell");
            }
        };
    }

    /// Move `atom_index` from `old_cell` to `new_cell`
    ///
    /// # Panics
    ///
    /// Panics if `atom_index`, `old_cell` or `new_cell` are invalid
    fn move_atom(&mut self, atom_index: usize, old_cell: usize, new_cell: usize) {
        if old_cell == new_cell { return };
        self.remove_atom(atom_index, old_cell);
        self.add_atom(atom_index, new_cell);
    }

    fn iter_cell_atoms(&self, cell: usize) -> CellIter<B> {
        CellIter {
            i: CellIterStage::Head(cell),
            cells: &self
        }
    }

    fn iter_cells(&self) -> impl Iterator<Item=usize> {
        let [i, j, k] = self.n_cells;

        (0..i).flat_map(move |x_idx| {
            (0..j).flat_map(move |y_idx| {
                (0..j).map(move |z_idx| z_idx + k * y_idx + j * k * x_idx)
            })
        })
    }
}

enum CellIterStage {
    Head(usize),
    Indices(usize),
    None
}

struct CellIter<'a, B: BoundaryConditions> {
    i: CellIterStage,
    cells: &'a CellList<B>
}

impl<B: BoundaryConditions> Iterator for CellIter<'_, B> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        match self.i {
            CellIterStage::Head(i) => match self.cells.head[i] {
                Some(i) => {
                    self.i = CellIterStage::Indices(i);
                    Some(i)
                },
                None =>  {
                    self.i = CellIterStage::None;
                    None
                }
            },
            CellIterStage::Indices(i) => match self.cells.indices[i] {
                Some(i) => {
                    self.i = CellIterStage::Indices(i);
                    Some(i)
                },
                None =>  {
                    self.i = CellIterStage::None;
                    None
                }
            },
            CellIterStage::None => None
        }
    }

}
