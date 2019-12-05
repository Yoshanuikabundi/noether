use crate::pairlist::Pairlist;
use crate::topology::Topology;
use crate::boundaries::BoundaryConditions;

pub struct State<'a, B: BoundaryConditions> {
    topology: Topology<'a>,
    boundary_conditions: B,
    pairlists: Vec<Box<dyn Pairlist<B> + 'a>>,
}
