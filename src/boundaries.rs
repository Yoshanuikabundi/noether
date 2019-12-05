use crate::pairlist::{Pairlist};
use crate::units::f64;
use crate::result::*;

use uom::typenum::consts::P2;

/// Boundary conditions of a simulation
///
/// Implements distance computations and probably more in the future
pub trait BoundaryConditions {
    /// Compute the distance between two points given the boundary conditions
    ///
    /// # Arguments
    ///
    /// * `a` - the first point
    /// * `b` - the second point
    ///
    /// # Examples
    ///
    /// The default implementation is just the square root of the squared distance.
    ///
    /// ```
    /// # use noether::units::f64;
    /// # let a = [1.2*f64::NM, 1.8*f64::NM, 1.1*f64::NM];
    /// # let b = [4.8*f64::NM, 1.2*f64::NM, 3.3*f64::NM];
    /// # let boundary_conditions = noether::boundaries::NoBounds;
    /// # use noether::boundaries::BoundaryConditions;
    /// #
    /// assert_eq!(boundary_conditions.dist2(a, b).sqrt(), boundary_conditions.dist(a, b))
    /// ```
    fn dist(
        &self,
        a: [f64::Length; 3],
        b: [f64::Length; 3],
    ) -> f64::Length {
        self.dist2(a, b).sqrt()
    }

    /// Compute the squared distance between two points given the boundary conditions
    ///
    /// # Arguments
    ///
    /// * `a` - the first point
    /// * `b` - the second point
    fn dist2(
        &self,
        a: [f64::Length; 3],
        b: [f64::Length; 3],
    ) -> f64::Area {
        let [x1, y1, z1] = a;
        let [x2, y2, z2] = b;

        (x2 - x1).powi(P2::new())
         + (y2 - y1).powi(P2::new())
         + (z2 - z1).powi(P2::new())
    }

    fn pairlist_checks<P: Pairlist<Self>>(
        &self,
        _pairlist: &P
    ) -> Result<()>  where Self: std::marker::Sized {
        Ok(())
    }
}

/// No boundary conditions
///
/// Particles are in an infinite box
pub struct NoBounds;

impl BoundaryConditions for NoBounds {}

/// Periodic boundary conditions
///
/// Particles are in a triclinic box that neighbours itself ad infinitum.
/// Distances are calculated with the minimum image convention. Boxes are
/// stored as box vectors.
pub struct Pbc(
    pub [f64::Length; 3],
    pub [f64::Length; 3],
    pub [f64::Length; 3]
);

impl Pbc {

    /// Create a cubic box with side length $d$
    ///
    /// The volume is $d^3$
    pub fn cubic(d: f64::Length) -> Self {
        Pbc(
            [          d, 0.0*f64::NM, 0.0*f64::NM],
            [0.0*f64::NM,           d, 0.0*f64::NM],
            [0.0*f64::NM, 0.0*f64::NM,           d])
    }

    /// Create a rectangular box with side lengths $l$, $w$ and $h$
    ///
    /// The volume is $l w h$
    pub fn rectangular(
        l: f64::Length,
        w: f64::Length,
        h: f64::Length
    ) -> Self {
        Pbc(
            [          l, 0.0*f64::NM, 0.0*f64::NM],
            [0.0*f64::NM,           w, 0.0*f64::NM],
            [0.0*f64::NM, 0.0*f64::NM,           h]
        )
    }

    /// Create a rhombic dodecahedral box with a square in the $xy$-plane
    ///
    /// The volume is $\frac{1}{2} \sqrt{2} d^3 \approx 0.707 d^3$
    pub fn rhombic_dodecahedral_xysquare(
        d: f64::Length
    ) -> Self {
        Pbc(
            [          d, 0.0*f64::NM,       0.0*f64::NM],
            [0.0*f64::NM,           d,       0.0*f64::NM],
            [      d/2.0,       d/2.0, d*2f64.sqrt()/2.0]
        )
    }

    /// Create a rhombic dodecahedral box with a hexagon in the $xy$-plane
    ///
    /// The volume is $\frac{1}{2} \sqrt{2} d^3 \approx 0.707 d^3$
    pub fn rhombic_dodecahedral_xyhex(
        d: f64::Length
    ) -> Self {
        Pbc(
            [     d,       0.0*f64::NM,       0.0*f64::NM],
            [ d/2.0, d*3f64.sqrt()/2.0,       0.0*f64::NM],
            [ d/2.0, d*3f64.sqrt()/6.0, d*6f64.sqrt()/3.0]
        )
    }
}

impl BoundaryConditions for Pbc {

    /// Compute the squared distance between two points with the minimum image convention
    fn dist2(
        &self,
        a: [f64::Length; 3],
        b: [f64::Length; 3],
    ) -> f64::Area {
        let [x1, y1, z1] = a;
        let [x2, y2, z2] = b;

        let Pbc([xi, yi, zi], [xj, yj, zj], [xk, yk, zk]) = *self;

        let mut min = std::f64::INFINITY * f64::NM * f64::NM;

        let dx = x2 - x1;
        let dy = y2 - y1;
        let dz = z2 - z1;

        for &i in &[-1.0, 0.0, 1.0] {
            for &j in &[-1.0, 0.0, 1.0] {
                for &k in &[-1.0, 0.0, 1.0] {
                    min = min.min(
                        (dx + i*xi + j*xj + k*xk).powi(P2::new())
                        + (dy + i*yi + j*yj + k*yk).powi(P2::new())
                        + (dz + i*zi + j*zj + k*zk).powi(P2::new())
                    );
                }
            }
        }

        min
    }

    fn pairlist_checks<P: Pairlist<Self>>(
        &self,
        pairlist: &P
    ) -> Result<()> {
        let cutoff = match pairlist.cutoff() {
            None => return Ok(()),
            Some(cutoff) => cutoff
        };

        let Pbc(a, b, c) = *self;

        let smallest_box_length = [a, b, c]
            .iter()
            .map(|[x, y, z]| (
                x.powi(P2::new())
                + y.powi(P2::new())
                + z.powi(P2::new())
            ).sqrt())
            .fold(
                std::f64::INFINITY * f64::NM,
                |min, x| if min <= x {min} else {x}
            );

        if cutoff < smallest_box_length {
            Ok(())
         } else {
            Err(MinimumImageConventionNotJustified)
         }
    }
}
