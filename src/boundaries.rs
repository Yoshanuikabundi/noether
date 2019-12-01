use crate::units::f64;

use uom::typenum::consts::P2;


pub trait BoundaryConditions {
    /// Compute the distance between two points given the boundary conditions
    fn dist(
        &self,
        a: [f64::Length; 3],
        b: [f64::Length; 3],
    ) -> f64::Length {
        self.dist2(a, b).sqrt()
    }

    /// Compute the squared distance between two points given the boundary conditions
    fn dist2(
        &self,
        a: [f64::Length; 3],
        b: [f64::Length; 3],
    ) -> f64::Area;

}

pub struct NoBounds;

impl BoundaryConditions for NoBounds {

    /// Compute the squared distance between two points
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
}

pub struct Pbc(
    [f64::Length; 3],
    [f64::Length; 3],
    [f64::Length; 3]
);

impl Pbc {

    /// Create a cubic box with side length d
    ///
    /// The volume is d^3
    pub fn cubic(d: f64::Length) -> Self {
        Pbc(
            [          d, 0.0*f64::NM, 0.0*f64::NM],
            [0.0*f64::NM,           d, 0.0*f64::NM],
            [0.0*f64::NM, 0.0*f64::NM,           d])
    }

    /// Create a cubic box with side lengths l, w and h
    ///
    /// The volume is lwh
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

    /// Create a rhombic dodecahedral box with a square in the xy-plane
    ///
    /// The volume is ~0.707 d^3
    pub fn rhombic_dodecahedral_xysquare(
        d: f64::Length
    ) -> Self {
        Pbc(
            [          d, 0.0*f64::NM,       0.0*f64::NM],
            [0.0*f64::NM,           d,       0.0*f64::NM],
            [      d/2.0,       d/2.0, d*2f64.sqrt()/2.0]
        )
    }

    /// Create a rhombic dodecahedral box with a hexagon in the xy-plane
    ///
    /// The volume is ~0.707 d^3
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

    /// Compute the minimum squared image distance between two points
    fn dist2(
        &self,
        a: [f64::Length; 3],
        b: [f64::Length; 3],
    ) -> f64::Area {
        let [x1, y1, z1] = a;
        let [x2, y2, z2] = b;

        let Pbc([xi, yi, zi], [xj, yj, zj], [xk, yk, zk]) = self;

        let mut min = std::f64::INFINITY * f64::NM * f64::NM;

        for &i in &[-1.0, 0.0, 1.0] {
            for &j in &[-1.0, 0.0, 1.0] {
                for &k in &[-1.0, 0.0, 1.0] {
                    min = min.min(
                        (x2 - x1 + i * *xi + j * *xj + k * *xk).powi(P2::new())
                        + (y2 - y1 + i * *yi + j * *yj + k * *yk).powi(P2::new())
                        + (z2 - z1 + i * *zi + j * *zj + k * *zk).powi(P2::new())
                    );
                }
            }
        }

        min
    }
}
