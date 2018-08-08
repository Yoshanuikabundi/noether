use std::ops::{
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Mul,
    MulAssign,
    Div,
    DivAssign,
    Rem,
    RemAssign,
    Neg,
    Index,
};

use units;
use units::GMXQ;
use units::GMXU;
use units::Dimension;

use typenum;
use typenum::operator_aliases::Square;

use std::marker::PhantomData;

// use typenum::consts::{
//     Z0,
//     P1,
//     P2,
//     N2,
//     N4,
// };
// type ForceDim = GMXQ<P1, P1, N2, Z0, Z0>;
// type Force2Dim = GMXQ<P2, P2, N4, Z0, Z0>;

type V = f32;
type Q<NM, DA, PS, E, K> = units::Quantity<GMXQ<NM, DA, PS, E, K>, GMXU<V>, V>;
// type Q = units::Quantity<ForceDim, GMXU<V>, V>;
// type Q2 = units::Quantity<Force2Dim, GMXU<V>, V>;

// const UNIT: Q = units::f32::KJPMNM;


/// A vector in 3 dimensions.
///
/// See trait implementations for implemented operations.
/// Note that the `%` (remainder) operator implements the
/// cross product, and not any sort of modulus arithmetic.
///
/// # Examples
///
/// ```
/// extern crate noether;
/// extern crate typenum;
///
/// use noether::geom::Vec3D;
/// use noether::units::f32::NM;
/// use noether::units::f32::KJPMNM;
/// use typenum::consts::{Z0, P1};
///
/// let origin: Vec3D<P1, Z0, Z0, Z0, Z0> = Vec3D::from(0.0, 0.0, 0.0);
/// let point = Vec3D::new(2.0 * NM, 1.0 * NM, 0.0 * NM);
///
/// // assert_eq!(point.clone(), point.clone() - origin);
/// assert_eq!(point.clone(), point.clone());
/// ```
#[derive(PartialEq, Debug, Clone)]
pub struct Vec3D<NM, DA, PS, E, K>
    where
        NM: typenum::Integer,
        DA: typenum::Integer,
        PS: typenum::Integer,
        E: typenum::Integer,
        K: typenum::Integer
{
    pub x: Q<NM, DA, PS, E, K>,
    pub y: Q<NM, DA, PS, E, K>,
    pub z: Q<NM, DA, PS, E, K>
}

impl<NM, DA, PS, E, K> Vec3D<NM, DA, PS, E, K>
    where
        NM: typenum::Integer,
        DA: typenum::Integer,
        PS: typenum::Integer,
        E: typenum::Integer,
        K: typenum::Integer
{
    pub fn from(x:V, y:V, z:V) -> Vec3D<NM, DA, PS, E, K> {
        let unit: Q<NM, DA, PS, E, K> = units::Quantity {
            dimension: PhantomData,
            units: PhantomData,
            value: 1.0
        };
        Vec3D {
            x: x * unit,
            y: y * unit,
            z: z * unit,
        }
    }

    pub fn new(
        x:Q<NM, DA, PS, E, K>,
        y:Q<NM, DA, PS, E, K>,
        z:Q<NM, DA, PS, E, K>
    ) -> Vec3D<NM, DA, PS, E, K> {
        Vec3D {
            x,
            y,
            z,
        }
    }
}

//     /// The squared cartesian norm of the vector.
//     ///
//     /// # Examples
//     ///
//     /// ```
//     /// use noether::geom::vec3::Vec3D;
//     /// use noether::units::f32::KJPMNM;
//     ///
//     /// let point: Vec3D<P1, P1, N2, Z0, Z0> = Vec3D::from(2.0, 1.0, 0.0);
//     /// assert_eq!(point.norm2(), 5.0 * KJPMNM * KJPMNM);
//     /// ```
//     pub fn norm2(&self) -> Q<Square<NM>, Square<DA>, Square<PS>, Square<E>, Square<K>> {
//         let unit: Q<Square<NM>, Square<DA>, Square<PS>, Square<E>, Square<K>> = units::Quantity {
//             dimension: PhantomData,
//             units: PhantomData,
//             value: 1.0
//         };
//         (self.x.value * self.x.value
//         + self.y.value * self.y.value
//         + self.z.value * self.z.value) * unit
//     }
// }

    // /// Normalize vector in place
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use noether::geom::vec3::Vec3D;
    // /// use noether::units::f32::KJPMNM;
    // ///
    // /// let mut point = Vec3D::from(2.0, -1.0, 2.0);
    // /// point.normalize();
    // /// assert_eq!(
    // ///     point,
    // ///     Vec3D::from(2.0/3.0, -1.0/3.0, 2.0/3.0)
    // /// );
    // /// assert_eq!(point.norm(), 1.0 * KJPMNM);
    // /// ```
    // pub fn normalize(&mut self) -> &Self {
    //     let n = self.norm();
    //     self.x /= n.value;
    //     self.y /= n.value;
    //     self.z /= n.value;
    //     self
    // }

    // /// Normalize vector in place after taking ownership
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use noether::geom::vec3::Vec3D;
    // /// use noether::units::f32::KJPMNM;
    // ///
    // /// let point = Vec3D::from(2.0, -1.0, 2.0);
    // /// let point = point.normalize_into();
    // /// assert_eq!(
    // ///     point,
    // ///     Vec3D::from(2.0/3.0, -1.0/3.0, 2.0/3.0)
    // /// );
    // /// assert_eq!(point.norm(), 1.0 * KJPMNM);
    // /// ```
    // pub fn normalize_into(mut self) -> Self {
    //     let n = self.norm();
    //     self.x /= n.value;
    //     self.y /= n.value;
    //     self.z /= n.value;
    //     self
    // }

    // /// Return a new, normalized vector
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use noether::geom::vec3::Vec3D;
    // /// use noether::units::f32::KJPMNM;
    // ///
    // /// let point = Vec3D::from(2.0, -1.0, 2.0);
    // /// assert_eq!(
    // ///     point.normalized(),
    // ///     Vec3D::from(2.0/3.0, -1.0/3.0, 2.0/3.0)
    // /// );
    // /// assert_eq!(point.normalized().norm(), 1.0 * KJPMNM);
    // /// ```
    // pub fn normalized(&self) -> Vec3D<D> {
    //     let n = self.norm();
    //     Vec3D {
    //         x: self.x / n.value,
    //         y: self.y / n.value,
    //         z: self.z / n.value,
    //     }
    // }

// //     /// Cross product
// //     ///
// //     /// See also the implementation of the `Rem` trait for
// //     /// cross products with the `%` operator.
// //     ///
// //     /// # Examples
// //     ///
// //     /// ```
// //     /// use noether::geom::vec3::Vec3D;
// //     ///
// //     /// let a = Vec3D::from(1.0, 0.0, 0.0);
// //     /// let b = Vec3D::from(1.0, 0.0, 0.0);
// //     /// assert_eq!(
// //     ///     a.cross(b),
// //     ///         Vec3D::from(0.0, 0.0, 0.0)
// //     /// );
// //     ///
// //     /// let a = Vec3D::from(1.0, 0.0, 0.0);
// //     /// let b = Vec3D::from(0.0, 1.0, 0.0);
// //     /// assert_eq!(
// //     ///     a.cross(b),
// //     ///     Vec3D::from(0.0, 0.0, 1.0)
// //     /// );
// //     ///
// //     /// let a = Vec3D::from(  3.0, -3.0,  1.0);
// //     /// let b = Vec3D::from(  4.0,  9.0,  2.0);
// //     /// assert_eq!(
// //     ///     a.cross(b),
// //     ///     Vec3D::from(-15.0, -2.0, 39.0)
// //     /// );
// //     /// ```
// //     pub fn cross(self, other:Vec3D) -> Vec3D {
// //         self % other
// //     }

//     /// Dot product
//     ///
//     /// See also the implementation of the `Mul` trait for
//     /// dot products with the `*` operator.
//     ///
//     /// # Examples
//     ///
//     /// ```
//     /// use noether::geom::vec3::Vec3D;
//     /// use noether::units::f32::KJPMNM;
//     ///
//     /// let a = Vec3D::from(1.0, 0.0, 0.0);
//     /// let b = Vec3D::from(1.0, 0.0, 0.0);
//     /// assert_eq!(a.dot(b), 1.0 * KJPMNM * KJPMNM);
//     ///
//     /// let a = Vec3D::from(1.0, 0.0, 0.0);
//     /// let b = Vec3D::from(0.0, 1.0, 0.0);
//     /// assert_eq!(a.dot(b), 0.0 * KJPMNM * KJPMNM);
//     ///
//     /// let a = Vec3D::from(  3.0, -3.0,  1.0);
//     /// let b = Vec3D::from(  4.0,  9.0,  2.0);
//     /// assert_eq!(a.dot(b), -13.0 * KJPMNM * KJPMNM);
//     /// ```
//     pub fn dot(self, other:Vec3D) -> Q2 {
//         self * other
//     }

//     /// The cartesian norm of the vector.
//     ///
//     /// # Examples
//     ///
//     /// ```
//     /// use noether::geom::vec3::Vec3D;
//     /// use noether::units::f32::KJPMNM;
//     ///
//     /// let point = Vec3D::from(2.0, -1.0, 2.0);
//     /// assert_eq!(point.norm2(), 9.0 * KJPMNM * KJPMNM);
//     /// assert_eq!(point.norm(), 3.0 * KJPMNM);
//     /// ```
//     pub fn norm(&self) -> Q {
//         self.norm2().value.sqrt() * UNIT
//     }
// }

/// Vector addition with the `+` operator.
///
/// # Examples
///
/// ```
/// extern crate noether;
/// extern crate typenum;
///
/// use noether::geom::vec3::Vec3D;
/// use typenum::consts::*;
///
/// assert_eq!(
///     Vec3D::from(  4.0,  3.9, 100204.23) as Vec3D<P1, Z0, Z0, Z0, Z0>
///     + Vec3D::from(2.6, -1.9,      0.00),
/// // -------------------------------------
///     Vec3D::from(  6.6,  2.0, 100204.23) as Vec3D<P1, Z0, Z0, Z0, Z0>
/// );
/// ```
impl<NM, DA, PS, E, K> Add for Vec3D<NM, DA, PS, E, K>
    where
        NM: typenum::Integer,
        DA: typenum::Integer,
        PS: typenum::Integer,
        E: typenum::Integer,
        K: typenum::Integer
{
    type Output = Vec3D<NM, DA, PS, E, K>;

    fn add(self, other:Vec3D<NM, DA, PS, E, K>) -> Vec3D<NM, DA, PS, E, K> {
        Vec3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

/// Vector addition assignment with the `+=` operator.
///
/// # Examples
///
/// ```
/// extern crate noether;
/// extern crate typenum;
///
/// use typenum::consts::*;
/// use noether::geom::vec3::Vec3D;
/// use noether::units::f32::NM;
///
/// let mut a = Vec3D::new(4.0 * NM, 3.9 * NM, 100204.23 * NM);
/// let b = Vec3D::from(2.6, -1.9, 0.0);
///
/// a += b;
///
/// assert_eq!(a, Vec3D::from(6.6, 2.0, 100204.23) as Vec3D<P1, Z0, Z0, Z0, Z0>);
/// ```
impl<NM, DA, PS, E, K> AddAssign for Vec3D<NM, DA, PS, E, K>
    where
        NM: typenum::Integer,
        DA: typenum::Integer,
        PS: typenum::Integer,
        E: typenum::Integer,
        K: typenum::Integer
{
    fn add_assign(&mut self, other:Vec3D<NM, DA, PS, E, K>) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

// /// Vector subtraction with the `-` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// assert_eq!(
// ///     Vec3D::from(   4.0, 3.9, 100204.23)
// ///     - Vec3D::from(-2.6, 1.9,      0.00),
// /// // -------------------------------------
// ///     Vec3D::from(   6.6, 2.0, 100204.23)
// /// );
// /// ```
// impl Sub for Vec3D {
//     type Output = Vec3D;

//     fn sub(self, other:Vec3D) -> Vec3D {
//         Vec3D {
//             x: self.x - other.x,
//             y: self.y - other.y,
//             z: self.z - other.z,
//         }
//     }
// }

// /// Vector subtraction assignment with the `-=` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// let mut a = Vec3D::from(4.0, 3.9, 100204.23);
// /// let b = Vec3D::from(-2.6, 1.9, -0.00);
// ///
// /// a -= b;
// ///
// /// assert_eq!(a, Vec3D::from(6.6, 2.0, 100204.23));
// /// ```
// impl SubAssign for Vec3D {
//     fn sub_assign(&mut self, other:Vec3D) {
//         self.x -= other.x;
//         self.y -= other.y;
//         self.z -= other.z;
//     }
// }

/// Multiplication of a vector by a dimensionless scalar with the `*` operator.
///
/// # Examples
///
/// ```
/// extern crate noether;
/// extern crate typenum;
///
/// use typenum::consts::*;
/// use noether::geom::vec3::Vec3D;
///
/// assert_eq!(
///     Vec3D::from( 4.0,  3.9,  -804.23) * 4.0,
///     Vec3D::from(16.0, 15.6, -3216.92) as Vec3D<P1, Z0, Z0, Z0, Z0>
/// );
/// ```
impl<NM, DA, PS, E, K> Mul<V> for Vec3D<NM, DA, PS, E, K>
    where
        NM: typenum::Integer,
        DA: typenum::Integer,
        PS: typenum::Integer,
        E: typenum::Integer,
        K: typenum::Integer
{
    type Output = Vec3D<NM, DA, PS, E, K>;

    fn mul(self, other:V) -> Vec3D<NM, DA, PS, E, K> {
        Vec3D {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

// /// Multiplication assignment of a vector by a dimensionless
// /// scalar with the `*=` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// let mut a = Vec3D::from(4.0, 3.9, -804.23);
// ///
// /// a *= 4.0;
// ///
// /// assert_eq!(a, Vec3D::from(16.0, 15.6, -3216.92));
// /// ```
// impl MulAssign<V> for Vec3D {
//     fn mul_assign(&mut self, other:V) {
//         self.x *= other;
//         self.y *= other;
//         self.z *= other;
//     }
// }

// /// Multiplication of a dimensionless scalar by a vector with the `*` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// assert_eq!(
// ///     4.0 * Vec3D::from( 4.0,  3.9,  -804.23),
// ///           Vec3D::from(16.0, 15.6, -3216.92)
// /// );
// /// ```
// impl Mul<Vec3D> for f32 {
//     type Output = Vec3D;

//     fn mul(self, other:Vec3D) -> Vec3D {
//         Vec3D {
//             x: self * other.x,
//             y: self * other.y,
//             z: self * other.z,
//         }
//     }
// }

// // TODO: Set up cross product when we have generics for Vec3Ds
// // /// Cross product of two vectors with the `%` operator.
// // ///
// // /// # Examples
// // ///
// // /// ```
// // /// use noether::geom::vec3::Vec3D;
// // ///
// // /// assert_eq!(
// // ///     Vec3D::from(  1.0, 0.0, 0.0)
// // ///     % Vec3D::from(1.0, 0.0, 0.0),
// // /// // -------------------------------
// // ///     Vec3D::from(  0.0, 0.0, 0.0)
// // /// );
// // /// assert_eq!(
// // ///     Vec3D::from(  1.0, 0.0, 0.0)
// // ///     % Vec3D::from(0.0, 1.0, 0.0),
// // /// // -------------------------------
// // ///     Vec3D::from(  0.0, 0.0, 1.0)
// // /// );
// // /// assert_eq!(
// // ///     Vec3D::from(  3.0, -3.0,  1.0)
// // ///     % Vec3D::from(4.0,  9.0,  2.0),
// // /// // --------------------------------
// // ///     Vec3D::from(-15.0, -2.0, 39.0)
// // /// );
// // /// ```
// // impl Rem<Vec3D> for Vec3D {
// //     type Output = Vec3D;
// //
// //     fn rem(self, other:Vec3D) -> Vec3D {
// //         Vec3D {
// //             x: (self.y * other.z) - (self.z * other.y),
// //             y: (self.z * other.x) - (self.x * other.z),
// //             z: (self.x * other.y) - (self.y * other.x),
// //         }
// //     }
// // }
// //
// // /// Cross product assignment of two vectors with the `%=` operator.
// // ///
// // /// # Examples
// // ///
// // /// ```
// // /// use noether::geom::vec3::Vec3D;
// // ///
// // /// let mut a = Vec3D::from(1.0, 0.0, 0.0);
// // /// let b = Vec3D::from(1.0, 0.0, 0.0);
// // /// a %= b;
// // /// assert_eq!(a, Vec3D::from(0.0, 0.0, 0.0));
// // ///
// // /// let mut a = Vec3D::from(1.0, 0.0, 0.0);
// // /// let b = Vec3D::from(0.0, 1.0, 0.0);
// // /// a %= b;
// // /// assert_eq!(a, Vec3D::from(0.0, 0.0, 1.0));
// // ///
// // /// let mut a = Vec3D::from(  3.0, -3.0,  1.0);
// // /// let b = Vec3D::from(  4.0,  9.0,  2.0);
// // /// a %= b;
// // /// assert_eq!(a, Vec3D::from(-15.0, -2.0, 39.0));
// // /// ```
// // impl RemAssign for Vec3D {
// //     fn rem_assign(&mut self, other:Vec3D) {
// //         let x = (self.y * other.z) - (self.z * other.y);
// //         let y = (self.z * other.x) - (self.x * other.z);
// //         let z = (self.x * other.y) - (self.y * other.x);
// //         self.x = x;
// //         self.y = y;
// //         self.z = z;
// //     }
// // }

// /// Dot product of two vectors with the `*` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// /// use noether::units::f32::KJPMNM;
// ///
// /// assert_eq!(
// ///     Vec3D::from(  1.0, 0.0, 0.0)
// ///     * Vec3D::from(1.0, 0.0, 0.0),
// /// // -------------------------------
// ///     1.0 * KJPMNM * KJPMNM
// /// );
// /// assert_eq!(
// ///     Vec3D::from(  1.0, 0.0, 0.0)
// ///     * Vec3D::from(0.0, 1.0, 0.0),
// /// // -------------------------------
// ///     0.0 * KJPMNM * KJPMNM
// /// );
// /// assert_eq!(
// ///     Vec3D::from(  3.0, -3.0, 1.0)
// ///     * Vec3D::from(4.0,  9.0, 2.0),
// /// // -------------------------------
// ///     -13.0 * KJPMNM * KJPMNM
// /// );
// /// ```
// impl Mul<Vec3D> for Vec3D {
//     type Output = Q2;

//     fn mul(self, other:Vec3D) -> Q2 {
//         self.x * other.x
//         + self.y * other.y
//         + self.z * other.z
//     }
// }

// /// Division of a vector by a dimensionless scalar with the `/` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// assert_eq!(
// ///     Vec3D::from(4.0, 3.9,   -804.23) / 4.0,
// ///     Vec3D::from(1.0, 0.975, -201.0575)
// /// );
// /// ```
// impl Div<V> for Vec3D {
//     type Output = Vec3D;

//     fn div(self, other:V) -> Vec3D {
//         Vec3D {
//             x: self.x / other,
//             y: self.y / other,
//             z: self.z / other,
//         }
//     }
// }

// /// Division assignment of a vector by a dimensionless
// /// scalar with the `/=` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// let mut a = Vec3D::from(4.0, 3.9, -804.23);
// ///
// /// a /= 4.0;
// ///
// /// assert_eq!(a, Vec3D::from(1.0, 0.975, -201.0575));
// /// ```
// impl DivAssign<V> for Vec3D {
//     fn div_assign(&mut self, other:V) {
//         self.x /= other;
//         self.y /= other;
//         self.z /= other;
//     }
// }

// // TODO: Implement division of scalar by vector with generics
// // /// Division of a dimensionless scalar by a vector with the `/` operator.
// // ///
// // /// # Examples
// // ///
// // /// ```
// // /// use noether::geom::vec3::Vec3D;
// // ///
// // /// assert_eq!(
// // ///     4.0 / Vec3D::from(4.0, 16.0, -0.0025),
// // ///           Vec3D::from(1.0, 0.25, -1600.0)
// // /// );
// // impl Div<Vec3D> for V {
// //     type Output = Vec3D;
// //
// //     fn div(self, other:Vec3D) -> Vec3D {
// //         Vec3D {
// //             x: self / other.x,
// //             y: self / other.y,
// //             z: self / other.z,
// //         }
// //     }
// // }

// /// Negation of a vector with the `-` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// ///
// /// assert_eq!(
// ///     - Vec3D::from( 4.0, 0.0, -0.0025),
// ///       Vec3D::from(-4.0, 0.0,  0.0025)
// /// );
// impl Neg for Vec3D {
//     type Output = Vec3D;

//     fn neg(self) -> Vec3D {
//         Vec3D {
//             x: -self.x,
//             y: -self.y,
//             z: -self.z,
//         }
//     }
// }

// /// Indexing with the `VecIdx` enum. Just in case it's useful.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// /// use noether::geom::force::VecIdx;
// /// use noether::units::f32::KJPMNM;
// ///
// /// let a_vec = Vec3D::from(4.0, 0.0, -0.0025);
// ///
// /// assert_eq!(a_vec[VecIdx::X], a_vec.x);
// /// assert_eq!(a_vec[VecIdx::Y], a_vec.y);
// /// assert_eq!(a_vec[VecIdx::Z], a_vec.z);
// /// assert_eq!(a_vec[VecIdx::X],  4.0 * KJPMNM);
// /// assert_eq!(a_vec[VecIdx::Y],  0.0 * KJPMNM);
// /// assert_eq!(a_vec[VecIdx::Z], -0.0025 * KJPMNM);
// /// ```
// impl Index<VecIdx> for Vec3D {
//     type Output = Q;
//     fn index(&self, index: VecIdx) -> &Q {
//         match index {
//             VecIdx::X => &self.x,
//             VecIdx::Y => &self.y,
//             VecIdx::Z => &self.z,
//         }
//     }
// }

// /// Index enum for accessing Vec3D with indexing syntax
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// /// use noether::geom::force::VecIdx;
// ///
// /// let a_vec = Vec3D::from( 4.0, 0.0, -0.0025);
// ///
// /// assert_eq!(a_vec[VecIdx::X], a_vec.x);
// /// assert_eq!(a_vec[VecIdx::Y], a_vec.y);
// /// assert_eq!(a_vec[VecIdx::Z], a_vec.z);
// /// ```
// ///
// /// ```
// /// use noether::geom::vec3::Vec3D;
// /// use noether::geom::force::VecIdx::*;
// /// use noether::units::f32::KJPMNM;
// ///
// /// let a_vec = Vec3D::from( 4.0, 0.0, -0.0025);
// ///
// /// assert_eq!(a_vec[X],  4.0 * KJPMNM);
// /// assert_eq!(a_vec[Y],  0.0 * KJPMNM);
// /// assert_eq!(a_vec[Z], -0.0025 * KJPMNM);
// /// ```
// pub enum VecIdx { X, Y, Z }
