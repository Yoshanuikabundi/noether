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

use crate::dim::{
    Sqrt,
    Dimensioned
};

/// A vector in 3 dimensions.
///
/// See trait implementations for implemented operations.
/// Note that the `%` (remainder) operator implements the
/// cross product, and not any sort of modulus arithmetic.
///
/// # Examples
///
/// ```
/// use noether::geom::vec3d::Vec3D;
/// use noether::units::KilojoulePerMolePerNanometer;
///
/// type Force = KilojoulePerMolePerNanometer<f32>;
/// type ForceVec = Vec3D<KilojoulePerMolePerNanometer<f32>>;
///
/// let origin = Vec3D::<Force>::zero();
/// let point: Vec3D<Force> = Vec3D::from(2.0, 1.0, 0.0);
/// let point2 = ForceVec::from(2.0, 1.0, 0.0);
///
/// // assert_eq!(point.clone(), point.clone() - origin);
/// assert_eq!(point, point2);
/// ```
#[derive(PartialEq, Clone, Debug)]
pub struct Vec3D<Q: Dimensioned>
{
    pub x: Q,
    pub y: Q,
    pub z: Q
}

impl<Q, V> Vec3D<Q>
    where
        Q: Copy + Dimensioned<Value=V>,
        V: Copy + From<f32>
{
    /// A vector in 3 dimensions.
    ///
    /// See trait implementations for implemented operations.
    /// Note that the `%` (remainder) operator implements the
    /// cross product, and not any sort of modulus arithmetic.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::f32consts::KJPMNM;
    ///
    /// let point = Vec3D::new(2.0 * KJPMNM, 1.0 * KJPMNM, 0.0 * KJPMNM);
    ///
    /// assert_eq!(point.clone(), point);
    /// ```
    pub fn new(x:Q, y:Q, z:Q) -> Vec3D<Q> {
        Vec3D {
            x,
            y,
            z,
        }
    }

    /// A vector in 3 dimensions.
    ///
    /// See trait implementations for implemented operations.
    /// Note that the `%` (remainder) operator implements the
    /// cross product, and not any sort of modulus arithmetic.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::KilojoulePerMolePerNanometer;
    ///
    /// type Force = KilojoulePerMolePerNanometer<f32>;
    ///
    /// let point: Vec3D<Force> = Vec3D::from(2.0, 1.0, 0.0);
    ///
    /// assert_eq!(point.clone(), point);
    /// ```
    pub fn from(x: V, y: V, z: V) -> Vec3D<Q> {
        Vec3D {
            x: Q::new(x),
            y: Q::new(y),
            z: Q::new(z)
        }
    }

    /// Create a vector of zeroes.
    pub fn zero() -> Vec3D<Q> {
        Vec3D::<Q>::from(
            V::from(0.0),
            V::from(0.0),
            V::from(0.0)
        )
    }
}

// impl<Q, QR, QQR> Vec3D<Q>
//     where
//         Q: Copy + Dimensioned + Mul<QR, Output=QQR>,
//         QR: Copy + Dimensioned,
//         QQR: Copy + Dimensioned + Add<QQR, Output=QQR> + Sub<QQR, Output=QQR> + Div<QR, Output=Q>
// {
//     /// Dot product
//     ///
//     /// See also the implementation of the `Mul` trait for
//     /// dot products with the `*` operator.
//     ///
//     /// # Examples
//     ///
//     /// ```
//     /// use noether::geom::vec3d::Vec3D;
//     /// use noether::units::Nanometer;
//     /// use noether::units::f32consts::NM2;
//     ///
//     /// type PosVec = Vec3D<Nanometer<f32>>;
//     ///
//     /// let a = PosVec::from(1.0, 0.0, 0.0);
//     /// let b = PosVec::from(1.0, 0.0, 0.0);
//     /// assert_eq!(a.dot(&b), 1.0 * NM2);
//     ///
//     /// let a = PosVec::from(1.0, 0.0, 0.0);
//     /// let b = PosVec::from(0.0, 1.0, 0.0);
//     /// assert_eq!(a.dot(&b), 0.0 * NM2);
//     ///
//     /// let a = PosVec::from(  3.0, -3.0,  1.0);
//     /// let b = PosVec::from(  4.0,  9.0,  2.0);
//     /// assert_eq!(a.dot(&b), -13.0 * NM2);
//     /// ```
//     pub fn dot(&self, other:&Vec3D<QR>) -> QQR {
//         self.x * other.x
//         + self.y * other.y
//         + self.z * other.z
//     }

//     /// Cross product
//     ///
//     /// See also the implementation of the `Rem` trait for
//     /// cross products with the `%` operator.
//     ///
//     /// # Examples
//     ///
//     /// ```
//     /// use noether::geom::vec3d::Vec3D;
//     /// use noether::units::Nanometer;
//     /// use noether::units::Nanometer2;
//     /// use noether::units::f32consts::NM2;
//     ///
//     /// type PosVec = Vec3D<Nanometer<f32>>;
//     /// type Pos2Vec = Vec3D<Nanometer2<f32>>;
//     ///
//     /// let a = PosVec::from(1.0, 0.0, 0.0);
//     /// let b = PosVec::from(1.0, 0.0, 0.0);
//     /// assert_eq!(
//     ///     a.cross(b),
//     ///     Pos2Vec::from(0.0, 0.0, 0.0)
//     /// );
//     ///
//     /// let a = PosVec::from(1.0, 0.0, 0.0);
//     /// let b = PosVec::from(0.0, 1.0, 0.0);
//     /// assert_eq!(
//     ///     a.cross(b),
//     ///     Pos2Vec::from(0.0, 0.0, 1.0)
//     /// );
//     ///
//     /// let a = PosVec::from(  3.0, -3.0,  1.0);
//     /// let b = PosVec::from(  4.0,  9.0,  2.0);
//     /// assert_eq!(
//     ///     a.cross(b),
//     ///     Pos2Vec::from(-15.0, -2.0, 39.0)
//     /// );
//     /// ```
//     pub fn cross(self, other:Vec3D<QR>) -> Vec3D<QQR> {
//         Vec3D {
//             x: (self.y * other.z) - (self.z * other.y),
//             y: (self.z * other.x) - (self.x * other.z),
//             z: (self.x * other.y) - (self.y * other.x),
//         }
//     }
// }

impl<Q, Q2, V> Vec3D<Q>
    where
        Q: Copy + Dimensioned<Value=V> + Mul<Q, Output=Q2> + DivAssign<V> + Div<V, Output=Q>,
        Q2: Copy + Dimensioned + Add<Q2, Output=Q2> + Sub<Q2, Output=Q2> + Sqrt<Output=Q>,
        V: Copy
{

    /// The squared cartesian norm of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::KilojoulePerMolePerNanometer;
    /// use noether::units::f32consts::KJPMNM;
    ///
    /// type ForceVec = Vec3D<KilojoulePerMolePerNanometer<f32>>;
    ///
    /// let point = ForceVec::from(2.0, 1.0, 0.0);
    /// assert!(point.norm2() == 5.0 * KJPMNM * KJPMNM);
    /// ```
    pub fn norm2(&self) -> Q2 {
        self.x * self.x
        + self.y * self.y
        + self.z * self.z
    }

    /// The cartesian norm of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::KilojoulePerMolePerNanometer;
    /// use noether::units::f32consts::KJPMNM;
    ///
    /// type ForceVec = Vec3D<KilojoulePerMolePerNanometer<f32>>;
    ///
    /// let point = ForceVec::from(2.0, -1.0, 2.0);
    /// assert_eq!(point.norm2(), 9.0 * KJPMNM * KJPMNM);
    /// assert_eq!(point.norm(), 3.0 * KJPMNM);
    /// ```
    pub fn norm(&self) -> Q {
        self.norm2().sqrt()
    }

    /// Normalize vector in place. Produces a unit vector
    /// in the same direction and _with the same units_ as
    /// self.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::KilojoulePerMolePerNanometer;
    /// use noether::units::f64consts::KJPMNM;
    ///
    /// type ForceVec = Vec3D<KilojoulePerMolePerNanometer<f64>>;
    ///
    /// let mut push = ForceVec::from(2.0, -1.0, 2.0);
    /// push.normalize();
    /// assert_eq!(
    ///     push,
    ///     ForceVec::from(2.0/3.0, -1.0/3.0, 2.0/3.0)
    /// );
    /// assert_eq!(push.norm(), 1.0 * KJPMNM);
    /// ```
    pub fn normalize(&mut self) -> &Self {
        let n = self.norm().value_unsafe().clone();
        self.x /= n;
        self.y /= n;
        self.z /= n;
        self
    }

    /// Normalize vector in place after taking ownership
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::KilojoulePerMolePerNanometer;
    /// use noether::units::f32consts::KJPMNM;
    ///
    /// type ForceVec = Vec3D<KilojoulePerMolePerNanometer<f32>>;
    ///
    /// let mut push = ForceVec::from(2.0, -1.0, 2.0);
    /// let push = push.normalize_into();
    /// assert_eq!(
    ///     push,
    ///     ForceVec::from(2.0/3.0, -1.0/3.0, 2.0/3.0)
    /// );
    /// assert_eq!(push.norm(), 1.0 * KJPMNM);
    /// ```
    pub fn normalize_into(mut self) -> Self {
        let n = self.norm().value_unsafe().clone();
        self.x /= n;
        self.y /= n;
        self.z /= n;
        self
    }

    /// Return a new, normalized vector
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// let point = PosVec::from(2.0, -1.0, 2.0);
    /// assert_eq!(
    ///     point.normalized(),
    ///     Vec3D::from(2.0/3.0, -1.0/3.0, 2.0/3.0)
    /// );
    /// assert_eq!(point.normalized().norm(), 1.0 * NM);
    /// ```
    pub fn normalized(&self) -> Vec3D<Q> {
        let n = self.norm().value_unsafe().clone();
        Vec3D {
            x: self.x / n,
            y: self.y / n,
            z: self.z / n,
        }
    }
}

impl<Q> Add for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Add<Q, Output=Q>
{
    type Output = Vec3D<Q>;

    /// Vector addition with the `+` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from(  4.0,  3.9, 100204.23)
    ///     + PosVec::from(2.6, -1.9,      0.00),
    /// // -------------------------------------
    ///     PosVec::from(  6.6,  2.0, 100204.23)
    /// );
    /// ```
    fn add(self, other: Vec3D<Q>) ->  Vec3D<Q> {
        Vec3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<Q> AddAssign for Vec3D<Q>
    where
        Q: Copy + Dimensioned + AddAssign<Q>
{
    /// Vector addition assignment with the `+=` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// let mut a = PosVec::from(4.0, 3.9, 100204.23);
    /// let b = PosVec::from(2.6, -1.9, 0.0);
    ///
    /// a += b;
    ///
    /// assert_eq!(a, PosVec::from(6.6, 2.0, 100204.23));
    /// ```
    fn add_assign(&mut self, other:Vec3D<Q>) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl<Q> Sub for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Sub<Q, Output=Q>
{
    type Output = Vec3D<Q>;

    /// Vector subtraction with the `-` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from(   4.0, 3.9, 100204.23)
    ///     - PosVec::from(-2.6, 1.9,      0.00),
    /// // -------------------------------------
    ///     PosVec::from(  6.6,  2.0, 100204.23)
    /// );
    /// ```
    fn sub(self, other: Vec3D<Q>) ->  Vec3D<Q> {
        Vec3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, Q> Sub for &'a Vec3D<Q>
    where
        Q: Copy + Dimensioned + Sub<Q, Output=Q>
{
    type Output = Vec3D<Q>;

    /// Vector subtraction with the `-` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from(   4.0, 3.9, 100204.23)
    ///     - PosVec::from(-2.6, 1.9,      0.00),
    /// // -------------------------------------
    ///     PosVec::from(  6.6,  2.0, 100204.23)
    /// );
    /// ```
    fn sub(self, other: &Vec3D<Q>) ->  Vec3D<Q> {
        Vec3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<Q> SubAssign for Vec3D<Q>
    where
        Q: Copy + Dimensioned + SubAssign<Q>
{
    /// Vector subtraction assignment with the `-=` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// let mut a = PosVec::from(4.0, 3.9, 100204.23);
    /// let b = PosVec::from(-2.6, 1.9, 0.0);
    ///
    /// a -= b;
    ///
    /// assert_eq!(a, PosVec::from(6.6, 2.0, 100204.23));
    /// ```
    fn sub_assign(&mut self, other:Vec3D<Q>) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl<Q, QR, QQR> Mul<QR> for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Mul<QR, Output=QQR>,
        QR: Copy,
        QQR: Copy + Dimensioned
{
    type Output = Vec3D<QQR>;

    /// Multiplication of a vector by a scalar with the `*` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::{
    ///     Nanometer,
    ///     ElemChargeNanometer
    /// };
    /// use noether::units::f32consts::E;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from( 4.0,  3.9,  -804.23) * 4.0,
    ///     PosVec::from(16.0, 15.6, -3216.92)
    /// );
    ///
    /// type DielectricVec = Vec3D<ElemChargeNanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from( 4.0,  3.9,  -804.23) * 4.0 * E,
    ///     DielectricVec::from(16.0, 15.6, -3216.92)
    /// );
    /// ```
    fn mul(self, other:QR) -> Vec3D<QQR> {
        Vec3D {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}
impl<Q, V> MulAssign<V> for Vec3D<Q>
    where
        Q: Copy + Dimensioned + MulAssign<V>,
        V: Copy
{

    /// Multiplication assignment of a vector by a dimensionless
    /// scalar with the `*=` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// let mut a = PosVec::from(4.0, 3.9, -804.23);
    ///
    /// a *= 4.0;
    ///
    /// assert_eq!(a, PosVec::from(16.0, 15.6, -3216.92));
    ///
    /// a *= (NM/NM);
    ///
    /// assert_eq!(a, PosVec::from(16.0, 15.6, -3216.92));
    /// ```
    fn mul_assign(&mut self, other:V) {
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

// TODO: Multiplication of a scalar by a vector
// /// Multiplication of a scalar by a vector with the `*` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3d::Vec3D;
// /// use noether::units::{
// ///     Nanometer,
// ///     ElemChargeNanometer
// /// };
// /// use noether::units::f32consts::E;
// ///
// /// type PosVec = Vec3D<Nanometer<f32>>;
// ///
// /// assert_eq!(
// ///     4 * PosVec::from( 4.0,  3.9,  -804.23),
// ///     PosVec::from(16.0, 15.6, -3216.92)
// /// );
// ///
// /// type DielectricVec = Vec3D<ElemChargeNanometer<f32>>;
// ///
// /// assert_eq!(
// ///     E * 4.0 * PosVec::from( 4.0,  3.9,  -804.23),
// ///     DielectricVec::from(16.0, 15.6, -3216.92)
// /// );
// /// ```
// impl<Q, QR, QQR> Mul<Vec3D<QR>> for Q
//     where
//         Q: Copy + Mul<QR, Output=QQR>,
//         QR: Copy + Dimensioned,
//         QQR: Copy + Dimensioned
// {
//     type Output = Vec3D<QQR>;

//     fn mul(self, other:Vec3D<QR>) -> Vec3D<QQR> {
//         Vec3D {
//             x: self * other.x,
//             y: self * other.y,
//             z: self * other.z,
//         }
//     }
// }

impl<Q, QR, QQR> Rem<Vec3D<QR>> for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Mul<QR, Output=QQR>,
        QR: Copy + Dimensioned,
        QQR: Copy + Dimensioned + Sub<QQR, Output=QQR>
{
    type Output = Vec3D<QQR>;

    /// Cross product of two vectors with the `%` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::{
    ///     Nanometer,
    ///     Nanometer2
    /// };
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    /// type AreaVec = Vec3D<Nanometer2<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from(  1.0, 0.0, 0.0)
    ///     % PosVec::from(1.0, 0.0, 0.0),
    /// // -------------------------------
    ///     AreaVec::from( 0.0, 0.0, 0.0)
    /// );
    /// assert_eq!(
    ///     PosVec::from(  1.0, 0.0, 0.0)
    ///     % PosVec::from(0.0, 1.0, 0.0),
    /// // -------------------------------
    ///     AreaVec::from( 0.0, 0.0, 1.0)
    /// );
    /// assert_eq!(
    ///     PosVec::from(   3.0, -3.0,  1.0)
    ///     % PosVec::from( 4.0,  9.0,  2.0),
    /// // --------------------------------
    ///     AreaVec::from(-15.0, -2.0, 39.0)
    /// );
    /// ```
    fn rem(self, other:Vec3D<QR>) -> Vec3D<QQR> {
        Vec3D {
            x: (self.y * other.z) - (self.z * other.y),
            y: (self.z * other.x) - (self.x * other.z),
            z: (self.x * other.y) - (self.y * other.x),
        }
    }
}

impl<Q, QR> RemAssign<Vec3D<QR>> for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Mul<QR, Output=Q> + Sub<Q, Output=Q>,
        QR: Copy + Dimensioned
{
    /// Cross product  assignment of two vectors
    /// with the `%=` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::{
    ///     Nanometer,
    ///     Unitless
    /// };
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    /// type NodimVec = Vec3D<Unitless<f32>>;
    ///
    /// let mut a = PosVec::from(1.0, 0.0, 0.0);
    /// let b = NodimVec::from(1.0, 0.0, 0.0);
    /// a %= b;
    /// assert_eq!(a, PosVec::from(0.0, 0.0, 0.0));
    ///
    /// let mut a = PosVec::from(1.0, 0.0, 0.0);
    /// let b = NodimVec::from(0.0, 1.0, 0.0);
    /// a %= b;
    /// assert_eq!(a, PosVec::from(0.0, 0.0, 1.0));
    ///
    /// let mut a = PosVec::from(  3.0, -3.0,  1.0);
    /// let b = NodimVec::from(  4.0,  9.0,  2.0);
    /// a %= b;
    /// assert_eq!(a, PosVec::from(-15.0, -2.0, 39.0));
    /// ```
    fn rem_assign(&mut self, other:Vec3D<QR>) {
        let x = (self.y * other.z) - (self.z * other.y);
        let y = (self.z * other.x) - (self.x * other.z);
        let z = (self.x * other.y) - (self.y * other.x);
        self.x = x;
        self.y = y;
        self.z = z;
    }
}

impl<Q, QR, QQR> Mul<Vec3D<QR>> for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Mul<QR, Output=QQR>,
        QR: Copy + Dimensioned,
        QQR: Copy + Dimensioned + Add<QQR, Output=QQR>
{
    type Output = QQR;

    /// Dot product of two vectors with the `*` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from(  1.0, 0.0, 0.0)
    ///     * PosVec::from(1.0, 0.0, 0.0),
    /// // -------------------------------
    ///     1.0 * NM * NM
    /// );
    /// assert_eq!(
    ///     PosVec::from(  1.0, 0.0, 0.0)
    ///     * PosVec::from(0.0, 1.0, 0.0),
    /// // -------------------------------
    ///     0.0 * NM * NM
    /// );
    /// assert_eq!(
    ///     PosVec::from(   3.0, -3.0,  1.0)
    ///     * PosVec::from( 4.0,  9.0,  2.0),
    /// // --------------------------------
    ///     -13.0 * NM * NM
    /// );
    /// ```
    fn mul(self, other:Vec3D<QR>) -> QQR {
        self.x * other.x
        + self.y * other.y
        + self.z * other.z
    }
}

impl<QQR, QR, Q> Div<QR> for Vec3D<QQR>
    where
        QQR: Copy + Dimensioned + Div<QR, Output=Q>,
        QR: Copy,
        Q: Copy + Dimensioned
{
    type Output = Vec3D<Q>;

    /// Division of a vector by a scalar with the `/` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::{
    ///     Nanometer,
    ///     ElemChargeNanometer
    /// };
    /// use noether::units::f32consts::E;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     PosVec::from(16.0, 15.6, -3216.92) / 4.0,
    ///     PosVec::from( 4.0,  3.9,  -804.23)
    /// );
    ///
    /// type DielectricVec = Vec3D<ElemChargeNanometer<f32>>;
    ///
    /// assert_eq!(
    ///     DielectricVec::from(16.0, 15.6, -3216.92) / 4.0 / E,
    ///     PosVec::from( 4.0,  3.9,  -804.23)
    /// );
    /// ```
    fn div(self, other:QR) -> Vec3D<Q> {
        Vec3D {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl<Q, V> DivAssign<V> for Vec3D<Q>
    where
        Q: Copy + Dimensioned + DivAssign<V>,
        V: Copy
{

    /// Division assignment of a vector by a dimensionless
    /// scalar with the `/=` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// let mut a = PosVec::from(16.0, 15.6, -3216.92);
    ///
    /// a /= 4.0;
    ///
    /// assert_eq!(a, PosVec::from(4.0, 3.9, -804.23));
    ///
    /// a /= (NM/NM);
    ///
    /// assert_eq!(a, PosVec::from(4.0, 3.9, -804.23));
    /// ```
    fn div_assign(&mut self, other:V) {
        self.x /= other;
        self.y /= other;
        self.z /= other;
    }
}

// TODO: Implement division of scalar by vector with generics
// /// Division of a dimensionless scalar by a vector with the `/` operator.
// ///
// /// # Examples
// ///
// /// ```
// /// use noether::geom::vec3d::Vec3D;
// ///
// /// assert_eq!(
// ///     4.0 / Vec3D::from(4.0, 16.0, -0.0025),
// ///           Vec3D::from(1.0, 0.25, -1600.0)
// /// );
// impl Div<Vec3D> for V {
//     type Output = Vec3D;
//
//     fn div(self, other:Vec3D) -> Vec3D {
//         Vec3D {
//             x: self / other.x,
//             y: self / other.y,
//             z: self / other.z,
//         }
//     }
// }

impl<Q> Neg for Vec3D<Q>
    where
        Q: Copy + Dimensioned + Neg<Output=Q>
{
    type Output = Vec3D<Q>;
    /// Negation of a vector with the `-` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::units::Nanometer;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// assert_eq!(
    ///     - PosVec::from( 4.0, 0.0, -0.0025),
    ///       PosVec::from(-4.0, 0.0,  0.0025)
    /// );
    fn neg(self) -> Vec3D<Q> {
        Vec3D {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
impl<Q> Index<VecIdx> for Vec3D<Q>
    where
        Q: Copy + Dimensioned
{
    type Output = Q;

    /// Indexing with the `VecIdx` enum. Just in case it's useful.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::vec3d::Vec3D;
    /// use noether::geom::vec3d::VecIdx;
    /// use noether::units::Nanometer;
    /// use noether::units::f32consts::NM;
    ///
    /// type PosVec = Vec3D<Nanometer<f32>>;
    ///
    /// let a_vec = PosVec::from(4.0, 0.0, -0.0025);
    ///
    /// assert_eq!(a_vec[VecIdx::X], a_vec.x);
    /// assert_eq!(a_vec[VecIdx::Y], a_vec.y);
    /// assert_eq!(a_vec[VecIdx::Z], a_vec.z);
    /// assert_eq!(a_vec[VecIdx::X],  4.0 * NM);
    /// assert_eq!(a_vec[VecIdx::Y],  0.0 * NM);
    /// assert_eq!(a_vec[VecIdx::Z], -0.0025 * NM);
    /// ```
    fn index(&self, index: VecIdx) -> &Q {
        match index {
            VecIdx::X => &self.x,
            VecIdx::Y => &self.y,
            VecIdx::Z => &self.z,
        }
    }
}

/// Index enum for accessing Vec3D with indexing syntax
pub enum VecIdx { X, Y, Z }
