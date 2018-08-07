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
use std::fmt;

/// A vector in 3 dimensions.
///
/// See trait implementations for implemented operations.
/// Note that the `%` (remainder) operator implements the
/// cross product, and not any sort of modulus arithmetic.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
/// use std::fmt::Display;
///
/// let origin = Vec3 (0.0, 0.0, 0.0);
/// let point = Vec3 (2.0, 1.0, 0.0);
/// println!("{}", point);
/// assert_eq!(point.clone(), point - origin);
/// ```
#[derive(PartialEq, Debug, Clone)]
pub struct Vec3 (pub f32, pub f32, pub f32);

impl Vec3 {
    /// The vector's x coordinate.
    pub fn x(&self) -> &f32 {
        &self.0
    }

    /// Mutable access to the vector's x coordinate.
    pub fn x_mut(&mut self) -> &mut f32 {
        &mut self.0
    }

    /// The vector's y coordinate.
    pub fn y(&self) -> &f32 {
        &self.1
    }

    /// Mutable access to the vector's y coordinate.
    pub fn y_mut(&mut self) -> &mut f32 {
        &mut self.1
    }

    /// The vector's z coordinate.
    pub fn z(&self) -> &f32 {
        &self.2
    }

    /// Mutable access to the vector's z coordinate.
    pub fn z_mut(&mut self) -> &mut f32 {
        &mut self.2
    }

    /// Cross product
    ///
    /// See also the implementation of the `Rem` trait for
    /// cross products with the `%` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let a = Vec3 (1.0, 0.0, 0.0);
    /// let b = Vec3 (1.0, 0.0, 0.0);
    /// assert_eq!(
    ///     a.cross(b),
    ///         Vec3 (0.0, 0.0, 0.0)
    /// );
    ///
    /// let a = Vec3 (1.0, 0.0, 0.0);
    /// let b = Vec3 (0.0, 1.0, 0.0);
    /// assert_eq!(
    ///     a.cross(b),
    ///     Vec3     (0.0, 0.0, 1.0)
    /// );
    ///
    /// let a = Vec3 (  3.0, -3.0,  1.0);
    /// let b = Vec3 (  4.0,  9.0,  2.0);
    /// assert_eq!(
    ///     a.cross(b),
    ///     Vec3     (-15.0, -2.0, 39.0)
    /// );
    /// ```
    pub fn cross(self, other:Vec3) -> Vec3 {
        self % other
    }

    /// Dot product
    ///
    /// See also the implementation of the `Mul` trait for
    /// dot products with the `*` operator.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let a = Vec3 (1.0, 0.0, 0.0);
    /// let b = Vec3 (1.0, 0.0, 0.0);
    /// assert_eq!(a.dot(b), 1.0);
    ///
    /// let a = Vec3 (1.0, 0.0, 0.0);
    /// let b = Vec3 (0.0, 1.0, 0.0);
    /// assert_eq!(a.dot(b), 0.0);
    ///
    /// let a = Vec3 (  3.0, -3.0,  1.0);
    /// let b = Vec3 (  4.0,  9.0,  2.0);
    /// assert_eq!(a.dot(b), -13.0);
    /// ```
    pub fn dot(self, other:Vec3) -> f32 {
        self * other
    }

    /// The squared cartesian norm of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let point = Vec3 (2.0, 1.0, 0.0);
    /// assert_eq!(point.norm2(), 5.0);
    /// ```
    pub fn norm2(&self) -> f32 {
        self.0 * self.0
        + self.1 * self.1
        + self.2 * self.2
    }

    /// The cartesian norm of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let point = Vec3 (2.0, -1.0, 2.0);
    /// assert_eq!(point.norm2(), 9.0);
    /// assert_eq!(point.norm(), 3.0);
    /// ```
    pub fn norm(&self) -> f32 {
        self.norm2().sqrt()
    }

    /// Normalize vector in place
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let mut point = Vec3 (2.0, -1.0, 2.0);
    /// point.normalize();
    /// assert_eq!(
    ///     point,
    ///     Vec3 (2.0/3.0, -1.0/3.0, 2.0/3.0)
    /// );
    /// assert_eq!(point.norm(), 1.0);
    /// ```
    pub fn normalize(&mut self) -> &Self {
        let n = self.norm();
        self.0 /= n;
        self.1 /= n;
        self.2 /= n;
        self
    }

    /// Normalize vector in place after taking ownership
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let point = Vec3 (2.0, -1.0, 2.0);
    /// let point = point.normalize_into();
    /// assert_eq!(
    ///     point,
    ///     Vec3 (2.0/3.0, -1.0/3.0, 2.0/3.0)
    /// );
    /// assert_eq!(point.norm(), 1.0);
    /// ```
    pub fn normalize_into(mut self) -> Self {
        let n = self.norm();
        self.0 /= n;
        self.1 /= n;
        self.2 /= n;
        self
    }

    /// Return a new, normalized vector
    ///
    /// # Examples
    ///
    /// ```
    /// use noether::geom::Vec3;
    ///
    /// let point = Vec3 (2.0, -1.0, 2.0);
    /// assert_eq!(
    ///     point.normalized(),
    ///     Vec3 (2.0/3.0, -1.0/3.0, 2.0/3.0)
    /// );
    /// assert_eq!(point.normalized().norm(), 1.0);
    /// ```
    pub fn normalized(&self) -> Vec3 {
        let n = self.norm();
        Vec3 (
            self.0 / n,
            self.1 / n,
            self.2 / n,
        )
    }

}

/// Vector addition with the `+` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     Vec3   (4.0,  3.9, 100204.23)
///     + Vec3 (2.6, -1.9,      0.00),
/// // -------------------------------
///     Vec3   (6.6,  2.0, 100204.23)
/// );
/// ```
impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other:Vec3) -> Vec3 {
        Vec3 (
            self.0 + other.0,
            self.1 + other.1,
            self.2 + other.2,
        )
    }
}

/// Vector addition assignment with the `+=` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// let mut a = Vec3 (4.0,  3.9, 100204.23);
/// let b = Vec3 (2.6, -1.9, 0.00);
///
/// a += b;
///
/// assert_eq!(a, Vec3 (6.6,  2.0, 100204.23));
/// ```
impl AddAssign for Vec3 {
    fn add_assign(&mut self, other:Vec3) {
        self.0 += other.0;
        self.1 += other.1;
        self.2 += other.2;
    }
}

/// Vector subtraction with the `-` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     Vec3   ( 4.0, 3.9, 100204.23)
///     - Vec3 (-2.6, 1.9,      0.00),
/// // -------------------------------
///     Vec3   (6.6,  2.0, 100204.23)
/// );
/// ```
impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, other:Vec3) -> Vec3 {
        Vec3 (
            self.0 - other.0,
            self.1 - other.1,
            self.2 - other.2,
        )
    }
}

/// Vector subtraction assignment with the `-=` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// let mut a = Vec3 (4.0,  3.9, 100204.23);
/// let b = Vec3 (-2.6, 1.9, -0.00);
///
/// a -= b;
///
/// assert_eq!(a, Vec3 (6.6,  2.0, 100204.23));
/// ```
impl SubAssign for Vec3 {
    fn sub_assign(&mut self, other:Vec3) {
        self.0 -= other.0;
        self.1 -= other.1;
        self.2 -= other.2;
    }
}

/// Multiplication of a vector by a scalar with the `*` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     Vec3 ( 4.0,  3.9,  -804.23) * 4.0,
///     Vec3 (16.0, 15.6, -3216.92)
/// );
/// ```
impl Mul<f32> for Vec3 {
    type Output = Vec3;

    fn mul(self, other:f32) -> Vec3 {
        Vec3 (
            self.0 * other,
            self.1 * other,
            self.2 * other,
        )
    }
}

/// Multiplication assignment of a vector by a scalar
/// with the `*=` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// let mut a = Vec3 (4.0, 3.9, -804.23);
///
/// a *= 4.0;
///
/// assert_eq!(a, Vec3 (16.0, 15.6, -3216.92));
/// ```
impl MulAssign<f32> for Vec3 {
    fn mul_assign(&mut self, other:f32) {
        self.0 *= other;
        self.1 *= other;
        self.2 *= other;
    }
}

/// Multiplication of a scalar by a vector with the `*` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     4.0 * Vec3 ( 4.0,  3.9,  -804.23),
///           Vec3 (16.0, 15.6, -3216.92)
/// );
/// ```
impl Mul<Vec3> for f32 {
    type Output = Vec3;

    fn mul(self, other:Vec3) -> Vec3 {
        Vec3 (
            self * other.0,
            self * other.1,
            self * other.2,
        )
    }
}

/// Cross product of two vectors with the `%` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     Vec3   ( 1.0, 0.0, 0.0)
///     % Vec3 ( 1.0, 0.0, 0.0),
/// // -------------------------------
///     Vec3   ( 0.0, 0.0, 0.0)
/// );
/// assert_eq!(
///     Vec3   ( 1.0, 0.0, 0.0)
///     % Vec3 ( 0.0, 1.0, 0.0),
/// // -------------------------------
///     Vec3   ( 0.0, 0.0, 1.0)
/// );
/// assert_eq!(
///     Vec3   (  3.0, -3.0,  1.0)
///     % Vec3 (  4.0,  9.0,  2.0),
/// // -------------------------------
///     Vec3   (-15.0, -2.0, 39.0)
/// );
/// ```
impl Rem<Vec3> for Vec3 {
    type Output = Vec3;

    fn rem(self, other:Vec3) -> Vec3 {
        Vec3 (
            (self.1 * other.2) - (self.2 * other.1),
            (self.2 * other.0) - (self.0 * other.2),
            (self.0 * other.1) - (self.1 * other.0),
        )
    }
}

/// Cross product assignment of two vectors with the `%=` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// let mut a = Vec3 (1.0, 0.0, 0.0);
/// let b = Vec3 (1.0, 0.0, 0.0);
/// a %= b;
/// assert_eq!(a, Vec3 (0.0, 0.0, 0.0));
///
/// let mut a = Vec3 (1.0, 0.0, 0.0);
/// let b = Vec3 (0.0, 1.0, 0.0);
/// a %= b;
/// assert_eq!(a, Vec3 (0.0, 0.0, 1.0));
///
/// let mut a = Vec3 (  3.0, -3.0,  1.0);
/// let b = Vec3 (  4.0,  9.0,  2.0);
/// a %= b;
/// assert_eq!(a, Vec3 (-15.0, -2.0, 39.0));
/// ```
impl RemAssign for Vec3 {
    fn rem_assign(&mut self, other:Vec3) {
        let x = (self.1 * other.2) - (self.2 * other.1);
        let y = (self.2 * other.0) - (self.0 * other.2);
        let z = (self.0 * other.1) - (self.1 * other.0);
        self.0 = x;
        self.1 = y;
        self.2 = z;
    }
}

/// Dot product of two vectors with the `*` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     Vec3   ( 1.0, 0.0, 0.0)
///     * Vec3 ( 1.0, 0.0, 0.0),
/// // -------------------------------
///     1.0
/// );
/// assert_eq!(
///     Vec3   ( 1.0, 0.0, 0.0)
///     * Vec3 ( 0.0, 1.0, 0.0),
/// // -------------------------------
///     0.0
/// );
/// assert_eq!(
///     Vec3   ( 3.0,  -3.0, 1.0)
///     * Vec3 ( 4.0,   9.0, 2.0),
/// // -------------------------------
///     -13.0
/// );
/// ```
impl Mul<Vec3> for Vec3 {
    type Output = f32;

    fn mul(self, other:Vec3) -> f32 {
        self.0 * other.0
        + self.1 * other.1
        + self.2 * other.2
    }
}

/// Division of a vector by a scalar with the `/` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     Vec3 (4.0, 3.9,   -804.23) / 4.0,
///     Vec3 (1.0, 0.975, -201.0575)
/// );
/// ```
impl Div<f32> for Vec3 {
    type Output = Vec3;

    fn div(self, other:f32) -> Vec3 {
        Vec3 (
            self.0 / other,
            self.1 / other,
            self.2 / other,
        )
    }
}

/// Division assignment of a vector by a scalar
/// with the `/=` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// let mut a = Vec3 (4.0, 3.9, -804.23);
///
/// a /= 4.0;
///
/// assert_eq!(a, Vec3 (1.0, 0.975, -201.0575));
/// ```
impl DivAssign<f32> for Vec3 {
    fn div_assign(&mut self, other:f32) {
        self.0 /= other;
        self.1 /= other;
        self.2 /= other;
    }
}

/// Division of a scalar by a vector with the `/` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     4.0 / Vec3 (4.0, 16.0, -0.0025),
///           Vec3 (1.0, 0.25, -1600.0)
/// );
impl Div<Vec3> for f32 {
    type Output = Vec3;

    fn div(self, other:Vec3) -> Vec3 {
        Vec3 (
            self / other.0,
            self / other.1,
            self / other.2,
        )
    }
}

/// Negation of a vector with the `-` operator.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
///
/// assert_eq!(
///     - Vec3 ( 4.0, 0.0, -0.0025),
///       Vec3 (-4.0, 0.0,  0.0025)
/// );
impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 (
            -self.0,
            -self.1,
            -self.2,
        )
    }
}

impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "<{}, {}, {}>", self.0, self.1, self.2)
    }
}

/// Indexing with the `VecIdx` enum. Just in case it's useful.
///
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
/// use noether::geom::vec3::VecIdx;
///
/// let a_vec = Vec3 ( 4.0, 0.0, -0.0025);
///
/// assert_eq!(a_vec[VecIdx::X], a_vec.0);
/// assert_eq!(a_vec[VecIdx::Y], a_vec.1);
/// assert_eq!(a_vec[VecIdx::Z], a_vec.2);
/// assert_eq!(a_vec[VecIdx::X],  4.0);
/// assert_eq!(a_vec[VecIdx::Y],  0.0);
/// assert_eq!(a_vec[VecIdx::Z], -0.0025);
/// ```
impl Index<VecIdx> for Vec3 {
    type Output = f32;
    fn index(&self, index: VecIdx) -> &f32 {
        match index {
            VecIdx::X => &self.0,
            VecIdx::Y => &self.1,
            VecIdx::Z => &self.2,
        }
    }
}

/// Index enum for accessing Vec3 with indexing syntax
/// # Examples
///
/// ```
/// use noether::geom::Vec3;
/// use noether::geom::vec3::VecIdx;
///
/// let a_vec = Vec3 ( 4.0, 0.0, -0.0025);
///
/// assert_eq!(a_vec[VecIdx::X], a_vec.0);
/// assert_eq!(a_vec[VecIdx::Y], a_vec.1);
/// assert_eq!(a_vec[VecIdx::Z], a_vec.2);
/// ```
///
/// ```
/// use noether::geom::Vec3;
/// use noether::geom::vec3::VecIdx::*;
///
/// let a_vec = Vec3 ( 4.0, 0.0, -0.0025);
///
/// assert_eq!(a_vec[X],  4.0);
/// assert_eq!(a_vec[Y],  0.0);
/// assert_eq!(a_vec[Z], -0.0025);
/// ```
pub enum VecIdx { X, Y, Z }
