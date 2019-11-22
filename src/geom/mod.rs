// pub mod vec3;
pub mod vec3d;

use crate::units;

pub use self::vec3d::Vec3D;

type V = f32;

/// 32 bit position in GMX units
pub type MDPos = units::Nanometer<V>;
/// 32 bit force in GMX units
pub type MDForce = units::KilojoulePerMolePerNanometer<V>;
/// 32 bit velocity in GMX units
pub type MDVeloc = units::NanometerPerPicosecond<V>;
/// 32 bit dimensionless quantity in GMX units
pub type MDNodim = units::Unitless<V>;

/// Vector holding 32 bit positions
pub type PosVec = Vec3D<MDPos>;
/// Vector holding 32 bit forces
pub type ForceVec = Vec3D<MDForce>;
/// Vector holding 32 bit velocities
pub type VelocVec = Vec3D<MDVeloc>;
/// Vector holding 32 bit dimensionless quantities
pub type NodimVec = Vec3D<MDNodim>;
