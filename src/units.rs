#![allow(clippy::excessive_precision, clippy::eq_op)]


system! {
    quantities: Q {
        length: nanometer, L;
        mass: dalton, M;
        time: picosecond, Ti;
        charge: elementary_charge, Q;
        temperature: kelvin, Te;
    }

    units: U {
        length::Length,
        mass::Mass,
        time::Time,
        charge::Charge,
        temperature::Temperature,

        energy::Energy,
        force::Force,
        entropy::Entropy,
        velocity::Velocity,
        area::Area,
    }
}

storage_types! {
    /// Type modules.
    pub types: Float; // Create storage types f32 and f64

    Q!(crate::units, V);

    use std::marker::PhantomData;

    pub const NM: Length = Length {value: 1.0, dimension: PhantomData, units: PhantomData};
    pub const DA: Mass = Mass {value: 1.0, dimension: PhantomData, units: PhantomData};
    pub const PS: Time = Time {value: 1.0, dimension: PhantomData, units: PhantomData};
    pub const E: Charge = Charge {value: 1.0, dimension: PhantomData, units: PhantomData};
    pub const K: Temperature = Temperature {value: 1.0, dimension: PhantomData, units: PhantomData};

    pub const KJPERMOL: Energy = Energy {value: 1.0, dimension: PhantomData, units: PhantomData};
    pub const KNPERMOL: Force = Force {value: 1.0, dimension: PhantomData, units: PhantomData};

    // Following values are truncated in f32
    pub const AVOGADROS_NUMBER: V = 6.022_140_76_E23;

    pub const BOLTZMANN_CONSTANT: Entropy = Entropy {value: 8.314_462_1_E-3, dimension: PhantomData, units: PhantomData};
}


#[cfg(test)]
mod unit_ratio_tests {
    use super::f64::*;

    #[test]
    fn check_length_ratios() {
        use super::length::{
            nanometer,
            angstrom,
            meter
        };

        assert_lt!((1.0 * NM - Length::new::<nanometer>(1.0)).abs(), 1e-200 * NM);
        assert_lt!((1.0 * NM - Length::new::<angstrom>(10.0)).abs(), 1e-15 * NM);
        assert_lt!((1.0 * NM - Length::new::<meter>(1e-9)).abs(), 1e-15 * NM);
    }

    #[test]
    fn check_mass_ratios() {
        use super::mass::{
            dalton,
            amu,
            gram,
            kilogram
        };

        assert_lt!((1.0 * DA - Mass::new::<dalton>(1.0)).abs(), 1e-200 * DA);
        assert_lt!((1.0 * DA - Mass::new::<amu>(1.0)).abs(), 1e-200 * DA);
        assert_lt!((1.0 * DA - (Mass::new::<gram>(1.0) / AVOGADROS_NUMBER)).abs(), 1e-200 * DA);
        assert_lt!((1.0 * DA - Mass::new::<kilogram>( 1.660_539_066_60_E-27)).abs(), 1e-9 * DA);

    }

    #[test]
    fn check_time_ratios() {
        use super::time::{
            picosecond,
            second,
            minute,
            hour,
            day,
            year
        };

        assert_lt!((1.0 * PS - Time::new::<picosecond>(1.0)).abs(), 1e-200 * PS);
        assert_lt!((1.0 * PS - Time::new::<second>(1e-12)).abs(), 1e-15 * PS);
        assert_lt!((Time::new::<second>(60.0) - Time::new::<minute>(1.0)).abs(), 1e-15 * PS);
        assert_lt!((Time::new::<minute>(60.0) - Time::new::<hour>(1.0)).abs(), 1e-15 * PS);
        assert_lt!((Time::new::<hour>(24.0) - Time::new::<day>(1.0)).abs(), 1e-15 * PS);
        assert_lt!((Time::new::<day>(365.0) - Time::new::<year>(1.0)).abs(), 1e-15 * PS);
    }

    #[test]
    fn check_charge_ratios() {
        use super::charge::{
            elementary_charge,
            coulomb
        };

        assert_lt!((1.0 * E - Charge::new::<elementary_charge>(1.0)).abs(), 1e-200 * E);
        assert_lt!((1.0 * E - Charge::new::<coulomb>(1.602_176_62_E-19)).abs(), 1e-200 * E);
    }

    #[test]
    fn check_temperature_ratios() {
        use super::temperature::{
            kelvin,
            microkelvin
        };

        assert_lt!((1.0 * K - Temperature::new::<kelvin>(1.0)).abs(), 1e-200 * K);
        assert_lt!((1.0 * K - Temperature::new::<microkelvin>(1.0E6)).abs(), 1e-200 * K);
    }

    #[test]
    fn check_energy_ratios() {
        use super::energy::{
            kilojoule_per_mole,
            kilojoule,
            joule_per_mole,
            kilocalorie_per_mole
        };

        assert_lt!((1.0 * KJPERMOL - Energy::new::<kilojoule_per_mole>(1.0)).abs(), 1e-200 * KJPERMOL);
        assert_lt!((AVOGADROS_NUMBER * KJPERMOL - Energy::new::<kilojoule>(1.0)).abs(), 1e-200 * KJPERMOL);
        assert_lt!((1.0 * KJPERMOL - Energy::new::<joule_per_mole>(1000.0)).abs(), 1e-200 * KJPERMOL);
        assert_lt!((4184.0 * KJPERMOL - Energy::new::<kilocalorie_per_mole>(1.0)).abs(), 1e-200 * KJPERMOL);
    }

    #[test]
    fn check_force_ratios() {
        use super::force::{
            kilonewton_per_mole,
            kilonewton,
            newton_per_mole
        };

        assert_lt!((1.0 * KNPERMOL - Force::new::<kilonewton_per_mole>(1.0)).abs(), 1e-200 * KNPERMOL);
        assert_lt!((AVOGADROS_NUMBER * KNPERMOL - Force::new::<kilonewton>(1.0)).abs(), 1e-200 * KNPERMOL);
        assert_lt!((1.0 * KNPERMOL - Force::new::<newton_per_mole>(1000.0)).abs(), 1e-200 * KNPERMOL);
    }
}
