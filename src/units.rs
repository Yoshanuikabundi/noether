/// Length (base unit nanometer, nm<sup>1</sup>).
#[macro_use]
pub mod length {
    quantity! {
        quantity: Length; "length";
        /// Length dimension, nm<sup>1</sup>.
        dimension: GMXQ<
            P1,  // length
            Z0,  // mass
            Z0,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Base length unit, equal to 1 × 10<sup>-9</sup> m
            @nanometer: prefix!(nano) / prefix!(nano); "nm", "nanometer", "nanometers";
            /// Length unit, equal to 0.1 nm (approximately the length of a bond or VdW radius of an atom)
            @angstrom: 1.0E-1; "Å", "ångström", "ångströms";
        }
    }
}

/// Mass (base unit Dalton, Da<sup>1</sup>).
#[macro_use]
pub mod mass {
    quantity! {
        quantity: Mass; "mass";
        /// Mass dimension, Da<sup>1</sup>.
        dimension: GMXQ<
            Z0,  // length
            P1,  // mass
            Z0,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Base mass unit, approximately 1.660 538 921 × 10<sup>−27</sup> kg (numerically 1 g mol<sup>-1</sup>)
            @dalton: prefix!(none); "Da", "Dalton", "Daltons";
            @gram_per_mole: prefix!(none); "g/mol", "gram per mole", "grams per mole";
            /// Mass unit, approximately 1.660 538 921 × 10<sup>−24</sup> kg
            @kilodalton: prefix!(kilo); "kDa", "kilodalton", "kilodaltons";
            /// Mass unit, approximately 1.660 538 921 × 10<sup>−21</sup> kg
            @megadalton: prefix!(mega); "MDa", "megadalton", "megadaltons";
            /// SI mass unit, approximately 6.022 140 857 × 10<sup>27</sup> Da
            @kilogram: 6.022_140_857E27; "kg", "kilogram", "kilograms";
        }
    }
}

/// Time (base unit picosecond, ps<sup>1</sup>).
#[macro_use]
pub mod time {
    quantity! {
        quantity: Time; "time";
        /// Time dimension, ps<sup>1</sup>.
        dimension: GMXQ<
            Z0,  // length
            Z0,  // mass
            P1,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Base time unit, equal to 1 × 10<sup>-12</sup> s
            @picosecond: prefix!(pico) / prefix!(pico); "ps", "picosecond", "picoseconds";
            /// SI time unit, equal to 1 × 10<sup>12</sup> ps
            @second: prefix!(none) / prefix!(pico); "s", "second", "seconds";
            /// SI time unit, equal to 1 × 10<sup>-3</sup> ps
            @femtosecond: prefix!(femto) / prefix!(pico); "fs", "femtosecond", "femtoseconds";
            /// SI time unit, equal to 1 × 10<sup>3</sup> ps
            @nanosecond: prefix!(nano) / prefix!(pico); "ns", "nanosecond", "nanoseconds";
            /// SI time unit, equal to 1 × 10<sup>6</sup> ps
            @microsecond: prefix!(micro) / prefix!(pico); "µs", "microsecond", "microseconds";
            /// SI time unit, equal to 1 × 10<sup>9</sup> ps
            @millisecond: prefix!(milli) / prefix!(pico); "ms", "millisecond", "milliseconds";
        }
    }
}

/// Charge (base unit elementary charge, *e*<sup>1</sup>).
#[macro_use]
pub mod charge {
    quantity! {
        quantity: Charge; "charge";
        /// Charge dimension, *e*<sup>1</sup>.
        dimension: GMXQ<
            Z0,  // length
            Z0,  // mass
            Z0,  // time
            P1,  // charge
            Z0>; // temperature
        units {
            // Base charge unit, charge of an electron (or proton)
            @elem_charge: 1.0; "e", "elementary charge", "elementary charges";
        }
    }
}

/// Temperature (base unit Kelvin, K<sup>1</sup>).
#[macro_use]
pub mod temperature {
    quantity! {
        quantity: Temperature; "temperature";
        /// Temperature dimension, K<sup>1</sup>.
        dimension: GMXQ<
            Z0,  // length
            Z0,  // mass
            Z0,  // time
            Z0,  // charge
            P1>; // temperature
        units {
            /// Base temperature unit, absolute
            @kelvin: 1.0; "K", "Kelvin", "Kelvins";
        }
    }
}

/// Energy (base unit kilojoule per mole, kJ<sup>1</sup>mol<sup>-1</sup>).
#[macro_use]
pub mod energy {
    quantity! {
        quantity: Energy; "energy";
        /// Energy dimension
        dimension: GMXQ<
            P2,  // length
            P1,  // mass
            N2,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Derived energy unit
            @kilojoule_per_mole: prefix!(kilo) / prefix!(kilo); "kJ/mol", "kilojoule per mole", "kilojoules per mole";
        }
    }
}

/// Force (base unit kilojoule per mole per nanometer, kJ<sup>1</sup>mol<sup>-1</sup>nm<sup>-1</sup>).
#[macro_use]
pub mod force {
    quantity! {
        quantity: Force; "force";
        /// Force dimension
        dimension: GMXQ<
            P1,  // length
            P1,  // mass
            N2,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Derived force unit
            @kilojoule_per_mole_per_nanometer: prefix!(kilo) / prefix!(kilo); "kJ/mol/nm", "kilojoule per mole per nanometer", "kilojoules per mole per nanometer";
            @teranewton_per_mole: prefix!(tera) / prefix!(tera); "TN/mol", "teranewton per mole", "teranewtons per mole";
        }
    }
}

/// Pressure (Base unit kilojoule per mol per nanometer cubed, kJ<sup>1</sup>mol<sup>-1</sup>nm<sup>-3</sup>).
#[macro_use]
pub mod pressure {
    quantity! {
        quantity: Pressure; "pressure";
        /// Pressure dimension
        dimension: GMXQ<
            N1,  // length
            P1,  // mass
            N2,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Derived pressure unit
            @kilojoule_per_mol_per_nanometer_cubed: 1.0; "kJ/mol/nm^3", "kilojoule per mole per cubic nanometer", "kilojoules per mole per cubic nanometer";
            @bar: 16.605_390_404; "bar", "bar", "bar";
        }
    }
}

/// Velocity (base unit nanometer per picosecond, nm<sup>1</sup>ps<sup>-1</sup>, 1000 m<sup>1</sup>s<sup>-1</sup>).
#[macro_use]
pub mod velocity {
    quantity! {
        quantity: Velocity; "velocity";
        /// Velocity dimension
        dimension: GMXQ<
            P1,  // length
            Z0,  // mass
            N1,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Derived velocity unit
            @nanometer_per_picosecond: 1.0; "nm/ps", "nanometer per picosecond", "nanometers per picosecond";
            @kilometer_per_second: prefix!(kilo) / prefix!(kilo); "km/s", "kilometer per second", "kilometers per second";
        }
    }
}

/// Dipole moment (base unit elementary charge nanometer, *e*<sup>1</sup>nm<sup>1</sup>).
#[macro_use]
pub mod dipole_moment {
    quantity! {
        quantity: DipoleMoment; "dipole moment";
        /// Dipole moment dimension
        dimension: GMXQ<
            P1,  // length
            Z0,  // mass
            Z0,  // time
            P1,  // charge
            Z0>; // temperature
        units {
            /// Derived dipole moment unit
            @elem_charge_nanometer: 1.0; "e ns", "elementary charge nanometer", "elementary charge nanometers";
        }
    }
}

/// Electric potential (base unit kilojoule per mole per elementary charge, kJ<sup>1</sup>mol<sup>-1</sup>*e*<sup>-1</sup>).
#[macro_use]
pub mod electric_potential {
    quantity! {
        quantity: ElectricPotential; "electric potential";
        /// Electric potential dimension
        dimension: GMXQ<
            P2,  // length
            P1,  // mass
            N2,  // time
            N1,  // charge
            Z0>; // temperature
        units {
            /// Derived electric potential unit
            @kilojoule_per_mole_per_elem_charge: 1.0; "kJ/mol/e", "kilojoule per mole per elementary charge", "kilojoules per mole per elementary charge";
        }
    }
}

/// Electric field (base unit kilojoule per mole per nanometer per elementary charge, kJ<sup>1</sup>mol<sup>-1</sup>nm<sup>-1</sup>*e*<sup>-1</sup>).
#[macro_use]
pub mod electric_field {
    quantity! {
        quantity: ElectricField; "electric field";
        /// Electric field dimension
        dimension: GMXQ<
            P1,  // length
            P1,  // mass
            N2,  // time
            N1,  // charge
            Z0>; // temperature
        units {
            /// Derived electric field unit
            @kilojoule_per_mole_per_nanometer_per_elem_charge: 1.0; "kJ/mol/nm/e", "kilojoule per mole per nanometer per elementary charge", "kilojoules per mole per nanometer per elementary charge";
        }
    }
}

/// Area (base unit nanometer, nm<sup>2</sup>).
#[macro_use]
pub mod area {
    quantity! {
        quantity: Area; "area";
        /// area dimension, nm<sup>2</sup>.
        dimension: GMXQ<
            P2,  // length
            Z0,  // mass
            Z0,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Base area unit, equal to 1 × 10<sup>-18</sup> m<sup>2</sup>
            @square_nanometer: 1.0; "nm^2", "square nanometer", "square nanometers";
            /// area unit
            @square_angstrom: 1.0E-2; "Å^2", "square ångström", "square ångströms";
        }
    }
}

/// Volume (base unit nanometer, nm<sup>3</sup>).
#[macro_use]
pub mod volume {
    quantity! {
        quantity: Volume; "volume";
        /// Volume dimension, nm<sup>3</sup>.
        dimension: GMXQ<
            P3,  // length
            Z0,  // mass
            Z0,  // time
            Z0,  // charge
            Z0>; // temperature
        units {
            /// Base volume unit, equal to 1 × 10<sup>-27</sup> m<sup>3</sup>
            @cubic_nanometer: 1.0; "nm^3", "cubic nanometer", "cubic nanometers";
            /// Volume unit
            @cubic_angstrom: 1.0E-3; "Å^3", "cubic ångström", "cubic ångströms";
        }
    }
}

#[macro_use]
system! {
    quantities: GMXQ {
        length: nanometer, NM;
        mass: dalton, DA;
        time: picosecond, PS;
        charge: elem_charge, E;
        temperature: kelvin, K;
    }
    units: GMXU {
        mod length::Length,
        mod mass::Mass,
        mod time::Time,
        mod charge::Charge,
        mod temperature::Temperature,
        mod energy::Energy,
        mod force::Force,
        mod pressure::Pressure,
        mod velocity::Velocity,
        mod dipole_moment::DipoleMoment,
        mod electric_potential::ElectricPotential,
        mod electric_field::ElectricField,
        mod area::Area,
        mod volume::Volume,
    }
}

storage_types! {
    pub types: Float;

    use std::marker::PhantomData;

    GMXQ!(units, V);

    /// Create a `const` length of 1 nm.
    pub const NM: Length = Length { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` mass of 1 Da.
    pub const DA: Mass = Mass { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` time of 1 ps.
    pub const PS: Time = Time { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` charge of 1 e.
    pub const E: Charge = Charge { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` temperature of 1 K.
    pub const K: Temperature = Temperature { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` energy of 1 kJ/mol.
    pub const KJPM: Energy = Energy { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` force of 1 kJ/mol.
    pub const KJPMNM: Force = Force { dimension: PhantomData, units: PhantomData, value: 1.0, };

    /// Create a `const` velocity of 1 nm/ps.
    pub const MPS: Velocity = Velocity { dimension: PhantomData, units: PhantomData, value: 1.0, };
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
