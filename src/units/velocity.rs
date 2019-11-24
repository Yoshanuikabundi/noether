quantity! {
    /// Velocity (base unit nanometer per picosecond, nm/ps).
    quantity: Velocity; "velocity";
    /// Velocity dimension
    dimension: Q<
        P1,  // length
        Z0,  // mass
        N1,  // time
        Z0,  // charge
        Z0>; // temperature
    units {
        @nanometer_per_picosecond: prefix!(nano) / prefix!(pico) / (prefix!(nano) / prefix!(pico)); "nm", "nanometer", "nanometers";

        @yottameter_per_picosecond: prefix!(yotta) / (prefix!(nano) / prefix!(pico)); "Ym", "yottameter", "yottameters";
        @zettameter_per_picosecond: prefix!(zetta) / (prefix!(nano) / prefix!(pico)); "Zm", "zettameter", "zettameters";
        @exameter_per_picosecond: prefix!(exa) / (prefix!(nano) / prefix!(pico)); "Em", "exameter", "exameters";
        @petameter_per_picosecond: prefix!(peta) / (prefix!(nano) / prefix!(pico)); "Pm", "petameter", "petameters";
        @terameter_per_picosecond: prefix!(tera) / (prefix!(nano) / prefix!(pico)); "Tm", "terameter", "terameters";
        @gigameter_per_picosecond: prefix!(giga) / (prefix!(nano) / prefix!(pico)); "Gm", "gigameter", "gigameters";
        @megameter_per_picosecond: prefix!(mega) / (prefix!(nano) / prefix!(pico)); "Mm", "megameter", "megameters";
        @kilometer_per_picosecond: prefix!(kilo) / (prefix!(nano) / prefix!(pico)); "km", "kilometer", "kilometers";
        @hectometer_per_picosecond: prefix!(hecto) / (prefix!(nano) / prefix!(pico)); "hm", "hectometer", "hectometers";
        @decameter_per_picosecond: prefix!(deca) / (prefix!(nano) / prefix!(pico)); "dam", "decameter", "decameters";
        /// The meter is the SI unit of length. It is defined by taking the fixed numerical value
        /// of the speed of light in vacuum *c* to be 299 792 458 when expressed in the unit m s⁻¹,
        /// where the second is defined in terms of the caesium frequency ∆*ν*<sub>Cs</sub>.
        @meter_per_picosecond: prefix!(none) / (prefix!(nano) / prefix!(pico)); "m", "meter", "meters";
        @decimeter_per_picosecond: prefix!(deci) / (prefix!(nano) / prefix!(pico)); "dm", "decimeter", "decimeters";
        @centimeter_per_picosecond: prefix!(centi) / (prefix!(nano) / prefix!(pico)); "cm", "centimeter", "centimeters";
        @millimeter_per_picosecond: prefix!(milli) / (prefix!(nano) / prefix!(pico)); "mm", "millimeter", "millimeters";
        @micrometer_per_picosecond: prefix!(micro) / (prefix!(nano) / prefix!(pico)); "µm", "micrometer", "micrometers";
        @picometer_per_picosecond: prefix!(pico) / (prefix!(nano) / prefix!(pico)); "pm", "picometer", "picometers";
        @femtometer_per_picosecond: prefix!(femto) / (prefix!(nano) / prefix!(pico)); "fm", "femtometer", "femtometers";
        @attometer_per_picosecond: prefix!(atto) / (prefix!(nano) / prefix!(pico)); "am", "attometer", "attometers";
        @zeptometer_per_picosecond: prefix!(zepto) / (prefix!(nano) / prefix!(pico)); "zm", "zeptometer", "zeptometers";
        @yoctometer_per_picosecond: prefix!(yocto) / (prefix!(nano) / prefix!(pico)); "ym", "yoctometer", "yoctometers";

        @angstrom_per_picosecond: 1.0_E-10 / (prefix!(nano) / prefix!(pico)); "Å", "Ångström", "Ångströms";

    }
}
