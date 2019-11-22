quantity! {
    /// Length (base unit nanometer, nm).
    quantity: Length; "length";
    /// Length dimension, nm.
    dimension: Q<
        P1,  // length
        Z0,  // mass
        Z0,  // time
        Z0,  // charge
        Z0>; // temperature
    units {
        @nanometer: prefix!(nano) / prefix!(nano); "nm", "nanometer", "nanometers";

        @yottameter: prefix!(yotta) / prefix!(nano); "Ym", "yottameter", "yottameters";
        @zettameter: prefix!(zetta) / prefix!(nano); "Zm", "zettameter", "zettameters";
        @exameter: prefix!(exa) / prefix!(nano); "Em", "exameter", "exameters";
        @petameter: prefix!(peta) / prefix!(nano); "Pm", "petameter", "petameters";
        @terameter: prefix!(tera) / prefix!(nano); "Tm", "terameter", "terameters";
        @gigameter: prefix!(giga) / prefix!(nano); "Gm", "gigameter", "gigameters";
        @megameter: prefix!(mega) / prefix!(nano); "Mm", "megameter", "megameters";
        @kilometer: prefix!(kilo) / prefix!(nano); "km", "kilometer", "kilometers";
        @hectometer: prefix!(hecto) / prefix!(nano); "hm", "hectometer", "hectometers";
        @decameter: prefix!(deca) / prefix!(nano); "dam", "decameter", "decameters";
        /// The meter is the SI unit of length. It is defined by taking the fixed numerical value
        /// of the speed of light in vacuum *c* to be 299 792 458 when expressed in the unit m s⁻¹,
        /// where the second is defined in terms of the caesium frequency ∆*ν*<sub>Cs</sub>.
        @meter: prefix!(none) / prefix!(nano); "m", "meter", "meters";
        @decimeter: prefix!(deci) / prefix!(nano); "dm", "decimeter", "decimeters";
        @centimeter: prefix!(centi) / prefix!(nano); "cm", "centimeter", "centimeters";
        @millimeter: prefix!(milli) / prefix!(nano); "mm", "millimeter", "millimeters";
        @micrometer: prefix!(micro) / prefix!(nano); "µm", "micrometer", "micrometers";
        @picometer: prefix!(pico) / prefix!(nano); "pm", "picometer", "picometers";
        @femtometer: prefix!(femto) / prefix!(nano); "fm", "femtometer", "femtometers";
        @attometer: prefix!(atto) / prefix!(nano); "am", "attometer", "attometers";
        @zeptometer: prefix!(zepto) / prefix!(nano); "zm", "zeptometer", "zeptometers";
        @yoctometer: prefix!(yocto) / prefix!(nano); "ym", "yoctometer", "yoctometers";

        @angstrom: 1.0_E-10 / prefix!(nano); "Å", "Ångström", "Ångströms";
        @astronomical_unit: 1.495_979_E11 / prefix!(nano); "ua", "astronomical unit", "astronomical units";
        @chain: 2.011_684_E1 / prefix!(nano); "ch", "chain", "chains";
        @fathom: 1.828_804_E0 / prefix!(nano); "fathom", "fathom", "fathoms";
        @fermi: 1.0_E-15 / prefix!(nano); "fermi", "fermi", "fermis";
        @foot: 3.048_E-1 / prefix!(nano); "ft", "foot", "feet";
        @foot_survey: 3.048_006_E-1 / prefix!(nano); "ft (U.S. survey)", "foot (U.S. survey)", "feet (U.S. survey)";
        @inch: 2.54_E-2 / prefix!(nano); "in", "inch", "inches";
        @light_year: 9.460_73_E15 / prefix!(nano); "l. y.", "light year", "light years";
        @microinch: 2.54_E-8 / prefix!(nano); "μin", "microinch", "microinches";
        @micron: 1.0_E-6 / prefix!(nano); "μ", "micron", "microns";
        @mil: 2.54_E-5 / prefix!(nano); "0.001 in", "mil", "mils";
        @mile: 1.609_344_E3 / prefix!(nano); "mi", "mile", "miles";
        @mile_survey: 1.609_347_E3 / prefix!(nano); "mi (U.S. survey)", "mile (U.S. survey)", "miles (U.S. survey)";
        @nautical_mile: 1.852_E3 / prefix!(nano); "M", "nautical mile", "nautical miles";
        @parsec: 3.085_678_E16 / prefix!(nano); "pc", "parsec", "parsecs";
        @pica_computer: 4.233_333_333_333_333_E-3 / prefix!(nano); "1/6 in (computer)", "pica (computer)",
            "picas (computer)";
        @pica_printers: 4.217_518_E-3 / prefix!(nano); "1/6 in", "pica (printer's)", "picas (printer's)";
        @point_computer: 3.527_778_E-4 / prefix!(nano); "1/72 in (computer)", "point (computer)",
            "points (computer)";
        @point_printers: 3.514_598_E-4 / prefix!(nano); "1/72 in", "point (printer's)", "points (printer's)";
        @rod: 5.029_21_E0 / prefix!(nano); "rd", "rod", "rods";
        @yard: 9.144_E-1 / prefix!(nano); "yd", "yard", "yards";

    }
}
