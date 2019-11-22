quantity! {
    /// Temperature (base unit kelvin, K).
    quantity: Temperature; "temperature";
    /// Temperature dimension, K.
    dimension: Q<
        Z0,  // length
        Z0,  // mass
        Z0,  // time
        Z0,  // charge
        P1>; // temperature
    units {
        @yottakelvin: prefix!(yotta); "YK", "yottakelvin", "yottakelvins";
        @zettakelvin: prefix!(zetta); "ZK", "zettakelvin", "zettakelvins";
        @exakelvin: prefix!(exa); "EK", "exakelvin", "exakelvins";
        @petakelvin: prefix!(peta); "PK", "petakelvin", "petakelvins";
        @terakelvin: prefix!(tera); "TK", "terakelvin", "terakelvins";
        @gigakelvin: prefix!(giga); "GK", "gigakelvin", "gigakelvins";
        @megakelvin: prefix!(mega); "MK", "megakelvin", "megakelvins";
        @kilokelvin: prefix!(kilo); "kK", "kilokelvin", "kilokelvins";
        @hectokelvin: prefix!(hecto); "hK", "hectokelvin", "hectokelvins";
        @decakelvin: prefix!(deca); "daK", "decakelvin", "decakelvins";
        /// The kelvin is the SI unit of thermodynamic temperature. It is defined by taking the
        /// fixed numerical value of the Boltzmann constant *k* to be 1.380 649 × 10⁻²³ when
        /// expressed in the unit J K⁻¹, which is equal to kg m² s⁻² K⁻¹, where the kilogram, meter,
        /// and second are defined in terms of *h*, *c*, and ∆*ν*<sub>Cs</sub>.
        @kelvin: prefix!(none); "K", "kelvin", "kelvins";
        @decikelvin: prefix!(deci); "dK", "decikelvin", "decikelvins";
        @centikelvin: prefix!(centi); "cK", "centikelvin", "centikelvins";
        @millikelvin: prefix!(milli); "mK", "millikelvin", "millikelvins";
        @microkelvin: prefix!(micro); "µK", "microkelvin", "microkelvins";
        @nanokelvin: prefix!(nano); "nK", "nanokelvin", "nanokelvins";
        @picokelvin: prefix!(pico); "pK", "picokelvin", "picokelvins";
        @femtokelvin: prefix!(femto); "fK", "femtokelvin", "femtokelvins";
        @attokelvin: prefix!(atto); "aK", "attokelvin", "attokelvins";
        @zeptokelvin: prefix!(zepto); "zK", "zeptokelvin", "zeptokelvins";
        @yoctokelvin: prefix!(yocto); "yK", "yoctokelvin", "yoctokelvins";
    }
}
