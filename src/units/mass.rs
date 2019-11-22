quantity! {
    /// Mass (base unit Dalton, Da).
    quantity: Mass; "mass";
    /// Mass dimension, Da.
    dimension: Q<
        Z0,  // length
        P1,  // mass
        Z0,  // time
        Z0,  // charge
        Z0>; // temperature
    units {
        @dalton: 6.022_140_76_E23 / 6.022_140_76_E23; "Da", "Dalton", "Daltons";
        @amu: 6.022_140_76_E23 / 6.022_140_76_E23; "u", "unified atomic mass unit", "unified atomic mass units";

        @yottagram: prefix!(yotta) / 6.022_140_76_E23; "Yg", "yottagram", "yottagrams";
        @zettagram: prefix!(zetta) / 6.022_140_76_E23; "Zg", "zettagram", "zettagrams";
        @exagram: prefix!(exa) / 6.022_140_76_E23; "Eg", "exagram", "exagrams";
        @petagram: prefix!(peta) / 6.022_140_76_E23; "Pg", "petagram", "petagrams";
        @teragram: prefix!(tera) / 6.022_140_76_E23; "Tg", "teragram", "teragrams";
        @gigagram: prefix!(giga) / 6.022_140_76_E23; "Gg", "gigagram", "gigagrams";
        @megagram: prefix!(mega) / 6.022_140_76_E23; "Mg", "megagram", "megagrams";
        /// The kilogram is the SI unit of mass. It is defined by taking the fixed numerical value
        /// of the Planck constant *h* to be 6.626 070 15 × 10⁻³⁴ when expressed in the unit J s,
        /// which is equal to kg m² s⁻¹, where the meter and the second are defined in terms of *c*
        /// and ∆*ν*<sub>Cs</sub>.
        @kilogram: prefix!(kilo) / 6.022_140_76_E23; "kg", "kilogram", "kilograms";
        @hectogram: prefix!(hecto) / 6.022_140_76_E23; "hg", "hectogram", "hectograms";
        @decagram: prefix!(deca) / 6.022_140_76_E23; "dag", "decagram", "decagrams";
        @gram: prefix!(none) / 6.022_140_76_E23; "g", "gram", "grams";
        @decigram: prefix!(deci) / 6.022_140_76_E23; "dg", "decigram", "decigrams";
        @centigram: prefix!(centi) / 6.022_140_76_E23; "cg", "centigram", "centigrams";
        @milligram: prefix!(milli) / 6.022_140_76_E23; "mg", "milligram", "milligrams";
        @microgram: prefix!(micro) / 6.022_140_76_E23; "µg", "microgram", "micrograms";
        @nanogram: prefix!(nano) / 6.022_140_76_E23; "ng", "nanogram", "nanograms";
        @picogram: prefix!(pico) / 6.022_140_76_E23; "pg", "picogram", "picograms";
        @femtogram: prefix!(femto) / 6.022_140_76_E23; "fg", "femtogram", "femtograms";
        @attogram: prefix!(atto) / 6.022_140_76_E23; "ag", "attogram", "attograms";
        @zeptogram: prefix!(zepto) / 6.022_140_76_E23; "zg", "zeptogram", "zeptograms";
        @yoctogram: prefix!(yocto) / 6.022_140_76_E23; "yg", "yoctogram", "yoctograms";

    }
}
