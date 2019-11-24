quantity! {
    /// Force (base unit kilonewton per mole, kg · m · s⁻² · mol⁻1).
    quantity: Force; "force";
    /// Dimension of Force
    dimension: Q<
        P1,  // length
        P1,  // mass
        N2,  // time
        Z0,  // charge
        Z0>; // temperature
    units {
        @kilonewton_per_mole: 1.0; "KN/mol", "kilonewton per mole", "kilonewtons per mole";

        @yottanewton_per_mole: prefix!(yotta) / prefix!(kilo); "YN", "yottanewton", "yottanewtons";
        @zettanewton_per_mole: prefix!(zetta) / prefix!(kilo); "ZN", "zettanewton", "zettanewtons";
        @exanewton_per_mole: prefix!(exa) / prefix!(kilo); "EN", "exanewton", "exanewtons";
        @petanewton_per_mole: prefix!(peta) / prefix!(kilo); "PN", "petanewton", "petanewtons";
        @teranewton_per_mole: prefix!(tera) / prefix!(kilo); "TN", "teranewton", "teranewtons";
        @giganewton_per_mole: prefix!(giga) / prefix!(kilo); "GN", "giganewton", "giganewtons";
        @meganewton_per_mole: prefix!(mega) / prefix!(kilo); "MN", "meganewton", "meganewtons";
        @hectonewton_per_mole: prefix!(hecto) / prefix!(kilo); "hN", "hectonewton", "hectonewtons";
        @decanewton_per_mole: prefix!(deca) / prefix!(kilo); "daN", "decanewton", "decanewtons";
        /// Derived unit of energy.
        @newton_per_mole: prefix!(none) / prefix!(kilo); "N", "newton", "newtons";
        @decinewton_per_mole: prefix!(deci) / prefix!(kilo); "dN", "decinewton", "decinewtons";
        @centinewton_per_mole: prefix!(centi) / prefix!(kilo); "cN", "centinewton", "centinewtons";
        @millinewton_per_mole: prefix!(milli) / prefix!(kilo); "mN", "millnewton", "millnewtons";
        @micronewton_per_mole: prefix!(micro) / prefix!(kilo); "µN", "micronewton", "micronewtons";
        @nanonewton_per_mole: prefix!(nano) / prefix!(kilo); "nN", "nanonewton", "nanonewtons";
        @piconewton_per_mole: prefix!(pico) / prefix!(kilo); "pN", "piconewton", "piconewtons";
        @femtonewton_per_mole: prefix!(femto) / prefix!(kilo); "fN", "femtonewton", "femtonewtons";
        @attonewton_per_mole: prefix!(atto) / prefix!(kilo); "aN", "attonewton", "attonewtons";
        @zeptonewton_per_mole: prefix!(zepto) / prefix!(kilo); "zN", "zeptonewton", "zeptonewtons";
        @yoctonewton_per_mole: prefix!(yocto) / prefix!(kilo); "yN", "yoctonewton", "yoctonewtons";

        @yottanewton: prefix!(yotta) / (prefix!(kilo) / 6.022_140_76_E23); "YN", "yottanewton", "yottanewtons";
        @zettanewton: prefix!(zetta) / (prefix!(kilo) / 6.022_140_76_E23); "ZN", "zettanewton", "zettanewtons";
        @exanewton: prefix!(exa) / (prefix!(kilo) / 6.022_140_76_E23); "EN", "exanewton", "exanewtons";
        @petanewton: prefix!(peta) / (prefix!(kilo) / 6.022_140_76_E23); "PN", "petanewton", "petanewtons";
        @teranewton: prefix!(tera) / (prefix!(kilo) / 6.022_140_76_E23); "TN", "teranewton", "teranewtons";
        @giganewton: prefix!(giga) / (prefix!(kilo) / 6.022_140_76_E23); "GN", "giganewton", "giganewtons";
        @meganewton: prefix!(mega) / (prefix!(kilo) / 6.022_140_76_E23); "MN", "meganewton", "meganewtons";
        @kilonewton: prefix!(kilo) / (prefix!(kilo) / 6.022_140_76_E23); "kN", "kilonewton", "kilonewtons";
        @hectonewton: prefix!(hecto) / (prefix!(kilo) / 6.022_140_76_E23); "hN", "hectonewton", "hectonewtons";
        @decanewton: prefix!(deca) / (prefix!(kilo) / 6.022_140_76_E23); "daN", "decanewton", "decanewtons";
        /// Derived unit of energy.
        @newton: prefix!(none) / (prefix!(kilo) / 6.022_140_76_E23); "N", "newton", "newtons";
        @decinewton: prefix!(deci) / (prefix!(kilo) / 6.022_140_76_E23); "dN", "decinewton", "decinewtons";
        @centinewton: prefix!(centi) / (prefix!(kilo) / 6.022_140_76_E23); "cN", "centinewton", "centinewtons";
        @millinewton: prefix!(milli) / (prefix!(kilo) / 6.022_140_76_E23); "mN", "millnewton", "millnewtons";
        @micronewton: prefix!(micro) / (prefix!(kilo) / 6.022_140_76_E23); "µN", "micronewton", "micronewtons";
        @nanonewton: prefix!(nano) / (prefix!(kilo) / 6.022_140_76_E23); "nN", "nanonewton", "nanonewtons";
        @piconewton: prefix!(pico) / (prefix!(kilo) / 6.022_140_76_E23); "pN", "piconewton", "piconewtons";
        @femtonewton: prefix!(femto) / (prefix!(kilo) / 6.022_140_76_E23); "fN", "femtonewton", "femtonewtons";
        @attonewton: prefix!(atto) / (prefix!(kilo) / 6.022_140_76_E23); "aN", "attonewton", "attonewtons";
        @zeptonewton: prefix!(zepto) / (prefix!(kilo) / 6.022_140_76_E23); "zN", "zeptonewton", "zeptonewtons";
        @yoctonewton: prefix!(yocto) / (prefix!(kilo) / 6.022_140_76_E23); "yN", "yoctonewton", "yoctonewtons";

        @calorie_nanometer: 4.184_E0 / (prefix!(kilo) / 6.022_140_76_E23); "cal/nm", "calorie per nanometer", "calories per nanometer";
        @kilocalorie_nanometer: prefix!(kilo) * 4.184_E0 / (prefix!(kilo) / 6.022_140_76_E23); "kcal/nm", "kilocalorie per nanometer", "kilocalories per nanometer";

        @calorie_per_mole_nanometer: 4.184_E0; "cal/mol/nm", "calorie per mole nanometer", "calories per mole nanometer";
        @kilocalorie_per_mole_nanometer: prefix!(kilo) * 4.184_E0; "kcal/mol/nm", "kilocalorie per mole nanometer", "kilocalories per mole nanometer";

        @electronvolt_per_nanometer: 1.602_177_E-19 / (prefix!(kilo) / 6.022_140_76_E23); "eV/nm", "electronvolt per nanometer", "electronvolts per nanometer";
        @erg_per_nanometer: 1.0_E-7 / (prefix!(kilo) / 6.022_140_76_E23); "erg/nm", "erg per nanometer", "ergs per nanometer";
        @hartree_per_nanometer: 2625.5; "Ha/nm", "hartree per nanometer", "hartrees per nanometer";
    }
}
