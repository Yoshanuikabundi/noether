quantity! {
    /// Energy (base unit kilojoule per mole, kg · m² · s⁻² · mol⁻1).
    quantity: Energy; "energy";
    /// Dimension of energy, (base unit joule, kg · m² · s⁻² · mol⁻1).
    dimension: Q<
        P2,  // length
        P1,  // mass
        N2,  // time
        Z0,  // charge
        Z0>; // temperature
    units {
        @kilojoule_per_mole: 1.0; "KJ/mol", "kilojoule per mole", "kilojoules per mole";

        @yottajoule_per_mole: prefix!(yotta) / prefix!(kilo); "YJ", "yottajoule", "yottajoules";
        @zettajoule_per_mole: prefix!(zetta) / prefix!(kilo); "ZJ", "zettajoule", "zettajoules";
        @exajoule_per_mole: prefix!(exa) / prefix!(kilo); "EJ", "exajoule", "exajoules";
        @petajoule_per_mole: prefix!(peta) / prefix!(kilo); "PJ", "petajoule", "petajoules";
        @terajoule_per_mole: prefix!(tera) / prefix!(kilo); "TJ", "terajoule", "terajoules";
        @gigajoule_per_mole: prefix!(giga) / prefix!(kilo); "GJ", "gigajoule", "gigajoules";
        @megajoule_per_mole: prefix!(mega) / prefix!(kilo); "MJ", "megajoule", "megajoules";
        @hectojoule_per_mole: prefix!(hecto) / prefix!(kilo); "hJ", "hectojoule", "hectojoules";
        @decajoule_per_mole: prefix!(deca) / prefix!(kilo); "daJ", "decajoule", "decajoules";
        /// Derived unit of energy.
        @joule_per_mole: prefix!(none) / prefix!(kilo); "J", "joule", "joules";
        @decijoule_per_mole: prefix!(deci) / prefix!(kilo); "dJ", "decijoule", "decijoules";
        @centijoule_per_mole: prefix!(centi) / prefix!(kilo); "cJ", "centijoule", "centijoules";
        @millijoule_per_mole: prefix!(milli) / prefix!(kilo); "mJ", "milljoule", "milljoules";
        @microjoule_per_mole: prefix!(micro) / prefix!(kilo); "µJ", "microjoule", "microjoules";
        @nanojoule_per_mole: prefix!(nano) / prefix!(kilo); "nJ", "nanojoule", "nanojoules";
        @picojoule_per_mole: prefix!(pico) / prefix!(kilo); "pJ", "picojoule", "picojoules";
        @femtojoule_per_mole: prefix!(femto) / prefix!(kilo); "fJ", "femtojoule", "femtojoules";
        @attojoule_per_mole: prefix!(atto) / prefix!(kilo); "aJ", "attojoule", "attojoules";
        @zeptojoule_per_mole: prefix!(zepto) / prefix!(kilo); "zJ", "zeptojoule", "zeptojoules";
        @yoctojoule_per_mole: prefix!(yocto) / prefix!(kilo); "yJ", "yoctojoule", "yoctojoules";

        @yottajoule: prefix!(yotta) / (prefix!(kilo) / 6.022_140_76_E23); "YJ", "yottajoule", "yottajoules";
        @zettajoule: prefix!(zetta) / (prefix!(kilo) / 6.022_140_76_E23); "ZJ", "zettajoule", "zettajoules";
        @exajoule: prefix!(exa) / (prefix!(kilo) / 6.022_140_76_E23); "EJ", "exajoule", "exajoules";
        @petajoule: prefix!(peta) / (prefix!(kilo) / 6.022_140_76_E23); "PJ", "petajoule", "petajoules";
        @terajoule: prefix!(tera) / (prefix!(kilo) / 6.022_140_76_E23); "TJ", "terajoule", "terajoules";
        @gigajoule: prefix!(giga) / (prefix!(kilo) / 6.022_140_76_E23); "GJ", "gigajoule", "gigajoules";
        @megajoule: prefix!(mega) / (prefix!(kilo) / 6.022_140_76_E23); "MJ", "megajoule", "megajoules";
        @kilojoule: prefix!(kilo) / (prefix!(kilo) / 6.022_140_76_E23); "kJ", "kilojoule", "kilojoules";
        @hectojoule: prefix!(hecto) / (prefix!(kilo) / 6.022_140_76_E23); "hJ", "hectojoule", "hectojoules";
        @decajoule: prefix!(deca) / (prefix!(kilo) / 6.022_140_76_E23); "daJ", "decajoule", "decajoules";
        /// Derived unit of energy.
        @joule: prefix!(none) / (prefix!(kilo) / 6.022_140_76_E23); "J", "joule", "joules";
        @decijoule: prefix!(deci) / (prefix!(kilo) / 6.022_140_76_E23); "dJ", "decijoule", "decijoules";
        @centijoule: prefix!(centi) / (prefix!(kilo) / 6.022_140_76_E23); "cJ", "centijoule", "centijoules";
        @millijoule: prefix!(milli) / (prefix!(kilo) / 6.022_140_76_E23); "mJ", "milljoule", "milljoules";
        @microjoule: prefix!(micro) / (prefix!(kilo) / 6.022_140_76_E23); "µJ", "microjoule", "microjoules";
        @nanojoule: prefix!(nano) / (prefix!(kilo) / 6.022_140_76_E23); "nJ", "nanojoule", "nanojoules";
        @picojoule: prefix!(pico) / (prefix!(kilo) / 6.022_140_76_E23); "pJ", "picojoule", "picojoules";
        @femtojoule: prefix!(femto) / (prefix!(kilo) / 6.022_140_76_E23); "fJ", "femtojoule", "femtojoules";
        @attojoule: prefix!(atto) / (prefix!(kilo) / 6.022_140_76_E23); "aJ", "attojoule", "attojoules";
        @zeptojoule: prefix!(zepto) / (prefix!(kilo) / 6.022_140_76_E23); "zJ", "zeptojoule", "zeptojoules";
        @yoctojoule: prefix!(yocto) / (prefix!(kilo) / 6.022_140_76_E23); "yJ", "yoctojoule", "yoctojoules";

        @petawatt_hour: 3.6_E18 / (prefix!(kilo) / 6.022_140_76_E23); "PW · h", "petawatt hour", "petawatt hours";
        @terawatt_hour: 3.6_E15 / (prefix!(kilo) / 6.022_140_76_E23); "TW · h", "terawatt hour", "terawatt hours";
        @gigawatt_hour: 3.6_E12 / (prefix!(kilo) / 6.022_140_76_E23); "GW · h", "gigawatt hour", "gigawatt hours";
        @megawatt_hour: 3.6_E9 / (prefix!(kilo) / 6.022_140_76_E23); "MW · h", "megawatt hour", "megawatt hours";
        @kilowatt_hour: 3.6_E6 / (prefix!(kilo) / 6.022_140_76_E23); "kW · h", "kilowatt hour", "kilowatt hours";
        @hectowatt_hour: 3.6_E5 / (prefix!(kilo) / 6.022_140_76_E23); "hW · h", "hectowatt hour", "hectowatt hours";
        @decawatt_hour: 3.6_E4 / (prefix!(kilo) / 6.022_140_76_E23); "daW · h", "decawatt hour", "decawatt hours";
        @watt_hour: 3.6_E3 / (prefix!(kilo) / 6.022_140_76_E23); "W · h", "watt hour", "watt hours";
        @milliwatt_hour: 3.6_E0 / (prefix!(kilo) / 6.022_140_76_E23); "mW · h", "milliwatt hour", "milliwatt hours";
        @microwatt_hour: 3.6_E-3 / (prefix!(kilo) / 6.022_140_76_E23); "µW · h", "microwatt hour", "microwatt hours";

        @calorie: 4.184_E0 / (prefix!(kilo) / 6.022_140_76_E23); "cal", "calorie", "calories";
        @kilocalorie: prefix!(kilo) * 4.184_E0 / (prefix!(kilo) / 6.022_140_76_E23); "kcal", "kilocalorie", "kilocalories";

        @calorie_per_mole: 4.184_E0; "cal", "calorie", "calories";
        @kilocalorie_per_mole: prefix!(kilo) * 4.184_E0; "kcal", "kilocalorie", "kilocalories";

        @electronvolt: 1.602_177_E-19 / (prefix!(kilo) / 6.022_140_76_E23); "eV", "electronvolt", "electronvolts";
        @erg: 1.0_E-7 / (prefix!(kilo) / 6.022_140_76_E23); "erg", "erg", "ergs";
        @hartree: 2625.5; "Ha", "hartree", "hartrees";
    }
}
