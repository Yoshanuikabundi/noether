quantity! {
    /// Charge (base unit elementary charge, e).
    quantity: Charge; "charge";
    /// Charge dimension, e.
    dimension: Q<
        Z0,  // length
        Z0,  // mass
        Z0,  // time
        P1,  // charge
        Z0>; // temperature
    units {
        @elementary_charge: 1.602_176_62_E-19 / 1.602_176_62_E-19; "e", "elementary charge", "elementary charges";

        @yottacoulomb: prefix!(yotta) / 1.602_176_62_E-19; "YC", "yottacoulomb", "yottacoulombs";
        @zettacoulomb: prefix!(zetta) / 1.602_176_62_E-19; "ZC", "zettacoulomb", "zettacoulombs";
        @exacoulomb: prefix!(exa) / 1.602_176_62_E-19; "EC", "exacoulomb", "exacoulombs";
        @petacoulomb: prefix!(peta) / 1.602_176_62_E-19; "PC", "petacoulomb", "petacoulombs";
        @teracoulomb: prefix!(tera) / 1.602_176_62_E-19; "TC", "teracoulomb", "teracoulombs";
        @gigacoulomb: prefix!(giga) / 1.602_176_62_E-19; "GC", "gigacoulomb", "gigacoulombs";
        @megacoulomb: prefix!(mega) / 1.602_176_62_E-19; "MC", "megacoulomb", "megacoulombs";
        @kilocoulomb: prefix!(kilo) / 1.602_176_62_E-19; "kC", "kilocoulomb", "kilocoulombs";
        @hectocoulomb: prefix!(hecto) / 1.602_176_62_E-19; "hC", "hectocoulomb", "hectocoulombs";
        @decacoulomb: prefix!(deca) / 1.602_176_62_E-19; "daC", "decacoulomb", "decacoulombs";
        /// Derived unit of electric charge.
        @coulomb: prefix!(none) / 1.602_176_62_E-19; "C", "coulomb", "coulombs";
        @decicoulomb: prefix!(deci) / 1.602_176_62_E-19; "dC", "decicoulomb", "decicoulombs";
        @centicoulomb: prefix!(centi) / 1.602_176_62_E-19; "cC", "centicoulomb", "centicoulombs";
        @millicoulomb: prefix!(milli) / 1.602_176_62_E-19; "mC", "millcoulomb", "millcoulombs";
        @microcoulomb: prefix!(micro) / 1.602_176_62_E-19; "µC", "microcoulomb", "microcoulombs";
        @nanocoulomb: prefix!(nano) / 1.602_176_62_E-19; "nC", "nanocoulomb", "nanocoulombs";
        @picocoulomb: prefix!(pico) / 1.602_176_62_E-19; "pC", "picocoulomb", "picocoulombs";
        @femtocoulomb: prefix!(femto) / 1.602_176_62_E-19; "fC", "femtocoulomb", "femtocoulombs";
        @attocoulomb: prefix!(atto) / 1.602_176_62_E-19; "aC", "attocoulomb", "attocoulombs";
        @zeptocoulomb: prefix!(zepto) / 1.602_176_62_E-19; "zC", "zeptocoulomb", "zeptocoulombs";
        @yoctocoulomb: prefix!(yocto) / 1.602_176_62_E-19; "yC", "yoctocoulomb", "yoctocoulombs";

        @petaampere_hour: 3.6_E18 / 1.602_176_62_E-19; "PA · h", "petaampere hour", "petaampere hours";
        @teraampere_hour: 3.6_E15 / 1.602_176_62_E-19; "TA · h", "teraampere hour", "teraampere hours";
        @gigaampere_hour: 3.6_E12 / 1.602_176_62_E-19; "GA · h", "gigaampere hour", "gigaampere hours";
        @megaampere_hour: 3.6_E9 / 1.602_176_62_E-19; "MA · h", "megaampere hour", "megaampere hours";
        @kiloampere_hour: 3.6_E6 / 1.602_176_62_E-19; "kA · h", "kiloampere hour", "kiloampere hours";
        @hectoampere_hour: 3.6_E5 / 1.602_176_62_E-19; "hA · h", "hectoampere hour", "hectoampere hours";
        @decaampere_hour: 3.6_E4 / 1.602_176_62_E-19; "daA · h", "decaampere hour", "decaampere hours";
        @ampere_hour: 3.6_E3 / 1.602_176_62_E-19; "A · h", "ampere hour", "ampere hours";
        @milliampere_hour: 3.6_E0 / 1.602_176_62_E-19; "mA · h", "milliampere hour", "milliampere hours";
        @microampere_hour: 3.6_E-3 / 1.602_176_62_E-19; "µA · h", "microampere hour", "microampere hours";

        @abcoulomb: 1.0_E1 / 1.602_176_62_E-19; "abC", "abcoulomb", "abcoulombs";
        @faraday: 9.648_531_E4 / 1.602_176_62_E-19; "F", "faraday", "faradays";
        @franklin: 3.335_641_E-10 / 1.602_176_62_E-19; "Fr", "franklin", "franklins";
        @statcoulomb: 3.335_641_E-10 / 1.602_176_62_E-19; "statC", "statcoulomb", "statcoulombs";
    }
}
