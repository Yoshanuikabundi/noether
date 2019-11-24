quantity! {
    /// Entropy (base unit kilojoule per mole per kelvin, kg · m² · s⁻² · mol⁻1 · K⁻1 ).
    quantity: Entropy; "entropy";
    /// Dimension of entropy
    dimension: Q<
        P2,  // length
        P1,  // mass
        N2,  // time
        Z0,  // charge
        N1>; // temperature
    units {
        @kilojoule_per_mole_kelvin: 1.0; "KJ/mol/K", "kilojoule per mole kelvin", "kilojoules per mole kelvin";

        @calorie_per_mole_kelvin: 4.184_E0; "cal/mol/K", "calorie per mole kelvin", "calories per mole kelvin";
        @kilocalorie_per_mole_kelvin: prefix!(kilo) * 4.184_E0; "kcal/mol/K", "kilocalorie per mole kelvin", "kilocalories per mole per kelvin";

        @electronvolt_per_kelvin: 1.602_177_E-19 / (prefix!(kilo) / 6.022_140_76_E23); "eV/K", "electronvolt per kelvin", "electronvolts per kelvin";
        @erg_per_kelvin: 1.0_E-7 / (prefix!(kilo) / 6.022_140_76_E23); "erg/K", "erg per kelvin", "ergs per kelvin";
        @hartree_per_kelvin: 2625.5; "Ha/K", "hartree per kelvin", "hartrees per kelvin";

        @boltzmann_constant: 8.314_462_1_E-3; "kB", "Boltzmann constant", "times the Boltzmann constant";
    }
}
