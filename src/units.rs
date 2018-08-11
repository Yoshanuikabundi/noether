make_units! {
    GMX;
    ONE: Unitless;

    base {
        NM: Nanometer, "nm", Length;
        DA: Dalton, "Da", Mass;
        PS: Picosecond, "ps", Time;
        E: ElemCharge, "e", Charge;
        K: Kelvin, "K", Temperature;
    }

    derived {
        NMPPS: NanometerPerPicosecond = (Nanometer / Picosecond), Velocity;
        NM2: Nanometer2 = (Nanometer * Nanometer), Area;
        NM3: Nanometer3 = (Nanometer2 * Nanometer), Volume;
        PS2: Picosecond2 = (Picosecond * Picosecond);
        KJPM: KilojoulePerMole = (Nanometer2 * Dalton / Picosecond2), Energy;
        KJPMNM: KilojoulePerMolePerNanometer = (KilojoulePerMole / Nanometer), Force;
        KJPMNM3: KilojoulePerMolePerNanometer3 = (KilojoulePerMole / Nanometer3), Pressure;
        ENM: ElemChargeNanometer = (ElemCharge * Nanometer);
        KJPME: KilojoulePerMolePerElemCharge = (KilojoulePerMole / ElemCharge), ElectricPotential;
        KJPMNME: KilojoulePerMolePerNanometerPerElemCharge = (KilojoulePerMolePerNanometer / ElemCharge);
    }

    constants {
        // Ångström
        A: Nanometer = 0.1;

        // KiloDalton
        KDA: Dalton = 1.0E3;
        // MegaDalton
        MDA: Dalton = 1.0E6;
        // Gram
        G: Dalton = 6.022_140_857E24;

        // Femtosecond
        FS: Picosecond = 1.0E-3;
        // Nanosecond
        NS: Picosecond = 1.0E3;
        // Microsecond
        US: Picosecond = 1.0E6;
        // Millisecond
        MS: Picosecond = 1.0E9;
        // Second
        S: Picosecond = 1.0E12;

        // Kilometers per second
        KMPS: NanometerPerPicosecond = 1.0;

        // Bar (Pressure unit, 100 kPA, ~1 atm)
        BAR: KilojoulePerMolePerNanometer3 = 16.605_390_404;

        PI: Unitless = consts::PI;
    }

    fmt = true;
}
