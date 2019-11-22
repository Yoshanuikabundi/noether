system! {
    quantities: Q {
        length: nanometer, L;
        mass: dalton, M;
        time: picosecond, Ti;
        charge: elementary_charge, Q;
        temperature: kelvin, Te;
    }

    units: U {
        length::Length,
        mass::Mass,
        time::Time,
        charge::Charge,
        temperature::Temperature,

        energy::Energy,
    }
}

mod f32 {
    Q!(crate::units, f32);
}
