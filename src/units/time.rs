#[cfg(feature = "try-from")]
use lib::time::Duration;
#[cfg(feature = "try-from")]
use num::{FromPrimitive, ToPrimitive, Zero};


quantity! {
    /// Mass (base unit picosecond, ps).
    quantity: Time; "time";
    /// Time dimension, ps.
    dimension: Q<
        Z0,  // length
        Z0,  // mass
        P1,  // time
        Z0,  // charge
        Z0>; // temperature
    units {
        @picosecond: prefix!(pico) / prefix!(pico); "ps", "picosecond", "picoseconds";

        @yottasecond: prefix!(yotta) / prefix!(pico); "Ys", "yottasecond", "yottaseconds";
        @zettasecond: prefix!(zetta) / prefix!(pico); "Zs", "zettasecond", "zettaseconds";
        @exasecond: prefix!(exa) / prefix!(pico); "Es", "exasecond", "exaseconds";
        @petasecond: prefix!(peta) / prefix!(pico); "Ps", "petasecond", "petaseconds";
        @terasecond: prefix!(tera) / prefix!(pico); "Ts", "terasecond", "teraseconds";
        @gigasecond: prefix!(giga) / prefix!(pico); "Gs", "gigasecond", "gigaseconds";
        @megasecond: prefix!(mega) / prefix!(pico); "Ms", "megasecond", "megaseconds";
        @kilosecond: prefix!(kilo) / prefix!(pico); "ks", "kilosecond", "kiloseconds";
        @hectosecond: prefix!(hecto) / prefix!(pico); "hs", "hectosecond", "hectoseconds";
        @decasecond: prefix!(deca) / prefix!(pico); "das", "decasecond", "decaseconds";
        /// The second is the SI unit of time. It is defined by taking the fixed numerical value of
        /// the caesium frequency ∆*ν*<sub>Cs</sub>, the unperturbed ground-state hyperfine
        /// transition frequency of the caesium 133 atom, to be 9 192 631 770 when expressed in the
        /// unit Hz, which is equal to s⁻¹.
        @second: prefix!(none) / prefix!(pico); "s", "second", "seconds";
        @decisecond: prefix!(deci) / prefix!(pico); "ds", "decisecond", "deciseconds";
        @centisecond: prefix!(centi) / prefix!(pico); "cs", "centisecond", "centiseconds";
        @millisecond: prefix!(milli) / prefix!(pico); "ms", "millisecond", "milliseconds";
        @microsecond: prefix!(micro) / prefix!(pico); "µs", "microsecond", "microseconds";
        @nanosecond: prefix!(nano) / prefix!(pico); "ns", "nanosecond", "nanoseconds";
        @femtosecond: prefix!(femto) / prefix!(pico); "fs", "femtosecond", "femtoseconds";
        @attosecond: prefix!(atto) / prefix!(pico); "as", "attosecond", "attoseconds";
        @zeptosecond: prefix!(zepto) / prefix!(pico); "zs", "zeptosecond", "zeptoseconds";
        @yoctosecond: prefix!(yocto) / prefix!(pico); "ys", "yoctosecond", "yoctoseconds";

        @day: 8.64_E4 / prefix!(pico); "d", "day", "days";
        @hour: 3.6_E3 / prefix!(pico); "h", "hour", "hours";
        @minute: 6.0_E1 / prefix!(pico); "min", "minute", "minutes";
        @shake: 1.0_E-8 / prefix!(pico); "10.0 ns", "shake", "shakes";
        @year: 3.1536_E7 / prefix!(pico); "a", "year", "years";

    }
}
