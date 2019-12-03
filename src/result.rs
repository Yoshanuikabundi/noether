/// Error type enum for this crate
#[derive(Debug)]
pub enum Error {
    /// Error variant when the minimum image convention is not
    /// justified for a combination of box size and cutoff.
    ///
    /// If the cutoff for an interaction is larger than the smallest
    /// dimension of the box, then it's possible that a pair of atoms
    /// should appear twice in the pairlist. Because constructing the
    /// pairlist is an $\mathcal{O}(N^2)$ operation already, we use
    /// the minimum image convention to construct the pairlist, which
    /// only permits a pair of atoms to appear once.
    MinimumImageConventionNotJustified,

    /// Error variant for illegal combinations of function arguments.
    ///
    /// The goal should be to eventually refactor to remove this form
    /// of error and make illegal combinations unrepresentable. I don't
    /// know if that's realistic yet though.
    ValueError(&'static str)
}

impl Error {
    /// Convert an error variant to its own name
    ///
    /// Used to point users to documentation
    fn get_variant_str(&self) -> &'static str{
        match self {
            MinimumImageConventionNotJustified => "MinimumImageConventionNotJustified",
            ValueError(_) => "ValueError",
        }
    }

    /// Get a brief description for an error
    fn get_description_str(&self) -> String {
        match self {
            MinimumImageConventionNotJustified => concat!(
                "Minimum image convention not justified. ",
                "Make the cutoff smaller or the box bigger."
            ).to_string(),
            ValueError(s) => format!(
                "Illegal combination of arguments: {}",
                s
            ),
        }
    }

    pub fn explain(self) {
        let variant = self.get_variant_str();
        let description = self.get_description_str();

        let url = match std::env::current_exe() {
            Ok(mut path) => {
                path.pop();
                path.pop();
                path.push("doc");
                format!("file://{}", path.display())
            },
            Err(_) => format!("https://docs.rs/noether/{}", crate::VERSION)
        };

        let url = format!(
            "{}/noether/result/enum.Error.html#variant.{}",
            url,
            variant
        );

        println!("{}", description);
        println!("For more information, see {}", url);
    }
}

pub use Error::*;

/// Result type for this module
pub type Result<T> = std::result::Result<T, Error>;

pub trait FriendlyResult {
    type Value;

    fn unwrap_nicely(self) -> Self::Value;
}

impl<T> FriendlyResult for Result<T> {
    type Value = T;
    fn unwrap_nicely(self) -> Self::Value {
        match self {
            Ok(v) => v,
            Err(e) => {
                println!("-----------------------");
                e.explain();
                println!("-----------------------");
                panic!();
            }
        }
    }
}
