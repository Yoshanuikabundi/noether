#[macro_use]
extern crate uom;
#[macro_use]
extern crate more_asserts;

pub mod units;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
