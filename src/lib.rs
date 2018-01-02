#![feature(test)]

#[macro_use]
extern crate itertools;
#[macro_use]
extern crate serde_derive;
extern crate serde;

#[cfg(test)]
extern crate test;

mod bitvec;
//mod pbwt;

pub use bitvec::BitVec;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
