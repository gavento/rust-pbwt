#![feature(test)]

//#[macro_use]
extern crate itertools;
#[macro_use]
extern crate serde_derive;
extern crate serde;

#[cfg(test)]
extern crate test;

#[allow(unused_macros)]
macro_rules! assert_panics {
    ($b:expr) => {
        let res = ::std::panic::catch_unwind(|| { $b });
        if !res.is_err() {
            panic!{"assertion failed: expression was expected to panic: {}",
                stringify!($b)}
        }        
    };
}

#[macro_use]
mod bitvec;
mod pbwt;

pub use bitvec::BitVec;


