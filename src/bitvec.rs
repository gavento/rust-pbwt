use std::vec::Vec;
use std::ops::Index;
use std::iter::FromIterator;
use std::fmt;

/// Type of BitVec chunk
pub(super) type BitVecType = u64;

/// Size of the BitVec chunk type
pub const CHUNK_SIZE: usize = 64;

#[derive(Clone, Serialize, Deserialize)]
pub struct BitVec {
    len: usize,
    data: Vec<BitVecType>,
}

impl BitVec {
    /// Create a new empty vector
    pub fn new() -> Self {
        BitVec {
            len: 0,
            data: Vec::new(),
        }
    }

    /// Create a new empty vector with pre-allocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        BitVec {
            len: 0,
            data: Vec::with_capacity((capacity + CHUNK_SIZE - 1) / CHUNK_SIZE),
        }
    }

    /// Create a vector with `len` repetitions of given `val`
    pub fn repeated(len: usize, val: bool) -> Self {
        let mut v = Vec::new();
        v.resize((len + CHUNK_SIZE - 1) / CHUNK_SIZE,
                 if val { BitVecType::max_value() } else { 0 });
        BitVec { len: len, data: v }
    }

    /// Return the vector length in bits
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Get a bit value by a bit-index.
    /// Panics on out-of-bounds access.
    #[inline]
    pub fn get(&self, index: usize) -> bool {
        self.try_get(index).expect("BitVec index out of bounds")
    }

    /// Get a bit value by a bit-index, None for out-of bounds
    #[inline]
    pub fn try_get(&self, index: usize) -> Option<bool> {
        if index >= self.len {
            None
        } else {
            Some((self.data[index / CHUNK_SIZE] & (1 << (index % CHUNK_SIZE))) != 0)
        }
    }

    /// Private set without checking the bit-index bounds
    #[inline]
    fn priv_set(&mut self, index: usize, val: bool) {
        if val {
            self.data[index / CHUNK_SIZE] |= 1 << (index % CHUNK_SIZE);
        } else {
            self.data[index / CHUNK_SIZE] &= !(1 << (index % CHUNK_SIZE));
        }
    }

    /// Set a bit in the current range to a given value.
    /// Panics on out-of-bounds access.
    #[inline]
    pub fn set(&mut self, index: usize, val: bool) {
        if index >= self.len {
            panic!("BitVec index {} out of bounds (len={})", index, self.len);
        }
        self.priv_set(index, val);
    }

    /// Extend the vector with a single value
    pub fn push(&mut self, val: bool) {
        if self.len % CHUNK_SIZE == 0 {
            self.data.push(0);
        }
        let idx = self.len;
        self.len += 1;
        self.priv_set(idx, val);
    }

    /// Return a iterator ove the bits
    pub fn iter<'a>(&'a self) -> BitVecIter<'a> {
        BitVecIter {
            offset: 0,
            bitvec: self,
        }
    }

    /// Return the last chunk masked with only valid bits.
    /// Returns 0 for empty `BitVec`.
    fn last_chunk_masked(&self) -> BitVecType {
        let last = *self.data.last().unwrap_or(&0);
        if self.len % CHUNK_SIZE == 0 {
            last
        } else {
            last & ((1 << (self.len % CHUNK_SIZE)) - 1)
        }
    }

    /// Retun a vec of number of ones up to chunk `0, 1, ...` (inclusive)
    /// including the partial chunk (if any). Returns `[]` for empty `BitVec`.
    pub fn chunk_ones(&self) -> Vec<usize> {
        let mut res = Vec::with_capacity(self.data.len());
        let mut s = 0;
        for (i, w) in self.data.iter().enumerate() {
            if i == self.data.len() - 1 {
                s += self.last_chunk_masked().count_ones() as usize;
            } else {
                s += w.count_ones() as usize;
            }
            res.push(s);
        }
        res
    }

    /// Get the indicated chunk, the last partial chunk is properly masked.
    /// Panics on index out of bounds.
    pub fn get_chunk(&self, index: usize) -> BitVecType {
        if index < self.data.len() - 1 {
            self.data[index]
        } else {
            self.last_chunk_masked()
        }
    }
}

impl fmt::Debug for BitVec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "BitVec {{ ")?;
        for i in self.iter() {
            write!(f, "{}", i as u8)?;
        }
        write!(f, " }}")?;
        Ok(())
    }
}

impl fmt::Display for BitVec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in self.iter() {
            write!(f, "{}", i as u8)?;
        }
        Ok(())
    }
}

impl Eq for BitVec {}

impl PartialEq for BitVec {
    fn eq(&self, other: &Self) -> bool {
        if self.len != other.len {
            return false;
        }
        assert_eq!(self.data.len(), other.data.len());
        for i in 0..self.data.len() {
            if self.get_chunk(i) != other.get_chunk(i) {
                return false;
            }
        }
        true
    }
}

impl Index<usize> for BitVec {
    type Output = bool;

    fn index(&self, index: usize) -> &bool {
        static ST: bool = true;
        static SF: bool = false;
        if self.get(index) { &ST } else { &SF }
    }
}

impl FromIterator<bool> for BitVec {
    fn from_iter<T: IntoIterator<Item = bool>>(iter: T) -> Self {
        let mut bv = BitVec::new();
        for i in iter {
            bv.push(i);
        }
        bv
    }
}

impl<'a> FromIterator<&'a bool> for BitVec {
    fn from_iter<T: IntoIterator<Item = &'a bool>>(iter: T) -> Self {
        let mut bv = BitVec::new();
        for i in iter {
            bv.push(*i);
        }
        bv
    }
}

#[derive(Clone, Debug)]
pub struct BitVecIter<'a> {
    offset: usize,
    bitvec: &'a BitVec,
}

impl<'a> Iterator for BitVecIter<'a> {
    type Item = bool;
    fn next(&mut self) -> Option<Self::Item> {
        if self.offset >= self.bitvec.len {
            None
        } else {
            self.offset += 1;
            Some(self.bitvec.get(self.offset - 1))
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let rem = self.bitvec.len - self.offset;
        (rem, Some(rem))
    }

    fn count(self) -> usize {
        return self.size_hint().0;
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base() {
        let mut v = BitVec::new();
        assert_eq!(v.len, 0);
        v.push(true);
        v.push(true);
        v.push(false);
        v.push(true);
        assert_eq!(v.len, 4);
        assert_eq!(v.get(0), true);
        assert_eq!(v.get(2), false);
        v.set(1, false);
        v.set(2, true);
        assert_eq!(v.get(1), false);
        assert_eq!(v[2], true);
        assert_eq!(v.iter().collect::<Vec<bool>>(),
                   vec![true, false, true, true]);
        let v2 = BitVec::from_iter(&[true, false, true, true]);
        assert_eq!(v, v2);
    }

    #[test]
    fn test_long() {
        let mut v = BitVec::repeated(62, true);
        assert_eq!(v.data.len(), 1);
        v.push(false);
        v.push(false);
        assert_eq!(v.data.len(), 1);
        v.push(false);
        assert_eq!(v.data.len(), 2);
        assert_eq!(v.get(62), false);
        assert_eq!(v.get(64), false);
        v.set(64, true);
        assert_eq!(v.get(64), true);
        assert_eq!(v.len(), 65);
    }

    #[test]
    fn test_eq() {
        let mut va = BitVec::repeated(62, true);
        let mut vb = BitVec::repeated(62, true);
        va.push(true);
        vb.push(false);
        assert_eq!(va, va);
        assert_ne!(va, vb);
        let mut vc = BitVec::repeated(64, true);
        let mut vd = BitVec::repeated(64, true);
        assert_eq!(vc, vd);
        vc.push(true);
        assert_ne!(vc, vd);
        vd.push(false);
        assert_ne!(vc, vd);
    }

    #[test]
    fn test_bounds() {
        let mut v = BitVec::repeated(162, true);
        assert_panics!( v.get(162) );
        assert!(v.try_get(162).is_none());
    }
}
