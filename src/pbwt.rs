use std::cmp::min;
use super::BitVec;
use std::iter::FromIterator;

pub type Range = ::std::ops::Range<usize>;

#[derive(Clone, Debug)]
pub struct PBWTState<L=usize> {
    /// Number of rows.
    pub n: usize,
    /// Row labels. Just cosmetic, may be any type.
    pub labels: Vec<L>,
    /// The number of columns already traversed.
    pub position: usize,
    /// The current permutaton of the rows.
    pub permutation: Vec<usize>,
    /// The next column, permuted by `permutation`.
    next_col: Option<BitVec>,
    /// Sum of ones in `next_bits` in `0..i` (not including `i`). Has length `N+1`.
    /// Only valid when `next_bits.is_some()`.
    true_count: Vec<usize>,
    /// Length of common suffix (match depth) of rows `i-1` and `i`.
    /// Has length `N` and `true_count[0]=0`.
    pub depths: Vec<usize>,
    // A static min-interval tree with the common suffix length of i and i+1 seq.
    // The first part is one smaller than length of `permutation`, the following parts
    // are the minimums of 2-tuples, 4-tuples, ... up to the entire tree min.
    // (Incomplete tuples are also minimized over.)
    //common_max_tree: Vec<Index>,
    // Indexes into `max_tree` with starts of windows of size 1, 2, 4, 8, ...
    // up to the entire tree.
    //max_tree_starts: Vec<usize>,
}

#[allow(dead_code)]
impl PBWTState<usize> {
    pub fn new(n: usize) -> Self {
        Self::with_labels((0..n).collect::<Vec<_>>())
    }
}

#[allow(dead_code)]
impl<L> PBWTState<L> {
    pub fn with_labels<T: Into<Vec<L>>>(labels: T) -> Self {
        let l: Vec<L> = labels.into();
        let n = l.len();
        assert!(n > 0);
        PBWTState {
            n: n,
            labels: l,
            position: 0,
            permutation: (0..n).collect(),
            next_col: None,
            true_count: vec![0; n + 1],
            depths: vec![0; n],
        }
        //s.common_max_tree.resize(s.permutation.len() - 1, 0);
        //s.compute_max_tree();
        //s
    }

    /// Sets the next column. The current `next_col` must be `None`.
    pub fn set_next_column(&mut self, next_col: BitVec) {
        assert!(self.next_col.is_none());
        assert_eq!(self.n, next_col.len());
        let next = BitVec::from_iter((0..self.n).map(
            |i| next_col.get(self.permutation[i])));
        self.true_count.truncate(0);
        self.true_count.push(0);
        self.true_count.extend(next.iter().scan(0,
            |st, x| { if x {*st += 1}; Some(*st) } ));
        self.next_col = Some(next);
    }

    /// Advance the internal state over the `next_col`, consuming it.
    /// Panics when `next_col` is `None`.
    pub fn advance(&mut self) {
        let next = self.next_col.take().expect("next_col must be set.");
        let mut perm = Vec::<usize>::with_capacity(self.n);
        let mut depth = Vec::<usize>::with_capacity(self.n);
        // First all zeros, then all ones
        for val in [false, true].iter() {
            let mut d = self.position;
            let mut first = true;
            for (i, v) in next.iter().enumerate() {
                d = min(d, self.depths[i]);
                if v == *val {
                    perm.push(self.permutation[i]);
                    depth.push(if first { 0 } else { d + 1 });
                    first = false;
                    d = self.position;
                }
            }
        }
        self.permutation = perm;
        self.depths = depth;
        self.position += 1;
    }

    pub fn extend_split(&self, r: &Range) -> [Range; 2] {
        [self.extend_with(r, false), self.extend_with(r, true)]
    }

    #[inline]
    pub fn extend_with(&self, r: &Range, val: bool) -> Range {
        self.extend_single_with(r.start, val) ..
            self.extend_single_with(r.end, val)
    }

    #[inline]
    pub fn extend_single_with(&self, i: usize, val: bool) -> usize {
        assert!(self.next_col.is_some());
        if val {
            self.n - self.true_count[self.n] + self.true_count[i]
        } else {
            i - self.true_count[i]
        }
    }

/*
    fn compute_max_tree(&mut self) {
        self.common_max_tree.truncate(self.len() - 1);
        self.max_tree_starts.clear();
        self.max_tree_starts.push(0);
        let mut wlen = 1;
        while wlen < self.len() {
            wlen *= 2;
            let start = *self.max_tree_starts.last().unwrap();
            let stop = self.common_max_tree.len();
            self.max_tree_starts.push(stop);
            for i in 0 .. (stop - start) / 2 {
                self.common_max_tree.push(min(
                    self.common_max_tree[start + 2 * i],
                    self.common_max_tree[start + 2 * i + 1]));
            }
            if (stop - start) % 2 != 0 {
                self.common_max_tree.push(self.common_max_tree[stop - 1]);
            }
        }
    }
    */
}

#[cfg(test)]
mod tests {
    use super::PBWTState;
    use super::super::BitVec;

    #[test]
    fn test_basic() {
        let mut p = PBWTState::new(6);
        assert_eq!(p.n, 6);
        assert_eq!(p.next_col, None);
        assert_eq!(p.depths, [0, 0, 0, 0, 0, 0]);

        p.set_next_column(bit_vec![0,1,1,0,1,0]);
        assert_eq!(p.extend_with(&(2..5), false), 1..2);
        assert_eq!(p.extend_with(&(2..5), true), 4..6);
        assert_eq!(p.extend_split(&(1..3)), [1..1, 3..5]);
        assert_eq!(p.extend_split(&(0..6)), [0..3, 3..6]);
        assert_eq!(p.extend_split(&(0..0)), [0..0, 3..3]);
        assert_eq!(p.extend_split(&(6..6)), [3..3, 6..6]);
        p.advance();
        assert_eq!(p.depths, [0, 1, 1, 0, 1, 1]);

        p.set_next_column(bit_vec![1,0,1,0,0,1]);
        assert_eq!(p.extend_split(&(1..5)), [0..2, 4..6]);
        p.advance();
        assert_eq!(p.depths, [0, 1, 2, 0, 2, 1]);
        assert_eq!(p.permutation, [3, 1, 4, 0, 5, 2]);
        assert_eq!(p.position, 2);
    }
}
