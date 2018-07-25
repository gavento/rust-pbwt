use super::BitVec;
use std::{ops::Range, cmp::max, iter::FromIterator, mem::size_of};

/// [1] Richard Durbin: 
/// Efficient haplotype matching and storage using the positional
/// Burrows–Wheeler transform (PBWT)

#[derive(Clone, Debug)]
pub struct PBWTState {

    /// Number of rows. Matches `M` in [1].
    pub m: usize,

    /// The number of columns already traversed.
    /// Called `k` in [1].
    pub position: usize,

    /// The current permutaton of the rows.
    /// Called `a_k` in [1].
    pub permutation: Vec<usize>,

    /// The next column, already permuted by `permutation`.
    next_col: Option<BitVec>,

    /// Sum of ones in `next_bits` in `0..i` (not including `i`). Has length `m+1`.
    /// Only valid when `next_col.is_some()`.
    true_count: Vec<usize>,

    /// Start of the common suffix of rows `i-1` and `i` (under the current permutation).
    /// Called `d_k` in [1]. Has length `m` and `match_start[0]=position`.
    pub match_start: Vec<usize>,
    
    /// Maxima interval tree for `match_start`. 
    /// `match_start_maxima[l][k]` contains the maximum of
    /// `match_start[(2 << l) * k .. (2 << l) * (k + 1) - 1]`.
    pub match_start_maxima: Vec<Vec<usize>>,
}

#[allow(dead_code)]
impl PBWTState {
    pub fn new(m: usize) -> Self {
        assert!(m > 0);
        let l = 1 + 8 * size_of::<usize>() - (m.leading_zeros() as usize);
        let mut s = PBWTState {
            m: m,
            position: 0,
            permutation: (0..m).collect(),
            next_col: None,
            true_count: Vec::new(),
            match_start: vec![0; m],
            match_start_maxima: vec![Vec::new(); l],
        };
        s.compute_max_tree();
        s
    }

    /// Sets the next column. The current `next_col` must be `None`.
    pub fn set_next_column(&mut self, next_col: BitVec) {
        assert!(self.next_col.is_none());
        assert_eq!(self.m, next_col.len());
        let next = BitVec::from_iter((0..self.m).map(
            |i| next_col.get(self.permutation[i])));
        self.true_count.truncate(0);
        self.true_count.push(0);
        self.true_count.extend(next.iter().scan(0,
            |st, x| { if x {*st += 1}; Some(*st) } ));
        self.next_col = Some(next);
    }

    /// Advance the internal state over the `next_col`.
    /// Returns a new structure. Panics when `next_col` is `None`.
    pub fn advance(&self) -> Self {
        let next = self.next_col.as_ref().expect("next_col must be set.");
        assert_eq!(self.true_count.len(), self.m + 1);
        let mut perm = Vec::<usize>::with_capacity(self.m);
        let mut match_start = Vec::<usize>::with_capacity(self.m);
        // First all zeros, then all ones
        for &val in &[false, true] {
            let mut d = self.position;
            let mut first = true;
            for (i, v) in next.iter().enumerate() {
                d = max(d, self.match_start[i]);
                if v == val {
                    perm.push(self.permutation[i]);
                    match_start.push(if first { self.position + 1 } else { d });
                    first = false;
                    d = 0;
                }
            }
        }
        let mut s = PBWTState {
            m: self.m,
            position: self.position + 1,
            permutation: perm,
            next_col: None,
            true_count: Vec::new(),
            match_start: match_start,
            match_start_maxima: vec![Vec::new(); self.match_start_maxima.len()],
        };
        s.compute_max_tree();
        s
    }

    pub fn extend_split(&self, r: &Range<usize>) -> [Range<usize>; 2] {
        [self.extend_with(r, false), self.extend_with(r, true)]
    }

    #[inline]
    pub fn extend_with(&self, r: &Range<usize>, val: bool) -> Range<usize> {
        self.extend_single_with(r.start, val) ..
            self.extend_single_with(r.end, val)
    }

    #[inline]
    pub fn extend_single_with(&self, i: usize, val: bool) -> usize {
        assert!(self.next_col.is_some());
        if val {
            self.m - self.true_count[self.m] + self.true_count[i]
        } else {
            i - self.true_count[i]
        }
    }

    /// Computes the max-tree in `match_start_maxima` for as many
    /// levels as `match_start_maxima` already has (NOP if it is empty).
    fn compute_max_tree(&mut self) {
        let l = self.match_start_maxima.len();
        self.match_start_maxima.clear();
        for d in 0 .. l {
            let mut ms = Vec::with_capacity(3 + (self.m >> (d + 1)));
            {
                let up = if d == 0 {
                    &self.match_start
                } else {
                    &self.match_start_maxima.last().unwrap()
                };
                for i in 0 .. (up.len() + 1) / 2 {
                    ms.push(max(up[2 * i], up[2 * i + 1]));
                }
            }
            self.match_start_maxima.push(ms);
        }
    }

    /// Return a (closed..open) range around `around_pos` that has match start <= `match_start`.
    /// A slow, linear-time version.
    fn match_start_range_slow(&self, match_start: usize, around_pos: usize) -> PBWTCursor {
        let mut a = around_pos;
        while a > 0 && self.match_start[a] <= match_start {
            a -= 1;
        }
        let mut b = around_pos + 1;
        while b < self.m && self.match_start[b] <= match_start {
            b += 1;
        }
        PBWTCursor {
            range: a .. b,
        }
    }

    /// Return a (closed..open) range around `around_pos` that has match start <= `match_start`.
    pub fn match_start_range(&self, match_start: usize, around_pos: usize) -> PBWTCursor {
        // TODO: log-time solution
        self.match_start_range_slow(match_start, around_pos)
    }
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct PBWTCursor {
    pub range: ::std::ops::Range<usize>,
    // Not working yet:
    //pub match_start: usize,
}

impl PBWTCursor {
    pub fn new_full_range(state: &PBWTState) -> Self {
        PBWTCursor {
            range: 0..state.m,
        }
    }

    pub fn extend_both(&self, state: &PBWTState) -> [Self; 2] {
        [self.extend_with(state, false), self.extend_with(state, true)]
    }

    pub fn extend_with(&self, state: &PBWTState, val: bool) -> Self {
        PBWTCursor {
            range: state.extend_single_with(self.range.start, val) .. 
                   state.extend_single_with(self.range.end, val),
        }
    }

    pub fn len(&self) -> usize {
        self.range.len()
    }
}


#[cfg(test)]
mod tests {
    use super::PBWTState;
    use super::super::BitVec;

    #[test]
    fn test_basic() {
        // Not permuted input:
        // 0 1
        // 1 0
        // 1 1
        // 0 0
        // 1 0
        // 0 1
        let mut p = PBWTState::new(6);
        assert_eq!(p.m, 6);
        assert_eq!(p.next_col, None);
        assert_eq!(p.match_start, [0, 0, 0, 0, 0, 0]);

        p.set_next_column(bit_vec![0,1,1,0,1,0]);
        assert_eq!(p.extend_with(&(2..5), false), 1..2);
        assert_eq!(p.extend_with(&(2..5), true), 4..6);
        assert_eq!(p.extend_split(&(1..3)), [1..1, 3..5]);
        assert_eq!(p.extend_split(&(0..6)), [0..3, 3..6]);
        assert_eq!(p.extend_split(&(0..0)), [0..0, 3..3]);
        assert_eq!(p.extend_split(&(6..6)), [3..3, 6..6]);
        let mut p = p.advance();
        assert_eq!(p.match_start, [1, 0, 0, 1, 0, 0]);

        p.set_next_column(bit_vec![1,0,1,0,0,1]);
        assert_eq!(p.extend_split(&(1..5)), [0..2, 4..6]);
        let p = p.advance();
        // Sorted at this point:
        // 3. 0 0
        // 1. 1 0
        // 4. 1 0
        // 0. 0 1
        // 5. 0 1
        // 2. 1 1
        assert_eq!(p.match_start, [2, 1, 0, 2, 0, 1]);
        assert_eq!(p.permutation, [3, 1, 4, 0, 5, 2]);
        assert_eq!(p.position, 2);
    }
}
