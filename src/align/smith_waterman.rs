//! Smith-Waterman algorithm
//! https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
//!
//! 1981 local alignment algorithm.
//! Perform local sequence alignment: find similar _regions_.
//! This is a variation of the Needleman-Wunsch algorithm.
//!
//! Guaranteed to find the optional local alignments relative
//! to the scoring matrix being used.
//!
//! The main difference to the Needleman–Wunsch algorithm is that
//! negative scoring matrix cells are set to zero, which renders the
//! (thus positively scoring) local alignments visible.
//!
//! Traceback procedure starts at the highest scoring matrix cell and
//! proceeds until a cell with score zero is encountered, yielding the
//! highest scoring local alignment.
//!
//! The Smith–Waterman algorithm finds the segments in two sequences
//! that have similarities while the Needleman–Wunsch algorithm aligns
//! two complete sequences.
//!
//! There is no one single type of gap penalty. The default is a penalty for
//! a gap, but often extending a gap is penalized less than the initial
//! opening of the gap: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty
//!
//! The default I'm familiar with is the linear penalty. A one-time penalty for
//! an adjacent gap in the alignment.

use super::{Align, Alignment, Scoring};

struct Aligner<'a> {
    grid: Vec<Vec<i32>>,
    a: &'a str,
    b: &'a str,
    scoring: Scoring,
}

impl<'a> Align for Aligner<'a> {
    fn align(&mut self) -> Alignment {
        self.init();
        self.fill();
        let (a, b) = self.backtrace();

        Alignment {
            grid: self.grid.to_owned(),
            a,
            b,
            a_orig: self.a.to_string(),
            b_orig: self.b.to_string(),
        }
    }
}

impl<'a> Aligner<'a> {
    pub fn new(a: &'a str, b: &'a str, scoring: Scoring) -> Aligner<'a> {
        Aligner {
            grid: Vec::new(),
            a,
            b,
            scoring,
        }
    }

    fn init(&mut self) {
        for _i in 0..self.b.len() + 1 {
            // unlike needleman-wunsch, these rows start as 0 (rather than -i for a global alignment)
            self.grid.push(vec![0; self.a.len() + 1]);
        }
    }

    fn fill(&mut self) {
        let a = self.a.as_bytes();
        let b = self.b.as_bytes();

        for i in 1..b.len() + 1 {
            for j in 1..a.len() + 1 {
                let mut opts: Vec<i32> = vec![
                    // negative values are ignored, 0 is as low as we'll go
                    0,
                    // gap
                    self.grid[i - 1][j] + self.scoring.gap,
                    self.grid[i][j - 1] + self.scoring.gap,
                ];

                // match or mismatch
                let match_val = self
                    .scoring
                    .rm
                    .get(&a[j - 1])
                    .unwrap()
                    .get(&b[i - 1])
                    .unwrap();
                opts.push(self.grid[i - 1][j - 1] + match_val);

                self.grid[i][j] = opts.iter().max().unwrap().to_owned();
            }
        }
    }

    fn backtrace(&self) -> (String, String) {
        // this finds the global maximum among the alignments
        let mut i = 0;
        let mut j = 0;
        for ti in 0..self.grid.len() {
            for tj in 0..self.grid[ti].len() {
                if self.grid[ti][tj] > self.grid[i][j] {
                    i = ti;
                    j = tj;
                }
            }
        }

        let a = self.a.as_bytes();
        let b = self.b.as_bytes();
        let mut a_row: Vec<u8> = Vec::new();
        let mut b_row: Vec<u8> = Vec::new();
        while i > 0 || j > 0 {
            let mut opts: Vec<i32> = Vec::new();
            if i > 0 && j > 0 {
                opts.push(self.grid[i - 1][j - 1]); // match or mismatch
            }
            if i > 0 {
                opts.push(self.grid[i - 1][j]); // gap
            }
            if j > 0 {
                opts.push(self.grid[i][j - 1]); // gap
            }

            let min = opts.iter().min().unwrap().to_owned();
            if min == 0 {
                a_row.push(a[j - 1]);
                b_row.push(b[i - 1]);
                break; // reached end of local alignment
            }

            let max = opts.iter().max().unwrap().to_owned();
            if i > 0 && j > 0 && self.grid[i - 1][j - 1] == max {
                a_row.push(a[j - 1]);
                b_row.push(b[i - 1]);
                i -= 1;
                j -= 1;
            } else if i > 0 && self.grid[i - 1][j] == max {
                a_row.push('-' as u8);
                b_row.push(b[i - 1]);
                i -= 1;
            } else if j > 0 && self.grid[i][j - 1] == max {
                a_row.push(a[j - 1]);
                b_row.push('-' as u8);
                j -= 1;
            }
        }

        let top: String = a_row.iter().rev().map(|c| *c as char).collect();
        let bottom: String = b_row.iter().rev().map(|c| *c as char).collect();
        (top, bottom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::matrix::MATRIX;

    #[test]
    fn test_aligner_align() {
        let mut a = Aligner::new(
            "TGTTACGG",
            "GGTTGACTA",
            Scoring {
                rm: MATRIX::NUC.read(),
                gap: -2,
                gap_extend: -1,
            },
        );
        let alignment = a.align();

        println!("{:?}", alignment);

        assert_eq!("GTT-AC", alignment.a);
        assert_eq!("GTTGAC", alignment.b);
    }
}
