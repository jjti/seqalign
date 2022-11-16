//! https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
//! 1970. Global optimal alignment, meaning total alignment between the two sequences.
//!
//! The original purpose of the algorithm described by Needleman and Wunsch was
//! to find similarities in the amino acid sequences of two proteins.
//!
//! The Needleman–Wunsch algorithm is still widely used for optimal global alignment,
//! particularly when the quality of the global alignment is of the utmost importance.
//! The algorithm assigns a score to every possible alignment, and the purpose of the
//! algorithm is to find all possible alignments having the highest score.
//!
//! TODO: add support for PAM and BLOSUM similarity matricies:
//! https://en.wikipedia.org/wiki/Point_accepted_mutation
//! https://en.wikipedia.org/wiki/BLOSUM
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
        }
    }
}

impl<'a> Aligner<'a> {
    pub fn new(a: &'a str, b: &'a str) -> Aligner<'a> {
        Aligner {
            grid: Vec::new(),
            a,
            b,
            scoring: Scoring {
                m: 1,
                mm: -1,
                indel: -1,
            },
        }
    }

    fn init(&mut self) {
        for i in 0..self.b.len() + 1 {
            self.grid.push(vec![0; self.a.len() + 1]);
            self.grid[i][0] = -(i as i32);
        }

        for j in 0..self.a.len() + 1 {
            self.grid[0][j] = -(j as i32);
        }
    }

    fn fill(&mut self) {
        let a = self.a.as_bytes();
        let b = self.b.as_bytes();

        for i in 1..b.len() + 1 {
            for j in 1..a.len() + 1 {
                let mut opts: Vec<i32> = vec![
                    // indel
                    self.grid[i - 1][j] + self.scoring.indel,
                    self.grid[i][j - 1] + self.scoring.indel,
                ];

                if a[j - 1] == b[i - 1] {
                    // match
                    opts.push(self.grid[i - 1][j - 1] + self.scoring.m);
                } else {
                    // mismatch
                    opts.push(self.grid[i - 1][j - 1] + self.scoring.mm);
                }

                self.grid[i][j] = opts.iter().max().unwrap().to_owned();
            }
        }
    }

    fn backtrace(&self) -> (String, String) {
        let a = self.a.as_bytes();
        let b = self.b.as_bytes();
        let mut a_row: Vec<u8> = Vec::new();
        let mut b_row: Vec<u8> = Vec::new();

        let mut i = self.b.len();
        let mut j = self.a.len();

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

    #[test]
    fn test_aligner_align() {
        let mut a = Aligner::new("GCATGCG", "GATTACA");
        let alignment = a.align();

        assert_eq!("GCAT-GCG", alignment.a);
        assert_eq!("G-ATTACA", alignment.b);
    }
}
