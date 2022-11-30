//! Smith-Waterman algorithm
//! https://dornsife.usc.edu/assets/sites/516/docs/papers/msw_papers/msw-042.pdf
//!
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

use super::{PWAlign, PWAlignment, Scoring, Step};

struct Aligner<'a> {
    grid: Vec<Vec<Step>>,
    a: &'a str,
    b: &'a str,
    scoring: Scoring,
}

impl<'a> PWAlign for Aligner<'a> {
    fn align(&mut self) -> PWAlignment {
        self.init();
        self.fill();
        let (a, b) = self.backtrace();

        PWAlignment {
            grid: self.grid.clone(),
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
        for _i in 0..=self.b.len() {
            // unlike needleman-wunsch, these rows start as 0 (rather than -i for a global alignment)
            self.grid.push(vec![Step::default(); self.a.len() + 1]);
        }

        for i in 0..=self.b.len() {
            self.grid[i][0] = Step::from(i, 0, 0);
        }

        for j in 0..=self.a.len() {
            self.grid[0][j] = Step::from(0, j, 0);
        }
    }

    fn fill(&mut self) {
        let a = self.a.as_bytes();
        let b = self.b.as_bytes();

        for i in 1..=b.len() {
            for j in 1..=a.len() {
                // negative values are ignored, 0 is as low as we'll go
                let mut opts: Vec<Step> = vec![Step {
                    val: 0,
                    i,
                    j,
                    next: None,
                }];

                // gaps
                let mut k = 1;
                while k < i {
                    if a[j - 1] == b[k - 1] {
                        opts.push(Step {
                            val: self.grid[k][j].val
                                + self.scoring.gap_opening
                                + self.scoring.gap_extension * (i - k - 1) as i32,
                            i,
                            j,
                            next: Some((k, j)),
                        });
                    }
                    k += 1;
                }

                let mut l = 1;
                while l < j {
                    if a[l - 1] == b[i - 1] {
                        opts.push(Step {
                            val: self.grid[i][l].val
                                + self.scoring.gap_opening
                                + self.scoring.gap_extension * (j - l - 1) as i32,
                            i,
                            j,
                            next: Some((i, l)),
                        });
                    }
                    l += 1;
                }

                // match or mismatch
                let match_val = self
                    .scoring
                    .replacement
                    .get(&a[j - 1])
                    .unwrap()
                    .get(&b[i - 1])
                    .unwrap();
                opts.push(Step {
                    val: self.grid[i - 1][j - 1].val + match_val,
                    i,
                    j,
                    next: Some((i - 1, j - 1)),
                });

                self.grid[i][j] = opts.iter().max().unwrap().clone();
            }
        }
    }

    fn backtrace(&self) -> (String, String) {
        // this finds the global maximum among the alignments
        let mut step = &Step::default();
        for i in 0..self.grid.len() {
            for j in 0..self.grid[i].len() {
                if &self.grid[i][j] > step {
                    step = &self.grid[i][j];
                }
            }
        }

        let a = self.a.as_bytes();
        let b = self.b.as_bytes();
        let mut a_row: Vec<u8> = Vec::new();
        let mut b_row: Vec<u8> = Vec::new();

        while let Some((next_i, next_j)) = step.next {
            let i_delta = step.i - next_i;
            let j_delta = step.j - next_j;

            if i_delta == 1 && j_delta == 1 {
                // match/mismatch
                a_row.push(a[step.j - 1]);
                b_row.push(b[step.i - 1]);
            } else if i_delta > 0 {
                // gap in seq a
                let mut i = step.i;
                while i > next_i {
                    a_row.push(b'-');
                    b_row.push(b[i - 1]);
                    i -= 1;
                }
            } else if j_delta > 0 {
                // gap in seq b
                let mut j = step.j;
                while j > next_j {
                    a_row.push(a[j - 1]);
                    b_row.push(b'-');
                    j -= 1;
                }
            } else {
                panic!("unexpected step");
            }

            // move to the next step in the alignment
            step = &self.grid[next_i][next_j];
        }

        let top: String = a_row.iter().rev().map(|c| *c as char).collect();
        let bottom: String = b_row.iter().rev().map(|c| *c as char).collect();
        (top, bottom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrices::NUC_4_4;

    #[test]
    fn test_aligner_align_small() {
        let mut a = Aligner::new(
            "GTT",
            "GAT",
            Scoring {
                replacement: NUC_4_4::MATRIX.clone(),
                gap_opening: -2,
                gap_extension: -2,
            },
        );
        let alignment = a.align();

        println!("{:?}", alignment);

        assert_eq!("G-T", alignment.a);
        assert_eq!("GAT", alignment.b);
    }

    /// same penalty for gap opening and extension
    #[test]
    fn test_aligner_align() {
        let mut a = Aligner::new(
            "TGTTACGG",
            "GGTTGACTA",
            Scoring {
                replacement: NUC_4_4::MATRIX.clone(),
                gap_opening: -2,
                gap_extension: -2,
            },
        );
        let alignment = a.align();

        println!("{:?}", alignment);

        assert_eq!("GTT-AC", alignment.a);
        assert_eq!("GTTGAC", alignment.b);
    }

    /// same gap extension as opening
    /// https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty_example
    #[test]
    fn test_aligner_align_same_gap_extension() {
        let mut a = Aligner::new(
            "TACGGGCCCGCTAC",
            "TAGCCCTATCGGTCA",
            Scoring {
                replacement: NUC_4_4::MATRIX.clone(),
                gap_opening: -1,
                gap_extension: -1,
            },
        );
        let alignment = a.align();

        println!("{:?}", alignment);

        assert_eq!("TACGGGCCCGCTA-C", alignment.a);
        assert_eq!("TA---G-CC-CTATC", alignment.b);
    }

    /// differing penalty for gap opening vs extension
    /// https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty_example
    #[test]
    fn test_aligner_align_small_gap_extension() {
        let mut a = Aligner::new(
            "TACGGGCCCGCTAC",
            "TAGCCCTATCGGTCA",
            Scoring {
                replacement: NUC_4_4::MATRIX.clone(),
                gap_opening: -5,
                gap_extension: -1,
            },
        );
        let alignment = a.align();

        println!("{:?}", alignment);

        assert_eq!("TACGGGCCCGCTA", alignment.a);
        assert_eq!("TA---GCC--CTA", alignment.b);
    }
}
