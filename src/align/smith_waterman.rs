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

use ordered_float::OrderedFloat;

use super::{Align, Alignment, Scoring, Step};

struct Aligner {
    scoring: Scoring,
}

impl Align for Aligner {
    fn align<I: IntoIterator<Item = String>>(&self, seqs: I) -> Alignment {
        let mut input = seqs.into_iter();
        let a = input.next().unwrap();
        let b = input.next().unwrap();

        // Initialize the alignment grid.
        let grid = &mut self.init_grid(a.len(), b.len());

        // Fill in the alignment grid.
        self.fill(grid, a.as_bytes(), b.as_bytes());

        // Backtrace the grid to get the final alignment.
        self.backtrace(grid, a.as_bytes(), b.as_bytes())
    }
}

impl Aligner {
    pub fn new(scoring: Scoring) -> Aligner {
        Aligner { scoring }
    }

    fn init_grid(&self, a_len: usize, b_len: usize) -> Vec<Vec<Step>> {
        let mut grid: Vec<Vec<Step>> = Vec::new();

        for _i in 0..=b_len {
            // unlike needleman-wunsch, these rows start as 0 (rather than -i for a global alignment)
            grid.push(vec![Step::default(); a_len + 1]);
        }

        for i in 0..=b_len {
            grid[i][0] = Step::from(i, 0, 0f32);
        }

        for j in 0..=a_len {
            grid[0][j] = Step::from(0, j, 0f32);
        }

        grid
    }

    fn fill(&self, grid: &mut [Vec<Step>], a: &[u8], b: &[u8]) {
        for i in 1..=b.len() {
            for j in 1..=a.len() {
                // negative values are ignored, 0 is as low as we'll go
                let mut options: Vec<Step> = vec![Step {
                    val: OrderedFloat(0f32),
                    i,
                    j,
                    next: None,
                }];

                // gaps
                let mut k = 1;
                while k < i {
                    if a[j - 1] == b[k - 1] {
                        options.push(Step {
                            val: grid[k][j].val
                                + self.scoring.gap_opening
                                + self.scoring.gap_extension * (i - k - 1) as f32,
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
                        options.push(Step {
                            val: grid[i][l].val
                                + self.scoring.gap_opening
                                + self.scoring.gap_extension * (j - l - 1) as f32,
                            i,
                            j,
                            next: Some((i, l)),
                        });
                    }
                    l += 1;
                }

                // match or mismatch
                let match_val =
                    self.scoring.replacement[a[j - 1] as usize][b[i - 1] as usize] as f32;
                options.push(Step {
                    val: grid[i - 1][j - 1].val + match_val,
                    i,
                    j,
                    next: Some((i - 1, j - 1)),
                });

                grid[i][j] = options.iter().max().unwrap().clone();
            }
        }
    }

    fn backtrace(&self, grid: &mut Vec<Vec<Step>>, a: &[u8], b: &[u8]) -> Alignment {
        // this finds the global maximum among the alignments
        let mut step = &Step::default();
        for i in 0..grid.len() {
            for j in 0..grid[i].len() {
                if &grid[i][j] > step {
                    step = &grid[i][j];
                }
            }
        }
        let score = step.val.0;

        let mut alignment: Vec<Vec<char>> = vec![Vec::new(), Vec::new()];
        while let Some((next_i, next_j)) = step.next {
            let i_delta = step.i - next_i;
            let j_delta = step.j - next_j;

            if i_delta == 1 && j_delta == 1 {
                // match/mismatch
                alignment[0].push(a[step.j - 1] as char);
                alignment[1].push(b[step.i - 1] as char);
            } else if i_delta > 0 {
                // gap in seq a
                let mut i = step.i;
                while i > next_i {
                    alignment[0].push('-');
                    alignment[1].push(b[i - 1] as char);
                    i -= 1;
                }
            } else if j_delta > 0 {
                // gap in seq b
                let mut j = step.j;
                while j > next_j {
                    alignment[0].push(a[j - 1] as char);
                    alignment[1].push('-');
                    j -= 1;
                }
            } else {
                panic!("unexpected step");
            }

            // move to the next step in the alignment
            step = &grid[next_i][next_j];
        }

        for line in alignment.iter_mut() {
            line.reverse();
        }

        Alignment::new(alignment, grid.to_vec(), score)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrices::NUC_4_4;

    #[test]
    fn test_aligner_align_small() {
        let mut a = Aligner::new(Scoring {
            replacement: NUC_4_4::MATRIX,
            gap_opening: -2f32,
            gap_extension: -2f32,
        });
        let alignment = a.align(vec!["GTT".to_string(), "GAT".to_string()]);

        println!("{:?}", alignment);

        assert_eq!("G-T", alignment.rows[0].iter().collect::<String>());
        assert_eq!("GAT", alignment.rows[1].iter().collect::<String>());
    }

    /// same penalty for gap opening and extension
    #[test]
    fn test_aligner_align() {
        let a = Aligner::new(Scoring {
            replacement: NUC_4_4::MATRIX,
            gap_opening: -2f32,
            gap_extension: -2f32,
        });
        let alignment = a.align(vec!["TGTTACGG".to_string(), "GGTTGACTA".to_string()]);

        println!("{:?}", alignment);

        assert_eq!("GTT-AC", alignment.rows[0].iter().collect::<String>());
        assert_eq!("GTTGAC", alignment.rows[1].iter().collect::<String>());
    }

    /// same gap extension as opening
    /// https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty_example
    #[test]
    fn test_aligner_align_same_gap_extension() {
        let a = Aligner::new(Scoring {
            replacement: NUC_4_4::MATRIX,
            gap_opening: -1f32,
            gap_extension: -1f32,
        });
        let alignment = a.align(vec![
            "TACGGGCCCGCTAC".to_string(),
            "TAGCCCTATCGGTCA".to_string(),
        ]);

        println!("{:?}", alignment);

        assert_eq!(
            "TACGGGCCCGCTA-C",
            alignment.rows[0].iter().collect::<String>()
        );
        assert_eq!(
            "TA---G-CC-CTATC",
            alignment.rows[1].iter().collect::<String>()
        );
    }

    /// differing penalty for gap opening vs extension
    /// https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty_example
    #[test]
    fn test_aligner_align_small_gap_extension() {
        let a = Aligner::new(Scoring {
            replacement: NUC_4_4::MATRIX,
            gap_opening: -5f32,
            gap_extension: -1f32,
        });
        let alignment = a.align(vec![
            "TAGCCCTATCGGTCA".to_string(),
            "TACGGGCCCGCTAC".to_string(),
        ]);

        println!("{:?}", alignment);

        assert_eq!(
            "TA---GCC--CTA",
            alignment.rows[0].iter().collect::<String>()
        );
        assert_eq!(
            "TACGGGCCCGCTA",
            alignment.rows[1].iter().collect::<String>()
        );
    }
}
