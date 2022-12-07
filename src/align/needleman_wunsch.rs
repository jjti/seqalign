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

use ordered_float::OrderedFloat;

use super::{Align, Alignment, Scoring, Step};

/// A Needleman-Wunsch sequence aligner.
pub struct Aligner {
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
        self.fill_grid(grid, a.as_bytes(), b.as_bytes());

        // Backtrace the grid to get the final alignment.
        self.backtrace(grid, a.as_bytes(), b.as_bytes())
    }
}

impl Aligner {
    pub fn new(scoring: Scoring) -> Aligner {
        Aligner { scoring }
    }

    /// init_grid creates a new 2D alignment grid with negatives values for each initial row.
    fn init_grid(&self, a_len: usize, b_len: usize) -> Vec<Vec<Step>> {
        let mut grid: Vec<Vec<Step>> = Vec::new();
        for i in 0..b_len + 1 {
            grid.push(vec![Step::default(); a_len + 1]);

            grid[i][0] = Step {
                i,
                j: 0,
                val: OrderedFloat(-(i as f32)),
                next: match i {
                    0 => None,
                    _ => Some((i - 1, 0)),
                },
            };
        }

        for j in 0..a_len + 1 {
            grid[0][j] = Step {
                i: 0,
                j,
                val: OrderedFloat(-(j as f32)),
                next: match j {
                    0 => None,
                    _ => Some((0, j - 1)),
                },
            };
        }

        grid
    }

    /// fill_grid in the alignment grid.
    fn fill_grid(&self, grid: &mut [Vec<Step>], a: &[u8], b: &[u8]) {
        for i in 1..b.len() + 1 {
            for j in 1..a.len() + 1 {
                let mut options: Vec<Step> = vec![
                    // indel
                    Step {
                        i,
                        j,
                        val: grid[i - 1][j].val + self.scoring.gap_opening,
                        next: Some((i - 1, j)),
                    },
                    Step {
                        i,
                        j,
                        val: grid[i][j - 1].val + self.scoring.gap_opening,
                        next: Some((i, j - 1)),
                    },
                ];

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

    /// backtrace builds up the the final alignment parsing it from the grid.
    fn backtrace(&self, grid: &mut [Vec<Step>], a: &[u8], b: &[u8]) -> Alignment {
        let mut alignment: Vec<Vec<char>> = vec![vec![], vec![]];

        let mut step = &grid[grid.len() - 1][a.len()];
        let score = step.val.0;
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
    use crate::matrices::MATCH;

    use super::*;

    #[test]
    fn test_aligner_align() {
        let a = Aligner::new(Scoring {
            replacement: MATCH::MATRIX,
            gap_opening: -1f32,
            gap_extension: -1f32,
        });
        let alignment = a.align(vec!["GCATGCG".to_string(), "GATTACA".to_string()]);

        println!("{:?}", alignment);

        assert_eq!(
            "GCA-TGCG",
            alignment.rows[0].clone().into_iter().collect::<String>()
        );
        assert_eq!(
            "G-ATTACA",
            alignment.rows[1].clone().into_iter().collect::<String>()
        );
    }
}
