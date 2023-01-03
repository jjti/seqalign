use std::{
    collections::HashSet,
    fmt::{Debug, Display},
};

use crate::matrices::Matrix;

use super::{step::Step, strategy::Strategy};

pub struct Alignment {
    /// the 2D-grid holding the alignment of the sequences.
    pub rows: Vec<Vec<char>>,

    /// steps are the steps from the full 2D alignment.
    pub steps: Vec<Vec<Step>>,

    /// score is the final alignment score.
    pub score: f32,

    /// distance is the ratio of characters that differ between the two sequences.
    pub distance: f32,
}

impl Alignment {
    pub fn new(alignment: Vec<Vec<char>>, steps: Vec<Vec<Step>>, score: f32) -> Self {
        if alignment.len() < 2 {
            panic!("Alignment must have at least one row")
        }

        let distance = Alignment::calc_distance(&alignment);

        Alignment {
            rows: alignment,
            steps,
            score,
            distance,
        }
    }

    /// returns the ratio of residues that have a difference between the
    /// sequences in the alignment. Residues that are entirely gap are ignored.
    fn calc_distance(alignment: &Vec<Vec<char>>) -> f32 {
        let mut residues = 0f32;
        let mut diffs = 0f32;
        for col in 0..alignment[0].len() {
            let chars = (0..alignment.len())
                .map(|i| {
                    if col < alignment[i].len() {
                        alignment[i][col]
                    } else {
                        '-'
                    }
                })
                .collect::<HashSet<_>>();
            if chars.len() == 1 {
                if chars.contains(&'-') {
                    continue; // it's all gap, skip
                }
            } else {
                diffs += 1f32;
            }
            residues += 1f32;
        }

        diffs / residues
    }
}

impl Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.rows
                .iter()
                .map(|a| a.iter().collect::<String>())
                .collect::<Vec<_>>()
                .join("\n")
        )
    }
}

// different formatting traits for different formatting types:
// https://doc.rust-lang.org/std/fmt/index.html#formatting-traits
//
// fmt::Debug implementations should be implemented for all public types.
impl Debug for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.rows.len() < 2 {
            panic!("Alignment must have at least two rows")
        }

        let header = format!(
            "{}\n{}\n",
            self.rows[0].iter().collect::<String>(),
            self.rows[1].iter().collect::<String>()
        );
        let mut result: Vec<String> = vec![header];

        let a: Vec<char> = self.rows[0]
            .clone()
            .into_iter()
            .filter(|c| *c != '-')
            .collect();
        let b: Vec<char> = self.rows[1]
            .clone()
            .into_iter()
            .filter(|c| *c != '-')
            .collect();

        for i in 0..b.len() + 1 {
            // Write seq A header row.
            if i == 0 {
                result.push("   |    |".to_string());
                a.iter().for_each(|c| result.push(format!(" {: <3}|", c)));
                result.push("\n".to_string());
            }

            // Write first column of grid.
            if i == 0 {
                result.push("   |".to_string());
            } else {
                result.push(format!("{: <3}|", b[i - 1]));
            }

            // Write each character of the grid.
            self.steps[i]
                .iter()
                .for_each(|f| result.push(format!(" {: <3}|", f.val)));
            result.push("\n".to_string());
        }

        write!(f, "{}", result.join(""))
    }
}

#[derive(Debug)]
pub struct Scoring {
    /// replacement matrix
    pub matrix: Matrix,

    /// penalty for a gap opening
    pub gap_opening: f32,

    /// penalty for a gap extension
    pub gap_extension: f32,
}

/// align sequences using a strategy and scoring scheme.
pub fn align<I: IntoIterator<Item = String>>(
    seqs: I,
    strategy: &Strategy,
    scoring: &Scoring,
) -> Alignment {
    let mut input = seqs.into_iter();
    let a = input.next().unwrap();
    let b = input.next().unwrap();

    // Initialize the alignment grid.
    let grid = &mut init_grid(strategy, a.len(), b.len());

    // Fill in the alignment grid.
    fill_grid(strategy, scoring, grid, a.as_bytes(), b.as_bytes());

    // Backtrace the grid to get the final alignment.
    backtrace(strategy, grid, a.as_bytes(), b.as_bytes())
}

/// init grid sets up the 2D alignment grid.
fn init_grid(strategy: &Strategy, a_len: usize, b_len: usize) -> Vec<Vec<Step>> {
    let mut grid: Vec<Vec<Step>> = Vec::new();
    for i in 0..b_len + 1 {
        grid.push(vec![Step::default(); a_len + 1]);

        grid[i][0] = Step {
            i,
            j: 0,
            val: (strategy.init_grid_value)(i),
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
            val: (strategy.init_grid_value)(j),
            next: match j {
                0 => None,
                _ => Some((0, j - 1)),
            },
        };
    }

    grid
}

/// fill_grid traverses the grid from top-left to bottom right, filling in each step.
fn fill_grid(strategy: &Strategy, scoring: &Scoring, grid: &mut [Vec<Step>], a: &[u8], b: &[u8]) {
    for i in 1..=b.len() {
        for j in 1..=a.len() {
            // negative values are ignored, 0 is as low as we'll go
            let mut options: Vec<Step> = (strategy.init_step_options)(i, j);

            // gaps
            let mut k = 1;
            while k < i {
                if a[j - 1] == b[k - 1] {
                    options.push(Step {
                        val: grid[k][j].val
                            + scoring.gap_opening
                            + scoring.gap_extension * (i - k - 1) as f32,
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
                            + scoring.gap_opening
                            + scoring.gap_extension * (j - l - 1) as f32,
                        i,
                        j,
                        next: Some((i, l)),
                    });
                }
                l += 1;
            }

            // match or mismatch
            let match_val = scoring.matrix[a[j - 1] as usize][b[i - 1] as usize] as f32;
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

/// backtrace walks to the result of alignment to find the optimal alignment.
fn backtrace(strategy: &Strategy, grid: &mut [Vec<Step>], a: &[u8], b: &[u8]) -> Alignment {
    let mut step = &(strategy.init_backtrace)(grid);
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

#[cfg(test)]
mod tests {
    use crate::align::step::Step;

    use super::*;

    #[test]
    fn test_alignment_debug() {
        let alignment = Alignment::new(
            vec![vec!['A', 'G', 'C'], vec!['C', 'G']],
            vec![
                vec![
                    Step::from(0, 0, 0f32),
                    Step::from(0, 1, -1f32),
                    Step::from(0, 2, -2f32),
                    Step::from(0, 3, -3f32),
                ],
                vec![
                    Step::from(1, 0, -1f32),
                    Step::from(1, 1, 0f32),
                    Step::from(1, 2, 0f32),
                    Step::from(1, 3, 1f32),
                ],
                vec![
                    Step::from(2, 0, -2f32),
                    Step::from(2, 1, 0f32),
                    Step::from(2, 2, 1f32),
                    Step::from(2, 3, 1f32),
                ],
            ],
            0f32,
        );

        assert_eq!(
            "AGC
CG
   |    | A  | G  | C  |
   | 0  | -1 | -2 | -3 |
C  | -1 | 0  | 0  | 1  |
G  | -2 | 0  | 1  | 1  |
",
            format!("{:?}", alignment)
        )
    }

    #[test]
    fn test_alignment_distance() {
        let alignment = Alignment::new(
            vec![vec!['A', 'C', 'C', 'G', 'T'], vec!['A', 'G', '-', 'C', 'T']],
            vec![vec![]],
            0f32,
        );

        assert_eq!(0.6, alignment.distance)
    }
}
