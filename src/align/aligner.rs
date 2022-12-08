use super::{step::Step, Alignment};
use crate::matrices::Matrix;
use ordered_float::OrderedFloat;

#[derive(Debug)]
pub struct Scoring {
    /// replacement matrix
    pub replacement: Matrix,

    /// penalty for a gap opening
    pub gap_opening: f32,

    /// penalty for a gap extension
    pub gap_extension: f32,
}

pub trait Aligner {
    /// get_scoring returns the alignment scoring to use.
    fn scoring(&self) -> &Scoring;

    /// init_grid_value defines the initial value of a cell in the alignment grid along the edge.
    ///
    /// For Needleman-Wunsch, this is a negative value proportional to the index.
    /// For Smith-Waterman, this is 0.
    fn default_grid_value(&self, _i: usize) -> OrderedFloat<f32>;

    /// default_step_options are those that are consistent across each set of options in an alignment.
    fn default_step_options(&self, _i: usize, _j: usize) -> Vec<Step>;

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

    /// init grid sets up the 2D alignment grid.
    fn init_grid(&self, a_len: usize, b_len: usize) -> Vec<Vec<Step>> {
        let mut grid: Vec<Vec<Step>> = Vec::new();
        for i in 0..b_len + 1 {
            grid.push(vec![Step::default(); a_len + 1]);

            grid[i][0] = Step {
                i,
                j: 0,
                val: self.default_grid_value(i),
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
                val: self.default_grid_value(j),
                next: match j {
                    0 => None,
                    _ => Some((0, j - 1)),
                },
            };
        }

        grid
    }

    /// fill_grid traverses the grid from top-left to bottom right, filling in each step.
    fn fill_grid(&self, grid: &mut [Vec<Step>], a: &[u8], b: &[u8]) {
        let scoring = self.scoring();

        for i in 1..=b.len() {
            for j in 1..=a.len() {
                // negative values are ignored, 0 is as low as we'll go
                let mut options: Vec<Step> = self.default_step_options(i, j);

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
                let match_val = scoring.replacement[a[j - 1] as usize][b[i - 1] as usize] as f32;
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

    /// backtrace_start finds the starting point for backtracing.
    fn backtrace_start(&self, grid: &[Vec<Step>]) -> Step {
        // this finds the global maximum among the alignments
        let mut step = Step::default();
        for i in 0..grid.len() {
            for j in 0..grid[i].len() {
                if grid[i][j] > step {
                    step = grid[i][j].clone();
                }
            }
        }
        step
    }

    /// backtrace walks to the result of alignment to find the optimal alignment.
    fn backtrace(&self, grid: &mut Vec<Vec<Step>>, a: &[u8], b: &[u8]) -> Alignment {
        let mut step = &self.backtrace_start(grid);
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
