use crate::io::matrix::Matrix;
use std::fmt::{Debug, Display};

pub mod needleman_wunsch;
pub mod smith_waterman;

pub trait Align {
    fn align(&mut self) -> Alignment;
}

/// Step is a single step through an alignment.
#[derive(Clone, Eq)]
pub struct Step {
    val: i32,
    i: usize,
    j: usize,
    next: Option<(usize, usize)>,
}

impl Step {
    fn from(i: usize, j: usize, val: i32) -> Self {
        Step {
            val,
            i,
            j,
            next: None,
        }
    }
}

impl Default for Step {
    fn default() -> Self {
        Step {
            val: i32::MIN,
            i: usize::MIN,
            j: usize::MIN,
            next: None,
        }
    }
}

impl Debug for Step {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {}): {}", self.i, self.j, self.val)
    }
}

impl PartialEq for Step {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j && self.val == other.val
    }
}

impl Ord for Step {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.val.cmp(&other.val)
    }
}

impl PartialOrd for Step {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.val.partial_cmp(&other.val)
    }
}

pub struct Alignment {
    /// grid holds the alignment of the two sequences
    grid: Vec<Vec<Step>>,

    /// a is the top sequence of the alignment
    a: String,

    /// b is the bottom sequence of the alignment
    b: String,

    /// a_orig is the original top input sequence
    a_orig: String,

    /// b_orig is the bottom input sequence
    b_orig: String,
}

impl Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\n{}", self.a, self.b)
    }
}

// different formatting traits for different formatting types:
// https://doc.rust-lang.org/std/fmt/index.html#formatting-traits
//
// fmt::Debug implementations should be implemented for all public types.
impl Debug for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let header = format!("{}\n{}\n", self.a, self.b);
        let mut result: Vec<String> = vec![header];

        for i in 0..self.b_orig.len() + 1 {
            // Write seq A header row.
            if i == 0 {
                result.push("   |    |".to_string());
                self.a_orig
                    .chars()
                    .for_each(|c| result.push(format!(" {: <3}|", c)));
                result.push("\n".to_string());
            }

            // Write first column of grid.
            if i == 0 {
                result.push("   |".to_string());
            } else {
                result.push(format!("{: <3}|", self.b_orig.chars().nth(i - 1).unwrap()));
            }

            // Write each character of the grid.
            self.grid[i]
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
    replacement: Matrix,

    /// penalty for a gap opening
    gap_opening: i32,

    /// penalty for a gap extension
    gap_extension: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_debug() {
        let a = Alignment {
            grid: vec![
                vec![
                    Step::from(0, 0, 0),
                    Step::from(0, 1, -1),
                    Step::from(0, 2, -2),
                    Step::from(0, 3, -3),
                ],
                vec![
                    Step::from(1, 0, -1),
                    Step::from(1, 1, 0),
                    Step::from(1, 2, 0),
                    Step::from(1, 3, 1),
                ],
                vec![
                    Step::from(2, 0, -2),
                    Step::from(2, 1, 0),
                    Step::from(2, 2, 1),
                    Step::from(2, 3, 1),
                ],
            ],
            a: "AGC".to_string(),
            b: "CG".to_string(),
            a_orig: "AGC".to_string(),
            b_orig: "CG".to_string(),
        };

        assert_eq!(
            "AGC
CG
   |    | A  | G  | C  |
   | 0  | -1 | -2 | -3 |
C  | -1 | 0  | 0  | 1  |
G  | -2 | 0  | 1  | 1  |
",
            format!("{:?}", a)
        )
    }
}
