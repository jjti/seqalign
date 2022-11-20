use crate::io::matrix::Matrix;
use std::fmt::{Debug, Display};

pub mod needleman_wunsch;
pub mod smith_waterman;

pub trait Align {
    fn align(&mut self) -> Alignment;
}

pub struct Alignment {
    grid: Vec<Vec<i32>>,
    a: String,
    b: String,
    a_orig: String,
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
                .for_each(|f| result.push(format!(" {: <3}|", f)));
            result.push("\n".to_string());
        }

        write!(f, "{}", result.join(""))
    }
}

pub struct Scoring {
    /// replacement matrix
    m: Matrix,

    /// insertion or deletion
    indel: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_debug() {
        let a = Alignment {
            grid: vec![vec![0, -1, -2, -3], vec![-1, 0, 0, 1], vec![-2, 0, 1, 1]],
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
