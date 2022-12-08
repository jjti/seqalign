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

use super::{Aligner, Scoring, Step};
use ordered_float::OrderedFloat;

/// A Needleman-Wunsch sequence aligner.
pub struct NeedlemanWunsch {
    scoring: Scoring,
}

impl Aligner for NeedlemanWunsch {
    fn scoring(&self) -> &Scoring {
        &self.scoring
    }

    fn default_grid_value(&self, index: usize) -> OrderedFloat<f32> {
        OrderedFloat(-(index as f32))
    }

    fn default_step_options(&self, _i: usize, _j: usize) -> Vec<Step> {
        vec![]
    }

    /// Needleman-Wunsch is a global alignment, so it starts at the very bottom right.
    fn backtrace_start(&self, grid: &[Vec<Step>]) -> Step {
        grid[grid.len() - 1][grid[0].len() - 1].clone()
    }
}

impl NeedlemanWunsch {
    pub fn new(scoring: Scoring) -> NeedlemanWunsch {
        NeedlemanWunsch { scoring }
    }
}

#[cfg(test)]
mod tests {
    use crate::matrices::MATCH;

    use super::*;

    #[test]
    fn test_aligner_align() {
        let a = NeedlemanWunsch::new(Scoring {
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
