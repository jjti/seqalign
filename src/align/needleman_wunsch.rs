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

use super::{step::Step, strategy::Strategy};
use ordered_float::OrderedFloat;

/// strategy for the Needleman-Wunsch algorithm.
pub const STRATEGY: Strategy = Strategy {
    init_grid_value: |index: usize| -> OrderedFloat<f32> { OrderedFloat(-(index as f32)) },

    init_step_options: |_i: usize, _j: usize| -> Vec<Step> { vec![] },

    init_backtrace: |grid: &[Vec<Step>]| -> Step {
        grid[grid.len() - 1][grid[0].len() - 1].clone()
    },
};

#[cfg(test)]
mod tests {
    use crate::{
        align::{align, Scoring},
        matrices::MATCH,
    };

    use super::*;

    #[test]
    fn test_aligner_align() {
        let alignment = align(
            vec!["GCATGCG".to_string(), "GATTACA".to_string()],
            &STRATEGY,
            &Scoring {
                matrix: MATCH::MATRIX,
                gap_opening: -1f32,
                gap_extension: -1f32,
            },
        );

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
