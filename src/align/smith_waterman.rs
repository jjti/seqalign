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

use super::{step::Step, strategy::Strategy};

/// strategy for the Smith-Waterman algorithm.
pub const STRATEGY: Strategy = Strategy {
    init_grid_value: |_i: usize| -> OrderedFloat<f32> { OrderedFloat(0f32) },

    init_step_options: |i: usize, j: usize| -> Vec<Step> {
        vec![Step {
            i,
            j,
            val: OrderedFloat(0f32),
            next: None,
        }]
    },

    init_backtrace: |grid: &[Vec<Step>]| -> Step {
        let mut step = grid[0][0].clone();
        for (i, _) in grid.iter().enumerate() {
            for j in 0..grid[i].len() {
                if grid[i][j] > step {
                    step = grid[i][j].clone();
                }
            }
        }
        step
    },
};

#[cfg(test)]
mod tests {
    use crate::{
        align::{align, smith_waterman::STRATEGY, Scoring},
        matrices::NUC_4_4,
    };

    #[test]
    fn test_aligner_align_small() {
        let alignment = align(
            vec!["GTT".to_string(), "GAT".to_string()],
            &STRATEGY,
            &Scoring {
                matrix: NUC_4_4::MATRIX,
                gap_opening: -2f32,
                gap_extension: -2f32,
            },
        );

        println!("{:?}", alignment);

        assert_eq!("G-T", alignment.rows[0].iter().collect::<String>());
        assert_eq!("GAT", alignment.rows[1].iter().collect::<String>());
    }

    /// same penalty for gap opening and extension
    #[test]
    fn test_aligner_align() {
        let alignment = align(
            vec!["TGTTACGG".to_string(), "GGTTGACTA".to_string()],
            &STRATEGY,
            &Scoring {
                matrix: NUC_4_4::MATRIX,
                gap_opening: -2f32,
                gap_extension: -2f32,
            },
        );

        println!("{:?}", alignment);

        assert_eq!("GTT-AC", alignment.rows[0].iter().collect::<String>());
        assert_eq!("GTTGAC", alignment.rows[1].iter().collect::<String>());
    }

    /// same gap extension as opening
    /// https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty_example
    #[test]
    fn test_aligner_align_same_gap_extension() {
        let alignment = align(
            vec!["TACGGGCCCGCTAC".to_string(), "TAGCCCTATCGGTCA".to_string()],
            &STRATEGY,
            &Scoring {
                matrix: NUC_4_4::MATRIX,
                gap_opening: -1f32,
                gap_extension: -1f32,
            },
        );

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
        let alignment = align(
            vec!["TAGCCCTATCGGTCA".to_string(), "TACGGGCCCGCTAC".to_string()],
            &STRATEGY,
            &Scoring {
                matrix: NUC_4_4::MATRIX,
                gap_opening: -5f32,
                gap_extension: -1f32,
            },
        );

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
