use ordered_float::OrderedFloat;

use super::{needleman_wunsch, smith_waterman, step::Step};

/// Strategy defines the strategy for aligning two sequences.
pub struct Strategy {
    /// init_grid_value defines the initial value of a cell in the alignment grid along the edge.
    ///
    /// For Needleman-Wunsch, this is a negative value proportional to the index.
    /// For Smith-Waterman, this is 0.
    pub init_grid_value: fn(_i: usize) -> OrderedFloat<f32>,

    /// init_step_options defines the default options for a cell in the alignment grid.
    pub init_step_options: fn(_i: usize, _j: usize) -> Vec<Step>,

    /// init_backtrace finds the initial Step for backtracing the alignment grid.
    pub init_backtrace: fn(grid: &[Vec<Step>]) -> Step,
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Method {
    NeedlemanWunsch,
    SmithWaterman,
}

impl Method {
    pub fn strategy(&self) -> &Strategy {
        match self {
            Method::NeedlemanWunsch => &needleman_wunsch::STRATEGY,
            Method::SmithWaterman => &smith_waterman::STRATEGY,
        }
    }
}
