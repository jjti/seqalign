pub use crate::align::alignment::align;
pub use crate::align::alignment::Alignment;
pub use crate::align::alignment::Scoring;
pub use crate::align::strategy::Method;

mod alignment;
mod clustal_w;
mod needleman_wunsch;
mod smith_waterman;
mod step;
mod strategy;
