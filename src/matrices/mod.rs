#![no_implicit_prelude]

/// Matrix is a single alignment matrix used in scoring an alignment.
///
/// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = [[i32; 128]; 128];

#[allow(non_snake_case)]
pub mod BLOSUM62;
#[allow(non_snake_case)]
pub mod MATCH;
#[allow(non_snake_case)]
pub mod NUC_4_4;
