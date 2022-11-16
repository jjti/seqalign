//! Smith-Waterman algorithm
//! https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
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



