//! Expressed with the big O notation commonly used to measure computational complexity,
//! a naïve MSA takes O(LengthNseqs) time to produce. To find the global optimum for n
//! sequences this way has been shown to be an NP-complete problem.
//!
//! In 1989, based on Carrillo-Lipman Algorithm,[8] Altschul introduced a practical method
//! that uses pairwise alignments to constrain the n-dimensional search space.[9] In this
//! approach pairwise dynamic programming alignments are performed on each pair of sequences
//! in the query set, and only the space near the n-dimensional intersection of these
//! alignments is searched for the n-way alignment. The MSA program optimizes the sum of
//! all of the pairs of characters at each position in the alignment (the so-called sum of
//! pair score) and has been implemented in a software program for constructing multiple
//! sequence alignments.
//!
//! The Multiple Sequence Alignment Problem in Biology:
//! https://zenodo.org/record/1236134
//!
//! A tool for multiple sequence alignment:
//! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC287279/?page=1
