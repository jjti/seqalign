//! Clustal-W paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC308517/pdf/nar00046-0131.pdf
//! cites neighbor joining method paper: https://pubmed.ncbi.nlm.nih.gov/3447015/
//!
//! unweighted pair group method with arithmetic mean
//! https://en.wikipedia.org/wiki/UPGMA
use crate::matrices::BLOSUM62;

use super::{needleman_wunsch::Aligner, PWAlign, PWAlignment};

fn clustal_w(seqs: &Vec<String>) {
    // create a matrix of initial pairwise alignments
    let mut alignments: Vec<Vec<PWAlignment>> = Vec::new();
    for i in 0..seqs.len() {
        alignments.push(Vec::new());

        for j in 0..i {
            alignments[i][j] = Aligner::new(
                seqs.get(i).unwrap(),
                seqs.get(j).unwrap(),
                super::Scoring {
                    replacement: BLOSUM62::MATRIX,
                    gap_opening: -2,
                    gap_extension: -1,
                },
            )
            .align();
        }
    }

    // create a distance matrix between each seq in the pairwise alignment
    // let mut distances: Vec<Vec<i32>> = Vec::new();
    // let mut min_pair: (usize, usize) = (0, 0);
    // for i in 0..alignments.len() {
    //     for j in 0..i {
    //         if alignments[i][j].distance() < alignments[min_pair.0][min_pair.1].distance() {
    //             min_val = alignments[i][j].distance();
    //         }
    //     }
    // }
}
