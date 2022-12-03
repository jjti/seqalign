//! Clustal-W paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC308517/pdf/nar00046-0131.pdf
//! cites neighbor joining method paper: https://pubmed.ncbi.nlm.nih.gov/3447015/
//!
//! unweighted pair group method with arithmetic mean
//! https://en.wikipedia.org/wiki/UPGMA
use crate::matrices::BLOSUM62;

use super::{needleman_wunsch::Aligner, PWAlign};

fn clustal_w(seqs: &Vec<String>) {}

fn upgma(seqs: &Vec<String>) {
    // create initial matrix of initial pairwise alignments
    let mut distances: Vec<Vec<f32>> = Vec::new();
    for i in 0..seqs.len() {
        distances.push(Vec::new());

        for j in 0..i {
            distances[i][j] = Aligner::new(
                seqs.get(i).unwrap(),
                seqs.get(j).unwrap(),
                super::Scoring {
                    replacement: BLOSUM62::MATRIX,
                    gap_opening: -2,
                    gap_extension: -1,
                },
            )
            .align()
            .distance();
        }
    }

    // start iterating over the matrix, finding, clustering the best pairs
    // a Cluster is a tuple of Node = (Node | String, Node | String)
    // distance is maintained in the distances matrix
    // while there are still distances in the matrix
    let clusters: Vec<Node> = Vec::new();
    while distances.len() > 2 {
        // find the minimum distance pair
        let mut min_pair: (usize, usize) = (0, 0);
        for i in 0..seqs.len() {
            for j in 0..1 {
                if distances[i][j] < distances[min_pair.0][min_pair.1] {
                    min_pair = (i, j);
                }
            }
        }

        //
    }
}

// Box::new(Nil)
struct Node {
    left_seq: Option<String>,
    left_node: Box<Node>,
    left_distance: f32,

    right_seq: Option<String>,
    right_node: Box<Node>,
    right_distance: f32,
}
