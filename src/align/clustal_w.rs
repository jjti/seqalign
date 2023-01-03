//! Clustal-W paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC308517/pdf/nar00046-0131.pdf
//! cites neighbor joining method paper: https://pubmed.ncbi.nlm.nih.gov/3447015/
//!
//! unweighted pair group method with arithmetic mean
//! https://en.wikipedia.org/wiki/UPGMA

use std::collections::{HashMap, HashSet};

use super::{align, needleman_wunsch, Alignment};
use crate::matrices::BLOSUM62;

type Distances = HashMap<(usize, usize), f32>;

#[derive(Clone)]
struct Node {
    index: usize,

    /// either seq is set OR this is truly a node
    seq: Option<String>,

    /// or left and right are set
    left: Option<Box<Node>>,
    left_branch_len: f32,
    right: Option<Box<Node>>,
    right_branch_len: f32,
}

impl Node {
    fn from_seq(index: usize, s: String) -> Node {
        Node {
            index,
            seq: Some(s),
            left: None,
            left_branch_len: 0f32,
            right: None,
            right_branch_len: 0f32,
        }
    }

    fn from_nodes(index: usize, left: &Node, right: &Node, dist: f32) -> Node {
        Node {
            index,
            seq: None,
            left: Some(Box::new(left.clone())),
            left_branch_len: (dist / 2f32) - left.left_branch_len,
            right: Some(Box::new(right.clone())),
            right_branch_len: (dist / 2f32) - right.right_branch_len,
        }
    }

    fn is_leaf(&self) -> bool {
        self.seq.is_some()
    }
}

/// upgma
///
/// Put each Seq in its own Node (by itself)
/// Maintain a clusters Vec that holds all remaining clusters
///   - Each sequence is initially stuck in its own cluster
/// Maintain a distances HashMap that maps pairs of clusters to their distance
///
/// On each iteration:
///  - Take the cluster in distances with the minimum value
///  - Merge the two clusters in that pair
///  - Remove the pair from distances
fn upgma(seqs: &mut [String]) -> Vec<Node> {
    let mut clusters: Vec<Node> = seqs
        .iter()
        .enumerate()
        .map(|(i, s)| Node::from_seq(i, s.to_string()))
        .collect();
    let mut unpaired: HashSet<usize> = (0..clusters.len()).collect(); // unpaired seqs

    // Initialize distances
    let mut distances: Distances = HashMap::new();
    for i in 0..clusters.len() {
        for j in i + 1..clusters.len() {
            let alignment = align(
                vec![seqs[i].clone(), seqs[j].clone()],
                &needleman_wunsch::STRATEGY,
                &super::Scoring {
                    matrix: BLOSUM62::MATRIX,
                    gap_opening: -1f32,
                    gap_extension: -0.5f32,
                },
            );

            distances.insert((i, j), alignment.distance);
        }
    }

    while !distances.is_empty() {
        // find minimum distance in distances
        let mut min = f32::MAX;
        let mut min_pair: (usize, usize) = (0, 0);
        for ((i, j), d) in distances.iter() {
            if *d < min {
                min = *d;
                min_pair = (*i, *j);
            }
        }

        // create a new cluster from the two in the minimum distance pair
        let left = &clusters[min_pair.0].clone();
        let right = &clusters[min_pair.1].clone();
        clusters.push(Node::from_nodes(clusters.len(), left, right, min));

        // remove the two clusters from unpaired
        unpaired.remove(&min_pair.0);
        unpaired.remove(&min_pair.1);

        // updated distances
        distances = new_distances(
            distances,
            &unpaired,
            min_pair.0,
            min_pair.1,
            clusters.len() - 1,
        );
    }

    clusters
}

/// new_distances removes distances from the right/left clustered sequences and
/// returns a new Distances map with updated distances to the new cluster.
fn new_distances(
    distances_old: Distances,
    unpaired: &HashSet<usize>,
    left: usize,
    right: usize,
    new: usize,
) -> Distances {
    // create a new distances map using the newly created cluster
    let mut distances_new: Distances = HashMap::new();

    // keep all unrelated distances
    distances_new.extend(
        distances_old
            .iter()
            .filter(|((i, j), _)| *i != left && *j != left && *i != right && *j != right),
    );

    // update the distances of all unpaired clusters to the new cluster
    for unpaired_cluster in unpaired.iter() {
        let left_key: (usize, usize) = (*unpaired_cluster.min(&left), *unpaired_cluster.max(&left));
        let right_key: (usize, usize) =
            (*unpaired_cluster.min(&right), *unpaired_cluster.max(&right));

        distances_new.insert(
            (*unpaired_cluster, new),
            (distances_old[&left_key] + distances_old[&right_key]) / 2f32,
        );
    }

    distances_new
}

impl Alignment {
    // align combines this alignment with another one. Gaps are preserved.
    //
    // This matches the process described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC308517/pdf/nar00046-0131.pdf
    // fn align(&self, other: Alignment) -> Alignment {
    //     let mut alignment = Alignment {
    //         seqs: vec![self.seqs[0].clone(), other.seqs[0].clone()],
    //         distance: 0f32,
    //     };

    //     for i in 0..self.seqs[0].len() {
    //         if self.seqs[0][i] == other.seqs[0][i] {
    //             alignment.seqs[0].push(self.seqs[0][i]);
    //             alignment.seqs[1].push(other.seqs[0][i]);
    //         } else {
    //             alignment.seqs[0].push('-');
    //             alignment.seqs[1].push('-');
    //             alignment.distance += 1f32;
    //         }
    //     }

    //     alignment
    // }

    // /// init_grid creates a new 2D alignment grid with negatives values for each initial row.
    // fn init_grid(&self, a_len: usize, b_len: usize) -> Vec<Vec<Step>> {
    //     let mut grid: Vec<Vec<Step>> = Vec::new();
    //     for i in 0..b_len + 1 {
    //         grid.push(vec![Step::default(); a_len + 1]);

    //         grid[i][0] = Step {
    //             i,
    //             j: 0,
    //             val: OrderedFloat(-(i as f32)),
    //             next: match i {
    //                 0 => None,
    //                 _ => Some((i - 1, 0)),
    //             },
    //         };
    //     }

    //     for j in 0..a_len + 1 {
    //         grid[0][j] = Step {
    //             i: 0,
    //             j,
    //             val: OrderedFloat(-(j as f32)),
    //             next: match j {
    //                 0 => None,
    //                 _ => Some((0, j - 1)),
    //             },
    //         };
    //     }

    //     grid
    // }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_upgma() {
        let clusters = upgma(&mut [
            "ACGTA".to_string(),
            "ACGCA".to_string(),
            "ACTTA".to_string(),
            "ACATA".to_string(),
        ]);

        // assert_eq!(7, clusters.len())
    }
}
