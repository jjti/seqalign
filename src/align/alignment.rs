use std::{
    collections::HashSet,
    fmt::{Debug, Display},
};

use super::Step;

pub struct Alignment {
    /// the 2D-grid holding the alignment of the sequences.
    pub rows: Vec<Vec<char>>,

    /// steps are the steps from the full 2D alignment.
    pub steps: Vec<Vec<Step>>,

    /// score is the final alignment score.
    pub score: f32,

    /// distance is the ratio of characters that differ between the two sequences.
    pub distance: f32,
}

impl Alignment {
    pub fn new(alignment: Vec<Vec<char>>, steps: Vec<Vec<Step>>, score: f32) -> Self {
        if alignment.len() < 2 {
            panic!("Alignment must have at least one row")
        }

        let distance = Alignment::calc_distance(&alignment);

        Alignment {
            rows: alignment,
            steps,
            score,
            distance,
        }
    }

    /// returns the ratio of residues that have a difference between the
    /// sequences in the alignment. Residues that are entirely gap are ignored.
    fn calc_distance(alignment: &Vec<Vec<char>>) -> f32 {
        let mut residues = 0f32;
        let mut diffs = 0f32;
        for col in 0..alignment[0].len() {
            let chars = (0..alignment.len())
                .map(|i| {
                    if col < alignment[i].len() {
                        alignment[i][col]
                    } else {
                        '-'
                    }
                })
                .collect::<HashSet<_>>();
            if chars.len() == 1 {
                if chars.contains(&'-') {
                    continue; // it's all gap, skip
                }
            } else {
                diffs += 1f32;
            }
            residues += 1f32;
        }

        diffs / residues
    }
}

impl Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.rows
                .iter()
                .map(|a| a.iter().collect::<String>())
                .collect::<Vec<_>>()
                .join("\n")
        )
    }
}

// different formatting traits for different formatting types:
// https://doc.rust-lang.org/std/fmt/index.html#formatting-traits
//
// fmt::Debug implementations should be implemented for all public types.
impl Debug for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.rows.len() < 2 {
            panic!("Alignment must have at least two rows")
        }

        let header = format!(
            "{}\n{}\n",
            self.rows[0].iter().collect::<String>(),
            self.rows[1].iter().collect::<String>()
        );
        let mut result: Vec<String> = vec![header];

        let a: Vec<char> = self.rows[0]
            .clone()
            .into_iter()
            .filter(|c| *c != '-')
            .collect();
        let b: Vec<char> = self.rows[1]
            .clone()
            .into_iter()
            .filter(|c| *c != '-')
            .collect();

        for i in 0..b.len() + 1 {
            // Write seq A header row.
            if i == 0 {
                result.push("   |    |".to_string());
                a.iter().for_each(|c| result.push(format!(" {: <3}|", c)));
                result.push("\n".to_string());
            }

            // Write first column of grid.
            if i == 0 {
                result.push("   |".to_string());
            } else {
                result.push(format!("{: <3}|", b[i - 1]));
            }

            // Write each character of the grid.
            self.steps[i]
                .iter()
                .for_each(|f| result.push(format!(" {: <3}|", f.val)));
            result.push("\n".to_string());
        }

        write!(f, "{}", result.join(""))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_debug() {
        let alignment = Alignment::new(
            vec![vec!['A', 'G', 'C'], vec!['C', 'G']],
            vec![
                vec![
                    Step::from(0, 0, 0f32),
                    Step::from(0, 1, -1f32),
                    Step::from(0, 2, -2f32),
                    Step::from(0, 3, -3f32),
                ],
                vec![
                    Step::from(1, 0, -1f32),
                    Step::from(1, 1, 0f32),
                    Step::from(1, 2, 0f32),
                    Step::from(1, 3, 1f32),
                ],
                vec![
                    Step::from(2, 0, -2f32),
                    Step::from(2, 1, 0f32),
                    Step::from(2, 2, 1f32),
                    Step::from(2, 3, 1f32),
                ],
            ],
            0f32,
        );

        assert_eq!(
            "AGC
CG
   |    | A  | G  | C  |
   | 0  | -1 | -2 | -3 |
C  | -1 | 0  | 0  | 1  |
G  | -2 | 0  | 1  | 1  |
",
            format!("{:?}", alignment)
        )
    }

    #[test]
    fn test_alignment_distance() {
        let alignment = Alignment::new(
            vec![vec!['A', 'C', 'C', 'G', 'T'], vec!['A', 'G', '-', 'C', 'T']],
            vec![vec![]],
            0f32,
        );

        assert_eq!(0.6, alignment.distance)
    }
}