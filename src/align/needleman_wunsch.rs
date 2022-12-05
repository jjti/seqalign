//! https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
//! 1970. Global optimal alignment, meaning total alignment between the two sequences.
//!
//! The original purpose of the algorithm described by Needleman and Wunsch was
//! to find similarities in the amino acid sequences of two proteins.
//!
//! The Needleman–Wunsch algorithm is still widely used for optimal global alignment,
//! particularly when the quality of the global alignment is of the utmost importance.
//! The algorithm assigns a score to every possible alignment, and the purpose of the
//! algorithm is to find all possible alignments having the highest score.

use ordered_float::OrderedFloat;

use super::{PWAlign, PWAlignment, Scoring, Step};

pub struct Aligner<'a> {
    grid: Vec<Vec<Step>>,
    a: &'a str,
    b: &'a str,
    scoring: Scoring,
}

impl<'a> PWAlign for Aligner<'a> {
    fn align(&mut self) -> PWAlignment {
        self.init();
        self.fill();
        let (a, b, score) = self.backtrace();

        PWAlignment {
            grid: self.grid.to_owned(),
            a,
            b,
            a_orig: self.a.to_string(),
            b_orig: self.b.to_string(),
            score,
        }
    }
}

impl<'a> Aligner<'a> {
    pub fn new(a: &'a str, b: &'a str, scoring: Scoring) -> Aligner<'a> {
        Aligner {
            grid: Vec::new(),
            a,
            b,
            scoring,
        }
    }

    fn init(&mut self) {
        for i in 0..self.b.len() + 1 {
            self.grid.push(vec![Step::default(); self.a.len() + 1]);

            self.grid[i][0] = Step {
                i,
                j: 0,
                val: OrderedFloat(-(i as f32)),
                next: match i {
                    0 => None,
                    _ => Some((i - 1, 0)),
                },
            };
        }

        for j in 0..self.a.len() + 1 {
            self.grid[0][j] = Step {
                i: 0,
                j,
                val: OrderedFloat(-(j as f32)),
                next: match j {
                    0 => None,
                    _ => Some((0, j - 1)),
                },
            };
        }
    }

    fn fill(&mut self) {
        let a = self.a.as_bytes();
        let b = self.b.as_bytes();

        for i in 1..b.len() + 1 {
            for j in 1..a.len() + 1 {
                let mut opts: Vec<Step> = vec![
                    // indel
                    Step {
                        i,
                        j,
                        val: self.grid[i - 1][j].val + self.scoring.gap_opening,
                        next: Some((i - 1, j)),
                    },
                    Step {
                        i,
                        j,
                        val: self.grid[i][j - 1].val + self.scoring.gap_opening,
                        next: Some((i, j - 1)),
                    },
                ];

                // match or mismatch
                let match_val =
                    self.scoring.replacement[a[j - 1] as usize][b[i - 1] as usize] as f32;
                opts.push(Step {
                    val: self.grid[i - 1][j - 1].val + match_val,
                    i,
                    j,
                    next: Some((i - 1, j - 1)),
                });

                self.grid[i][j] = opts.iter().max().unwrap().clone();
            }
        }
    }

    fn backtrace(&self) -> (String, String, f32) {
        let a = self.a.as_bytes();
        let b = self.b.as_bytes();
        let mut a_row: Vec<u8> = Vec::new();
        let mut b_row: Vec<u8> = Vec::new();

        let mut step = &self.grid[self.grid.len() - 1][self.a.len()];
        let score = step.val.0;
        while let Some((next_i, next_j)) = step.next {
            let i_delta = step.i - next_i;
            let j_delta = step.j - next_j;

            if i_delta == 1 && j_delta == 1 {
                // match/mismatch
                a_row.push(a[step.j - 1]);
                b_row.push(b[step.i - 1]);
            } else if i_delta > 0 {
                // gap in seq a
                let mut i = step.i;
                while i > next_i {
                    a_row.push(b'-');
                    b_row.push(b[i - 1]);
                    i -= 1;
                }
            } else if j_delta > 0 {
                // gap in seq b
                let mut j = step.j;
                while j > next_j {
                    a_row.push(a[j - 1]);
                    b_row.push(b'-');
                    j -= 1;
                }
            } else {
                panic!("unexpected step");
            }

            // move to the next step in the alignment
            step = &self.grid[next_i][next_j];
        }

        let top: String = a_row.iter().rev().map(|c| *c as char).collect();
        let bottom: String = b_row.iter().rev().map(|c| *c as char).collect();
        (top, bottom, score)
    }
}

#[cfg(test)]
mod tests {
    use crate::{io::fasta::Reader, matrices::BLOSUM62, matrices::MATCH};

    use super::*;

    #[test]
    fn test_aligner_align() {
        let mut a = Aligner::new(
            "GCATGCG",
            "GATTACA",
            Scoring {
                replacement: MATCH::MATRIX,
                gap_opening: -1f32,
                gap_extension: -1f32,
            },
        );
        let alignment = a.align();

        println!("{:?}", alignment);

        assert_eq!("GCA-TGCG", alignment.a);
        assert_eq!("G-ATTACA", alignment.b);
    }

    /// read in tests/data/reper.pep to fasta and align the first two records
    #[test]
    fn test_aligner_align_fasta() {
        let mut r =
            Reader::new(std::fs::File::open("./tests/data/reper.pep").expect("can't open file"));

        let first = r.next().unwrap().unwrap();
        let second = r.next().unwrap().unwrap();

        let a = Aligner::new(
            first.seq.as_str(),
            second.seq.as_str(),
            Scoring {
                replacement: BLOSUM62::MATRIX,
                gap_opening: -1f32,
                gap_extension: -0.5f32,
            },
        )
        .align();

        // EMBOSS Needle alignment results
        // assert_eq!("METKNLTIGERIRY-RRKNLKH---TQRSLAKALKIS--HVSVSQWERGDSEPTGKNLFALSKVLQCSP-TWILFGDEDKQPTPPVEKPVALSPKELEL-LELFNALPESEQDTQLAEMRARVKNFNKLFEELLKAR--QR---TN-KR---", a.a);
        // assert_eq!("LDGKKL--GALIK-DKRKE-KHLKQTE--MAKALGMSRTYLS-------DIE-NGR--Y-L-------PST--------K--T--LSR-IAI----L-INLDL-NVL---KM-T---EI--QV-----V-EE---G-GYDRAAGT-CRRQAL", a.b);
    }
}
