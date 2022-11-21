//! https://en.wikipedia.org/wiki/Point_accepted_mutation
//! A PAM matrix is a matrix where each column and row represents one of the twenty standard
//! amino acids. In bioinformatics, PAM matrices are sometimes used as substitution matrices
//! to score sequence alignments for proteins. Each entry in a PAM matrix indicates
//! the likelihood of the amino acid of that row being replaced with the amino acid of
//! that column through a series of one or more point accepted mutations during a
//! specified evolutionary interval, rather than these two amino acids being aligned
//! due to chance. Different PAM matrices correspond to different lengths of time in the
//! evolution of the protein sequence.
//!
//! open question: how are the PAM matricies used in BLAST for DNA? Does the DNA
//! have to correspond to a translation?

use std::collections::HashMap;
use std::fs;

/// Matrix is a single alignment matrix used in scoring an alignment.
///
/// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = HashMap<u8, HashMap<u8, i32>>;

#[derive(Debug)]
pub enum MATRIX {
    NUC,
    PAM10,
    PAM50,
    PAM100,
    PAM250,
    BLOSUM30,
    BLOSUM50,
    BLOSUM100,
}

impl MATRIX {
    pub fn read(&self) -> Matrix {
        let matrix = format!("{:?}", self);

        if matrix == "NUC" {
            read("NUC.4.4".to_string())
        } else {
            read(matrix)
        }
    }
}

/// read returns a new PAM matrix parsed from a given file.
fn read(f: String) -> Matrix {
    let mut aas: Vec<u8> = Vec::new();
    let mut m: Matrix = HashMap::new();

    let contents = fs::read_to_string(format!("matricies/{}", f.to_string())).unwrap();
    for l in contents.lines() {
        if l.starts_with("#") {
            continue;
        }
        if l.starts_with(" ") {
            for c in l.split_whitespace() {
                aas.push(c.chars().nth(0).unwrap() as u8);
            }
            continue;
        }

        let mut row = l.split_whitespace();
        let aa: u8;
        if let Some(cols) = row.next() {
            aa = cols.chars().next().unwrap() as u8;
        } else {
            continue;
        }

        let mut mm: HashMap<u8, i32> = HashMap::new();
        for (i, v) in row.map(|f| f.parse::<i32>().unwrap()).enumerate() {
            let sub = aas[i];
            mm.insert(sub, v);
        }
        m.insert(aa as u8, mm);
    }

    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read() {
        let m = MATRIX::NUC.read();
        assert_eq!(-4, *m.get(&('A' as u8)).unwrap().get(&('G' as u8)).unwrap());

        MATRIX::PAM10.read();
        MATRIX::PAM250.read();
        MATRIX::BLOSUM50.read();
    }

    #[test]
    fn test_read_pam50() {
        let m = MATRIX::PAM50.read();

        for (l, r, v) in vec![('A', 'A', 5), ('C', 'R', -6), ('A', 'D', -2)] {
            assert_eq!(v, *m.get(&(l as u8)).unwrap().get(&(r as u8)).unwrap());
        }
    }
}
