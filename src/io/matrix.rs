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

// Matrix is a single alignment matrix used in scoring an alignment.
//
// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = HashMap<char, HashMap<char, i8>>;

#[derive(Debug)]
pub enum MATRIX {
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
        read(format!("{:?}", self))
    }
}

// read returns a new PAM matrix parsed from a given file.
fn read(f: String) -> Matrix {
    let aas = "ARNDCQEGHILKMFPSTWYVBZX*".as_bytes();
    let mut m: Matrix = HashMap::new();

    let contents = fs::read_to_string(format!("matricies/{}", f.to_string())).unwrap();
    for l in contents.lines() {
        if l.starts_with("#") || l.starts_with(" ") {
            continue;
        }

        let mut row = l.split_whitespace();
        let aa = row.next().unwrap().chars().next().unwrap();
        let mut mm: HashMap<char, i8> = HashMap::new();
        for (i, v) in row.map(|f| f.parse::<i8>().unwrap()).enumerate() {
            let sub = aas[i];
            mm.insert(sub as char, v);
        }
        m.insert(aa, mm);
    }

    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read() {
        MATRIX::PAM10.read();
        MATRIX::PAM250.read();
        MATRIX::BLOSUM50.read();
    }

    #[test]
    fn test_read_pam50() {
        let m = MATRIX::PAM50.read();

        for (l, r, v) in vec![('A', 'A', 5), ('C', 'R', -6), ('A', 'D', -2)] {
            assert_eq!(v, *m.get(&l).unwrap().get(&r).unwrap());
        }
    }
}