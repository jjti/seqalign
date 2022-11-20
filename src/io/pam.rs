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

// Matrix is a single alignment matrix used in scoring an alignment.
//
// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = HashMap<char, HashMap<char, i8>>;

pub enum PAM {
    PAM10,
    PAM100,
}

impl PAM {
    pub fn read(&self) -> Matrix {
        match self {
            PAM::PAM10 => read("PAM10"),
            _ => HashMap::new(),
        }
    }
}

// read returns a new PAM matrix parsed from a given file.
fn read(f: &str) -> Matrix {
    let aas = "ARNDCQEGHILKMFPSTWYVBZX*";

    let mut m: Matrix = HashMap::new();
    for aa in aas.chars() {
        let mut mm: HashMap<char, i8> = HashMap::new();
        for aa in aas.chars() {
            mm.insert(aa, 0);
        }
        m.insert(aa, mm);
    }
    m
}
