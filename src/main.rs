use std::{fs::File, io::BufReader};

use align::Aligner;
use clap::Parser;

pub mod align;
pub mod io;
pub mod matrices;

/// CLI for the `seqalign` crate
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the file to read in
    #[arg(
        short,
        long,
        value_name = "FILE",
        help = "Path to FASTA file with sequences to align"
    )]
    file: String,

    /// Method to use for the alignment
    #[arg(short, long, default_value = "smith-waterman")]
    method: String,

    /// Name of the replacement matrix
    #[arg(short, long, default_value = "blosum62")]
    matrix: String,

    /// Cost of a gap opening
    #[arg(long, default_value_t = -1.0)]
    gap_opening_cost: f32,

    /// Name of the replacement matrix
    #[arg(long, default_value_t = -1.0)]
    gap_extension_cost: f32,
}

fn main() {
    let args = Args::parse();

    // Read seqs
    let file = File::open(&args.file).expect("Unable to open file");
    let reader = BufReader::new(file);
    let mut reader_fasta = io::fasta::Reader::new(reader);
    let seq1 = reader_fasta
        .next()
        .expect("Missing first seq to align")
        .unwrap();
    let seq2 = reader_fasta
        .next()
        .expect("Missing second seq to align")
        .unwrap();

    // Create scoring
    let scoring = align::Scoring {
        matrix: match args.matrix.as_str() {
            "blosum62" => matrices::BLOSUM62::MATRIX,
            "nuc" => matrices::NUC_4_4::MATRIX,
            _ => panic!("Unknown matrix"),
        },
        gap_opening: args.gap_opening_cost,
        gap_extension: args.gap_extension_cost,
    };

    // Align
    let alignment = match args.method.as_str() {
        "needleman-wunsch" => {
            let aligner = align::NeedlemanWunsch::new();
            aligner.align(vec![seq1.seq, seq2.seq], scoring)
        }
        "smith-waterman" => {
            let aligner = align::SmithWaterman::new();
            aligner.align(vec![seq1.seq, seq2.seq], scoring)
        }
        _ => panic!("Unknown method"),
    };

    println!("{}", alignment)
}
