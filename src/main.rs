use std::{fs::File, io::BufReader};

use align::{align, Method};
use clap::Parser;

pub mod align;
pub mod io;
pub mod matrices;

/// Align sequences in a FASTA using one of the supported sequence alignment algorithms
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// FASTA file with sequences to align
    #[arg(index = 1, value_name = "FILE", required = true)]
    file: String,

    /// Algorithm to use
    #[arg(value_enum, short, long, default_value_t = Method::NeedlemanWunsch)]
    algo: Method,

    /// Name of the replacement matrix
    #[arg(short, long, default_value = "BLOSUM62")]
    replacement_matrix: String,

    /// Cost of a gap opening
    #[arg(long, default_value_t = -1.0)]
    gap_opening_cost: f32,

    /// Cost of extending a gap
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
    let scoring = &align::Scoring {
        matrix: match args.replacement_matrix.as_str() {
            "blosum62" => matrices::BLOSUM62::MATRIX,
            "nuc" => matrices::NUC_4_4::MATRIX,
            _ => panic!("Unknown matrix"),
        },
        gap_opening: args.gap_opening_cost,
        gap_extension: args.gap_extension_cost,
    };

    // Get the strategy
    let strategy = args.algo.strategy();

    // Align a couple sequences
    let alignment = align(vec![seq1.seq, seq2.seq], strategy, scoring);

    println!("{}", alignment)
}
