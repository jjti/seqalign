use std::{
    io::{self, BufRead},
    path::PathBuf,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("expected '@' at record start")]
    MissingAt,

    #[error("can't open {path} file: {source}")]
    FileOpen { path: PathBuf, source: io::Error },

    #[error("can't read input")]
    ReadError(#[from] io::Error),
}
pub type Result<T, E = Error> = std::result::Result<T, E>;

// A FASTA record.
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
}

impl Record {
    fn new() -> Record {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
        }
    }

    fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    fn read<R: io::Read>(&mut self, r: &mut Reader<R>) -> io::Result<()> {
        if r.line_buffer.is_empty() || !r.line_buffer.starts_with(">") {
            r.reader.read_line(&mut r.line_buffer)?;

            if r.line_buffer.is_empty() {
                return Ok(());
            }
        }

        if !r.line_buffer.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }

        let mut headers = r.line_buffer[1..].trim_end().splitn(2, " ");
        self.id = headers.next().unwrap_or_default().to_string();
        self.desc = headers.next().map(str::to_string);

        loop {
            r.line_buffer.clear();
            r.reader.read_line(&mut r.line_buffer)?;
            if r.line_buffer.is_empty() || r.line_buffer.starts_with(">") {
                return Ok(());
            }

            r.line_buffer
                .split_whitespace()
                .for_each(|s| self.seq.push_str(s));
        }
    }
}

// A FASTA Reader.
pub struct Reader<R> {
    reader: io::BufReader<R>,
    line_buffer: String,
}

impl<R: io::Read> Reader<R> {
    /// Read from a given [`io::Read`](https://doc.rust-lang.org/std/io/trait.Read.html).
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line_buffer: String::new(),
        }
    }
}

impl<R: io::Read> Iterator for Reader<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        let mut r = Record::new();

        match r.read(self) {
            Ok(()) if r.is_empty() => None,
            Ok(()) => Some(Ok(r)),
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reader_fasta() {
        let mut r = Reader::new(
            ">MCHU - Calmodulin
MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVM
DIDGDGQVNYEEFVQMMTAK*
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANG"
                .as_bytes(),
        );

        let first = r.next().unwrap().unwrap();
        assert_eq!("MCHU", first.id);
        assert_eq!("- Calmodulin", first.desc.unwrap());
        assert_eq!(
            "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMDIDGDGQVNYEEFVQMMTAK*",
            first.seq
        );

        let second = r.next().unwrap().unwrap();
        assert_eq!("SEQUENCE_2", second.id);
        assert_eq!(
            "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANG",
            second.seq
        );
    }
}
