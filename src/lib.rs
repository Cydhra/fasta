extern crate core;

use memchr::memchr;

/// A Multi FASTA file containing zero, one, or more [FastaSequence]s.
#[derive(Clone, Debug, PartialOrd, PartialEq)]
pub struct Fasta<'a> {
    pub sequences: Vec<FastaSequence<'a>>,
}

/// A FASTA sequence with a description from a FASTA file.
/// The sequence is not processed in any way, meaning accessing it will perform further parsing.
#[derive(Clone, Debug, PartialOrd, PartialEq)]
pub struct FastaSequence<'a> {
    pub description: &'a [u8],
    sequence: &'a [u8],
}

/// FASTA parsing error
#[derive(Clone, Debug, PartialOrd, PartialEq)]
pub enum ParseError {
    /// Invalid descriptor start character.
    /// The parser expects any FASTA description line to start with '>'.
    /// The invalid character is returned in the error.
    InvalidDescription { invalid: u8 },

    /// A valid descriptor was parsed, but no sequence is following
    EmptySequence,
}

impl<'a> FastaSequence<'a> {

    /// Iterate through the FASTA sequence characters without newlines.
    pub fn iter(&self) -> impl Iterator<Item = &u8> {
        self.sequence.iter().filter(|b| **b != b'\n')
    }
}

/// Parse a FASTA or Multi FASTA file.
/// Sequence descriptions are expected to start with '>'.
/// The deprecated comment character ';' is not accepted, neither for sequence descriptors nor for
/// additional comment lines.
/// Parsing is done lazily: Sequence descriptions and sequences are identified, but are not further
/// processed.
///
/// # Returns
/// A [Fasta] instance or a [ParseError] if a sequence description didn't start with a `>` character.
pub fn parse_fasta<'a, T: AsRef<[u8]>>(data: &'a T) -> Result<Fasta<'a>, ParseError> {
    let data: &'a [u8] = data.as_ref();

    let mut sequences = Vec::new();

    if data.is_empty() {
        return Ok(Fasta { sequences });
    }

    let mut cursor = 0usize;

    loop {
        if !expect(&data, b'>', &mut cursor) {
            return Err(ParseError::InvalidDescription { invalid: data[cursor] });
        }

        let header_end = memchr(b'\n', &data[cursor..]).unwrap_or(data.len() - cursor);
        let description = &data[cursor..cursor + header_end];
        cursor += header_end + 1;

        if cursor >= data.len() {
            return Err(ParseError::EmptySequence);
        }

        let sequence_end = memchr(b'>', &data[cursor..]).unwrap_or(data.len() - cursor);
        // may contain trailing white space
        let sequence = &data[cursor..cursor + sequence_end];
        cursor += sequence_end;

        sequences.push(FastaSequence {
            description,
            sequence,
        });

        if cursor >= data.len() {
            break;
        }
    }

    Ok(Fasta {
        sequences,
    })
}

/// Expect that the byte at [cursor] is equal to [expected]. If it is, advance the cursor by one.
/// Returns false, if the byte is not equal to the expected byte.
#[inline]
fn expect(data: &[u8], expected: u8, cursor: &mut usize) -> bool {
    if data[*cursor] == expected {
        *cursor += 1;
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests;