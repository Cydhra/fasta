use memchr::memchr;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// A Multi FASTA file containing zero, one, or more [FastaSequence]s.
#[derive(Clone, Debug)]
pub struct Fasta<'a> {
    pub sequences: Vec<FastaSequence<'a>>,
}

/// A FASTA sequence with a description from a FASTA file.
/// The sequence is not processed in any way, meaning accessing it will perform further parsing.
#[derive(Clone, Debug)]
pub struct FastaSequence<'a> {
    pub description: &'a [u8],
    sequence: &'a [u8],
}

/// FASTA parsing error
#[derive(Clone, Debug)]
pub enum ParseError {
    /// Invalid descriptor start character.
    /// The parser expects any FASTA description line to start with '>'.
    /// The invalid character is returned in the error.
    InvalidDescription { invalid: u8 },

    /// A valid descriptor was parsed, but no sequence is following
    EmptySequence,
}

impl Display for ParseError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Error for ParseError {}

impl<'a> FastaSequence<'a> {
    /// Returns an iterator over the FASTA sequence characters, excluding newlines.
    /// Note that the parser expects unix-style line breaks, thus CR-characters are preserved.
    pub fn iter(&self) -> impl Iterator<Item = &u8> {
        self.sequence.iter().filter(|&x| *x != b'\n')
    }

    /// Copy the sequence into a consecutive memory region.
    /// This method allocates a buffer and copies the sequence into it without newline symbols.
    /// Note that any other symbol (including whitespace) gets preserved.
    /// The returned allocation capacity may be larger than the actual sequence.
    /// It is guaranteed, however, that only one allocation is performed.
    pub fn copy_sequential(&self) -> Box<[u8]> {
        let mut buffer = vec![0u8; self.sequence.len()];
        let mut target = 0;
        let mut pos = 0;
        loop {
            let pivot = memchr(b'\n', &self.sequence[pos..]);
            if let Some(q) = pivot {
                buffer[target..target + q].copy_from_slice(&self.sequence[pos..pos + q]);
                pos += q + 1;
                target += q;

                if pos >= self.sequence.len() {
                    break;
                }
            }
        }
        buffer.truncate(target);
        buffer.into_boxed_slice()
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
            return Err(ParseError::InvalidDescription {
                invalid: data[cursor],
            });
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

    Ok(Fasta { sequences })
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
