# Fire Fasta
Ultra-fast, lightweight, zero-copy, lazy Multi-FASTA parser.

The parser is intended for high performance applications where the input is expected to be well-formed.
Therefore, it sacrifices input validation and deprecated features for parsing performance.

### Sequence Characters
The parser makes no assumptions about the sequence alphabet:
The parser is explicitly intended for custom sequences with characters that do not conform to NCBI specifications.
The only illegal characters in sequences are unix-style newlines (LF), which are ignored, and the greater-than sign,
which starts a new sequence descriptor in Multi-FASTA files.
Note, that the parser does not validate whether a sequence description starts at the beginning of a new line.

The parser expects input data that is compatible with ASCII.
Multibyte UTF-8 codepoints are processed as separate ASCII characters.

Windows-style newlines (CRLF) are not supported.
Instead, the parser will treat the LF as a unix-style newline and preserve the CR as a valid sequence character.
Old FASTA comments starting with `;` are also not supported, they are treated as part of the sequence.

### Usage and Lazy Parsing
Calling the parser will do one pass over the entire input, separating individual fasta sequences from each other.
No further processing is done and no data is copied.
```rust
let seq = ">example\nMSTIL\nAATIL\n\n";
let fasta = parse_fasta_str(&seq).expect("Failed to parse FASTA");
// or parse_fasta(&data) for &[u8] slices

assert_eq!(fasta.sequences.len(), 1);
```

Iterating over a sequence will remove newlines from the iterator on the fly:
```rust
assert_eq!(
    String::from_utf8(fasta.sequences[0].iter().copied().collect::<Vec<_>>()).unwrap(),
    "MSTILAATIL"
);
```

If you want to iterate over a sequence multiple times, it may be faster to first copy the full sequence into its own buffer:
```rust
let copied: Box<[u8]> = fasta.sequences[0].copy_sequential();
assert_eq!(copied.as_ref(), b"MSTILAATIL");
```

Parsing and copying use the [memchr](https://crates.io/crates/memchr) crate, 
and thus operations use SIMD instructions when available.

### Validation and Convenience
If you require input validation or features like Windows-style newlines, have a look at [seq_io](https://crates.io/crates/seq_io).

### License
Licensed under either of

* Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license
  ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Contribution
Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.