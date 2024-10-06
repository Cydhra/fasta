use crate::{parse_fasta, Fasta};

#[test]
fn empty_fasta() {
    let empty = "";
    assert_eq!(parse_fasta(&empty), Ok(Fasta { sequences: vec![] }));
}

#[test]
fn one_sequence() {
    let seq = r#"
>P32234 1-368
MSTILEKISAIESEMARTQKNKATSAHLGLLKAKLAKLRRELISPKGGGGGTGEAGFEVAKTGDARVGFVGFPSVGKSTL
    "#.trim();

    let fasta = parse_fasta(&seq).expect("Failed to parse FASTA");
    assert_eq!(fasta.sequences.len(), 1);
    assert_eq!(fasta.sequences[0].description, b"P32234 1-368");
    assert_eq!(String::from_utf8(fasta.sequences[0].sequence.to_vec()).unwrap().trim(), "MSTILEKISAIESEMARTQKNKATSAHLGLLKAKLAKLRRELISPKGGGGGTGEAGFEVAKTGDARVGFVGFPSVGKSTL")
}

#[test]
fn multi_sequence() {
    let seq = r"
>P32234 1-368
MSTILEKISAIESEMARTQKNKATSAHLGLLKAKLAKLRRELISPKGGGGGTGEAGFEVAKTGDARVGFVGFPSVGKSTL

>O77448 1-1117
MQKINNINNNKQMLTRKEDLLTVLKQISALKYVSNLYEFLLATEKIVQTSELDTQFQEFLTTTIIASEQNLVENYKQKYN
    ".trim();

    let fasta = parse_fasta(&seq).expect("Failed to parse FASTA");
    assert_eq!(fasta.sequences.len(), 2);

    assert_eq!(fasta.sequences[0].description, b"P32234 1-368");
    assert_eq!(String::from_utf8(fasta.sequences[0].sequence.to_vec()).unwrap().trim(), "MSTILEKISAIESEMARTQKNKATSAHLGLLKAKLAKLRRELISPKGGGGGTGEAGFEVAKTGDARVGFVGFPSVGKSTL");

    assert_eq!(fasta.sequences[1].description, b"O77448 1-1117");
    assert_eq!(String::from_utf8(fasta.sequences[1].sequence.to_vec()).unwrap().trim(), "MQKINNINNNKQMLTRKEDLLTVLKQISALKYVSNLYEFLLATEKIVQTSELDTQFQEFLTTTIIASEQNLVENYKQKYN");
}