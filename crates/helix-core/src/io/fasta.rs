//! FASTA sequence and alignment parser/writer with automatic UniProt header parsing.

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use crate::sequence::{Alphabet, Sequence};

/// Parses one or more FASTA records into an [`Alignment`].
///
/// Automatically recognizes standard UniProt FASTA headers of the form:
/// `>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53`
/// and populates `Sequence::uniprot_accession` and `Sequence::organism`.
pub fn parse_fasta(name: impl Into<String>, input: &str) -> Result<Alignment> {
    let mut sequences: Vec<Sequence> = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_buf: Vec<u8> = Vec::new();

    for (line_idx, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(header) = current_header.take() {
                sequences.push(build_sequence_from_header(&header, &current_buf));
                current_buf.clear();
            }
            current_header = Some(rest.trim().to_owned());
        } else {
            if current_header.is_none() {
                return Err(HelixError::parse(
                    "FASTA",
                    Some(line_idx + 1),
                    "sequence data encountered before first '>' header line",
                ));
            }
            for b in line.bytes() {
                if !b.is_ascii_whitespace() {
                    current_buf.push(b.to_ascii_uppercase());
                }
            }
        }
    }

    if let Some(header) = current_header {
        sequences.push(build_sequence_from_header(&header, &current_buf));
    }

    if sequences.is_empty() {
        return Err(HelixError::parse("FASTA", None, "no FASTA records found"));
    }

    let alphabet = Alphabet::infer_from_slices(sequences.iter().map(|s| s.residues.as_slice()));
    Alignment::new(name, alphabet, sequences)
}

/// Serializes an [`Alignment`] into line-wrapped FASTA format (80 columns per line).
#[must_use]
pub fn write_fasta(alignment: &Alignment) -> String {
    let mut out = String::with_capacity(alignment.row_count() * (alignment.width() + 64));
    for seq in alignment.sequences() {
        out.push('>');
        if let Some(acc) = &seq.uniprot_accession {
            out.push_str("sp|");
            out.push_str(acc);
            out.push('|');
        }
        out.push_str(&seq.name);
        if !seq.description.is_empty() {
            out.push(' ');
            out.push_str(&seq.description);
        }
        if let Some(org) = &seq.organism {
            if !seq.description.contains("OS=") {
                out.push_str(" OS=");
                out.push_str(org);
            }
        }
        out.push('\n');

        for chunk in seq.residues.chunks(80) {
            for &b in chunk {
                out.push(char::from(b));
            }
            out.push('\n');
        }
    }
    out
}

/// Extracts name, description, UniProt accession, and organism from a FASTA header line.
fn build_sequence_from_header(header: &str, residues: &[u8]) -> Sequence {
    let mut parts = header.splitn(2, char::is_whitespace);
    let raw_id = parts.next().unwrap_or("Unnamed");
    let desc = parts.next().unwrap_or("").trim();

    let mut uniprot_acc = None;
    let display_name = if let Some((acc, entry_name)) = parse_uniprot_pipe_id(raw_id) {
        uniprot_acc = Some(acc.to_owned());
        entry_name.to_owned()
    } else {
        raw_id.to_owned()
    };

    let organism = extract_uniprot_field(desc, "OS=");

    let mut seq = Sequence::new(display_name, residues);
    if !desc.is_empty() {
        seq.description = desc.to_owned();
    }
    seq.uniprot_accession = uniprot_acc;
    seq.organism = organism;
    seq
}

/// Recognizes UniProt FASTA identifiers like `sp|P04637|P53_HUMAN` or `tr|A0A024R4M0|...`.
fn parse_uniprot_pipe_id(raw_id: &str) -> Option<(&str, &str)> {
    let mut segs = raw_id.split('|');
    let db = segs.next()?;
    let acc = segs.next()?;
    let entry = segs.next()?;
    if matches!(db, "sp" | "tr" | "UP") && acc.len() >= 6 && !entry.is_empty() {
        Some((acc, entry))
    } else {
        None
    }
}

/// Extracts a UniProt key=value metadata field such as `OS=Homo sapiens` from a FASTA header.
fn extract_uniprot_field(desc: &str, key: &str) -> Option<String> {
    let start = desc.find(key)? + key.len();
    let rest = &desc[start..];
    // UniProt fields end at the next two-letter uppercase key followed by '=' (e.g. ` OX=`).
    let mut end = rest.len();
    let bytes = rest.as_bytes();
    for i in 0..bytes.len().saturating_sub(3) {
        if bytes[i] == b' '
            && bytes[i + 1].is_ascii_uppercase()
            && bytes[i + 2].is_ascii_uppercase()
            && bytes[i + 3] == b'='
        {
            end = i;
            break;
        }
    }
    let val = rest[..end].trim();
    if val.is_empty() {
        None
    } else {
        Some(val.to_owned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_parses_uniprot_headers_and_roundtrips() {
        let raw = ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53\n\
                   MEEPQSDPSV\n\
                   >MOUSE_P53 Mouse ortholog\n\
                   MTAMEESQSD\n";
        let aln = parse_fasta("p53", raw).expect("parse fasta");
        assert_eq!(aln.row_count(), 2);
        assert_eq!(aln.alphabet, Alphabet::Protein);
        assert_eq!(aln.sequences()[0].name, "P53_HUMAN");
        assert_eq!(
            aln.sequences()[0].uniprot_accession.as_deref(),
            Some("P04637")
        );
        assert_eq!(
            aln.sequences()[0].organism.as_deref(),
            Some("Homo sapiens")
        );

        let exported = write_fasta(&aln);
        let reparsed = parse_fasta("roundtrip", &exported).expect("reparse");
        assert_eq!(reparsed.sequences()[0].residues, aln.sequences()[0].residues);
    }
}
