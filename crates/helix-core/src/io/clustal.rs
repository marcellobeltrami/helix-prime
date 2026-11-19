//! ClustalW / Clustal Omega (`.aln`) alignment parser and writer.

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use crate::sequence::{Alphabet, Sequence};

/// Parses a Clustal-format multiple sequence alignment (`.aln`).
pub fn parse_clustal(name: impl Into<String>, input: &str) -> Result<Alignment> {
    let mut order: Vec<String> = Vec::new();
    let mut seq_data: Vec<Vec<u8>> = Vec::new();
    let mut saw_header = false;

    for (line_idx, raw_line) in input.lines().enumerate() {
        let trimmed = raw_line.trim();
        if !saw_header {
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.to_ascii_uppercase().starts_with("CLUSTAL")
                || trimmed.to_ascii_uppercase().starts_with("MUSCLE")
            {
                saw_header = true;
                continue;
            }
            return Err(HelixError::parse(
                "Clustal",
                Some(line_idx + 1),
                "expected CLUSTAL or MUSCLE header on first non-empty line",
            ));
        }

        // Blank lines and consensus symbol lines (`*`, `:`, `.`) begin with whitespace or are empty.
        if trimmed.is_empty()
            || raw_line.starts_with(' ')
            || raw_line.starts_with('\t')
            || trimmed
                .bytes()
                .all(|b| matches!(b, b'*' | b':' | b'.' | b' '))
        {
            continue;
        }

        let mut tokens = trimmed.split_whitespace();
        let seq_name = match tokens.next() {
            Some(s) => s,
            None => continue,
        };
        let chunk = match tokens.next() {
            Some(c) => c,
            None => {
                return Err(HelixError::parse(
                    "Clustal",
                    Some(line_idx + 1),
                    format!("missing aligned sequence chunk for '{seq_name}'"),
                ));
            }
        };

        if let Some(pos) = order.iter().position(|n| n == seq_name) {
            seq_data[pos].extend(chunk.bytes().map(|b| b.to_ascii_uppercase()));
        } else {
            order.push(seq_name.to_owned());
            seq_data.push(chunk.bytes().map(|b| b.to_ascii_uppercase()).collect());
        }
    }

    if order.is_empty() {
        return Err(HelixError::parse(
            "Clustal",
            None,
            "no aligned sequence blocks found",
        ));
    }

    let sequences: Vec<Sequence> = order
        .into_iter()
        .zip(seq_data)
        .map(|(seq_name, residues)| Sequence::new(seq_name, &residues))
        .collect();

    let alphabet = Alphabet::infer_from_slices(sequences.iter().map(|s| s.residues.as_slice()));
    Alignment::new(name, alphabet, sequences)
}

/// Serializes an [`Alignment`] into standard ClustalW (`.aln`) format with 60 columns per block.
#[must_use]
pub fn write_clustal(alignment: &Alignment) -> String {
    let mut out = String::from("CLUSTAL W (HelixPrime) multiple sequence alignment\n\n");
    let max_name = alignment
        .sequences()
        .iter()
        .map(|s| s.name.len())
        .max()
        .unwrap_or(10)
        .clamp(10, 30);

    let block_size = 60usize;
    let width = alignment.width();
    let mut pos = 0usize;

    while pos < width {
        let end = (pos + block_size).min(width);
        for seq in alignment.sequences() {
            let short_name: String = seq.name.chars().take(max_name).collect();
            out.push_str(&short_name);
            for _ in short_name.len()..(max_name + 4) {
                out.push(' ');
            }
            for &b in &seq.residues[pos..end] {
                out.push(char::from(b));
            }
            out.push('\n');
        }
        // Conservation line.
        for _ in 0..(max_name + 4) {
            out.push(' ');
        }
        for col in pos..end {
            let id = alignment.columns()[col].identity_pct;
            if id >= 99.9 {
                out.push('*');
            } else if id >= 75.0 {
                out.push(':');
            } else if id >= 50.0 {
                out.push('.');
            } else {
                out.push(' ');
            }
        }
        out.push_str("\n\n");
        pos = end;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clustal_roundtrip_preserves_sequences() {
        let sample = "CLUSTAL W (1.81) multiple sequence alignment\n\n\
                      HUMAN_HBB      MVHLTPEEKSAVTALWGKVNVDEVGGEALGR 31\n\
                      MOUSE_HBB      MVHLTDAEKAAVSCLWGKVNSDEVGGEALGR 31\n\
                                     *****  **:**: ******.**********\n\n\
                      HUMAN_HBB      LLVVYPWTQRFFESFGDLSTPDAVMGNPKVK 62\n\
                      MOUSE_HBB      LLVVYPWTQRYFDSFGDLSSASAIMGNAKVK 62\n\
                                     **********:*:******...*:***.***\n";
        let aln = parse_clustal("HBB", sample).expect("parse clustal");
        assert_eq!(aln.row_count(), 2);
        assert_eq!(aln.width(), 62);
        let serialized = write_clustal(&aln);
        let reparsed = parse_clustal("HBB2", &serialized).expect("reparse clustal");
        assert_eq!(reparsed.sequences()[0].residues, aln.sequences()[0].residues);
    }
}
