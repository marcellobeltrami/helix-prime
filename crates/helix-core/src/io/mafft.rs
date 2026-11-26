//! Pre-aligned MAFFT file importer supporting all three standard MAFFT output formats:
//! 1. **Aligned FASTA** (default `mafft` stdout, including lowercase sequences and `_R_`
//!    reverse-complemented headers produced by `mafft --adjustdirection`),
//! 2. **MAFFT Clustal output** (`mafft --clustalout`, e.g. `CLUSTAL format alignment by MAFFT L-INS-i`),
//! 3. **MAFFT PHYLIP output** (`mafft --phylipout`, both interleaved and sequential blocks).
//!
//! Pre-aligned columns and gap characters (`-`) are preserved verbatim without re-aligning.

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use crate::io::{clustal::parse_clustal, fasta::parse_fasta};
use crate::sequence::{Alphabet, Sequence};
use std::path::Path;

/// Detected serialization format of a pre-aligned MAFFT file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MafftFlavor {
    /// Standard aligned multi-FASTA output (`mafft in.fa > out.mafft`).
    AlignedFasta,
    /// Clustal-formatted MAFFT output (`mafft --clustalout`).
    ClustalOut,
    /// PHYLIP-formatted MAFFT output (`mafft --phylipout`).
    PhylipOut,
}

impl MafftFlavor {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::AlignedFasta => "MAFFT Aligned FASTA",
            Self::ClustalOut => "MAFFT Clustal (--clustalout)",
            Self::PhylipOut => "MAFFT PHYLIP (--phylipout)",
        }
    }
}

/// Result of importing a pre-aligned MAFFT file, preserving exact columns and MAFFT metadata.
#[derive(Debug, Clone)]
pub struct MafftImportReport {
    pub alignment: Alignment,
    pub flavor: MafftFlavor,
    /// Optional header banner captured from `--clustalout` (e.g. `CLUSTAL format alignment by MAFFT L-INS-i (v7.520)`).
    pub algorithm_banner: Option<String>,
    /// Number of sequences prefixed with `_R_` by `mafft --adjustdirection`.
    pub adjusted_direction_count: usize,
}

/// Loads a pre-aligned MAFFT file from `path`, auto-detecting FASTA, `--clustalout`, or `--phylipout`.
pub fn load_mafft_alignment_from_path(path: &Path) -> Result<MafftImportReport> {
    let content = std::fs::read_to_string(path)?;
    let doc_name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("MAFFT Alignment")
        .to_owned();
    parse_mafft(doc_name, &content)
}

/// Parses an already-aligned MAFFT text payload (aligned FASTA, `--clustalout`, or `--phylipout`),
/// preserving all gap columns verbatim without running any re-alignment.
pub fn parse_mafft(name: impl Into<String>, input: &str) -> Result<MafftImportReport> {
    let doc_name = name.into();
    let first_line = input
        .lines()
        .map(str::trim)
        .find(|l| !l.is_empty())
        .ok_or_else(|| HelixError::parse("MAFFT", None, "file is empty"))?;

    let (flavor, banner, mut alignment) = if first_line.starts_with('>') {
        let aln = parse_fasta(doc_name, input)?;
        (MafftFlavor::AlignedFasta, None, aln)
    } else if first_line.to_ascii_uppercase().starts_with("CLUSTAL")
        || first_line.to_ascii_uppercase().starts_with("MUSCLE")
    {
        let aln = parse_clustal(doc_name, input)?;
        (MafftFlavor::ClustalOut, Some(first_line.to_owned()), aln)
    } else if parse_phylip_header(first_line).is_some() {
        let aln = parse_phylip(doc_name, input)?;
        (MafftFlavor::PhylipOut, None, aln)
    } else {
        return Err(HelixError::parse(
            "MAFFT",
            Some(1),
            "unrecognized MAFFT alignment format (expected aligned FASTA '>', MAFFT '--clustalout', or MAFFT '--phylipout')",
        ));
    };

    let adjusted_direction_count = alignment
        .sequences()
        .iter()
        .filter(|s| s.name.starts_with("_R_"))
        .count();

    let mut desc = match &banner {
        Some(b) => format!("Pre-aligned {b}"),
        None => format!("Pre-aligned {} ({} cols)", flavor.label(), alignment.width()),
    };
    if adjusted_direction_count > 0 {
        desc.push_str(&format!(
            " | {adjusted_direction_count} seq(s) reverse-adjusted (_R_)"
        ));
    }
    alignment.description = desc;

    Ok(MafftImportReport {
        alignment,
        flavor,
        algorithm_banner: banner,
        adjusted_direction_count,
    })
}

/// Returns `Some((num_taxa, alignment_len))` if `line` is a valid PHYLIP header (`<N> <L>`).
#[must_use]
pub fn parse_phylip_header(line: &str) -> Option<(usize, usize)> {
    let mut parts = line.split_whitespace();
    let n: usize = parts.next()?.parse().ok()?;
    let len: usize = parts.next()?.parse().ok()?;
    // Reject extra non-flag trailing tokens so arbitrary text doesn't false-match.
    if let Some(extra) = parts.next() {
        if !matches!(extra, "I" | "S" | "i" | "s") {
            return None;
        }
    }
    if n > 0 && len > 0 {
        Some((n, len))
    } else {
        None
    }
}

/// Parses a PHYLIP multiple sequence alignment (supporting both interleaved `--phylipout`
/// produced by MAFFT and sequential/relaxed PHYLIP files).
pub fn parse_phylip(name: impl Into<String>, input: &str) -> Result<Alignment> {
    let mut lines = input.lines().enumerate().peekable();

    let mut header: Option<(usize, usize)> = None;
    for (_, raw) in lines.by_ref() {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            continue;
        }
        header = parse_phylip_header(trimmed);
        break;
    }

    let Some((num_taxa, expected_len)) = header else {
        return Err(HelixError::parse(
            "PHYLIP",
            Some(1),
            "expected PHYLIP header '<num_sequences> <alignment_length>' on first line",
        ));
    };

    let mut names: Vec<String> = Vec::with_capacity(num_taxa);
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(num_taxa);

    // Collect remaining non-empty lines with 1-based line numbers.
    let body: Vec<(usize, &str)> = lines
        .filter_map(|(idx, l)| {
            let t = l.trim();
            if t.is_empty() {
                None
            } else {
                Some((idx + 1, l))
            }
        })
        .collect();

    if body.len() < num_taxa {
        return Err(HelixError::parse(
            "PHYLIP",
            None,
            format!("header declares {num_taxa} sequences, but found fewer data lines"),
        ));
    }

    // Parse first block (one line per taxon with taxon name + initial sequence chunk).
    // If sequential PHYLIP is used, each taxon reaches `expected_len` before the next taxon starts.
    let mut cursor = 0usize;
    let mut sequential_mode = false;

    for taxon_idx in 0..num_taxa {
        if cursor >= body.len() {
            return Err(HelixError::parse(
                "PHYLIP",
                None,
                format!("unexpected end of PHYLIP file before sequence {}", taxon_idx + 1),
            ));
        }
        let (line_no, raw_line) = body[cursor];
        cursor += 1;

        let (taxon_name, first_chunk) = split_phylip_name_and_seq(raw_line).ok_or_else(|| {
            HelixError::parse(
                "PHYLIP",
                Some(line_no),
                "could not split taxon name and sequence data",
            )
        })?;

        let mut buf = Vec::with_capacity(expected_len);
        append_clean_residues(&mut buf, first_chunk);

        // Detect sequential PHYLIP on the first taxon: if continuation lines belong to taxon 0
        // (leading whitespace or reaching expected_len before taxon 1), keep reading until full.
        if taxon_idx == 0 && buf.len() < expected_len && cursor < body.len() {
            let next_raw = body[cursor].1;
            let next_starts_with_ws = next_raw.starts_with(' ') || next_raw.starts_with('\t');
            // In interleaved MAFFT --phylipout, the first `num_taxa` lines all start with taxon names
            // (non-whitespace at col 0), whereas sequential PHYLIP continuation lines often start with spaces
            // OR have total lines == num_taxa * ceil(expected_len / block_width) without blank separators.
            if next_starts_with_ws {
                sequential_mode = true;
            }
        }

        if sequential_mode {
            while buf.len() < expected_len && cursor < body.len() {
                append_clean_residues(&mut buf, body[cursor].1);
                cursor += 1;
            }
        }

        names.push(taxon_name);
        seqs.push(buf);
    }

    // Interleaved continuation blocks: cycle across taxa `0..num_taxa` until all remaining lines are read.
    if !sequential_mode {
        let mut idx = 0usize;
        while cursor < body.len() {
            let (_, raw_line) = body[cursor];
            cursor += 1;
            append_clean_residues(&mut seqs[idx % num_taxa], raw_line);
            idx += 1;
        }
    }

    let sequences: Vec<Sequence> = names
        .into_iter()
        .zip(seqs)
        .map(|(seq_name, residues)| Sequence::new(seq_name, &residues))
        .collect();

    if sequences.is_empty() || sequences.iter().all(Sequence::is_empty) {
        return Err(HelixError::parse(
            "PHYLIP",
            None,
            "no aligned residues found in PHYLIP file",
        ));
    }

    let alphabet = Alphabet::infer_from_slices(sequences.iter().map(|s| s.residues.as_slice()));
    Alignment::new(name, alphabet, sequences)
}

/// Splits a PHYLIP taxon line into `(taxon_name, sequence_chunk)`.
/// Supports both strict 10-character PHYLIP names (`mafft --phylipout`) and relaxed whitespace-separated names.
fn split_phylip_name_and_seq(raw_line: &str) -> Option<(String, &str)> {
    let trimmed_start = raw_line.trim_start();
    if trimmed_start.is_empty() {
        return None;
    }

    // Standard MAFFT `--phylipout` emits a 10-character fixed name column followed by sequence
    // blocks separated by spaces (e.g. `seq_01    atgcgtacgt acgt--acgt`).
    if let Some(ws_idx) = trimmed_start.find(char::is_whitespace) {
        let name = trimmed_start[..ws_idx].trim();
        let rest = trimmed_start[ws_idx..].trim();
        if !name.is_empty() && !rest.is_empty() {
            return Some((name.to_owned(), rest));
        }
    }

    // Strict 10-column fallback when name occupies all 10 chars with no separating space.
    if raw_line.len() > 10 {
        let name = raw_line[..10].trim();
        let rest = raw_line[10..].trim();
        if !name.is_empty() && !rest.is_empty() {
            return Some((name.to_owned(), rest));
        }
    }

    None
}

fn append_clean_residues(buf: &mut Vec<u8>, chunk: &str) {
    for b in chunk.bytes() {
        if b.is_ascii_alphabetic() || matches!(b, b'-' | b'.' | b'*' | b'?') {
            let norm = if b == b'.' { b'-' } else { b.to_ascii_uppercase() };
            buf.push(norm);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_mafft_aligned_fasta_with_adjustdirection_and_lowercase() {
        let input = ">seq_alpha\natgcgtacg--tacg\n>_R_seq_beta\natgcg-acgtttacg\n>seq_gamma\natgcgtacg--tacc\n";
        let rep = parse_mafft("MyMafft", input).expect("parse mafft fasta");
        assert_eq!(rep.flavor, MafftFlavor::AlignedFasta);
        assert_eq!(rep.adjusted_direction_count, 1);
        assert_eq!(rep.alignment.row_count(), 3);
        assert_eq!(rep.alignment.width(), 15);
        assert_eq!(rep.alignment.sequences()[0].residues, b"ATGCGTACG--TACG");
        assert_eq!(rep.alignment.sequences()[1].residues, b"ATGCG-ACGTTTACG");
    }

    #[test]
    fn parses_mafft_clustalout_format() {
        let input = "\
CLUSTAL format alignment by MAFFT L-INS-i (v7.520)

seq1            atgcgtacgt--acgtacgt
seq2            atgcg-acgtttacgtacgt
                ***** ****  ********

seq1            ggccaatt
seq2            ggccaatc
                *******
";
        let rep = parse_mafft("MafftClustal", input).expect("parse mafft clustalout");
        assert_eq!(rep.flavor, MafftFlavor::ClustalOut);
        assert!(rep
            .algorithm_banner
            .as_deref()
            .is_some_and(|b| b.contains("MAFFT L-INS-i")));
        assert_eq!(rep.alignment.row_count(), 2);
        assert_eq!(rep.alignment.width(), 28);
        assert_eq!(rep.alignment.sequences()[0].residues, b"ATGCGTACGT--ACGTACGTGGCCAATT");
    }

    #[test]
    fn parses_mafft_phylipout_interleaved_format() {
        let input = "\
 3 24
taxon_a   atgcgtacgt --acgtacgt
taxon_b   atgcg-acgt ttacgtacgt
taxon_c   atgcgtacgt --acgtacca

          ggcc
          ggct
          ggcc
";
        let rep = parse_mafft("MafftPhylip", input).expect("parse mafft phylipout");
        assert_eq!(rep.flavor, MafftFlavor::PhylipOut);
        assert_eq!(rep.alignment.row_count(), 3);
        assert_eq!(rep.alignment.width(), 24);
        assert_eq!(rep.alignment.sequences()[0].name, "taxon_a");
        assert_eq!(rep.alignment.sequences()[0].residues, b"ATGCGTACGT--ACGTACGTGGCC");
        assert_eq!(rep.alignment.sequences()[1].residues, b"ATGCG-ACGTTTACGTACGTGGCT");
    }
}
