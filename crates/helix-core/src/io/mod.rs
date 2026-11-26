//! Unified bioinformatics file format detection, loading, and export.

pub mod clustal;
pub mod fasta;
pub mod genbank;
pub mod gff3;
pub mod mafft;

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use std::path::Path;

/// Supported alignment and sequence input formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentFormat {
    Fasta,
    Clustal,
    GenBank,
    Mafft,
    Phylip,
}

/// Auto-detects format from file extension or leading content and loads an [`Alignment`].
pub fn load_alignment_from_path(path: &Path) -> Result<Alignment> {
    let content = std::fs::read_to_string(path)?;
    let doc_name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("Imported Alignment")
        .to_owned();

    let format = detect_format(path, &content)?;
    match format {
        AlignmentFormat::Fasta => fasta::parse_fasta(doc_name, &content),
        AlignmentFormat::Clustal => clustal::parse_clustal(doc_name, &content),
        AlignmentFormat::GenBank => genbank::parse_genbank(doc_name, &content),
        AlignmentFormat::Mafft => mafft::parse_mafft(doc_name, &content).map(|r| r.alignment),
        AlignmentFormat::Phylip => mafft::parse_phylip(doc_name, &content),
    }
}

/// Detects sequence/alignment format using content-aware sniffing paired with file extensions
/// (so e.g. a MAFFT aligned FASTA saved as `.aln` or `.mafft` still parses cleanly).
pub fn detect_format(path: &Path, content: &str) -> Result<AlignmentFormat> {
    let first_line = content.lines().map(str::trim).find(|l| !l.is_empty());

    // First check unambiguous content headers (critical when MAFFT FASTA output is saved as `.aln`).
    if let Some(line) = first_line {
        let upper = line.to_ascii_uppercase();
        if upper.starts_with("CLUSTAL") && upper.contains("MAFFT") {
            return Ok(AlignmentFormat::Mafft);
        }
        if upper.starts_with("LOCUS") {
            return Ok(AlignmentFormat::GenBank);
        }
    }

    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "mafft" | "afa" | "mfa" | "msa" => return Ok(AlignmentFormat::Mafft),
            "phy" | "phylip" => return Ok(AlignmentFormat::Phylip),
            "fa" | "fasta" | "fna" | "faa" | "fas" => return Ok(AlignmentFormat::Fasta),
            "aln" | "clustal" | "clw" => {
                // Researchers often redirect default MAFFT FASTA output to `*.aln`; sniff leading `>`.
                if first_line.is_some_and(|l| l.starts_with('>')) {
                    return Ok(AlignmentFormat::Mafft);
                }
                if first_line.is_some_and(|l| mafft::parse_phylip_header(l).is_some()) {
                    return Ok(AlignmentFormat::Phylip);
                }
                return Ok(AlignmentFormat::Clustal);
            }
            "gb" | "gbk" | "genbank" => return Ok(AlignmentFormat::GenBank),
            _ => {}
        }
    }

    if let Some(trimmed) = first_line {
        if trimmed.starts_with('>') {
            return Ok(AlignmentFormat::Fasta);
        }
        let upper = trimmed.to_ascii_uppercase();
        if upper.starts_with("CLUSTAL") || upper.starts_with("MUSCLE") {
            return Ok(AlignmentFormat::Clustal);
        }
        if upper.starts_with("LOCUS") {
            return Ok(AlignmentFormat::GenBank);
        }
        if mafft::parse_phylip_header(trimmed).is_some() {
            return Ok(AlignmentFormat::Phylip);
        }
    }

    Err(HelixError::parse(
        "AutoDetect",
        None,
        "unrecognized alignment format (expected MAFFT .mafft/.afa/.aln/.phy, FASTA, Clustal, PHYLIP, or GenBank)",
    ))
}
