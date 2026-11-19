//! Unified bioinformatics file format detection, loading, and export.

pub mod clustal;
pub mod fasta;
pub mod genbank;
pub mod gff3;

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use std::path::Path;

/// Supported alignment and sequence input formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentFormat {
    Fasta,
    Clustal,
    GenBank,
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
    }
}

/// Detects sequence/alignment format from extension first, falling back to content sniffing.
pub fn detect_format(path: &Path, content: &str) -> Result<AlignmentFormat> {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "fa" | "fasta" | "fna" | "faa" | "fas" => return Ok(AlignmentFormat::Fasta),
            "aln" | "clustal" => return Ok(AlignmentFormat::Clustal),
            "gb" | "gbk" | "genbank" => return Ok(AlignmentFormat::GenBank),
            _ => {}
        }
    }

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
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
        break;
    }

    Err(HelixError::parse(
        "AutoDetect",
        None,
        "unrecognized alignment format (expected FASTA, Clustal .aln, or GenBank .gb/.gbk)",
    ))
}
