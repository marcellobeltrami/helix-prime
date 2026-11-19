//! Biological sequence representation, alphabet inference, and genetic code translation.
//!
//! Sequences are stored internally as contiguous uppercase ASCII byte buffers (`Vec<u8>`)
//! so residue indexing is O(1), cache-friendly, and free of UTF-8 decoding overhead.

use serde::{Deserialize, Serialize};

/// Biological alphabet of a sequence or multiple sequence alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Alphabet {
    /// Deoxyribonucleic acid (`A`, `C`, `G`, `T`, plus IUPAC ambiguity codes).
    Dna,
    /// Ribonucleic acid (`A`, `C`, `G`, `U`, plus IUPAC ambiguity codes).
    Rna,
    /// Amino acid polypeptide (`A`..`Z`, `*` for stop).
    Protein,
}

impl Alphabet {
    /// Short display label used in the GUI document table and status bar.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
            Self::Protein => "Protein",
        }
    }

    /// Returns `true` if this alphabet represents nucleotides (DNA or RNA).
    #[must_use]
    pub const fn is_nucleotide(self) -> bool {
        matches!(self, Self::Dna | Self::Rna)
    }

    /// Infers the most likely biological alphabet from a slice of raw residue slices.
    ///
    /// Inspects non-gap ASCII letters: if unambiguously protein-only amino acids appear
    /// (`E`, `F`, `I`, `L`, `P`, `Q`, `Z`), classifies as [`Alphabet::Protein`]. Otherwise
    /// distinguishes RNA from DNA by comparing `U` vs `T` counts.
    #[must_use]
    pub fn infer_from_slices<'a>(slices: impl IntoIterator<Item = &'a [u8]>) -> Self {
        let mut total_letters: usize = 0;
        let mut standard_nuc: usize = 0;
        let mut protein_only: usize = 0;
        let mut u_count: usize = 0;
        let mut t_count: usize = 0;

        for bytes in slices {
            for &raw in bytes {
                let ch = raw.to_ascii_uppercase();
                if is_gap(ch) || !ch.is_ascii_alphabetic() {
                    continue;
                }
                total_letters += 1;
                match ch {
                    b'A' | b'C' | b'G' | b'N' => standard_nuc += 1,
                    b'T' => {
                        standard_nuc += 1;
                        t_count += 1;
                    }
                    b'U' => {
                        standard_nuc += 1;
                        u_count += 1;
                    }
                    b'E' | b'F' | b'I' | b'L' | b'P' | b'Q' | b'J' | b'O' => {
                        protein_only += 1;
                    }
                    _ => {}
                }
            }
        }

        if total_letters == 0 {
            return Self::Dna;
        }

        // Even a small fraction of strictly amino-acid letters rules out DNA/RNA.
        if protein_only > 0 || (standard_nuc * 100 / total_letters) < 80 {
            Self::Protein
        } else if u_count > t_count {
            Self::Rna
        } else {
            Self::Dna
        }
    }
}

/// Strand orientation for nucleotide annotations and translations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum Strand {
    /// Forward (5' -> 3' / sense) strand.
    #[default]
    Forward,
    /// Reverse (3' -> 5' / antisense) strand.
    Reverse,
}

impl Strand {
    /// Single-character GFF3/GenBank strand symbol.
    #[must_use]
    pub const fn symbol(self) -> char {
        match self {
            Self::Forward => '+',
            Self::Reverse => '-',
        }
    }
}

/// Returns `true` if `byte` represents an alignment gap (`-` or `.`).
#[inline]
#[must_use]
pub const fn is_gap(byte: u8) -> bool {
    matches!(byte, b'-' | b'.')
}

/// A single biological sequence (aligned or unaligned) with optional UniProt metadata.
///
/// # Invariants
/// - `residues` contains uppercase ASCII symbols (`A`–`Z`, `-`, `.`, `*`, `?`).
/// - `name` is non-empty after normalization (callers fall back to `"Unnamed"`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sequence {
    /// Short identifier displayed in the alignment row header (e.g., `P53_HUMAN`).
    pub name: String,
    /// Optional longer description from FASTA/GenBank headers.
    pub description: String,
    /// Optional organism scientific name (e.g., `Homo sapiens`).
    pub organism: Option<String>,
    /// Optional UniProtKB accession extracted from header or linked manually (e.g., `P04637`).
    pub uniprot_accession: Option<String>,
    /// Uppercase ASCII residue bytes, including gap symbols when aligned.
    pub residues: Vec<u8>,
}

impl Sequence {
    /// Creates a new normalized sequence, stripping whitespace and uppercasing ASCII residues.
    #[must_use]
    pub fn new(name: impl Into<String>, raw_residues: &[u8]) -> Self {
        let mut name_str = name.into();
        if name_str.trim().is_empty() {
            name_str = String::from("Unnamed");
        }
        let mut residues = Vec::with_capacity(raw_residues.len());
        for &b in raw_residues {
            if !b.is_ascii_whitespace() {
                residues.push(b.to_ascii_uppercase());
            }
        }
        Self {
            name: name_str,
            description: String::new(),
            organism: None,
            uniprot_accession: None,
            residues,
        }
    }

    /// Builder helper to attach a description string.
    #[must_use]
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = description.into();
        self
    }

    /// Builder helper to attach an organism name.
    #[must_use]
    pub fn with_organism(mut self, organism: impl Into<String>) -> Self {
        self.organism = Some(organism.into());
        self
    }

    /// Builder helper to attach a UniProt accession.
    #[must_use]
    pub fn with_uniprot(mut self, accession: impl Into<String>) -> Self {
        self.uniprot_accession = Some(accession.into());
        self
    }

    /// Total length including alignment gap characters.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// Returns `true` if the sequence has zero columns.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Counts non-gap residues in this sequence.
    #[must_use]
    pub fn ungapped_len(&self) -> usize {
        self.residues.iter().filter(|&&b| !is_gap(b)).count()
    }

    /// Returns a new byte vector containing only the non-gap residues.
    #[must_use]
    pub fn ungapped_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.residues.len());
        for &b in &self.residues {
            if !is_gap(b) {
                out.push(b);
            }
        }
        out
    }

    /// Returns the 1-based ungapped residue number at alignment column `col`,
    /// or `None` if `col` is out of bounds or holds a gap character.
    #[must_use]
    pub fn ungapped_position_at(&self, col: usize) -> Option<usize> {
        let byte = *self.residues.get(col)?;
        if is_gap(byte) {
            return None;
        }
        let mut count = 0usize;
        for &b in &self.residues[..=col] {
            if !is_gap(b) {
                count += 1;
            }
        }
        Some(count)
    }

    /// Computes GC percentage (`0.0..=100.0`) over non-gap nucleotide residues.
    /// Returns `None` if there are no non-gap bases.
    #[must_use]
    pub fn gc_content_pct(&self) -> Option<f32> {
        let mut total = 0usize;
        let mut gc = 0usize;
        for &b in &self.residues {
            if is_gap(b) {
                continue;
            }
            total += 1;
            if matches!(b, b'G' | b'C') {
                gc += 1;
            }
        }
        if total == 0 {
            None
        } else {
            Some((gc as f32 / total as f32) * 100.0)
        }
    }

    /// Returns the residue byte at `col`, if within bounds.
    #[inline]
    #[must_use]
    pub fn residue_at(&self, col: usize) -> Option<u8> {
        self.residues.get(col).copied()
    }
}

/// Returns the Watson-Crick complement of a single IUPAC nucleotide byte.
#[must_use]
pub const fn complement_nucleotide(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

/// Computes the reverse complement of an ungapped or gapped nucleotide slice.
#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(complement_nucleotide(b));
    }
    out
}

/// Translates a single 3-base DNA/RNA codon into a 1-letter amino acid (NCBI Standard Code 1).
///
/// Returns `b'*'` for stop codons and `b'X'` for ambiguous or gapped codons.
#[must_use]
pub fn translate_codon(b0: u8, b1: u8, b2: u8) -> u8 {
    let idx = |b: u8| -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'T' | b'U' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    };

    // Standard NCBI Translation Table 1 indexed by (b0 * 16 + b1 * 4 + b2), order T, C, A, G.
    const TABLE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    match (idx(b0), idx(b1), idx(b2)) {
        (Some(i0), Some(i1), Some(i2)) => TABLE[i0 * 16 + i1 * 4 + i2],
        _ => b'X',
    }
}

/// Translates an ungapped nucleotide slice starting at reading frame offset `frame` (`0..=2`).
///
/// Any trailing incomplete codon (< 3 bases) is ignored.
#[must_use]
pub fn translate_ungapped(seq: &[u8], frame: usize) -> Vec<u8> {
    let frame = frame % 3;
    if seq.len() < frame + 3 {
        return Vec::new();
    }
    let slice = &seq[frame..];
    let codon_count = slice.len() / 3;
    let mut protein = Vec::with_capacity(codon_count);
    for chunk in slice.chunks_exact(3) {
        protein.push(translate_codon(chunk[0], chunk[1], chunk[2]));
    }
    protein
}

/// A translated codon anchored back to its alignment column span `[start_col, end_col)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignedCodon {
    /// 0-based inclusive alignment column of the first base of this codon.
    pub start_col: usize,
    /// 0-based exclusive alignment column past the third base of this codon.
    pub end_col: usize,
    /// 0-based alignment column of the middle base (ideal for centering the amino acid label).
    pub center_col: usize,
    /// Translated 1-letter amino acid (`A`–`Z`, `*`, or `X`).
    pub amino_acid: u8,
    /// 1-based translated residue index within the reading frame.
    pub residue_index: usize,
}

/// Translates an aligned nucleotide slice over `[range_start, range_end)` while skipping
/// alignment gaps so frameshifts caused by indels in other sequences do not corrupt codons.
///
/// Respects `strand`: for [`Strand::Reverse`], codons are read from `range_end` down to
/// `range_start` using complementary bases.
#[must_use]
pub fn translate_aligned_region(
    aligned: &[u8],
    range_start: usize,
    range_end: usize,
    strand: Strand,
    frame: usize,
) -> Vec<AlignedCodon> {
    let end = range_end.min(aligned.len());
    let start = range_start.min(end);
    if start >= end {
        return Vec::new();
    }

    // Collect non-gap (column_index, base) pairs in transcription order.
    let mut bases: Vec<(usize, u8)> = Vec::with_capacity(end - start);
    match strand {
        Strand::Forward => {
            for (col, &b) in aligned.iter().enumerate().take(end).skip(start) {
                if !is_gap(b) {
                    bases.push((col, b));
                }
            }
        }
        Strand::Reverse => {
            for col in (start..end).rev() {
                let b = aligned[col];
                if !is_gap(b) {
                    bases.push((col, complement_nucleotide(b)));
                }
            }
        }
    }

    let offset = frame % 3;
    if bases.len() < offset + 3 {
        return Vec::new();
    }

    let usable = &bases[offset..];
    let mut codons = Vec::with_capacity(usable.len() / 3);
    for (idx, chunk) in usable.chunks_exact(3).enumerate() {
        let (c0, b0) = chunk[0];
        let (c1, b1) = chunk[1];
        let (c2, b2) = chunk[2];
        let min_col = c0.min(c2);
        let max_col = c0.max(c2);
        codons.push(AlignedCodon {
            start_col: min_col,
            end_col: max_col.saturating_add(1),
            center_col: c1,
            amino_acid: translate_codon(b0, b1, b2),
            residue_index: idx + 1,
        });
    }
    codons
}

/// Converts a 3-letter PDB/mmCIF amino acid code (e.g. `"ALA"`) into a 1-letter IUPAC code.
#[must_use]
pub fn three_to_one_aa(code: &str) -> u8 {
    match code.trim().to_ascii_uppercase().as_str() {
        "ALA" => b'A',
        "ARG" => b'R',
        "ASN" => b'N',
        "ASP" => b'D',
        "CYS" => b'C',
        "GLN" => b'Q',
        "GLU" => b'E',
        "GLY" => b'G',
        "HIS" => b'H',
        "ILE" => b'I',
        "LEU" => b'L',
        "LYS" => b'K',
        "MET" | "MSE" => b'M',
        "PHE" => b'F',
        "PRO" => b'P',
        "SER" => b'S',
        "THR" => b'T',
        "TRP" => b'W',
        "TYR" => b'Y',
        "VAL" => b'V',
        "SEC" => b'U',
        "PYL" => b'O',
        _ => b'X',
    }
}

/// Converts a 1-letter amino acid code into its standard 3-letter abbreviation.
#[must_use]
pub const fn one_to_three_aa(code: u8) -> &'static str {
    match code.to_ascii_uppercase() {
        b'A' => "ALA",
        b'R' => "ARG",
        b'N' => "ASN",
        b'D' => "ASP",
        b'C' => "CYS",
        b'Q' => "GLN",
        b'E' => "GLU",
        b'G' => "GLY",
        b'H' => "HIS",
        b'I' => "ILE",
        b'L' => "LEU",
        b'K' => "LYS",
        b'M' => "MET",
        b'F' => "PHE",
        b'P' => "PRO",
        b'S' => "SER",
        b'T' => "THR",
        b'W' => "TRP",
        b'Y' => "TYR",
        b'V' => "VAL",
        b'*' => "TER",
        _ => "UNK",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alphabet_inference_distinguishes_dna_rna_protein() {
        assert_eq!(
            Alphabet::infer_from_slices([b"ATGCGATCGAT--CG".as_slice()]),
            Alphabet::Dna
        );
        assert_eq!(
            Alphabet::infer_from_slices([b"AUGCGAUCGAU--CG".as_slice()]),
            Alphabet::Rna
        );
        assert_eq!(
            Alphabet::infer_from_slices([b"MEEPQSDPSVEPPLSQETFSDLWKLLPENN".as_slice()]),
            Alphabet::Protein
        );
    }

    #[test]
    fn translation_handles_standard_codons_and_reverse_strand() {
        assert_eq!(translate_ungapped(b"ATGGTGCACCTGACTCCTGAGGAGAAG", 0), b"MVHLTPEEK");
        // Reverse complement of "ATGGTGCAC" is "GTGCACCAT" -> VHH
        let aligned = b"ATG--GTGCAC";
        let fwd = translate_aligned_region(aligned, 0, aligned.len(), Strand::Forward, 0);
        assert_eq!(fwd.len(), 3);
        assert_eq!(fwd[0].amino_acid, b'M');
        assert_eq!(fwd[1].amino_acid, b'V');
        assert_eq!(fwd[2].amino_acid, b'H');

        let rev = translate_aligned_region(aligned, 0, aligned.len(), Strand::Reverse, 0);
        assert_eq!(rev.len(), 3);
        assert_eq!(rev[0].amino_acid, b'V');
        assert_eq!(rev[1].amino_acid, b'H');
        assert_eq!(rev[2].amino_acid, b'H');
    }

    #[test]
    fn ungapped_position_skips_gaps_safely() {
        let seq = Sequence::new("test", b"--AT-G");
        assert_eq!(seq.ungapped_position_at(0), None);
        assert_eq!(seq.ungapped_position_at(2), Some(1));
        assert_eq!(seq.ungapped_position_at(3), Some(2));
        assert_eq!(seq.ungapped_position_at(4), None);
        assert_eq!(seq.ungapped_position_at(5), Some(3));
        assert_eq!(seq.ungapped_position_at(99), None);
    }
}
