//! Restriction enzyme library, circular-aware cleavage site finder, and in-silico digestion.

use crate::sequence::{complement_nucleotide, is_gap};
use serde::{Deserialize, Serialize};

/// Type of DNA overhang left by a restriction endonuclease cleavage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum OverhangKind {
    FivePrime,
    ThreePrime,
    Blunt,
}

/// Static definition of a Type II restriction enzyme recognition site.
///
/// # Invariants
/// - `recognition` contains uppercase IUPAC DNA bases (`A`, `C`, `G`, `T`, `N`, `R`, `Y`, ...).
/// - `top_cut <= recognition.len()` and `bottom_cut <= recognition.len()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RestrictionEnzyme {
    pub name: &'static str,
    pub recognition: &'static [u8],
    /// 0-based cleavage offset on the forward (5'->3') strand relative to site start.
    pub top_cut: usize,
    /// 0-based cleavage offset on the reverse strand relative to site start.
    pub bottom_cut: usize,
}

impl RestrictionEnzyme {
    /// Classifies the sticky/blunt end produced by cutting with this enzyme.
    #[must_use]
    pub const fn overhang(&self) -> OverhangKind {
        if self.top_cut < self.bottom_cut {
            OverhangKind::FivePrime
        } else if self.top_cut > self.bottom_cut {
            OverhangKind::ThreePrime
        } else {
            OverhangKind::Blunt
        }
    }

    /// Recognition sequence formatted as UTF-8 string (e.g. `"GAATTC"`).
    #[must_use]
    pub fn recognition_str(&self) -> &'static str {
        std::str::from_utf8(self.recognition).unwrap_or("")
    }
}

/// Curated library of 24 commonly used molecular cloning & BioBrick restriction enzymes.
pub const COMMON_ENZYMES: &[RestrictionEnzyme] = &[
    RestrictionEnzyme { name: "EcoRI",  recognition: b"GAATTC",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "BamHI",  recognition: b"GGATCC",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "HindIII",recognition: b"AAGCTT",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "XbaI",   recognition: b"TCTAGA",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "SpeI",   recognition: b"ACTAGT",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "PstI",   recognition: b"CTGCAG",   top_cut: 5, bottom_cut: 1 },
    RestrictionEnzyme { name: "NotI",   recognition: b"GCGGCCGC", top_cut: 2, bottom_cut: 6 },
    RestrictionEnzyme { name: "XhoI",   recognition: b"CTCGAG",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "SalI",   recognition: b"GTCGAC",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "KpnI",   recognition: b"GGTACC",   top_cut: 5, bottom_cut: 1 },
    RestrictionEnzyme { name: "SacI",   recognition: b"GAGCTC",   top_cut: 5, bottom_cut: 1 },
    RestrictionEnzyme { name: "SmaI",   recognition: b"CCCGGG",   top_cut: 3, bottom_cut: 3 },
    RestrictionEnzyme { name: "EcoRV",  recognition: b"GATATC",   top_cut: 3, bottom_cut: 3 },
    RestrictionEnzyme { name: "NcoI",   recognition: b"CCATGG",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "NdeI",   recognition: b"CATATG",   top_cut: 2, bottom_cut: 4 },
    RestrictionEnzyme { name: "BglII",  recognition: b"AGATCT",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "NheI",   recognition: b"GCTAGC",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "AgeI",   recognition: b"ACCGGT",   top_cut: 1, bottom_cut: 5 },
    RestrictionEnzyme { name: "ClaI",   recognition: b"ATCGAT",   top_cut: 2, bottom_cut: 4 },
    RestrictionEnzyme { name: "SphI",   recognition: b"GCATGC",   top_cut: 5, bottom_cut: 1 },
    RestrictionEnzyme { name: "ApaI",   recognition: b"GGGCCC",   top_cut: 5, bottom_cut: 1 },
    RestrictionEnzyme { name: "BsaI",   recognition: b"GGTCTC",   top_cut: 7, bottom_cut: 11 },
    RestrictionEnzyme { name: "PacI",   recognition: b"TTAATTAA", top_cut: 5, bottom_cut: 3 },
    RestrictionEnzyme { name: "AscI",   recognition: b"GGCGCGCC", top_cut: 2, bottom_cut: 6 },
];

/// A discovered restriction enzyme site anchored to alignment column coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RestrictionSiteHit {
    pub enzyme_name: &'static str,
    /// 0-based inclusive start column in the alignment.
    pub start_col: usize,
    /// 0-based exclusive end column in the alignment.
    pub end_col: usize,
    /// 0-based alignment column immediately after the top-strand cleavage point.
    pub cut_col: usize,
    pub overhang: OverhangKind,
}

/// Looks up a built-in enzyme by name (case-insensitive).
#[must_use]
pub fn find_enzyme_by_name(name: &str) -> Option<&'static RestrictionEnzyme> {
    COMMON_ENZYMES
        .iter()
        .find(|e| e.name.eq_ignore_ascii_case(name.trim()))
}

/// Scans an aligned nucleotide sequence for all recognition sites of `enzymes`, skipping gaps.
///
/// When `circular` is `true`, sites spanning the circular junction are also detected.
#[must_use]
pub fn find_restriction_sites(
    aligned: &[u8],
    enzymes: &[RestrictionEnzyme],
    circular: bool,
) -> Vec<RestrictionSiteHit> {
    let mut ungapped: Vec<(usize, u8)> = Vec::with_capacity(aligned.len());
    for (col, &b) in aligned.iter().enumerate() {
        if !is_gap(b) {
            ungapped.push((col, b.to_ascii_uppercase()));
        }
    }
    let n = ungapped.len();
    if n < 4 {
        return Vec::new();
    }

    let mut hits = Vec::new();

    for enz in enzymes {
        let pat = enz.recognition;
        let plen = pat.len();
        if plen == 0 || plen > n {
            continue;
        }
        let max_start = if circular { n } else { n - plen + 1 };

        // Precompute reverse complement pattern for non-palindromic enzymes (e.g. BsaI).
        let rev_pat: Vec<u8> = pat.iter().rev().map(|&b| complement_nucleotide(b)).collect();
        let is_palindrome = rev_pat.as_slice() == pat;

        for start_i in 0..max_start {
            let fwd_match = (0..plen).all(|k| {
                let (_, base) = ungapped[(start_i + k) % n];
                iupac_match(pat[k], base)
            });

            let rev_match = !fwd_match
                && !is_palindrome
                && (0..plen).all(|k| {
                    let (_, base) = ungapped[(start_i + k) % n];
                    iupac_match(rev_pat[k], base)
                });

            if fwd_match || rev_match {
                let start_col = ungapped[start_i].0;
                let end_col = ungapped[(start_i + plen - 1) % n].0.saturating_add(1);
                let cut_offset = enz.top_cut.min(plen.saturating_sub(1));
                let cut_col = ungapped[(start_i + cut_offset) % n].0;
                hits.push(RestrictionSiteHit {
                    enzyme_name: enz.name,
                    start_col,
                    end_col: end_col.max(start_col + 1),
                    cut_col,
                    overhang: enz.overhang(),
                });
            }
        }
    }

    hits.sort_by_key(|h| (h.cut_col, h.enzyme_name));
    hits
}

/// Checks if an IUPAC pattern byte matches a concrete DNA/RNA base byte.
#[inline]
fn iupac_match(pat: u8, base: u8) -> bool {
    let b = if base == b'U' { b'T' } else { base };
    match pat {
        b'A' | b'C' | b'G' | b'T' => pat == b,
        b'N' => matches!(b, b'A' | b'C' | b'G' | b'T'),
        b'R' => matches!(b, b'A' | b'G'),
        b'Y' => matches!(b, b'C' | b'T'),
        b'S' => matches!(b, b'G' | b'C'),
        b'W' => matches!(b, b'A' | b'T'),
        b'K' => matches!(b, b'G' | b'T'),
        b'M' => matches!(b, b'A' | b'C'),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restriction_scanner_finds_linear_and_circular_sites() {
        // EcoRI (GAATTC) at index 2, BamHI (GGATCC) at index 10.
        let seq = b"ATGAATTCGGGGATCCTT";
        let hits = find_restriction_sites(seq, COMMON_ENZYMES, false);
        assert!(hits.iter().any(|h| h.enzyme_name == "EcoRI" && h.start_col == 2));
        assert!(hits.iter().any(|h| h.enzyme_name == "BamHI" && h.start_col == 10));

        // Circular EcoRI spanning end->start ("...GAA" + "TTC...").
        let circ = b"TTCAAAAAAAAAAAGAA";
        let circ_hits = find_restriction_sites(circ, COMMON_ENZYMES, true);
        assert!(circ_hits.iter().any(|h| h.enzyme_name == "EcoRI"));
    }
}
