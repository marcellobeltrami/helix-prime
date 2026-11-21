//! PCR & sequencing primer thermodynamic analysis, automated primer pair design,
//! and template binding-site testing.

use crate::sequence::{complement_nucleotide, is_gap, reverse_complement, Strand};
use serde::{Deserialize, Serialize};

/// Thermodynamic & structural quality metrics for a single oligonucleotide primer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrimerMetrics {
    pub sequence: Vec<u8>,
    pub length: usize,
    /// Estimated melting temperature ($T_m$) in degrees Celsius.
    pub tm_celsius: f32,
    /// GC percentage (`0.0..=100.0`).
    pub gc_pct: f32,
    /// `true` if the 3' terminal base is `G` or `C` (improves polymerase anchoring).
    pub has_gc_clamp: bool,
    /// Longest internal self-complementary stem (heuristic hairpin/self-dimer risk indicator).
    pub max_self_comp_run: usize,
    /// Overall heuristic quality score (`0.0..=100.0`, higher is better).
    pub quality_score: f32,
}

/// A designed forward + reverse PCR primer pair targeting an amplicon span.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DesignedPrimerPair {
    pub forward: PrimerMetrics,
    pub reverse: PrimerMetrics,
    /// 0-based inclusive start column of the forward primer on the template.
    pub fwd_start: usize,
    /// 0-based exclusive end column of the forward primer on the template.
    pub fwd_end: usize,
    /// 0-based inclusive start column of the reverse primer on the template.
    pub rev_start: usize,
    /// 0-based exclusive end column of the reverse primer on the template.
    pub rev_end: usize,
    /// Expected PCR product length in base pairs.
    pub amplicon_len: usize,
    /// Absolute $T_m$ difference between forward and reverse primers.
    pub tm_diff: f32,
}

/// Constraints governing automated PCR primer design.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrimerDesignConfig {
    pub min_len: usize,
    pub max_len: usize,
    pub target_tm: f32,
    pub min_gc_pct: f32,
    pub max_gc_pct: f32,
}

impl Default for PrimerDesignConfig {
    fn default() -> Self {
        Self {
            min_len: 18,
            max_len: 24,
            target_tm: 60.0,
            min_gc_pct: 40.0,
            max_gc_pct: 62.0,
        }
    }
}

/// A primer binding match found when testing an existing oligo against a DNA template.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerBindingHit {
    pub start_col: usize,
    pub end_col: usize,
    pub strand: Strand,
    pub mismatches: usize,
}

/// Computes salt-adjusted $T_m$, GC %, 3' GC clamp, self-complementarity, and quality score
/// for a DNA primer sequence.
#[must_use]
pub fn analyze_primer(raw: &[u8]) -> PrimerMetrics {
    let seq: Vec<u8> = raw
        .iter()
        .copied()
        .filter(|b| !b.is_ascii_whitespace() && !is_gap(*b))
        .map(|b| b.to_ascii_uppercase())
        .collect();

    let len = seq.len();
    if len == 0 {
        return PrimerMetrics {
            sequence: seq,
            length: 0,
            tm_celsius: 0.0,
            gc_pct: 0.0,
            has_gc_clamp: false,
            max_self_comp_run: 0,
            quality_score: 0.0,
        };
    }

    let gc = seq.iter().filter(|&&b| matches!(b, b'G' | b'C')).count();
    let at = len.saturating_sub(gc);
    let gc_pct = (gc as f32 / len as f32) * 100.0;

    // Wallace rule for short oligos (< 14 bp); standard salt-adjusted GC formula for >= 14 bp
    // assuming 50 mM monovalent Na+ (`64.9 + 41 * (G+C - 16.4) / len`).
    let tm_celsius = if len < 14 {
        (2 * at + 4 * gc) as f32
    } else {
        64.9 + 41.0 * (gc as f32 - 16.4) / (len as f32)
    };

    let has_gc_clamp = seq.last().is_some_and(|&b| matches!(b, b'G' | b'C'));
    let max_self_comp_run = longest_self_complement_run(&seq);

    // Heuristic 0..=100 quality score rewarding ~60 C Tm, ~50% GC, 18-24 bp, GC clamp, low dimer.
    let mut score = 100.0f32;
    score -= (tm_celsius - 60.0).abs() * 2.4;
    score -= (gc_pct - 50.0).abs() * 0.8;
    if !(18..=25).contains(&len) {
        score -= 12.0;
    }
    if !has_gc_clamp {
        score -= 6.0;
    }
    if max_self_comp_run >= 5 {
        score -= ((max_self_comp_run - 4) as f32) * 9.0;
    }

    PrimerMetrics {
        sequence: seq,
        length: len,
        tm_celsius,
        gc_pct,
        has_gc_clamp,
        max_self_comp_run,
        quality_score: score.clamp(0.0, 100.0),
    }
}

/// Designs top-scoring Forward & Reverse PCR primer pairs flanking `[region_start, region_end)`
/// on an aligned DNA template sequence.
#[must_use]
pub fn design_primer_pairs(
    aligned_template: &[u8],
    region_start: usize,
    region_end: usize,
    config: PrimerDesignConfig,
    max_results: usize,
) -> Vec<DesignedPrimerPair> {
    let mut ungapped: Vec<(usize, u8)> = Vec::with_capacity(aligned_template.len());
    for (col, &b) in aligned_template.iter().enumerate() {
        if !is_gap(b) && matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U') {
            let base = if b.eq_ignore_ascii_case(&b'U') {
                b'T'
            } else {
                b.to_ascii_uppercase()
            };
            ungapped.push((col, base));
        }
    }

    let n = ungapped.len();
    let min_l = config.min_len.clamp(12, 36);
    let max_l = config.max_len.clamp(min_l, 40);
    if n < min_l * 2 + 10 {
        return Vec::new();
    }

    // Locate ungapped indices overlapping [region_start, region_end).
    let u_start = ungapped
        .iter()
        .position(|&(col, _)| col >= region_start)
        .unwrap_or(0);
    let u_end = ungapped
        .iter()
        .rposition(|&(col, _)| col < region_end)
        .map_or(n, |p| (p + 1).min(n));

    if u_end.saturating_sub(u_start) < min_l * 2 + 4 {
        return Vec::new();
    }

    // Sample candidate forward primers near `u_start` and reverse primers near `u_end`.
    let scan_window = ((u_end - u_start) / 4).clamp(1, 60);
    let mut fwd_candidates: Vec<(usize, usize, PrimerMetrics)> = Vec::new();
    let mut rev_candidates: Vec<(usize, usize, PrimerMetrics)> = Vec::new();

    for offset in 0..scan_window {
        let f_pos = u_start + offset;
        for len in min_l..=max_l {
            if f_pos + len > u_end {
                break;
            }
            let oligo: Vec<u8> = ungapped[f_pos..(f_pos + len)]
                .iter()
                .map(|&(_, b)| b)
                .collect();
            let metrics = analyze_primer(&oligo);
            if metrics.gc_pct >= config.min_gc_pct
                && metrics.gc_pct <= config.max_gc_pct
                && (metrics.tm_celsius - config.target_tm).abs() <= 10.0
            {
                fwd_candidates.push((f_pos, f_pos + len, metrics));
            }
        }

        let r_end = u_end.saturating_sub(offset);
        for len in min_l..=max_l {
            if r_end < u_start + len {
                break;
            }
            let r_start = r_end - len;
            let sense: Vec<u8> = ungapped[r_start..r_end].iter().map(|&(_, b)| b).collect();
            let rev_oligo = reverse_complement(&sense);
            let metrics = analyze_primer(&rev_oligo);
            if metrics.gc_pct >= config.min_gc_pct
                && metrics.gc_pct <= config.max_gc_pct
                && (metrics.tm_celsius - config.target_tm).abs() <= 10.0
            {
                rev_candidates.push((r_start, r_end, metrics));
            }
        }
    }

    // Fallback: relax GC/Tm bounds if strict filter yielded no candidates (e.g., AT-rich phage loci).
    if fwd_candidates.is_empty() && u_start + min_l <= u_end {
        let oligo: Vec<u8> = ungapped[u_start..(u_start + min_l)]
            .iter()
            .map(|&(_, b)| b)
            .collect();
        fwd_candidates.push((u_start, u_start + min_l, analyze_primer(&oligo)));
    }
    if rev_candidates.is_empty() && u_end >= u_start + min_l {
        let r_start = u_end - min_l;
        let sense: Vec<u8> = ungapped[r_start..u_end].iter().map(|&(_, b)| b).collect();
        rev_candidates.push((r_start, u_end, analyze_primer(&reverse_complement(&sense))));
    }

    fwd_candidates.sort_by(|a, b| b.2.quality_score.total_cmp(&a.2.quality_score));
    rev_candidates.sort_by(|a, b| b.2.quality_score.total_cmp(&a.2.quality_score));
    fwd_candidates.truncate(12);
    rev_candidates.truncate(12);

    let mut pairs = Vec::new();
    for (fs, fe, fm) in &fwd_candidates {
        for (rs, re, rm) in &rev_candidates {
            if *fe >= *rs {
                continue;
            }
            let tm_diff = (fm.tm_celsius - rm.tm_celsius).abs();
            let amplicon_len = re.saturating_sub(*fs);
            pairs.push(DesignedPrimerPair {
                forward: fm.clone(),
                reverse: rm.clone(),
                fwd_start: ungapped[*fs].0,
                fwd_end: ungapped[fe - 1].0.saturating_add(1),
                rev_start: ungapped[*rs].0,
                rev_end: ungapped[re - 1].0.saturating_add(1),
                amplicon_len,
                tm_diff,
            });
        }
    }

    // Rank pairs by combined primer quality minus Tm mismatch penalty.
    pairs.sort_by(|a, b| {
        let score_a = a.forward.quality_score + a.reverse.quality_score - a.tm_diff * 5.0;
        let score_b = b.forward.quality_score + b.reverse.quality_score - b.tm_diff * 5.0;
        score_b.total_cmp(&score_a)
    });
    pairs.truncate(max_results.max(1));
    pairs
}

/// Tests where a given primer binds on a DNA template (both Forward and Reverse strands),
/// allowing up to `max_mismatches` substitutions.
#[must_use]
pub fn test_primer_binding(
    aligned_template: &[u8],
    primer: &[u8],
    max_mismatches: usize,
) -> Vec<PrimerBindingHit> {
    let clean_primer: Vec<u8> = primer
        .iter()
        .copied()
        .filter(|b| !b.is_ascii_whitespace() && !is_gap(*b))
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if clean_primer.len() < 6 {
        return Vec::new();
    }
    let rev_primer = reverse_complement(&clean_primer);

    let mut ungapped: Vec<(usize, u8)> = Vec::with_capacity(aligned_template.len());
    for (col, &b) in aligned_template.iter().enumerate() {
        if !is_gap(b) {
            ungapped.push((col, b.to_ascii_uppercase()));
        }
    }
    let plen = clean_primer.len();
    if ungapped.len() < plen {
        return Vec::new();
    }

    let mut hits = Vec::new();
    for i in 0..=(ungapped.len() - plen) {
        let mut fwd_mm = 0usize;
        let mut rev_mm = 0usize;
        for k in 0..plen {
            let tb = ungapped[i + k].1;
            if tb != clean_primer[k] {
                fwd_mm += 1;
            }
            if tb != rev_primer[k] {
                rev_mm += 1;
            }
        }
        if fwd_mm <= max_mismatches {
            hits.push(PrimerBindingHit {
                start_col: ungapped[i].0,
                end_col: ungapped[i + plen - 1].0.saturating_add(1),
                strand: Strand::Forward,
                mismatches: fwd_mm,
            });
        }
        if rev_mm <= max_mismatches {
            hits.push(PrimerBindingHit {
                start_col: ungapped[i].0,
                end_col: ungapped[i + plen - 1].0.saturating_add(1),
                strand: Strand::Reverse,
                mismatches: rev_mm,
            });
        }
    }

    hits
}

/// Computes the longest contiguous reverse-complementary run inside `seq` (self-dimer / hairpin check).
fn longest_self_complement_run(seq: &[u8]) -> usize {
    let n = seq.len();
    if n < 4 {
        return 0;
    }
    let rc: Vec<u8> = seq.iter().rev().map(|&b| complement_nucleotide(b)).collect();
    let mut best = 0usize;
    for shift in -(n as isize - 1)..(n as isize) {
        let mut run = 0usize;
        for (i, &sb) in seq.iter().enumerate() {
            let j = i as isize + shift;
            if j >= 0 && (j as usize) < n && sb == rc[j as usize] {
                run += 1;
                if run > best {
                    best = run;
                }
            } else {
                run = 0;
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primer_design_and_binding_find_valid_pcr_pair() {
        let template = concat!(
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC",
            "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG",
            "AGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAG"
        )
        .as_bytes();
        let pairs = design_primer_pairs(template, 0, template.len(), PrimerDesignConfig::default(), 5);
        assert!(!pairs.is_empty());
        let top = &pairs[0];
        assert!(top.amplicon_len >= 80);
        assert!(top.forward.tm_celsius > 45.0);

        let hits = test_primer_binding(template, &top.forward.sequence, 0);
        assert!(hits.iter().any(|h| h.strand == Strand::Forward && h.start_col == top.fwd_start));
        let rev_hits = test_primer_binding(template, &top.reverse.sequence, 0);
        assert!(rev_hits.iter().any(|h| h.strand == Strand::Reverse && h.start_col == top.rev_start));
    }
}
