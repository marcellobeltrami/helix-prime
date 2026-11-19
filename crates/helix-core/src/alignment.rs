//! Multiple sequence alignment data model, column statistics, motif search, and progressive aligner.

use crate::annotation::Annotation;
use crate::error::{HelixError, Result};
use crate::sequence::{is_gap, Alphabet, Sequence};
use serde::{Deserialize, Serialize};

/// Configuration controlling how the alignment consensus track is computed.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConsensusConfig {
    /// Minimum fraction (`0.0..=1.0`) of sequences that must agree for a residue to be called.
    /// Below this threshold, an ambiguity symbol (`N` for DNA/RNA, `X` for Protein) is emitted.
    pub threshold: f32,
    /// When `true`, gap characters are ignored when calculating residue frequencies.
    pub ignore_gaps: bool,
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            threshold: 0.5,
            ignore_gaps: true,
        }
    }
}

/// Precomputed per-column statistics used by the alignment renderer and 3D conservation mapper.
///
/// Fields are ordered from largest to smallest alignment to keep the struct compact (`16` bytes).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ColumnStats {
    /// Mean pairwise residue identity percentage (`0.0..=100.0`), matching Geneious Prime's
    /// signature green/yellow/red identity histogram track.
    pub identity_pct: f32,
    /// Shannon information entropy in bits (`0.0` = fully conserved).
    pub shannon_entropy: f32,
    /// Number of sequences having a non-gap symbol in this column.
    pub non_gap_count: u32,
    /// Called consensus byte (`A`–`Z`, `-`, `N`, or `X`).
    pub consensus: u8,
    /// Modal (most frequent non-gap) residue byte, or `b'-'` if the column is all gaps.
    pub modal_residue: u8,
}

/// A motif or subsequence match located within a specific aligned sequence row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MotifHit {
    /// 0-based sequence row index.
    pub row: usize,
    /// 0-based inclusive start column in the alignment.
    pub start_col: usize,
    /// 0-based exclusive end column in the alignment.
    pub end_col: usize,
}

/// A multiple sequence alignment (MSA) with annotations and precomputed column statistics.
///
/// # Invariants
/// - Every `Sequence` in `sequences` has `seq.len() == self.width`.
/// - `columns.len() == self.width`.
/// - `annotations` have `ann.start < ann.end` and `ann.end <= self.width` (when `width > 0`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Alignment {
    /// Document title displayed in the Sources panel and tab bar.
    pub name: String,
    /// Optional description shown in the document table.
    pub description: String,
    /// Biological alphabet (DNA, RNA, or Protein).
    pub alphabet: Alphabet,
    /// Uniform alignment column count across all sequences.
    width: usize,
    /// Aligned sequences (all guaranteed to have length `width`).
    sequences: Vec<Sequence>,
    /// Annotations attached to the alignment or individual sequences.
    annotations: Vec<Annotation>,
    /// Next unique annotation ID counter.
    next_annotation_id: u64,
    /// Consensus calculation parameters.
    consensus_config: ConsensusConfig,
    /// Cached per-column statistics (`columns.len() == width`).
    columns: Box<[ColumnStats]>,
}

impl Alignment {
    /// Builds a new alignment from one or more sequences.
    ///
    /// Shorter sequences are right-padded with gap characters (`'-'`) so the rectangular
    /// alignment invariant (`seq.len() == width`) always holds.
    pub fn new(
        name: impl Into<String>,
        alphabet: Alphabet,
        mut sequences: Vec<Sequence>,
    ) -> Result<Self> {
        if sequences.is_empty() {
            return Err(HelixError::Alignment(
                "cannot create an alignment with zero sequences".to_owned(),
            ));
        }

        let width = sequences.iter().map(Sequence::len).max().unwrap_or(0);
        for seq in &mut sequences {
            if seq.residues.len() < width {
                seq.residues.resize(width, b'-');
            }
        }

        let mut aln = Self {
            name: name.into(),
            description: String::new(),
            alphabet,
            width,
            sequences,
            annotations: Vec::new(),
            next_annotation_id: 1,
            consensus_config: ConsensusConfig::default(),
            columns: Box::new([]),
        };
        aln.recompute_stats();
        Ok(aln)
    }

    /// Builder helper to set document description.
    #[must_use]
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = description.into();
        self
    }

    /// Number of columns in the alignment.
    #[inline]
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// Number of sequences (rows) in the alignment.
    #[inline]
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.sequences.len()
    }

    /// Read-only slice of aligned sequences.
    #[inline]
    #[must_use]
    pub fn sequences(&self) -> &[Sequence] {
        &self.sequences
    }

    /// Mutable slice of aligned sequences (caller should invoke [`Self::recompute_stats`] after
    /// modifying residues).
    #[inline]
    pub fn sequences_mut(&mut self) -> &mut [Sequence] {
        &mut self.sequences
    }

    /// Read-only slice of precomputed per-column statistics.
    #[inline]
    #[must_use]
    pub fn columns(&self) -> &[ColumnStats] {
        &self.columns
    }

    /// Read-only slice of annotations attached to this alignment.
    #[inline]
    #[must_use]
    pub fn annotations(&self) -> &[Annotation] {
        &self.annotations
    }

    /// Mutable access to annotations (e.g., to attach UniProt metadata or edit notes).
    #[inline]
    pub fn annotations_mut(&mut self) -> &mut Vec<Annotation> {
        &mut self.annotations
    }

    /// Returns the current consensus configuration.
    #[inline]
    #[must_use]
    pub const fn consensus_config(&self) -> ConsensusConfig {
        self.consensus_config
    }

    /// Updates consensus settings and recomputes column statistics.
    pub fn set_consensus_config(&mut self, config: ConsensusConfig) {
        self.consensus_config = ConsensusConfig {
            threshold: config.threshold.clamp(0.0, 1.0),
            ignore_gaps: config.ignore_gaps,
        };
        self.recompute_stats();
    }

    /// Adds an annotation, clamping its span to the alignment width and assigning a fresh ID.
    pub fn add_annotation(&mut self, mut ann: Annotation) -> u64 {
        let id = self.next_annotation_id;
        self.next_annotation_id = self.next_annotation_id.saturating_add(1);
        ann.id = id;
        if self.width > 0 {
            ann.start = ann.start.min(self.width.saturating_sub(1));
            ann.end = ann.end.clamp(ann.start + 1, self.width);
        }
        self.annotations.push(ann);
        id
    }

    /// Removes an annotation by its unique ID. Returns `true` if removed.
    pub fn remove_annotation(&mut self, id: u64) -> bool {
        let before = self.annotations.len();
        self.annotations.retain(|a| a.id != id);
        self.annotations.len() < before
    }

    /// Looks up an annotation by ID.
    #[must_use]
    pub fn annotation_by_id(&self, id: u64) -> Option<&Annotation> {
        self.annotations.iter().find(|a| a.id == id)
    }

    /// Looks up a mutable annotation by ID.
    pub fn annotation_by_id_mut(&mut self, id: u64) -> Option<&mut Annotation> {
        self.annotations.iter_mut().find(|a| a.id == id)
    }

    /// Computes mean pairwise identity across the entire alignment (`0.0..=100.0`).
    #[must_use]
    pub fn mean_pairwise_identity_pct(&self) -> f32 {
        if self.columns.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.columns.iter().map(|c| c.identity_pct).sum();
        sum / self.columns.len() as f32
    }

    /// Fraction of columns (`0.0..=100.0`) that are 100% identical across all non-gap rows.
    #[must_use]
    pub fn identical_sites_pct(&self) -> f32 {
        if self.columns.is_empty() {
            return 0.0;
        }
        let identical = self
            .columns
            .iter()
            .filter(|c| c.identity_pct >= 99.99 && c.non_gap_count > 0)
            .count();
        (identical as f32 / self.columns.len() as f32) * 100.0
    }

    /// Computes a circular sliding-window `(gc_fraction, gc_skew)` profile with `bins` angular
    /// buckets around the genome (wrapping seamlessly across the circular origin).
    ///
    /// - `gc_fraction` is in `0.0..=1.0` (proportion of G/C among A/C/G/T/U bases in the window).
    /// - `gc_skew` is `(G - C) / (G + C)` in `-1.0..=1.0` (standard phage/bacterial `ori`/`ter` signal).
    ///
    /// Returns `None` for protein alignments or empty alignments.
    #[must_use]
    pub fn circular_gc_profile(&self, bins: usize) -> Option<Vec<(f32, f32)>> {
        if !self.alphabet.is_nucleotide() || self.width == 0 || bins == 0 {
            return None;
        }
        let ref_seq = self.sequences.first()?;
        let len = ref_seq.residues.len();
        if len == 0 {
            return None;
        }

        let safe_bins = bins.clamp(12, 360);
        // Use a modest sliding window (~1/18th of genome, at least 9 bp) wrapped circularly.
        let half_win = (len / 36).max(4);
        let mut profile = Vec::with_capacity(safe_bins);

        for bin in 0..safe_bins {
            let center = (bin * len) / safe_bins;
            let mut g_cnt = 0u32;
            let mut c_cnt = 0u32;
            let mut total_atgc = 0u32;

            for offset in 0..=(half_win * 2) {
                // Circular index wrap: (center + len - half_win + offset) % len.
                let idx = (center + len - (half_win % len) + offset) % len;
                match ref_seq.residues[idx] {
                    b'G' => {
                        g_cnt += 1;
                        total_atgc += 1;
                    }
                    b'C' => {
                        c_cnt += 1;
                        total_atgc += 1;
                    }
                    b'A' | b'T' | b'U' => {
                        total_atgc += 1;
                    }
                    _ => {}
                }
            }

            let gc_frac = if total_atgc == 0 {
                0.5
            } else {
                (g_cnt + c_cnt) as f32 / total_atgc as f32
            };
            let gc_sum = g_cnt + c_cnt;
            let gc_skew = if gc_sum == 0 {
                0.0
            } else {
                (g_cnt as f32 - c_cnt as f32) / gc_sum as f32
            };
            profile.push((gc_frac, gc_skew));
        }

        Some(profile)
    }

    /// Recomputes consensus, Geneious-style pairwise column identity, and Shannon entropy.
    ///
    /// Uses a fixed 256-entry stack frequency table per column so recomputation performs
    /// only a single heap allocation for the output boxed slice.
    pub fn recompute_stats(&mut self) {
        let mut stats = Vec::with_capacity(self.width);
        let row_count = self.sequences.len();
        let ambiguity = if self.alphabet.is_nucleotide() {
            b'N'
        } else {
            b'X'
        };

        for col in 0..self.width {
            // Stack-allocated ASCII byte frequency histogram: zero heap allocations per column.
            let mut counts = [0u32; 256];
            let mut non_gap_count = 0u32;
            let mut gap_count = 0u32;

            for seq in &self.sequences {
                let b = seq.residues.get(col).copied().unwrap_or(b'-');
                if is_gap(b) {
                    gap_count += 1;
                } else {
                    non_gap_count += 1;
                    counts[usize::from(b)] += 1;
                }
            }

            // Find modal non-gap residue.
            let mut modal_residue = b'-';
            let mut modal_count = 0u32;
            for (byte_val, &c) in counts.iter().enumerate() {
                if c > modal_count {
                    modal_count = c;
                    modal_residue = byte_val as u8;
                }
            }

            // Consensus calling according to ConsensusConfig.
            let denom = if self.consensus_config.ignore_gaps {
                non_gap_count
            } else {
                row_count as u32
            };

            let consensus = if denom == 0 {
                b'-'
            } else if !self.consensus_config.ignore_gaps && gap_count > modal_count {
                let gap_frac = gap_count as f32 / denom as f32;
                if gap_frac >= self.consensus_config.threshold {
                    b'-'
                } else {
                    ambiguity
                }
            } else {
                let frac = modal_count as f32 / denom as f32;
                if frac >= self.consensus_config.threshold {
                    modal_residue
                } else {
                    ambiguity
                }
            };

            // Mean pairwise identity across all pairs of rows in this column:
            // Two residues match iff they are equal (and gaps are excluded from matching pairs
            // when `non_gap_count >= 2`, or counted consistently against total pairs).
            let identity_pct = if row_count <= 1 {
                if non_gap_count == 1 {
                    100.0
                } else {
                    0.0
                }
            } else {
                let total_pairs = (row_count as u64 * (row_count as u64 - 1)) / 2;
                let mut matching_pairs = 0u64;
                for &c in &counts {
                    if c >= 2 {
                        let cu = u64::from(c);
                        matching_pairs += (cu * (cu - 1)) / 2;
                    }
                }
                if total_pairs == 0 {
                    0.0
                } else {
                    ((matching_pairs as f64 / total_pairs as f64) * 100.0) as f32
                }
            };

            // Shannon entropy over non-gap residue distribution (in bits).
            let shannon_entropy = if non_gap_count <= 1 {
                0.0
            } else {
                let n = f64::from(non_gap_count);
                let mut h = 0.0f64;
                for &c in &counts {
                    if c > 0 {
                        let p = f64::from(c) / n;
                        h -= p * p.log2();
                    }
                }
                h as f32
            };

            stats.push(ColumnStats {
                identity_pct,
                shannon_entropy,
                non_gap_count,
                consensus,
                modal_residue,
            });
        }

        self.columns = stats.into_boxed_slice();
    }

    /// Computes an `N x N` row-major pairwise sequence identity percentage matrix (`0.0..=100.0`).
    ///
    /// Pairs where both sequences have a gap are excluded from the denominator, matching
    /// standard bioinformatics pairwise distance/identity tables.
    #[must_use]
    pub fn pairwise_identity_matrix(&self) -> Vec<f32> {
        let n = self.sequences.len();
        let mut matrix = vec![100.0f32; n.saturating_mul(n)];

        for i in 0..n {
            for j in (i + 1)..n {
                let s1 = &self.sequences[i].residues;
                let s2 = &self.sequences[j].residues;
                let mut compared = 0usize;
                let mut matches = 0usize;

                for (&a, &b) in s1.iter().zip(s2.iter()) {
                    if is_gap(a) && is_gap(b) {
                        continue;
                    }
                    compared += 1;
                    if a == b {
                        matches += 1;
                    }
                }

                let pct = if compared == 0 {
                    0.0
                } else {
                    (matches as f32 / compared as f32) * 100.0
                };
                matrix[i * n + j] = pct;
                matrix[j * n + i] = pct;
            }
        }
        matrix
    }

    /// Searches all sequences for a motif query, skipping alignment gaps inside each row.
    ///
    /// Supports IUPAC ambiguity wildcards when `allow_iupac` is `true` (for nucleotide
    /// alignments) or `X` wildcard for proteins.
    #[must_use]
    pub fn find_motif(&self, query: &str, allow_iupac: bool) -> Vec<MotifHit> {
        let pattern: Vec<u8> = query
            .bytes()
            .filter(|b| !b.is_ascii_whitespace() && !is_gap(*b))
            .map(|b| b.to_ascii_uppercase())
            .collect();

        if pattern.is_empty() {
            return Vec::new();
        }

        let mut hits = Vec::new();
        let is_nuc = self.alphabet.is_nucleotide();

        for (row, seq) in self.sequences.iter().enumerate() {
            // Map each non-gap position back to its alignment column so matches spanning
            // alignment gaps are still discovered accurately.
            let mut ungapped_cols: Vec<(usize, u8)> = Vec::with_capacity(seq.residues.len());
            for (col, &b) in seq.residues.iter().enumerate() {
                if !is_gap(b) {
                    ungapped_cols.push((col, b));
                }
            }

            if ungapped_cols.len() < pattern.len() {
                continue;
            }

            let last_start = ungapped_cols.len() - pattern.len();
            for start_idx in 0..=last_start {
                let mut matched = true;
                for (offset, &pat_b) in pattern.iter().enumerate() {
                    let seq_b = ungapped_cols[start_idx + offset].1;
                    if !residue_matches(pat_b, seq_b, is_nuc, allow_iupac) {
                        matched = false;
                        break;
                    }
                }
                if matched {
                    let start_col = ungapped_cols[start_idx].0;
                    let end_col = ungapped_cols[start_idx + pattern.len() - 1]
                        .0
                        .saturating_add(1);
                    hits.push(MotifHit {
                        row,
                        start_col,
                        end_col,
                    });
                }
            }
        }

        hits
    }

    /// Runs a pure-Rust center-star progressive global multiple sequence alignment
    /// (Needleman-Wunsch against the longest/reference sequence, then merges inserted gap columns).
    ///
    /// Useful when importing unaligned FASTA sequences or after adding a new sequence.
    /// Rejects sequences exceeding a safe DP matrix memory cap (`2_500` residues each) so
    /// interactive GUI execution cannot exhaust memory.
    pub fn realign_progressive(&mut self) -> Result<()> {
        const MAX_SEQ_LEN_FOR_NW: usize = 2_500;

        if self.sequences.len() <= 1 {
            return Ok(());
        }

        let ungapped: Vec<Vec<u8>> = self
            .sequences
            .iter()
            .map(Sequence::ungapped_bytes)
            .collect();

        for (idx, seq) in ungapped.iter().enumerate() {
            if seq.len() > MAX_SEQ_LEN_FOR_NW {
                return Err(HelixError::Alignment(format!(
                    "sequence {} length ({}) exceeds interactive Needleman-Wunsch limit ({MAX_SEQ_LEN_FOR_NW})",
                    self.sequences[idx].name,
                    seq.len()
                )));
            }
        }

        // Use row 0 as the reference anchor (standard Geneious "Align to Reference" / star MSA).
        let ref_seq = &ungapped[0];
        let mut pairwise_aligned: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(ungapped.len() - 1);

        // Track maximum gap run inserted before each reference position `0..=ref_seq.len()`.
        let mut max_gaps_before = vec![0usize; ref_seq.len() + 1];

        for other in ungapped.iter().skip(1) {
            let (aln_ref, aln_other) = needleman_wunsch(ref_seq, other, self.alphabet);
            let mut ref_pos = 0usize;
            let mut gap_run = 0usize;
            for &b in &aln_ref {
                if is_gap(b) {
                    gap_run += 1;
                } else {
                    if gap_run > max_gaps_before[ref_pos] {
                        max_gaps_before[ref_pos] = gap_run;
                    }
                    gap_run = 0;
                    ref_pos += 1;
                }
            }
            if gap_run > max_gaps_before[ref_pos] {
                max_gaps_before[ref_pos] = gap_run;
            }
            pairwise_aligned.push((aln_ref, aln_other));
        }

        // Compute total unified alignment width.
        let total_width = ref_seq.len() + max_gaps_before.iter().sum::<usize>();

        // Rebuild reference sequence with unified gap columns.
        let mut new_ref = Vec::with_capacity(total_width);
        for (ref_pos, &base) in ref_seq.iter().enumerate() {
            new_ref.resize(new_ref.len() + max_gaps_before[ref_pos], b'-');
            new_ref.push(base);
        }
        new_ref.resize(new_ref.len() + max_gaps_before[ref_seq.len()], b'-');
        self.sequences[0].residues = new_ref;

        // Project each non-reference pairwise alignment into the unified column layout.
        for (pair_idx, (aln_ref, aln_other)) in pairwise_aligned.into_iter().enumerate() {
            let mut projected = Vec::with_capacity(total_width);
            let mut ref_pos = 0usize;
            let mut insertion_buf: Vec<u8> = Vec::new();

            for (&r_b, &o_b) in aln_ref.iter().zip(aln_other.iter()) {
                if is_gap(r_b) {
                    insertion_buf.push(o_b);
                } else {
                    let slot = max_gaps_before[ref_pos];
                    projected.extend_from_slice(&insertion_buf);
                    projected.resize(
                        projected.len() + slot.saturating_sub(insertion_buf.len()),
                        b'-',
                    );
                    insertion_buf.clear();
                    projected.push(o_b);
                    ref_pos += 1;
                }
            }
            let trailing_slot = max_gaps_before[ref_seq.len()];
            projected.extend_from_slice(&insertion_buf);
            projected.resize(
                projected.len() + trailing_slot.saturating_sub(insertion_buf.len()),
                b'-',
            );

            self.sequences[pair_idx + 1].residues = projected;
        }

        self.width = total_width;
        // Clamp any existing annotations to the new alignment width.
        if self.width > 0 {
            for ann in &mut self.annotations {
                ann.start = ann.start.min(self.width - 1);
                ann.end = ann.end.clamp(ann.start + 1, self.width);
            }
        }
        self.recompute_stats();
        Ok(())
    }
}

/// Checks whether a query residue `pat` matches a sequence residue `target`.
fn residue_matches(pat: u8, target: u8, is_nuc: bool, allow_iupac: bool) -> bool {
    if pat == target {
        return true;
    }
    if !allow_iupac {
        return false;
    }
    if !is_nuc {
        return pat == b'X' || target == b'X';
    }
    let t = if target == b'U' { b'T' } else { target };
    match pat {
        b'N' => matches!(t, b'A' | b'C' | b'G' | b'T' | b'N'),
        b'R' => matches!(t, b'A' | b'G'),
        b'Y' => matches!(t, b'C' | b'T'),
        b'S' => matches!(t, b'G' | b'C'),
        b'W' => matches!(t, b'A' | b'T'),
        b'K' => matches!(t, b'G' | b'T'),
        b'M' => matches!(t, b'A' | b'C'),
        b'B' => matches!(t, b'C' | b'G' | b'T'),
        b'D' => matches!(t, b'A' | b'G' | b'T'),
        b'H' => matches!(t, b'A' | b'C' | b'T'),
        b'V' => matches!(t, b'A' | b'C' | b'G'),
        b'U' => t == b'T',
        _ => false,
    }
}

/// Global pairwise Needleman-Wunsch alignment using compact flat DP + traceback buffers.
fn needleman_wunsch(seq_a: &[u8], seq_b: &[u8], alphabet: Alphabet) -> (Vec<u8>, Vec<u8>) {
    let m = seq_a.len();
    let n = seq_b.len();
    let cols = n + 1;
    let gap_penalty: i32 = if alphabet.is_nucleotide() { -4 } else { -6 };

    // Use two rolling score rows for lower memory footprint, plus a compact 1-byte
    // traceback direction matrix (0 = diag, 1 = up/gap in B, 2 = left/gap in A).
    let mut prev = vec![0i32; cols];
    let mut curr = vec![0i32; cols];
    let mut trace = vec![0u8; (m + 1) * cols];

    for j in 1..=n {
        prev[j] = gap_penalty * j as i32;
        trace[j] = 2;
    }

    for i in 1..=m {
        curr[0] = gap_penalty * i as i32;
        trace[i * cols] = 1;
        let a = seq_a[i - 1];

        for j in 1..=n {
            let b = seq_b[j - 1];
            let diag = prev[j - 1] + substitution_score(a, b, alphabet);
            let up = prev[j] + gap_penalty;
            let left = curr[j - 1] + gap_penalty;

            if diag >= up && diag >= left {
                curr[j] = diag;
                trace[i * cols + j] = 0;
            } else if up >= left {
                curr[j] = up;
                trace[i * cols + j] = 1;
            } else {
                curr[j] = left;
                trace[i * cols + j] = 2;
            }
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    // Backtrack from (m, n) to (0, 0).
    let mut out_a = Vec::with_capacity(m.max(n) + 16);
    let mut out_b = Vec::with_capacity(m.max(n) + 16);
    let mut i = m;
    let mut j = n;

    while i > 0 || j > 0 {
        let dir = trace[i * cols + j];
        if i > 0 && j > 0 && dir == 0 {
            out_a.push(seq_a[i - 1]);
            out_b.push(seq_b[j - 1]);
            i -= 1;
            j -= 1;
        } else if i > 0 && (j == 0 || dir == 1) {
            out_a.push(seq_a[i - 1]);
            out_b.push(b'-');
            i -= 1;
        } else {
            out_a.push(b'-');
            out_b.push(seq_b[j - 1]);
            j -= 1;
        }
    }

    out_a.reverse();
    out_b.reverse();
    (out_a, out_b)
}

/// Substitution score for Needleman-Wunsch (biochemical similarity groups for proteins,
/// transition/transversion weighting for nucleotides).
fn substitution_score(a: u8, b: u8, alphabet: Alphabet) -> i32 {
    if a == b {
        return if alphabet.is_nucleotide() { 4 } else { 5 };
    }
    if alphabet.is_nucleotide() {
        // Transitions (A<->G, C<->T/U) penalized less harshly than transversions.
        let is_transition = matches!(
            (a, b),
            (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') | (b'C', b'U') | (b'U', b'C')
        );
        return if is_transition { -1 } else { -3 };
    }

    // Physicochemical conservative amino acid substitution groups (BLOSUM-inspired).
    let group = |aa: u8| -> u8 {
        match aa {
            b'I' | b'L' | b'V' | b'M' => 1, // Aliphatic hydrophobic
            b'F' | b'Y' | b'W' => 2,        // Aromatic
            b'K' | b'R' | b'H' => 3,        // Positively charged
            b'D' | b'E' => 4,               // Negatively charged
            b'S' | b'T' | b'N' | b'Q' => 5, // Polar uncharged
            b'A' | b'G' => 6,               // Small
            _ => 0,
        }
    };
    let ga = group(a);
    if ga != 0 && ga == group(b) {
        2
    } else {
        -2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_computes_consensus_and_pairwise_identity() {
        let seqs = vec![
            Sequence::new("S1", b"ACGT"),
            Sequence::new("S2", b"ACGA"),
            Sequence::new("S3", b"ACGT"),
        ];
        let aln = Alignment::new("Test", Alphabet::Dna, seqs).expect("valid alignment");
        assert_eq!(aln.width(), 4);
        assert_eq!(aln.columns()[0].consensus, b'A');
        assert!((aln.columns()[0].identity_pct - 100.0).abs() < 0.01);
        // Column 3 has T, A, T -> 1 matching pair out of 3 = 33.33%
        assert!((aln.columns()[3].identity_pct - 33.333).abs() < 0.1);
        assert_eq!(aln.columns()[3].consensus, b'T');
    }

    #[test]
    fn motif_search_skips_gaps_and_supports_iupac() {
        let seqs = vec![
            Sequence::new("S1", b"AT--GCA"),
            Sequence::new("S2", b"ATGTA"),
        ];
        let aln = Alignment::new("MotifTest", Alphabet::Dna, seqs).expect("valid");
        let hits = aln.find_motif("ATGYA", true);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0], MotifHit { row: 0, start_col: 0, end_col: 7 });
        assert_eq!(hits[1], MotifHit { row: 1, start_col: 0, end_col: 5 });
    }

    #[test]
    fn progressive_realign_aligns_unaligned_sequences() {
        let seqs = vec![
            Sequence::new("Ref", b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFF"),
            Sequence::new("Var", b"MVHLTPVEKSAVTALWGKVNDEVGGEALGRLLVVYPWTQRFF"),
        ];
        let mut aln = Alignment::new("HbB", Alphabet::Protein, seqs).expect("valid");
        aln.realign_progressive().expect("realign succeeds");
        assert_eq!(aln.sequences()[0].len(), aln.sequences()[1].len());
        assert!(aln.mean_pairwise_identity_pct() > 90.0);
    }

    #[test]
    fn circular_gc_profile_wraps_cleanly_for_dna_and_returns_none_for_protein() {
        let dna = Alignment::new(
            "PhageMini",
            Alphabet::Dna,
            vec![Sequence::new("P1", b"GGGGGGGGCCCCAAAATTTTGGGGCCCC")],
        )
        .expect("dna");
        let prof = dna.circular_gc_profile(24).expect("gc profile");
        assert_eq!(prof.len(), 24);
        for (gc, skew) in prof {
            assert!((0.0..=1.0).contains(&gc));
            assert!((-1.0..=1.0).contains(&skew));
        }

        let prot = Alignment::new(
            "Prot",
            Alphabet::Protein,
            vec![Sequence::new("S1", b"MVHLTPEEK")],
        )
        .expect("prot");
        assert!(prot.circular_gc_profile(24).is_none());
    }
}
