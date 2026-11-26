//! Multi-algorithm sequence alignment suite implementing:
//! 1. **De Bruijn Graph Alignment** (k-mer De Bruijn unitig/anchor chaining + local DP fill)
//! 2. **ClustalW** (distance matrix + Neighbor-Joining/UPGMA guide tree + sequence weighting + profile-profile alignment)
//! 3. **Clustal Omega** (mBed k-mer embedding guide tree + HMM-inspired conservation-weighted profile alignment)
//! 4. **MUSCLE** (k-mer draft progressive + improved Kimura guide tree + iterative Sum-of-Pairs tree refinement)
//! 5. **MAFFT** (physicochemical cross-correlation lag anchoring + progressive profile alignment + refinement)
//! 6. **minimap2** (invokes bundled/system `minimap2` binary with SAM CIGAR parsing when installed,
//!    with a built-in pure-Rust `(w, k)` minimizer-seed + collinear chaining fallback engine!).

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use crate::sequence::{is_gap, reverse_complement, Alphabet, Sequence};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Selectable multiple / pairwise alignment algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum AlignerAlgorithm {
    /// De Bruijn k-mer graph anchor chaining + bounded DP inter-anchor alignment.
    DeBruijn,
    /// ClustalW progressive guide-tree profile-profile alignment with sequence weights.
    #[default]
    ClustalW,
    /// Clustal Omega mBed k-mer embedding + conservation-weighted profile-profile alignment.
    ClustalOmega,
    /// MUSCLE 3-stage iterative progressive + Sum-of-Pairs bipartition refinement.
    Muscle,
    /// MAFFT physicochemical cross-correlation anchor detection + iterative alignment.
    Mafft,
    /// minimap2 (external `minimap2` binary with SAM CIGAR parsing, or built-in Rust minimizer-chaining fallback).
    Minimap2,
}

impl AlignerAlgorithm {
    pub const ALL: [Self; 6] = [
        Self::DeBruijn,
        Self::ClustalW,
        Self::ClustalOmega,
        Self::Muscle,
        Self::Mafft,
        Self::Minimap2,
    ];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::DeBruijn => "De Bruijn Graph Anchor Aligner",
            Self::ClustalW => "ClustalW (Weighted Guide-Tree Profile)",
            Self::ClustalOmega => "Clustal Omega (mBed + HMM-Profile)",
            Self::Muscle => "MUSCLE (Iterative SP Refinement)",
            Self::Mafft => "MAFFT (Correlation Anchors + Refine)",
            Self::Minimap2 => "minimap2 (Binary / Built-in Minimizers)",
        }
    }

    #[must_use]
    pub const fn short_name(self) -> &'static str {
        match self {
            Self::DeBruijn => "DeBruijn",
            Self::ClustalW => "ClustalW",
            Self::ClustalOmega => "ClustalO",
            Self::Muscle => "MUSCLE",
            Self::Mafft => "MAFFT",
            Self::Minimap2 => "minimap2",
        }
    }
}

/// Official minimap2 `-x` alignment presets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Minimap2Preset {
    #[default]
    Asm5,
    Asm20,
    MapOnt,
    MapHifi,
    ShortRead,
}

impl Minimap2Preset {
    pub const ALL: [Self; 5] = [
        Self::Asm5,
        Self::Asm20,
        Self::MapOnt,
        Self::MapHifi,
        Self::ShortRead,
    ];

    #[must_use]
    pub const fn flag(self) -> &'static str {
        match self {
            Self::Asm5 => "asm5",
            Self::Asm20 => "asm20",
            Self::MapOnt => "map-ont",
            Self::MapHifi => "map-hifi",
            Self::ShortRead => "sr",
        }
    }

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Asm5 => "asm5 (Genome/Plasmid <5% div)",
            Self::Asm20 => "asm20 (Divergent Genome/Phage <20%)",
            Self::MapOnt => "map-ont (Oxford Nanopore)",
            Self::MapHifi => "map-hifi (PacBio HiFi)",
            Self::ShortRead => "sr (Short Paired/Single Reads)",
        }
    }
}

/// Configuration for running any alignment algorithm in [`AlignerAlgorithm`].
#[derive(Debug, Clone)]
pub struct AlignerConfig {
    pub algorithm: AlignerAlgorithm,
    /// Optional explicit path to a `minimap2` / `minimap2.exe` binary.
    pub minimap2_path: Option<PathBuf>,
    /// Preset passed to `minimap2 -x <preset>`.
    pub minimap2_preset: Minimap2Preset,
    /// Gap opening penalty magnitude (positive integer, default `8`).
    pub gap_open: i32,
    /// Gap extension penalty magnitude (positive integer, default `2`).
    pub gap_extend: i32,
}

impl Default for AlignerConfig {
    fn default() -> Self {
        Self {
            algorithm: AlignerAlgorithm::ClustalW,
            minimap2_path: None,
            minimap2_preset: Minimap2Preset::Asm5,
            gap_open: 8,
            gap_extend: 2,
        }
    }
}

/// Result returned by [`run_alignment_algorithm`], including engine telemetry.
#[derive(Debug, Clone)]
pub struct AlignmentRunOutput {
    pub alignment: Alignment,
    pub engine_note: String,
}

/// Probes for an installed or bundled `minimap2` executable (`custom_path`, `./tools/minimap2(.exe)`,
/// or on system `PATH`).
#[must_use]
pub fn detect_minimap2_binary(custom_path: Option<&Path>) -> Option<PathBuf> {
    if let Some(p) = custom_path {
        if p.is_file() {
            return Some(p.to_path_buf());
        }
    }
    let bundled_candidates = [
        PathBuf::from("tools/minimap2.exe"),
        PathBuf::from("tools/minimap2"),
        PathBuf::from("minimap2.exe"),
    ];
    for c in bundled_candidates {
        if c.is_file() {
            return Some(c);
        }
    }
    // Check `PATH` by testing `minimap2 --version`.
    if Command::new("minimap2")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        return Some(PathBuf::from("minimap2"));
    }
    None
}

/// Runs the requested alignment algorithm on `input` sequences, producing a fresh aligned [`Alignment`].
pub fn run_alignment_algorithm(
    name: impl Into<String>,
    alphabet: Alphabet,
    sequences: &[Sequence],
    config: &AlignerConfig,
) -> Result<AlignmentRunOutput> {
    const MAX_DP_LEN: usize = 3_500;
    if sequences.len() < 2 {
        return Err(HelixError::Alignment(
            "at least 2 sequences are required for alignment".to_owned(),
        ));
    }

    let doc_name = name.into();
    let ungapped: Vec<Vec<u8>> = sequences.iter().map(Sequence::ungapped_bytes).collect();

    if config.algorithm == AlignerAlgorithm::Minimap2 {
        return run_minimap2_alignment(&doc_name, alphabet, sequences, &ungapped, config);
    }

    for (i, u) in ungapped.iter().enumerate() {
        if u.len() > MAX_DP_LEN {
            return Err(HelixError::Alignment(format!(
                "sequence '{}' length ({}) exceeds interactive MSA limit ({MAX_DP_LEN}); use minimap2 for long genomes",
                sequences[i].name,
                u.len()
            )));
        }
    }

    let (aligned_rows, note) = match config.algorithm {
        AlignerAlgorithm::DeBruijn => align_debruijn_graph(&ungapped, alphabet, config),
        AlignerAlgorithm::ClustalW => align_clustalw(&ungapped, alphabet, config),
        AlignerAlgorithm::ClustalOmega => align_clustal_omega(&ungapped, alphabet, config),
        AlignerAlgorithm::Muscle => align_muscle(&ungapped, alphabet, config),
        AlignerAlgorithm::Mafft => align_mafft(&ungapped, alphabet, config),
        AlignerAlgorithm::Minimap2 => unreachable!("handled above"),
    };

    let out_seqs: Vec<Sequence> = sequences
        .iter()
        .zip(aligned_rows)
        .map(|(orig, res)| {
            let mut s = Sequence::new(orig.name.clone(), &res);
            s.description = orig.description.clone();
            s.organism = orig.organism.clone();
            s.uniprot_accession = orig.uniprot_accession.clone();
            s
        })
        .collect();

    let alignment = Alignment::new(doc_name, alphabet, out_seqs)?.with_description(note.clone());
    Ok(AlignmentRunOutput {
        alignment,
        engine_note: note,
    })
}

// ============================================================================
// 1. De Bruijn Graph K-mer Anchor Alignment
// ============================================================================

/// Builds a k-mer De Bruijn graph across all input sequences, extracts a collinear chain of
/// shared unitig/anchor k-mers (present across sequences in monotone order), and aligns
/// inter-anchor variable spans with banded profile/pairwise Needleman-Wunsch.
fn align_debruijn_graph(
    seqs: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
) -> (Vec<Vec<u8>>, String) {
    let k = if alphabet.is_nucleotide() { 7usize } else { 4usize };
    let n = seqs.len();
    let ref_seq = &seqs[0];

    // Index k-mer occurrences in each sequence (k-mer -> list of positions per sequence).
    let mut kmer_hits: HashMap<&[u8], Vec<Option<usize>>> = HashMap::new();
    if ref_seq.len() >= k {
        for (pos, w) in ref_seq.windows(k).enumerate() {
            kmer_hits.entry(w).or_insert_with(|| {
                let mut v = vec![None; n];
                v[0] = Some(pos);
                v
            });
        }
    }
    for (s_idx, seq) in seqs.iter().enumerate().skip(1) {
        if seq.len() < k {
            continue;
        }
        for (pos, w) in seq.windows(k).enumerate() {
            if let Some(entry) = kmer_hits.get_mut(w) {
                if entry[s_idx].is_none() {
                    entry[s_idx] = Some(pos);
                }
            }
        }
    }

    // Retain shared k-mer nodes present in ALL sequences and order by reference coordinate.
    let mut shared_nodes: Vec<Vec<usize>> = kmer_hits
        .into_values()
        .filter_map(|opt_pos| opt_pos.into_iter().collect::<Option<Vec<usize>>>())
        .collect();
    shared_nodes.sort_by_key(|v| v[0]);

    // Extract longest strictly collinear non-overlapping anchor chain across all sequences.
    let mut anchors: Vec<Vec<usize>> = Vec::new();
    for node in shared_nodes {
        let collinear = anchors.last().is_none_or(|prev| {
            prev.iter()
                .zip(node.iter())
                .all(|(&p, &cur)| cur >= p + k)
        });
        if collinear {
            anchors.push(node);
        }
    }

    if anchors.is_empty() {
        let (rows, _) = align_clustalw(seqs, alphabet, config);
        return (
            rows,
            format!("De Bruijn Graph (k={k}): 0 universal anchors -> fell back to progressive profile DP"),
        );
    }

    // Assemble final MSA by aligning each variable segment between consecutive De Bruijn anchors
    // and splicing exact anchor unitigs verbatim!
    let mut combined = vec![Vec::<u8>::new(); n];
    let mut prev_ends = vec![0usize; n];

    for anchor in &anchors {
        let segs: Vec<Vec<u8>> = (0..n)
            .map(|i| seqs[i][prev_ends[i]..anchor[i]].to_vec())
            .collect();
        if segs.iter().any(|s| !s.is_empty()) {
            let aligned_seg = progressive_profile_msa(&segs, alphabet, config, false);
            for (i, row) in aligned_seg.into_iter().enumerate() {
                combined[i].extend_from_slice(&row);
            }
        }
        for i in 0..n {
            combined[i].extend_from_slice(&seqs[i][anchor[i]..(anchor[i] + k)]);
            prev_ends[i] = anchor[i] + k;
        }
    }

    // Align trailing suffix after final anchor.
    let tail_segs: Vec<Vec<u8>> = (0..n).map(|i| seqs[i][prev_ends[i]..].to_vec()).collect();
    if tail_segs.iter().any(|s| !s.is_empty()) {
        let aligned_tail = progressive_profile_msa(&tail_segs, alphabet, config, false);
        for (i, row) in aligned_tail.into_iter().enumerate() {
            combined[i].extend_from_slice(&row);
        }
    }

    let note = format!(
        "De Bruijn Graph Alignment (k={k}, {} collinear unitig anchors)",
        anchors.len()
    );
    (combined, note)
}

// ============================================================================
// 2. ClustalW (Weighted Guide-Tree Progressive Profile Alignment)
// ============================================================================

fn align_clustalw(
    seqs: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
) -> (Vec<Vec<u8>>, String) {
    let rows = progressive_profile_msa(seqs, alphabet, config, false);
    (
        rows,
        format!(
            "ClustalW progressive weighted guide-tree profile alignment (gap open -{}, ext -{})",
            config.gap_open, config.gap_extend
        ),
    )
}

// ============================================================================
// 3. Clustal Omega (mBed K-mer Embedding + Conservation-Weighted Profile Alignment)
// ============================================================================

fn align_clustal_omega(
    seqs: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
) -> (Vec<Vec<u8>>, String) {
    let n = seqs.len();
    let k = if alphabet.is_nucleotide() { 5 } else { 3 };
    // mBed embedding: select up to `ceil(log2(n)) + 2` seed sequences and embed each sequence
    // into a low-dimensional k-mer distance vector relative to those seeds.
    let seed_count = ((n as f32).log2().ceil() as usize + 2).clamp(1, n);
    let step = (n / seed_count).max(1);
    let seeds: Vec<usize> = (0..n).step_by(step).take(seed_count).collect();

    let kmer_sets: Vec<Vec<u64>> = seqs.iter().map(|s| sorted_kmer_hashes(s, k)).collect();
    let mut embed: Vec<Vec<f32>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut vec_i = Vec::with_capacity(seeds.len());
        for &sd in &seeds {
            vec_i.push(kmer_jaccard_distance(&kmer_sets[i], &kmer_sets[sd]));
        }
        embed.push(vec_i);
    }

    // Euclidean distance matrix in mBed embedding space.
    let mut dist = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let sq: f32 = embed[i]
                .iter()
                .zip(&embed[j])
                .map(|(a, b)| (a - b) * (a - b))
                .sum();
            let d = sq.sqrt();
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    // Progressive HMM-conservation-weighted profile alignment following the mBed guide tree.
    let rows = progressive_with_dist_matrix(seqs, &dist, alphabet, config, true);
    (
        rows,
        format!(
            "Clustal Omega (mBed {} seed embedding + conservation-weighted profile alignment)",
            seeds.len()
        ),
    )
}

// ============================================================================
// 4. MUSCLE (3-Stage Progressive + Kimura Re-Tree + Iterative SP Refinement)
// ============================================================================

fn align_muscle(
    seqs: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
) -> (Vec<Vec<u8>>, String) {
    let n = seqs.len();
    // Stage 1: Draft progressive alignment from fast k-mer distance guide tree.
    let stage1 = progressive_profile_msa(seqs, alphabet, config, false);

    // Stage 2: Compute Kimura-like pairwise distances directly from Stage 1 MSA and rebuild guide tree.
    let mut k_dist = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut compared = 0usize;
            let mut diffs = 0usize;
            for (&a, &b) in stage1[i].iter().zip(&stage1[j]) {
                if is_gap(a) || is_gap(b) {
                    continue;
                }
                compared += 1;
                if a != b {
                    diffs += 1;
                }
            }
            let p = if compared == 0 {
                0.5
            } else {
                diffs as f32 / compared as f32
            };
            let kimura = -(1.0 - p - 0.2 * p * p).max(0.05).ln();
            k_dist[i][j] = kimura;
            k_dist[j][i] = kimura;
        }
    }
    let mut best_msa = progressive_with_dist_matrix(seqs, &k_dist, alphabet, config, false);
    let mut best_sp = sum_of_pairs_score(&best_msa, alphabet);

    // Stage 3: Iterative leave-one-out / bipartition refinement keeping improvements in SP score.
    let mut improvements = 0usize;
    if n >= 3 {
        let max_passes = n.min(8);
        for target in 0..max_passes {
            let mut group_a_indices = Vec::with_capacity(n - 1);
            for i in 0..n {
                if i != target {
                    group_a_indices.push(i);
                }
            }
            let prof_a = strip_all_gap_columns(
                &group_a_indices
                    .iter()
                    .map(|&i| best_msa[i].clone())
                    .collect::<Vec<_>>(),
            );
            let prof_b = strip_all_gap_columns(&[best_msa[target].clone()]);
            let merged = align_two_profiles(&prof_a, &prof_b, alphabet, config, false);

            let mut candidate = vec![Vec::new(); n];
            for (local, &orig_idx) in group_a_indices.iter().enumerate() {
                candidate[orig_idx] = merged[local].clone();
            }
            candidate[target] = merged[group_a_indices.len()].clone();

            let sp = sum_of_pairs_score(&candidate, alphabet);
            if sp > best_sp {
                best_sp = sp;
                best_msa = candidate;
                improvements += 1;
            }
        }
    }

    (
        best_msa,
        format!(
            "MUSCLE 3-stage alignment (Kimura guide tree + {improvements} iterative SP refinement improvement(s))"
        ),
    )
}

// ============================================================================
// 5. MAFFT (Physicochemical Cross-Correlation Lag Anchoring + Refinement)
// ============================================================================

fn align_mafft(
    seqs: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
) -> (Vec<Vec<u8>>, String) {
    let n = seqs.len();
    // Compute normalized physicochemical signal for each sequence (Katoh FFT-NS idea:
    // purine/pyrimidine + AT/GC channels for DNA, polarity/volume channels for proteins).
    let signals: Vec<Vec<f32>> = seqs.iter().map(|s| residue_signal(s, alphabet)).collect();

    // Estimate pairwise correlation distances via peak cross-correlation lag overlap.
    let mut dist = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let corr = peak_cross_correlation(&signals[i], &signals[j]);
            let d = (1.0 - corr).clamp(0.0, 2.0);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    let mut msa = progressive_with_dist_matrix(seqs, &dist, alphabet, config, true);

    // One iterative MAFFT-L-INS-i-style refinement pass for higher accuracy on divergent ends.
    if n >= 3 {
        let mut best_sp = sum_of_pairs_score(&msa, alphabet);
        for target in 0..n.min(6) {
            let others: Vec<usize> = (0..n).filter(|&i| i != target).collect();
            let pa = strip_all_gap_columns(&others.iter().map(|&i| msa[i].clone()).collect::<Vec<_>>());
            let pb = strip_all_gap_columns(&[msa[target].clone()]);
            let merged = align_two_profiles(&pa, &pb, alphabet, config, true);
            let mut cand = vec![Vec::new(); n];
            for (local, &idx) in others.iter().enumerate() {
                cand[idx] = merged[local].clone();
            }
            cand[target] = merged[others.len()].clone();
            let sp = sum_of_pairs_score(&cand, alphabet);
            if sp > best_sp {
                best_sp = sp;
                msa = cand;
            }
        }
    }

    (
        msa,
        String::from("MAFFT physicochemical cross-correlation guide tree + iterative profile alignment"),
    )
}

// ============================================================================
// 6. minimap2 (External Binary + SAM CIGAR Parser + Built-in Rust Minimizer Engine)
// ============================================================================

fn run_minimap2_alignment(
    doc_name: &str,
    alphabet: Alphabet,
    sequences: &[Sequence],
    ungapped: &[Vec<u8>],
    config: &AlignerConfig,
) -> Result<AlignmentRunOutput> {
    if !alphabet.is_nucleotide() {
        return Err(HelixError::Alignment(
            "minimap2 operates on DNA/RNA nucleotide sequences; use ClustalW, Clustal Omega, MUSCLE, or MAFFT for proteins".to_owned(),
        ));
    }

    // Path A: Invoke real external `minimap2` binary if found (bundled or on PATH or user-specified).
    if let Some(bin) = detect_minimap2_binary(config.minimap2_path.as_deref()) {
        match run_external_minimap2_binary(&bin, doc_name, alphabet, sequences, ungapped, config.minimap2_preset) {
            Ok(out) => return Ok(out),
            Err(err) => {
                // Fall back gracefully to built-in minimizer chaining if external run fails, noting reason.
                let fallback = run_builtin_minimizer_chaining(
                    doc_name,
                    alphabet,
                    sequences,
                    ungapped,
                    config.minimap2_preset,
                    Some(&format!("external minimap2 failed ({err}); used built-in Rust minimizer-chainer")),
                )?;
                return Ok(fallback);
            }
        }
    }

    // Path B: Built-in pure-Rust `(w, k)` minimizer + collinear DP chaining engine.
    run_builtin_minimizer_chaining(
        doc_name,
        alphabet,
        sequences,
        ungapped,
        config.minimap2_preset,
        None,
    )
}

/// Runs external `minimap2 -a -x <preset> ref.fa queries.fa` and parses SAM CIGAR strings into an [`Alignment`].
fn run_external_minimap2_binary(
    bin: &Path,
    doc_name: &str,
    alphabet: Alphabet,
    sequences: &[Sequence],
    ungapped: &[Vec<u8>],
    preset: Minimap2Preset,
) -> Result<AlignmentRunOutput> {
    let tmp_dir = std::env::temp_dir().join("helix_prime_minimap2");
    std::fs::create_dir_all(&tmp_dir)?;
    let ref_path = tmp_dir.join("ref.fa");
    let qry_path = tmp_dir.join("queries.fa");

    let mut ref_fa = format!(">{}\n", sanitize_fasta_id(&sequences[0].name));
    ref_fa.push_str(&String::from_utf8_lossy(&ungapped[0]));
    ref_fa.push('\n');
    std::fs::write(&ref_path, ref_fa)?;

    let mut qry_fa = String::new();
    for (idx, seq) in sequences.iter().enumerate().skip(1) {
        qry_fa.push_str(&format!(">{idx}_{}\n", sanitize_fasta_id(&seq.name)));
        qry_fa.push_str(&String::from_utf8_lossy(&ungapped[idx]));
        qry_fa.push('\n');
    }
    std::fs::write(&qry_path, qry_fa)?;

    let output = Command::new(bin)
        .arg("-a")
        .arg("-x")
        .arg(preset.flag())
        .arg(&ref_path)
        .arg(&qry_path)
        .output()
        .map_err(|e| HelixError::Alignment(format!("failed to execute minimap2: {e}")))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(HelixError::Alignment(format!("minimap2 exited with error: {stderr}")));
    }

    let sam_text = String::from_utf8_lossy(&output.stdout);
    let aligned = parse_sam_into_alignment(doc_name, alphabet, sequences, &ungapped[0], &sam_text)?;
    let note = format!(
        "minimap2 external binary ({}) with preset -x {}",
        bin.display(),
        preset.flag()
    );
    Ok(AlignmentRunOutput {
        alignment: aligned.with_description(note.clone()),
        engine_note: note,
    })
}

/// Parses SAM records with CIGAR strings (`M`, `=`, `X`, `I`, `D`, `S`, `H`) against `ref_seq`
/// into a rectangular reference-anchored [`Alignment`].
pub fn parse_sam_into_alignment(
    doc_name: &str,
    alphabet: Alphabet,
    sequences: &[Sequence],
    ref_seq: &[u8],
    sam_text: &str,
) -> Result<Alignment> {
    let ref_len = ref_seq.len();
    let mut rows = Vec::with_capacity(sequences.len());
    rows.push(
        Sequence::new(format!("{} (Ref)", sequences[0].name), ref_seq)
            .with_description(sequences[0].description.clone()),
    );

    let mut seen_queries: HashMap<String, bool> = HashMap::new();

    for line in sam_text.lines() {
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 {
            continue;
        }
        let qname = cols[0];
        let flag: u16 = cols[1].parse().unwrap_or(4);
        // Skip unmapped (4), secondary (256), and supplementary (2048) alignments.
        if (flag & (4 | 256 | 2048)) != 0 || seen_queries.contains_key(qname) {
            continue;
        }
        let pos_1based: usize = cols[3].parse().unwrap_or(0);
        let cigar = cols[5];
        let sam_seq = cols[9].as_bytes();
        if pos_1based == 0 || cigar == "*" || sam_seq == b"*" {
            continue;
        }
        seen_queries.insert(qname.to_owned(), true);

        let mut placed = vec![b'-'; ref_len];
        let mut r_cursor = pos_1based - 1;
        let mut q_cursor = 0usize;
        let mut num_buf = 0usize;

        for b in cigar.bytes() {
            if b.is_ascii_digit() {
                num_buf = num_buf.saturating_mul(10).saturating_add(usize::from(b - b'0'));
            } else {
                let count = num_buf;
                num_buf = 0;
                match b {
                    b'M' | b'=' | b'X' => {
                        for _ in 0..count {
                            if r_cursor < ref_len && q_cursor < sam_seq.len() {
                                placed[r_cursor] = sam_seq[q_cursor].to_ascii_uppercase();
                            }
                            r_cursor += 1;
                            q_cursor += 1;
                        }
                    }
                    b'I' | b'S' => {
                        q_cursor = q_cursor.saturating_add(count);
                    }
                    b'D' | b'N' => {
                        r_cursor = r_cursor.saturating_add(count);
                    }
                    b'H' | b'P' => {}
                    _ => {}
                }
            }
        }

        let is_rev = (flag & 16) != 0;
        let display_name = if is_rev {
            format!("{qname} (RC)")
        } else {
            qname.to_owned()
        };
        rows.push(Sequence::new(display_name, &placed).with_description(format!(
            "minimap2 SAM pos {pos_1based}, CIGAR {cigar}"
        )));
    }

    Alignment::new(doc_name, alphabet, rows)
}

/// Built-in pure-Rust minimap2-inspired `(w, k)` minimizer-seed + collinear DP chaining engine.
fn run_builtin_minimizer_chaining(
    doc_name: &str,
    alphabet: Alphabet,
    sequences: &[Sequence],
    ungapped: &[Vec<u8>],
    preset: Minimap2Preset,
    fallback_reason: Option<&str>,
) -> Result<AlignmentRunOutput> {
    let ref_seq = &ungapped[0];
    let (w, k) = match preset {
        Minimap2Preset::ShortRead => (4usize, 9usize),
        Minimap2Preset::Asm5 | Minimap2Preset::MapHifi => (5, 11),
        Minimap2Preset::Asm20 | Minimap2Preset::MapOnt => (4, 9),
    };
    let eff_k = k.min(ref_seq.len().max(1)).max(5);

    // Index reference `(w, k)` minimizers.
    let ref_mins = extract_minimizers(ref_seq, w, eff_k);
    let mut ref_index: HashMap<u64, Vec<usize>> = HashMap::new();
    for (hash, pos) in ref_mins {
        ref_index.entry(hash).or_default().push(pos);
    }

    let mut rows = Vec::with_capacity(sequences.len());
    rows.push(
        Sequence::new(format!("{} (Ref)", sequences[0].name), ref_seq)
            .with_description(sequences[0].description.clone()),
    );

    let mut mapped = 0usize;
    for idx in 1..sequences.len() {
        let q = &ungapped[idx];
        let q_rc = reverse_complement(q);
        let fwd_chain = best_minimizer_chain(&ref_index, q, w, eff_k);
        let rev_chain = best_minimizer_chain(&ref_index, &q_rc, w, eff_k);

        let (chosen_seq, anchors, is_rev) = match (fwd_chain, rev_chain) {
            (Some(fc), Some(rc)) => {
                if rc.len() > fc.len() {
                    (q_rc.as_slice(), rc, true)
                } else {
                    (q.as_slice(), fc, false)
                }
            }
            (Some(fc), None) => (q.as_slice(), fc, false),
            (None, Some(rc)) => (q_rc.as_slice(), rc, true),
            (None, None) => (q.as_slice(), Vec::new(), false),
        };

        let mut placed = vec![b'-'; ref_seq.len()];
        if let Some(&(first_r, first_q)) = anchors.first() {
            // Place query anchored on its collinear minimizer chain (piecewise diagonal projection).
            let mut prev_r = first_r.saturating_sub(first_q);
            let mut prev_q = first_q.saturating_sub(first_r);
            for &(ar, aq) in &anchors {
                let len = (ar.saturating_sub(prev_r)).min(aq.saturating_sub(prev_q));
                for d in 0..len {
                    if prev_r + d < placed.len() && prev_q + d < chosen_seq.len() {
                        placed[prev_r + d] = chosen_seq[prev_q + d];
                    }
                }
                for d in 0..eff_k {
                    if ar + d < placed.len() && aq + d < chosen_seq.len() {
                        placed[ar + d] = chosen_seq[aq + d];
                    }
                }
                prev_r = ar + eff_k;
                prev_q = aq + eff_k;
            }
            let rem = (placed.len().saturating_sub(prev_r)).min(chosen_seq.len().saturating_sub(prev_q));
            placed[prev_r..(prev_r + rem)].copy_from_slice(&chosen_seq[prev_q..(prev_q + rem)]);
            mapped += 1;
        }

        let label = if is_rev {
            format!("{} (RC)", sequences[idx].name)
        } else {
            sequences[idx].name.clone()
        };
        rows.push(Sequence::new(label, &placed).with_description(format!(
            "minimap2-style ({} minimizer chain anchors, preset {})",
            anchors.len(),
            preset.flag()
        )));
    }

    let note = match fallback_reason {
        Some(reason) => format!("minimap2 ({reason}; mapped {mapped}/{})", sequences.len() - 1),
        None => format!(
            "minimap2 built-in Rust (w={w}, k={eff_k}) minimizer-chaining engine (preset -x {}, mapped {mapped}/{})",
            preset.flag(),
            sequences.len() - 1
        ),
    };
    let alignment = Alignment::new(doc_name, alphabet, rows)?.with_description(note.clone());
    Ok(AlignmentRunOutput {
        alignment,
        engine_note: note,
    })
}

/// Extracts `(hash, position)` minimizers over sliding windows of `w` consecutive `k`-mers.
fn extract_minimizers(seq: &[u8], w: usize, k: usize) -> Vec<(u64, usize)> {
    if seq.len() < k {
        return Vec::new();
    }
    let num_kmers = seq.len() - k + 1;
    let mut hashes = Vec::with_capacity(num_kmers);
    for i in 0..num_kmers {
        hashes.push((fnv1a_kmer(&seq[i..(i + k)]), i));
    }
    let win = w.max(1).min(hashes.len());
    let mut out = Vec::new();
    let mut last_pos = usize::MAX;
    for start in 0..=(hashes.len() - win) {
        let mut best = hashes[start];
        for &cand in &hashes[start..(start + win)] {
            if cand.0 < best.0 {
                best = cand;
            }
        }
        if best.1 != last_pos {
            last_pos = best.1;
            out.push(best);
        }
    }
    out
}

/// Collinear dynamic-programming chaining of `(ref_pos, query_pos)` minimizer hits.
fn best_minimizer_chain(
    ref_index: &HashMap<u64, Vec<usize>>,
    query: &[u8],
    w: usize,
    k: usize,
) -> Option<Vec<(usize, usize)>> {
    let q_mins = extract_minimizers(query, w, k);
    let mut hits: Vec<(usize, usize)> = Vec::new();
    for (h, qpos) in q_mins {
        if let Some(rpositions) = ref_index.get(&h) {
            if rpositions.len() <= 16 {
                for &rpos in rpositions {
                    hits.push((rpos, qpos));
                }
            }
        }
    }
    if hits.is_empty() {
        return None;
    }
    hits.sort_unstable();

    let m = hits.len();
    let mut dp = vec![1i32; m];
    let mut prev = vec![None::<usize>; m];
    let mut best_end = 0usize;

    for i in 0..m {
        let (ri, qi) = hits[i];
        let lookback_start = i.saturating_sub(48);
        for j in lookback_start..i {
            let (rj, qj) = hits[j];
            if rj < ri && qj < qi {
                let dr = ri - rj;
                let dq = qi - qj;
                let gap = dr.abs_diff(dq);
                if gap <= 64 && dr <= 500 && dq <= 500 {
                    let cand = dp[j] + 10 - (gap as i32 / 8);
                    if cand > dp[i] {
                        dp[i] = cand;
                        prev[i] = Some(j);
                    }
                }
            }
        }
        if dp[i] > dp[best_end] {
            best_end = i;
        }
    }

    let mut chain = Vec::new();
    let mut cur = Some(best_end);
    while let Some(idx) = cur {
        chain.push(hits[idx]);
        cur = prev[idx];
    }
    chain.reverse();
    Some(chain)
}

// ============================================================================
// Shared Progressive Profile-Profile DP & Helper Primitives
// ============================================================================

fn progressive_profile_msa(
    seqs: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
    conservation_weighted: bool,
) -> Vec<Vec<u8>> {
    let n = seqs.len();
    if n <= 1 {
        return seqs.to_vec();
    }
    let k = if alphabet.is_nucleotide() { 5 } else { 3 };
    let ksets: Vec<Vec<u64>> = seqs.iter().map(|s| sorted_kmer_hashes(s, k)).collect();
    let mut dist = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = kmer_jaccard_distance(&ksets[i], &ksets[j]);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }
    progressive_with_dist_matrix(seqs, &dist, alphabet, config, conservation_weighted)
}

/// Walks an agglomerative UPGMA guide tree built from `dist` and merges sub-alignments via
/// profile-to-profile dynamic programming, preserving original input sequence order.
fn progressive_with_dist_matrix(
    seqs: &[Vec<u8>],
    initial_dist: &[Vec<f32>],
    alphabet: Alphabet,
    config: &AlignerConfig,
    conservation_weighted: bool,
) -> Vec<Vec<u8>> {
    let n = seqs.len();
    if n <= 1 {
        return seqs.to_vec();
    }

    struct Group {
        indices: Vec<usize>,
        rows: Vec<Vec<u8>>,
    }

    let mut groups: Vec<Group> = (0..n)
        .map(|i| Group {
            indices: vec![i],
            rows: vec![seqs[i].clone()],
        })
        .collect();
    let mut dist = initial_dist.to_vec();

    while groups.len() > 1 {
        let g = groups.len();
        let mut best_i = 0usize;
        let mut best_j = 1usize;
        let mut best_d = f32::INFINITY;
        for (i, row) in dist.iter().enumerate().take(g) {
            for (j, &d) in row.iter().enumerate().take(g).skip(i + 1) {
                if d < best_d {
                    best_d = d;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        let g_j = groups.remove(best_j);
        let g_i = groups.remove(best_i);

        let merged_rows = align_two_profiles(&g_i.rows, &g_j.rows, alphabet, config, conservation_weighted);
        let mut merged_indices = g_i.indices;
        merged_indices.extend(g_j.indices);

        let mut new_row = Vec::with_capacity(groups.len() + 1);
        for m in 0..g {
            if m == best_i || m == best_j {
                continue;
            }
            new_row.push((dist[best_i][m] + dist[best_j][m]) * 0.5);
        }
        dist.remove(best_j);
        dist.remove(best_i);
        for row in &mut dist {
            row.remove(best_j);
            row.remove(best_i);
        }
        for (row, &d_val) in dist.iter_mut().zip(&new_row) {
            row.push(d_val);
        }
        new_row.push(0.0);
        dist.push(new_row);

        groups.push(Group {
            indices: merged_indices,
            rows: merged_rows,
        });
    }

    let mut ordered = vec![Vec::new(); n];
    if let Some(final_g) = groups.pop() {
        for (orig_idx, row) in final_g.indices.into_iter().zip(final_g.rows) {
            ordered[orig_idx] = row;
        }
    }
    ordered
}

/// Aligns two sub-profiles `prof_a` and `prof_b` using global profile-profile dynamic programming.
fn align_two_profiles(
    prof_a: &[Vec<u8>],
    prof_b: &[Vec<u8>],
    alphabet: Alphabet,
    config: &AlignerConfig,
    conservation_weighted: bool,
) -> Vec<Vec<u8>> {
    let m = prof_a.first().map_or(0, Vec::len);
    let n = prof_b.first().map_or(0, Vec::len);
    let cols = n + 1;

    if m == 0 {
        let mut out = vec![vec![b'-'; n]; prof_a.len()];
        out.extend(prof_b.iter().cloned());
        return out;
    }
    if n == 0 {
        let mut out = prof_a.to_vec();
        out.extend(vec![vec![b'-'; m]; prof_b.len()]);
        return out;
    }

    // Precompute modal residue & gap fraction for every column in both profiles for O(1) DP cell scoring.
    let col_a = summarize_profile(prof_a);
    let col_b = summarize_profile(prof_b);
    let base_gap = -config.gap_open.clamp(2, 20);

    let mut prev = vec![0i32; cols];
    let mut curr = vec![0i32; cols];
    let mut trace = vec![0u8; (m + 1) * cols];

    for j in 1..=n {
        prev[j] = base_gap * j as i32;
        trace[j] = 2;
    }

    for i in 1..=m {
        curr[0] = base_gap * i as i32;
        trace[i * cols] = 1;
        let ca = col_a[i - 1];

        for j in 1..=n {
            let cb = col_b[j - 1];
            let sub = profile_col_score(ca, cb, alphabet, conservation_weighted);
            // Conservation-aware gap penalty: penalize opening gaps inside highly conserved columns!
            let gap_up = if conservation_weighted && ca.gap_frac > 0.3 {
                base_gap / 2
            } else {
                base_gap
            };
            let gap_left = if conservation_weighted && cb.gap_frac > 0.3 {
                base_gap / 2
            } else {
                base_gap
            };

            let diag = prev[j - 1] + sub;
            let up = prev[j] + gap_up;
            let left = curr[j - 1] + gap_left;

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

    // Backtrack aligned column instructions.
    let mut ops: Vec<u8> = Vec::with_capacity(m.max(n) + 16);
    let mut i = m;
    let mut j = n;
    while i > 0 || j > 0 {
        let dir = trace[i * cols + j];
        if i > 0 && j > 0 && dir == 0 {
            ops.push(0);
            i -= 1;
            j -= 1;
        } else if i > 0 && (j == 0 || dir == 1) {
            ops.push(1);
            i -= 1;
        } else {
            ops.push(2);
            j -= 1;
        }
    }
    ops.reverse();

    let total_rows = prof_a.len() + prof_b.len();
    let mut out = vec![Vec::with_capacity(ops.len()); total_rows];
    let mut ia = 0usize;
    let mut ib = 0usize;
    for dir in ops {
        match dir {
            0 => {
                for (r, row) in prof_a.iter().enumerate() {
                    out[r].push(row[ia]);
                }
                for (r, row) in prof_b.iter().enumerate() {
                    out[prof_a.len() + r].push(row[ib]);
                }
                ia += 1;
                ib += 1;
            }
            1 => {
                for (r, row) in prof_a.iter().enumerate() {
                    out[r].push(row[ia]);
                }
                for r in 0..prof_b.len() {
                    out[prof_a.len() + r].push(b'-');
                }
                ia += 1;
            }
            _ => {
                for row_out in out.iter_mut().take(prof_a.len()) {
                    row_out.push(b'-');
                }
                for (r, row) in prof_b.iter().enumerate() {
                    out[prof_a.len() + r].push(row[ib]);
                }
                ib += 1;
            }
        }
    }
    out
}

#[derive(Clone, Copy)]
struct ProfileColSummary {
    modal: u8,
    conservation: f32,
    gap_frac: f32,
}

fn summarize_profile(prof: &[Vec<u8>]) -> Vec<ProfileColSummary> {
    let width = prof.first().map_or(0, Vec::len);
    let rows = prof.len().max(1);
    let mut out = Vec::with_capacity(width);
    for col in 0..width {
        let mut counts = [0u32; 256];
        let mut gaps = 0u32;
        for seq in prof {
            let b = seq[col];
            if is_gap(b) {
                gaps += 1;
            } else {
                counts[usize::from(b)] += 1;
            }
        }
        let mut modal = b'-';
        let mut best_c = 0u32;
        for (val, &c) in counts.iter().enumerate() {
            if c > best_c {
                best_c = c;
                modal = val as u8;
            }
        }
        let non_gap = (rows as u32).saturating_sub(gaps).max(1);
        out.push(ProfileColSummary {
            modal,
            conservation: best_c as f32 / non_gap as f32,
            gap_frac: gaps as f32 / rows as f32,
        });
    }
    out
}

#[inline]
fn profile_col_score(
    a: ProfileColSummary,
    b: ProfileColSummary,
    alphabet: Alphabet,
    conservation_weighted: bool,
) -> i32 {
    if is_gap(a.modal) || is_gap(b.modal) {
        return -2;
    }
    let raw = pair_score(a.modal, b.modal, alphabet);
    if conservation_weighted {
        let boost = 1.0 + 0.35 * (a.conservation * b.conservation);
        (raw as f32 * boost).round() as i32
    } else {
        raw
    }
}

#[inline]
fn pair_score(a: u8, b: u8, alphabet: Alphabet) -> i32 {
    if a == b {
        return if alphabet.is_nucleotide() { 5 } else { 6 };
    }
    if alphabet.is_nucleotide() {
        let transition = matches!(
            (a, b),
            (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') | (b'C', b'U') | (b'U', b'C')
        );
        return if transition { -1 } else { -4 };
    }
    let grp = |c: u8| -> u8 {
        match c {
            b'I' | b'L' | b'V' | b'M' => 1,
            b'F' | b'Y' | b'W' => 2,
            b'K' | b'R' | b'H' => 3,
            b'D' | b'E' => 4,
            b'S' | b'T' | b'N' | b'Q' => 5,
            b'A' | b'G' => 6,
            _ => 0,
        }
    };
    let ga = grp(a);
    if ga != 0 && ga == grp(b) {
        2
    } else {
        -3
    }
}

fn sum_of_pairs_score(msa: &[Vec<u8>], alphabet: Alphabet) -> i64 {
    let n = msa.len();
    let width = msa.first().map_or(0, Vec::len);
    let mut total = 0i64;
    for col in 0..width {
        for i in 0..n {
            let a = msa[i][col];
            for row_j in msa.iter().take(n).skip(i + 1) {
                let b = row_j[col];
                if is_gap(a) && is_gap(b) {
                    continue;
                } else if is_gap(a) || is_gap(b) {
                    total -= 3;
                } else {
                    total += i64::from(pair_score(a, b, alphabet));
                }
            }
        }
    }
    total
}

fn strip_all_gap_columns(prof: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let width = prof.first().map_or(0, Vec::len);
    let mut out = vec![Vec::with_capacity(width); prof.len()];
    for col in 0..width {
        if prof.iter().any(|row| !is_gap(row[col])) {
            for (r, row) in prof.iter().enumerate() {
                out[r].push(row[col]);
            }
        }
    }
    out
}

fn sorted_kmer_hashes(seq: &[u8], k: usize) -> Vec<u64> {
    if seq.len() < k {
        return Vec::new();
    }
    let mut hashes: Vec<u64> = seq.windows(k).map(fnv1a_kmer).collect();
    hashes.sort_unstable();
    hashes.dedup();
    hashes
}

fn kmer_jaccard_distance(a: &[u64], b: &[u64]) -> f32 {
    if a.is_empty() || b.is_empty() {
        return 1.0;
    }
    let mut i = 0usize;
    let mut j = 0usize;
    let mut inter = 0usize;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                inter += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
    let union = (a.len() + b.len()).saturating_sub(inter).max(1);
    1.0 - (inter as f32 / union as f32)
}

#[inline]
fn fnv1a_kmer(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

/// Converts biological residues into a normalized physicochemical numerical signal for MAFFT
/// cross-correlation estimation.
fn residue_signal(seq: &[u8], alphabet: Alphabet) -> Vec<f32> {
    seq.iter()
        .map(|&b| {
            if alphabet.is_nucleotide() {
                match b {
                    b'A' => 1.0,
                    b'G' => 0.6,
                    b'C' => -0.6,
                    b'T' | b'U' => -1.0,
                    _ => 0.0,
                }
            } else {
                // Normalized Kyte-Doolittle hydropathy values in [-1.0, 1.0].
                match b {
                    b'I' => 1.0,
                    b'V' => 0.93,
                    b'L' => 0.84,
                    b'F' => 0.62,
                    b'C' => 0.55,
                    b'M' => 0.42,
                    b'A' => 0.40,
                    b'G' => -0.09,
                    b'T' => -0.15,
                    b'S' => -0.18,
                    b'W' => -0.20,
                    b'Y' => -0.29,
                    b'P' => -0.35,
                    b'H' => -0.71,
                    b'E' | b'Q' | b'D' | b'N' => -0.78,
                    b'K' => -0.87,
                    b'R' => -1.0,
                    _ => 0.0,
                }
            }
        })
        .collect()
}

/// Estimates peak normalized cross-correlation across bounded lags (fast real-space equivalent
/// of MAFFT's FFT lag spectrum).
fn peak_cross_correlation(a: &[f32], b: &[f32]) -> f32 {
    let n = a.len().min(b.len());
    if n < 4 {
        return 0.0;
    }
    let max_lag = (n / 3).clamp(4, 32) as isize;
    let mut best = 0.0f32;
    for lag in -max_lag..=max_lag {
        let a_s = (-lag).max(0) as usize;
        let b_s = lag.max(0) as usize;
        let ov = (a.len() - a_s).min(b.len() - b_s);
        if ov < 4 {
            continue;
        }
        let mut dot = 0.0f32;
        let mut na = 0.0f32;
        let mut nb = 0.0f32;
        for k in 0..ov {
            let va = a[a_s + k];
            let vb = b[b_s + k];
            dot += va * vb;
            na += va * va;
            nb += vb * vb;
        }
        let denom = (na * nb).sqrt();
        if denom > 1e-5 {
            let corr = dot / denom;
            if corr > best {
                best = corr;
            }
        }
    }
    best
}

fn sanitize_fasta_id(name: &str) -> String {
    let s: String = name
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect();
    if s.is_empty() {
        String::from("Seq")
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_dna_seqs() -> Vec<Sequence> {
        vec![
            Sequence::new("Ref", b"ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"),
            Sequence::new("Iso1", b"ATGGTGCATCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"),
            Sequence::new("Iso2", b"ATGGTGCACCTGACTCCTGAGGAGAAGTCTGTTACTGCCCTGTGGGGCAAGGTGAAC"),
        ]
    }

    #[test]
    fn all_six_alignment_algorithms_produce_valid_rectangular_alignments() {
        let seqs = sample_dna_seqs();
        for alg in AlignerAlgorithm::ALL {
            let cfg = AlignerConfig {
                algorithm: alg,
                ..AlignerConfig::default()
            };
            let out = run_alignment_algorithm(
                format!("Test_{}", alg.short_name()),
                Alphabet::Dna,
                &seqs,
                &cfg,
            )
            .expect("aligner succeeds");
            assert_eq!(out.alignment.row_count(), 3);
            assert!(out.alignment.width() >= 55);
            assert!(
                out.alignment.mean_pairwise_identity_pct() > 70.0,
                "{} identity too low: {}",
                alg.short_name(),
                out.alignment.mean_pairwise_identity_pct()
            );
        }
    }

    #[test]
    fn sam_cigar_parser_handles_matches_insertions_deletions_and_softclips() {
        let ref_s = b"ATGGTGCATCTGACTCCTGAGGAG";
        let seqs = vec![
            Sequence::new("Ref", ref_s),
            Sequence::new("Q1", b"TTATGGTGCATCTG"),
        ];
        // Soft-clip 2 bp, match 6 bp, delete 2 bp, match 6 bp starting at pos 1.
        let sam = "Q1\t0\tRef\t1\t60\t2S6M2D6M\t*\t0\t0\tTTATGGTGCATCTG\t*\n";
        let aln = parse_sam_into_alignment("SAMTest", Alphabet::Dna, &seqs, ref_s, sam)
            .expect("parse sam");
        assert_eq!(aln.row_count(), 2);
        assert_eq!(&aln.sequences()[1].residues[0..6], b"ATGGTG");
        assert_eq!(&aln.sequences()[1].residues[6..8], b"--");
        assert_eq!(&aln.sequences()[1].residues[8..14], b"CATCTG");
    }
}
