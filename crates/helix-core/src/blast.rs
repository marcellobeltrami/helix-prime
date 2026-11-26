//! Multi-program BLAST similarity search engine supporting:
//! - **totalBLAST**: Unified cross-alphabet search (runs nucleotide + 6-frame translated protein
//!   searches simultaneously across all workspace documents!).
//! - **blastn**: Nucleotide query vs Nucleotide database (both forward & reverse strands).
//! - **blastp**: Protein query vs Protein database (BLOSUM-inspired substitution scoring).
//! - **blastx**: 6-frame translated nucleotide query vs Protein database.
//! - **tblastn**: Protein query vs 6-frame translated nucleotide database.
//! - **tblastx**: 6-frame translated nucleotide query vs 6-frame translated nucleotide database.
//! - **megablast**: Fast long-word (`k = 11`) nucleotide search for closely related strains/phages.

use crate::alignment::Alignment;
use crate::sequence::{is_gap, reverse_complement, translate_ungapped, Alphabet};
use serde::{Deserialize, Serialize};

/// Supported BLAST search program variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum BlastProgram {
    /// Unified all-vs-all search: automatically runs blastn + blastx (for DNA/RNA query)
    /// or blastp + tblastn (for Protein query) across every document in the workspace.
    #[default]
    TotalBlast,
    /// Nucleotide query vs Nucleotide database (forward & reverse-complement).
    BlastN,
    /// Protein query vs Protein database.
    BlastP,
    /// 6-frame translated nucleotide query vs Protein database.
    BlastX,
    /// Protein query vs 6-frame translated nucleotide database.
    TBlastN,
    /// 6-frame translated nucleotide query vs 6-frame translated nucleotide database.
    TBlastX,
    /// Fast exact long-word (`k = 11`) nucleotide search for high-identity homologs.
    MegaBlast,
}

impl BlastProgram {
    pub const ALL: [Self; 7] = [
        Self::TotalBlast,
        Self::BlastN,
        Self::BlastP,
        Self::BlastX,
        Self::TBlastN,
        Self::TBlastX,
        Self::MegaBlast,
    ];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::TotalBlast => "totalBLAST (All DNA + 6-Frame Protein)",
            Self::BlastN => "blastn (Nucleotide -> Nucleotide)",
            Self::BlastP => "blastp (Protein -> Protein)",
            Self::BlastX => "blastx (Translated DNA -> Protein)",
            Self::TBlastN => "tblastn (Protein -> Translated DNA)",
            Self::TBlastX => "tblastx (Translated DNA -> Translated DNA)",
            Self::MegaBlast => "MegaBLAST (Fast High-Identity DNA)",
        }
    }

    #[must_use]
    pub const fn short_tag(self) -> &'static str {
        match self {
            Self::TotalBlast => "totalBLAST",
            Self::BlastN => "blastn",
            Self::BlastP => "blastp",
            Self::BlastX => "blastx",
            Self::TBlastN => "tblastn",
            Self::TBlastX => "tblastx",
            Self::MegaBlast => "megablast",
        }
    }
}

/// Configurable parameters for a BLAST run.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BlastConfig {
    pub program: BlastProgram,
    /// Minimum percent identity (`0.0..=100.0`) required to report a hit.
    pub min_identity_pct: f32,
    /// Maximum E-value threshold (hits with higher E-value are filtered out).
    pub max_e_value: f64,
    /// Maximum number of top-ranked hits to return.
    pub max_hits: usize,
}

impl Default for BlastConfig {
    fn default() -> Self {
        Self {
            program: BlastProgram::TotalBlast,
            min_identity_pct: 45.0,
            max_e_value: 10.0,
            max_hits: 30,
        }
    }
}

/// A single BLAST High-Scoring Segment Pair (HSP) hit.
#[derive(Debug, Clone)]
pub struct BlastHit {
    /// Specific BLAST mode & strand/frame that produced this hit (e.g. `"blastn (+)"`, `"blastx (+2)"`).
    pub program_tag: String,
    pub doc_index: usize,
    pub doc_name: String,
    pub subject_name: String,
    pub subject_row: usize,
    pub query_start: usize,
    pub query_end: usize,
    /// 0-based start coordinate in the original subject sequence (mapped back to DNA bp for tblastn/tblastx!).
    pub subject_start: usize,
    /// 0-based exclusive end coordinate in the original subject sequence.
    pub subject_end: usize,
    pub align_len: usize,
    pub identity_pct: f32,
    pub positives_pct: f32,
    pub bit_score: f32,
    pub e_value: f64,
    pub query_snippet: String,
    pub match_midline: String,
    pub subject_snippet: String,
}

/// Runs a BLAST search using the selected [`BlastProgram`] (`totalBLAST`, `blastn`, `blastp`,
/// `blastx`, `tblastn`, `tblastx`, or `megablast`) across all provided workspace alignments.
#[must_use]
pub fn run_blast(
    raw_query: &str,
    subjects: &[(usize, &Alignment)],
    config: BlastConfig,
) -> Vec<BlastHit> {
    let query: Vec<u8> = raw_query
        .bytes()
        .filter(|b| !b.is_ascii_whitespace() && !is_gap(*b))
        .map(|b| b.to_ascii_uppercase())
        .collect();

    if query.len() < 4 {
        return Vec::new();
    }

    let query_alpha = Alphabet::infer_from_slices([query.as_slice()]);
    let total_db_len: usize = subjects
        .iter()
        .flat_map(|(_, a)| a.sequences().iter())
        .map(crate::sequence::Sequence::ungapped_len)
        .sum::<usize>()
        .max(100);

    let mut hits: Vec<BlastHit> = Vec::new();

    match config.program {
        BlastProgram::TotalBlast => {
            if query_alpha.is_nucleotide() {
                // Run both blastn (+/-) against DNA/RNA and 6-frame blastx against Protein!
                collect_blastn(&query, subjects, 6, "blastn", total_db_len, config, &mut hits);
                collect_blastx(&query, subjects, total_db_len, config, &mut hits);
            } else {
                // Run both blastp against Protein and 6-frame tblastn against DNA/RNA!
                collect_blastp(&query, subjects, total_db_len, config, &mut hits);
                collect_tblastn(&query, subjects, total_db_len, config, &mut hits);
            }
        }
        BlastProgram::BlastN => {
            collect_blastn(&query, subjects, 6, "blastn", total_db_len, config, &mut hits);
        }
        BlastProgram::MegaBlast => {
            collect_blastn(
                &query,
                subjects,
                11,
                "megablast",
                total_db_len,
                config,
                &mut hits,
            );
        }
        BlastProgram::BlastP => {
            collect_blastp(&query, subjects, total_db_len, config, &mut hits);
        }
        BlastProgram::BlastX => {
            collect_blastx(&query, subjects, total_db_len, config, &mut hits);
        }
        BlastProgram::TBlastN => {
            collect_tblastn(&query, subjects, total_db_len, config, &mut hits);
        }
        BlastProgram::TBlastX => {
            collect_tblastx(&query, subjects, total_db_len, config, &mut hits);
        }
    }

    hits.sort_by(|a, b| b.bit_score.total_cmp(&a.bit_score));
    hits.truncate(config.max_hits.max(1));
    hits
}

/// Backward-compatible convenience entry point running `totalBLAST` with default settings.
#[must_use]
pub fn run_local_blast(
    raw_query: &str,
    subjects: &[(usize, &Alignment)],
    max_hits: usize,
) -> Vec<BlastHit> {
    run_blast(
        raw_query,
        subjects,
        BlastConfig {
            max_hits,
            ..BlastConfig::default()
        },
    )
}

/// Nucleotide vs Nucleotide (`blastn` / `megablast`) scanning both forward and reverse-complement strands.
fn collect_blastn(
    query: &[u8],
    subjects: &[(usize, &Alignment)],
    k: usize,
    tag_prefix: &str,
    total_db_len: usize,
    config: BlastConfig,
    hits: &mut Vec<BlastHit>,
) {
    if query.len() < k {
        return;
    }
    let rc_query = reverse_complement(query);

    for &(doc_idx, aln) in subjects {
        if !aln.alphabet.is_nucleotide() {
            continue;
        }
        for (row_idx, seq) in aln.sequences().iter().enumerate() {
            let subj = seq.ungapped_bytes();
            if subj.len() < k {
                continue;
            }
            let fwd = find_best_hsp(query, &subj, k, true);
            let rev = find_best_hsp(&rc_query, &subj, k, true);

            let best = match (fwd, rev) {
                (Some(f), Some(r)) => {
                    if r.score > f.score {
                        Some((r, &rc_query[..], format!("{tag_prefix} (-)")))
                    } else {
                        Some((f, query, format!("{tag_prefix} (+)")))
                    }
                }
                (Some(f), None) => Some((f, query, format!("{tag_prefix} (+)"))) ,
                (None, Some(r)) => Some((r, &rc_query[..], format!("{tag_prefix} (-)"))) ,
                (None, None) => None,
            };

            if let Some((hsp, used_q, tag)) = best {
                if let Some(hit) = build_hit(
                    tag,
                    doc_idx,
                    aln,
                    row_idx,
                    &seq.name,
                    used_q,
                    &subj,
                    &hsp,
                    hsp.s_start,
                    hsp.s_end,
                    true,
                    total_db_len,
                    config,
                ) {
                    hits.push(hit);
                }
            }
        }
    }
}

/// Protein vs Protein (`blastp`) using physicochemical substitution scoring.
fn collect_blastp(
    query: &[u8],
    subjects: &[(usize, &Alignment)],
    total_db_len: usize,
    config: BlastConfig,
    hits: &mut Vec<BlastHit>,
) {
    let k = 3usize;
    if query.len() < k {
        return;
    }
    for &(doc_idx, aln) in subjects {
        if aln.alphabet.is_nucleotide() {
            continue;
        }
        for (row_idx, seq) in aln.sequences().iter().enumerate() {
            let subj = seq.ungapped_bytes();
            if subj.len() < k {
                continue;
            }
            if let Some(hsp) = find_best_hsp(query, &subj, k, false) {
                if let Some(hit) = build_hit(
                    String::from("blastp"),
                    doc_idx,
                    aln,
                    row_idx,
                    &seq.name,
                    query,
                    &subj,
                    &hsp,
                    hsp.s_start,
                    hsp.s_end,
                    false,
                    total_db_len,
                    config,
                ) {
                    hits.push(hit);
                }
            }
        }
    }
}

/// 6-frame translated DNA query vs Protein database (`blastx`).
fn collect_blastx(
    dna_query: &[u8],
    subjects: &[(usize, &Alignment)],
    total_db_len: usize,
    config: BlastConfig,
    hits: &mut Vec<BlastHit>,
) {
    let frames = translate_six_frames(dna_query);
    let k = 3usize;

    for &(doc_idx, aln) in subjects {
        if aln.alphabet.is_nucleotide() {
            continue;
        }
        for (row_idx, seq) in aln.sequences().iter().enumerate() {
            let subj = seq.ungapped_bytes();
            if subj.len() < k {
                continue;
            }
            let mut best: Option<(Hsp, usize)> = None;
            for (f_idx, (_, prot_q)) in frames.iter().enumerate() {
                if prot_q.len() < k {
                    continue;
                }
                if let Some(hsp) = find_best_hsp(prot_q, &subj, k, false) {
                    if best.as_ref().is_none_or(|(b, _)| hsp.score > b.score) {
                        best = Some((hsp, f_idx));
                    }
                }
            }
            if let Some((hsp, f_idx)) = best {
                let (frame_lbl, ref prot_q) = frames[f_idx];
                if let Some(hit) = build_hit(
                    format!("blastx ({frame_lbl})"),
                    doc_idx,
                    aln,
                    row_idx,
                    &seq.name,
                    prot_q,
                    &subj,
                    &hsp,
                    hsp.s_start,
                    hsp.s_end,
                    false,
                    total_db_len,
                    config,
                ) {
                    hits.push(hit);
                }
            }
        }
    }
}

/// Protein query vs 6-frame translated DNA database (`tblastn`).
fn collect_tblastn(
    prot_query: &[u8],
    subjects: &[(usize, &Alignment)],
    total_db_len: usize,
    config: BlastConfig,
    hits: &mut Vec<BlastHit>,
) {
    let k = 3usize;
    if prot_query.len() < k {
        return;
    }

    for &(doc_idx, aln) in subjects {
        if !aln.alphabet.is_nucleotide() {
            continue;
        }
        for (row_idx, seq) in aln.sequences().iter().enumerate() {
            let dna = seq.ungapped_bytes();
            if dna.len() < 12 {
                continue;
            }
            let subj_frames = translate_six_frames(&dna);
            let mut best: Option<(Hsp, usize)> = None;
            for (f_idx, (_, prot_s)) in subj_frames.iter().enumerate() {
                if prot_s.len() < k {
                    continue;
                }
                if let Some(hsp) = find_best_hsp(prot_query, prot_s, k, false) {
                    if best.as_ref().is_none_or(|(b, _)| hsp.score > b.score) {
                        best = Some((hsp, f_idx));
                    }
                }
            }
            if let Some((hsp, f_idx)) = best {
                let (frame_lbl, ref prot_s) = subj_frames[f_idx];
                let (dna_s, dna_e) = map_frame_to_dna_span(f_idx, hsp.s_start, hsp.s_end, dna.len());
                if let Some(hit) = build_hit(
                    format!("tblastn ({frame_lbl})"),
                    doc_idx,
                    aln,
                    row_idx,
                    &seq.name,
                    prot_query,
                    prot_s,
                    &hsp,
                    dna_s,
                    dna_e,
                    false,
                    total_db_len,
                    config,
                ) {
                    hits.push(hit);
                }
            }
        }
    }
}

/// 6-frame translated DNA query vs 6-frame translated DNA database (`tblastx`).
fn collect_tblastx(
    dna_query: &[u8],
    subjects: &[(usize, &Alignment)],
    total_db_len: usize,
    config: BlastConfig,
    hits: &mut Vec<BlastHit>,
) {
    let q_frames = translate_six_frames(dna_query);
    let k = 3usize;

    for &(doc_idx, aln) in subjects {
        if !aln.alphabet.is_nucleotide() {
            continue;
        }
        for (row_idx, seq) in aln.sequences().iter().enumerate() {
            let dna = seq.ungapped_bytes();
            if dna.len() < 12 {
                continue;
            }
            let s_frames = translate_six_frames(&dna);
            let mut best: Option<(Hsp, usize, usize)> = None;
            for (qi, (_, pq)) in q_frames.iter().enumerate() {
                if pq.len() < k {
                    continue;
                }
                for (si, (_, ps)) in s_frames.iter().enumerate() {
                    if ps.len() < k {
                        continue;
                    }
                    if let Some(hsp) = find_best_hsp(pq, ps, k, false) {
                        if best.as_ref().is_none_or(|(b, _, _)| hsp.score > b.score) {
                            best = Some((hsp, qi, si));
                        }
                    }
                }
            }
            if let Some((hsp, qi, si)) = best {
                let (q_lbl, ref pq) = q_frames[qi];
                let (s_lbl, ref ps) = s_frames[si];
                let (dna_s, dna_e) = map_frame_to_dna_span(si, hsp.s_start, hsp.s_end, dna.len());
                if let Some(hit) = build_hit(
                    format!("tblastx ({q_lbl}/{s_lbl})"),
                    doc_idx,
                    aln,
                    row_idx,
                    &seq.name,
                    pq,
                    ps,
                    &hsp,
                    dna_s,
                    dna_e,
                    false,
                    total_db_len,
                    config,
                ) {
                    hits.push(hit);
                }
            }
        }
    }
}

/// Translates a DNA sequence in all 6 reading frames (`+1`, `+2`, `+3`, `-1`, `-2`, `-3`).
fn translate_six_frames(dna: &[u8]) -> [(&'static str, Vec<u8>); 6] {
    let rc = reverse_complement(dna);
    [
        ("+1", translate_ungapped(dna, 0)),
        ("+2", translate_ungapped(dna, 1)),
        ("+3", translate_ungapped(dna, 2)),
        ("-1", translate_ungapped(&rc, 0)),
        ("-2", translate_ungapped(&rc, 1)),
        ("-3", translate_ungapped(&rc, 2)),
    ]
}

/// Maps amino-acid coordinates in reading frame `frame_idx` (`0..=5`) back to 0-based DNA bp coordinates.
fn map_frame_to_dna_span(frame_idx: usize, aa_start: usize, aa_end: usize, dna_len: usize) -> (usize, usize) {
    let offset = frame_idx % 3;
    let raw_start = offset + aa_start * 3;
    let raw_end = (offset + aa_end * 3).min(dna_len);
    if frame_idx < 3 {
        (raw_start.min(dna_len), raw_end.max(raw_start.min(dna_len)))
    } else {
        let rev_start = dna_len.saturating_sub(raw_end);
        let rev_end = dna_len.saturating_sub(raw_start).max(rev_start);
        (rev_start, rev_end)
    }
}

#[allow(clippy::too_many_arguments)]
fn build_hit(
    program_tag: String,
    doc_index: usize,
    aln: &Alignment,
    subject_row: usize,
    subject_name: &str,
    query: &[u8],
    subject: &[u8],
    hsp: &Hsp,
    mapped_subj_start: usize,
    mapped_subj_end: usize,
    is_nuc: bool,
    total_db_len: usize,
    config: BlastConfig,
) -> Option<BlastHit> {
    let alen = hsp.q_end.saturating_sub(hsp.q_start);
    if alen < 3 {
        return None;
    }

    let q_sub = &query[hsp.q_start..hsp.q_end];
    let s_sub = &subject[hsp.s_start..hsp.s_end];

    let mut exact = 0usize;
    let mut positives = 0usize;
    let mut midline = String::with_capacity(alen);
    for (&qb, &sb) in q_sub.iter().zip(s_sub.iter()) {
        if qb == sb {
            exact += 1;
            positives += 1;
            midline.push('|');
        } else if !is_nuc && residue_score(qb, sb, false) > 0 {
            positives += 1;
            midline.push('+');
        } else {
            midline.push(' ');
        }
    }

    let identity_pct = (exact as f32 / alen as f32) * 100.0;
    let positives_pct = (positives as f32 / alen as f32) * 100.0;
    if identity_pct < config.min_identity_pct {
        return None;
    }

    // Karlin-Altschul bit-score scaling (lambda / ln(2) ≈ 0.3176).
    let bit_score = (hsp.score as f32 * 0.3176 - 0.14).max(1.0);
    let search_space = (query.len() as f64) * (total_db_len as f64);
    let e_value = search_space * 2.0f64.powf(-f64::from(bit_score));
    if e_value > config.max_e_value {
        return None;
    }

    Some(BlastHit {
        program_tag,
        doc_index,
        doc_name: aln.name.clone(),
        subject_name: subject_name.to_owned(),
        subject_row,
        query_start: hsp.q_start,
        query_end: hsp.q_end,
        subject_start: mapped_subj_start,
        subject_end: mapped_subj_end,
        align_len: alen,
        identity_pct,
        positives_pct,
        bit_score,
        e_value,
        query_snippet: String::from_utf8_lossy(q_sub).into_owned(),
        match_midline: midline,
        subject_snippet: String::from_utf8_lossy(s_sub).into_owned(),
    })
}

struct Hsp {
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    score: i32,
}

fn find_best_hsp(query: &[u8], subject: &[u8], k: usize, is_nuc: bool) -> Option<Hsp> {
    if query.len() < k || subject.len() < k {
        return None;
    }
    let x_drop = 16i32;
    let mut best: Option<Hsp> = None;

    let step = (query.len() / 32).max(1);
    let mut q_i = 0usize;
    while q_i + k <= query.len() {
        let seed = &query[q_i..(q_i + k)];
        // Skip low-complexity homopolymer seeds (e.g., AAAAAA) for cleaner hits.
        if seed.iter().all(|&b| b == seed[0]) {
            q_i += step;
            continue;
        }
        for s_i in 0..=(subject.len() - k) {
            if &subject[s_i..(s_i + k)] != seed {
                continue;
            }
            let mut seed_score = 0i32;
            for d in 0..k {
                seed_score += residue_score(query[q_i + d], subject[s_i + d], is_nuc);
            }

            // Extend right with X-drop.
            let mut best_r = k;
            let mut cur_r_score = seed_score;
            let mut max_r_score = cur_r_score;
            let max_ext_r = (query.len() - q_i).min(subject.len() - s_i);
            for d in k..max_ext_r {
                cur_r_score += residue_score(query[q_i + d], subject[s_i + d], is_nuc);
                if cur_r_score > max_r_score {
                    max_r_score = cur_r_score;
                    best_r = d + 1;
                } else if max_r_score - cur_r_score > x_drop {
                    break;
                }
            }

            // Extend left with X-drop.
            let mut best_l = 0usize;
            let mut cur_l_score = 0i32;
            let mut max_l_score = 0i32;
            let max_ext_l = q_i.min(s_i);
            for d in 1..=max_ext_l {
                cur_l_score += residue_score(query[q_i - d], subject[s_i - d], is_nuc);
                if cur_l_score > max_l_score {
                    max_l_score = cur_l_score;
                    best_l = d;
                } else if max_l_score - cur_l_score > x_drop {
                    break;
                }
            }

            let total_score = max_r_score + max_l_score;
            if best.as_ref().is_none_or(|b| total_score > b.score) {
                best = Some(Hsp {
                    q_start: q_i - best_l,
                    q_end: q_i + best_r,
                    s_start: s_i - best_l,
                    s_end: s_i + best_r,
                    score: total_score,
                });
            }
        }
        q_i += step;
    }

    best
}

#[inline]
fn residue_score(a: u8, b: u8, is_nuc: bool) -> i32 {
    if a == b {
        return if is_nuc { 2 } else { 5 };
    }
    if is_nuc {
        return -3;
    }
    // Conservative amino-acid similarity groups (BLOSUM62-inspired positive substitutions).
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
        -2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::demo::build_default_workspace;

    #[test]
    fn total_blast_finds_both_protein_and_translated_dna_hits() {
        let docs = build_default_workspace();
        let subjects: Vec<(usize, &Alignment)> = docs
            .iter()
            .enumerate()
            .map(|(i, d)| (i, &d.alignment))
            .collect();

        // Protein query should match BOTH Hemoglobin protein docs (blastp) AND HBB DNA CDS (tblastn)!
        let hits = run_blast(
            "LSELHCDKLHVDPENFRLLGNVLV",
            &subjects,
            BlastConfig {
                program: BlastProgram::TotalBlast,
                ..BlastConfig::default()
            },
        );
        assert!(hits.iter().any(|h| h.program_tag == "blastp"));
        assert!(hits.iter().any(|h| h.program_tag.starts_with("tblastn")));
    }

    #[test]
    fn blastx_and_megablast_work_on_dna_query() {
        let docs = build_default_workspace();
        let subjects: Vec<(usize, &Alignment)> = docs
            .iter()
            .enumerate()
            .map(|(i, d)| (i, &d.alignment))
            .collect();

        let dna_query = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC";
        let bx = run_blast(
            dna_query,
            &subjects,
            BlastConfig {
                program: BlastProgram::BlastX,
                ..BlastConfig::default()
            },
        );
        assert!(bx.iter().any(|h| h.program_tag.starts_with("blastx")));

        let mb = run_blast(
            dna_query,
            &subjects,
            BlastConfig {
                program: BlastProgram::MegaBlast,
                ..BlastConfig::default()
            },
        );
        assert!(mb.iter().any(|h| h.program_tag.starts_with("megablast")));
    }
}
