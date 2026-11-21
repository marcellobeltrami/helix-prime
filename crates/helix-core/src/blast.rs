//! Local BLAST-style similarity search engine (k-mer seeding + ungapped HSP extension +
//! bit score & E-value estimation) across workspace sequences.

use crate::alignment::Alignment;
use crate::sequence::{is_gap, Alphabet};

/// A single local BLAST High-Scoring Segment Pair (HSP) match.
#[derive(Debug, Clone)]
pub struct BlastHit {
    pub doc_index: usize,
    pub doc_name: String,
    pub subject_name: String,
    pub subject_row: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub align_len: usize,
    pub identity_pct: f32,
    pub bit_score: f32,
    pub e_value: f64,
    pub query_snippet: String,
    pub match_midline: String,
    pub subject_snippet: String,
}

/// Runs a fast pure-Rust local BLAST search of `raw_query` across all documents in `subjects`.
///
/// Uses k-mer seed indexing (`k = 6` for DNA/RNA, `k = 3` for proteins) followed by bidirectional
/// ungapped HSP extension (X-drop), ranking results by bit score descending.
#[must_use]
pub fn run_local_blast(
    raw_query: &str,
    subjects: &[(usize, &Alignment)],
    max_hits: usize,
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
    let k: usize = if query_alpha.is_nucleotide() { 6 } else { 3 };
    if query.len() < k {
        return Vec::new();
    }

    let total_db_len: usize = subjects
        .iter()
        .flat_map(|(_, a)| a.sequences().iter())
        .map(crate::sequence::Sequence::ungapped_len)
        .sum::<usize>()
        .max(100);

    let mut hits: Vec<BlastHit> = Vec::new();

    for &(doc_idx, aln) in subjects {
        // Match compatible alphabets (DNA/RNA together, Protein with Protein).
        if aln.alphabet.is_nucleotide() != query_alpha.is_nucleotide() {
            continue;
        }
        for (row_idx, seq) in aln.sequences().iter().enumerate() {
            let subj = seq.ungapped_bytes();
            if subj.len() < k {
                continue;
            }
            if let Some(best_hsp) = find_best_hsp(&query, &subj, k, query_alpha.is_nucleotide()) {
                let raw_score = best_hsp.score as f32;
                // Karlin-Altschul bit-score scaling (lambda / ln(2) ≈ 0.3176).
                let bit_score = (raw_score * 0.3176 - 0.14).max(1.0);
                let search_space = (query.len() as f64) * (total_db_len as f64);
                let e_value = search_space * 2.0f64.powf(-f64::from(bit_score));

                let alen = best_hsp.q_end.saturating_sub(best_hsp.q_start);
                if alen == 0 {
                    continue;
                }
                let q_sub = &query[best_hsp.q_start..best_hsp.q_end];
                let s_sub = &subj[best_hsp.s_start..best_hsp.s_end];

                let mut matches = 0usize;
                let mut midline = String::with_capacity(alen);
                for (&qb, &sb) in q_sub.iter().zip(s_sub.iter()) {
                    if qb == sb {
                        matches += 1;
                        midline.push('|');
                    } else {
                        midline.push(' ');
                    }
                }
                let identity_pct = (matches as f32 / alen as f32) * 100.0;

                if identity_pct >= 55.0 && alen >= k {
                    hits.push(BlastHit {
                        doc_index: doc_idx,
                        doc_name: aln.name.clone(),
                        subject_name: seq.name.clone(),
                        subject_row: row_idx,
                        query_start: best_hsp.q_start,
                        query_end: best_hsp.q_end,
                        subject_start: best_hsp.s_start,
                        subject_end: best_hsp.s_end,
                        align_len: alen,
                        identity_pct,
                        bit_score,
                        e_value,
                        query_snippet: String::from_utf8_lossy(q_sub).into_owned(),
                        match_midline: midline,
                        subject_snippet: String::from_utf8_lossy(s_sub).into_owned(),
                    });
                }
            }
        }
    }

    hits.sort_by(|a, b| b.bit_score.total_cmp(&a.bit_score));
    hits.truncate(max_hits.max(1));
    hits
}

struct Hsp {
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    score: i32,
}

fn find_best_hsp(query: &[u8], subject: &[u8], k: usize, is_nuc: bool) -> Option<Hsp> {
    let match_sc = if is_nuc { 2 } else { 4 };
    let mismatch_sc = if is_nuc { -3 } else { -2 };
    let x_drop = 14i32;

    let mut best: Option<Hsp> = None;

    // Sample seeds across query to locate matching diagonals rapidly.
    let step = (query.len() / 28).max(1);
    let mut q_i = 0usize;
    while q_i + k <= query.len() {
        let seed = &query[q_i..(q_i + k)];
        for s_i in 0..=(subject.len() - k) {
            if &subject[s_i..(s_i + k)] != seed {
                continue;
            }
            // Extend right with X-drop.
            let mut best_r = k;
            let mut cur_r_score = (k as i32) * match_sc;
            let mut max_r_score = cur_r_score;
            let max_ext_r = (query.len() - q_i).min(subject.len() - s_i);
            for d in k..max_ext_r {
                cur_r_score += if query[q_i + d] == subject[s_i + d] {
                    match_sc
                } else {
                    mismatch_sc
                };
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
                cur_l_score += if query[q_i - d] == subject[s_i - d] {
                    match_sc
                } else {
                    mismatch_sc
                };
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::demo::build_default_workspace;

    #[test]
    fn local_blast_finds_exact_and_homologous_hits() {
        let docs = build_default_workspace();
        let subjects: Vec<(usize, &Alignment)> = docs
            .iter()
            .enumerate()
            .map(|(i, d)| (i, &d.alignment))
            .collect();
        let hits = run_local_blast("LSELHCDKLHVDPENFRLLGNVLV", &subjects, 10);
        assert!(!hits.is_empty());
        assert!(hits[0].identity_pct > 95.0);
    }
}
