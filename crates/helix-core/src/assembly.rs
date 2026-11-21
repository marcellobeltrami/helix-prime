//! Reference read mapping and greedy Overlap-Layout-Consensus (OLC) de novo contig assembly.

use crate::alignment::Alignment;
use crate::error::{HelixError, Result};
use crate::sequence::{reverse_complement, Alphabet, Sequence};

/// Maps a collection of query reads onto a reference sequence using k-mer seeding and
/// ungapped placement (checking both forward and reverse-complement orientations).
///
/// Produces a rectangular [`Alignment`] with the reference as row 0 and all mapped reads
/// padded cleanly to their mapped reference coordinates.
pub fn map_reads_to_reference(
    name: impl Into<String>,
    reference: &Sequence,
    reads: &[Sequence],
    alphabet: Alphabet,
) -> Result<(Alignment, usize)> {
    let ref_bytes = reference.ungapped_bytes();
    if ref_bytes.len() < 12 {
        return Err(HelixError::Alignment(
            "reference sequence is too short for read mapping (minimum 12 residues)".to_owned(),
        ));
    }

    let k = if ref_bytes.len() >= 100 { 8 } else { 5 };
    let mut aligned_rows = Vec::with_capacity(reads.len() + 1);
    aligned_rows.push(
        Sequence::new(format!("{} (Ref)", reference.name), &ref_bytes)
            .with_description(reference.description.clone()),
    );

    let mut mapped_count = 0usize;

    for read in reads {
        let raw = read.ungapped_bytes();
        if raw.len() < k {
            continue;
        }
        let fwd_best = best_placement(&ref_bytes, &raw, k);
        let rev_raw = if alphabet.is_nucleotide() {
            Some(reverse_complement(&raw))
        } else {
            None
        };
        let rev_best = rev_raw.as_deref().and_then(|r| best_placement(&ref_bytes, r, k));

        let chosen = match (fwd_best, rev_best) {
            (Some((p1, id1)), Some((p2, id2))) => {
                if id2 > id1 {
                    Some((p2, id2, true))
                } else {
                    Some((p1, id1, false))
                }
            }
            (Some((p1, id1)), None) => Some((p1, id1, false)),
            (None, Some((p2, id2))) => Some((p2, id2, true)),
            (None, None) => None,
        };

        if let Some((start_pos, identity, is_rev)) = chosen {
            if identity >= 0.65 {
                let placed_seq = if is_rev {
                    rev_raw.as_deref().unwrap_or(&raw)
                } else {
                    &raw
                };
                let mut padded = vec![b'-'; ref_bytes.len()];
                for (i, &b) in placed_seq.iter().enumerate() {
                    let target = start_pos + i;
                    if target < padded.len() {
                        padded[target] = b;
                    }
                }
                let label = if is_rev {
                    format!("{} (RC)", read.name)
                } else {
                    read.name.clone()
                };
                aligned_rows.push(Sequence::new(label, &padded).with_description(format!(
                    "Mapped at bp {}..{} ({:.0}% identity)",
                    start_pos + 1,
                    (start_pos + placed_seq.len()).min(ref_bytes.len()),
                    identity * 100.0
                )));
                mapped_count += 1;
            }
        }
    }

    let doc_name = name.into();
    let aln = Alignment::new(doc_name, alphabet, aligned_rows)?.with_description(format!(
        "Reference mapping: {mapped_count} of {} reads mapped onto {}",
        reads.len(),
        reference.name
    ));
    Ok((aln, mapped_count))
}

/// Performs greedy Overlap-Layout-Consensus (OLC) de novo contig assembly on a set of
/// overlapping sequences, returning an aligned contig view (consensus contig + contributing reads).
pub fn assemble_contigs_denovo(
    name: impl Into<String>,
    reads: &[Sequence],
    alphabet: Alphabet,
    min_overlap: usize,
) -> Result<Alignment> {
    if reads.len() < 2 {
        return Err(HelixError::Alignment(
            "de novo assembly requires at least 2 reads".to_owned(),
        ));
    }

    let min_ov = min_overlap.max(6);
    // Start with the longest read as the initial contig seed.
    let mut remaining: Vec<(String, Vec<u8>)> = reads
        .iter()
        .map(|s| (s.name.clone(), s.ungapped_bytes()))
        .filter(|(_, b)| !b.is_empty())
        .collect();
    remaining.sort_by_key(|(_, b)| std::cmp::Reverse(b.len()));

    let (seed_name, mut contig) = remaining.remove(0);
    // Track placed reads as (name, offset_relative_to_contig_0, bytes).
    let mut placed: Vec<(String, isize, Vec<u8>)> = vec![(seed_name, 0, contig.clone())];

    loop {
        let mut best_idx: Option<usize> = None;
        let mut best_offset: isize = 0;
        let mut best_ov_len = 0usize;
        let mut best_read_bytes: Vec<u8> = Vec::new();

        for (idx, (_, rbytes)) in remaining.iter().enumerate() {
            if let Some((off, ov)) = find_best_overlap(&contig, rbytes, min_ov) {
                if ov > best_ov_len {
                    best_ov_len = ov;
                    best_offset = off;
                    best_idx = Some(idx);
                    best_read_bytes = rbytes.clone();
                }
            }
            if alphabet.is_nucleotide() {
                let rc = reverse_complement(rbytes);
                if let Some((off, ov)) = find_best_overlap(&contig, &rc, min_ov) {
                    if ov > best_ov_len {
                        best_ov_len = ov;
                        best_offset = off;
                        best_idx = Some(idx);
                        best_read_bytes = rc;
                    }
                }
            }
        }

        let Some(chosen) = best_idx else {
            break;
        };
        let (rname, _) = remaining.remove(chosen);

        // Extend contig left or right as needed.
        if best_offset < 0 {
            let prepend = (-best_offset) as usize;
            let mut new_contig = Vec::with_capacity(prepend + contig.len());
            new_contig.extend_from_slice(&best_read_bytes[..prepend]);
            new_contig.extend_from_slice(&contig);
            contig = new_contig;
            for entry in &mut placed {
                entry.1 += prepend as isize;
            }
            placed.push((rname, 0, best_read_bytes));
        } else {
            let start = best_offset as usize;
            let end = start + best_read_bytes.len();
            if end > contig.len() {
                let extra_from = contig.len().saturating_sub(start);
                contig.extend_from_slice(&best_read_bytes[extra_from..]);
            }
            placed.push((rname, best_offset, best_read_bytes));
        }
    }

    let width = contig.len();
    let mut rows = Vec::with_capacity(placed.len() + 1);
    rows.push(
        Sequence::new("Contig_1 (Consensus)", &contig)
            .with_description(format!("De novo OLC contig ({} reads assembled)", placed.len())),
    );
    for (rname, offset, bytes) in placed {
        let mut padded = vec![b'-'; width];
        let start = offset.max(0) as usize;
        for (i, &b) in bytes.iter().enumerate() {
            if start + i < width {
                padded[start + i] = b;
            }
        }
        rows.push(Sequence::new(rname, &padded));
    }

    Alignment::new(name, alphabet, rows).map(|a| {
        a.with_description(format!(
            "De novo contig assembly ({} bp consensus)",
            width
        ))
    })
}

/// Finds the best start index of `read` inside `reference` using k-mer seeds.
fn best_placement(reference: &[u8], read: &[u8], k: usize) -> Option<(usize, f32)> {
    if read.len() > reference.len() || read.len() < k {
        return None;
    }

    // Collect candidate offsets from up to 4 evenly spaced k-mer seeds inside `read`.
    let mut candidates = Vec::new();
    let seed_positions = [0, read.len() / 3, (read.len() * 2) / 3, read.len() - k];
    for &sp in &seed_positions {
        let seed = &read[sp..(sp + k)];
        for ref_i in 0..=(reference.len() - k) {
            if &reference[ref_i..(ref_i + k)] == seed && ref_i >= sp {
                let start = ref_i - sp;
                if start + read.len() <= reference.len() && !candidates.contains(&start) {
                    candidates.push(start);
                }
            }
        }
    }

    let mut best: Option<(usize, f32)> = None;
    for start in candidates {
        let matches = reference[start..(start + read.len())]
            .iter()
            .zip(read.iter())
            .filter(|(a, b)| a == b)
            .count();
        let id = matches as f32 / read.len() as f32;
        if best.is_none_or(|(_, b_id)| id > b_id) {
            best = Some((start, id));
        }
    }
    best
}

/// Finds the longest suffix/prefix or containment overlap between `contig` and `read` with >= 88% identity.
fn find_best_overlap(contig: &[u8], read: &[u8], min_ov: usize) -> Option<(isize, usize)> {
    let clen = contig.len() as isize;
    let rlen = read.len() as isize;
    let mut best: Option<(isize, usize)> = None;

    for offset in -(rlen - min_ov as isize)..=(clen - min_ov as isize) {
        let c_start = offset.max(0) as usize;
        let r_start = (-offset).max(0) as usize;
        let ov_len = (contig.len() - c_start).min(read.len() - r_start);
        if ov_len < min_ov {
            continue;
        }
        let matches = contig[c_start..(c_start + ov_len)]
            .iter()
            .zip(&read[r_start..(r_start + ov_len)])
            .filter(|(a, b)| a == b)
            .count();
        let id = matches as f32 / ov_len as f32;
        if id >= 0.88 && best.is_none_or(|(_, b_ov)| ov_len > b_ov) {
            best = Some((offset, ov_len));
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn denovo_assembly_merges_overlapping_reads_into_single_contig() {
        let reads = vec![
            Sequence::new("R1", b"ATGGTGCATCTGACTCCTGAGGAG"),
            Sequence::new("R2", b"CTGACTCCTGAGGAGAAGTCTGCCGTT"),
            Sequence::new("R3", b"GAGAAGTCTGCCGTTACTGCCCTGT"),
        ];
        let contig = assemble_contigs_denovo("ContigTest", &reads, Alphabet::Dna, 10)
            .expect("assemble");
        assert_eq!(contig.row_count(), 4);
        assert!(contig.width() >= 36);
    }

    #[test]
    fn map_reads_aligns_forward_and_reverse_complement_reads() {
        let ref_seq = Sequence::new("Ref", b"ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC");
        let reads = vec![
            Sequence::new("Read_Fwd", b"CTGACTCCTGAGGAGAAGTC"),
            Sequence::new("Read_Rev", &reverse_complement(b"TGCCCTGTGGGGCAAGGTGA")),
        ];
        let (aln, mapped) =
            map_reads_to_reference("MapTest", &ref_seq, &reads, Alphabet::Dna).expect("map");
        assert_eq!(mapped, 2);
        assert_eq!(aln.row_count(), 3);
    }
}
