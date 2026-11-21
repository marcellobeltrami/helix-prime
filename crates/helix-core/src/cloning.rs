//! In-silico molecular cloning simulation: Restriction Enzyme, TOPO-TA, Gateway, and BioBrick
//! assembly with automatic annotation transfer and construct lineage tracking.

use crate::alignment::Alignment;
use crate::annotation::{Annotation, FeatureKind};
use crate::error::{HelixError, Result};
use crate::restriction::{find_enzyme_by_name, find_restriction_sites};
use crate::sequence::{reverse_complement, Alphabet, Sequence, Strand};
use serde::{Deserialize, Serialize};

/// Supported in-silico cloning methodologies matching Geneious Prime's Cloning menu.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum CloningMethod {
    /// Traditional Type II restriction digest & ligation (1 or 2 enzymes).
    #[default]
    RestrictionDigest,
    /// TOPO-TA linearized vector insertion with 3'-A overhang ligation.
    TopoTa,
    /// Site-specific Gateway BP/LR recombination (`attB`/`attP`/`attL`/`attR`).
    Gateway,
    /// BioBrick RFC 10 standard assembly (creates 8 bp `TACTAGAG` mixed scar).
    BioBrick,
}

impl CloningMethod {
    pub const ALL: [Self; 4] = [
        Self::RestrictionDigest,
        Self::TopoTa,
        Self::Gateway,
        Self::BioBrick,
    ];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::RestrictionDigest => "Restriction Enzyme Cloning",
            Self::TopoTa => "TOPO-TA Cloning",
            Self::Gateway => "Gateway Recombination",
            Self::BioBrick => "BioBrick (RFC 10) Assembly",
        }
    }
}

/// Parameters controlling an in-silico cloning simulation.
#[derive(Debug, Clone)]
pub struct CloningRequest {
    pub method: CloningMethod,
    pub construct_name: String,
    /// Primary restriction enzyme name (used by RestrictionDigest).
    pub enzyme_5prime: String,
    /// Optional secondary enzyme name (if empty or equal to 5', single-enzyme digest).
    pub enzyme_3prime: String,
    /// When `true`, inserts the reverse complement of the insert fragment.
    pub reverse_insert: bool,
}

/// Outcome of an in-silico cloning simulation ready to be added to the workspace.
#[derive(Debug, Clone)]
pub struct CloningOutcome {
    pub construct: Alignment,
    pub summary: String,
    pub insertion_site: usize,
    pub insert_len: usize,
}

/// Simulates cloning `insert_doc` into `vector_doc` using `req`, transferring all intact
/// vector annotations and insert annotations onto the resulting circular construct.
pub fn simulate_cloning(
    vector_doc: &Alignment,
    insert_doc: &Alignment,
    req: &CloningRequest,
) -> Result<CloningOutcome> {
    if !vector_doc.alphabet.is_nucleotide() || !insert_doc.alphabet.is_nucleotide() {
        return Err(HelixError::Alignment(
            "cloning requires DNA/RNA vector and insert sequences".to_owned(),
        ));
    }

    let vector_seq = vector_doc
        .sequences()
        .first()
        .ok_or_else(|| HelixError::Alignment("vector document has no sequences".to_owned()))?
        .ungapped_bytes();
    let mut insert_seq = insert_doc
        .sequences()
        .first()
        .ok_or_else(|| HelixError::Alignment("insert document has no sequences".to_owned()))?
        .ungapped_bytes();

    if vector_seq.is_empty() || insert_seq.is_empty() {
        return Err(HelixError::Alignment(
            "vector and insert sequences must be non-empty".to_owned(),
        ));
    }

    if req.reverse_insert {
        insert_seq = reverse_complement(&insert_seq);
    }

    let (backbone_left, backbone_right, insert_payload, scar_note) = match req.method {
        CloningMethod::RestrictionDigest => {
            let enz1 = find_enzyme_by_name(&req.enzyme_5prime).ok_or_else(|| {
                HelixError::Alignment(format!("unknown restriction enzyme '{}'", req.enzyme_5prime))
            })?;
            let sites1 = find_restriction_sites(&vector_seq, &[*enz1], true);
            let cut1 = sites1
                .first()
                .map(|h| h.cut_col)
                .ok_or_else(|| {
                    HelixError::Alignment(format!(
                        "vector '{}' contains no {} site",
                        vector_doc.name, enz1.name
                    ))
                })?;

            let enz2_opt = if req.enzyme_3prime.trim().is_empty()
                || req.enzyme_3prime.eq_ignore_ascii_case(enz1.name)
            {
                None
            } else {
                Some(find_enzyme_by_name(&req.enzyme_3prime).ok_or_else(|| {
                    HelixError::Alignment(format!(
                        "unknown second restriction enzyme '{}'",
                        req.enzyme_3prime
                    ))
                })?)
            };

            let (drop_start, drop_end, desc) = if let Some(enz2) = enz2_opt {
                let sites2 = find_restriction_sites(&vector_seq, &[*enz2], true);
                let cut2 = sites2.first().map(|h| h.cut_col).ok_or_else(|| {
                    HelixError::Alignment(format!(
                        "vector '{}' contains no {} site",
                        vector_doc.name, enz2.name
                    ))
                })?;
                let (lo, hi) = if cut1 <= cut2 {
                    (cut1, cut2)
                } else {
                    (cut2, cut1)
                };
                (
                    lo,
                    hi,
                    format!(
                        "Double digest ({} / {}) replacing {} bp polylinker with {} bp insert",
                        enz1.name,
                        enz2.name,
                        hi.saturating_sub(lo),
                        insert_seq.len()
                    ),
                )
            } else {
                (
                    cut1,
                    cut1,
                    format!(
                        "Single-site {} linearization at bp {} + {} bp insert ligation",
                        enz1.name,
                        cut1 + 1,
                        insert_seq.len()
                    ),
                )
            };

            (
                vector_seq[..drop_start].to_vec(),
                vector_seq[drop_end..].to_vec(),
                insert_seq,
                desc,
            )
        }
        CloningMethod::TopoTa => {
            // Locate EcoRI / TOPO cloning site in LacZ-alpha if present, otherwise insert near
            // the LacZ-alpha annotation midpoint or 15% around the vector.
            let topo_cut = find_enzyme_by_name("EcoRI")
                .and_then(|e| find_restriction_sites(&vector_seq, &[*e], true).first().map(|h| h.cut_col))
                .or_else(|| {
                    vector_doc
                        .annotations()
                        .iter()
                        .find(|a| a.name.to_ascii_lowercase().contains("lacz"))
                        .map(|a| (a.start + a.end) / 2)
                })
                .unwrap_or(vector_seq.len() / 6);

            let mut payload = Vec::with_capacity(insert_seq.len() + 2);
            payload.push(b'A'); // 3'-A Taq polymerase overhang
            payload.extend_from_slice(&insert_seq);
            payload.push(b'T');

            (
                vector_seq[..topo_cut].to_vec(),
                vector_seq[topo_cut..].to_vec(),
                payload,
                format!(
                    "TOPO-TA isomerase ligation at bp {} ({} bp PCR product + A/T overhangs)",
                    topo_cut + 1,
                    insert_seq.len()
                ),
            )
        }
        CloningMethod::Gateway => {
            // Recombine at attP1/attR1/attB1 site if present, replacing any ccdB cassette between att sites.
            let att_hits: Vec<(usize, usize)> = vector_doc
                .annotations()
                .iter()
                .filter(|a| a.name.to_ascii_lowercase().contains("att"))
                .map(|a| (a.start, a.end))
                .collect();

            let (site_start, site_end) = if att_hits.len() >= 2 {
                (att_hits[0].1.min(vector_seq.len()), att_hits[1].0.min(vector_seq.len()))
            } else if let Some(&(s, e)) = att_hits.first() {
                (s.min(vector_seq.len()), e.min(vector_seq.len()))
            } else {
                let mid = vector_seq.len() / 4;
                (mid, mid)
            };

            let (lo, hi) = if site_start <= site_end {
                (site_start, site_end)
            } else {
                (site_end, site_start)
            };

            // Add canonical attB1/attB2 recombination junctions around the inserted gene.
            let mut payload = Vec::with_capacity(insert_seq.len() + 50);
            payload.extend_from_slice(b"ACAAGTTTGTACAAAAAAGCAGGCT");
            payload.extend_from_slice(&insert_seq);
            payload.extend_from_slice(b"ACCCAGCTTTCTTGTACAAAGTGGT");

            (
                vector_seq[..lo].to_vec(),
                vector_seq[hi..].to_vec(),
                payload,
                format!(
                    "Gateway LR/BP recombination replacing bp {}..{} with attB-flanked {} bp insert",
                    lo + 1,
                    hi,
                    insert_seq.len()
                ),
            )
        }
        CloningMethod::BioBrick => {
            // BioBrick RFC 10: ligate downstream of SpeI/XbaI creating the canonical 8 bp scar `TACTAGAG`.
            let cut = find_enzyme_by_name("SpeI")
                .or_else(|| find_enzyme_by_name("XbaI"))
                .and_then(|e| find_restriction_sites(&vector_seq, &[*e], true).first().map(|h| h.cut_col))
                .unwrap_or(vector_seq.len() / 5);

            let mut payload = Vec::with_capacity(insert_seq.len() + 16);
            payload.extend_from_slice(b"TACTAGAG"); // RFC 10 mixed scar
            payload.extend_from_slice(&insert_seq);
            payload.extend_from_slice(b"ACTAGT"); // Regenerated downstream SpeI site

            (
                vector_seq[..cut].to_vec(),
                vector_seq[cut..].to_vec(),
                payload,
                format!(
                    "BioBrick RFC 10 assembly at bp {} with TACTAGAG mixed scar + {} bp part",
                    cut + 1,
                    insert_seq.len()
                ),
            )
        }
    };

    let insertion_site = backbone_left.len();
    let dropped_len = vector_seq
        .len()
        .saturating_sub(backbone_left.len() + backbone_right.len());
    let insert_len = insert_payload.len();

    let mut combined = Vec::with_capacity(backbone_left.len() + insert_len + backbone_right.len());
    combined.extend_from_slice(&backbone_left);
    combined.extend_from_slice(&insert_payload);
    combined.extend_from_slice(&backbone_right);

    let name = if req.construct_name.trim().is_empty() {
        format!("{}::{}", vector_doc.name, insert_doc.name)
    } else {
        req.construct_name.trim().to_owned()
    };

    let mut construct = Alignment::new(
        name.clone(),
        Alphabet::Dna,
        vec![Sequence::new(name, &combined).with_description(scar_note.clone())],
    )?
    .with_description(scar_note.clone());

    // Transfer intact vector annotations adjusted for the inserted payload offset.
    let drop_end = insertion_site + dropped_len;
    let delta: isize = insert_len as isize - dropped_len as isize;
    for ann in vector_doc.annotations() {
        if ann.end <= insertion_site {
            construct.add_annotation(ann.clone());
        } else if ann.start >= drop_end {
            let mut shifted = ann.clone();
            shifted.start = (ann.start as isize + delta).max(0) as usize;
            shifted.end = (ann.end as isize + delta).max(shifted.start as isize + 1) as usize;
            construct.add_annotation(shifted);
        }
    }

    // Add an explicit Insert annotation and transfer inner insert annotations.
    construct.add_annotation(
        Annotation::new(
            0,
            format!("Insert: {}", insert_doc.name),
            FeatureKind::Cds,
            insertion_site,
            insertion_site + insert_len,
        )
        .with_strand(if req.reverse_insert {
            Strand::Reverse
        } else {
            Strand::Forward
        })
        .with_notes(format!("Cloned via {}: {}", req.method.label(), scar_note)),
    );

    if !req.reverse_insert {
        for ann in insert_doc.annotations() {
            let mut moved = ann.clone();
            moved.start = insertion_site + ann.start;
            moved.end = (insertion_site + ann.end).min(combined.len());
            if moved.start < moved.end {
                construct.add_annotation(moved);
            }
        }
    }

    Ok(CloningOutcome {
        construct,
        summary: scar_note,
        insertion_site,
        insert_len,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restriction_and_biobrick_cloning_produce_valid_constructs() {
        let vec_aln = Alignment::new(
            "pTestVec",
            Alphabet::Dna,
            vec![Sequence::new("pTestVec", b"AAAAAGAATTCGGGGGGATCCTTTTT")],
        )
        .expect("vec");
        let ins_aln = Alignment::new(
            "GFP_mini",
            Alphabet::Dna,
            vec![Sequence::new("GFP_mini", b"ATGGTGAGCAAGGGCGAGGAG")],
        )
        .expect("ins");

        let res = simulate_cloning(
            &vec_aln,
            &ins_aln,
            &CloningRequest {
                method: CloningMethod::RestrictionDigest,
                construct_name: String::from("pTestVec-GFP"),
                enzyme_5prime: String::from("EcoRI"),
                enzyme_3prime: String::from("BamHI"),
                reverse_insert: false,
            },
        )
        .expect("cloning");
        assert!(res.construct.width() > ins_aln.width());
        assert!(res.construct.annotations().iter().any(|a| a.name.contains("GFP_mini")));

        let bb = simulate_cloning(
            &vec_aln,
            &ins_aln,
            &CloningRequest {
                method: CloningMethod::BioBrick,
                construct_name: String::from("pBB-GFP"),
                enzyme_5prime: String::new(),
                enzyme_3prime: String::new(),
                reverse_insert: false,
            },
        )
        .expect("biobrick");
        let seq = &bb.construct.sequences()[0].residues;
        assert!(seq.windows(8).any(|w| w == b"TACTAGAG"));
    }
}
