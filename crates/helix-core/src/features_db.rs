//! Built-in Reference Features database (Helix Plasmid & Phage Features) and
//! automated sequence annotation scanner.

use crate::alignment::Alignment;
use crate::annotation::{Annotation, FeatureKind, UniProtLink};
use crate::sequence::{is_gap, reverse_complement, Strand};

/// A curated reference genetic element stored in the built-in Reference Features database.
#[derive(Debug, Clone)]
pub struct ReferenceFeature {
    pub name: &'static str,
    pub kind: FeatureKind,
    pub signature: &'static [u8],
    pub description: &'static str,
    pub uniprot_accession: Option<&'static str>,
}

/// Built-in Reference Feature library inspired by "Geneious Plasmid Features".
///
/// Includes classic plasmid origins, selection markers, promoters, tags, Gateway `att` sites,
/// BioBrick flanking scars, and bacteriophage cassettes.
pub const PLASMID_REFERENCE_FEATURES: &[ReferenceFeature] = &[
    ReferenceFeature {
        name: "pUC origin",
        kind: FeatureKind::Motif,
        signature: b"TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCC",
        description: "High-copy-number ColE1/pMB1/pBR322/pUC bacterial origin of replication",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "f1 origin",
        kind: FeatureKind::Motif,
        signature: b"ACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCC",
        description: "f1 bacteriophage single-stranded DNA packaging origin",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "Ampicillin (AmpR)",
        kind: FeatureKind::Gene,
        signature: b"ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAA",
        description: "Beta-lactamase resistance marker conferring ampicillin/carbenicillin resistance",
        uniprot_accession: Some("P62593"),
    },
    ReferenceFeature {
        name: "Kanamycin (KanR)",
        kind: FeatureKind::Gene,
        signature: b"ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCC",
        description: "Aminoglycoside phosphotransferase conferring kanamycin/neomycin resistance",
        uniprot_accession: Some("P00552"),
    },
    ReferenceFeature {
        name: "LacZ-alpha",
        kind: FeatureKind::Cds,
        signature: b"ATGACCATGATTACGCCAAGCTTGGTACCGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTGCTGGAATTCGCCCTT",
        description: "LacZ alpha-complementation peptide for blue/white colony screening",
        uniprot_accession: Some("P00722"),
    },
    ReferenceFeature {
        name: "T7 Promoter",
        kind: FeatureKind::Motif,
        signature: b"TAATACGACTCACTATAGGG",
        description: "Bacteriophage T7 RNA polymerase promoter",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "SP6 Promoter",
        kind: FeatureKind::Motif,
        signature: b"ATTTAGGTGACACTATAG",
        description: "Bacteriophage SP6 RNA polymerase promoter",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "T7 Terminator",
        kind: FeatureKind::Motif,
        signature: b"CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        description: "Stem-loop transcription terminator for T7 RNA polymerase",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "6xHis Tag",
        kind: FeatureKind::Domain,
        signature: b"CATCACCATCACCATCAC",
        description: "Polyhistidine metal-affinity purification tag",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "Gateway attB1",
        kind: FeatureKind::BindingSite,
        signature: b"ACAAGTTTGTACAAAAAAGCAGGCT",
        description: "Gateway site-specific recombination attB1 site",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "Gateway attB2",
        kind: FeatureKind::BindingSite,
        signature: b"ACCCAGCTTTCTTGTACAAAGTGGT",
        description: "Gateway site-specific recombination attB2 site",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "Gateway attP1",
        kind: FeatureKind::BindingSite,
        signature: b"AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAAC",
        description: "Gateway donor vector attP1 recombination site",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "BioBrick Prefix (RFC10)",
        kind: FeatureKind::Primer,
        signature: b"GAATTCGCGGCCGCTTCTAGAG",
        description: "BioBrick RFC 10 standard prefix (EcoRI-NotI-XbaI)",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "BioBrick Suffix (RFC10)",
        kind: FeatureKind::Primer,
        signature: b"TACTAGTAGCGGCCGCTGCAG",
        description: "BioBrick RFC 10 standard suffix (SpeI-NotI-PstI)",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "BioBrick Mixed Scar",
        kind: FeatureKind::Motif,
        signature: b"TACTAGAG",
        description: "Standard 8 bp XbaI/SpeI BioBrick ligation scar",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "Kozak Consensus",
        kind: FeatureKind::Motif,
        signature: b"GCCACCATGG",
        description: "Eukaryotic translation initiation Kozak consensus sequence",
        uniprot_accession: None,
    },
    ReferenceFeature {
        name: "Phage cos/ori",
        kind: FeatureKind::Motif,
        signature: b"AATTTAAATGGGCGGCGACCTTAATTAATTTATATGCTTAAGCGATCGATAATTTATAAA",
        description: "Bacteriophage AT-rich packaging site (cos) and origin of replication",
        uniprot_accession: None,
    },
];

/// Automatically scans the reference sequence of `alignment` against [`PLASMID_REFERENCE_FEATURES`]
/// (both forward and reverse strands, allowing up to `min_identity_pct` similarity) and adds any
/// newly discovered annotations. Returns the count of added annotations.
pub fn auto_annotate_from_reference_db(
    alignment: &mut Alignment,
    min_identity_pct: f32,
) -> usize {
    if !alignment.alphabet.is_nucleotide() || alignment.row_count() == 0 {
        return 0;
    }

    let ref_residues = &alignment.sequences()[0].residues;
    let mut ungapped: Vec<(usize, u8)> = Vec::with_capacity(ref_residues.len());
    for (col, &b) in ref_residues.iter().enumerate() {
        if !is_gap(b) {
            ungapped.push((col, b.to_ascii_uppercase()));
        }
    }
    let n = ungapped.len();
    if n < 8 {
        return 0;
    }

    let mut discovered: Vec<Annotation> = Vec::new();
    let min_frac = (min_identity_pct / 100.0).clamp(0.70, 1.0);

    for feat in PLASMID_REFERENCE_FEATURES {
        let sig = feat.signature;
        let slen = sig.len();
        if slen > n {
            continue;
        }
        let rc_sig = reverse_complement(sig);
        let max_mm = ((1.0 - min_frac) * slen as f32).floor() as usize;

        for (strand, pattern) in [(Strand::Forward, sig), (Strand::Reverse, rc_sig.as_slice())] {
            let mut i = 0usize;
            while i + slen <= n {
                let mut mm = 0usize;
                for k in 0..slen {
                    if ungapped[i + k].1 != pattern[k] {
                        mm += 1;
                        if mm > max_mm {
                            break;
                        }
                    }
                }
                if mm <= max_mm {
                    let col_start = ungapped[i].0;
                    let col_end = ungapped[i + slen - 1].0.saturating_add(1);
                    // Avoid duplicating an identical feature already present at the same locus.
                    let already_exists = alignment.annotations().iter().any(|existing| {
                        existing.name.eq_ignore_ascii_case(feat.name)
                            && existing.start.abs_diff(col_start) <= 3
                    }) || discovered.iter().any(|existing| {
                        existing.name.eq_ignore_ascii_case(feat.name)
                            && existing.start.abs_diff(col_start) <= 3
                    });

                    if !already_exists {
                        let id_pct = ((slen - mm) as f32 / slen as f32) * 100.0;
                        let mut ann = Annotation::new(0, feat.name, feat.kind, col_start, col_end)
                            .with_strand(strand)
                            .with_notes(format!(
                                "{} (Auto-annotated: {:.0}% match)",
                                feat.description, id_pct
                            ));
                        if let Some(acc) = feat.uniprot_accession {
                            ann.uniprot = Some(UniProtLink::new(acc));
                        }
                        discovered.push(ann);
                    }
                    i += slen;
                } else {
                    i += 1;
                }
            }
        }
    }

    let added = discovered.len();
    for ann in discovered {
        alignment.add_annotation(ann);
    }
    added
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{Alphabet, Sequence};

    #[test]
    fn auto_annotate_detects_t7_promoter_and_his_tag() {
        let seq = b"GGGGTAATACGACTCACTATAGGGGGGCCACCATGGCATCACCATCACCATCACTAA";
        let mut aln =
            Alignment::new("Construct", Alphabet::Dna, vec![Sequence::new("C1", seq)]).expect("ok");
        let added = auto_annotate_from_reference_db(&mut aln, 95.0);
        assert!(added >= 3);
        assert!(aln.annotations().iter().any(|a| a.name == "T7 Promoter"));
        assert!(aln.annotations().iter().any(|a| a.name == "Kozak Consensus"));
        assert!(aln.annotations().iter().any(|a| a.name == "6xHis Tag"));
    }
}
