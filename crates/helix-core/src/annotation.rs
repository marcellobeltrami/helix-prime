//! Sequence & alignment annotations with UniProt/AlphaFold cross-references and track packing.

use crate::sequence::Strand;
use serde::{Deserialize, Serialize};

/// Biological classification of an annotation feature, matching GenBank/GFF3/UniProt conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FeatureKind {
    Gene,
    Cds,
    Domain,
    ActiveSite,
    BindingSite,
    Variant,
    Motif,
    Primer,
    Helix,
    BetaStrand,
    Misc,
}

impl FeatureKind {
    /// Human-readable label for the GUI inspector and annotation badges.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Gene => "Gene",
            Self::Cds => "CDS",
            Self::Domain => "Domain",
            Self::ActiveSite => "Active Site",
            Self::BindingSite => "Binding Site",
            Self::Variant => "Variant",
            Self::Motif => "Motif",
            Self::Primer => "Primer",
            Self::Helix => "Alpha Helix",
            Self::BetaStrand => "Beta Strand",
            Self::Misc => "Misc Feature",
        }
    }

    /// Parses a GenBank, GFF3, or UniProt feature type string into a [`FeatureKind`].
    #[must_use]
    pub fn from_type_str(raw: &str) -> Self {
        let lower = raw.trim().to_ascii_lowercase();
        match lower.as_str() {
            "gene" => Self::Gene,
            "cds" | "coding_sequence" => Self::Cds,
            "domain" | "region" | "repeat" | "zinc finger" | "dna binding" | "transmembrane"
            | "topological domain" | "coiled coil" => Self::Domain,
            "active site" | "active_site" => Self::ActiveSite,
            "binding site" | "binding_site" | "metal binding" | "nucleotide binding"
            | "site" => Self::BindingSite,
            "variant" | "natural variant" | "mutagenesis" | "mutation" | "variation"
            | "misc_difference" => Self::Variant,
            "motif" | "promoter" | "regulatory" => Self::Motif,
            "primer" | "primer_bind" => Self::Primer,
            "helix" => Self::Helix,
            "beta strand" | "strand" | "sheet" => Self::BetaStrand,
            _ => Self::Misc,
        }
    }

    /// Default RGB color palette inspired by Geneious Prime's annotation tracks.
    #[must_use]
    pub const fn default_rgb(self) -> [u8; 3] {
        match self {
            Self::Gene => [66, 165, 118],
            Self::Cds => [242, 194, 58],
            Self::Domain => [72, 138, 224],
            Self::ActiveSite => [230, 82, 82],
            Self::BindingSite => [235, 134, 52],
            Self::Variant => [198, 78, 178],
            Self::Motif => [56, 178, 172],
            Self::Primer => [102, 187, 106],
            Self::Helix => [224, 92, 128],
            Self::BetaStrand => [235, 182, 52],
            Self::Misc => [138, 149, 165],
        }
    }
}

/// Link connecting an annotation to a UniProtKB protein entry and AlphaFold 3D model.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct UniProtLink {
    /// UniProtKB primary accession (e.g., `"P04637"` for human TP53).
    pub accession: String,
    /// Optional UniProt entry name (e.g., `"P53_HUMAN"`).
    pub entry_name: Option<String>,
    /// Recommended protein name from UniProt.
    pub protein_name: Option<String>,
    /// Organism scientific name.
    pub organism: Option<String>,
    /// 1-based residue number in the UniProt canonical sequence corresponding to
    /// the start of this annotation (default: `1`).
    pub uniprot_start_residue: usize,
    /// Preferred chain identifier when displaying multi-chain structures (default: `'A'`).
    pub chain_id: char,
}

impl UniProtLink {
    /// Creates a minimal UniProt link from an accession string.
    #[must_use]
    pub fn new(accession: impl Into<String>) -> Self {
        Self {
            accession: accession.into().trim().to_ascii_uppercase(),
            entry_name: None,
            protein_name: None,
            organism: None,
            uniprot_start_residue: 1,
            chain_id: 'A',
        }
    }

    /// Builder helper to set the starting 1-based UniProt residue index.
    #[must_use]
    pub const fn with_start_residue(mut self, residue_1based: usize) -> Self {
        self.uniprot_start_residue = if residue_1based == 0 {
            1
        } else {
            residue_1based
        };
        self
    }

    /// Builder helper to attach entry/protein metadata.
    #[must_use]
    pub fn with_metadata(
        mut self,
        entry_name: impl Into<String>,
        protein_name: impl Into<String>,
        organism: impl Into<String>,
    ) -> Self {
        self.entry_name = Some(entry_name.into());
        self.protein_name = Some(protein_name.into());
        self.organism = Some(organism.into());
        self
    }
}

/// A biological annotation anchored to a column span `[start, end)` in an alignment.
///
/// # Invariants
/// - `start < end` (half-open 0-based alignment column interval).
/// - If `target_sequence` is `None`, the annotation applies to the consensus / whole alignment;
///   otherwise it references a specific sequence row index.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Annotation {
    /// Unique numeric identifier within its owning alignment.
    pub id: u64,
    /// Display name drawn inside or beside the annotation arrow.
    pub name: String,
    /// Biological feature category.
    pub kind: FeatureKind,
    /// 0-based inclusive start column in alignment coordinates.
    pub start: usize,
    /// 0-based exclusive end column in alignment coordinates (`start < end`).
    pub end: usize,
    /// Strand orientation (`+` renders a right-pointing arrow, `-` left-pointing).
    pub strand: Strand,
    /// Optional row index if this annotation belongs to a specific sequence.
    pub target_sequence: Option<usize>,
    /// Optional RGB override; falls back to [`FeatureKind::default_rgb`] when `None`.
    pub color_override: Option<[u8; 3]>,
    /// Free-form notes or functional description.
    pub notes: String,
    /// Optional UniProt / AlphaFold 3D protein link.
    pub uniprot: Option<UniProtLink>,
    /// Additional GenBank/GFF3 key-value qualifiers.
    pub qualifiers: Vec<(String, String)>,
}

impl Annotation {
    /// Creates a validated annotation, clamping `end` to be at least `start + 1`.
    #[must_use]
    pub fn new(
        id: u64,
        name: impl Into<String>,
        kind: FeatureKind,
        start: usize,
        end: usize,
    ) -> Self {
        let safe_end = if end <= start {
            start.saturating_add(1)
        } else {
            end
        };
        Self {
            id,
            name: name.into(),
            kind,
            start,
            end: safe_end,
            strand: Strand::Forward,
            target_sequence: None,
            color_override: None,
            notes: String::new(),
            uniprot: None,
            qualifiers: Vec::new(),
        }
    }

    /// Length of the annotation in alignment columns.
    #[inline]
    #[must_use]
    pub const fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// Returns `true` if `len() == 0` (always false for valid annotations).
    #[inline]
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns `true` if the alignment column `col` lies inside `[start, end)`.
    #[inline]
    #[must_use]
    pub const fn contains_column(&self, col: usize) -> bool {
        col >= self.start && col < self.end
    }

    /// Effective RGB color for rendering this annotation track item.
    #[inline]
    #[must_use]
    pub const fn rgb(&self) -> [u8; 3] {
        match self.color_override {
            Some(rgb) => rgb,
            None => self.kind.default_rgb(),
        }
    }

    /// Builder helper to set strand orientation.
    #[must_use]
    pub const fn with_strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Builder helper to set target sequence row index.
    #[must_use]
    pub const fn with_target_sequence(mut self, row: usize) -> Self {
        self.target_sequence = Some(row);
        self
    }

    /// Builder helper to set notes text.
    #[must_use]
    pub fn with_notes(mut self, notes: impl Into<String>) -> Self {
        self.notes = notes.into();
        self
    }

    /// Builder helper to attach a UniProt/AlphaFold link.
    #[must_use]
    pub fn with_uniprot(mut self, link: UniProtLink) -> Self {
        self.uniprot = Some(link);
        self
    }
}

/// Packs annotations into non-overlapping horizontal visual lanes (tracks) using a greedy
/// interval scheduling sweep.
///
/// Returns `(lane_count, lane_assignment)` where `lane_assignment[i]` is the 0-based track lane
/// for `annotations[i]`. A small `padding_cols` gap prevents adjacent arrow tips from touching.
#[must_use]
pub fn pack_annotation_lanes(
    annotations: &[Annotation],
    padding_cols: usize,
) -> (usize, Vec<usize>) {
    if annotations.is_empty() {
        return (0, Vec::new());
    }

    let mut order: Vec<usize> = (0..annotations.len()).collect();
    order.sort_by_key(|&idx| (annotations[idx].start, std::cmp::Reverse(annotations[idx].end)));

    let mut lane_ends: Vec<usize> = Vec::new();
    let mut assignments = vec![0usize; annotations.len()];

    for idx in order {
        let ann = &annotations[idx];
        let mut chosen_lane = None;
        for (lane_idx, lane_end) in lane_ends.iter_mut().enumerate() {
            if ann.start >= *lane_end {
                *lane_end = ann.end.saturating_add(padding_cols);
                chosen_lane = Some(lane_idx);
                break;
            }
        }
        let lane = match chosen_lane {
            Some(l) => l,
            None => {
                lane_ends.push(ann.end.saturating_add(padding_cols));
                lane_ends.len() - 1
            }
        };
        assignments[idx] = lane;
    }

    (lane_ends.len(), assignments)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lane_packing_separates_overlapping_features() {
        let anns = vec![
            Annotation::new(1, "Full", FeatureKind::Gene, 0, 50),
            Annotation::new(2, "DomA", FeatureKind::Domain, 5, 20),
            Annotation::new(3, "DomB", FeatureKind::Domain, 22, 40),
            Annotation::new(4, "Site", FeatureKind::ActiveSite, 10, 11),
        ];
        let (lanes, assigned) = pack_annotation_lanes(&anns, 1);
        assert_eq!(lanes, 3);
        // DomA and DomB don't overlap, so they share lane 1.
        assert_eq!(assigned[1], assigned[2]);
        assert_ne!(assigned[0], assigned[1]);
        assert_ne!(assigned[1], assigned[3]);
    }
}
