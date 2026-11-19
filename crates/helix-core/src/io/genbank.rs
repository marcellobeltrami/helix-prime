//! GenBank (`.gb` / `.gbk`) parser supporting sequence extraction, feature annotations,
//! strand/complement locations, and `/db_xref="UniProtKB:..."` cross-references.

use crate::alignment::Alignment;
use crate::annotation::{Annotation, FeatureKind, UniProtLink};
use crate::error::{HelixError, Result};
use crate::sequence::{Alphabet, Sequence, Strand};

/// Parses a GenBank flatfile into an [`Alignment`] populated with all extracted annotations.
pub fn parse_genbank(fallback_name: impl Into<String>, input: &str) -> Result<Alignment> {
    let mut locus_name = fallback_name.into();
    let mut definition = String::new();
    let mut organism: Option<String> = None;
    let mut seq_bytes: Vec<u8> = Vec::new();
    let mut raw_features: Vec<RawFeature> = Vec::new();

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Section {
        Header,
        Features,
        Origin,
    }

    let mut section = Section::Header;
    let mut current_feature: Option<RawFeature> = None;

    for raw_line in input.lines() {
        let trimmed = raw_line.trim();
        if trimmed == "//" {
            break;
        }

        if raw_line.starts_with("LOCUS") {
            if let Some(token) = trimmed.split_whitespace().nth(1) {
                locus_name = token.to_owned();
            }
            section = Section::Header;
            continue;
        }
        if raw_line.starts_with("DEFINITION") {
            definition = raw_line
                .get(10..)
                .unwrap_or("")
                .trim()
                .trim_end_matches('.')
                .to_owned();
            continue;
        }
        if trimmed.starts_with("ORGANISM") {
            let org = trimmed.trim_start_matches("ORGANISM").trim();
            if !org.is_empty() {
                organism = Some(org.to_owned());
            }
            continue;
        }
        if raw_line.starts_with("FEATURES") {
            section = Section::Features;
            continue;
        }
        if raw_line.starts_with("ORIGIN") {
            if let Some(feat) = current_feature.take() {
                raw_features.push(feat);
            }
            section = Section::Origin;
            continue;
        }

        match section {
            Section::Header => {}
            Section::Features => {
                // Standard GenBank feature keys start at column 5 (`     CDS             1..300`),
                // whereas qualifiers start at column 21 (`                     /gene="HBB"`).
                let indent = raw_line.bytes().take_while(|&b| b == b' ').count();
                if (2..18).contains(&indent) && !trimmed.starts_with('/') {
                    if let Some(feat) = current_feature.take() {
                        raw_features.push(feat);
                    }
                    let mut parts = trimmed.split_whitespace();
                    if let (Some(key), Some(loc)) = (parts.next(), parts.next()) {
                        current_feature = Some(RawFeature {
                            key: key.to_owned(),
                            location: loc.to_owned(),
                            qualifiers: Vec::new(),
                        });
                    }
                } else if let Some(feat) = current_feature.as_mut() {
                    if let Some(stripped) = trimmed.strip_prefix('/') {
                        if let Some((k, v)) = stripped.split_once('=') {
                            let clean_v = v.trim().trim_matches('"').to_owned();
                            feat.qualifiers.push((k.trim().to_owned(), clean_v));
                        } else {
                            feat.qualifiers
                                .push((stripped.trim().to_owned(), String::new()));
                        }
                    } else if let Some(last) = feat.qualifiers.last_mut() {
                        // Multi-line qualifier continuation.
                        let cont = trimmed.trim_end_matches('"');
                        if !cont.is_empty() {
                            if !last.1.is_empty() {
                                last.1.push(' ');
                            }
                            last.1.push_str(cont);
                        }
                    } else {
                        feat.location.push_str(trimmed);
                    }
                }
            }
            Section::Origin => {
                for b in raw_line.bytes() {
                    if b.is_ascii_alphabetic() || b == b'-' {
                        seq_bytes.push(b.to_ascii_uppercase());
                    }
                }
            }
        }
    }

    if let Some(feat) = current_feature.take() {
        raw_features.push(feat);
    }

    if seq_bytes.is_empty() {
        return Err(HelixError::parse(
            "GenBank",
            None,
            "no sequence residues found under ORIGIN section",
        ));
    }

    let alphabet = Alphabet::infer_from_slices([seq_bytes.as_slice()]);
    let mut seq = Sequence::new(locus_name.clone(), &seq_bytes).with_description(definition.clone());
    seq.organism = organism.clone();

    let mut aln = Alignment::new(locus_name, alphabet, vec![seq])?.with_description(definition);

    for feat in raw_features {
        if feat.key.eq_ignore_ascii_case("source") {
            continue;
        }
        if let Some((start, end, strand)) = parse_genbank_location(&feat.location) {
            let kind = FeatureKind::from_type_str(&feat.key);
            let label = feat
                .qualifier("label")
                .or_else(|| feat.qualifier("gene"))
                .or_else(|| feat.qualifier("product"))
                .unwrap_or(&feat.key)
                .to_owned();

            let mut ann = Annotation::new(0, label, kind, start, end)
                .with_strand(strand)
                .with_target_sequence(0);

            if let Some(note) = feat.qualifier("note").or_else(|| feat.qualifier("product")) {
                ann.notes = note.to_owned();
            }

            // Detect UniProtKB db_xref qualifiers automatically.
            for (k, v) in &feat.qualifiers {
                if k == "db_xref" {
                    if let Some(acc) = extract_uniprot_db_xref(v) {
                        let mut link = UniProtLink::new(acc);
                        if let Some(prod) = feat.qualifier("product") {
                            link.protein_name = Some(prod.to_owned());
                        }
                        link.organism = organism.clone();
                        ann.uniprot = Some(link);
                    }
                }
            }

            ann.qualifiers = feat.qualifiers;
            aln.add_annotation(ann);
        }
    }

    Ok(aln)
}

struct RawFeature {
    key: String,
    location: String,
    qualifiers: Vec<(String, String)>,
}

impl RawFeature {
    fn qualifier(&self, key: &str) -> Option<&str> {
        self.qualifiers
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(key))
            .map(|(_, v)| v.as_str())
    }
}

/// Parses GenBank location expressions such as `1..444`, `complement(10..200)`, or
/// `join(1..50,100..250)` into a 0-based half-open span `(start, end, strand)`.
fn parse_genbank_location(raw: &str) -> Option<(usize, usize, Strand)> {
    let strand = if raw.contains("complement(") {
        Strand::Reverse
    } else {
        Strand::Forward
    };

    // Extract all positive integer coordinates in the expression and take the outer span.
    let mut min_coord: Option<usize> = None;
    let mut max_coord: Option<usize> = None;
    let mut num_buf = String::new();

    for ch in raw.chars() {
        if ch.is_ascii_digit() {
            num_buf.push(ch);
        } else if !num_buf.is_empty() {
            if let Ok(val) = num_buf.parse::<usize>() {
                min_coord = Some(min_coord.map_or(val, |m| m.min(val)));
                max_coord = Some(max_coord.map_or(val, |m| m.max(val)));
            }
            num_buf.clear();
        }
    }
    if !num_buf.is_empty() {
        if let Ok(val) = num_buf.parse::<usize>() {
            min_coord = Some(min_coord.map_or(val, |m| m.min(val)));
            max_coord = Some(max_coord.map_or(val, |m| m.max(val)));
        }
    }

    let start_1 = min_coord?;
    let end_1 = max_coord?;
    let start_0 = start_1.saturating_sub(1);
    let end_0 = end_1.max(start_0 + 1);
    Some((start_0, end_0, strand))
}

fn extract_uniprot_db_xref(val: &str) -> Option<&str> {
    let lower = val.to_ascii_lowercase();
    if lower.starts_with("uniprot") {
        val.split(':').nth(1).map(str::trim).filter(|s| !s.is_empty())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genbank_parser_extracts_sequence_features_and_uniprot_xref() {
        let gb = "\
LOCUS       HBB_MINI                  27 bp    DNA     linear   PRI 01-JAN-2026
DEFINITION  Human hemoglobin beta mini CDS.
  ORGANISM  Homo sapiens
FEATURES             Location/Qualifiers
     CDS             1..27
                     /gene=\"HBB\"
                     /product=\"Hemoglobin subunit beta\"
                     /db_xref=\"UniProtKB/Swiss-Prot:P68871\"
ORIGIN
        1 atggtgcacc tgactcctga ggagaag
//
";
        let aln = parse_genbank("fallback", gb).expect("parse genbank");
        assert_eq!(aln.name, "HBB_MINI");
        assert_eq!(aln.width(), 27);
        assert_eq!(aln.annotations().len(), 1);
        let ann = &aln.annotations()[0];
        assert_eq!(ann.name, "HBB");
        assert_eq!(ann.kind, FeatureKind::Cds);
        assert_eq!(
            ann.uniprot.as_ref().map(|u| u.accession.as_str()),
            Some("P68871")
        );
    }
}
