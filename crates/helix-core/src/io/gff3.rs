//! GFF3 annotation parser for overlaying genomic/protein features onto an alignment.

use crate::alignment::Alignment;
use crate::annotation::{Annotation, FeatureKind, UniProtLink};
use crate::error::{HelixError, Result};
use crate::sequence::Strand;

/// Parses GFF3 records and attaches them to `alignment`, returning the count of imported features.
pub fn import_gff3_annotations(alignment: &mut Alignment, input: &str) -> Result<usize> {
    let mut imported = 0usize;

    for (line_idx, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(HelixError::parse(
                "GFF3",
                Some(line_idx + 1),
                format!("expected 9 tab-separated columns, found {}", cols.len()),
            ));
        }

        let seq_id = cols[0];
        let feat_type = cols[2];
        let start_1: usize = cols[3].parse().map_err(|_| {
            HelixError::parse("GFF3", Some(line_idx + 1), "invalid start coordinate")
        })?;
        let end_1: usize = cols[4].parse().map_err(|_| {
            HelixError::parse("GFF3", Some(line_idx + 1), "invalid end coordinate")
        })?;
        let strand = if cols[6] == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        let attrs = parse_gff3_attributes(cols[8]);

        let label = attrs
            .iter()
            .find(|(k, _)| matches!(k.as_str(), "Name" | "gene" | "ID" | "product"))
            .map(|(_, v)| v.clone())
            .unwrap_or_else(|| feat_type.to_owned());

        let start_0 = start_1.saturating_sub(1);
        let end_0 = end_1.max(start_0 + 1);

        let mut ann =
            Annotation::new(0, label, FeatureKind::from_type_str(feat_type), start_0, end_0)
                .with_strand(strand);

        if let Some(row_idx) = alignment
            .sequences()
            .iter()
            .position(|s| s.name.eq_ignore_ascii_case(seq_id))
        {
            ann.target_sequence = Some(row_idx);
        }

        for (k, v) in &attrs {
            if k.eq_ignore_ascii_case("Note") || k.eq_ignore_ascii_case("product") {
                ann.notes = v.clone();
            }
            if k.eq_ignore_ascii_case("Dbxref") || k.eq_ignore_ascii_case("UniProt") {
                for entry in v.split(',') {
                    if let Some(acc) = entry
                        .strip_prefix("UniProtKB:")
                        .or_else(|| entry.strip_prefix("UniProt:"))
                    {
                        ann.uniprot = Some(UniProtLink::new(acc));
                    }
                }
            }
        }
        ann.qualifiers = attrs;
        alignment.add_annotation(ann);
        imported += 1;
    }

    Ok(imported)
}

fn parse_gff3_attributes(raw: &str) -> Vec<(String, String)> {
    raw.split(';')
        .filter_map(|part| {
            let trimmed = part.trim();
            if trimmed.is_empty() {
                return None;
            }
            let (k, v) = trimmed.split_once('=')?;
            Some((k.trim().to_owned(), v.trim().replace("%20", " ")))
        })
        .collect()
}
