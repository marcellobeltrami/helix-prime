//! Live NCBI E-utilities REST client supporting Gene, Genome, Nucleotide, PopSet,
//! Protein, PubMed, Structure, and Taxonomy searches and FASTA/GenBank downloads.

use crate::alignment::Alignment;
use crate::database::build_http_agent;
use crate::error::{HelixError, Result};
use crate::io::{fasta::parse_fasta, genbank::parse_genbank};
use serde_json::Value;

/// Public NCBI databases accessible directly from the Sources sidebar (matching Geneious Prime).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NcbiDatabase {
    Gene,
    Genome,
    Nucleotide,
    PopSet,
    Protein,
    PubMed,
    Structure,
    Taxonomy,
}

impl NcbiDatabase {
    pub const ALL: [Self; 8] = [
        Self::Gene,
        Self::Genome,
        Self::Nucleotide,
        Self::PopSet,
        Self::Protein,
        Self::PubMed,
        Self::Structure,
        Self::Taxonomy,
    ];

    /// Sidebar display name.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Gene => "Gene",
            Self::Genome => "Genome",
            Self::Nucleotide => "Nucleotide",
            Self::PopSet => "PopSet",
            Self::Protein => "Protein",
            Self::PubMed => "PubMed",
            Self::Structure => "Structure",
            Self::Taxonomy => "Taxonomy",
        }
    }

    /// Official NCBI E-utilities `db=` parameter value.
    #[must_use]
    pub const fn api_name(self) -> &'static str {
        match self {
            Self::Gene => "gene",
            Self::Genome => "genome",
            Self::Nucleotide => "nuccore",
            Self::PopSet => "popset",
            Self::Protein => "protein",
            Self::PubMed => "pubmed",
            Self::Structure => "structure",
            Self::Taxonomy => "taxonomy",
        }
    }

    /// Returns `true` if records from this database can be downloaded directly as sequences.
    #[must_use]
    pub const fn supports_sequence_download(self) -> bool {
        matches!(self, Self::Nucleotide | Self::Protein | Self::PopSet)
    }
}

/// Summary row returned by an NCBI E-utilities `esearch` + `esummary` query.
#[derive(Debug, Clone)]
pub struct NcbiSearchHit {
    pub database: NcbiDatabase,
    pub uid: String,
    pub accession: String,
    pub title: String,
    pub organism: String,
    pub length: usize,
    pub updated: String,
}

/// Searches an NCBI public database using `esearch.fcgi` + `esummary.fcgi`.
pub fn search_ncbi(db: NcbiDatabase, query: &str, limit: usize) -> Result<Vec<NcbiSearchHit>> {
    let trimmed = query.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }
    let retmax = limit.clamp(1, 20);
    let agent = build_http_agent()?;

    let search_url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={}&term={}&retmax={retmax}&retmode=json",
        db.api_name(),
        encode_query(trimmed)
    );

    let search_json: Value = agent
        .get(&search_url)
        .call()
        .map_err(|e| HelixError::database("NCBI", format!("esearch failed: {e}")))?
        .into_json()
        .map_err(|e| HelixError::database("NCBI", format!("invalid esearch JSON: {e}")))?;

    let id_list: Vec<String> = search_json
        .pointer("/esearchresult/idlist")
        .and_then(Value::as_array)
        .map(|arr| {
            arr.iter()
                .filter_map(Value::as_str)
                .map(str::to_owned)
                .collect()
        })
        .unwrap_or_default();

    if id_list.is_empty() {
        return Ok(Vec::new());
    }

    let joined_ids = id_list.join(",");
    let summary_url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={}&id={joined_ids}&retmode=json",
        db.api_name()
    );

    let sum_json: Value = agent
        .get(&summary_url)
        .call()
        .map_err(|e| HelixError::database("NCBI", format!("esummary failed: {e}")))?
        .into_json()
        .map_err(|e| HelixError::database("NCBI", format!("invalid esummary JSON: {e}")))?;

    let mut hits = Vec::with_capacity(id_list.len());
    for uid in id_list {
        let Some(rec) = sum_json.pointer(&format!("/result/{uid}")) else {
            continue;
        };
        let accession = rec
            .get("accessionversion")
            .or_else(|| rec.get("caption"))
            .or_else(|| rec.get("name"))
            .or_else(|| rec.get("pdbacc"))
            .and_then(Value::as_str)
            .unwrap_or(&uid)
            .to_owned();
        let title = rec
            .get("title")
            .or_else(|| rec.get("description"))
            .or_else(|| rec.get("scientificname"))
            .or_else(|| rec.get("pdbdescr"))
            .and_then(Value::as_str)
            .unwrap_or("NCBI Record")
            .to_owned();
        let organism = rec
            .get("organism")
            .and_then(Value::as_str)
            .or_else(|| rec.pointer("/organism/scientificname").and_then(Value::as_str))
            .unwrap_or("-")
            .to_owned();
        let length = rec
            .get("slen")
            .and_then(Value::as_u64)
            .unwrap_or(0) as usize;
        let updated = rec
            .get("updatedate")
            .or_else(|| rec.get("pubdate"))
            .or_else(|| rec.get("createdate"))
            .and_then(Value::as_str)
            .unwrap_or("-")
            .to_owned();

        hits.push(NcbiSearchHit {
            database: db,
            uid,
            accession,
            title,
            organism,
            length,
            updated,
        });
    }

    Ok(hits)
}

/// Downloads a sequence record from NCBI (`Nucleotide`, `Protein`, or `PopSet`) using `efetch.fcgi`
/// and parses it into an [`Alignment`] (GenBank with annotations for Nucleotide, FASTA for Protein/PopSet).
pub fn download_ncbi_record(db: NcbiDatabase, uid: &str, fallback_name: &str) -> Result<Alignment> {
    let clean_uid: String = uid
        .trim()
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.'))
        .collect();
    if clean_uid.is_empty() {
        return Err(HelixError::database("NCBI", "empty record identifier"));
    }

    let agent = build_http_agent()?;
    let (rettype, use_gb) = match db {
        NcbiDatabase::Nucleotide => ("gb", true),
        _ => ("fasta", false),
    };

    let url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&id={clean_uid}&rettype={rettype}&retmode=text",
        db.api_name()
    );
    let body = agent
        .get(&url)
        .call()
        .map_err(|e| HelixError::database("NCBI", format!("efetch failed: {e}")))?
        .into_string()
        .map_err(|e| HelixError::database("NCBI", format!("failed reading efetch body: {e}")))?;

    if use_gb {
        parse_genbank(fallback_name, &body).or_else(|_| parse_fasta(fallback_name, &body))
    } else {
        parse_fasta(fallback_name, &body)
    }
}

fn encode_query(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 16);
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(char::from(b));
        } else if b == b' ' {
            out.push('+');
        } else {
            out.push('%');
            let hi = (b >> 4) & 0xF;
            let lo = b & 0xF;
            out.push(if hi < 10 { (b'0' + hi) as char } else { (b'A' + hi - 10) as char });
            out.push(if lo < 10 { (b'0' + lo) as char } else { (b'A' + lo - 10) as char });
        }
    }
    out
}
