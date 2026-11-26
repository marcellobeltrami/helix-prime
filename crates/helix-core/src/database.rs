//! Live REST clients for UniProtKB, AlphaFold Protein Structure Database (EBI), and RCSB PDB.
//!
//! All network calls are synchronous and designed to execute on a background worker thread
//! in the GUI so the interface stays completely responsive. Downloaded 3D structure files are
//! cached on disk inside the OS temporary directory to make repeated lookups instantaneous.

use crate::alignment::Alignment;
use crate::annotation::{Annotation, FeatureKind, UniProtLink};
use crate::error::{HelixError, Result};
use crate::sequence::is_gap;
use crate::structure::{parse_mmcif, parse_pdb, ProteinStructure};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Duration;

const USER_AGENT: &str = "HelixPrime/0.1 (Open-Source Alignment Inspector)";

/// Remote database source selector for fetching 3D protein structures and metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum ProteinStructureSource {
    /// Fetch 3D model and pLDDT metadata directly from AlphaFold Protein Structure Database (EBI).
    #[default]
    AlphaFoldDb,
    /// Fetch metadata, functional features, and linked experimental PDB (or AFDB fallback) via UniProtKB.
    UniProtDb,
    /// Try AlphaFold DB first, automatically falling back to UniProtKB + experimental PDB if unavailable.
    AutoFallback,
}

impl ProteinStructureSource {
    pub const ALL: [Self; 3] = [Self::AlphaFoldDb, Self::UniProtDb, Self::AutoFallback];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::AlphaFoldDb => "AlphaFold DB",
            Self::UniProtDb => "UniProtKB (PDB / Meta)",
            Self::AutoFallback => "Auto (AlphaFold / UniProt)",
        }
    }
}

/// Combined result returned when fetching a 3D protein structure + metadata from either
/// AlphaFold DB or UniProtKB.
#[derive(Debug, Clone)]
pub struct FetchedProteinBundle {
    pub structure: ProteinStructure,
    pub alphafold_meta: Option<AlphaFoldMeta>,
    pub uniprot_entry: Option<UniProtEntry>,
    pub resolved_source: ProteinStructureSource,
}

/// Summary metadata returned by a UniProtKB search query.
#[derive(Debug, Clone)]
pub struct UniProtSearchHit {
    pub accession: String,
    pub entry_name: String,
    pub protein_name: String,
    pub gene_name: String,
    pub organism: String,
    pub length: usize,
    pub reviewed: bool,
}

/// A curated sequence feature returned by UniProtKB (domain, active site, binding site, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniProtFeature {
    pub feature_type: String,
    pub description: String,
    /// 1-based inclusive start residue in the canonical UniProt sequence.
    pub start: usize,
    /// 1-based inclusive end residue in the canonical UniProt sequence.
    pub end: usize,
}

/// Full UniProtKB protein record with functional annotations and structural cross-references.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniProtEntry {
    pub accession: String,
    pub entry_name: String,
    pub protein_name: String,
    pub gene_names: Vec<String>,
    pub organism: String,
    pub function_summary: String,
    pub sequence: String,
    pub alphafold_id: Option<String>,
    pub pdb_ids: Vec<String>,
    pub features: Vec<UniProtFeature>,
}

impl UniProtEntry {
    /// Converts curated UniProt functional features (domains, active sites, binding sites,
    /// motifs, regions) into [`Annotation`] objects mapped onto `alignment` row `target_row`.
    #[must_use]
    pub fn build_alignment_annotations(
        &self,
        alignment: &Alignment,
        target_row: usize,
    ) -> Vec<Annotation> {
        let Some(seq) = alignment.sequences().get(target_row) else {
            return Vec::new();
        };

        // Build map from 1-based ungapped residue number -> 0-based alignment column.
        let mut res_to_col: Vec<usize> = Vec::with_capacity(seq.len());
        for (col, &b) in seq.residues.iter().enumerate() {
            if !is_gap(b) {
                res_to_col.push(col);
            }
        }
        if res_to_col.is_empty() {
            return Vec::new();
        }

        // Only import biologically informative tracks (skip hundreds of minor variants/conflicts).
        let mut out = Vec::new();
        for feat in &self.features {
            let kind = match feat.feature_type.as_str() {
                "Domain" | "Region" | "DNA binding" | "Zinc finger" | "Motif" | "Repeat" => {
                    FeatureKind::Domain
                }
                "Active site" => FeatureKind::ActiveSite,
                "Binding site" | "Metal binding" | "Site" => FeatureKind::BindingSite,
                _ => continue,
            };

            if feat.start == 0 || feat.end < feat.start {
                continue;
            }
            let idx_start = feat.start - 1;
            let idx_end = feat.end - 1;
            if idx_start >= res_to_col.len() {
                continue;
            }

            let col_start = res_to_col[idx_start];
            let col_end = res_to_col[idx_end.min(res_to_col.len() - 1)].saturating_add(1);
            let label = if feat.description.is_empty() {
                feat.feature_type.clone()
            } else {
                feat.description.clone()
            };

            let link = UniProtLink::new(&self.accession)
                .with_start_residue(feat.start)
                .with_metadata(&self.entry_name, &self.protein_name, &self.organism);

            let ann = Annotation::new(0, label, kind, col_start, col_end)
                .with_target_sequence(target_row)
                .with_notes(format!(
                    "UniProt {} ({}..{}): {}",
                    self.accession, feat.start, feat.end, feat.feature_type
                ))
                .with_uniprot(link);
            out.push(ann);
        }

        out
    }
}

/// Metadata summary returned by the AlphaFold Protein Structure Database API.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlphaFoldMeta {
    pub entry_id: String,
    pub uniprot_accession: String,
    pub uniprot_description: String,
    pub organism: String,
    pub latest_version: u32,
    pub global_plddt: Option<f32>,
    pub pdb_url: Option<String>,
    pub cif_url: Option<String>,
}

pub(crate) fn build_http_agent() -> Result<ureq::Agent> {
    // `ureq` with `native-tls` (and `default-features = false`) requires explicitly
    // attaching a `native_tls::TlsConnector` so HTTPS (Windows SChannel / macOS / OpenSSL)
    // requests are enabled. `.try_proxy_from_env(true)` respects system HTTP(S)_PROXY settings.
    let tls = native_tls::TlsConnector::new()
        .map_err(|e| HelixError::database("TLS", format!("failed to initialize native TLS: {e}")))?;

    Ok(ureq::AgentBuilder::new()
        .tls_connector(Arc::new(tls))
        .try_proxy_from_env(true)
        .timeout_connect(Duration::from_secs(15))
        .timeout_read(Duration::from_secs(35))
        .user_agent(USER_AGENT)
        .build())
}

/// Executes an HTTP GET request with automatic retries + short backoff on transient DNS
/// (`os error 11002` / `WSATRY_AGAIN`), connection resets, or rate-limit (`429`/`503`) responses.
pub(crate) fn http_get_with_retry(
    agent: &ureq::Agent,
    url: &str,
    accept: Option<&str>,
    service: &'static str,
) -> Result<ureq::Response> {
    const MAX_ATTEMPTS: usize = 3;
    let mut last_err = String::new();

    for attempt in 0..MAX_ATTEMPTS {
        let mut req = agent.get(url);
        if let Some(hdr) = accept {
            req = req.set("Accept", hdr);
        }
        match req.call() {
            Ok(resp) => return Ok(resp),
            Err(ureq::Error::Status(code, resp)) => {
                // Retry only on transient server/rate-limit statuses.
                let retryable = matches!(code, 429 | 500 | 502 | 503 | 504);
                let body_hint = resp.status_text().to_owned();
                last_err = format!("HTTP {code} ({body_hint}) for {url}");
                if !retryable || attempt + 1 == MAX_ATTEMPTS {
                    return Err(HelixError::database(service, last_err));
                }
            }
            Err(ureq::Error::Transport(t)) => {
                last_err = format!("{url}: {t}");
                if attempt + 1 == MAX_ATTEMPTS {
                    return Err(HelixError::database(service, last_err));
                }
            }
        }
        // Exponential backoff: 350ms, then 900ms.
        let delay_ms = if attempt == 0 { 350 } else { 900 };
        std::thread::sleep(Duration::from_millis(delay_ms));
    }

    Err(HelixError::database(service, last_err))
}

/// Searches UniProtKB by gene name, protein name, or free-text query.
pub fn search_uniprot(query: &str, limit: usize) -> Result<Vec<UniProtSearchHit>> {
    let trimmed = query.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }
    let capped = limit.clamp(1, 25);
    let url = format!(
        "https://rest.uniprot.org/uniprotkb/search?query={}&format=json&size={capped}",
        percent_encode(trimmed)
    );

    let agent = build_http_agent()?;
    let resp = http_get_with_retry(&agent, &url, Some("application/json"), "UniProt")?;

    let json: Value = resp
        .into_json()
        .map_err(|e| HelixError::database("UniProt", format!("invalid JSON response: {e}")))?;

    let results = json
        .get("results")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();

    let mut hits = Vec::with_capacity(results.len());
    for item in &results {
        let accession = item
            .get("primaryAccession")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_owned();
        if accession.is_empty() {
            continue;
        }
        let entry_name = item
            .get("uniProtkbId")
            .and_then(Value::as_str)
            .unwrap_or(&accession)
            .to_owned();
        let protein_name = extract_protein_name(item);
        let organism = item
            .pointer("/organism/scientificName")
            .and_then(Value::as_str)
            .unwrap_or("Unknown organism")
            .to_owned();
        let gene_name = item
            .pointer("/genes/0/geneName/value")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_owned();
        let length = item
            .pointer("/sequence/length")
            .and_then(Value::as_u64)
            .unwrap_or(0) as usize;
        let reviewed = item
            .get("entryType")
            .and_then(Value::as_str)
            .is_some_and(|t| t.contains("reviewed"));

        hits.push(UniProtSearchHit {
            accession,
            entry_name,
            protein_name,
            gene_name,
            organism,
            length,
            reviewed,
        });
    }

    Ok(hits)
}

/// Fetches a full UniProtKB entry (metadata, sequence, curated features, PDB/AlphaFold refs).
pub fn fetch_uniprot_entry(accession: &str) -> Result<UniProtEntry> {
    let clean_acc = sanitize_accession(accession)?;
    let url = format!("https://rest.uniprot.org/uniprotkb/{clean_acc}.json");

    let agent = build_http_agent()?;
    let resp = http_get_with_retry(&agent, &url, Some("application/json"), "UniProt")?;

    let json: Value = resp
        .into_json()
        .map_err(|e| HelixError::database("UniProt", format!("invalid JSON payload: {e}")))?;

    let primary_acc = json
        .get("primaryAccession")
        .and_then(Value::as_str)
        .unwrap_or(&clean_acc)
        .to_owned();
    let entry_name = json
        .get("uniProtkbId")
        .and_then(Value::as_str)
        .unwrap_or(&primary_acc)
        .to_owned();
    let protein_name = extract_protein_name(&json);
    let organism = json
        .pointer("/organism/scientificName")
        .and_then(Value::as_str)
        .unwrap_or("Unknown organism")
        .to_owned();
    let sequence = json
        .pointer("/sequence/value")
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_owned();

    let mut gene_names = Vec::new();
    if let Some(genes) = json.get("genes").and_then(Value::as_array) {
        for g in genes {
            if let Some(name) = g.pointer("/geneName/value").and_then(Value::as_str) {
                gene_names.push(name.to_owned());
            }
        }
    }

    // Extract primary FUNCTION comment if present.
    let mut function_summary = String::new();
    if let Some(comments) = json.get("comments").and_then(Value::as_array) {
        for c in comments {
            if c.get("commentType").and_then(Value::as_str) == Some("FUNCTION") {
                if let Some(text) = c.pointer("/texts/0/value").and_then(Value::as_str) {
                    function_summary = text.to_owned();
                    break;
                }
            }
        }
    }

    // Extract AlphaFoldDB and experimental PDB cross-references.
    let mut alphafold_id = None;
    let mut pdb_ids = Vec::new();
    if let Some(xrefs) = json
        .get("uniProtKBCrossReferences")
        .and_then(Value::as_array)
    {
        for xr in xrefs {
            let db = xr.get("database").and_then(Value::as_str).unwrap_or("");
            let id = xr.get("id").and_then(Value::as_str).unwrap_or("");
            if db == "AlphaFoldDB" && !id.is_empty() {
                alphafold_id = Some(id.to_owned());
            } else if db == "PDB" && !id.is_empty() && pdb_ids.len() < 12 {
                pdb_ids.push(id.to_owned());
            }
        }
    }

    // Extract curated sequence features.
    let mut features = Vec::new();
    if let Some(raw_feats) = json.get("features").and_then(Value::as_array) {
        for f in raw_feats {
            let ftype = f
                .get("type")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_owned();
            let desc = f
                .get("description")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_owned();
            let start = f
                .pointer("/location/start/value")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize;
            let end = f
                .pointer("/location/end/value")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize;

            if !ftype.is_empty() && start > 0 && end >= start {
                features.push(UniProtFeature {
                    feature_type: ftype,
                    description: desc,
                    start,
                    end,
                });
            }
        }
    }

    Ok(UniProtEntry {
        accession: primary_acc,
        entry_name,
        protein_name,
        gene_names,
        organism,
        function_summary,
        sequence,
        alphafold_id,
        pdb_ids,
        features,
    })
}

/// Fetches an AlphaFold DB 3D structure by UniProt accession (with automatic local disk caching).
///
/// Queries the EBI AlphaFold API (`/api/prediction/{accession}`) to discover the latest model
/// version URL, then downloads and parses the PDB (falling back to mmCIF if PDB is omitted).
pub fn fetch_alphafold_structure(accession: &str) -> Result<(ProteinStructure, AlphaFoldMeta)> {
    let clean_acc = sanitize_accession(accession)?;
    let api_url = format!("https://alphafold.ebi.ac.uk/api/prediction/{clean_acc}");
    let agent = build_http_agent()?;
    let resp = http_get_with_retry(&agent, &api_url, Some("application/json"), "AlphaFold DB")?;

    let json: Value = resp
        .into_json()
        .map_err(|e| HelixError::database("AlphaFold DB", format!("invalid JSON metadata: {e}")))?;

    let entries = json.as_array().ok_or_else(|| {
        HelixError::database(
            "AlphaFold DB",
            format!("no AlphaFold prediction available for '{clean_acc}'"),
        )
    })?;

    let first = entries.first().ok_or_else(|| {
        HelixError::database(
            "AlphaFold DB",
            format!("empty AlphaFold prediction list for '{clean_acc}'"),
        )
    })?;

    let entry_id = first
        .get("entryId")
        .and_then(Value::as_str)
        .unwrap_or(&clean_acc)
        .to_owned();
    let description = first
        .get("uniprotDescription")
        .and_then(Value::as_str)
        .unwrap_or(&entry_id)
        .to_owned();
    let organism = first
        .get("organismScientificName")
        .and_then(Value::as_str)
        .unwrap_or("Unknown organism")
        .to_owned();
    let latest_version = first
        .get("latestVersion")
        .and_then(Value::as_u64)
        .unwrap_or(4) as u32;
    let global_plddt = first
        .get("globalMetricValue")
        .and_then(Value::as_f64)
        .map(|v| v as f32);
    let pdb_url = first
        .get("pdbUrl")
        .and_then(Value::as_str)
        .map(str::to_owned);
    let cif_url = first
        .get("cifUrl")
        .and_then(Value::as_str)
        .map(str::to_owned);

    let meta = AlphaFoldMeta {
        entry_id: entry_id.clone(),
        uniprot_accession: clean_acc.clone(),
        uniprot_description: description.clone(),
        organism: organism.clone(),
        latest_version,
        global_plddt,
        pdb_url: pdb_url.clone(),
        cif_url: cif_url.clone(),
    };

    let source_label = format!("AlphaFold DB v{latest_version}");

    let mut structure = if let Some(url) = &pdb_url {
        let cache_name = format!("{entry_id}_v{latest_version}.pdb");
        let pdb_text = fetch_text_cached(&agent, url, &cache_name, "AlphaFold DB")?;
        parse_pdb(
            &pdb_text,
            &entry_id,
            Some(clean_acc.clone()),
            &source_label,
            true,
        )?
    } else if let Some(url) = &cif_url {
        let cache_name = format!("{entry_id}_v{latest_version}.cif");
        let cif_text = fetch_text_cached(&agent, url, &cache_name, "AlphaFold DB")?;
        parse_mmcif(
            &cif_text,
            &entry_id,
            Some(clean_acc.clone()),
            &source_label,
            true,
        )?
    } else {
        return Err(HelixError::database(
            "AlphaFold DB",
            format!("no PDB or mmCIF download URL in AlphaFold metadata for '{clean_acc}'"),
        ));
    };

    structure.title = description;
    structure.organism = Some(organism);
    Ok((structure, meta))
}

/// Fetches a 3D protein structure and metadata from either **AlphaFold DB** or **UniProtKB**
/// (or automatic fallback between them).
///
/// - [`ProteinStructureSource::AlphaFoldDb`]: Downloads the AlphaFold DB v4 model & pLDDT
///   metadata directly from `alphafold.ebi.ac.uk`, enriching with UniProtKB metadata if reachable.
/// - [`ProteinStructureSource::UniProtDb`]: Queries UniProtKB first for curated metadata,
///   functional annotations, and structural cross-references; loads the primary experimental
///   PDB linked in UniProt (falling back to UniProt's AlphaFoldDB cross-reference if no
///   experimental PDB exists) and stamps UniProt's curated protein name & organism onto the model.
/// - [`ProteinStructureSource::AutoFallback`]: Attempts AlphaFold DB first, transparently falling
///   back to UniProtKB + experimental PDB if AlphaFold DB has no prediction for `accession`.
pub fn fetch_protein_bundle(
    accession: &str,
    source: ProteinStructureSource,
) -> Result<FetchedProteinBundle> {
    let clean_acc = sanitize_accession(accession)?;

    match source {
        ProteinStructureSource::AlphaFoldDb => {
            let (mut structure, meta) = fetch_alphafold_structure(&clean_acc)?;
            // Best-effort metadata enrichment from UniProtKB so function & PDB xrefs populate.
            let uniprot_entry = fetch_uniprot_entry(&clean_acc).ok();
            if let Some(ref up) = uniprot_entry {
                if structure.title.is_empty() || structure.title == structure.entry_id {
                    structure.title = up.protein_name.clone();
                }
                if structure.organism.is_none() {
                    structure.organism = Some(up.organism.clone());
                }
            }
            Ok(FetchedProteinBundle {
                structure,
                alphafold_meta: Some(meta),
                uniprot_entry,
                resolved_source: ProteinStructureSource::AlphaFoldDb,
            })
        }
        ProteinStructureSource::UniProtDb => {
            let up = fetch_uniprot_entry(&clean_acc)?;
            // Prioritize experimental PDB cross-references curated by UniProtKB.
            let mut loaded_structure: Option<ProteinStructure> = None;
            let mut af_meta: Option<AlphaFoldMeta> = None;

            for pdb_id in up.pdb_ids.iter().take(3) {
                if let Ok(mut st) = fetch_rcsb_pdb(pdb_id, Some(up.accession.clone())) {
                    st.title = format!("{} (UniProt {} / PDB {pdb_id})", up.protein_name, up.accession);
                    st.organism = Some(up.organism.clone());
                    st.source_label = format!("UniProtKB -> RCSB PDB ({pdb_id})");
                    loaded_structure = Some(st);
                    break;
                }
            }

            // If UniProt has no experimental PDB (or download failed), follow UniProt's AlphaFoldDB xref.
            if loaded_structure.is_none() {
                let af_acc = up.alphafold_id.as_deref().unwrap_or(&up.accession);
                let (mut st, meta) = fetch_alphafold_structure(af_acc)?;
                st.title = format!("{} (UniProt {})", up.protein_name, up.accession);
                st.organism = Some(up.organism.clone());
                st.source_label = format!("UniProtKB -> AlphaFold DB v{}", meta.latest_version);
                loaded_structure = Some(st);
                af_meta = Some(meta);
            }

            let structure = loaded_structure.ok_or_else(|| {
                HelixError::database(
                    "UniProt",
                    format!("UniProt entry '{}' has no downloadable PDB or AlphaFold 3D structure", up.accession),
                )
            })?;

            Ok(FetchedProteinBundle {
                structure,
                alphafold_meta: af_meta,
                uniprot_entry: Some(up),
                resolved_source: ProteinStructureSource::UniProtDb,
            })
        }
        ProteinStructureSource::AutoFallback => {
            match fetch_protein_bundle(&clean_acc, ProteinStructureSource::AlphaFoldDb) {
                Ok(bundle) => Ok(bundle),
                Err(_) => fetch_protein_bundle(&clean_acc, ProteinStructureSource::UniProtDb),
            }
        }
    }
}

/// Fetches an experimental 3D structure directly from RCSB PDB by 4-character PDB ID.
pub fn fetch_rcsb_pdb(pdb_id: &str, accession: Option<String>) -> Result<ProteinStructure> {
    let clean = pdb_id.trim().to_ascii_uppercase();
    if clean.len() != 4 || !clean.chars().all(|c| c.is_ascii_alphanumeric()) {
        return Err(HelixError::database(
            "RCSB PDB",
            format!("invalid 4-character PDB ID '{pdb_id}'"),
        ));
    }
    let url = format!("https://files.rcsb.org/download/{clean}.pdb");
    let agent = build_http_agent()?;
    let cache_name = format!("rcsb_{clean}.pdb");
    let text = fetch_text_cached(&agent, &url, &cache_name, "RCSB PDB")?;
    parse_pdb(&text, &clean, accession, format!("RCSB PDB ({clean})"), false)
}

/// Downloads a text resource with best-effort disk caching inside `<temp_dir>/helix_prime_cache`.
fn fetch_text_cached(
    agent: &ureq::Agent,
    url: &str,
    cache_filename: &str,
    service: &'static str,
) -> Result<String> {
    let cache_path = cache_file_path(cache_filename);
    if let Some(ref path) = cache_path {
        if let Ok(cached) = std::fs::read_to_string(path) {
            if !cached.is_empty() {
                return Ok(cached);
            }
        }
    }

    let resp = http_get_with_retry(agent, url, None, service)?;

    let text = resp
        .into_string()
        .map_err(|e| HelixError::database(service, format!("failed reading response body: {e}")))?;

    if let Some(path) = cache_path {
        if let Some(parent) = path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        let _ = std::fs::write(path, &text);
    }

    Ok(text)
}

fn cache_file_path(filename: &str) -> Option<PathBuf> {
    // Sanitize filename to prevent path traversal.
    let safe: String = filename
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.'))
        .collect();
    if safe.is_empty() {
        return None;
    }
    Some(std::env::temp_dir().join("helix_prime_cache").join(safe))
}

fn sanitize_accession(raw: &str) -> Result<String> {
    let clean = raw.trim().to_ascii_uppercase();
    if clean.is_empty()
        || !clean
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_')
    {
        return Err(HelixError::database(
            "UniProt",
            format!("invalid accession identifier '{raw}'"),
        ));
    }
    Ok(clean)
}

fn extract_protein_name(json: &Value) -> String {
    json.pointer("/proteinDescription/recommendedName/fullName/value")
        .or_else(|| json.pointer("/proteinDescription/submissionNames/0/fullName/value"))
        .or_else(|| json.pointer("/proteinDescription/alternativeNames/0/fullName/value"))
        .and_then(Value::as_str)
        .unwrap_or("Unnamed protein")
        .to_owned()
}

/// Minimal RFC 3986 percent-encoder for query parameters.
fn percent_encode(input: &str) -> String {
    let mut out = String::with_capacity(input.len() + 16);
    for b in input.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(char::from(b));
        } else {
            out.push('%');
            let hi = (b >> 4) & 0xF;
            let lo = b & 0xF;
            out.push(hex_nibble(hi));
            out.push(hex_nibble(lo));
        }
    }
    out
}

const fn hex_nibble(n: u8) -> char {
    match n {
        0..=9 => (b'0' + n) as char,
        _ => (b'A' + (n - 10)) as char,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percent_encoding_handles_spaces_and_special_chars() {
        assert_eq!(percent_encode("TP53 human"), "TP53%20human");
        assert_eq!(percent_encode("gene:HBB"), "gene%3AHBB");
    }

    #[test]
    fn sanitize_accession_rejects_invalid_characters() {
        assert_eq!(sanitize_accession(" p04637 ").expect("ok"), "P04637");
        assert!(sanitize_accession("P04637/../bad").is_err());
        assert!(sanitize_accession("").is_err());
    }

    #[test]
    fn http_agent_initializes_native_tls_connector() {
        assert!(build_http_agent().is_ok());
    }
}
