//! Atomic, per-user workspace session persistence (`session.json`).
//!
//! Saves the complete HelixPrime workspace (documents, edits, annotations, trees, lineage,
//! active selection, cursor & scroll position, view/display toggles, circular camera, tool
//! preferences, theme, and 3D protein viewer state) into a standard user-writable directory
//! so the app reopens right where the user left off — requiring zero administrator privileges.

use crate::app::{
    CircularViewState, DisplaySettings, DocTab, HelixApp, InspectorTab, LineageEntry,
    OperationRecord, OperationStatus, SourceSelection, WorkspaceDoc,
};
use crate::renderer3d::{Camera3D, ColorBy3D, RenderStyle};
use crate::theme::ThemeMode;
use helix_core::aligners::{AlignerAlgorithm, Minimap2Preset};
use helix_core::blast::BlastConfig;
use helix_core::cloning::CloningMethod;
use helix_core::database::{AlphaFoldMeta, ProteinStructureSource, UniProtEntry};
use helix_core::phylo::{DistanceModel, PhyloTree, TreeMethod};
use helix_core::{Alignment, ProteinStructure};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

/// Current on-disk session schema version. Incompatible versions fall back cleanly to defaults.
pub const SESSION_SCHEMA_VERSION: u32 = 1;

/// Maximum atom count embedded inline in `session.json` to keep close/open latency sub-second.
const MAX_PERSISTED_3D_ATOMS: usize = 60_000;

/// Serializable subset of [`WorkspaceDoc`] excluding lazily recomputed caches.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PersistedDoc {
    pub id: u64,
    pub folder: String,
    pub modified: String,
    pub is_circular: bool,
    #[serde(default)]
    pub deleted: bool,
    #[serde(default)]
    pub shared: bool,
    #[serde(default)]
    pub checked: bool,
    pub alignment: Alignment,
    #[serde(default)]
    pub tree: Option<PhyloTree>,
    #[serde(default)]
    pub lineage: Vec<LineageEntry>,
}

impl From<&WorkspaceDoc> for PersistedDoc {
    fn from(doc: &WorkspaceDoc) -> Self {
        Self {
            id: doc.id,
            folder: doc.folder.clone(),
            modified: doc.modified.clone(),
            is_circular: doc.is_circular,
            deleted: doc.deleted,
            shared: doc.shared,
            checked: doc.checked,
            alignment: doc.alignment.clone(),
            tree: doc.tree.clone(),
            lineage: doc.lineage.clone(),
        }
    }
}

impl From<PersistedDoc> for WorkspaceDoc {
    fn from(doc: PersistedDoc) -> Self {
        Self {
            id: doc.id,
            folder: doc.folder,
            modified: doc.modified,
            is_circular: doc.is_circular,
            deleted: doc.deleted,
            shared: doc.shared,
            checked: doc.checked,
            alignment: doc.alignment,
            tree: doc.tree,
            lineage: doc.lineage,
            pairwise_cache: None,
            restriction_cache: None,
        }
    }
}

/// Persisted preferences for the tool dialogs so repeated analyses keep user parameters.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct PersistedToolPrefs {
    pub assemble_mode: usize,
    pub assemble_algorithm: AlignerAlgorithm,
    pub minimap2_preset: Minimap2Preset,
    pub minimap2_custom_path: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub min_overlap: usize,
    pub blast_config: BlastConfig,
    pub blast_query: String,
    pub cloning_method: CloningMethod,
    pub cloning_enzyme_5prime: String,
    pub cloning_enzyme_3prime: String,
    pub cloning_reverse_insert: bool,
    pub tree_method: TreeMethod,
    pub tree_model: DistanceModel,
    pub workflow_auto_annotate: bool,
    pub workflow_realign: bool,
    pub workflow_translate: bool,
    pub workflow_build_tree: bool,
    pub workflow_design_primers: bool,
}

/// Persisted state for the separate OS 3D Protein Viewer window.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct PersistedProteinViewer {
    pub open: bool,
    pub accession_input: String,
    pub source: ProteinStructureSource,
    pub linked_annotation_id: Option<u64>,
    pub camera: Camera3D,
    pub style: RenderStyle,
    pub color_by: ColorBy3D,
    pub selected_residue: Option<usize>,
    pub structure: Option<ProteinStructure>,
    pub alphafold_meta: Option<AlphaFoldMeta>,
    pub uniprot_entry: Option<UniProtEntry>,
}

/// Versioned snapshot of the full HelixPrime application session.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct PersistedSession {
    pub schema_version: u32,
    pub theme: ThemeMode,
    pub docs: Vec<PersistedDoc>,
    pub next_doc_id: u64,
    pub active_doc: usize,
    pub source_selection: SourceSelection,
    pub active_tab: DocTab,
    pub global_search: String,
    pub history: Vec<usize>,
    pub history_pos: usize,
    pub show_sources_panel: bool,
    pub show_inspector_panel: bool,
    pub show_doc_table: bool,
    pub inspector_tab: InspectorTab,
    pub display: DisplaySettings,
    pub circular: CircularViewState,
    pub horizontal_scroll_px: f32,
    pub cursor_col: usize,
    pub cursor_row: usize,
    pub selected_cols: Option<(usize, usize)>,
    pub selected_annotation_id: Option<u64>,
    pub motif_query: String,
    pub motif_allow_iupac: bool,
    pub annotation_filter: String,
    pub uniprot_search_query: String,
    pub ncbi_search_query: String,
    pub operations: Vec<OperationRecord>,
    pub next_op_id: u64,
    pub tool_prefs: PersistedToolPrefs,
    pub viewer3d: PersistedProteinViewer,
}

impl Default for PersistedSession {
    fn default() -> Self {
        Self {
            schema_version: SESSION_SCHEMA_VERSION,
            theme: ThemeMode::Dark,
            docs: Vec::new(),
            next_doc_id: 1,
            active_doc: 0,
            source_selection: SourceSelection::AllLocal,
            active_tab: DocTab::SequenceView,
            global_search: String::new(),
            history: Vec::new(),
            history_pos: 0,
            show_sources_panel: true,
            show_inspector_panel: true,
            show_doc_table: true,
            inspector_tab: InspectorTab::Display,
            display: DisplaySettings::default(),
            circular: CircularViewState::default(),
            horizontal_scroll_px: 0.0,
            cursor_col: 0,
            cursor_row: 0,
            selected_cols: None,
            selected_annotation_id: None,
            motif_query: String::new(),
            motif_allow_iupac: true,
            annotation_filter: String::new(),
            uniprot_search_query: String::from("TP53 human"),
            ncbi_search_query: String::from("Escherichia phage lambda"),
            operations: Vec::new(),
            next_op_id: 1,
            tool_prefs: PersistedToolPrefs::default(),
            viewer3d: PersistedProteinViewer::default(),
        }
    }
}

impl PersistedSession {
    /// Captures a serializable snapshot from the live [`HelixApp`] state.
    #[must_use]
    pub fn from_app(app: &HelixApp) -> Self {
        let structure = app.viewer3d.structure.as_ref().and_then(|st| {
            if st.atoms.len() <= MAX_PERSISTED_3D_ATOMS {
                Some(st.clone())
            } else {
                None
            }
        });

        // Cap persisted operations history to keep session files lean.
        let operations: Vec<OperationRecord> = app.operations.iter().take(40).cloned().collect();

        Self {
            schema_version: SESSION_SCHEMA_VERSION,
            theme: app.theme,
            docs: app.docs.iter().map(PersistedDoc::from).collect(),
            next_doc_id: app.next_doc_id,
            active_doc: app.active_doc,
            source_selection: app.source_selection.clone(),
            active_tab: app.active_tab,
            global_search: app.global_search.clone(),
            history: app.history.clone(),
            history_pos: app.history_pos,
            show_sources_panel: app.show_sources_panel,
            show_inspector_panel: app.show_inspector_panel,
            show_doc_table: app.show_doc_table,
            inspector_tab: app.inspector_tab,
            display: app.display.clone(),
            circular: app.circular.clone(),
            horizontal_scroll_px: app.horizontal_scroll_px,
            cursor_col: app.cursor_col,
            cursor_row: app.cursor_row,
            selected_cols: app.selected_cols,
            selected_annotation_id: app.selected_annotation_id,
            motif_query: app.motif_query.clone(),
            motif_allow_iupac: app.motif_allow_iupac,
            annotation_filter: app.annotation_filter.clone(),
            uniprot_search_query: app.uniprot_search_query.clone(),
            ncbi_search_query: app.ncbi_search_query.clone(),
            operations,
            next_op_id: app.next_op_id,
            tool_prefs: PersistedToolPrefs {
                assemble_mode: app.assemble_modal.mode,
                assemble_algorithm: app.assemble_modal.algorithm,
                minimap2_preset: app.assemble_modal.minimap2_preset,
                minimap2_custom_path: app.assemble_modal.minimap2_custom_path.clone(),
                gap_open: app.assemble_modal.gap_open,
                gap_extend: app.assemble_modal.gap_extend,
                min_overlap: app.assemble_modal.min_overlap,
                blast_config: app.blast_modal.config,
                blast_query: app.blast_modal.query.clone(),
                cloning_method: app.cloning_modal.method,
                cloning_enzyme_5prime: app.cloning_modal.enzyme_5prime.clone(),
                cloning_enzyme_3prime: app.cloning_modal.enzyme_3prime.clone(),
                cloning_reverse_insert: app.cloning_modal.reverse_insert,
                tree_method: app.tree_modal.method,
                tree_model: app.tree_modal.model,
                workflow_auto_annotate: app.workflow_modal.step_auto_annotate,
                workflow_realign: app.workflow_modal.step_realign,
                workflow_translate: app.workflow_modal.step_translate,
                workflow_build_tree: app.workflow_modal.step_build_tree,
                workflow_design_primers: app.workflow_modal.step_design_primers,
            },
            viewer3d: PersistedProteinViewer {
                open: app.viewer3d.open,
                accession_input: app.viewer3d.accession_input.clone(),
                source: app.viewer3d.source,
                linked_annotation_id: app.viewer3d.linked_annotation_id,
                camera: app.viewer3d.camera.clone(),
                style: app.viewer3d.style,
                color_by: app.viewer3d.color_by,
                selected_residue: app.viewer3d.selected_residue,
                structure,
                alphafold_meta: app.viewer3d.alphafold_meta.clone(),
                uniprot_entry: app.viewer3d.uniprot_entry.clone(),
            },
        }
    }

    /// Applies this persisted snapshot onto `app`, sanitizing every index and floating-point
    /// bound so a corrupt or partially edited file can never cause an out-of-bounds panic.
    pub fn apply_to(self, app: &mut HelixApp) -> bool {
        if self.schema_version != SESSION_SCHEMA_VERSION || self.docs.is_empty() {
            return false;
        }

        let docs: Vec<WorkspaceDoc> = self.docs.into_iter().map(WorkspaceDoc::from).collect();
        let doc_count = docs.len();
        let max_doc_id = docs.iter().map(|d| d.id).max().unwrap_or(0);

        app.theme = self.theme;
        app.docs = docs;
        app.next_doc_id = self.next_doc_id.max(max_doc_id.saturating_add(1));
        app.active_doc = self.active_doc.min(doc_count.saturating_sub(1));
        app.source_selection = self.source_selection;
        app.active_tab = self.active_tab;
        app.global_search = self.global_search;

        // Sanitize navigation history against the loaded document count.
        let mut history: Vec<usize> = self
            .history
            .into_iter()
            .filter(|&idx| idx < doc_count)
            .collect();
        if history.is_empty() {
            history.push(app.active_doc);
        }
        app.history_pos = self.history_pos.min(history.len().saturating_sub(1));
        app.history = history;

        app.show_sources_panel = self.show_sources_panel;
        app.show_inspector_panel = self.show_inspector_panel;
        app.show_doc_table = self.show_doc_table;
        app.inspector_tab = self.inspector_tab;

        let mut display = self.display;
        display.col_width = sanitize_f32(display.col_width, 14.0, 4.0, 26.0);
        display.row_height = sanitize_f32(display.row_height, 22.0, 16.0, 30.0);
        display.translation_frame %= 3;
        app.display = display;

        let mut circular = self.circular;
        circular.rotation_rad = if circular.rotation_rad.is_finite() {
            circular.rotation_rad
        } else {
            0.0
        };
        circular.zoom = sanitize_f32(circular.zoom, 1.0, 0.4, 3.0);
        app.circular = circular;

        // Clamp cursor, scroll, and selection against the active document.
        let (width, rows) = app
            .current_alignment()
            .map(|a| (a.width(), a.row_count()))
            .unwrap_or((0, 0));
        app.horizontal_scroll_px = if self.horizontal_scroll_px.is_finite() {
            self.horizontal_scroll_px.max(0.0)
        } else {
            0.0
        };
        app.scroll_to_column = None;
        app.cursor_col = if width > 0 {
            self.cursor_col.min(width.saturating_sub(1))
        } else {
            0
        };
        app.cursor_row = if rows > 0 {
            self.cursor_row.min(rows.saturating_sub(1))
        } else {
            0
        };
        app.selected_cols = self.selected_cols.and_then(|(s, e)| {
            if width == 0 {
                return None;
            }
            let s = s.min(width.saturating_sub(1));
            let e = e.clamp(s + 1, width);
            Some((s, e))
        });
        app.selected_annotation_id = self.selected_annotation_id;

        app.motif_query = self.motif_query;
        app.motif_allow_iupac = self.motif_allow_iupac;
        app.annotation_filter = self.annotation_filter;
        if !self.uniprot_search_query.trim().is_empty() {
            app.uniprot_search_query = self.uniprot_search_query;
        }
        if !self.ncbi_search_query.trim().is_empty() {
            app.ncbi_search_query = self.ncbi_search_query;
        }

        // Convert any in-flight operations from a previous abrupt exit so spinners never stick.
        let mut ops: Vec<OperationRecord> = self
            .operations
            .into_iter()
            .map(|mut op| {
                if op.status == OperationStatus::Running {
                    op.status = OperationStatus::Failed;
                    op.detail = format!("{} (interrupted on previous exit)", op.detail);
                }
                op
            })
            .collect();
        let max_op_id = ops.iter().map(|o| o.id).max().unwrap_or(0);
        let next_op_id = self.next_op_id.max(max_op_id.saturating_add(1));
        ops.insert(
            0,
            OperationRecord {
                id: next_op_id,
                name: String::from("Session Restored"),
                status: OperationStatus::Completed,
                detail: format!("Restored {doc_count} workspace documents and view settings"),
                started_label: String::from("Startup"),
            },
        );
        if ops.len() > 80 {
            ops.truncate(80);
        }
        app.operations = ops;
        app.next_op_id = next_op_id.saturating_add(1);

        // Restore tool preferences.
        app.assemble_modal.mode = self.tool_prefs.assemble_mode.min(2);
        app.assemble_modal.algorithm = self.tool_prefs.assemble_algorithm;
        app.assemble_modal.minimap2_preset = self.tool_prefs.minimap2_preset;
        app.assemble_modal.minimap2_custom_path = self.tool_prefs.minimap2_custom_path;
        if self.tool_prefs.gap_open > 0 {
            app.assemble_modal.gap_open = self.tool_prefs.gap_open;
        }
        if self.tool_prefs.gap_extend > 0 {
            app.assemble_modal.gap_extend = self.tool_prefs.gap_extend;
        }
        if self.tool_prefs.min_overlap > 0 {
            app.assemble_modal.min_overlap = self.tool_prefs.min_overlap;
        }
        app.blast_modal.config = self.tool_prefs.blast_config;
        app.blast_modal.query = self.tool_prefs.blast_query;
        app.cloning_modal.method = self.tool_prefs.cloning_method;
        if !self.tool_prefs.cloning_enzyme_5prime.is_empty() {
            app.cloning_modal.enzyme_5prime = self.tool_prefs.cloning_enzyme_5prime;
        }
        if !self.tool_prefs.cloning_enzyme_3prime.is_empty() {
            app.cloning_modal.enzyme_3prime = self.tool_prefs.cloning_enzyme_3prime;
        }
        app.cloning_modal.reverse_insert = self.tool_prefs.cloning_reverse_insert;
        app.tree_modal.method = self.tool_prefs.tree_method;
        app.tree_modal.model = self.tool_prefs.tree_model;
        app.workflow_modal.step_auto_annotate = self.tool_prefs.workflow_auto_annotate;
        app.workflow_modal.step_realign = self.tool_prefs.workflow_realign;
        app.workflow_modal.step_translate = self.tool_prefs.workflow_translate;
        app.workflow_modal.step_build_tree = self.tool_prefs.workflow_build_tree;
        app.workflow_modal.step_design_primers = self.tool_prefs.workflow_design_primers;

        // Restore 3D Protein Viewer state (including separate OS window open state).
        let mut cam = self.viewer3d.camera;
        cam.yaw = if cam.yaw.is_finite() { cam.yaw } else { 0.55 };
        cam.pitch = sanitize_f32(cam.pitch, 0.25, -1.5, 1.5);
        cam.zoom = sanitize_f32(cam.zoom, 1.0, 0.2, 8.0);
        cam.pan_x = if cam.pan_x.is_finite() { cam.pan_x } else { 0.0 };
        cam.pan_y = if cam.pan_y.is_finite() { cam.pan_y } else { 0.0 };

        app.viewer3d.open = self.viewer3d.open;
        if !self.viewer3d.accession_input.trim().is_empty() {
            app.viewer3d.accession_input = self.viewer3d.accession_input;
        }
        app.viewer3d.source = self.viewer3d.source;
        app.viewer3d.linked_annotation_id = self.viewer3d.linked_annotation_id;
        app.viewer3d.camera = cam;
        app.viewer3d.style = self.viewer3d.style;
        app.viewer3d.color_by = self.viewer3d.color_by;
        app.viewer3d.structure = self.viewer3d.structure;
        app.viewer3d.alphafold_meta = self.viewer3d.alphafold_meta;
        app.viewer3d.uniprot_entry = self.viewer3d.uniprot_entry;
        app.viewer3d.selected_residue = self.viewer3d.selected_residue.filter(|&idx| {
            app.viewer3d
                .structure
                .as_ref()
                .is_some_and(|st| idx < st.residues.len())
        });
        app.rebuild_3d_alignment_map();

        if !app.motif_query.trim().is_empty() {
            app.run_motif_search();
        }

        let active_name = app
            .current_alignment()
            .map(|a| a.name.clone())
            .unwrap_or_else(|| String::from("Workspace"));
        app.status_message = format!("Restored previous session — {active_name}.");
        true
    }
}

/// Returns the per-user writable path where HelixPrime stores its session snapshot.
///
/// Uses `%LOCALAPPDATA%\HelixPrime\session.json` on Windows (falling back to `%APPDATA%`
/// and the user temp dir) and `$XDG_DATA_HOME/HelixPrime/session.json` on Unix/macOS,
/// so no administrator privileges are ever needed.
#[must_use]
pub fn default_session_path() -> PathBuf {
    if cfg!(target_os = "windows") {
        if let Some(local) = std::env::var_os("LOCALAPPDATA") {
            if !local.is_empty() {
                return PathBuf::from(local).join("HelixPrime").join("session.json");
            }
        }
        if let Some(roaming) = std::env::var_os("APPDATA") {
            if !roaming.is_empty() {
                return PathBuf::from(roaming).join("HelixPrime").join("session.json");
            }
        }
    } else {
        if let Some(xdg) = std::env::var_os("XDG_DATA_HOME") {
            if !xdg.is_empty() {
                return PathBuf::from(xdg).join("HelixPrime").join("session.json");
            }
        }
        if let Some(home) = std::env::var_os("HOME") {
            if !home.is_empty() {
                return PathBuf::from(home)
                    .join(".local")
                    .join("share")
                    .join("HelixPrime")
                    .join("session.json");
            }
        }
    }
    std::env::temp_dir().join("HelixPrime").join("session.json")
}

/// Loads a previously persisted session from disk, returning `None` if no file exists or
/// if it cannot be parsed.
#[must_use]
pub fn load_session_from_path(path: &Path) -> Option<PersistedSession> {
    let bytes = fs::read(path).ok()?;
    let session: PersistedSession = serde_json::from_slice(&bytes).ok()?;
    if session.schema_version != SESSION_SCHEMA_VERSION || session.docs.is_empty() {
        return None;
    }
    Some(session)
}

/// Loads the default per-user session file if present.
#[must_use]
pub fn load_default_session() -> Option<PersistedSession> {
    load_session_from_path(&default_session_path())
}

/// Atomically writes `session` to `path` via a sibling `.tmp` file so abrupt termination
/// cannot leave a truncated JSON file on disk.
pub fn save_session_to_path(session: &PersistedSession, path: &Path) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create session directory {}: {e}", parent.display()))?;
    }

    let payload = serde_json::to_vec(session)
        .map_err(|e| format!("failed to serialize session state: {e}"))?;

    let tmp_path = path.with_extension("json.tmp");
    fs::write(&tmp_path, &payload)
        .map_err(|e| format!("failed to write temporary session {}: {e}", tmp_path.display()))?;

    // Replace existing session file atomically; if antivirus/Windows indexing briefly locks
    // the rename target, fall back to a direct overwrite and clean up `.tmp`.
    if fs::rename(&tmp_path, path).is_err() {
        fs::write(path, &payload)
            .map_err(|e| format!("failed to save session {}: {e}", path.display()))?;
        let _ = fs::remove_file(&tmp_path);
    }

    Ok(())
}

/// Saves `session` into the default per-user session path.
pub fn save_default_session(session: &PersistedSession) -> Result<(), String> {
    save_session_to_path(session, &default_session_path())
}

/// Deletes the saved session file so the next startup (or manual reset) starts fresh.
pub fn clear_default_session() {
    let path = default_session_path();
    let _ = fs::remove_file(&path);
    let _ = fs::remove_file(path.with_extension("json.tmp"));
}

fn sanitize_f32(val: f32, fallback: f32, min: f32, max: f32) -> f32 {
    if val.is_finite() {
        val.clamp(min, max)
    } else {
        fallback
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use helix_core::{Alphabet, Sequence};

    #[test]
    fn session_round_trip_preserves_state_and_theme() {
        let aln = Alignment::new(
            "Test Plasmid",
            Alphabet::Dna,
            vec![Sequence::new("seq1", b"ATGCGATCGATCGATCGA")],
        )
        .expect("valid test alignment");

        let session = PersistedSession {
            theme: ThemeMode::Light,
            docs: vec![PersistedDoc {
                id: 42,
                folder: String::from("Cloning"),
                modified: String::from("Today"),
                is_circular: true,
                deleted: false,
                shared: true,
                checked: false,
                alignment: aln,
                tree: None,
                lineage: vec![LineageEntry {
                    operation: String::from("TOPO Cloning"),
                    parent_ids: vec![1, 2],
                    parent_names: vec![String::from("Vec"), String::from("Ins")],
                    details: String::from("Inserted fragment"),
                }],
            }],
            next_doc_id: 43,
            active_doc: 0,
            cursor_col: 7,
            horizontal_scroll_px: 128.5,
            ..PersistedSession::default()
        };

        let tmp_dir = std::env::temp_dir().join(format!(
            "helix_prime_persist_test_{}",
            std::process::id()
        ));
        let path = tmp_dir.join("session.json");

        save_session_to_path(&session, &path).expect("save succeeds");
        let loaded = load_session_from_path(&path).expect("load succeeds");

        assert_eq!(loaded.schema_version, SESSION_SCHEMA_VERSION);
        assert_eq!(loaded.theme, ThemeMode::Light);
        assert_eq!(loaded.docs.len(), 1);
        assert_eq!(loaded.docs[0].id, 42);
        assert_eq!(loaded.docs[0].alignment.name, "Test Plasmid");
        assert_eq!(loaded.docs[0].lineage.len(), 1);
        assert_eq!(loaded.cursor_col, 7);
        assert!((loaded.horizontal_scroll_px - 128.5).abs() < 1e-3);

        let _ = fs::remove_dir_all(&tmp_dir);
    }

    #[test]
    fn default_theme_is_dark() {
        assert_eq!(ThemeMode::default(), ThemeMode::Dark);
        assert_eq!(PersistedSession::default().theme, ThemeMode::Dark);
    }
}
