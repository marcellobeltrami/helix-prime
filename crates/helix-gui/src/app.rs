//! Top-level application state, workspace document database, lineage tracking,
//! Operations monitor, background NCBI/UniProt/AlphaFold worker orchestration,
//! and main `eframe::App` update loop including the separate OS 3D viewport window.

use crate::persistence;
use crate::renderer3d::{Camera3D, ColorBy3D, RenderStyle, Renderer3D};
use crate::sys_mem::{sample_memory_mb, MemorySnapshot};
use crate::theme::{apply_theme, ColorScheme, HighlightMode, Palette, ThemeMode};
use crate::widgets;
use egui::{Context, ViewportBuilder, ViewportId};
use serde::{Deserialize, Serialize};
use helix_core::aligners::{AlignerAlgorithm, Minimap2Preset};
use helix_core::blast::{BlastConfig, BlastHit};
use helix_core::cloning::CloningMethod;
use helix_core::database::{
    fetch_protein_bundle, fetch_rcsb_pdb, fetch_uniprot_entry, search_uniprot, AlphaFoldMeta,
    FetchedProteinBundle, ProteinStructureSource, UniProtEntry, UniProtSearchHit,
};
use helix_core::demo::{build_default_workspace, build_offline_preview_structure};
use helix_core::io::{
    clustal::write_clustal, fasta::write_fasta, gff3::import_gff3_annotations,
    load_alignment_from_path, mafft::load_mafft_alignment_from_path,
};
use helix_core::ncbi::{download_ncbi_record, search_ncbi, NcbiDatabase, NcbiSearchHit};
use helix_core::phylo::{build_phylogenetic_tree, DistanceModel, PhyloTree, TreeMethod};
use helix_core::primers::{DesignedPrimerPair, PrimerBindingHit, PrimerDesignConfig};
use helix_core::restriction::{find_restriction_sites, RestrictionSiteHit, COMMON_ENZYMES};
use helix_core::{
    parse_mmcif, parse_pdb, Alignment, Annotation, FeatureKind, MotifHit, ProteinStructure, Strand,
    StructureAlignmentMap, UniProtLink,
};
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver, Sender};
use std::time::Instant;

/// Tracks how a derived document was created from one or more parent documents (Geneious Lineage tab).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageEntry {
    pub operation: String,
    pub parent_ids: Vec<u64>,
    pub parent_names: Vec<String>,
    pub details: String,
}

/// A single biological document inside the HelixPrime workspace database.
pub struct WorkspaceDoc {
    pub id: u64,
    pub folder: String,
    pub modified: String,
    pub is_circular: bool,
    pub deleted: bool,
    pub shared: bool,
    pub checked: bool,
    pub alignment: Alignment,
    /// Optional phylogenetic tree when this document is a Tree Document.
    pub tree: Option<PhyloTree>,
    /// Derivation history shown in the Lineage tab.
    pub lineage: Vec<LineageEntry>,
    /// Cached N x N pairwise identity matrix (recomputed lazily).
    pub pairwise_cache: Option<Vec<f32>>,
    /// Cached restriction enzyme cut sites for the reference sequence.
    pub restriction_cache: Option<Vec<RestrictionSiteHit>>,
}

impl WorkspaceDoc {
    pub fn new(
        id: u64,
        folder: impl Into<String>,
        modified: impl Into<String>,
        is_circular: bool,
        alignment: Alignment,
    ) -> Self {
        Self {
            id,
            folder: folder.into(),
            modified: modified.into(),
            is_circular,
            deleted: false,
            shared: false,
            checked: false,
            alignment,
            tree: None,
            lineage: Vec::new(),
            pairwise_cache: None,
            restriction_cache: None,
        }
    }

    pub fn invalidate_caches(&mut self) {
        self.pairwise_cache = None;
        self.restriction_cache = None;
    }

    /// Ensures the restriction site cache is populated and returns a shared slice.
    pub fn ensure_restriction_sites(&mut self) -> &[RestrictionSiteHit] {
        if self.restriction_cache.is_none() {
            let sites = if self.alignment.alphabet.is_nucleotide() {
                self.alignment
                    .sequences()
                    .first()
                    .map(|s| find_restriction_sites(&s.residues, COMMON_ENZYMES, self.is_circular))
                    .unwrap_or_default()
            } else {
                Vec::new()
            };
            self.restriction_cache = Some(sites);
        }
        self.restriction_cache.as_deref().unwrap_or(&[])
    }
}

/// Active selection in the left Geneious Prime Sources tree.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SourceSelection {
    AllLocal,
    Folder(String),
    ReferenceFeatures,
    DeletedItems,
    SharedDatabases,
    Operations,
    Ncbi(NcbiDatabase),
    UniProt,
}

/// Active tab inside the lower Document Viewer pane (matching Geneious Prime).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum DocTab {
    #[default]
    SequenceView,
    Annotations,
    TextView,
    Lineage,
    Info,
    TreeView,
}

/// Right-side inspector drawer tabs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum InspectorTab {
    #[default]
    Display,
    Annotations,
    FindMotif,
    Statistics,
}

/// Interactive camera and track settings for the Circular Overview / Phage Genome ring.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct CircularViewState {
    pub rotation_rad: f32,
    pub zoom: f32,
    pub show_labels: bool,
    pub show_sequence_rings: bool,
    pub show_identity_ring: bool,
    pub show_gc_ring: bool,
}

impl Default for CircularViewState {
    fn default() -> Self {
        Self {
            rotation_rad: 0.0,
            zoom: 1.0,
            show_labels: true,
            show_sequence_rings: true,
            show_identity_ring: true,
            show_gc_ring: true,
        }
    }
}

/// Visual configuration for the main Sequence View (matching Geneious Prime's General sidebar).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct DisplaySettings {
    pub color_scheme: ColorScheme,
    pub highlight_mode: HighlightMode,
    pub dots_for_matches: bool,
    pub show_overview: bool,
    pub show_identity_graph: bool,
    pub show_consensus: bool,
    pub show_annotations: bool,
    pub show_complement: bool,
    pub show_translation: bool,
    pub show_restriction_sites: bool,
    pub translation_frame: usize,
    /// Show the Circular Overview / Phage Genome ring.
    pub show_circular_overview: bool,
    /// Show the zoomed Linear sequence/alignment canvas.
    pub show_linear_view: bool,
    /// Wrap the Linear view into multi-row blocks (60–100 columns per row).
    pub wrap_sequence: bool,
    /// Allow direct keyboard base editing in Sequence View.
    pub allow_editing: bool,
    /// Width of one residue/base column in logical pixels (`4.0..=26.0`).
    pub col_width: f32,
    /// Height of one sequence row in logical pixels (`16.0..=30.0`).
    pub row_height: f32,
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            color_scheme: ColorScheme::Geneious,
            highlight_mode: HighlightMode::AllResidues,
            dots_for_matches: false,
            show_overview: true,
            show_identity_graph: true,
            show_consensus: true,
            show_annotations: true,
            show_complement: false,
            show_translation: false,
            show_restriction_sites: false,
            translation_frame: 0,
            show_circular_overview: true,
            show_linear_view: true,
            wrap_sequence: false,
            allow_editing: false,
            col_width: 14.0,
            row_height: 22.0,
        }
    }
}

/// Background operation lifecycle status shown in the Operations panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OperationStatus {
    Running,
    Completed,
    Failed,
}

/// A single tracked operation in the Operations log.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OperationRecord {
    pub id: u64,
    pub name: String,
    pub status: OperationStatus,
    pub detail: String,
    pub started_label: String,
}

/// Draft state for the "Add Annotation" modal dialog.
#[derive(Debug, Clone)]
pub struct NewAnnotationDraft {
    pub open: bool,
    pub name: String,
    pub kind: FeatureKind,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Strand,
    pub uniprot_accession: String,
    pub notes: String,
}

impl Default for NewAnnotationDraft {
    fn default() -> Self {
        Self {
            open: false,
            name: String::from("New Feature"),
            kind: FeatureKind::Domain,
            start_1based: 1,
            end_1based: 20,
            strand: Strand::Forward,
            uniprot_accession: String::new(),
            notes: String::new(),
        }
    }
}

/// Modal state for the In-Silico Cloning simulator.
#[derive(Debug, Clone)]
pub struct CloningModalState {
    pub open: bool,
    pub method: CloningMethod,
    pub vector_doc_idx: usize,
    pub insert_doc_idx: usize,
    pub enzyme_5prime: String,
    pub enzyme_3prime: String,
    pub reverse_insert: bool,
    pub construct_name: String,
    pub error: Option<String>,
}

impl Default for CloningModalState {
    fn default() -> Self {
        Self {
            open: false,
            method: CloningMethod::TopoTa,
            vector_doc_idx: 0,
            insert_doc_idx: 9,
            enzyme_5prime: String::from("EcoRI"),
            enzyme_3prime: String::from("HindIII"),
            reverse_insert: false,
            construct_name: String::from("pCRII-TOPO::eGFP"),
            error: None,
        }
    }
}

/// Modal state for PCR Primer Design & Binding Site Testing.
#[derive(Debug, Clone)]
pub struct PrimersModalState {
    pub open: bool,
    pub config: PrimerDesignConfig,
    pub designed_pairs: Vec<DesignedPrimerPair>,
    pub test_primer_seq: String,
    pub test_max_mismatches: usize,
    pub binding_hits: Vec<PrimerBindingHit>,
}

impl Default for PrimersModalState {
    fn default() -> Self {
        Self {
            open: false,
            config: PrimerDesignConfig::default(),
            designed_pairs: Vec::new(),
            test_primer_seq: String::from("ATGGTGCATCTGACTCCTGAGG"),
            test_max_mismatches: 1,
            binding_hits: Vec::new(),
        }
    }
}

/// Modal state for Align / Map to Reference / De Novo Contig Assembly.
#[derive(Debug, Clone)]
pub struct AssembleModalState {
    pub open: bool,
    /// 0 = Multiple Sequence Alignment (De Bruijn / ClustalW / Clustal Omega / MUSCLE / MAFFT / minimap2),
    /// 1 = De Novo Contig Assembly, 2 = Map Reads to Reference.
    pub mode: usize,
    pub algorithm: AlignerAlgorithm,
    pub minimap2_preset: Minimap2Preset,
    pub minimap2_custom_path: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub reference_doc_idx: usize,
    pub reads_doc_idx: usize,
    pub min_overlap: usize,
    pub output_name: String,
    pub last_note: Option<String>,
    pub error: Option<String>,
}

impl Default for AssembleModalState {
    fn default() -> Self {
        Self {
            open: false,
            mode: 0,
            algorithm: AlignerAlgorithm::ClustalW,
            minimap2_preset: Minimap2Preset::Asm5,
            minimap2_custom_path: String::new(),
            gap_open: 8,
            gap_extend: 2,
            reference_doc_idx: 6,
            reads_doc_idx: 8,
            min_overlap: 12,
            output_name: String::from("Aligned_Result"),
            last_note: None,
            error: None,
        }
    }
}

/// Modal state for Phylogenetic Tree builder.
#[derive(Debug, Clone)]
pub struct TreeModalState {
    pub open: bool,
    pub method: TreeMethod,
    pub model: DistanceModel,
    pub error: Option<String>,
}

impl Default for TreeModalState {
    fn default() -> Self {
        Self {
            open: false,
            method: TreeMethod::NeighborJoining,
            model: DistanceModel::JukesCantor,
            error: None,
        }
    }
}

/// Modal state for Multi-Program BLAST (`totalBLAST`, `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `megablast`).
#[derive(Debug, Clone, Default)]
pub struct BlastModalState {
    pub open: bool,
    pub config: BlastConfig,
    pub query: String,
    pub hits: Vec<BlastHit>,
    pub selected_hit: Option<usize>,
}

/// Modal state for Automated Workflows pipeline builder.
#[derive(Debug, Clone)]
pub struct WorkflowModalState {
    pub open: bool,
    pub step_auto_annotate: bool,
    pub step_realign: bool,
    pub step_translate: bool,
    pub step_build_tree: bool,
    pub step_design_primers: bool,
    pub last_report: Option<String>,
}

impl Default for WorkflowModalState {
    fn default() -> Self {
        Self {
            open: false,
            step_auto_annotate: true,
            step_realign: false,
            step_translate: false,
            step_build_tree: false,
            step_design_primers: true,
            last_report: None,
        }
    }
}

/// State backing the separate native OS 3D Protein Viewer window.
pub struct ProteinViewerState {
    pub open: bool,
    pub accession_input: String,
    /// Selected remote backend for 3D structure & metadata fetching (AlphaFold DB vs UniProtKB).
    pub source: ProteinStructureSource,
    pub structure: Option<ProteinStructure>,
    pub alphafold_meta: Option<AlphaFoldMeta>,
    pub uniprot_entry: Option<UniProtEntry>,
    pub alignment_map: Option<StructureAlignmentMap>,
    pub linked_annotation_id: Option<u64>,
    pub camera: Camera3D,
    pub style: RenderStyle,
    pub color_by: ColorBy3D,
    pub hovered_residue: Option<usize>,
    pub selected_residue: Option<usize>,
    pub loading_message: Option<String>,
    pub error_message: Option<String>,
    pub renderer: Renderer3D,
}

impl Default for ProteinViewerState {
    fn default() -> Self {
        Self {
            open: false,
            accession_input: String::from("P68871"),
            source: ProteinStructureSource::AutoFallback,
            structure: None,
            alphafold_meta: None,
            uniprot_entry: None,
            alignment_map: None,
            linked_annotation_id: None,
            camera: Camera3D::default(),
            style: RenderStyle::Cartoon,
            color_by: ColorBy3D::AlphaFoldPlddt,
            hovered_residue: None,
            selected_residue: None,
            loading_message: None,
            error_message: None,
            renderer: Renderer3D::default(),
        }
    }
}

/// Asynchronous events delivered from background network worker threads to the GUI thread.
pub enum AsyncEvent {
    UniProtSearchDone {
        op_id: u64,
        result: Result<Vec<UniProtSearchHit>, String>,
    },
    UniProtEntryLoaded {
        op_id: u64,
        entry: Result<UniProtEntry, String>,
        import_annotations_for_doc: Option<usize>,
        attach_to_annotation: Option<(usize, u64)>,
        save_as_new_doc: bool,
    },
    ProteinBundleLoaded {
        op_id: u64,
        // Boxed so `AsyncEvent` stays compact on the channel stack (avoids 600+ byte enum padding).
        result: Result<Box<FetchedProteinBundle>, String>,
        linked_annotation_id: Option<u64>,
    },
    RcsbPdbLoaded {
        op_id: u64,
        result: Result<ProteinStructure, String>,
    },
    NcbiSearchDone {
        op_id: u64,
        db: NcbiDatabase,
        result: Result<Vec<NcbiSearchHit>, String>,
    },
    NcbiRecordDownloaded {
        op_id: u64,
        db: NcbiDatabase,
        result: Result<Alignment, String>,
    },
}

/// Main HelixPrime desktop application state.
pub struct HelixApp {
    pub theme: ThemeMode,
    pub docs: Vec<WorkspaceDoc>,
    pub next_doc_id: u64,
    pub active_doc: usize,
    pub source_selection: SourceSelection,
    pub active_tab: DocTab,
    pub global_search: String,

    // Back / Forward document history.
    pub history: Vec<usize>,
    pub history_pos: usize,

    // Panel visibility & inspector.
    pub show_sources_panel: bool,
    pub show_inspector_panel: bool,
    pub show_doc_table: bool,
    pub inspector_tab: InspectorTab,
    pub display: DisplaySettings,
    pub circular: CircularViewState,

    // Sequence/alignment cursor, selection & cross-window sync.
    pub horizontal_scroll_px: f32,
    pub scroll_to_column: Option<usize>,
    pub cursor_col: usize,
    pub cursor_row: usize,
    pub selected_cols: Option<(usize, usize)>,
    pub selected_annotation_id: Option<u64>,
    pub hovered_column: Option<usize>,
    pub hovered_row: Option<usize>,

    // Motif search state.
    pub motif_query: String,
    pub motif_allow_iupac: bool,
    pub motif_hits: Vec<MotifHit>,
    pub active_motif_hit: usize,

    // Annotations tab filter string.
    pub annotation_filter: String,

    // UniProt & NCBI live database search state.
    pub uniprot_search_query: String,
    pub uniprot_search_results: Vec<UniProtSearchHit>,
    pub uniprot_searching: bool,
    pub uniprot_search_error: Option<String>,

    pub ncbi_search_query: String,
    pub ncbi_search_results: Vec<NcbiSearchHit>,
    pub ncbi_searching: bool,
    pub ncbi_search_error: Option<String>,

    // Operations & Memory telemetry.
    pub operations: Vec<OperationRecord>,
    pub next_op_id: u64,
    pub memory: MemorySnapshot,
    pub last_mem_sample: Instant,

    // Tool dialogs.
    pub new_ann_draft: NewAnnotationDraft,
    pub cloning_modal: CloningModalState,
    pub primers_modal: PrimersModalState,
    pub assemble_modal: AssembleModalState,
    pub tree_modal: TreeModalState,
    pub blast_modal: BlastModalState,
    pub workflow_modal: WorkflowModalState,

    // Separate OS window 3D protein viewer.
    pub viewer3d: ProteinViewerState,

    // Background worker channel.
    pub event_tx: Sender<AsyncEvent>,
    pub event_rx: Receiver<AsyncEvent>,

    /// Timestamp of the last disk session flush to prevent duplicate writes during shutdown.
    pub last_session_save: Option<Instant>,

    pub status_message: String,
}

impl HelixApp {
    /// Initializes HelixPrime, restoring the previous user session if present or starting
    /// fresh in Dark Mode with the built-in Geneious Prime-inspired workspace.
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self::build_fresh_workspace();

        if let Some(saved) = persistence::load_default_session() {
            let restored = saved.apply_to(&mut app);
            if restored && app.viewer3d.open && app.viewer3d.structure.is_none() {
                // Re-populate the 3D viewport when a very large structure wasn't embedded inline.
                let acc = app.viewer3d.accession_input.clone();
                if !acc.trim().is_empty() {
                    app.open_3d_for_accession(&cc.egui_ctx, &acc, app.viewer3d.linked_annotation_id);
                } else {
                    app.open_offline_3d_preview();
                }
            }
        }

        apply_theme(&cc.egui_ctx, app.theme);
        app
    }

    /// Builds a clean default workspace starting in Dark Mode (pCRII-TOPO, pDEST15, pDONR201,
    /// BioBrick, HBB/TP53 alignments, circular Phage genome, primers, Sanger reads, and tree).
    fn build_fresh_workspace() -> Self {
        let (event_tx, event_rx) = mpsc::channel();
        let mut next_doc_id = 1u64;
        let mut docs: Vec<WorkspaceDoc> = build_default_workspace()
            .into_iter()
            .map(|d| {
                let id = next_doc_id;
                next_doc_id += 1;
                WorkspaceDoc::new(id, d.folder, d.modified, d.is_circular, d.alignment)
            })
            .collect();

        // Pre-build a sample Phylogenetic Tree document from the HBB Orthologs alignment so
        // `Tree Documents` is populated immediately!
        if let Some(hbb_doc) = docs
            .iter()
            .find(|d| d.alignment.name.contains("Hemoglobin Beta Orthologs"))
        {
            let parent_id = hbb_doc.id;
            let parent_name = hbb_doc.alignment.name.clone();
            if let Ok(tree) = build_phylogenetic_tree(
                &hbb_doc.alignment,
                TreeMethod::NeighborJoining,
                DistanceModel::JukesCantor,
            ) {
                let mut tree_doc = WorkspaceDoc::new(
                    next_doc_id,
                    "Tree Documents",
                    "19 Nov 2026",
                    false,
                    hbb_doc.alignment.clone().with_description(
                        "Neighbor-Joining phylogenetic tree built from Hemoglobin Beta Orthologs",
                    ),
                );
                tree_doc.alignment.name = String::from("HBB Orthologs Tree (NJ)");
                tree_doc.tree = Some(tree);
                tree_doc.lineage.push(LineageEntry {
                    operation: String::from("Build Phylogenetic Tree (Neighbor-Joining / Jukes-Cantor)"),
                    parent_ids: vec![parent_id],
                    parent_names: vec![parent_name],
                    details: String::from("Inferred additive tree from pairwise distance matrix"),
                });
                next_doc_id += 1;
                docs.push(tree_doc);
            }
        }

        let mut app = Self {
            theme: ThemeMode::Dark,
            docs,
            next_doc_id,
            active_doc: 0, // Start on pCRII-TOPO in Cloning (split Circular + Linear view!)
            source_selection: SourceSelection::Folder(String::from("Cloning")),
            active_tab: DocTab::SequenceView,
            global_search: String::new(),
            history: vec![0],
            history_pos: 0,
            show_sources_panel: true,
            show_inspector_panel: true,
            show_doc_table: true,
            inspector_tab: InspectorTab::Display,
            display: DisplaySettings::default(),
            circular: CircularViewState::default(),
            horizontal_scroll_px: 0.0,
            scroll_to_column: Some(1178), // Focus around bp 1,179 like the reference screenshot!
            cursor_col: 1178,
            cursor_row: 0,
            selected_cols: None,
            selected_annotation_id: None,
            hovered_column: None,
            hovered_row: None,
            motif_query: String::new(),
            motif_allow_iupac: true,
            motif_hits: Vec::new(),
            active_motif_hit: 0,
            annotation_filter: String::new(),
            uniprot_search_query: String::from("TP53 human"),
            uniprot_search_results: Vec::new(),
            uniprot_searching: false,
            uniprot_search_error: None,
            ncbi_search_query: String::from("Escherichia phage lambda"),
            ncbi_search_results: Vec::new(),
            ncbi_searching: false,
            ncbi_search_error: None,
            operations: vec![OperationRecord {
                id: 1,
                name: String::from("Workspace Initialized"),
                status: OperationStatus::Completed,
                detail: String::from("Loaded Sample Documents, Reference Plasmid Features, and 3D AlphaFold Engine"),
                started_label: String::from("Startup"),
            }],
            next_op_id: 2,
            memory: sample_memory_mb(),
            last_mem_sample: Instant::now(),
            new_ann_draft: NewAnnotationDraft::default(),
            cloning_modal: CloningModalState::default(),
            primers_modal: PrimersModalState::default(),
            assemble_modal: AssembleModalState::default(),
            tree_modal: TreeModalState::default(),
            blast_modal: BlastModalState::default(),
            workflow_modal: WorkflowModalState::default(),
            viewer3d: ProteinViewerState::default(),
            event_tx,
            event_rx,
            last_session_save: None,
            status_message: String::from(
                "Ready — Viewing pCRII-TOPO (Circular Overview + Linear Zoom). Session auto-saves on exit.",
            ),
        };

        if let Some(first_ann) = app
            .current_alignment()
            .and_then(|a| a.annotations().first().map(|ann| ann.id))
        {
            app.selected_annotation_id = Some(first_ann);
        }

        app
    }

    /// Persists the current workspace state to disk, updating `last_session_save`.
    pub fn persist_session(&mut self, force: bool) -> Result<(), String> {
        if !force {
            if let Some(last) = self.last_session_save {
                if last.elapsed().as_millis() < 500 {
                    return Ok(());
                }
            }
        }
        let snapshot = persistence::PersistedSession::from_app(self);
        persistence::save_default_session(&snapshot)?;
        self.last_session_save = Some(Instant::now());
        Ok(())
    }

    /// Explicitly saves the current session and reports status in the footer bar.
    pub fn save_session_now(&mut self) {
        match self.persist_session(true) {
            Ok(()) => {
                self.status_message = format!(
                    "Saved workspace session ({} documents).",
                    self.docs.len()
                );
            }
            Err(err) => {
                self.status_message = format!("Session save failed: {err}");
            }
        }
    }

    /// Clears saved session state on disk and resets the workspace to factory defaults (Dark Mode).
    pub fn reset_workspace_to_defaults(&mut self, ctx: &Context) {
        persistence::clear_default_session();
        *self = Self::build_fresh_workspace();
        apply_theme(ctx, self.theme);
        self.status_message = String::from("Reset workspace to default sample documents (Dark Mode).");
    }

    /// Updates Light/Dark theme and applies styling immediately.
    pub fn set_theme(&mut self, ctx: &Context, mode: ThemeMode) {
        self.theme = mode;
        apply_theme(ctx, mode);
    }

    /// Registers a new tracked operation in the Operations log and returns its unique ID.
    pub fn start_operation(&mut self, name: impl Into<String>, detail: impl Into<String>) -> u64 {
        let id = self.next_op_id;
        self.next_op_id = self.next_op_id.saturating_add(1);
        self.operations.insert(
            0,
            OperationRecord {
                id,
                name: name.into(),
                status: OperationStatus::Running,
                detail: detail.into(),
                started_label: format!("Job #{id}"),
            },
        );
        if self.operations.len() > 80 {
            self.operations.truncate(80);
        }
        id
    }

    /// Marks a tracked operation completed or failed.
    pub fn finish_operation(&mut self, id: u64, status: OperationStatus, detail: impl Into<String>) {
        if let Some(op) = self.operations.iter_mut().find(|o| o.id == id) {
            op.status = status;
            op.detail = detail.into();
        }
    }

    /// Appends a new document to the workspace, records optional Lineage, selects it, and returns its index.
    pub fn add_workspace_doc(
        &mut self,
        folder: impl Into<String>,
        is_circular: bool,
        alignment: Alignment,
        tree: Option<PhyloTree>,
        lineage: Option<LineageEntry>,
    ) -> usize {
        let id = self.next_doc_id;
        self.next_doc_id = self.next_doc_id.saturating_add(1);
        let folder_str = folder.into();
        let mut doc = WorkspaceDoc::new(id, folder_str.clone(), "Today", is_circular, alignment);
        doc.tree = tree;
        if let Some(lin) = lineage {
            doc.lineage.push(lin);
        }
        self.docs.push(doc);
        let idx = self.docs.len() - 1;
        self.source_selection = SourceSelection::Folder(folder_str);
        self.select_document(idx);
        idx
    }

    #[must_use]
    pub fn current_doc(&self) -> Option<&WorkspaceDoc> {
        self.docs.get(self.active_doc)
    }

    #[must_use]
    pub fn current_alignment(&self) -> Option<&Alignment> {
        self.docs.get(self.active_doc).map(|d| &d.alignment)
    }

    pub fn current_doc_mut(&mut self) -> Option<&mut WorkspaceDoc> {
        self.docs.get_mut(self.active_doc)
    }

    /// Selects a document by index and records it in Back/Forward navigation history.
    pub fn select_document(&mut self, idx: usize) {
        self.select_document_internal(idx, true);
    }

    fn select_document_internal(&mut self, idx: usize, push_history: bool) {
        if idx >= self.docs.len() {
            return;
        }
        self.active_doc = idx;
        if push_history {
            self.history.truncate(self.history_pos + 1);
            if self.history.last().copied() != Some(idx) {
                self.history.push(idx);
                self.history_pos = self.history.len() - 1;
            }
        }

        let is_tree = self.docs[idx].tree.is_some();
        let is_circular = self.docs[idx].is_circular;
        let is_phage = self.docs[idx].alignment.name.contains("Bacteriophage");

        if is_tree {
            self.active_tab = DocTab::TreeView;
        } else if self.active_tab == DocTab::TreeView {
            self.active_tab = DocTab::SequenceView;
        }

        // Adapt default view layout naturally: circular plasmids open Split Circular+Linear,
        // bacteriophage circular genome opens Circular Phage overview, linear alignments open Linear.
        if is_phage {
            self.display.show_circular_overview = true;
            self.display.show_linear_view = false;
        } else if is_circular {
            self.display.show_circular_overview = true;
            self.display.show_linear_view = true;
        } else {
            self.display.show_circular_overview = false;
            self.display.show_linear_view = true;
        }

        self.horizontal_scroll_px = 0.0;
        self.cursor_col = 0;
        self.cursor_row = 0;
        self.selected_cols = None;
        self.hovered_column = None;
        self.motif_hits.clear();
        self.active_motif_hit = 0;
        self.selected_annotation_id = self
            .current_alignment()
            .and_then(|a| a.annotations().first().map(|ann| ann.id));

        if let Some(aln) = self.current_alignment() {
            let unit = if aln.alphabet.is_nucleotide() { "bp" } else { "aa" };
            self.status_message = format!(
                "Opened '{}' ({} seqs, {} {unit}, {:.1}% pairwise ID)",
                aln.name,
                aln.row_count(),
                aln.width(),
                aln.mean_pairwise_identity_pct()
            );
        }
        self.rebuild_3d_alignment_map();
    }

    pub fn navigate_back(&mut self) {
        if self.history_pos > 0 {
            self.history_pos -= 1;
            let idx = self.history[self.history_pos];
            self.select_document_internal(idx, false);
        }
    }

    pub fn navigate_forward(&mut self) {
        if self.history_pos + 1 < self.history.len() {
            self.history_pos += 1;
            let idx = self.history[self.history_pos];
            self.select_document_internal(idx, false);
        }
    }

    /// Recomputes the bidirectional coordinate map between the loaded 3D structure and the active alignment.
    pub fn rebuild_3d_alignment_map(&mut self) {
        let Some(ref structure) = self.viewer3d.structure else {
            self.viewer3d.alignment_map = None;
            return;
        };
        let Some(aln) = self.current_alignment() else {
            self.viewer3d.alignment_map = None;
            return;
        };
        let linked_ann = self
            .viewer3d
            .linked_annotation_id
            .and_then(|id| aln.annotation_by_id(id))
            .or_else(|| {
                self.selected_annotation_id
                    .and_then(|id| aln.annotation_by_id(id))
            });
        self.viewer3d.alignment_map =
            Some(StructureAlignmentMap::build(structure, aln, linked_ann));
    }

    /// Triggers an asynchronous 3D protein structure & metadata fetch from the selected source
    /// (`AlphaFold DB`, `UniProtKB`, or `AutoFallback`) and opens the separate 3D OS window.
    pub fn open_3d_for_accession(
        &mut self,
        ctx: &Context,
        accession: &str,
        linked_annotation_id: Option<u64>,
    ) {
        let source = self.viewer3d.source;
        self.open_3d_with_source(ctx, accession, source, linked_annotation_id);
    }

    /// Fetches a 3D protein structure + metadata explicitly from `source` (`AlphaFoldDb`,
    /// `UniProtDb`, or `AutoFallback`).
    pub fn open_3d_with_source(
        &mut self,
        ctx: &Context,
        accession: &str,
        source: ProteinStructureSource,
        linked_annotation_id: Option<u64>,
    ) {
        let clean = accession.trim().to_ascii_uppercase();
        if clean.is_empty() {
            self.status_message = String::from("Enter a valid UniProt accession (e.g. P68871 or P04637).");
            return;
        }

        let op_id = self.start_operation(
            format!("3D Fetch ({clean})"),
            format!("Source: {}", source.label()),
        );

        self.viewer3d.open = true;
        self.viewer3d.source = source;
        self.viewer3d.accession_input = clean.clone();
        self.viewer3d.linked_annotation_id = linked_annotation_id;
        self.viewer3d.error_message = None;
        self.viewer3d.loading_message = Some(format!(
            "Fetching 3D model & metadata for {clean} via {}...",
            source.label()
        ));
        self.status_message = format!("Querying {} for {clean}...", source.label());

        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let result = fetch_protein_bundle(&clean, source)
                .map(Box::new)
                .map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::ProteinBundleLoaded {
                op_id,
                result,
                linked_annotation_id,
            });
            ctx_clone.request_repaint();
        });
    }

    pub fn open_offline_3d_preview(&mut self) {
        let structure = build_offline_preview_structure();
        self.viewer3d.open = true;
        self.viewer3d.accession_input = String::from("P68871");
        self.viewer3d.structure = Some(structure);
        self.viewer3d.alphafold_meta = None;
        self.viewer3d.loading_message = None;
        self.viewer3d.error_message = None;
        self.viewer3d.camera.reset();
        self.rebuild_3d_alignment_map();
        self.status_message = String::from("Loaded offline HBB 3D schematic preview in separate window.");
    }

    pub fn trigger_uniprot_search(&mut self, ctx: &Context) {
        let q = self.uniprot_search_query.trim().to_owned();
        if q.is_empty() {
            return;
        }
        let op_id = self.start_operation("UniProt Search", format!("Query: '{q}'"));
        self.uniprot_searching = true;
        self.uniprot_search_error = None;
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let result = search_uniprot(&q, 15).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::UniProtSearchDone { op_id, result });
            ctx_clone.request_repaint();
        });
    }

    pub fn trigger_ncbi_search(&mut self, ctx: &Context, db: NcbiDatabase) {
        let q = self.ncbi_search_query.trim().to_owned();
        if q.is_empty() {
            return;
        }
        let op_id = self.start_operation(
            format!("NCBI {} Search", db.label()),
            format!("Query: '{q}'"),
        );
        self.ncbi_searching = true;
        self.ncbi_search_error = None;
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let result = search_ncbi(db, &q, 15).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::NcbiSearchDone { op_id, db, result });
            ctx_clone.request_repaint();
        });
    }

    pub fn download_ncbi_hit(&mut self, ctx: &Context, hit: &NcbiSearchHit) {
        let db = hit.database;
        let uid = hit.uid.clone();
        let name = hit.accession.clone();
        let op_id = self.start_operation(
            format!("NCBI Download ({name})"),
            format!("Fetching from NCBI {}...", db.label()),
        );
        self.status_message = format!("Downloading {name} from NCBI {}...", db.label());
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let result = download_ncbi_record(db, &uid, &name).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::NcbiRecordDownloaded { op_id, db, result });
            ctx_clone.request_repaint();
        });
    }

    pub fn download_uniprot_as_doc(&mut self, ctx: &Context, accession: &str) {
        let clean = accession.trim().to_ascii_uppercase();
        if clean.is_empty() {
            return;
        }
        let op_id = self.start_operation(
            format!("UniProt Download ({clean})"),
            "Fetching protein sequence & curated features",
        );
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let entry = fetch_uniprot_entry(&clean).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::UniProtEntryLoaded {
                op_id,
                entry,
                import_annotations_for_doc: None,
                attach_to_annotation: None,
                save_as_new_doc: true,
            });
            ctx_clone.request_repaint();
        });
    }

    pub fn fetch_and_import_uniprot_annotations(
        &mut self,
        ctx: &Context,
        accession: &str,
        attach_to_ann_id: Option<u64>,
    ) {
        let clean = accession.trim().to_ascii_uppercase();
        if clean.is_empty() {
            return;
        }
        let doc_idx = self.active_doc;
        let op_id = self.start_operation(
            format!("UniProt Features ({clean})"),
            "Importing domains & active sites",
        );
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let entry = fetch_uniprot_entry(&clean).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::UniProtEntryLoaded {
                op_id,
                entry,
                import_annotations_for_doc: Some(doc_idx),
                attach_to_annotation: attach_to_ann_id.map(|id| (doc_idx, id)),
                save_as_new_doc: false,
            });
            ctx_clone.request_repaint();
        });
    }

    pub fn fetch_rcsb_structure(&mut self, ctx: &Context, pdb_id: &str) {
        let clean = pdb_id.trim().to_ascii_uppercase();
        let acc = self
            .viewer3d
            .uniprot_entry
            .as_ref()
            .map(|e| e.accession.clone());
        let op_id = self.start_operation(format!("RCSB PDB ({clean})"), "Downloading experimental PDB");
        self.viewer3d.open = true;
        self.viewer3d.error_message = None;
        self.viewer3d.loading_message = Some(format!("Downloading experimental RCSB PDB {clean}..."));
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let result = fetch_rcsb_pdb(&clean, acc).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::RcsbPdbLoaded { op_id, result });
            ctx_clone.request_repaint();
        });
    }

    pub fn import_alignment_dialog(&mut self) {
        let picked: Option<PathBuf> = rfd::FileDialog::new()
            .add_filter(
                "Bioinformatics Sequences & Alignments (FASTA, MAFFT, Clustal, PHYLIP, GenBank)",
                &[
                    "fasta", "fa", "fna", "faa", "fas", "mafft", "afa", "mfa", "msa", "aln",
                    "clustal", "clw", "phy", "phylip", "gb", "gbk", "genbank",
                ],
            )
            .add_filter(
                "Pre-Aligned MAFFT Files",
                &["mafft", "afa", "mfa", "msa", "aln", "clw", "phy", "phylip", "fasta", "fa"],
            )
            .pick_file();

        let Some(path) = picked else {
            return;
        };

        match load_alignment_from_path(&path) {
            Ok(aln) => {
                let name = aln.name.clone();
                let folder = if aln.row_count() > 1 {
                    "Alignments"
                } else {
                    "Imported Files"
                };
                self.add_workspace_doc(folder, false, aln, None, None);
                self.status_message = format!("Imported '{name}' into {folder} from {}", path.display());
            }
            Err(err) => {
                self.status_message = format!("Import failed: {err}");
            }
        }
    }

    /// Imports an already-aligned MAFFT file (aligned FASTA, `--clustalout`, or `--phylipout`)
    /// directly into the `Alignments` folder while preserving all gap columns verbatim.
    pub fn import_mafft_dialog(&mut self) {
        let picked: Option<PathBuf> = rfd::FileDialog::new()
            .set_title("Import Pre-Aligned MAFFT Alignment")
            .add_filter(
                "MAFFT Pre-Aligned Files (.mafft, .afa, .aln, .phy, .fasta)",
                &[
                    "mafft", "afa", "mfa", "msa", "aln", "clustal", "clw", "phy", "phylip",
                    "fasta", "fa", "fna", "faa", "fas",
                ],
            )
            .pick_file();

        let Some(path) = picked else {
            return;
        };

        let op_id = self.start_operation(
            "Import Pre-Aligned MAFFT",
            format!("Loading {}", path.display()),
        );

        match load_mafft_alignment_from_path(&path) {
            Ok(report) => {
                let name = report.alignment.name.clone();
                let rows = report.alignment.row_count();
                let cols = report.alignment.width();
                let flavor_label = report.flavor.label();
                let mut detail = format!(
                    "Imported {flavor_label} ({rows} sequences, {cols} aligned columns; gaps preserved verbatim)"
                );
                if report.adjusted_direction_count > 0 {
                    detail.push_str(&format!(
                        "; {} sequence(s) reverse-adjusted (_R_)",
                        report.adjusted_direction_count
                    ));
                }
                let lineage = LineageEntry {
                    operation: format!("Import Pre-Aligned MAFFT ({flavor_label})"),
                    parent_ids: Vec::new(),
                    parent_names: vec![path.display().to_string()],
                    details: detail.clone(),
                };
                self.add_workspace_doc("Alignments", false, report.alignment, None, Some(lineage));
                self.finish_operation(op_id, OperationStatus::Completed, detail);
                self.status_message = format!(
                    "Imported pre-aligned MAFFT '{name}' ({flavor_label}: {rows} seqs x {cols} cols)."
                );
            }
            Err(err) => {
                let msg = err.to_string();
                self.finish_operation(op_id, OperationStatus::Failed, msg.clone());
                self.status_message = format!("MAFFT import failed: {msg}");
            }
        }
    }

    pub fn import_gff3_dialog(&mut self) {
        let picked = rfd::FileDialog::new()
            .add_filter("GFF3 Annotations", &["gff", "gff3"])
            .pick_file();
        let Some(path) = picked else {
            return;
        };
        let Ok(text) = std::fs::read_to_string(&path) else {
            self.status_message = format!("Could not read {}", path.display());
            return;
        };
        if let Some(doc) = self.current_doc_mut() {
            match import_gff3_annotations(&mut doc.alignment, &text) {
                Ok(count) => {
                    self.status_message = format!("Imported {count} GFF3 annotations.");
                }
                Err(err) => {
                    self.status_message = format!("GFF3 import failed: {err}");
                }
            }
        }
    }

    pub fn export_alignment_dialog(&mut self) {
        let Some(doc) = self.current_doc() else {
            return;
        };
        if let Some(ref tree) = doc.tree {
            let default_name = format!(
                "{}.nwk",
                doc.alignment.name.replace([' ', '/', '\\', ':', '(', ')'], "_")
            );
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name(&default_name)
                .add_filter("Newick Phylogenetic Tree", &["nwk", "newick"])
                .save_file()
            {
                match std::fs::write(&path, tree.to_newick()) {
                    Ok(()) => self.status_message = format!("Exported Newick tree to {}", path.display()),
                    Err(e) => self.status_message = format!("Export failed: {e}"),
                }
            }
            return;
        }

        let aln = &doc.alignment;
        let default_name = format!(
            "{}.fasta",
            aln.name.replace([' ', '/', '\\', ':', '(', ')'], "_")
        );
        let picked = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("FASTA Alignment", &["fasta", "fa"])
            .add_filter("Clustal Alignment", &["aln"])
            .save_file();
        let Some(path) = picked else {
            return;
        };
        let is_clustal = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("aln"));
        let content = if is_clustal {
            write_clustal(aln)
        } else {
            write_fasta(aln)
        };
        match std::fs::write(&path, content) {
            Ok(()) => self.status_message = format!("Exported alignment to {}", path.display()),
            Err(err) => self.status_message = format!("Export failed: {err}"),
        }
    }

    pub fn load_local_structure_dialog(&mut self) {
        let picked = rfd::FileDialog::new()
            .add_filter("Protein 3D Structures", &["pdb", "cif", "mmcif"])
            .pick_file();
        let Some(path) = picked else {
            return;
        };
        let Ok(text) = std::fs::read_to_string(&path) else {
            self.viewer3d.error_message = Some(format!("Failed reading {}", path.display()));
            return;
        };
        let entry_id = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("LocalStructure")
            .to_owned();
        let is_cif = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("cif") || e.eq_ignore_ascii_case("mmcif"));

        let parsed = if is_cif {
            parse_mmcif(&text, &entry_id, None, "Local mmCIF", false)
        } else {
            parse_pdb(&text, &entry_id, None, "Local PDB", false)
        };

        match parsed {
            Ok(structure) => {
                self.viewer3d.open = true;
                self.viewer3d.structure = Some(structure);
                self.viewer3d.alphafold_meta = None;
                self.viewer3d.error_message = None;
                self.viewer3d.loading_message = None;
                self.viewer3d.camera.reset();
                self.rebuild_3d_alignment_map();
                self.status_message = format!("Loaded local 3D structure {}", path.display());
            }
            Err(err) => {
                self.viewer3d.error_message = Some(format!("3D parse failed: {err}"));
            }
        }
    }

    pub fn run_motif_search(&mut self) {
        let query = self.motif_query.clone();
        let iupac = self.motif_allow_iupac;
        if let Some(aln) = self.current_alignment() {
            self.motif_hits = aln.find_motif(&query, iupac);
            self.active_motif_hit = 0;
            if let Some(first) = self.motif_hits.first() {
                self.scroll_to_column = Some(first.start_col);
                self.selected_cols = Some((first.start_col, first.end_col));
            }
        }
    }

    fn poll_async_events(&mut self) {
        while let Ok(ev) = self.event_rx.try_recv() {
            match ev {
                AsyncEvent::UniProtSearchDone { op_id, result } => {
                    self.uniprot_searching = false;
                    match result {
                        Ok(hits) => {
                            let n = hits.len();
                            self.uniprot_search_results = hits;
                            self.uniprot_search_error = None;
                            self.finish_operation(op_id, OperationStatus::Completed, format!("Found {n} UniProt entries"));
                        }
                        Err(err) => {
                            self.uniprot_search_error = Some(err.clone());
                            self.finish_operation(op_id, OperationStatus::Failed, err);
                        }
                    }
                }
                AsyncEvent::UniProtEntryLoaded {
                    op_id,
                    entry,
                    import_annotations_for_doc,
                    attach_to_annotation,
                    save_as_new_doc,
                } => match entry {
                    Ok(up_entry) => {
                        if save_as_new_doc && !up_entry.sequence.is_empty() {
                            let seq = helix_core::Sequence::new(&up_entry.entry_name, up_entry.sequence.as_bytes())
                                .with_description(&up_entry.protein_name)
                                .with_organism(&up_entry.organism)
                                .with_uniprot(&up_entry.accession);
                            if let Ok(mut aln) = Alignment::new(
                                format!("{} ({})", up_entry.entry_name, up_entry.accession),
                                helix_core::Alphabet::Protein,
                                vec![seq],
                            ) {
                                aln.description = up_entry.protein_name.clone();
                                let feats = up_entry.build_alignment_annotations(&aln, 0);
                                for f in feats {
                                    aln.add_annotation(f);
                                }
                                let acc = up_entry.accession.clone();
                                self.add_workspace_doc("Protein Documents", false, aln, None, None);
                                self.finish_operation(
                                    op_id,
                                    OperationStatus::Completed,
                                    format!("Saved UniProt {acc} into Protein Documents"),
                                );
                            }
                        }

                        if let Some((doc_idx, ann_id)) = attach_to_annotation {
                            if let Some(doc) = self.docs.get_mut(doc_idx) {
                                if let Some(ann) = doc.alignment.annotation_by_id_mut(ann_id) {
                                    let start_res = ann
                                        .uniprot
                                        .as_ref()
                                        .map_or(1, |u| u.uniprot_start_residue);
                                    ann.uniprot = Some(
                                        UniProtLink::new(&up_entry.accession)
                                            .with_start_residue(start_res)
                                            .with_metadata(
                                                &up_entry.entry_name,
                                                &up_entry.protein_name,
                                                &up_entry.organism,
                                            ),
                                    );
                                }
                            }
                        }

                        if let Some(doc_idx) = import_annotations_for_doc {
                            if let Some(doc) = self.docs.get_mut(doc_idx) {
                                let new_anns = up_entry.build_alignment_annotations(&doc.alignment, 0);
                                let added = new_anns.len();
                                for ann in new_anns {
                                    doc.alignment.add_annotation(ann);
                                }
                                self.finish_operation(
                                    op_id,
                                    OperationStatus::Completed,
                                    format!("Imported {added} curated annotations from UniProt {}", up_entry.accession),
                                );
                                self.status_message = format!(
                                    "UniProt {} loaded — imported {added} domains/sites.",
                                    up_entry.accession
                                );
                            }
                        }
                        self.viewer3d.uniprot_entry = Some(up_entry);
                    }
                    Err(err) => {
                        if op_id > 0 {
                            self.finish_operation(op_id, OperationStatus::Failed, err.clone());
                        }
                        self.status_message = format!("UniProt fetch failed: {err}");
                    }
                },
                AsyncEvent::ProteinBundleLoaded {
                    op_id,
                    result,
                    linked_annotation_id,
                } => {
                    self.viewer3d.loading_message = None;
                    match result {
                        Ok(boxed_bundle) => {
                            let bundle = *boxed_bundle;
                            let msg = format!(
                                "Loaded {} via {} ({} residues, {})",
                                bundle.structure.entry_id,
                                bundle.structure.source_label,
                                bundle.structure.residues.len(),
                                if bundle.structure.is_alphafold {
                                    format!("mean pLDDT {:.1}", bundle.structure.mean_plddt)
                                } else {
                                    format!("{} atoms", bundle.structure.atoms.len())
                                }
                            );
                            // Automatically switch coloring to Secondary Structure for experimental X-ray PDBs,
                            // or AlphaFold pLDDT for AlphaFold models!
                            if bundle.structure.is_alphafold {
                                self.viewer3d.color_by = ColorBy3D::AlphaFoldPlddt;
                            } else if self.viewer3d.color_by == ColorBy3D::AlphaFoldPlddt {
                                self.viewer3d.color_by = ColorBy3D::SecondaryStructure;
                            }
                            self.finish_operation(op_id, OperationStatus::Completed, msg.clone());
                            self.viewer3d.structure = Some(bundle.structure);
                            self.viewer3d.alphafold_meta = bundle.alphafold_meta;
                            if let Some(up) = bundle.uniprot_entry {
                                self.viewer3d.uniprot_entry = Some(up);
                            }
                            self.viewer3d.linked_annotation_id = linked_annotation_id;
                            self.viewer3d.error_message = None;
                            self.viewer3d.camera.reset();
                            self.rebuild_3d_alignment_map();
                            self.status_message = msg;
                        }
                        Err(err) => {
                            self.finish_operation(op_id, OperationStatus::Failed, err.clone());
                            self.viewer3d.error_message = Some(err.clone());
                            self.status_message = format!("3D structure fetch failed: {err}");
                        }
                    }
                }
                AsyncEvent::RcsbPdbLoaded { op_id, result } => {
                    self.viewer3d.loading_message = None;
                    match result {
                        Ok(structure) => {
                            let msg = format!(
                                "Loaded experimental PDB {} ({} residues)",
                                structure.entry_id,
                                structure.residues.len()
                            );
                            self.finish_operation(op_id, OperationStatus::Completed, msg.clone());
                            self.viewer3d.structure = Some(structure);
                            self.viewer3d.alphafold_meta = None;
                            self.viewer3d.error_message = None;
                            self.viewer3d.camera.reset();
                            self.rebuild_3d_alignment_map();
                            self.status_message = msg;
                        }
                        Err(err) => {
                            self.finish_operation(op_id, OperationStatus::Failed, err.clone());
                            self.viewer3d.error_message = Some(err);
                        }
                    }
                }
                AsyncEvent::NcbiSearchDone { op_id, db, result } => {
                    self.ncbi_searching = false;
                    match result {
                        Ok(hits) => {
                            let n = hits.len();
                            self.ncbi_search_results = hits;
                            self.ncbi_search_error = None;
                            self.finish_operation(
                                op_id,
                                OperationStatus::Completed,
                                format!("NCBI {} returned {n} records", db.label()),
                            );
                        }
                        Err(err) => {
                            self.ncbi_search_error = Some(err.clone());
                            self.finish_operation(op_id, OperationStatus::Failed, err);
                        }
                    }
                }
                AsyncEvent::NcbiRecordDownloaded { op_id, db, result } => match result {
                    Ok(aln) => {
                        let name = aln.name.clone();
                        let folder = match db {
                            NcbiDatabase::Protein => "Protein Documents",
                            NcbiDatabase::PopSet => "Alignments",
                            _ => "Genomes",
                        };
                        self.add_workspace_doc(folder, false, aln, None, None);
                        self.finish_operation(
                            op_id,
                            OperationStatus::Completed,
                            format!("Downloaded '{name}' into {folder}"),
                        );
                        self.status_message = format!("Downloaded NCBI record '{name}' into {folder}.");
                    }
                    Err(err) => {
                        self.finish_operation(op_id, OperationStatus::Failed, err.clone());
                        self.status_message = format!("NCBI download failed: {err}");
                    }
                },
            }
        }
    }

    fn show_new_annotation_modal(&mut self, ctx: &Context) {
        if !self.new_ann_draft.open {
            return;
        }
        let mut open = self.new_ann_draft.open;
        let mut commit = false;
        let mut cancel = false;

        egui::Window::new("Add / Edit Annotation")
            .open(&mut open)
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                egui::Grid::new("new_ann_grid")
                    .num_columns(2)
                    .spacing([12.0, 8.0])
                    .show(ui, |ui| {
                        ui.label("Name:");
                        ui.text_edit_singleline(&mut self.new_ann_draft.name);
                        ui.end_row();

                        ui.label("Type:");
                        egui::ComboBox::from_id_salt("new_ann_kind")
                            .selected_text(self.new_ann_draft.kind.label())
                            .show_ui(ui, |ui| {
                                for kind in [
                                    FeatureKind::Domain,
                                    FeatureKind::Cds,
                                    FeatureKind::Gene,
                                    FeatureKind::ActiveSite,
                                    FeatureKind::BindingSite,
                                    FeatureKind::Variant,
                                    FeatureKind::Motif,
                                    FeatureKind::Primer,
                                    FeatureKind::Helix,
                                    FeatureKind::BetaStrand,
                                    FeatureKind::Misc,
                                ] {
                                    ui.selectable_value(
                                        &mut self.new_ann_draft.kind,
                                        kind,
                                        kind.label(),
                                    );
                                }
                            });
                        ui.end_row();

                        ui.label("Start (1-based):");
                        ui.add(egui::DragValue::new(&mut self.new_ann_draft.start_1based).range(1..=100_000));
                        ui.end_row();

                        ui.label("End (1-based):");
                        ui.add(egui::DragValue::new(&mut self.new_ann_draft.end_1based).range(1..=100_000));
                        ui.end_row();

                        ui.label("Strand:");
                        ui.horizontal(|ui| {
                            ui.radio_value(&mut self.new_ann_draft.strand, Strand::Forward, "Forward (+)");
                            ui.radio_value(&mut self.new_ann_draft.strand, Strand::Reverse, "Reverse (-)");
                        });
                        ui.end_row();

                        ui.label("UniProt Accession:");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.new_ann_draft.uniprot_accession)
                                .hint_text("Optional, e.g. P04637 or P68871"),
                        );
                        ui.end_row();

                        ui.label("Notes:");
                        ui.text_edit_multiline(&mut self.new_ann_draft.notes);
                        ui.end_row();
                    });

                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    if ui.button("Create Annotation").clicked() {
                        commit = true;
                    }
                    if ui.button("Cancel").clicked() {
                        cancel = true;
                    }
                });
            });

        if cancel {
            open = false;
        }
        if commit {
            let start = self.new_ann_draft.start_1based.saturating_sub(1);
            let end = self.new_ann_draft.end_1based.max(start + 1);
            let mut ann = Annotation::new(
                0,
                self.new_ann_draft.name.trim(),
                self.new_ann_draft.kind,
                start,
                end,
            )
            .with_strand(self.new_ann_draft.strand)
            .with_notes(self.new_ann_draft.notes.clone());

            let acc = self.new_ann_draft.uniprot_accession.trim();
            if !acc.is_empty() {
                ann.uniprot = Some(UniProtLink::new(acc));
            }

            if let Some(doc) = self.current_doc_mut() {
                let id = doc.alignment.add_annotation(ann);
                self.selected_annotation_id = Some(id);
                self.status_message = format!("Added annotation '{}'.", self.new_ann_draft.name);
            }
            open = false;
        }
        self.new_ann_draft.open = open;
    }
}

impl eframe::App for HelixApp {
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        self.poll_async_events();

        // Refresh process & system RAM monitor every 2 seconds (low overhead).
        if self.last_mem_sample.elapsed().as_secs() >= 2 {
            self.memory = sample_memory_mb();
            self.last_mem_sample = Instant::now();
        }

        let pal = Palette::for_mode(self.theme);

        // 1. Top Geneious Prime-style ribbon toolbar.
        egui::TopBottomPanel::top("helix_toolbar")
            .exact_height(48.0)
            .frame(
                egui::Frame::none()
                    .fill(pal.header_bg)
                    .inner_margin(egui::Margin::symmetric(10.0, 6.0))
                    .stroke(egui::Stroke::new(1.0, pal.border)),
            )
            .show(ctx, |ui| {
                widgets::toolbar::show_toolbar(self, ctx, ui);
            });

        // 2. Bottom status bar with RAM Memory Monitor (left) & Cursor Base Position (middle/right).
        egui::TopBottomPanel::bottom("helix_status_bar")
            .exact_height(26.0)
            .frame(
                egui::Frame::none()
                    .fill(pal.header_bg)
                    .inner_margin(egui::Margin::symmetric(10.0, 4.0))
                    .stroke(egui::Stroke::new(1.0, pal.border)),
            )
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    // Left RAM monitor matching Geneious Prime (`574 / 14,380 MB Memory`).
                    if self.memory.total_mb > 0 {
                        ui.label(
                            egui::RichText::new(format!(
                                "{} / {} MB Memory",
                                self.memory.used_mb, self.memory.total_mb
                            ))
                            .small()
                            .monospace()
                            .color(pal.text_muted),
                        );
                        ui.separator();
                    }

                    // Cursor base coordinate indicator (`Cursor before base 1,179.`).
                    let pos_col = self.hovered_column.unwrap_or(self.cursor_col);
                    ui.label(
                        egui::RichText::new(format!("Cursor before base {}.", pos_col + 1))
                            .small()
                            .color(pal.text_primary),
                    );
                    ui.separator();

                    ui.label(
                        egui::RichText::new(&self.status_message)
                            .small()
                            .color(pal.text_muted),
                    );

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let running_ops = self
                            .operations
                            .iter()
                            .filter(|o| o.status == OperationStatus::Running)
                            .count();
                        if running_ops > 0 {
                            ui.label(
                                egui::RichText::new(format!("{running_ops} running ops"))
                                    .small()
                                    .strong()
                                    .color(pal.accent),
                            );
                            ui.spinner();
                        } else if let Some(aln) = self.current_alignment() {
                            if let Some(col) = self.hovered_column {
                                if let Some(stats) = aln.columns().get(col) {
                                    ui.label(
                                        egui::RichText::new(format!(
                                            "Col {} | Cons: {} | ID: {:.0}% | H: {:.2}b",
                                            col + 1,
                                            char::from(stats.consensus),
                                            stats.identity_pct,
                                            stats.shannon_entropy
                                        ))
                                        .small()
                                        .monospace()
                                        .color(pal.text_primary),
                                    );
                                }
                            }
                        }
                    });
                });
            });

        // 3. Left Geneious Prime Sources Tree (Local, Reference Features, Deleted, Shared, Operations, NCBI, UniProt).
        if self.show_sources_panel {
            egui::SidePanel::left("helix_sources_panel")
                .default_width(235.0)
                .min_width(190.0)
                .max_width(340.0)
                .frame(
                    egui::Frame::none()
                        .fill(pal.panel_bg)
                        .inner_margin(egui::Margin::same(8.0))
                        .stroke(egui::Stroke::new(1.0, pal.border)),
                )
                .show(ctx, |ui| {
                    widgets::sources_panel::show_sources_panel(self, ctx, ui);
                });
        }

        // 4. Right Geneious Prime General & Controls Drawer.
        if self.show_inspector_panel {
            egui::SidePanel::right("helix_inspector_panel")
                .default_width(290.0)
                .min_width(240.0)
                .max_width(420.0)
                .frame(
                    egui::Frame::none()
                        .fill(pal.panel_bg)
                        .inner_margin(egui::Margin::same(8.0))
                        .stroke(egui::Stroke::new(1.0, pal.border)),
                )
                .show(ctx, |ui| {
                    widgets::inspector_panel::show_inspector_panel(self, ctx, ui);
                });
        }

        // 5. Center Workspace: Upper Document Table + Lower Tabbed Document Viewer (exact Geneious layout!).
        egui::CentralPanel::default()
            .frame(egui::Frame::none().fill(pal.canvas_bg).inner_margin(egui::Margin::same(0.0)))
            .show(ctx, |ui| {
                if self.show_doc_table {
                    let table_h = match self.source_selection {
                        SourceSelection::Operations
                        | SourceSelection::ReferenceFeatures
                        | SourceSelection::Ncbi(_)
                        | SourceSelection::UniProt => (ui.available_height() * 0.44).clamp(180.0, 340.0),
                        _ => (ui.available_height() * 0.25).clamp(125.0, 220.0),
                    };
                    egui::TopBottomPanel::top("helix_doc_table_pane")
                        .resizable(true)
                        .default_height(table_h)
                        .min_height(90.0)
                        .max_height(460.0)
                        .frame(
                            egui::Frame::none()
                                .fill(pal.canvas_bg)
                                .stroke(egui::Stroke::new(1.0, pal.border)),
                        )
                        .show_inside(ui, |ui| {
                            widgets::document_table::show_document_table(self, ctx, ui);
                        });
                }

                widgets::doc_viewer::show_doc_viewer(self, ctx, ui);
            });

        // 6. Modals for Annotation Creation & Molecular Biology Tools.
        self.show_new_annotation_modal(ctx);
        widgets::tools_modals::show_all_tool_modals(self, ctx);

        // 7. Separate native OS Window for the 3D Protein Viewer!
        if self.viewer3d.open {
            let window_title = match &self.viewer3d.structure {
                Some(st) => format!("HelixPrime 3D Protein Viewer — {}", st.entry_id),
                None => String::from("HelixPrime 3D Protein Viewer — AlphaFold / UniProt"),
            };

            let viewport_id = ViewportId::from_hash_of("helix_prime_3d_protein_viewport");
            let builder = ViewportBuilder::default()
                .with_title(window_title)
                .with_icon(crate::app_window_icon())
                .with_inner_size([920.0, 680.0])
                .with_min_inner_size([560.0, 420.0]);

            ctx.show_viewport_immediate(viewport_id, builder, |vp_ctx, _class| {
                apply_theme(vp_ctx, self.theme);
                widgets::protein_window::show_protein_window(self, vp_ctx);
                if vp_ctx.input(|i| i.viewport().close_requested()) {
                    self.viewer3d.open = false;
                }
            });
        }

        // Flush state immediately when the main OS window close is requested.
        if ctx.input(|i| i.viewport().close_requested()) {
            let _ = self.persist_session(true);
        }
    }

    fn save(&mut self, _storage: &mut dyn eframe::Storage) {
        // Called periodically by eframe (~30s) and during shutdown.
        let _ = self.persist_session(false);
    }

    fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {
        // Guarantee a final session flush on close even if eframe's built-in storage is unavailable.
        let _ = self.persist_session(false);
    }
}
