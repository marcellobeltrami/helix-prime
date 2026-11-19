//! Top-level application state, background UniProt/AlphaFold worker orchestration,
//! and main `eframe::App` update loop including the separate OS 3D viewport window.

use crate::renderer3d::{Camera3D, ColorBy3D, RenderStyle, Renderer3D};
use crate::theme::{apply_theme, ColorScheme, HighlightMode, Palette, ThemeMode};
use crate::widgets;
use egui::{Context, ViewportBuilder, ViewportId};
use helix_core::database::{
    fetch_alphafold_structure, fetch_rcsb_pdb, fetch_uniprot_entry, search_uniprot, AlphaFoldMeta,
    UniProtEntry, UniProtSearchHit,
};
use helix_core::demo::{build_default_workspace, build_offline_preview_structure};
use helix_core::io::{
    clustal::write_clustal, fasta::write_fasta, gff3::import_gff3_annotations,
    load_alignment_from_path,
};
use helix_core::{
    parse_mmcif, parse_pdb, Alignment, Annotation, FeatureKind, MotifHit, ProteinStructure, Strand,
    StructureAlignmentMap, UniProtLink,
};
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver, Sender};

/// A single alignment document inside the Geneious-style Sources workspace.
pub struct WorkspaceDoc {
    pub folder: String,
    pub alignment: Alignment,
    /// Cached N x N pairwise identity matrix (recomputed lazily when alignment changes).
    pub pairwise_cache: Option<Vec<f32>>,
}

impl WorkspaceDoc {
    pub fn new(folder: impl Into<String>, alignment: Alignment) -> Self {
        Self {
            folder: folder.into(),
            alignment,
            pairwise_cache: None,
        }
    }

    pub fn invalidate_caches(&mut self) {
        self.pairwise_cache = None;
    }
}

/// Right-side inspector drawer tabs inspired by Geneious Prime's right control sidebar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InspectorTab {
    #[default]
    Display,
    Annotations,
    FindMotif,
    Statistics,
}

/// Primary alignment visualization layout mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ViewMode {
    /// Traditional Geneious Prime horizontal linear alignment grid.
    #[default]
    Linear,
    /// Circular phage / plasmid genome ring visualization.
    CircularPhage,
}

/// Interactive camera and track settings for [`ViewMode::CircularPhage`].
#[derive(Debug, Clone)]
pub struct CircularViewState {
    /// Angular rotation offset in radians (`0.0` places origin column 1 at 12 o'clock).
    pub rotation_rad: f32,
    /// Radial zoom multiplier (`0.55..=2.4`).
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

/// Visual configuration for the main alignment viewport.
#[derive(Debug, Clone)]
pub struct DisplaySettings {
    pub color_scheme: ColorScheme,
    pub highlight_mode: HighlightMode,
    pub dots_for_matches: bool,
    pub show_overview: bool,
    pub show_identity_graph: bool,
    pub show_consensus: bool,
    pub show_annotations: bool,
    pub show_translation: bool,
    pub translation_frame: usize,
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
            show_translation: true,
            translation_frame: 0,
            col_width: 14.0,
            row_height: 22.0,
        }
    }
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

/// State backing the separate native OS 3D Protein Viewer window.
pub struct ProteinViewerState {
    pub open: bool,
    pub accession_input: String,
    pub structure: Option<ProteinStructure>,
    pub alphafold_meta: Option<AlphaFoldMeta>,
    pub uniprot_entry: Option<UniProtEntry>,
    pub alignment_map: Option<StructureAlignmentMap>,
    /// Optional ID of the alignment annotation that opened this structure.
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

/// Asynchronous messages delivered from background network threads back to the GUI thread.
pub enum AsyncEvent {
    UniProtSearchDone(Result<Vec<UniProtSearchHit>, String>),
    UniProtEntryLoaded {
        entry: Result<UniProtEntry, String>,
        import_annotations_for_doc: Option<usize>,
        attach_to_annotation: Option<(usize, u64)>,
    },
    AlphaFoldLoaded {
        result: Result<(ProteinStructure, AlphaFoldMeta), String>,
        linked_annotation_id: Option<u64>,
    },
    RcsbPdbLoaded {
        result: Result<ProteinStructure, String>,
    },
}

/// Main HelixPrime desktop application state.
pub struct HelixApp {
    pub theme: ThemeMode,
    pub view_mode: ViewMode,
    pub circular: CircularViewState,
    pub docs: Vec<WorkspaceDoc>,
    pub active_doc: usize,
    pub show_sources_panel: bool,
    pub show_inspector_panel: bool,
    pub inspector_tab: InspectorTab,
    pub display: DisplaySettings,

    // Alignment interaction & cross-window synchronization state.
    pub horizontal_scroll_px: f32,
    pub scroll_to_column: Option<usize>,
    pub selected_cols: Option<(usize, usize)>,
    pub selected_annotation_id: Option<u64>,
    pub hovered_column: Option<usize>,
    pub hovered_row: Option<usize>,

    // Motif search state.
    pub motif_query: String,
    pub motif_allow_iupac: bool,
    pub motif_hits: Vec<MotifHit>,
    pub active_motif_hit: usize,

    // UniProt search state inside the Inspector.
    pub uniprot_search_query: String,
    pub uniprot_search_results: Vec<UniProtSearchHit>,
    pub uniprot_searching: bool,
    pub uniprot_search_error: Option<String>,

    // Annotation creation modal.
    pub new_ann_draft: NewAnnotationDraft,

    // Separate OS window 3D protein viewer.
    pub viewer3d: ProteinViewerState,

    // Background worker channel.
    pub event_tx: Sender<AsyncEvent>,
    pub event_rx: Receiver<AsyncEvent>,

    // Status bar feedback message.
    pub status_message: String,
}

impl HelixApp {
    /// Initializes HelixPrime with curated biological demo alignments and ready-to-use 3D preview.
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        apply_theme(&cc.egui_ctx, ThemeMode::Light);

        let (event_tx, event_rx) = mpsc::channel();
        let docs: Vec<WorkspaceDoc> = build_default_workspace()
            .into_iter()
            .map(|d| WorkspaceDoc::new(d.folder, d.alignment))
            .collect();

        let mut app = Self {
            theme: ThemeMode::Light,
            view_mode: ViewMode::Linear,
            circular: CircularViewState::default(),
            docs,
            active_doc: 0,
            show_sources_panel: true,
            show_inspector_panel: true,
            inspector_tab: InspectorTab::Annotations,
            display: DisplaySettings::default(),
            horizontal_scroll_px: 0.0,
            scroll_to_column: None,
            selected_cols: None,
            selected_annotation_id: None,
            hovered_column: None,
            hovered_row: None,
            motif_query: String::new(),
            motif_allow_iupac: true,
            motif_hits: Vec::new(),
            active_motif_hit: 0,
            uniprot_search_query: String::from("TP53 human"),
            uniprot_search_results: Vec::new(),
            uniprot_searching: false,
            uniprot_search_error: None,
            new_ann_draft: NewAnnotationDraft::default(),
            viewer3d: ProteinViewerState::default(),
            event_tx,
            event_rx,
            status_message: String::from(
                "Ready — Click any UniProt-linked annotation or '3D Protein Window' to inspect AlphaFold 3D structures.",
            ),
        };

        // Pre-select the first annotation so the Inspector immediately displays rich UniProt controls.
        if let Some(first_ann_id) = app
            .current_alignment()
            .and_then(|a| a.annotations().first().map(|ann| ann.id))
        {
            app.selected_annotation_id = Some(first_ann_id);
        }

        app
    }

    /// Returns a shared reference to the currently active alignment, if any.
    #[must_use]
    pub fn current_alignment(&self) -> Option<&Alignment> {
        self.docs.get(self.active_doc).map(|d| &d.alignment)
    }

    /// Returns a mutable reference to the currently active workspace document.
    pub fn current_doc_mut(&mut self) -> Option<&mut WorkspaceDoc> {
        self.docs.get_mut(self.active_doc)
    }

    /// Updates the active Light/Dark theme and applies styling immediately.
    pub fn set_theme(&mut self, ctx: &Context, mode: ThemeMode) {
        self.theme = mode;
        apply_theme(ctx, mode);
    }

    /// Switches the active alignment document and refreshes dependent caches/mappings.
    pub fn select_document(&mut self, idx: usize) {
        if idx >= self.docs.len() {
            return;
        }
        self.active_doc = idx;
        if self.docs[idx].folder.eq_ignore_ascii_case("Phage Genomes") {
            self.view_mode = ViewMode::CircularPhage;
        }
        self.horizontal_scroll_px = 0.0;
        self.selected_cols = None;
        self.hovered_column = None;
        self.motif_hits.clear();
        self.active_motif_hit = 0;
        self.selected_annotation_id = self
            .current_alignment()
            .and_then(|a| a.annotations().first().map(|ann| ann.id));

        if let Some(aln) = self.current_alignment() {
            self.status_message = format!(
                "Opened '{}' ({} sequences, {} {}, {:.1}% pairwise identity)",
                aln.name,
                aln.row_count(),
                aln.width(),
                if aln.alphabet.is_nucleotide() {
                    "bp"
                } else {
                    "aa"
                },
                aln.mean_pairwise_identity_pct()
            );
        }
        self.rebuild_3d_alignment_map();
    }

    /// Recomputes the bidirectional coordinate map between the loaded 3D structure and the
    /// active alignment/annotation.
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

    /// Triggers an asynchronous AlphaFold DB 3D structure fetch + UniProt metadata fetch
    /// and opens the separate 3D Protein Viewer OS window.
    pub fn open_3d_for_accession(
        &mut self,
        ctx: &Context,
        accession: &str,
        linked_annotation_id: Option<u64>,
    ) {
        let clean = accession.trim().to_ascii_uppercase();
        if clean.is_empty() {
            self.status_message = String::from("Enter a valid UniProt accession (e.g. P68871 or P04637).");
            return;
        }

        self.viewer3d.open = true;
        self.viewer3d.accession_input = clean.clone();
        self.viewer3d.linked_annotation_id = linked_annotation_id;
        self.viewer3d.error_message = None;
        self.viewer3d.loading_message = Some(format!(
            "Fetching AlphaFold DB 3D model & UniProt record for {clean}..."
        ));
        self.status_message = format!("Querying AlphaFold DB & UniProtKB for {clean}...");

        // Spawn background fetch for AlphaFold 3D structure.
        let tx_af = self.event_tx.clone();
        let ctx_af = ctx.clone();
        let acc_af = clean.clone();
        std::thread::spawn(move || {
            let result = fetch_alphafold_structure(&acc_af).map_err(|e| e.to_string());
            let _ = tx_af.send(AsyncEvent::AlphaFoldLoaded {
                result,
                linked_annotation_id,
            });
            ctx_af.request_repaint();
        });

        // Also fetch UniProt metadata in parallel so functional summary & PDB xrefs populate.
        let tx_up = self.event_tx.clone();
        let ctx_up = ctx.clone();
        let acc_up = clean;
        std::thread::spawn(move || {
            let entry = fetch_uniprot_entry(&acc_up).map_err(|e| e.to_string());
            let _ = tx_up.send(AsyncEvent::UniProtEntryLoaded {
                entry,
                import_annotations_for_doc: None,
                attach_to_annotation: None,
            });
            ctx_up.request_repaint();
        });
    }

    /// Loads the built-in offline schematic 3D protein preview immediately in the separate window.
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
        self.status_message = String::from(
            "Loaded offline HBB 3D schematic preview in separate window. Click 'Fetch AlphaFold DB' for live full-atom model.",
        );
    }

    /// Triggers an asynchronous UniProtKB text search.
    pub fn trigger_uniprot_search(&mut self, ctx: &Context) {
        let q = self.uniprot_search_query.trim().to_owned();
        if q.is_empty() {
            return;
        }
        self.uniprot_searching = true;
        self.uniprot_search_error = None;
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let res = search_uniprot(&q, 10).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::UniProtSearchDone(res));
            ctx_clone.request_repaint();
        });
    }

    /// Fetches a UniProt entry and imports its curated domains/active sites onto the active alignment.
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
        self.status_message = format!("Fetching UniProtKB features for {clean}...");
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let entry = fetch_uniprot_entry(&clean).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::UniProtEntryLoaded {
                entry,
                import_annotations_for_doc: Some(doc_idx),
                attach_to_annotation: attach_to_ann_id.map(|id| (doc_idx, id)),
            });
            ctx_clone.request_repaint();
        });
    }

    /// Fetches an experimental X-ray / cryo-EM structure from RCSB PDB by 4-letter PDB ID.
    pub fn fetch_rcsb_structure(&mut self, ctx: &Context, pdb_id: &str) {
        let clean = pdb_id.trim().to_ascii_uppercase();
        let acc = self
            .viewer3d
            .uniprot_entry
            .as_ref()
            .map(|e| e.accession.clone());
        self.viewer3d.open = true;
        self.viewer3d.error_message = None;
        self.viewer3d.loading_message = Some(format!("Downloading experimental RCSB PDB {clean}..."));
        let tx = self.event_tx.clone();
        let ctx_clone = ctx.clone();
        std::thread::spawn(move || {
            let result = fetch_rcsb_pdb(&clean, acc).map_err(|e| e.to_string());
            let _ = tx.send(AsyncEvent::RcsbPdbLoaded { result });
            ctx_clone.request_repaint();
        });
    }

    /// Opens a native file dialog to import a FASTA, Clustal, or GenBank alignment.
    pub fn import_alignment_dialog(&mut self) {
        let picked: Option<PathBuf> = rfd::FileDialog::new()
            .add_filter(
                "Bioinformatics Alignments",
                &["fasta", "fa", "fna", "faa", "aln", "clustal", "gb", "gbk"],
            )
            .pick_file();

        let Some(path) = picked else {
            return;
        };

        match load_alignment_from_path(&path) {
            Ok(aln) => {
                let name = aln.name.clone();
                self.docs.push(WorkspaceDoc::new("Imported Files", aln));
                let new_idx = self.docs.len() - 1;
                self.select_document(new_idx);
                self.status_message = format!("Imported '{name}' from {}", path.display());
            }
            Err(err) => {
                self.status_message = format!("Import failed: {err}");
            }
        }
    }

    /// Opens a native file dialog to overlay GFF3 annotations onto the active alignment.
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

    /// Opens a native file dialog to export the active alignment in FASTA or Clustal format.
    pub fn export_alignment_dialog(&mut self) {
        let Some(aln) = self.current_alignment() else {
            return;
        };
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
            Ok(()) => {
                self.status_message = format!("Exported alignment to {}", path.display());
            }
            Err(err) => {
                self.status_message = format!("Export failed: {err}");
            }
        }
    }

    /// Opens a native file dialog to load a local `.pdb` or `.cif` 3D structure into the 3D window.
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

    /// Re-runs the motif search over the active alignment.
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

    /// Drains all completed background worker events without blocking the UI thread.
    fn poll_async_events(&mut self) {
        while let Ok(ev) = self.event_rx.try_recv() {
            match ev {
                AsyncEvent::UniProtSearchDone(res) => {
                    self.uniprot_searching = false;
                    match res {
                        Ok(hits) => {
                            self.uniprot_search_results = hits;
                            self.uniprot_search_error = None;
                        }
                        Err(err) => {
                            self.uniprot_search_error = Some(err);
                        }
                    }
                }
                AsyncEvent::UniProtEntryLoaded {
                    entry,
                    import_annotations_for_doc,
                    attach_to_annotation,
                } => match entry {
                    Ok(up_entry) => {
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
                                    if ann.notes.is_empty() && !up_entry.function_summary.is_empty()
                                    {
                                        ann.notes = up_entry.function_summary.clone();
                                    }
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
                                self.status_message = format!(
                                    "UniProt {} ({}) loaded — imported {added} curated domain/site annotations.",
                                    up_entry.accession, up_entry.entry_name
                                );
                            }
                        }
                        self.viewer3d.uniprot_entry = Some(up_entry);
                    }
                    Err(err) => {
                        self.status_message = format!("UniProt fetch failed: {err}");
                    }
                },
                AsyncEvent::AlphaFoldLoaded {
                    result,
                    linked_annotation_id,
                } => {
                    self.viewer3d.loading_message = None;
                    match result {
                        Ok((structure, meta)) => {
                            let msg = format!(
                                "Loaded AlphaFold 3D model {} ({} residues, mean pLDDT {:.1})",
                                meta.entry_id,
                                structure.residues.len(),
                                structure.mean_plddt
                            );
                            self.viewer3d.structure = Some(structure);
                            self.viewer3d.alphafold_meta = Some(meta);
                            self.viewer3d.linked_annotation_id = linked_annotation_id;
                            self.viewer3d.error_message = None;
                            self.viewer3d.camera.reset();
                            self.rebuild_3d_alignment_map();
                            self.status_message = msg;
                        }
                        Err(err) => {
                            self.viewer3d.error_message = Some(err.clone());
                            self.status_message = format!("AlphaFold fetch failed: {err}");
                        }
                    }
                }
                AsyncEvent::RcsbPdbLoaded { result } => {
                    self.viewer3d.loading_message = None;
                    match result {
                        Ok(structure) => {
                            self.status_message = format!(
                                "Loaded experimental PDB {} ({} residues, {} atoms)",
                                structure.entry_id,
                                structure.residues.len(),
                                structure.atoms.len()
                            );
                            self.viewer3d.structure = Some(structure);
                            self.viewer3d.alphafold_meta = None;
                            self.viewer3d.error_message = None;
                            self.viewer3d.camera.reset();
                            self.rebuild_3d_alignment_map();
                        }
                        Err(err) => {
                            self.viewer3d.error_message = Some(err);
                        }
                    }
                }
            }
        }
    }

    /// Renders the "Add Annotation" modal window if open.
    fn show_new_annotation_modal(&mut self, ctx: &Context) {
        if !self.new_ann_draft.open {
            return;
        }
        let mut open = self.new_ann_draft.open;
        let mut commit = false;
        let mut cancel = false;

        egui::Window::new("Add Alignment Annotation")
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

        let pal = Palette::for_mode(self.theme);

        // 1. Top Geneious Prime-style modern toolbar.
        egui::TopBottomPanel::top("helix_toolbar")
            .exact_height(46.0)
            .frame(
                egui::Frame::none()
                    .fill(pal.header_bg)
                    .inner_margin(egui::Margin::symmetric(10.0, 6.0))
                    .stroke(egui::Stroke::new(1.0, pal.border)),
            )
            .show(ctx, |ui| {
                widgets::toolbar::show_toolbar(self, ctx, ui);
            });

        // 2. Bottom status bar showing coordinates, conservation, and network status.
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
                    ui.label(
                        egui::RichText::new(&self.status_message)
                            .small()
                            .color(pal.text_muted),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if let Some(aln) = self.current_alignment() {
                            if let Some(col) = self.hovered_column {
                                if let Some(stats) = aln.columns().get(col) {
                                    ui.label(
                                        egui::RichText::new(format!(
                                            "Col {}  |  Consensus: {}  |  Identity: {:.0}%  |  Entropy: {:.2} bits",
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
                            } else if let Some((s, e)) = self.selected_cols {
                                ui.label(
                                    egui::RichText::new(format!(
                                        "Selected: {}..{} ({} cols)",
                                        s + 1,
                                        e,
                                        e.saturating_sub(s)
                                    ))
                                    .small()
                                    .monospace()
                                    .color(pal.accent),
                                );
                            }
                        }
                    });
                });
            });

        // 3. Left Geneious Prime Sources & Documents panel.
        if self.show_sources_panel {
            egui::SidePanel::left("helix_sources_panel")
                .default_width(240.0)
                .min_width(180.0)
                .max_width(360.0)
                .frame(
                    egui::Frame::none()
                        .fill(pal.panel_bg)
                        .inner_margin(egui::Margin::same(8.0))
                        .stroke(egui::Stroke::new(1.0, pal.border)),
                )
                .show(ctx, |ui| {
                    widgets::sources_panel::show_sources_panel(self, ui);
                });
        }

        // 4. Right Geneious Prime Inspector & Controls drawer.
        if self.show_inspector_panel {
            egui::SidePanel::right("helix_inspector_panel")
                .default_width(305.0)
                .min_width(250.0)
                .max_width(440.0)
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

        // 5. Center Alignment Inspection Canvas (Linear MSA or Circular Phage Map).
        egui::CentralPanel::default()
            .frame(egui::Frame::none().fill(pal.canvas_bg).inner_margin(egui::Margin::same(0.0)))
            .show(ctx, |ui| match self.view_mode {
                ViewMode::Linear => widgets::alignment_view::show_alignment_view(self, ctx, ui),
                ViewMode::CircularPhage => {
                    widgets::circular_view::show_circular_view(self, ctx, ui)
                }
            });

        // 6. Optional Add-Annotation modal.
        self.show_new_annotation_modal(ctx);

        // 7. Separate native OS Window for the 3D Protein Viewer!
        if self.viewer3d.open {
            let window_title = match &self.viewer3d.structure {
                Some(st) => format!("HelixPrime 3D Protein Viewer — {}", st.entry_id),
                None => String::from("HelixPrime 3D Protein Viewer — AlphaFold / UniProt"),
            };

            let viewport_id = ViewportId::from_hash_of("helix_prime_3d_protein_viewport");
            let builder = ViewportBuilder::default()
                .with_title(window_title)
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
    }
}
