//! Right-click context menus providing fast mouse access to commonly used
//! HelixPrime bioinformatics tools across the Sequence View, Circular Phage/Plasmid Map,
//! Document Table, and 3D Protein Viewer.

use crate::app::{DocTab, HelixApp, InspectorTab, LineageEntry};
use crate::renderer3d::{ColorBy3D, RenderStyle};
use crate::theme::ThemeMode;
use egui::{Context, RichText, Ui};
use helix_core::database::ProteinStructureSource;
use helix_core::features_db::{auto_annotate_from_reference_db, PLASMID_REFERENCE_FEATURES};

/// Metadata about an annotation under the mouse pointer in the Circular View.
#[derive(Debug, Clone)]
pub struct HoveredAnnotationContext {
    pub id: u64,
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub uniprot_acc: Option<String>,
}

/// Renders the right-click context menu for the Linear and Wrapped Sequence/Alignment views.
pub fn show_sequence_context_menu(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    ui.set_min_width(235.0);

    let (is_nuc, doc_width, active_acc) = match app.current_alignment() {
        Some(aln) => (
            aln.alphabet.is_nucleotide(),
            aln.width(),
            resolve_active_uniprot(app),
        ),
        None => return,
    };

    let header = match app.selected_cols {
        Some((s, e)) if e > s + 1 => format!("Selection {}..{} ({} cols)", s + 1, e, e - s),
        _ => format!("Cursor at base {}", app.cursor_col + 1),
    };
    ui.label(RichText::new(header).strong().small());
    ui.separator();

    // --- Sequence & Selection Actions ---
    if ui.button("BLAST Selection / Sequence...").clicked() {
        open_blast_with_selection(app);
        ui.close_menu();
    }

    if ui.button("Copy Selection / Sequence").clicked() {
        if let Some(seq_text) = extract_selection_text(app, 120) {
            ctx.copy_text(seq_text);
            app.status_message = String::from("Copied sequence text to clipboard.");
        }
        ui.close_menu();
    }

    if ui.button("Extract Region to New Document").clicked() {
        run_extract_region(app);
        ui.close_menu();
    }

    if ui
        .add_enabled(is_nuc, egui::Button::new("Reverse Complement (R.C.)"))
        .clicked()
    {
        run_reverse_complement(app);
        ui.close_menu();
    }

    if ui
        .add_enabled(is_nuc, egui::Button::new("Translate DNA/RNA to Protein"))
        .clicked()
    {
        run_translate(app);
        ui.close_menu();
    }

    ui.separator();

    // --- Annotation & Motif Actions ---
    if ui.button("Add Annotation on Selection...").clicked() {
        let (s, e) = app.selected_cols.unwrap_or_else(|| {
            let start = app.cursor_col.min(doc_width.saturating_sub(1));
            let end = (start + 30).min(doc_width.max(1));
            (start, end)
        });
        app.new_ann_draft.start_1based = s + 1;
        app.new_ann_draft.end_1based = e.max(s + 1);
        app.new_ann_draft.open = true;
        ui.close_menu();
    }

    if ui
        .add_enabled(is_nuc, egui::Button::new("Auto-Annotate Reference Features"))
        .clicked()
    {
        run_auto_annotate(app);
        ui.close_menu();
    }

    if ui.button("Find Motif in Sequence...").clicked() {
        if let Some((s, e)) = app.selected_cols {
            if e > s && e - s <= 40 {
                if let Some(txt) = extract_selection_text(app, 40) {
                    let ungapped: String = txt.chars().filter(|&c| c != '-' && c != '.').collect();
                    if !ungapped.is_empty() {
                        app.motif_query = ungapped;
                        app.run_motif_search();
                    }
                }
            }
        }
        app.show_inspector_panel = true;
        app.inspector_tab = InspectorTab::FindMotif;
        ui.close_menu();
    }

    ui.separator();

    // --- Core Molecular Biology & Alignment Tools ---
    if ui.button("Align / Assemble (MAFFT, MUSCLE, Clustal, minimap2)...").clicked() {
        app.assemble_modal.open = true;
        ui.close_menu();
    }

    if ui.button("Import Pre-Aligned MAFFT File...").clicked() {
        app.import_mafft_dialog();
        ui.close_menu();
    }

    if ui.button("Build Phylogenetic Tree...").clicked() {
        app.tree_modal.open = true;
        ui.close_menu();
    }

    if ui.button("Design PCR Primers...").clicked() {
        app.primers_modal.open = true;
        ui.close_menu();
    }

    if ui.button("In-Silico Cloning...").clicked() {
        app.cloning_modal.open = true;
        ui.close_menu();
    }

    if ui.button("Run Automated Workflow...").clicked() {
        app.workflow_modal.open = true;
        ui.close_menu();
    }

    ui.separator();

    // --- 3D Protein Viewer & View Overlays ---
    let label_3d = match &active_acc {
        Some(acc) => format!("Open 3D Protein Viewer ({acc})"),
        None => String::from("Open 3D Protein Viewer (Preview)"),
    };
    if ui.button(label_3d).clicked() {
        if let Some(acc) = active_acc {
            app.open_3d_for_accession(ctx, &acc, app.selected_annotation_id);
        } else {
            app.open_offline_3d_preview();
        }
        ui.close_menu();
    }

    ui.menu_button("View & Overlays", |ui| {
        ui.checkbox(&mut app.display.show_circular_overview, "Circular Phage / Plasmid Overview");
        ui.checkbox(&mut app.display.show_linear_view, "Linear Sequence View");
        ui.checkbox(&mut app.display.wrap_sequence, "Wrap Sequence");
        ui.separator();
        ui.checkbox(&mut app.display.allow_editing, "Allow Direct Base Editing");
        ui.checkbox(&mut app.display.show_annotations, "Annotations Track");
        ui.checkbox(&mut app.display.show_identity_graph, "Identity Graph");
        ui.checkbox(&mut app.display.show_consensus, "Consensus Track");
        if is_nuc {
            ui.checkbox(&mut app.display.show_complement, "Complement Strand (3'->5')");
            ui.checkbox(&mut app.display.show_translation, "Translation Track");
            ui.checkbox(&mut app.display.show_restriction_sites, "Restriction Sites");
        }
        ui.separator();
        let next_theme = if app.theme == ThemeMode::Dark {
            ThemeMode::Light
        } else {
            ThemeMode::Dark
        };
        let theme_label = if next_theme == ThemeMode::Dark {
            "Switch to Dark Mode"
        } else {
            "Switch to Light Mode"
        };
        if ui.button(theme_label).clicked() {
            app.set_theme(ctx, next_theme);
            ui.close_menu();
        }
    });
}

/// Renders the right-click context menu for the Circular Phage / Plasmid Genome Map.
pub fn show_circular_context_menu(
    app: &mut HelixApp,
    ctx: &Context,
    ui: &mut Ui,
    hovered_col: Option<usize>,
    hovered_ann: Option<&HoveredAnnotationContext>,
) {
    ui.set_min_width(235.0);

    let is_nuc = app
        .current_alignment()
        .is_some_and(|a| a.alphabet.is_nucleotide());

    if let Some(ann) = hovered_ann {
        ui.label(
            RichText::new(format!("{} ({}..{})", ann.name, ann.start + 1, ann.end))
                .strong()
                .small(),
        );
        ui.separator();

        if let Some(ref acc) = ann.uniprot_acc {
            if ui.button(format!("Open 3D Protein Viewer ({acc})")).clicked() {
                app.selected_annotation_id = Some(ann.id);
                app.open_3d_for_accession(ctx, acc, Some(ann.id));
                ui.close_menu();
            }
        }
        if ui.button("Zoom Linear View to Feature").clicked() {
            app.selected_annotation_id = Some(ann.id);
            app.selected_cols = Some((ann.start, ann.end));
            app.cursor_col = ann.start;
            app.scroll_to_column = Some(ann.start);
            app.display.show_linear_view = true;
            ui.close_menu();
        }
        if ui.button("Extract Feature to New Document").clicked() {
            app.selected_cols = Some((ann.start, ann.end));
            run_extract_region(app);
            ui.close_menu();
        }
        if ui.button("BLAST This Feature...").clicked() {
            app.selected_cols = Some((ann.start, ann.end));
            open_blast_with_selection(app);
            ui.close_menu();
        }
        ui.separator();
    } else {
        let base_num = hovered_col.unwrap_or(app.cursor_col) + 1;
        ui.label(RichText::new(format!("Circular Map — Base {base_num}")).strong().small());
        ui.separator();

        if ui.button("Zoom Linear View to Base").clicked() {
            let col = hovered_col.unwrap_or(app.cursor_col);
            app.cursor_col = col;
            app.scroll_to_column = Some(col);
            app.display.show_linear_view = true;
            ui.close_menu();
        }
    }

    // Core tools
    if ui.button("Add Annotation Here...").clicked() {
        let col = hovered_col.unwrap_or(app.cursor_col);
        let width = app.current_alignment().map_or(1, |a| a.width().max(1));
        app.new_ann_draft.start_1based = col + 1;
        app.new_ann_draft.end_1based = (col + 60).min(width);
        app.new_ann_draft.open = true;
        ui.close_menu();
    }
    if ui
        .add_enabled(is_nuc, egui::Button::new("Auto-Annotate Reference Features"))
        .clicked()
    {
        run_auto_annotate(app);
        ui.close_menu();
    }
    if ui.button("BLAST Genome / Selection...").clicked() {
        open_blast_with_selection(app);
        ui.close_menu();
    }
    if ui.button("Align / Assemble...").clicked() {
        app.assemble_modal.open = true;
        ui.close_menu();
    }
    if ui.button("Import Pre-Aligned MAFFT File...").clicked() {
        app.import_mafft_dialog();
        ui.close_menu();
    }
    if ui.button("Design PCR Primers...").clicked() {
        app.primers_modal.open = true;
        ui.close_menu();
    }
    if ui.button("In-Silico Cloning...").clicked() {
        app.cloning_modal.open = true;
        ui.close_menu();
    }
    if ui
        .add_enabled(is_nuc, egui::Button::new("Reverse Complement (R.C.)"))
        .clicked()
    {
        run_reverse_complement(app);
        ui.close_menu();
    }
    if ui
        .add_enabled(is_nuc, egui::Button::new("Translate to Protein"))
        .clicked()
    {
        run_translate(app);
        ui.close_menu();
    }

    ui.separator();

    ui.menu_button("Circular Display Options", |ui| {
        ui.checkbox(&mut app.circular.show_labels, "Outer Feature Labels");
        ui.checkbox(&mut app.circular.show_sequence_rings, "Aligned Strain Rings");
        ui.checkbox(&mut app.circular.show_identity_ring, "Pairwise Identity Ring");
        if is_nuc {
            ui.checkbox(&mut app.circular.show_gc_ring, "GC Content Ring");
        }
        ui.separator();
        if ui.button("Split Circular + Linear View").clicked() {
            app.display.show_circular_overview = true;
            app.display.show_linear_view = true;
            ui.close_menu();
        }
        if ui.button("Full Circular Phage Map").clicked() {
            app.display.show_circular_overview = true;
            app.display.show_linear_view = false;
            ui.close_menu();
        }
        if ui.button("Linear View Only").clicked() {
            app.display.show_circular_overview = false;
            app.display.show_linear_view = true;
            ui.close_menu();
        }
        ui.separator();
        if ui.button("Reset Circle Rotation & Zoom").clicked() {
            app.circular.rotation_rad = 0.0;
            app.circular.zoom = 1.0;
            ui.close_menu();
        }
    });
}

/// Renders the right-click context menu for a row in the upper Document Table.
pub fn show_document_row_context_menu(
    app: &mut HelixApp,
    ctx: &Context,
    ui: &mut Ui,
    doc_idx: usize,
) {
    let Some(doc) = app.docs.get(doc_idx) else {
        return;
    };
    let doc_name = doc.alignment.name.clone();
    let is_tree = doc.tree.is_some();
    let is_nuc = doc.alignment.alphabet.is_nucleotide();
    let is_deleted = doc.deleted;
    let is_shared = doc.shared;
    let uniprot_acc = doc
        .alignment
        .annotations()
        .iter()
        .find_map(|a| a.uniprot.as_ref().map(|u| u.accession.clone()))
        .or_else(|| {
            doc.alignment
                .sequences()
                .iter()
                .find_map(|s| s.uniprot_accession.clone())
        });

    ui.set_min_width(225.0);
    ui.label(RichText::new(&doc_name).strong().small());
    ui.separator();

    if ui.button(if is_tree { "Open in Tree View" } else { "Open in Sequence View" }).clicked() {
        app.select_document(doc_idx);
        app.active_tab = if is_tree {
            DocTab::TreeView
        } else {
            DocTab::SequenceView
        };
        ui.close_menu();
    }

    if ui.button("BLAST Document...").clicked() {
        app.select_document(doc_idx);
        open_blast_with_selection(app);
        ui.close_menu();
    }

    if ui.button("Align / Assemble...").clicked() {
        app.select_document(doc_idx);
        app.assemble_modal.open = true;
        ui.close_menu();
    }

    if ui.button("Build Phylogenetic Tree...").clicked() {
        app.select_document(doc_idx);
        app.tree_modal.open = true;
        ui.close_menu();
    }

    if ui.button("Design PCR Primers...").clicked() {
        app.select_document(doc_idx);
        app.primers_modal.open = true;
        ui.close_menu();
    }

    if ui.button("In-Silico Cloning...").clicked() {
        app.select_document(doc_idx);
        app.cloning_modal.open = true;
        ui.close_menu();
    }

    if ui
        .add_enabled(is_nuc, egui::Button::new("Auto-Annotate Reference Features"))
        .clicked()
    {
        app.select_document(doc_idx);
        run_auto_annotate(app);
        ui.close_menu();
    }

    if let Some(acc) = uniprot_acc {
        if ui.button(format!("Open 3D Protein Viewer ({acc})")).clicked() {
            app.select_document(doc_idx);
            app.open_3d_for_accession(ctx, &acc, None);
            ui.close_menu();
        }
    }

    ui.separator();

    if ui.button("Duplicate Document").clicked() {
        let copy_payload = app.docs.get(doc_idx).map(|src| {
            (
                src.folder.clone(),
                src.is_circular,
                src.alignment.clone(),
                src.tree.clone(),
                src.id,
                src.alignment.name.clone(),
            )
        });
        if let Some((folder, is_circular, mut aln, tree, parent_id, parent_name)) = copy_payload {
            aln.name = format!("{parent_name} (Copy)");
            app.add_workspace_doc(
                folder,
                is_circular,
                aln,
                tree,
                Some(LineageEntry {
                    operation: String::from("Duplicate Document"),
                    parent_ids: vec![parent_id],
                    parent_names: vec![parent_name],
                    details: String::from("Created editable copy via context menu"),
                }),
            );
        }
        ui.close_menu();
    }

    if ui.button("Export Document...").clicked() {
        app.select_document(doc_idx);
        app.export_alignment_dialog();
        ui.close_menu();
    }

    if ui.button("Import Pre-Aligned MAFFT File...").clicked() {
        app.import_mafft_dialog();
        ui.close_menu();
    }

    ui.separator();

    if !is_deleted {
        let share_txt = if is_shared {
            "Move to Local Workspace"
        } else {
            "Share to Team DB"
        };
        if ui.button(share_txt).clicked() {
            if let Some(d) = app.docs.get_mut(doc_idx) {
                d.shared = !is_shared;
            }
            ui.close_menu();
        }
        if ui.button("Move to Deleted Items").clicked() {
            if let Some(d) = app.docs.get_mut(doc_idx) {
                d.deleted = true;
            }
            app.status_message = format!("Moved '{doc_name}' to Deleted Items.");
            ui.close_menu();
        }
    } else if ui.button("Restore from Deleted Items").clicked() {
        if let Some(d) = app.docs.get_mut(doc_idx) {
            d.deleted = false;
        }
        app.status_message = format!("Restored '{doc_name}' from Deleted Items.");
        ui.close_menu();
    }
}

/// Renders the right-click context menu on the separate 3D Protein Viewer canvas.
pub fn show_3d_viewport_context_menu(app: &mut HelixApp, vp_ctx: &Context, ui: &mut Ui) {
    ui.set_min_width(215.0);

    let title = app
        .viewer3d
        .structure
        .as_ref()
        .map_or_else(|| String::from("3D Protein Viewer"), |s| s.entry_id.clone());
    ui.label(RichText::new(title).strong().small());
    ui.separator();

    // Jump to mapped alignment column if a residue is focused.
    let focus_res = app.viewer3d.hovered_residue.or(app.viewer3d.selected_residue);
    if let Some(res_idx) = focus_res {
        if let Some(Some((c0, c1))) = app
            .viewer3d
            .alignment_map
            .as_ref()
            .and_then(|m| m.residue_to_cols.get(res_idx))
            .copied()
        {
            if ui
                .button(format!("Jump Alignment to Col {}..{}", c0 + 1, c1))
                .clicked()
            {
                app.selected_cols = Some((c0, c1));
                app.cursor_col = c0;
                app.scroll_to_column = Some(c0);
                ui.close_menu();
            }
            ui.separator();
        }
    }

    ui.menu_button("Render Style", |ui| {
        for st in RenderStyle::ALL {
            if ui.selectable_label(app.viewer3d.style == st, st.label()).clicked() {
                app.viewer3d.style = st;
                ui.close_menu();
            }
        }
    });

    ui.menu_button("Color Scheme", |ui| {
        for cb in ColorBy3D::ALL {
            if ui
                .selectable_label(app.viewer3d.color_by == cb, cb.label())
                .clicked()
            {
                app.viewer3d.color_by = cb;
                ui.close_menu();
            }
        }
    });

    ui.menu_button("Fetch / Load Structure", |ui| {
        let acc = app.viewer3d.accession_input.clone();
        if ui.button("Fetch from AlphaFold DB").clicked() {
            app.open_3d_with_source(
                vp_ctx,
                &acc,
                ProteinStructureSource::AlphaFoldDb,
                app.viewer3d.linked_annotation_id,
            );
            ui.close_menu();
        }
        if ui.button("Fetch from UniProtKB (PDB / Meta)").clicked() {
            app.open_3d_with_source(
                vp_ctx,
                &acc,
                ProteinStructureSource::UniProtDb,
                app.viewer3d.linked_annotation_id,
            );
            ui.close_menu();
        }
        if ui.button("Auto Fetch (AlphaFold -> UniProt)").clicked() {
            app.open_3d_with_source(
                vp_ctx,
                &acc,
                ProteinStructureSource::AutoFallback,
                app.viewer3d.linked_annotation_id,
            );
            ui.close_menu();
        }
        ui.separator();
        if ui.button("Load Local PDB / mmCIF File...").clicked() {
            app.load_local_structure_dialog();
            ui.close_menu();
        }
        if ui.button("Load Offline Schematic Preview").clicked() {
            app.open_offline_3d_preview();
            ui.close_menu();
        }
    });

    ui.separator();

    ui.checkbox(&mut app.viewer3d.camera.auto_rotate, "Auto-Spin Structure");
    if ui.button("Reset Camera View").clicked() {
        app.viewer3d.camera.reset();
        ui.close_menu();
    }
}

// --- Shared action helpers used by context menus ---

fn open_blast_with_selection(app: &mut HelixApp) {
    if let Some(seq_text) = extract_selection_text(app, 180) {
        let ungapped: String = seq_text.chars().filter(|&c| c != '-' && c != '.').collect();
        if !ungapped.is_empty() {
            app.blast_modal.query = ungapped;
        }
    }
    app.blast_modal.open = true;
}

fn extract_selection_text(app: &HelixApp, default_window: usize) -> Option<String> {
    let aln = app.current_alignment()?;
    let seq = aln
        .sequences()
        .get(app.cursor_row)
        .or_else(|| aln.sequences().first())?;
    let len = seq.len();
    if len == 0 {
        return None;
    }
    let (s, e) = match app.selected_cols {
        Some((s, e)) if e > s + 1 => (s.min(len), e.min(len)),
        _ => {
            let start = app.cursor_col.min(len.saturating_sub(1));
            let end = (start + default_window).min(len);
            (start, end)
        }
    };
    Some(String::from_utf8_lossy(&seq.residues[s..e]).into_owned())
}

fn resolve_active_uniprot(app: &HelixApp) -> Option<String> {
    let aln = app.current_alignment()?;
    if let Some(ann_id) = app.selected_annotation_id {
        if let Some(acc) = aln
            .annotation_by_id(ann_id)
            .and_then(|a| a.uniprot.as_ref())
            .map(|u| u.accession.clone())
        {
            return Some(acc);
        }
    }
    aln.annotations()
        .iter()
        .find_map(|a| a.uniprot.as_ref().map(|u| u.accession.clone()))
        .or_else(|| {
            aln.sequences()
                .iter()
                .find_map(|s| s.uniprot_accession.clone())
        })
}

fn run_extract_region(app: &mut HelixApp) {
    let Some(doc) = app.current_doc() else {
        return;
    };
    let parent_id = doc.id;
    let parent_name = doc.alignment.name.clone();
    let width = doc.alignment.width();
    let (s, e) = match app.selected_cols {
        Some((s, e)) if e > s + 1 => (s, e),
        _ => {
            let start = app.cursor_col.min(width.saturating_sub(1));
            let end = (start + 100).min(width);
            (start, end)
        }
    };
    let new_name = format!("{parent_name}_{}..{}", s + 1, e);
    match doc.alignment.extract_region(new_name, s, e) {
        Ok(extracted) => {
            app.add_workspace_doc(
                "Alignments",
                false,
                extracted,
                None,
                Some(LineageEntry {
                    operation: String::from("Extract Region"),
                    parent_ids: vec![parent_id],
                    parent_names: vec![parent_name],
                    details: format!("Extracted bp {}..{}", s + 1, e),
                }),
            );
        }
        Err(err) => app.status_message = format!("Extract failed: {err}"),
    }
}

fn run_reverse_complement(app: &mut HelixApp) {
    let Some(doc) = app.current_doc() else {
        return;
    };
    let parent_id = doc.id;
    let parent_name = doc.alignment.name.clone();
    let folder = doc.folder.clone();
    let is_circ = doc.is_circular;
    match doc.alignment.reverse_complement_copy() {
        Ok(rc) => {
            app.add_workspace_doc(
                folder,
                is_circ,
                rc,
                None,
                Some(LineageEntry {
                    operation: String::from("Reverse Complement"),
                    parent_ids: vec![parent_id],
                    parent_names: vec![parent_name],
                    details: String::from("Generated 5'->3' reverse complement with flipped annotations"),
                }),
            );
        }
        Err(err) => app.status_message = format!("Reverse complement failed: {err}"),
    }
}

fn run_translate(app: &mut HelixApp) {
    let Some(doc) = app.current_doc() else {
        return;
    };
    let parent_id = doc.id;
    let parent_name = doc.alignment.name.clone();
    let frame = app.display.translation_frame;
    match doc.alignment.translate_to_protein_alignment(frame) {
        Ok(prot) => {
            app.add_workspace_doc(
                "Protein Documents",
                false,
                prot,
                None,
                Some(LineageEntry {
                    operation: format!("Translate DNA/RNA (Frame +{})", frame + 1),
                    parent_ids: vec![parent_id],
                    parent_names: vec![parent_name],
                    details: String::from("Translated using standard genetic code"),
                }),
            );
        }
        Err(err) => app.status_message = format!("Translation failed: {err}"),
    }
}

fn run_auto_annotate(app: &mut HelixApp) {
    if let Some(doc) = app.current_doc_mut() {
        let added = auto_annotate_from_reference_db(&mut doc.alignment, 90.0);
        app.status_message = format!(
            "Auto-Annotate scanned {} reference features and added {added} new annotation(s).",
            PLASMID_REFERENCE_FEATURES.len()
        );
    }
}
