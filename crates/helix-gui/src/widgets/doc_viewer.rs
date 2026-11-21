//! Lower Geneious Prime Document Viewer pane with tabs:
//! `Sequence View` | `Annotations` | `Text View` | `Lineage` | `Info` | `Tree View`.
//!
//! Includes the Geneious Sequence View manipulation sub-toolbar (`Extract`, `R.C.`, `Translate`,
//! `Add/Edit Annotation`, `Auto-Annotate`, `Allow Editing`) and supports simultaneous
//! Circular Overview + Linear Zoom split rendering (matching the Geneious Prime screenshot).

use crate::app::{DocTab, HelixApp, LineageEntry};
use crate::theme::Palette;
use crate::widgets::{alignment_view, circular_view, tree_view};
use egui::{Button, Color32, Context, Frame, RichText, ScrollArea, Stroke, Ui};
use helix_core::features_db::auto_annotate_from_reference_db;
use helix_core::io::fasta::write_fasta;

/// Renders the lower tabbed Document Viewer.
pub fn show_doc_viewer(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    if app.current_doc().is_none() {
        ui.centered_and_justified(|ui| {
            ui.label(RichText::new("Select a document from the table above.").color(pal.text_muted));
        });
        return;
    }

    let has_tree = app.current_doc().is_some_and(|d| d.tree.is_some());

    // 1. Document Tab Strip (`Sequence View`, `Annotations`, `Text View`, `Lineage`, `Info`, `Tree View`).
    Frame::none()
        .fill(pal.header_bg)
        .stroke(Stroke::new(1.0, pal.border))
        .inner_margin(egui::Margin::symmetric(10.0, 4.0))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut app.active_tab, DocTab::SequenceView, "Sequence View");
                ui.selectable_value(&mut app.active_tab, DocTab::Annotations, "Annotations");
                ui.selectable_value(&mut app.active_tab, DocTab::TextView, "Text View");
                ui.selectable_value(&mut app.active_tab, DocTab::Lineage, "Lineage");
                ui.selectable_value(&mut app.active_tab, DocTab::Info, "Info");
                if has_tree {
                    ui.selectable_value(&mut app.active_tab, DocTab::TreeView, "Tree View");
                }
            });
        });

    // 2. Active Tab Body.
    match app.active_tab {
        DocTab::SequenceView => show_sequence_view_tab(app, ctx, ui, &pal),
        DocTab::Annotations => show_annotations_table_tab(app, ctx, ui, &pal),
        DocTab::TextView => show_text_view_tab(app, ui, &pal),
        DocTab::Lineage => show_lineage_tab(app, ui, &pal),
        DocTab::Info => show_info_tab(app, ui, &pal),
        DocTab::TreeView => tree_view::show_tree_view(app, ui),
    }
}

fn show_sequence_view_tab(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    // Sequence manipulation sub-toolbar matching Geneious Prime (`Extract`, `R.C.`, `Translate`,
    // `Add/Edit Annotation`, `Auto-Annotate`, `Allow Editing`).
    show_sequence_subtoolbar(app, ctx, ui, pal);

    // Ensure at least one of Circular Overview or Linear View is visible.
    if !app.display.show_circular_overview && !app.display.show_linear_view {
        app.display.show_linear_view = true;
    }

    match (app.display.show_circular_overview, app.display.show_linear_view) {
        // Iconic Geneious Prime Split View: Circular Overview on the left, Linear Zoom on the right!
        (true, true) => {
            let split_w = (ui.available_width() * 0.40).clamp(240.0, 520.0);
            egui::SidePanel::left("helix_circular_split_left")
                .resizable(true)
                .default_width(split_w)
                .min_width(210.0)
                .max_width(680.0)
                .frame(
                    Frame::none()
                        .fill(pal.canvas_bg)
                        .stroke(Stroke::new(1.0, pal.border)),
                )
                .show_inside(ui, |ui| {
                    circular_view::show_circular_canvas_embedded(app, ctx, ui, pal);
                });

            alignment_view::show_alignment_view_body(app, ctx, ui);
        }
        // Full Circular Phage / Plasmid Map.
        (true, false) => {
            circular_view::show_circular_view(app, ctx, ui);
        }
        // Full Linear or Wrapped Sequence Alignment View.
        (false, true) => {
            alignment_view::show_alignment_view_body(app, ctx, ui);
        }
        (false, false) => {}
    }
}

fn show_sequence_subtoolbar(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let is_nuc = app
        .current_alignment()
        .is_some_and(|a| a.alphabet.is_nucleotide());
    let mut launch_3d_req: Option<(String, Option<u64>)> = None;

    Frame::none()
        .fill(pal.panel_bg)
        .stroke(Stroke::new(1.0, pal.border))
        .inner_margin(egui::Margin::symmetric(10.0, 4.0))
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                // Step between annotations / features.
                if ui
                    .small_button("←")
                    .on_hover_text("Jump to previous annotation")
                    .clicked()
                {
                    jump_relative_annotation(app, false);
                }
                if ui
                    .small_button("→")
                    .on_hover_text("Jump to next annotation")
                    .clicked()
                {
                    jump_relative_annotation(app, true);
                }
                ui.separator();

                // Extract selected region into a new workspace document with Lineage!
                if ui
                    .small_button("Extract")
                    .on_hover_text("Extract selected column region into a new sequence document")
                    .clicked()
                {
                    if let Some(doc) = app.current_doc() {
                        let parent_id = doc.id;
                        let parent_name = doc.alignment.name.clone();
                        let (s, e) = app
                            .selected_cols
                            .unwrap_or((0, doc.alignment.width().min(100)));
                        let new_name = format!("{}_{}..{}", parent_name, s + 1, e);
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
                }

                // Reverse Complement (R.C.)
                if ui
                    .add_enabled(is_nuc, Button::new(RichText::new("R.C.").small()))
                    .on_hover_text("Create reverse-complement document with flipped annotations & lineage")
                    .clicked()
                {
                    if let Some(doc) = app.current_doc() {
                        let parent_id = doc.id;
                        let parent_name = doc.alignment.name.clone();
                        let folder = doc.folder.clone();
                        let circ = doc.is_circular;
                        match doc.alignment.reverse_complement_copy() {
                            Ok(rc) => {
                                app.add_workspace_doc(
                                    folder,
                                    circ,
                                    rc,
                                    None,
                                    Some(LineageEntry {
                                        operation: String::from("Reverse Complement (R.C.)"),
                                        parent_ids: vec![parent_id],
                                        parent_names: vec![parent_name],
                                        details: String::from("Generated 3'->5' reverse complement"),
                                    }),
                                );
                            }
                            Err(err) => app.status_message = format!("R.C. failed: {err}"),
                        }
                    }
                }

                // Translate DNA/RNA -> Protein document with Lineage!
                if ui
                    .add_enabled(is_nuc, Button::new(RichText::new("Translate").small()))
                    .on_hover_text("Translate sequence into a new Protein document in Protein Documents")
                    .clicked()
                {
                    if let Some(doc) = app.current_doc() {
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
                                        operation: format!("Translate (Frame +{})", frame + 1),
                                        parent_ids: vec![parent_id],
                                        parent_names: vec![parent_name],
                                        details: String::from("Standard NCBI Genetic Code 1"),
                                    }),
                                );
                            }
                            Err(err) => app.status_message = format!("Translate failed: {err}"),
                        }
                    }
                }

                // Add/Edit Annotation
                if ui
                    .small_button("Add/Edit Annotation")
                    .on_hover_text("Add a new feature annotation at the selected range")
                    .clicked()
                {
                    if let Some((s, e)) = app.selected_cols {
                        app.new_ann_draft.start_1based = s + 1;
                        app.new_ann_draft.end_1based = e;
                    } else {
                        app.new_ann_draft.start_1based = app.cursor_col + 1;
                        app.new_ann_draft.end_1based = app.cursor_col + 30;
                    }
                    app.new_ann_draft.open = true;
                }

                // Auto-Annotate from built-in Reference Plasmid/Phage Features!
                if ui
                    .add_enabled(is_nuc, Button::new(RichText::new("Auto-Annotate").small()))
                    .on_hover_text("Automatically identify promoters, origins, markers, tags & att sites from Reference Features")
                    .clicked()
                {
                    if let Some(doc) = app.current_doc_mut() {
                        let added = auto_annotate_from_reference_db(&mut doc.alignment, 90.0);
                        app.status_message = format!("Auto-Annotate added {added} feature(s) from Helix Plasmid Features.");
                    }
                }

                ui.separator();

                // Allow Editing toggle (matches Geneious Prime's Allow Editing button!).
                let edit_btn = if app.display.allow_editing {
                    Button::new(RichText::new("Editing: ON").small().strong().color(Color32::WHITE))
                        .fill(Color32::from_rgb(212, 108, 42))
                } else {
                    Button::new(RichText::new("Allow Editing").small())
                };
                if ui
                    .add(edit_btn)
                    .on_hover_text("Enable direct keyboard editing in Linear View (type A/C/G/T, Backspace, Delete)")
                    .clicked()
                {
                    app.display.allow_editing = !app.display.allow_editing;
                    app.status_message = if app.display.allow_editing {
                        String::from("Sequence Editing ENABLED — click a base in Linear View and type to edit or press Backspace/Delete.")
                    } else {
                        String::from("Sequence Editing locked (read-only).")
                    };
                }

                // Right side: selected annotation AlphaFold 3D quick button.
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if let Some(ann_id) = app.selected_annotation_id {
                        if let Some(ann) = app
                            .current_alignment()
                            .and_then(|a| a.annotation_by_id(ann_id))
                        {
                            if let Some(up) = &ann.uniprot {
                                let btn = Button::new(
                                    RichText::new(format!("3D AlphaFold ({})", up.accession))
                                        .small()
                                        .strong()
                                        .color(Color32::WHITE),
                                )
                                .fill(Color32::from_rgb(16, 138, 134));
                                if ui.add(btn).clicked() {
                                    launch_3d_req = Some((up.accession.clone(), Some(ann.id)));
                                }
                            }
                        }
                    }
                });
            });
        });

    if let Some((acc, ann_id)) = launch_3d_req {
        app.open_3d_for_accession(ctx, &acc, ann_id);
    }
}

fn jump_relative_annotation(app: &mut HelixApp, forward: bool) {
    let Some(aln) = app.current_alignment() else {
        return;
    };
    let anns = aln.annotations();
    if anns.is_empty() {
        return;
    }
    let cur_idx = app
        .selected_annotation_id
        .and_then(|id| anns.iter().position(|a| a.id == id))
        .unwrap_or(0);
    let next_idx = if forward {
        (cur_idx + 1) % anns.len()
    } else if cur_idx == 0 {
        anns.len() - 1
    } else {
        cur_idx - 1
    };
    let ann = &anns[next_idx];
    let (id, start, end) = (ann.id, ann.start, ann.end);
    app.selected_annotation_id = Some(id);
    app.selected_cols = Some((start, end));
    app.cursor_col = start;
    app.scroll_to_column = Some(start);
}

/// Full searchable Annotations table tab (matching Geneious Prime's Annotations tab).
fn show_annotations_table_tab(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        ui.add_space(10.0);
        ui.label("Filter:");
        ui.add(
            egui::TextEdit::singleline(&mut app.annotation_filter)
                .desired_width(220.0)
                .hint_text("Filter by name, type, UniProt, or notes..."),
        );
        if ui.button("+ Add Annotation").clicked() {
            app.new_ann_draft.open = true;
        }
    });
    ui.separator();

    let filter = app.annotation_filter.trim().to_ascii_lowercase();
    let mut jump_ann: Option<(u64, usize, usize)> = None;
    let mut open_3d: Option<(String, u64)> = None;
    let mut delete_id: Option<u64> = None;

    if let Some(aln) = app.current_alignment() {
        ScrollArea::vertical().id_salt("ann_tab_scroll").show(ui, |ui| {
            ui.add_space(4.0);
            egui::Grid::new("ann_full_table")
                .num_columns(8)
                .striped(true)
                .spacing([14.0, 5.0])
                .show(ui, |ui| {
                    for hdr in ["Name", "Type", "Min", "Max", "Length", "Strand", "UniProt / 3D", "Notes"] {
                        ui.label(RichText::new(hdr).strong().small().color(pal.text_muted));
                    }
                    ui.end_row();

                    for ann in aln.annotations() {
                        if !filter.is_empty()
                            && !ann.name.to_ascii_lowercase().contains(&filter)
                            && !ann.kind.label().to_ascii_lowercase().contains(&filter)
                            && !ann.notes.to_ascii_lowercase().contains(&filter)
                            && !ann
                                .uniprot
                                .as_ref()
                                .is_some_and(|u| u.accession.to_ascii_lowercase().contains(&filter))
                        {
                            continue;
                        }
                        let is_sel = app.selected_annotation_id == Some(ann.id);
                        if ui.selectable_label(is_sel, RichText::new(&ann.name).strong()).clicked() {
                            jump_ann = Some((ann.id, ann.start, ann.end));
                        }
                        ui.label(ann.kind.label());
                        ui.label(RichText::new((ann.start + 1).to_string()).monospace());
                        ui.label(RichText::new(ann.end.to_string()).monospace());
                        ui.label(RichText::new(ann.len().to_string()).monospace());
                        ui.label(ann.strand.symbol().to_string());
                        if let Some(up) = &ann.uniprot {
                            if ui.small_button(format!("3D: {}", up.accession)).clicked() {
                                open_3d = Some((up.accession.clone(), ann.id));
                            }
                        } else {
                            ui.label("-");
                        }
                        ui.horizontal(|ui| {
                            let short: String = ann.notes.chars().take(55).collect();
                            ui.label(RichText::new(short).small().color(pal.text_muted));
                            if ui.small_button("X").on_hover_text("Delete annotation").clicked() {
                                delete_id = Some(ann.id);
                            }
                        });
                        ui.end_row();
                    }
                });
        });
    }

    if let Some((id, s, e)) = jump_ann {
        app.selected_annotation_id = Some(id);
        app.selected_cols = Some((s, e));
        app.cursor_col = s;
        app.scroll_to_column = Some(s);
        app.active_tab = DocTab::SequenceView;
    }
    if let Some((acc, id)) = open_3d {
        app.open_3d_for_accession(ctx, &acc, Some(id));
    }
    if let Some(id) = delete_id {
        if let Some(doc) = app.current_doc_mut() {
            doc.alignment.remove_annotation(id);
        }
    }
}

/// Raw FASTA Text View tab (matching Geneious Prime's Text View tab).
fn show_text_view_tab(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    let Some(aln) = app.current_alignment() else {
        return;
    };
    let mut fasta_text = write_fasta(aln);

    ui.add_space(6.0);
    ui.horizontal(|ui| {
        ui.add_space(10.0);
        ui.label(RichText::new("FASTA Text View (Read-only preview)").strong().small().color(pal.text_muted));
        if ui.small_button("Copy FASTA").clicked() {
            ui.ctx().copy_text(fasta_text.clone());
            app.status_message = String::from("Copied FASTA text to clipboard.");
        }
        if ui.small_button("Export File...").clicked() {
            app.export_alignment_dialog();
        }
    });
    ui.separator();

    ScrollArea::vertical().id_salt("text_view_scroll").show(ui, |ui| {
        ui.add(
            egui::TextEdit::multiline(&mut fasta_text)
                .font(egui::TextStyle::Monospace)
                .desired_width(f32::INFINITY)
                .interactive(false),
        );
    });
}

/// Interactive Lineage derivation tab showing parent constructs & child descendants!
fn show_lineage_tab(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    let Some(doc) = app.current_doc() else {
        return;
    };
    let doc_id = doc.id;
    let doc_name = doc.alignment.name.clone();
    let lineage = doc.lineage.clone();

    // Find any child documents in the workspace derived from `doc_id`.
    let children: Vec<(usize, String, String)> = app
        .docs
        .iter()
        .enumerate()
        .filter(|(_, d)| !d.deleted && d.lineage.iter().any(|l| l.parent_ids.contains(&doc_id)))
        .map(|(i, d)| {
            let op = d
                .lineage
                .first()
                .map(|l| l.operation.clone())
                .unwrap_or_default();
            (i, d.alignment.name.clone(), op)
        })
        .collect();

    let mut jump_to_doc: Option<usize> = None;

    ScrollArea::vertical().show(ui, |ui| {
        ui.add_space(10.0);
        ui.horizontal(|ui| {
            ui.add_space(12.0);
            ui.vertical(|ui| {
                ui.label(RichText::new("DOCUMENT LINEAGE & PROVENANCE").strong().color(pal.accent));
                ui.add_space(8.0);

                // Ancestors
                if lineage.is_empty() {
                    ui.label(
                        RichText::new("Original source document (no parent derivation recorded). Perform Cloning, Extract, R.C., Translate, Assembly, or Tree Building to create tracked lineage!")
                            .color(pal.text_muted),
                    );
                } else {
                    for entry in &lineage {
                        Frame::none()
                            .fill(pal.panel_bg)
                            .stroke(Stroke::new(1.0, pal.border))
                            .rounding(egui::Rounding::same(6.0))
                            .inner_margin(egui::Margin::same(10.0))
                            .show(ui, |ui| {
                                ui.label(RichText::new(&entry.operation).strong().color(pal.accent));
                                ui.label(RichText::new(&entry.details).small().color(pal.text_muted));
                                ui.add_space(4.0);
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(RichText::new("Derived from:").small().strong());
                                    for (pid, pname) in entry.parent_ids.iter().zip(&entry.parent_names) {
                                        if ui.small_button(pname).clicked() {
                                            if let Some(pidx) = app.docs.iter().position(|d| d.id == *pid) {
                                                jump_to_doc = Some(pidx);
                                            }
                                        }
                                    }
                                });
                            });
                        ui.add_space(6.0);
                    }
                }

                // Current Document Node
                ui.add_space(6.0);
                Frame::none()
                    .fill(pal.accent_soft)
                    .stroke(Stroke::new(1.5, pal.accent))
                    .rounding(egui::Rounding::same(6.0))
                    .inner_margin(egui::Margin::same(10.0))
                    .show(ui, |ui| {
                        ui.label(RichText::new(format!("Current Document: {doc_name}")).strong());
                    });

                // Derived Descendants
                if !children.is_empty() {
                    ui.add_space(10.0);
                    ui.label(RichText::new("Derived Descendant Documents:").strong().color(pal.text_muted));
                    ui.add_space(4.0);
                    for (cidx, cname, op) in children {
                        ui.horizontal(|ui| {
                            if ui.button(&cname).clicked() {
                                jump_to_doc = Some(cidx);
                            }
                            ui.label(RichText::new(format!("via {op}")).small().color(pal.text_muted));
                        });
                    }
                }
            });
        });
    });

    if let Some(idx) = jump_to_doc {
        app.select_document(idx);
    }
}

/// Editable document metadata & physicochemical Info tab.
fn show_info_tab(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    let Some(doc) = app.current_doc_mut() else {
        return;
    };

    ScrollArea::vertical().show(ui, |ui| {
        ui.add_space(10.0);
        ui.horizontal(|ui| {
            ui.add_space(14.0);
            egui::Grid::new("doc_info_grid")
                .num_columns(2)
                .spacing([18.0, 8.0])
                .show(ui, |ui| {
                    ui.label(RichText::new("Document Name").strong().color(pal.text_muted));
                    ui.text_edit_singleline(&mut doc.alignment.name);
                    ui.end_row();

                    ui.label(RichText::new("Description").strong().color(pal.text_muted));
                    ui.add(
                        egui::TextEdit::multiline(&mut doc.alignment.description)
                            .desired_rows(2)
                            .desired_width(420.0),
                    );
                    ui.end_row();

                    ui.label(RichText::new("Folder").strong().color(pal.text_muted));
                    ui.label(&doc.folder);
                    ui.end_row();

                    ui.label(RichText::new("Topology").strong().color(pal.text_muted));
                    ui.horizontal(|ui| {
                        if ui.radio_value(&mut doc.is_circular, false, "Linear").changed() {
                            doc.invalidate_caches();
                        }
                        if ui
                            .radio_value(&mut doc.is_circular, true, "Circular (Plasmid / Phage)")
                            .changed()
                        {
                            doc.invalidate_caches();
                        }
                    });
                    ui.end_row();

                    ui.label(RichText::new("Alphabet").strong().color(pal.text_muted));
                    ui.label(doc.alignment.alphabet.label());
                    ui.end_row();

                    ui.label(RichText::new("Sequences / Width").strong().color(pal.text_muted));
                    ui.label(format!(
                        "{} sequence(s), {} columns",
                        doc.alignment.row_count(),
                        doc.alignment.width()
                    ));
                    ui.end_row();

                    ui.label(RichText::new("Mean Pairwise ID").strong().color(pal.text_muted));
                    ui.label(format!("{:.1}%", doc.alignment.mean_pairwise_identity_pct()));
                    ui.end_row();

                    if doc.alignment.alphabet.is_nucleotide() {
                        if let Some(gc) = doc.alignment.sequences().first().and_then(|s| s.gc_content_pct()) {
                            ui.label(RichText::new("GC Content").strong().color(pal.text_muted));
                            ui.label(format!("{gc:.1}%"));
                            ui.end_row();
                        }
                        let mw_kda = doc.alignment.width() as f32 * 0.65;
                        ui.label(RichText::new("Est. dsDNA MW").strong().color(pal.text_muted));
                        ui.label(format!("{mw_kda:.1} kDa"));
                        ui.end_row();
                    } else {
                        let mw_kda = doc.alignment.width() as f32 * 0.11;
                        ui.label(RichText::new("Est. Protein MW").strong().color(pal.text_muted));
                        ui.label(format!("{mw_kda:.2} kDa"));
                        ui.end_row();
                    }
                });
        });
    });
}
