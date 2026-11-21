//! Upper Geneious Prime-style Document Table pane displaying local folder items,
//! Search Everywhere matches, Reference Features, Deleted Items, Shared Databases,
//! Operations telemetry, and live NCBI / UniProtKB database search results.

use crate::app::{HelixApp, OperationStatus, SourceSelection};
use crate::theme::Palette;
use egui::{Color32, Context, RichText, ScrollArea, Ui};
use helix_core::features_db::{auto_annotate_from_reference_db, PLASMID_REFERENCE_FEATURES};
use helix_core::ncbi::NcbiDatabase;

/// Renders the upper document/database table pane corresponding to `app.source_selection`.
pub fn show_document_table(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    match app.source_selection.clone() {
        SourceSelection::ReferenceFeatures => show_reference_features_table(app, ui, &pal),
        SourceSelection::Operations => show_operations_table(app, ui, &pal),
        SourceSelection::Ncbi(db) => show_ncbi_table(app, ctx, ui, db, &pal),
        SourceSelection::UniProt => show_uniprot_table(app, ctx, ui, &pal),
        SourceSelection::DeletedItems => show_local_docs_table(app, ui, &pal, true, false, None),
        SourceSelection::SharedDatabases => {
            show_local_docs_table(app, ui, &pal, false, true, None);
        }
        SourceSelection::AllLocal => show_local_docs_table(app, ui, &pal, false, false, None),
        SourceSelection::Folder(folder) => {
            show_local_docs_table(app, ui, &pal, false, false, Some(folder));
        }
    }
}

fn show_local_docs_table(
    app: &mut HelixApp,
    ui: &mut Ui,
    pal: &Palette,
    deleted_only: bool,
    shared_only: bool,
    folder_filter: Option<String>,
) {
    let query = app.global_search.trim().to_ascii_lowercase();
    let matching_indices: Vec<usize> = app
        .docs
        .iter()
        .enumerate()
        .filter(|(_, d)| {
            if d.deleted != deleted_only {
                return false;
            }
            if !deleted_only && shared_only && !d.shared {
                return false;
            }
            if !deleted_only && !shared_only && d.shared {
                return false;
            }
            if query.is_empty() {
                if let Some(ref folder) = folder_filter {
                    if &d.folder != folder {
                        return false;
                    }
                }
                true
            } else {
                // Global Search Everywhere scans document name, description, folder, organism, and annotations!
                d.alignment.name.to_ascii_lowercase().contains(&query)
                    || d.alignment.description.to_ascii_lowercase().contains(&query)
                    || d.folder.to_ascii_lowercase().contains(&query)
                    || d.alignment
                        .annotations()
                        .iter()
                        .any(|a| a.name.to_ascii_lowercase().contains(&query))
            }
        })
        .map(|(i, _)| i)
        .collect();

    // Header strip showing selection count and quick actions (Delete / Restore / Share).
    ui.horizontal(|ui| {
        let heading = if !query.is_empty() {
            format!("Search Everywhere: '{}' ({} matches)", app.global_search, matching_indices.len())
        } else if deleted_only {
            format!("Deleted Items ({})", matching_indices.len())
        } else if shared_only {
            format!("Shared Databases ({})", matching_indices.len())
        } else if let Some(ref f) = folder_filter {
            format!("{f} ({})", matching_indices.len())
        } else {
            format!("All Local Documents ({})", matching_indices.len())
        };
        ui.label(RichText::new(heading).strong().small().color(pal.text_primary));

        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if deleted_only {
                if ui.small_button("Restore Selected").clicked() {
                    for &idx in &matching_indices {
                        if app.docs[idx].checked || idx == app.active_doc {
                            app.docs[idx].deleted = false;
                            app.docs[idx].checked = false;
                        }
                    }
                    app.status_message = String::from("Restored selected documents from Deleted Items.");
                }
                if ui.small_button("Empty Trash").clicked() {
                    app.docs.retain(|d| !d.deleted);
                    if app.active_doc >= app.docs.len() && !app.docs.is_empty() {
                        app.active_doc = 0;
                    }
                    app.status_message = String::from("Emptied Deleted Items.");
                }
            } else {
                if ui.small_button("Move to Deleted Items").clicked() {
                    let mut moved = 0usize;
                    for &idx in &matching_indices {
                        if app.docs[idx].checked || (idx == app.active_doc && matching_indices.len() == 1) {
                            app.docs[idx].deleted = true;
                            app.docs[idx].checked = false;
                            moved += 1;
                        }
                    }
                    if moved == 0 && app.active_doc < app.docs.len() {
                        app.docs[app.active_doc].deleted = true;
                        moved = 1;
                    }
                    app.status_message = format!("Moved {moved} item(s) to Deleted Items (recoverable).");
                }
                let share_label = if shared_only {
                    "Move to Local"
                } else {
                    "Share to Team DB"
                };
                if ui.small_button(share_label).clicked() {
                    for &idx in &matching_indices {
                        if app.docs[idx].checked || idx == app.active_doc {
                            app.docs[idx].shared = !shared_only;
                            app.docs[idx].checked = false;
                        }
                    }
                    app.status_message = if shared_only {
                        String::from("Moved selected item(s) back to Local.")
                    } else {
                        String::from("Shared selected item(s) to Shared Databases.")
                    };
                }
            }
        });
    });
    ui.separator();

    let mut select_idx: Option<usize> = None;

    ScrollArea::vertical().id_salt("local_doc_table_scroll").show(ui, |ui| {
        egui::Grid::new("local_docs_grid")
            .num_columns(7)
            .striped(true)
            .spacing([12.0, 4.0])
            .show(ui, |ui| {
                ui.label(RichText::new("").small());
                ui.label(RichText::new("Type").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Name").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Description").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Modified").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Organism").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Length").strong().small().color(pal.text_muted));
                ui.end_row();

                for idx in matching_indices {
                    let is_active = idx == app.active_doc;
                    let doc = &mut app.docs[idx];
                    ui.checkbox(&mut doc.checked, "");

                    let (type_badge, badge_col) = if doc.tree.is_some() {
                        ("TREE", Color32::from_rgb(186, 120, 48))
                    } else if doc.is_circular {
                        ("CIRC", Color32::from_rgb(38, 152, 92))
                    } else if doc.alignment.row_count() > 1 {
                        ("MSA", Color32::from_rgb(46, 124, 208))
                    } else {
                        (doc.alignment.alphabet.label(), Color32::from_rgb(132, 94, 212))
                    };
                    ui.label(RichText::new(type_badge).strong().size(10.0).color(badge_col));

                    if ui
                        .selectable_label(is_active, RichText::new(&doc.alignment.name).strong().size(12.0))
                        .clicked()
                    {
                        select_idx = Some(idx);
                    }

                    let short_desc: String = doc.alignment.description.chars().take(54).collect();
                    ui.label(RichText::new(short_desc).size(11.5).color(pal.text_muted));
                    ui.label(RichText::new(&doc.modified).size(11.0).color(pal.text_muted));

                    let org = doc
                        .alignment
                        .sequences()
                        .first()
                        .and_then(|s| s.organism.as_deref())
                        .unwrap_or("-");
                    ui.label(RichText::new(org).size(11.0).color(pal.text_muted));

                    let len_str = if doc.tree.is_some() {
                        format!("{} taxa", doc.alignment.row_count())
                    } else {
                        format!("{}", doc.alignment.width())
                    };
                    ui.label(RichText::new(len_str).monospace().size(11.0));
                    ui.end_row();
                }
            });
    });

    if let Some(idx) = select_idx {
        app.select_document(idx);
    }
}

fn show_reference_features_table(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    ui.horizontal(|ui| {
        ui.label(
            RichText::new(format!(
                "Helix Plasmid & Phage Reference Features ({} curated elements)",
                PLASMID_REFERENCE_FEATURES.len()
            ))
            .strong()
            .small()
            .color(Color32::from_rgb(38, 148, 86)),
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui
                .button("Auto-Annotate Active Document Now")
                .on_hover_text("Scan active nucleotide sequence against all built-in reference features")
                .clicked()
            {
                if let Some(doc) = app.current_doc_mut() {
                    let added = auto_annotate_from_reference_db(&mut doc.alignment, 90.0);
                    app.status_message = format!(
                        "Auto-Annotate scanned {} reference features and added {added} new annotation(s).",
                        PLASMID_REFERENCE_FEATURES.len()
                    );
                }
            }
        });
    });
    ui.separator();

    ScrollArea::vertical().id_salt("ref_features_scroll").show(ui, |ui| {
        egui::Grid::new("ref_features_grid")
            .num_columns(5)
            .striped(true)
            .spacing([12.0, 4.0])
            .show(ui, |ui| {
                ui.label(RichText::new("Feature").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Type").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Length").strong().small().color(pal.text_muted));
                ui.label(RichText::new("UniProt").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Description").strong().small().color(pal.text_muted));
                ui.end_row();

                for feat in PLASMID_REFERENCE_FEATURES {
                    ui.label(RichText::new(feat.name).strong().size(12.0));
                    ui.label(RichText::new(feat.kind.label()).size(11.0).color(pal.accent));
                    ui.label(RichText::new(format!("{} bp", feat.signature.len())).monospace().size(11.0));
                    ui.label(
                        RichText::new(feat.uniprot_accession.unwrap_or("-"))
                            .monospace()
                            .size(11.0),
                    );
                    ui.label(RichText::new(feat.description).size(11.0).color(pal.text_muted));
                    ui.end_row();
                }
            });
    });
}

fn show_operations_table(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    ui.horizontal(|ui| {
        ui.label(
            RichText::new(format!("Operations & Background Tasks ({})", app.operations.len()))
                .strong()
                .small(),
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.small_button("Clear Completed").clicked() {
                app.operations
                    .retain(|o| o.status == OperationStatus::Running);
            }
        });
    });
    ui.separator();

    ScrollArea::vertical().id_salt("operations_scroll").show(ui, |ui| {
        egui::Grid::new("ops_grid")
            .num_columns(4)
            .striped(true)
            .spacing([14.0, 4.0])
            .show(ui, |ui| {
                ui.label(RichText::new("Job").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Operation").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Status").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Details").strong().small().color(pal.text_muted));
                ui.end_row();

                for op in &app.operations {
                    ui.label(RichText::new(&op.started_label).monospace().size(11.0));
                    ui.label(RichText::new(&op.name).strong().size(12.0));
                    let (st_txt, st_col) = match op.status {
                        OperationStatus::Running => ("Running...", pal.accent),
                        OperationStatus::Completed => ("Completed", pal.identity_high),
                        OperationStatus::Failed => ("Failed", pal.identity_low),
                    };
                    ui.label(RichText::new(st_txt).strong().size(11.0).color(st_col));
                    ui.label(RichText::new(&op.detail).size(11.0).color(pal.text_muted));
                    ui.end_row();
                }
            });
    });
}

fn show_ncbi_table(
    app: &mut HelixApp,
    ctx: &Context,
    ui: &mut Ui,
    db: NcbiDatabase,
    pal: &Palette,
) {
    let mut trigger = false;
    ui.horizontal(|ui| {
        ui.label(
            RichText::new(format!("NCBI {}", db.label()))
                .strong()
                .color(pal.accent),
        );
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.ncbi_search_query)
                .desired_width(280.0)
                .hint_text("Search NCBI (e.g. Escherichia phage lambda, TP53, pUC19)..."),
        );
        if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            trigger = true;
        }
        if ui.button("Search NCBI").clicked() {
            trigger = true;
        }
        if app.ncbi_searching {
            ui.spinner();
            ui.label(RichText::new("Querying E-utilities...").small());
        }
    });
    if trigger {
        app.trigger_ncbi_search(ctx, db);
    }

    if let Some(err) = &app.ncbi_search_error {
        ui.label(RichText::new(err).small().color(pal.identity_low));
    }
    ui.separator();

    let mut download_hit: Option<helix_core::NcbiSearchHit> = None;
    let mut open_pdb_3d: Option<String> = None;

    ScrollArea::vertical().id_salt("ncbi_results_scroll").show(ui, |ui| {
        egui::Grid::new("ncbi_grid")
            .num_columns(6)
            .striped(true)
            .spacing([10.0, 4.0])
            .show(ui, |ui| {
                ui.label(RichText::new("Accession / ID").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Title").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Organism").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Length").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Updated").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Action").strong().small().color(pal.text_muted));
                ui.end_row();

                for hit in &app.ncbi_search_results {
                    ui.label(RichText::new(&hit.accession).monospace().strong().size(11.5));
                    let title_short: String = hit.title.chars().take(65).collect();
                    ui.label(RichText::new(title_short).size(11.5));
                    ui.label(RichText::new(&hit.organism).size(11.0).color(pal.text_muted));
                    ui.label(
                        RichText::new(if hit.length > 0 {
                            hit.length.to_string()
                        } else {
                            "-".to_owned()
                        })
                        .monospace()
                        .size(11.0),
                    );
                    ui.label(RichText::new(&hit.updated).size(10.5).color(pal.text_muted));
                    if hit.database.supports_sequence_download() {
                        if ui.small_button("Download to Local").clicked() {
                            download_hit = Some(hit.clone());
                        }
                    } else if hit.database == NcbiDatabase::Structure && hit.accession.len() == 4 {
                        if ui.small_button("View 3D (RCSB)").clicked() {
                            open_pdb_3d = Some(hit.accession.clone());
                        }
                    } else {
                        ui.label(RichText::new("Metadata").small().color(pal.text_muted));
                    }
                    ui.end_row();
                }
            });
    });

    if let Some(hit) = download_hit {
        app.download_ncbi_hit(ctx, &hit);
    }
    if let Some(pdb) = open_pdb_3d {
        app.fetch_rcsb_structure(ctx, &pdb);
    }
}

fn show_uniprot_table(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let mut trigger = false;
    ui.horizontal(|ui| {
        ui.label(RichText::new("UniProtKB").strong().color(pal.accent));
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.uniprot_search_query)
                .desired_width(280.0)
                .hint_text("Search UniProtKB (e.g. TP53 human, T4 lysozyme, GFP)..."),
        );
        if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            trigger = true;
        }
        if ui.button("Search UniProt").clicked() {
            trigger = true;
        }
        if app.uniprot_searching {
            ui.spinner();
            ui.label(RichText::new("Searching UniProtKB...").small());
        }
    });
    if trigger {
        app.trigger_uniprot_search(ctx);
    }
    if let Some(err) = &app.uniprot_search_error {
        ui.label(RichText::new(err).small().color(pal.identity_low));
    }
    ui.separator();

    let mut download_acc: Option<String> = None;
    let mut view_3d_acc: Option<String> = None;

    ScrollArea::vertical().id_salt("uniprot_table_scroll").show(ui, |ui| {
        egui::Grid::new("uniprot_tbl_grid")
            .num_columns(6)
            .striped(true)
            .spacing([10.0, 4.0])
            .show(ui, |ui| {
                ui.label(RichText::new("Accession").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Entry").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Protein Name").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Organism").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Length").strong().small().color(pal.text_muted));
                ui.label(RichText::new("Actions").strong().small().color(pal.text_muted));
                ui.end_row();

                for hit in &app.uniprot_search_results {
                    ui.label(RichText::new(&hit.accession).monospace().strong().color(pal.accent));
                    ui.label(RichText::new(&hit.entry_name).monospace().size(11.0));
                    let pname: String = hit.protein_name.chars().take(52).collect();
                    ui.label(RichText::new(pname).size(11.5));
                    ui.label(RichText::new(&hit.organism).size(11.0).color(pal.text_muted));
                    ui.label(RichText::new(format!("{} aa", hit.length)).monospace().size(11.0));
                    ui.horizontal(|ui| {
                        if ui.small_button("Download").clicked() {
                            download_acc = Some(hit.accession.clone());
                        }
                        if ui.small_button("3D AlphaFold").clicked() {
                            view_3d_acc = Some(hit.accession.clone());
                        }
                    });
                    ui.end_row();
                }
            });
    });

    if let Some(acc) = download_acc {
        app.download_uniprot_as_doc(ctx, &acc);
    }
    if let Some(acc) = view_3d_acc {
        app.open_3d_for_accession(ctx, &acc, None);
    }
}
