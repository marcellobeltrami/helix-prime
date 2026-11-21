//! Left Geneious Prime-inspired Sources tree: Local folders, Reference Features,
//! Deleted Items recovery, Shared Databases, Operations monitor, NCBI (8 sub-databases),
//! and UniProtKB.

use crate::app::{HelixApp, OperationStatus, SourceSelection};
use crate::theme::Palette;
use egui::{Color32, Context, RichText, ScrollArea, Ui};
use helix_core::features_db::PLASMID_REFERENCE_FEATURES;
use helix_core::ncbi::NcbiDatabase;

const STANDARD_FOLDERS: &[&str] = &[
    "Alignments",
    "Cloning",
    "Contig Assembly",
    "Genomes",
    "Plasmids",
    "Primers",
    "Protein Documents",
    "Tree Documents",
];

/// Renders the left Sources hierarchy sidebar matching the Geneious Prime reference layout.
pub fn show_sources_panel(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    let local_total = app.docs.iter().filter(|d| !d.deleted && !d.shared).count();
    let deleted_total = app.docs.iter().filter(|d| d.deleted).count();
    let shared_total = app.docs.iter().filter(|d| !d.deleted && d.shared).count();
    let running_ops = app
        .operations
        .iter()
        .filter(|o| o.status == OperationStatus::Running)
        .count();

    // Collect any extra custom/imported folders beyond the standard list.
    let mut extra_folders: Vec<String> = Vec::new();
    for doc in &app.docs {
        if !doc.deleted
            && !doc.shared
            && !STANDARD_FOLDERS.iter().any(|&f| f == doc.folder)
            && !extra_folders.iter().any(|f| f == &doc.folder)
        {
            extra_folders.push(doc.folder.clone());
        }
    }

    let mut next_source: Option<SourceSelection> = None;

    ScrollArea::vertical().show(ui, |ui| {
        // --- 1. Local Database Tree ---
        egui::CollapsingHeader::new(
            RichText::new(format!("Local  ({local_total})"))
                .strong()
                .color(pal.text_primary),
        )
        .default_open(true)
        .show(ui, |ui| {
            if tree_row(
                ui,
                app.source_selection == SourceSelection::AllLocal,
                "All Local Documents",
                Some(local_total),
                pal.text_primary,
                &pal,
            ) {
                next_source = Some(SourceSelection::AllLocal);
            }

            egui::CollapsingHeader::new(
                RichText::new("Sample Documents")
                    .strong()
                    .color(pal.text_primary),
            )
            .default_open(true)
            .show(ui, |ui| {
                for &folder in STANDARD_FOLDERS {
                    let count = app
                        .docs
                        .iter()
                        .filter(|d| !d.deleted && !d.shared && d.folder == folder)
                        .count();
                    let sel = app.source_selection == SourceSelection::Folder(folder.to_owned());
                    if tree_row(ui, sel, folder, Some(count), pal.text_primary, &pal) {
                        next_source = Some(SourceSelection::Folder(folder.to_owned()));
                    }
                }
                for folder in &extra_folders {
                    let count = app
                        .docs
                        .iter()
                        .filter(|d| !d.deleted && !d.shared && &d.folder == folder)
                        .count();
                    let sel = app.source_selection == SourceSelection::Folder(folder.clone());
                    if tree_row(ui, sel, folder, Some(count), pal.text_primary, &pal) {
                        next_source = Some(SourceSelection::Folder(folder.clone()));
                    }
                }
            });

            // Reference Features (Automated Plasmid/Phage Annotation Library).
            egui::CollapsingHeader::new(
                RichText::new("Reference Features")
                    .strong()
                    .color(Color32::from_rgb(38, 148, 86)),
            )
            .default_open(true)
            .show(ui, |ui| {
                let sel = app.source_selection == SourceSelection::ReferenceFeatures;
                if tree_row(
                    ui,
                    sel,
                    "Helix Plasmid Features",
                    Some(PLASMID_REFERENCE_FEATURES.len()),
                    Color32::from_rgb(38, 148, 86),
                    &pal,
                ) {
                    next_source = Some(SourceSelection::ReferenceFeatures);
                }
            });

            // Deleted Items recovery bin.
            let del_sel = app.source_selection == SourceSelection::DeletedItems;
            if tree_row(
                ui,
                del_sel,
                "Deleted Items",
                Some(deleted_total),
                pal.text_muted,
                &pal,
            ) {
                next_source = Some(SourceSelection::DeletedItems);
            }
        });

        ui.add_space(4.0);

        // --- 2. Shared Databases ---
        let shared_sel = app.source_selection == SourceSelection::SharedDatabases;
        if tree_row(
            ui,
            shared_sel,
            "Shared Databases",
            Some(shared_total),
            pal.text_primary,
            &pal,
        ) {
            next_source = Some(SourceSelection::SharedDatabases);
        }

        // --- 3. Operations Monitor ---
        let ops_sel = app.source_selection == SourceSelection::Operations;
        let ops_label = if running_ops > 0 {
            format!("Operations ({running_ops} active)")
        } else {
            String::from("Operations")
        };
        if tree_row(
            ui,
            ops_sel,
            &ops_label,
            Some(app.operations.len()),
            if running_ops > 0 {
                pal.accent
            } else {
                pal.text_primary
            },
            &pal,
        ) {
            next_source = Some(SourceSelection::Operations);
        }

        ui.add_space(4.0);
        ui.separator();

        // --- 4. NCBI Public Databases (all 8 databases from screenshot) ---
        egui::CollapsingHeader::new(
            RichText::new("NCBI")
                .strong()
                .color(pal.text_primary),
        )
        .default_open(true)
        .show(ui, |ui| {
            for db in NcbiDatabase::ALL {
                let sel = app.source_selection == SourceSelection::Ncbi(db);
                if tree_row(ui, sel, db.label(), None, pal.text_primary, &pal) {
                    next_source = Some(SourceSelection::Ncbi(db));
                }
            }
        });

        // --- 5. UniProtKB Public Database ---
        let up_sel = app.source_selection == SourceSelection::UniProt;
        if tree_row(ui, up_sel, "UniProt", None, pal.accent, &pal) {
            next_source = Some(SourceSelection::UniProt);
        }
    });

    if let Some(sel) = next_source {
        app.source_selection = sel.clone();
        app.show_doc_table = true;
        // If switching into a folder with documents, focus the first document in that folder.
        if let SourceSelection::Folder(ref folder) = sel {
            if let Some(first_idx) = app
                .docs
                .iter()
                .position(|d| !d.deleted && !d.shared && &d.folder == folder)
            {
                app.select_document(first_idx);
            }
        } else if let SourceSelection::Ncbi(db) = sel {
            if app.ncbi_search_results.is_empty() {
                app.trigger_ncbi_search(ctx, db);
            }
        } else if sel == SourceSelection::UniProt && app.uniprot_search_results.is_empty() {
            app.trigger_uniprot_search(ctx);
        }
    }
}

fn tree_row(
    ui: &mut Ui,
    selected: bool,
    label: &str,
    count: Option<usize>,
    text_color: Color32,
    pal: &Palette,
) -> bool {
    let mut clicked = false;
    ui.horizontal(|ui| {
        let resp = ui.selectable_label(
            selected,
            RichText::new(label)
                .size(12.0)
                .color(if selected { pal.accent } else { text_color }),
        );
        if resp.clicked() {
            clicked = true;
        }
        if let Some(n) = count {
            ui.label(
                RichText::new(n.to_string())
                    .size(10.5)
                    .color(pal.text_muted),
            );
        }
    });
    clicked
}
