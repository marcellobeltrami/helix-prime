//! Top Geneious Prime-inspired toolbar with Back/Forward history, Add, Export, BLAST,
//! Workflows, Align/Assemble, Tree, Primers, Cloning, Phage Circular, 3D Viewer, Dark Mode,
//! and right-aligned Search Everywhere.

use crate::app::{HelixApp, SourceSelection};
use crate::theme::{Palette, ThemeMode};
use egui::{Button, Color32, Context, RichText, Ui};

/// Renders the main application ribbon toolbar matching the Geneious Prime interface.
pub fn show_toolbar(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    ui.horizontal_centered(|ui| {
        // Back / Forward document history navigation.
        let can_back = app.history_pos > 0;
        let can_fwd = app.history_pos + 1 < app.history.len();
        if ui
            .add_enabled(can_back, Button::new("←"))
            .on_hover_text("Back to previous document")
            .clicked()
        {
            app.navigate_back();
        }
        if ui
            .add_enabled(can_fwd, Button::new("→"))
            .on_hover_text("Forward to next document")
            .clicked()
        {
            app.navigate_forward();
        }
        ui.separator();

        // Add / Export dropdown-style quick actions.
        ui.menu_button(RichText::new("+ Add").strong(), |ui| {
            if ui.button("Import Sequence / Alignment (FASTA, Clustal, GenBank)...").clicked() {
                app.import_alignment_dialog();
                ui.close_menu();
            }
            if ui.button("Import GFF3 Annotations onto Active Document...").clicked() {
                app.import_gff3_dialog();
                ui.close_menu();
            }
            if ui.button("Load Local 3D Structure (PDB / mmCIF)...").clicked() {
                app.load_local_structure_dialog();
                ui.close_menu();
            }
            ui.separator();
            if ui.button("New Annotation on Current Selection...").clicked() {
                if let Some((s, e)) = app.selected_cols {
                    app.new_ann_draft.start_1based = s + 1;
                    app.new_ann_draft.end_1based = e;
                }
                app.new_ann_draft.open = true;
                ui.close_menu();
            }
        });

        if ui
            .button("Export")
            .on_hover_text("Export active alignment (FASTA/Clustal) or tree (Newick)")
            .clicked()
        {
            app.export_alignment_dialog();
        }

        ui.separator();

        // Core Molecular Biology & Analysis Tools (matching the Geneious Prime toolbar).
        if ui
            .button(RichText::new("BLAST").strong())
            .on_hover_text("Run local workspace BLAST or search NCBI/UniProt")
            .clicked()
        {
            if app.blast_modal.query.trim().is_empty() {
                if let Some(aln) = app.current_alignment() {
                    if let Some(seq) = aln.sequences().first() {
                        let (s, e) = app
                            .selected_cols
                            .unwrap_or((0, aln.width().min(60)));
                        let frag = &seq.residues[s.min(seq.len())..e.min(seq.len())];
                        app.blast_modal.query = String::from_utf8_lossy(frag).into_owned();
                    }
                }
            }
            app.blast_modal.open = true;
        }

        if ui
            .button("Workflows")
            .on_hover_text("Run automated multi-step bioinformatics pipelines")
            .clicked()
        {
            app.workflow_modal.open = true;
        }

        if ui
            .button("Align/Assemble")
            .on_hover_text("Multiple sequence alignment, de novo contig assembly, or read mapping")
            .clicked()
        {
            app.assemble_modal.open = true;
        }

        if ui
            .button("Tree")
            .on_hover_text("Build UPGMA or Neighbor-Joining phylogenetic tree")
            .clicked()
        {
            app.tree_modal.open = true;
        }

        if ui
            .button("Primers")
            .on_hover_text("Design PCR primer pairs and test primer binding sites")
            .clicked()
        {
            app.primers_modal.open = true;
        }

        if ui
            .button(RichText::new("Cloning").strong())
            .on_hover_text("In-silico Cloning: Restriction Digest, TOPO-TA, Gateway, and BioBrick")
            .clicked()
        {
            app.cloning_modal.open = true;
        }

        ui.separator();

        // Phage Circular View quick toggle.
        let is_circular = app.display.show_circular_overview;
        let phage_btn = if is_circular && !app.display.show_linear_view {
            Button::new(RichText::new("Phage Circle").strong().color(Color32::WHITE))
                .fill(Color32::from_rgb(16, 138, 134))
        } else {
            Button::new(RichText::new("Phage Circle"))
        };
        if ui
            .add(phage_btn)
            .on_hover_text("Toggle full Circular Phage / Plasmid genome ring view")
            .clicked()
        {
            if app.display.show_circular_overview && !app.display.show_linear_view {
                app.display.show_circular_overview = false;
                app.display.show_linear_view = true;
            } else {
                app.display.show_circular_overview = true;
                app.display.show_linear_view = false;
            }
        }

        // Separate 3D Protein Window CTA.
        let active_acc = find_active_uniprot_accession(app);
        let btn_label = match &app.viewer3d.structure {
            Some(st) if app.viewer3d.open => format!("3D ({})", st.entry_id),
            _ => match &active_acc {
                Some(acc) => format!("3D ({acc})"),
                None => String::from("3D Viewer"),
            },
        };
        let btn_3d =
            Button::new(RichText::new(btn_label).strong().color(Color32::WHITE)).fill(pal.accent);
        if ui
            .add(btn_3d)
            .on_hover_text("Open interactive AlphaFold / UniProt 3D protein viewer in a separate window")
            .clicked()
        {
            if app.viewer3d.structure.is_some() {
                app.viewer3d.open = true;
            } else if let Some(acc) = active_acc {
                app.open_3d_for_accession(ctx, &acc, app.selected_annotation_id);
            } else {
                app.open_offline_3d_preview();
            }
        }

        // Right side: Search Everywhere + Dark Mode toggle + Panel toggles.
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let is_dark = app.theme == ThemeMode::Dark;
            let theme_btn = if is_dark {
                Button::new(RichText::new("Light").strong())
            } else {
                Button::new(RichText::new("Dark").strong().color(Color32::WHITE))
                    .fill(Color32::from_rgb(42, 52, 68))
            };
            if ui
                .add(theme_btn)
                .on_hover_text("Toggle Dark Mode / Light Mode")
                .clicked()
            {
                app.set_theme(
                    ctx,
                    if is_dark {
                        ThemeMode::Light
                    } else {
                        ThemeMode::Dark
                    },
                );
            }

            if ui
                .selectable_label(app.show_doc_table, "Table")
                .on_hover_text("Toggle upper Document Table browser")
                .clicked()
            {
                app.show_doc_table = !app.show_doc_table;
            }

            // Global Search Everywhere input (matching top-right of screenshot).
            let search_resp = ui.add(
                egui::TextEdit::singleline(&mut app.global_search)
                    .desired_width(165.0)
                    .hint_text("Search Everywhere..."),
            );
            if search_resp.changed() && !app.global_search.trim().is_empty() {
                app.show_doc_table = true;
                if matches!(
                    app.source_selection,
                    SourceSelection::Ncbi(_) | SourceSelection::UniProt | SourceSelection::Operations
                ) {
                    app.source_selection = SourceSelection::AllLocal;
                }
            }
        });
    });
}

fn find_active_uniprot_accession(app: &HelixApp) -> Option<String> {
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
    if let Some(acc) = aln
        .annotations()
        .iter()
        .find_map(|a| a.uniprot.as_ref().map(|u| u.accession.clone()))
    {
        return Some(acc);
    }
    aln.sequences()
        .iter()
        .find_map(|s| s.uniprot_accession.clone())
}
