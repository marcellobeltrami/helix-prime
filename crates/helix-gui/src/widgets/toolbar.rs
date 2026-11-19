//! Top Geneious Prime-inspired ribbon toolbar with file import/export, progressive alignment,
//! annotation creation, zoom controls, theme toggle, and separate 3D Protein Window launcher.

use crate::app::{HelixApp, InspectorTab, ViewMode};
use crate::theme::{Palette, ThemeMode};
use egui::{Button, Color32, Context, RichText, Ui};

/// Renders the main application toolbar.
pub fn show_toolbar(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    ui.horizontal_centered(|ui| {
        // Brand badge.
        ui.label(
            RichText::new("HelixPrime")
                .strong()
                .size(16.0)
                .color(pal.accent),
        );
        ui.separator();

        // Sidebar visibility toggles.
        if ui
            .selectable_label(app.show_sources_panel, "Sources")
            .on_hover_text("Toggle left Sources & Documents browser")
            .clicked()
        {
            app.show_sources_panel = !app.show_sources_panel;
        }
        if ui
            .selectable_label(app.show_inspector_panel, "Inspector")
            .on_hover_text("Toggle right Display, Annotations & UniProt inspector")
            .clicked()
        {
            app.show_inspector_panel = !app.show_inspector_panel;
        }
        ui.separator();

        // File I/O actions.
        if ui
            .button("Import...")
            .on_hover_text("Open FASTA (.fa/.fasta), Clustal (.aln), or GenBank (.gb/.gbk)")
            .clicked()
        {
            app.import_alignment_dialog();
        }
        if ui
            .button("+ GFF3")
            .on_hover_text("Overlay GFF3 annotations onto the active alignment")
            .clicked()
        {
            app.import_gff3_dialog();
        }
        if ui
            .button("Export...")
            .on_hover_text("Save active alignment as FASTA or Clustal (.aln)")
            .clicked()
        {
            app.export_alignment_dialog();
        }
        ui.separator();

        // Alignment & Annotation actions.
        if ui
            .button("Align (NW)")
            .on_hover_text("Run progressive Needleman-Wunsch multiple sequence alignment")
            .clicked()
        {
            let res = app
                .current_doc_mut()
                .map(|doc| {
                    let r = doc.alignment.realign_progressive();
                    if r.is_ok() {
                        doc.invalidate_caches();
                    }
                    r
                });
            match res {
                Some(Ok(())) => {
                    app.rebuild_3d_alignment_map();
                    app.status_message = String::from("Re-aligned sequences using progressive Needleman-Wunsch.");
                }
                Some(Err(err)) => {
                    app.status_message = format!("Alignment failed: {err}");
                }
                None => {}
            }
        }

        if ui
            .button("+ Annotation")
            .on_hover_text("Create a new annotation from the selected column range")
            .clicked()
        {
            if let Some((s, e)) = app.selected_cols {
                app.new_ann_draft.start_1based = s + 1;
                app.new_ann_draft.end_1based = e;
            }
            app.new_ann_draft.open = true;
        }

        if ui
            .button("UniProt Search")
            .on_hover_text("Open UniProtKB lookup drawer to search and link protein annotations")
            .clicked()
        {
            app.show_inspector_panel = true;
            app.inspector_tab = InspectorTab::Annotations;
        }

        ui.separator();

        // Prominent Phage / Circular Genome View toggle button.
        let is_circular = app.view_mode == ViewMode::CircularPhage;
        let phage_btn = if is_circular {
            Button::new(
                RichText::new("Phage Circular: ON")
                    .strong()
                    .color(Color32::WHITE),
            )
            .fill(Color32::from_rgb(16, 138, 134))
        } else {
            Button::new(RichText::new("Phage Circular").strong())
        };
        if ui
            .add(phage_btn)
            .on_hover_text("Toggle circular phage / plasmid genome ring visualization")
            .clicked()
        {
            app.view_mode = if is_circular {
                ViewMode::Linear
            } else {
                ViewMode::CircularPhage
            };
            app.status_message = match app.view_mode {
                ViewMode::Linear => String::from("Switched to Linear Alignment View."),
                ViewMode::CircularPhage => {
                    String::from("Switched to Circular Phage Genome Visualization (drag to rotate, scroll to zoom).")
                }
            };
        }

        // Highlighted primary CTA: Launch/Focus the separate native OS 3D Protein Viewer window.
        let active_acc = find_active_uniprot_accession(app);
        let btn_label = match &app.viewer3d.structure {
            Some(st) if app.viewer3d.open => format!("3D Window ({})", st.entry_id),
            _ => match &active_acc {
                Some(acc) => format!("View 3D ({acc})"),
                None => String::from("3D Protein Window"),
            },
        };

        let btn_3d = Button::new(RichText::new(btn_label).strong().color(Color32::WHITE))
            .fill(pal.accent);
        if ui
            .add(btn_3d)
            .on_hover_text("Open interactive 3D protein viewer in a separate OS window via AlphaFold DB & UniProt")
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

        // Right-aligned zoom and Dark Mode controls.
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let is_dark = app.theme == ThemeMode::Dark;
            let theme_btn = if is_dark {
                Button::new(RichText::new("Light Mode").strong())
            } else {
                Button::new(RichText::new("Dark Mode").strong().color(Color32::WHITE))
                    .fill(Color32::from_rgb(42, 52, 68))
            };
            if ui
                .add(theme_btn)
                .on_hover_text("Switch between modern Light Mode and Dark Mode")
                .clicked()
            {
                let next = if is_dark {
                    ThemeMode::Light
                } else {
                    ThemeMode::Dark
                };
                app.set_theme(ctx, next);
            }

            ui.separator();

            if ui.button("+").on_hover_text("Zoom in").clicked() {
                match app.view_mode {
                    ViewMode::Linear => {
                        app.display.col_width = (app.display.col_width + 2.0).min(26.0);
                    }
                    ViewMode::CircularPhage => {
                        app.circular.zoom = (app.circular.zoom * 1.15).min(2.4);
                    }
                }
            }
            let zoom_txt = match app.view_mode {
                ViewMode::Linear => format!("{:.0}px", app.display.col_width),
                ViewMode::CircularPhage => format!("{:.0}%", app.circular.zoom * 100.0),
            };
            ui.label(
                RichText::new(zoom_txt)
                    .small()
                    .monospace()
                    .color(pal.text_muted),
            );
            if ui.button("-").on_hover_text("Zoom out").clicked() {
                match app.view_mode {
                    ViewMode::Linear => {
                        app.display.col_width = (app.display.col_width - 2.0).max(4.0);
                    }
                    ViewMode::CircularPhage => {
                        app.circular.zoom = (app.circular.zoom * 0.87).max(0.55);
                    }
                }
            }
        });
    });
}

/// Resolves the best UniProt accession for the current selection (selected annotation -> first
/// UniProt-linked annotation -> first sequence UniProt accession).
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
