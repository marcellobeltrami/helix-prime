//! Right-hand Geneious Prime-inspired Inspector Drawer: Display controls, Annotations &
//! UniProtKB lookup/3D launcher, Motif search, and Pairwise Identity Statistics.

use crate::app::{HelixApp, InspectorTab};
use crate::theme::{ColorScheme, HighlightMode, Palette, ThemeMode};
use egui::{Button, Color32, Context, Frame, RichText, Rounding, ScrollArea, Stroke, Ui};
use helix_core::{ConsensusConfig, UniProtLink};

/// Renders the right-hand inspector sidebar with 4 specialized tabs (General/Display,
/// Annotations & UniProt, Motif Search, and Statistics).
pub fn show_inspector_panel(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    ui.horizontal(|ui| {
        ui.selectable_value(&mut app.inspector_tab, InspectorTab::Display, "General");
        ui.selectable_value(
            &mut app.inspector_tab,
            InspectorTab::Annotations,
            "Annotations",
        );
        ui.selectable_value(&mut app.inspector_tab, InspectorTab::FindMotif, "Motif");
        ui.selectable_value(&mut app.inspector_tab, InspectorTab::Statistics, "Stats");
    });
    ui.separator();

    ScrollArea::vertical().show(ui, |ui| match app.inspector_tab {
        InspectorTab::Display => show_display_tab(app, ctx, ui, &pal),
        InspectorTab::Annotations => show_annotations_uniprot_tab(app, ctx, ui, &pal),
        InspectorTab::FindMotif => show_motif_tab(app, ui, &pal),
        InspectorTab::Statistics => show_statistics_tab(app, ui, &pal),
    });
}

fn show_display_tab(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    // Zoom row matching the top of Geneious Prime's General sidebar (`-`, `+`, `100 %`).
    ui.horizontal(|ui| {
        if ui.small_button("-").clicked() {
            app.display.col_width = (app.display.col_width - 2.0).max(4.0);
        }
        let zoom_pct = ((app.display.col_width / 14.0) * 100.0).round();
        ui.label(RichText::new(format!("{zoom_pct:.0}%")).monospace().strong());
        if ui.small_button("+").clicked() {
            app.display.col_width = (app.display.col_width + 2.0).min(26.0);
        }
        ui.separator();
        let mut dark_mode = app.theme == ThemeMode::Dark;
        if ui.checkbox(&mut dark_mode, "Dark Mode").changed() {
            app.set_theme(
                ctx,
                if dark_mode {
                    ThemeMode::Dark
                } else {
                    ThemeMode::Light
                },
            );
        }
    });
    ui.separator();

    section_heading(ui, "COLORS & HIGHLIGHTING", pal);

    egui::ComboBox::from_label("Colors")
        .selected_text(app.display.color_scheme.label())
        .show_ui(ui, |ui| {
            for scheme in ColorScheme::ALL {
                ui.selectable_value(&mut app.display.color_scheme, scheme, scheme.label());
            }
        });

    egui::ComboBox::from_label("Highlighting")
        .selected_text(app.display.highlight_mode.label())
        .show_ui(ui, |ui| {
            for mode in HighlightMode::ALL {
                ui.selectable_value(&mut app.display.highlight_mode, mode, mode.label());
            }
        });

    ui.checkbox(
        &mut app.display.dots_for_matches,
        "Dots (.) for consensus matches",
    );

    ui.add_space(8.0);
    section_heading(ui, "OVERLAYS", pal);

    ui.checkbox(&mut app.display.show_identity_graph, "Graphs (Identity / GC)");
    ui.checkbox(&mut app.display.show_annotations, "Annotations");
    ui.checkbox(&mut app.display.show_overview, "Overview minimap");
    ui.checkbox(&mut app.display.show_consensus, "Consensus track");

    let is_nuc = app
        .current_alignment()
        .is_some_and(|a| a.alphabet.is_nucleotide());
    if is_nuc {
        ui.checkbox(&mut app.display.show_complement, "Complement strand (3'->5')");
        ui.checkbox(&mut app.display.show_translation, "Translation");
        if app.display.show_translation {
            ui.indent("trans_frame_opts", |ui| {
                ui.horizontal(|ui| {
                    ui.label("Frame:");
                    ui.radio_value(&mut app.display.translation_frame, 0, "+1");
                    ui.radio_value(&mut app.display.translation_frame, 1, "+2");
                    ui.radio_value(&mut app.display.translation_frame, 2, "+3");
                });
            });
        }
        ui.checkbox(&mut app.display.show_restriction_sites, "Restriction Sites");
    }

    ui.add_space(8.0);
    section_heading(ui, "VIEW LAYOUT (PHAGE / LINEAR / WRAP)", pal);

    ui.checkbox(&mut app.display.show_circular_overview, "Circular Overview (Phage / Plasmid)");
    ui.checkbox(&mut app.display.show_linear_view, "Linear View");
    ui.checkbox(&mut app.display.wrap_sequence, "Wrap Sequence");

    if app.display.show_circular_overview {
        ui.indent("circular_opts", |ui| {
            ui.checkbox(&mut app.circular.show_labels, "Outer feature labels");
            ui.checkbox(&mut app.circular.show_sequence_rings, "Aligned strain rings");
            ui.checkbox(&mut app.circular.show_identity_ring, "Pairwise identity ring");
            ui.checkbox(&mut app.circular.show_gc_ring, "GC content ring");
            ui.add(egui::Slider::new(&mut app.circular.zoom, 0.55..=2.4).text("Circle zoom"));
        });
    }

    ui.add_space(8.0);
    section_heading(ui, "ZOOM & GEOMETRY", pal);

    ui.add(
        egui::Slider::new(&mut app.display.col_width, 4.0..=26.0)
            .text("Column width (px)"),
    );
    ui.add(
        egui::Slider::new(&mut app.display.row_height, 16.0..=30.0)
            .text("Row height (px)"),
    );

    ui.add_space(8.0);
    section_heading(ui, "CONSENSUS CALLING", pal);

    if let Some(doc) = app.current_doc_mut() {
        let mut cfg: ConsensusConfig = doc.alignment.consensus_config();
        let mut pct = cfg.threshold * 100.0;
        let changed_t = ui
            .add(egui::Slider::new(&mut pct, 25.0..=100.0).text("Threshold %"))
            .changed();
        let changed_g = ui
            .checkbox(&mut cfg.ignore_gaps, "Ignore gaps in consensus")
            .changed();
        if changed_t || changed_g {
            cfg.threshold = pct / 100.0;
            doc.alignment.set_consensus_config(cfg);
        }
    }
}

fn show_annotations_uniprot_tab(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    section_heading(ui, "ALIGNMENT ANNOTATIONS", pal);

    let mut select_ann: Option<(u64, usize, usize)> = None;
    let mut launch_3d: Option<(String, Option<u64>)> = None;
    let mut fetch_uniprot_feats: Option<(String, Option<u64>)> = None;
    let mut delete_ann_id: Option<u64> = None;

    if let Some(aln) = app.current_alignment() {
        if aln.annotations().is_empty() {
            ui.label(
                RichText::new("No annotations yet. Select columns and click '+ Annotation' or import GFF3/UniProt features.")
                    .small()
                    .color(pal.text_muted),
            );
        } else {
            for ann in aln.annotations() {
                let is_sel = app.selected_annotation_id == Some(ann.id);
                let [r, g, b] = ann.rgb();
                let resp = Frame::none()
                    .fill(if is_sel {
                        pal.accent_soft
                    } else {
                        pal.canvas_bg
                    })
                    .stroke(Stroke::new(
                        1.0,
                        if is_sel { pal.accent } else { pal.border },
                    ))
                    .rounding(Rounding::same(4.0))
                    .inner_margin(egui::Margin::symmetric(7.0, 5.0))
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            let (dot_rect, _) =
                                ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                            ui.painter()
                                .rect_filled(dot_rect, 2.0, Color32::from_rgb(r, g, b));
                            ui.label(RichText::new(&ann.name).strong().size(12.0));
                            ui.with_layout(
                                egui::Layout::right_to_left(egui::Align::Center),
                                |ui| {
                                    if let Some(up) = &ann.uniprot {
                                        Frame::none()
                                            .fill(Color32::from_rgb(16, 138, 134))
                                            .rounding(Rounding::same(3.0))
                                            .inner_margin(egui::Margin::symmetric(4.0, 1.0))
                                            .show(ui, |ui| {
                                                ui.label(
                                                    RichText::new(&up.accession)
                                                        .size(9.5)
                                                        .strong()
                                                        .color(Color32::WHITE),
                                                );
                                            });
                                    }
                                },
                            );
                        });
                        ui.label(
                            RichText::new(format!(
                                "{}  |  {}..{} ({} cols)",
                                ann.kind.label(),
                                ann.start + 1,
                                ann.end,
                                ann.len()
                            ))
                            .size(10.5)
                            .color(pal.text_muted),
                        );
                    })
                    .response
                    .interact(egui::Sense::click());

                if resp.clicked() {
                    select_ann = Some((ann.id, ann.start, ann.end));
                }
                ui.add_space(3.0);
            }
        }
    }

    if let Some((id, s, e)) = select_ann {
        app.selected_annotation_id = Some(id);
        app.selected_cols = Some((s, e));
        app.scroll_to_column = Some(s);
        app.rebuild_3d_alignment_map();
    }

    // Selected annotation detail editor & UniProt / AlphaFold 3D action card.
    if let Some(sel_id) = app.selected_annotation_id {
        ui.add_space(8.0);
        section_heading(ui, "SELECTED ANNOTATION & 3D LINK", pal);

        if let Some(doc) = app.current_doc_mut() {
            if let Some(ann) = doc.alignment.annotation_by_id_mut(sel_id) {
                Frame::none()
                    .fill(pal.canvas_bg)
                    .stroke(Stroke::new(1.0, pal.border))
                    .rounding(Rounding::same(6.0))
                    .inner_margin(egui::Margin::same(8.0))
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(RichText::new("Name:").small().color(pal.text_muted));
                            ui.text_edit_singleline(&mut ann.name);
                        });

                        ui.label(RichText::new("Notes:").small().color(pal.text_muted));
                        ui.add(
                            egui::TextEdit::multiline(&mut ann.notes)
                                .desired_rows(2)
                                .desired_width(f32::INFINITY),
                        );

                        ui.add_space(4.0);
                        let mut acc_buf = ann
                            .uniprot
                            .as_ref()
                            .map(|u| u.accession.clone())
                            .unwrap_or_default();
                        ui.horizontal(|ui| {
                            ui.label(RichText::new("UniProt ID:").small().strong());
                            let resp = ui.add(
                                egui::TextEdit::singleline(&mut acc_buf)
                                    .desired_width(110.0)
                                    .hint_text("e.g. P68871"),
                            );
                            if resp.changed() {
                                let trimmed = acc_buf.trim().to_ascii_uppercase();
                                if trimmed.is_empty() {
                                    ann.uniprot = None;
                                } else if let Some(link) = ann.uniprot.as_mut() {
                                    link.accession = trimmed;
                                } else {
                                    ann.uniprot = Some(UniProtLink::new(trimmed));
                                }
                            }
                        });

                        if let Some(up) = &ann.uniprot {
                            if let Some(pname) = &up.protein_name {
                                ui.label(
                                    RichText::new(pname).small().color(pal.text_primary),
                                );
                            }
                            if let Some(org) = &up.organism {
                                ui.label(
                                    RichText::new(org).small().italics().color(pal.text_muted),
                                );
                            }

                            ui.add_space(6.0);
                            let open_btn = Button::new(
                                RichText::new(format!("View 3D in AlphaFold ({})", up.accession))
                                    .strong()
                                    .color(Color32::WHITE),
                            )
                            .fill(Color32::from_rgb(16, 138, 134));
                            if ui.add_sized([ui.available_width(), 28.0], open_btn).clicked() {
                                launch_3d = Some((up.accession.clone(), Some(ann.id)));
                            }

                            ui.add_space(4.0);
                            if ui
                                .button("Import UniProt Domains & Sites")
                                .on_hover_text("Fetch curated UniProtKB functional annotations and overlay onto this alignment")
                                .clicked()
                            {
                                fetch_uniprot_feats = Some((up.accession.clone(), Some(ann.id)));
                            }
                        }

                        ui.add_space(4.0);
                        if ui
                            .button(RichText::new("Delete Annotation").small().color(pal.identity_low))
                            .clicked()
                        {
                            delete_ann_id = Some(ann.id);
                        }
                    });
            }
        }
    }

    if let Some(del_id) = delete_ann_id {
        if let Some(doc) = app.current_doc_mut() {
            doc.alignment.remove_annotation(del_id);
        }
        app.selected_annotation_id = None;
    }

    if let Some((acc, ann_id)) = launch_3d {
        app.open_3d_for_accession(ctx, &acc, ann_id);
    }

    if let Some((acc, ann_id)) = fetch_uniprot_feats {
        app.fetch_and_import_uniprot_annotations(ctx, &acc, ann_id);
    }

    // Live UniProtKB Search Engine right inside the Inspector!
    ui.add_space(12.0);
    section_heading(ui, "SEARCH UNIPROTKB DATABASE", pal);
    ui.label(
        RichText::new("Search by gene, protein name, or organism to link an annotation and open its AlphaFold 3D model:")
            .small()
            .color(pal.text_muted),
    );
    ui.add_space(3.0);

    let mut do_search = false;
    ui.horizontal(|ui| {
        let te = ui.add(
            egui::TextEdit::singleline(&mut app.uniprot_search_query)
                .desired_width(ui.available_width() - 68.0)
                .hint_text("e.g. BRCA1 human, insulin, KRAS"),
        );
        if te.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            do_search = true;
        }
        if ui.button("Search").clicked() {
            do_search = true;
        }
    });
    if do_search {
        app.trigger_uniprot_search(ctx);
    }

    if app.uniprot_searching {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label(RichText::new("Searching UniProtKB REST API...").small());
        });
    }
    if let Some(err) = &app.uniprot_search_error {
        ui.label(RichText::new(err).small().color(pal.identity_low));
    }

    let mut attach_hit: Option<(String, String, String, String)> = None;
    let mut view_hit_3d: Option<String> = None;

    for hit in &app.uniprot_search_results {
        Frame::none()
            .fill(pal.canvas_bg)
            .stroke(Stroke::new(1.0, pal.border))
            .rounding(Rounding::same(4.0))
            .inner_margin(egui::Margin::same(6.0))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new(&hit.accession)
                            .strong()
                            .monospace()
                            .color(pal.accent),
                    );
                    if hit.reviewed {
                        ui.label(
                            RichText::new("Swiss-Prot")
                                .size(9.5)
                                .color(pal.identity_high),
                        );
                    }
                    ui.label(
                        RichText::new(format!("{} aa", hit.length))
                            .size(10.0)
                            .color(pal.text_muted),
                    );
                });
                ui.label(RichText::new(&hit.protein_name).size(11.5));
                ui.label(
                    RichText::new(&hit.organism)
                        .size(10.5)
                        .italics()
                        .color(pal.text_muted),
                );
                ui.horizontal(|ui| {
                    if ui.small_button("Link to Selected Annotation").clicked() {
                        attach_hit = Some((
                            hit.accession.clone(),
                            hit.entry_name.clone(),
                            hit.protein_name.clone(),
                            hit.organism.clone(),
                        ));
                    }
                    if ui.small_button("View 3D").clicked() {
                        view_hit_3d = Some(hit.accession.clone());
                    }
                });
            });
        ui.add_space(3.0);
    }

    if let Some((acc, ename, pname, org)) = attach_hit {
        if let Some(sel_id) = app.selected_annotation_id {
            if let Some(doc) = app.current_doc_mut() {
                if let Some(ann) = doc.alignment.annotation_by_id_mut(sel_id) {
                    ann.uniprot =
                        Some(UniProtLink::new(&acc).with_metadata(ename, pname, org));
                }
            }
            app.status_message = format!("Linked annotation to UniProt {acc}.");
        } else {
            app.status_message = format!("Select an annotation first to link UniProt {acc}.");
        }
    }

    if let Some(acc) = view_hit_3d {
        app.open_3d_for_accession(ctx, &acc, app.selected_annotation_id);
    }
}

fn show_motif_tab(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    section_heading(ui, "FIND MOTIF / SUBSEQUENCE", pal);
    ui.label(
        RichText::new("Searches across alignment gaps. Supports IUPAC ambiguity codes (R, Y, N, X).")
            .small()
            .color(pal.text_muted),
    );
    ui.add_space(4.0);

    let mut run = false;
    ui.horizontal(|ui| {
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.motif_query)
                .desired_width(ui.available_width() - 60.0)
                .hint_text("e.g. LSELHCD or ATGGTG"),
        );
        if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            run = true;
        }
        if ui.button("Find").clicked() {
            run = true;
        }
    });
    if ui
        .checkbox(&mut app.motif_allow_iupac, "Allow IUPAC / X wildcards")
        .changed()
    {
        run = true;
    }

    if run {
        app.run_motif_search();
    }

    ui.add_space(8.0);
    if app.motif_hits.is_empty() {
        if !app.motif_query.trim().is_empty() {
            ui.label(RichText::new("No matches found.").small().color(pal.text_muted));
        }
    } else {
        ui.horizontal(|ui| {
            ui.label(
                RichText::new(format!(
                    "Match {} of {}",
                    app.active_motif_hit + 1,
                    app.motif_hits.len()
                ))
                .strong(),
            );
            if ui.small_button("Prev").clicked() {
                if app.active_motif_hit == 0 {
                    app.active_motif_hit = app.motif_hits.len() - 1;
                } else {
                    app.active_motif_hit -= 1;
                }
                let hit = app.motif_hits[app.active_motif_hit];
                app.scroll_to_column = Some(hit.start_col);
                app.selected_cols = Some((hit.start_col, hit.end_col));
            }
            if ui.small_button("Next").clicked() {
                app.active_motif_hit = (app.active_motif_hit + 1) % app.motif_hits.len();
                let hit = app.motif_hits[app.active_motif_hit];
                app.scroll_to_column = Some(hit.start_col);
                app.selected_cols = Some((hit.start_col, hit.end_col));
            }
        });

        ui.add_space(6.0);
        let seq_names: Vec<String> = app
            .current_alignment()
            .map(|a| a.sequences().iter().map(|s| s.name.clone()).collect())
            .unwrap_or_default();

        for (idx, hit) in app.motif_hits.iter().enumerate().take(60) {
            let sname = seq_names.get(hit.row).map(String::as_str).unwrap_or("Seq");
            let selected = idx == app.active_motif_hit;
            if ui
                .selectable_label(
                    selected,
                    format!("{sname}  ({}..{})", hit.start_col + 1, hit.end_col),
                )
                .clicked()
            {
                app.active_motif_hit = idx;
                app.scroll_to_column = Some(hit.start_col);
                app.selected_cols = Some((hit.start_col, hit.end_col));
            }
        }
    }
}

fn show_statistics_tab(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    section_heading(ui, "ALIGNMENT SUMMARY", pal);

    // Lazily populate pairwise identity matrix cache.
    if let Some(doc) = app.current_doc_mut() {
        if doc.pairwise_cache.is_none() {
            doc.pairwise_cache = Some(doc.alignment.pairwise_identity_matrix());
        }
    }

    let Some(doc) = app.docs.get(app.active_doc) else {
        return;
    };
    let aln = &doc.alignment;

    Frame::none()
        .fill(pal.canvas_bg)
        .stroke(Stroke::new(1.0, pal.border))
        .rounding(Rounding::same(6.0))
        .inner_margin(egui::Margin::same(8.0))
        .show(ui, |ui| {
            egui::Grid::new("stats_summary_grid")
                .num_columns(2)
                .spacing([12.0, 4.0])
                .show(ui, |ui| {
                    ui.label(RichText::new("Alphabet").color(pal.text_muted));
                    ui.label(RichText::new(aln.alphabet.label()).strong());
                    ui.end_row();

                    ui.label(RichText::new("Sequences").color(pal.text_muted));
                    ui.label(RichText::new(aln.row_count().to_string()).strong());
                    ui.end_row();

                    ui.label(RichText::new("Alignment Width").color(pal.text_muted));
                    ui.label(RichText::new(aln.width().to_string()).strong());
                    ui.end_row();

                    ui.label(RichText::new("Mean Pairwise ID").color(pal.text_muted));
                    ui.label(
                        RichText::new(format!("{:.1}%", aln.mean_pairwise_identity_pct()))
                            .strong()
                            .color(pal.identity_high),
                    );
                    ui.end_row();

                    ui.label(RichText::new("Identical Sites").color(pal.text_muted));
                    ui.label(
                        RichText::new(format!("{:.1}%", aln.identical_sites_pct())).strong(),
                    );
                    ui.end_row();

                    if aln.alphabet.is_nucleotide() {
                        if let Some(gc) = aln.sequences().first().and_then(|s| s.gc_content_pct()) {
                            ui.label(RichText::new("Ref GC Content").color(pal.text_muted));
                            ui.label(RichText::new(format!("{gc:.1}%")).strong());
                            ui.end_row();
                        }
                    }
                });
        });

    ui.add_space(10.0);
    section_heading(ui, "PAIRWISE IDENTITY MATRIX (%)", pal);

    if let Some(matrix) = &doc.pairwise_cache {
        let n = aln.row_count().min(12);
        ScrollArea::horizontal().show(ui, |ui| {
            egui::Grid::new("pairwise_matrix_grid")
                .spacing([3.0, 3.0])
                .show(ui, |ui| {
                    ui.label("");
                    for j in 0..n {
                        let short: String = aln.sequences()[j].name.chars().take(6).collect();
                        ui.label(RichText::new(short).size(9.5).color(pal.text_muted));
                    }
                    ui.end_row();

                    for i in 0..n {
                        let short: String = aln.sequences()[i].name.chars().take(10).collect();
                        ui.label(RichText::new(short).size(10.0).strong());
                        for j in 0..n {
                            let val = matrix[i * aln.row_count() + j];
                            let fill = if i == j {
                                pal.header_bg
                            } else if val >= 90.0 {
                                Color32::from_rgba_premultiplied(46, 160, 92, 65)
                            } else if val >= 70.0 {
                                Color32::from_rgba_premultiplied(214, 174, 42, 65)
                            } else {
                                Color32::from_rgba_premultiplied(212, 76, 76, 65)
                            };
                            Frame::none()
                                .fill(fill)
                                .rounding(Rounding::same(2.0))
                                .inner_margin(egui::Margin::symmetric(4.0, 2.0))
                                .show(ui, |ui| {
                                    ui.label(
                                        RichText::new(format!("{val:.0}"))
                                            .size(10.0)
                                            .monospace(),
                                    );
                                });
                        }
                        ui.end_row();
                    }
                });
        });
    }
}

fn section_heading(ui: &mut Ui, title: &str, pal: &Palette) {
    ui.label(
        RichText::new(title)
            .strong()
            .small()
            .color(pal.text_muted),
    );
    ui.add_space(3.0);
}
