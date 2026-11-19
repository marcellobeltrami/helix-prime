//! Separate native OS window (`egui::ViewportId`) hosting the interactive 3D Protein Viewer
//! connected to UniProtKB, AlphaFold Protein Structure Database (EBI), and RCSB PDB.

use crate::app::HelixApp;
use crate::renderer3d::{ColorBy3D, RenderHighlights, RenderStyle};
use crate::theme::{plddt_color, Palette};
use egui::{
    Align2, Button, Color32, Context, FontId, Frame, Pos2, RichText, Rounding, ScrollArea, Sense,
    Stroke, Ui,
};

/// Renders the complete UI of the separate 3D Protein Viewer OS window.
pub fn show_protein_window(app: &mut HelixApp, vp_ctx: &Context) {
    let pal = Palette::for_mode(app.theme);

    // 1. Top controls bar inside the 3D window.
    egui::TopBottomPanel::top("protein3d_toolbar")
        .frame(
            Frame::none()
                .fill(pal.header_bg)
                .inner_margin(egui::Margin::symmetric(10.0, 7.0))
                .stroke(Stroke::new(1.0, pal.border)),
        )
        .show(vp_ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("UniProt:").strong());
                let acc_edit = ui.add(
                    egui::TextEdit::singleline(&mut app.viewer3d.accession_input)
                        .desired_width(86.0)
                        .hint_text("P68871"),
                );
                let mut trigger_fetch =
                    acc_edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));

                let fetch_btn = Button::new(
                    RichText::new("Fetch AlphaFold DB")
                        .strong()
                        .color(Color32::WHITE),
                )
                .fill(Color32::from_rgb(16, 138, 134));
                if ui.add(fetch_btn).clicked() {
                    trigger_fetch = true;
                }

                if trigger_fetch {
                    let acc = app.viewer3d.accession_input.clone();
                    app.open_3d_for_accession(vp_ctx, &acc, app.selected_annotation_id);
                }

                if ui
                    .button("Local PDB/mmCIF...")
                    .on_hover_text("Load a local .pdb or .cif protein structure file")
                    .clicked()
                {
                    app.load_local_structure_dialog();
                }

                ui.separator();

                egui::ComboBox::from_id_salt("style_3d_combo")
                    .selected_text(app.viewer3d.style.label())
                    .width(120.0)
                    .show_ui(ui, |ui| {
                        for st in RenderStyle::ALL {
                            ui.selectable_value(&mut app.viewer3d.style, st, st.label());
                        }
                    });

                egui::ComboBox::from_id_salt("color_3d_combo")
                    .selected_text(app.viewer3d.color_by.label())
                    .width(136.0)
                    .show_ui(ui, |ui| {
                        for cb in ColorBy3D::ALL {
                            ui.selectable_value(&mut app.viewer3d.color_by, cb, cb.label());
                        }
                    });

                ui.separator();

                ui.checkbox(&mut app.viewer3d.camera.auto_rotate, "Spin");
                if ui.button("Reset View").clicked() {
                    app.viewer3d.camera.reset();
                }
            });
        });

    // 2. Right metadata & UniProt/AlphaFold inspector sidebar inside the 3D window.
    egui::SidePanel::right("protein3d_info_panel")
        .default_width(250.0)
        .min_width(210.0)
        .frame(
            Frame::none()
                .fill(pal.panel_bg)
                .inner_margin(egui::Margin::same(10.0))
                .stroke(Stroke::new(1.0, pal.border)),
        )
        .show(vp_ctx, |ui| {
            show_3d_metadata_sidebar(app, vp_ctx, ui, &pal);
        });

    // 3. Center interactive 3D molecular viewport.
    egui::CentralPanel::default()
        .frame(Frame::none().fill(pal.viewport_3d_bg_top).inner_margin(egui::Margin::same(0.0)))
        .show(vp_ctx, |ui| {
            let rect = ui.available_rect_before_wrap();
            let response = ui.allocate_rect(rect, Sense::click_and_drag());
            let painter = ui.painter_at(rect);

            // Subtle vertical gradient laboratory backdrop.
            let mut bg_mesh = egui::Mesh::default();
            bg_mesh.vertices.push(egui::epaint::Vertex {
                pos: rect.left_top(),
                uv: egui::epaint::WHITE_UV,
                color: pal.viewport_3d_bg_top,
            });
            bg_mesh.vertices.push(egui::epaint::Vertex {
                pos: rect.right_top(),
                uv: egui::epaint::WHITE_UV,
                color: pal.viewport_3d_bg_top,
            });
            bg_mesh.vertices.push(egui::epaint::Vertex {
                pos: rect.right_bottom(),
                uv: egui::epaint::WHITE_UV,
                color: pal.viewport_3d_bg_bottom,
            });
            bg_mesh.vertices.push(egui::epaint::Vertex {
                pos: rect.left_bottom(),
                uv: egui::epaint::WHITE_UV,
                color: pal.viewport_3d_bg_bottom,
            });
            bg_mesh.indices.extend_from_slice(&[0, 1, 2, 0, 2, 3]);
            painter.add(egui::Shape::mesh(bg_mesh));

            // Orbit (left-drag), Pan (right-drag or Shift+drag), Zoom (scroll), Reset (double-click).
            if response.double_clicked() {
                app.viewer3d.camera.reset();
            }
            if response.dragged_by(egui::PointerButton::Primary) {
                let delta = response.drag_delta();
                let shift = vp_ctx.input(|i| i.modifiers.shift);
                if shift {
                    app.viewer3d.camera.pan_x += delta.x;
                    app.viewer3d.camera.pan_y += delta.y;
                } else {
                    app.viewer3d.camera.yaw += delta.x * 0.011;
                    app.viewer3d.camera.pitch = (app.viewer3d.camera.pitch + delta.y * 0.011)
                        .clamp(-std::f32::consts::FRAC_PI_2 + 0.05, std::f32::consts::FRAC_PI_2 - 0.05);
                }
            } else if response.dragged_by(egui::PointerButton::Secondary)
                || response.dragged_by(egui::PointerButton::Middle)
            {
                let delta = response.drag_delta();
                app.viewer3d.camera.pan_x += delta.x;
                app.viewer3d.camera.pan_y += delta.y;
            }

            if response.hovered() {
                let scroll_y = vp_ctx.input(|i| i.smooth_scroll_delta.y);
                if scroll_y.abs() > 0.1 {
                    let zoom_factor = (1.0 + scroll_y * 0.0025).clamp(0.75, 1.35);
                    app.viewer3d.camera.zoom = (app.viewer3d.camera.zoom * zoom_factor).clamp(0.3, 6.0);
                }
            }

            if app.viewer3d.camera.auto_rotate {
                app.viewer3d.camera.yaw += 0.008;
                vp_ctx.request_repaint();
            }

            // Compute active annotation residue span (if an annotation sub-region is linked).
            let ann_res_span = compute_annotation_residue_span(app);

            // Sync hover from main alignment window if not currently hovering directly in 3D.
            let aln_synced_hover = app
                .hovered_column
                .and_then(|c| app.viewer3d.alignment_map.as_ref()?.residue_for_column(c));

            let highlights = RenderHighlights {
                hovered_residue: app.viewer3d.hovered_residue.or(aln_synced_hover),
                selected_residue: app.viewer3d.selected_residue,
                annotation_residue_span: ann_res_span,
            };

            if let Some(ref structure) = app.viewer3d.structure {
                let pointer_pos = response.hover_pos();
                let pick = app.viewer3d.renderer.render(
                    &painter,
                    rect,
                    structure,
                    app.viewer3d.alignment_map.as_ref(),
                    &app.viewer3d.camera,
                    app.viewer3d.style,
                    app.viewer3d.color_by,
                    &highlights,
                    pointer_pos,
                );
                app.viewer3d.hovered_residue = pick.hovered_residue;

                // Clicking a residue in 3D selects it and jumps the main alignment viewport right to it!
                if response.clicked() {
                    if let Some(res_idx) = pick.hovered_residue {
                        app.viewer3d.selected_residue = Some(res_idx);
                        if let Some(Some((col_start, col_end))) = app
                            .viewer3d
                            .alignment_map
                            .as_ref()
                            .and_then(|m| m.residue_to_cols.get(res_idx))
                        {
                            app.selected_cols = Some((*col_start, *col_end));
                            app.scroll_to_column = Some(*col_start);
                        }
                    }
                }

                // Floating HUD overlay in top-left of 3D canvas for hovered/selected residue.
                let inspect_idx = app.viewer3d.hovered_residue.or(app.viewer3d.selected_residue);
                if let Some(res_idx) = inspect_idx {
                    if let Some(res) = structure.residues.get(res_idx) {
                        let msa_info = app
                            .viewer3d
                            .alignment_map
                            .as_ref()
                            .and_then(|m| {
                                let cols = m.residue_to_cols.get(res_idx).copied().flatten()?;
                                let cons = m.residue_conservation.get(res_idx).copied().flatten();
                                Some((cols, cons))
                            });

                        let mut hud_text = format!(
                            "{} {} (Chain {})  |  {}  |  {}: {:.1}",
                            res.three_letter(),
                            res.seq_id,
                            res.chain_id,
                            res.secondary_structure.label(),
                            if structure.is_alphafold { "pLDDT" } else { "B-factor" },
                            res.plddt
                        );
                        if let Some(((c0, c1), Some(c))) = msa_info {
                            hud_text.push_str(&format!(
                                "  |  MSA Col {}..{} ({:.0}% id)",
                                c0 + 1,
                                c1,
                                c
                            ));
                        }

                        let hud_pos = Pos2::new(rect.min.x + 12.0, rect.min.y + 12.0);
                        let galley = painter.layout_no_wrap(
                            hud_text,
                            FontId::monospace(12.0),
                            Color32::from_rgb(240, 246, 255),
                        );
                        let pill = egui::Rect::from_min_size(
                            hud_pos,
                            galley.size() + egui::vec2(16.0, 10.0),
                        );
                        painter.rect_filled(pill, 6.0, Color32::from_rgba_premultiplied(15, 23, 42, 215));
                        painter.rect_stroke(
                            pill,
                            6.0,
                            Stroke::new(1.0, Color32::from_rgb(56, 189, 248)),
                        );
                        painter.galley(hud_pos + egui::vec2(8.0, 5.0), galley, Color32::WHITE);
                    }
                }

                // Subtle controls hint in bottom-left.
                painter.text(
                    Pos2::new(rect.min.x + 12.0, rect.max.y - 12.0),
                    Align2::LEFT_BOTTOM,
                    "Left-Drag: Orbit  |  Right/Shift-Drag: Pan  |  Scroll: Zoom  |  Click Residue: Jump in Alignment",
                    FontId::proportional(11.0),
                    Color32::from_rgba_unmultiplied(200, 214, 230, 175),
                );
            } else {
                // Empty state inside 3D canvas while loading or awaiting fetch.
                let msg = app
                    .viewer3d
                    .loading_message
                    .as_deref()
                    .unwrap_or("Enter a UniProt accession above and click 'Fetch AlphaFold DB' — or open the Offline Preview.");
                painter.text(
                    rect.center(),
                    Align2::CENTER_CENTER,
                    msg,
                    FontId::proportional(14.0),
                    Color32::from_rgb(203, 213, 225),
                );
            }

            // Loading overlay pill when fetching a new structure over an existing one.
            if let Some(loading) = &app.viewer3d.loading_message {
                let banner = egui::Rect::from_center_size(
                    Pos2::new(rect.center().x, rect.max.y - 44.0),
                    egui::vec2(380.0, 30.0),
                );
                painter.rect_filled(banner, 6.0, Color32::from_rgba_premultiplied(16, 138, 134, 230));
                painter.text(
                    banner.center(),
                    Align2::CENTER_CENTER,
                    loading,
                    FontId::proportional(12.0),
                    Color32::WHITE,
                );
            }
        });
}

fn show_3d_metadata_sidebar(app: &mut HelixApp, vp_ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let mut fetch_pdb: Option<String> = None;
    let mut load_offline_preview = false;

    ScrollArea::vertical().show(ui, |ui| {
        ui.label(
            RichText::new("3D STRUCTURE & DATABASE")
                .strong()
                .small()
                .color(pal.text_muted),
        );
        ui.add_space(4.0);

        if let Some(err) = &app.viewer3d.error_message {
            Frame::none()
                .fill(Color32::from_rgba_premultiplied(212, 76, 76, 35))
                .stroke(Stroke::new(1.0, pal.identity_low))
                .rounding(Rounding::same(5.0))
                .inner_margin(egui::Margin::same(8.0))
                .show(ui, |ui| {
                    ui.label(RichText::new("Fetch Error").strong().color(pal.identity_low));
                    ui.label(RichText::new(err).small());
                    ui.add_space(4.0);
                    if ui.small_button("Load Offline Preview Fold").clicked() {
                        load_offline_preview = true;
                    }
                });
            ui.add_space(8.0);
        }

        if let Some(st) = &app.viewer3d.structure {
            Frame::none()
                .fill(pal.canvas_bg)
                .stroke(Stroke::new(1.0, pal.border))
                .rounding(Rounding::same(6.0))
                .inner_margin(egui::Margin::same(8.0))
                .show(ui, |ui| {
                    ui.label(RichText::new(&st.title).strong().size(12.5));
                    if let Some(org) = &st.organism {
                        ui.label(RichText::new(org).small().italics().color(pal.text_muted));
                    }
                    ui.add_space(4.0);
                    egui::Grid::new("meta_3d_grid")
                        .num_columns(2)
                        .spacing([8.0, 3.0])
                        .show(ui, |ui| {
                            ui.label(RichText::new("Entry:").small().color(pal.text_muted));
                            ui.label(RichText::new(&st.entry_id).small().monospace().strong());
                            ui.end_row();

                            ui.label(RichText::new("Source:").small().color(pal.text_muted));
                            ui.label(RichText::new(&st.source_label).small());
                            ui.end_row();

                            ui.label(RichText::new("Residues:").small().color(pal.text_muted));
                            ui.label(RichText::new(st.residues.len().to_string()).small().strong());
                            ui.end_row();

                            ui.label(RichText::new("Heavy Atoms:").small().color(pal.text_muted));
                            ui.label(RichText::new(st.atoms.len().to_string()).small());
                            ui.end_row();

                            if st.is_alphafold {
                                ui.label(RichText::new("Mean pLDDT:").small().color(pal.text_muted));
                                ui.label(
                                    RichText::new(format!("{:.1} / 100", st.mean_plddt))
                                        .small()
                                        .strong()
                                        .color(plddt_color(st.mean_plddt)),
                                );
                                ui.end_row();
                            }

                            if let Some(map) = &app.viewer3d.alignment_map {
                                ui.label(RichText::new("MSA Mapped:").small().color(pal.text_muted));
                                ui.label(
                                    RichText::new(format!("{} residues", map.mapped_count))
                                        .small()
                                        .strong()
                                        .color(pal.sync_3d_highlight),
                                );
                                ui.end_row();
                            }
                        });
                });
        }

        // Color legend matching active 3D coloring mode.
        ui.add_space(10.0);
        ui.label(
            RichText::new("COLOR LEGEND")
                .strong()
                .small()
                .color(pal.text_muted),
        );
        ui.add_space(4.0);
        Frame::none()
            .fill(pal.canvas_bg)
            .stroke(Stroke::new(1.0, pal.border))
            .rounding(Rounding::same(6.0))
            .inner_margin(egui::Margin::same(8.0))
            .show(ui, |ui| match app.viewer3d.color_by {
                ColorBy3D::AlphaFoldPlddt => {
                    legend_row(ui, plddt_color(95.0), "Very high (pLDDT > 90)");
                    legend_row(ui, plddt_color(80.0), "Confident (70 – 90)");
                    legend_row(ui, plddt_color(60.0), "Low (50 – 70)");
                    legend_row(ui, plddt_color(35.0), "Very low (< 50)");
                }
                ColorBy3D::SecondaryStructure => {
                    legend_row(ui, Color32::from_rgb(230, 86, 134), "Alpha Helix");
                    legend_row(ui, Color32::from_rgb(240, 186, 46), "Beta Strand");
                    legend_row(ui, Color32::from_rgb(92, 168, 224), "Coil / Loop");
                }
                ColorBy3D::AlignmentConservation => {
                    legend_row(ui, Color32::from_rgb(195, 30, 95), "100% Conserved");
                    legend_row(ui, Color32::from_rgb(238, 204, 84), "Moderate (~60%)");
                    legend_row(ui, Color32::from_rgb(72, 168, 224), "Variable (<40%)");
                    legend_row(ui, Color32::from_rgb(122, 132, 148), "Unmapped region");
                }
                ColorBy3D::Rainbow => {
                    legend_row(ui, Color32::from_rgb(60, 120, 240), "N-Terminus");
                    legend_row(ui, Color32::from_rgb(80, 210, 120), "Mid-Chain");
                    legend_row(ui, Color32::from_rgb(235, 70, 70), "C-Terminus");
                }
                ColorBy3D::ChainOrElement => {
                    legend_row(ui, Color32::from_rgb(72, 152, 232), "Chain A / N (blue)");
                    legend_row(ui, Color32::from_rgb(236, 78, 78), "Oxygen (red)");
                    legend_row(ui, Color32::from_rgb(238, 206, 62), "Sulfur (yellow)");
                    legend_row(ui, Color32::from_rgb(148, 142, 212), "Metal / Ligand");
                }
            });

        // UniProt function summary & experimental RCSB PDB quick-switches.
        if let Some(up) = &app.viewer3d.uniprot_entry {
            ui.add_space(10.0);
            ui.label(
                RichText::new(format!("UNIPROTKB ({})", up.accession))
                    .strong()
                    .small()
                    .color(pal.text_muted),
            );
            ui.add_space(4.0);
            Frame::none()
                .fill(pal.canvas_bg)
                .stroke(Stroke::new(1.0, pal.border))
                .rounding(Rounding::same(6.0))
                .inner_margin(egui::Margin::same(8.0))
                .show(ui, |ui| {
                    if !up.function_summary.is_empty() {
                        let snippet: String = up.function_summary.chars().take(280).collect();
                        let display = if up.function_summary.len() > 280 {
                            format!("{snippet}...")
                        } else {
                            snippet
                        };
                        ui.label(RichText::new(display).size(11.0).color(pal.text_primary));
                    }

                    if !up.pdb_ids.is_empty() {
                        ui.add_space(6.0);
                        ui.label(
                            RichText::new("Experimental RCSB PDB Structures:")
                                .small()
                                .strong()
                                .color(pal.text_muted),
                        );
                        ui.horizontal_wrapped(|ui| {
                            for pdb in up.pdb_ids.iter().take(6) {
                                if ui.small_button(pdb).clicked() {
                                    fetch_pdb = Some(pdb.clone());
                                }
                            }
                        });
                    }
                });
        }
    });

    if load_offline_preview {
        app.open_offline_3d_preview();
    }
    if let Some(pdb_id) = fetch_pdb {
        app.fetch_rcsb_structure(vp_ctx, &pdb_id);
    }
}

fn legend_row(ui: &mut Ui, swatch: Color32, text: &str) {
    ui.horizontal(|ui| {
        let (r, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), Sense::hover());
        ui.painter().rect_filled(r, 3.0, swatch);
        ui.label(RichText::new(text).size(11.0));
    });
}

/// Computes the 3D residue index range corresponding to the currently linked alignment annotation.
fn compute_annotation_residue_span(app: &HelixApp) -> Option<(usize, usize)> {
    let ann_id = app.viewer3d.linked_annotation_id.or(app.selected_annotation_id)?;
    let aln = app.current_alignment()?;
    let ann = aln.annotation_by_id(ann_id)?;
    let map = app.viewer3d.alignment_map.as_ref()?;

    let mut min_res: Option<usize> = None;
    let mut max_res: Option<usize> = None;
    for (res_idx, entry) in map.residue_to_cols.iter().enumerate() {
        if let Some((c0, c1)) = entry {
            if *c1 > ann.start && *c0 < ann.end {
                min_res = Some(min_res.map_or(res_idx, |m| m.min(res_idx)));
                max_res = Some(max_res.map_or(res_idx + 1, |m| m.max(res_idx + 1)));
            }
        }
    }
    Some((min_res?, max_res?))
}
