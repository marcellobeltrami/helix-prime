//! High-performance Geneious Prime-inspired multiple sequence alignment inspection canvas.
//!
//! Renders an interactive Overview minimap, Identity % histogram, Consensus track, Coordinate
//! ruler, packed directional Annotation tracks (with UniProt 3D launch badges), optional
//! nucleotide CDS translation track, and horizontally virtualized sequence rows.

use crate::app::{HelixApp, InspectorTab};
use crate::theme::{residue_colors, Palette};
use egui::{
    Align2, Color32, Context, FontId, Frame, Pos2, Rect, RichText, Sense, Shape, Stroke, Ui, Vec2,
};
use helix_core::sequence::{complement_nucleotide, is_gap, translate_aligned_region};
use helix_core::{pack_annotation_lanes, Strand};

const GUTTER_WIDTH: f32 = 165.0;
const OVERVIEW_HEIGHT: f32 = 30.0;
const IDENTITY_HEIGHT: f32 = 26.0;
const CONSENSUS_HEIGHT: f32 = 22.0;
const RULER_HEIGHT: f32 = 20.0;
const RESTRICTION_HEIGHT: f32 = 22.0;
const ANN_LANE_HEIGHT: f32 = 22.0;
const TRANSLATION_HEIGHT: f32 = 20.0;

/// Renders the alignment body (Overview + Linear or Wrapped canvas) inside the Document Viewer.
pub fn show_alignment_view_body(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    if app.current_alignment().is_none() {
        return;
    }

    // Handle direct keyboard editing when `Allow Editing` is enabled and no text box has focus.
    handle_direct_sequence_editing(app, ctx);

    // Ensure restriction sites cache is warm when overlay is toggled on.
    if app.display.show_restriction_sites {
        if let Some(doc) = app.current_doc_mut() {
            let _ = doc.ensure_restriction_sites();
        }
    }

    if app.display.wrap_sequence {
        show_wrapped_view(app, ui, &pal);
        return;
    }

    if app.display.show_overview {
        show_overview_bar(app, ui, &pal);
    }

    show_main_grid(app, ctx, ui, &pal);
}

/// Processes keyboard input for direct sequence editing (`Allow Editing`) and cursor movement.
fn handle_direct_sequence_editing(app: &mut HelixApp, ctx: &Context) {
    if ctx.wants_keyboard_input() {
        return;
    }
    let Some(doc) = app.docs.get(app.active_doc) else {
        return;
    };
    let width = doc.alignment.width();
    if width == 0 {
        return;
    }
    let is_nuc = doc.alignment.alphabet.is_nucleotide();

    // Left / Right arrow cursor navigation always works nicely.
    let (left, right, backspace, delete, typed_chars) = ctx.input(|i| {
        let mut chars = Vec::new();
        for ev in &i.events {
            if let egui::Event::Text(txt) = ev {
                for ch in txt.chars() {
                    if ch.is_ascii_alphabetic() || ch == '-' {
                        chars.push(ch.to_ascii_uppercase() as u8);
                    }
                }
            }
        }
        (
            i.key_pressed(egui::Key::ArrowLeft),
            i.key_pressed(egui::Key::ArrowRight),
            i.key_pressed(egui::Key::Backspace),
            i.key_pressed(egui::Key::Delete),
            chars,
        )
    });

    if left && app.cursor_col > 0 {
        app.cursor_col -= 1;
        app.scroll_to_column = Some(app.cursor_col);
    }
    if right && app.cursor_col + 1 < width {
        app.cursor_col += 1;
        app.scroll_to_column = Some(app.cursor_col);
    }

    if !app.display.allow_editing {
        return;
    }

    let row = app.cursor_row;
    let mut mutated = false;
    if backspace && app.cursor_col > 0 {
        let target = app.cursor_col - 1;
        if let Some(doc) = app.current_doc_mut() {
            if doc.alignment.delete_column_at(target) {
                doc.invalidate_caches();
                app.cursor_col = target;
                mutated = true;
            }
        }
    } else if delete {
        let target = app.cursor_col;
        let new_width = app.current_doc_mut().and_then(|doc| {
            if doc.alignment.delete_column_at(target) {
                doc.invalidate_caches();
                Some(doc.alignment.width())
            } else {
                None
            }
        });
        if let Some(w) = new_width {
            app.cursor_col = target.min(w.saturating_sub(1));
            mutated = true;
        }
    } else if !typed_chars.is_empty() {
        for b in typed_chars {
            let valid = if is_nuc {
                matches!(b, b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'-')
            } else {
                b.is_ascii_alphabetic() || b == b'-'
            };
            if valid {
                let col = app.cursor_col;
                if let Some(doc) = app.current_doc_mut() {
                    if doc.alignment.insert_residue_at(row, col, b) {
                        doc.invalidate_caches();
                        app.cursor_col = col.saturating_add(1);
                        mutated = true;
                    }
                }
            }
        }
    }

    if mutated {
        app.rebuild_3d_alignment_map();
        app.status_message = format!("Edited sequence at base {}.", app.cursor_col + 1);
    }
}

/// Renders the multi-line Wrapped Sequence View (`Wrap` checkbox in Geneious Prime).
fn show_wrapped_view(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    let Some(doc) = app.docs.get(app.active_doc) else {
        return;
    };
    let aln = &doc.alignment;
    let width = aln.width();
    let is_nuc = aln.alphabet.is_nucleotide();
    let show_comp = app.display.show_complement && is_nuc;
    let col_w = app.display.col_width.clamp(8.0, 22.0);
    let avail_w = (ui.available_width() - 120.0).max(160.0);
    let cols_per_line = ((avail_w / col_w).floor() as usize).clamp(20, 120);
    let blocks = width.div_ceil(cols_per_line);

    let mut clicked_col: Option<usize> = None;

    egui::ScrollArea::vertical()
        .id_salt("wrapped_seq_scroll")
        .show(ui, |ui| {
            ui.add_space(6.0);
            for block in 0..blocks {
                let start = block * cols_per_line;
                let end = (start + cols_per_line).min(width);
                let len = end - start;

                let row_h = app.display.row_height;
                let seq_rows = aln.row_count();
                let per_seq_h = if show_comp { row_h * 1.85 } else { row_h };
                let block_h = RULER_HEIGHT + seq_rows as f32 * per_seq_h + 10.0;

                let (rect, resp) = ui.allocate_exact_size(
                    Vec2::new(ui.available_width(), block_h),
                    Sense::click(),
                );
                let painter = ui.painter_at(rect);
                let left_x = rect.min.x + 110.0;

                // Coordinate ruler for this wrapped block.
                for local in 0..len {
                    let col = start + local;
                    if (col + 1) % 10 == 0 || local == 0 {
                        let x = left_x + local as f32 * col_w + col_w * 0.5;
                        painter.text(
                            Pos2::new(x, rect.min.y + 8.0),
                            Align2::CENTER_CENTER,
                            (col + 1).to_string(),
                            FontId::monospace(9.5),
                            pal.text_muted,
                        );
                    }
                }

                let mono = FontId::monospace((col_w * 0.78).clamp(8.5, 14.0));
                for (r_idx, seq) in aln.sequences().iter().enumerate() {
                    let ry = rect.min.y + RULER_HEIGHT + r_idx as f32 * per_seq_h;
                    let short_name: String = seq.name.chars().take(13).collect();
                    painter.text(
                        Pos2::new(rect.min.x + 8.0, ry + row_h * 0.5),
                        Align2::LEFT_CENTER,
                        short_name,
                        FontId::proportional(11.0),
                        pal.text_primary,
                    );

                    for local in 0..len {
                        let col = start + local;
                        let x0 = left_x + local as f32 * col_w;
                        let cell = Rect::from_min_max(
                            Pos2::new(x0, ry + 1.0),
                            Pos2::new(x0 + col_w, ry + row_h - 1.0),
                        );
                        let res = seq.residues[col];
                        let cons = aln.columns()[col].consensus;
                        let (bg, fg) = residue_colors(
                            res,
                            cons,
                            aln.alphabet,
                            app.display.color_scheme,
                            app.display.highlight_mode,
                            app.theme,
                        );
                        if let Some(bg_c) = bg {
                            painter.rect_filled(cell, 1.0, bg_c);
                        }
                        let b = [res];
                        painter.text(
                            cell.center(),
                            Align2::CENTER_CENTER,
                            std::str::from_utf8(&b).unwrap_or("?"),
                            mono.clone(),
                            fg,
                        );

                        if show_comp {
                            let comp = [complement_nucleotide(res)];
                            painter.text(
                                Pos2::new(cell.center().x, ry + row_h + row_h * 0.38),
                                Align2::CENTER_CENTER,
                                std::str::from_utf8(&comp).unwrap_or(" "),
                                FontId::monospace((col_w * 0.68).clamp(8.0, 12.0)),
                                pal.text_muted,
                            );
                        }

                        if col == app.cursor_col {
                            painter.line_segment(
                                [Pos2::new(x0, ry), Pos2::new(x0, ry + per_seq_h - 2.0)],
                                Stroke::new(2.0, pal.accent),
                            );
                        }
                    }
                }

                if resp.clicked() {
                    if let Some(ptr) = resp.interact_pointer_pos() {
                        if ptr.x >= left_x {
                            let local = ((ptr.x - left_x) / col_w) as usize;
                            if local < len {
                                clicked_col = Some(start + local);
                            }
                        }
                    }
                }
            }
        });

    if let Some(col) = clicked_col {
        app.cursor_col = col;
        app.selected_cols = Some((col, col + 1));
    }
}

/// Renders the central alignment document header, overview minimap, and main sequence canvas.
#[allow(dead_code)]
pub fn show_alignment_view(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    if app.current_alignment().is_none() {
        ui.centered_and_justified(|ui| {
            ui.label(
                RichText::new("No alignment loaded. Click 'Import...' to open a FASTA, Clustal, or GenBank file.")
                    .color(pal.text_muted),
            );
        });
        return;
    }

    show_document_header(app, ctx, ui, &pal);
    show_alignment_view_body(app, ctx, ui);
}

fn show_document_header(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let mut launch_3d_req: Option<(String, Option<u64>)> = None;

    Frame::none()
        .fill(pal.panel_bg)
        .stroke(Stroke::new(1.0, pal.border))
        .inner_margin(egui::Margin::symmetric(12.0, 6.0))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                if let Some(aln) = app.current_alignment() {
                    ui.label(
                        RichText::new(&aln.name)
                            .strong()
                            .size(14.0)
                            .color(pal.text_primary),
                    );
                    let unit = if aln.alphabet.is_nucleotide() {
                        "bp"
                    } else {
                        "aa"
                    };
                    ui.label(
                        RichText::new(format!(
                            "{}  |  {} seqs  |  {} {unit}  |  Pairwise ID: {:.1}%  |  Identical sites: {:.1}%",
                            aln.alphabet.label(),
                            aln.row_count(),
                            aln.width(),
                            aln.mean_pairwise_identity_pct(),
                            aln.identical_sites_pct()
                        ))
                        .small()
                        .color(pal.text_muted),
                    );
                }

                // Right side: quick 3D CTA when a UniProt-linked annotation is selected.
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if let Some(ann_id) = app.selected_annotation_id {
                        if let Some(ann) = app
                            .current_alignment()
                            .and_then(|a| a.annotation_by_id(ann_id))
                        {
                            if let Some(up) = &ann.uniprot {
                                let btn = egui::Button::new(
                                    RichText::new(format!("Open 3D AlphaFold ({})", up.accession))
                                        .small()
                                        .strong()
                                        .color(Color32::WHITE),
                                )
                                .fill(Color32::from_rgb(16, 138, 134));
                                if ui.add(btn).clicked() {
                                    launch_3d_req = Some((up.accession.clone(), Some(ann.id)));
                                }
                            }
                            ui.label(
                                RichText::new(format!(
                                    "Annotation: {} ({}..{})",
                                    ann.name,
                                    ann.start + 1,
                                    ann.end
                                ))
                                .small()
                                .strong()
                                .color(pal.accent),
                            );
                        }
                    }
                });
            });
        });

    if let Some((acc, ann_id)) = launch_3d_req {
        app.open_3d_for_accession(ctx, &acc, ann_id);
    }
}

fn show_overview_bar(app: &mut HelixApp, ui: &mut Ui, pal: &Palette) {
    let Some(aln) = app.current_alignment() else {
        return;
    };
    let width_cols = aln.width().max(1);
    let col_w = app.display.col_width;

    let avail_w = ui.available_width();
    let (rect, response) =
        ui.allocate_exact_size(Vec2::new(avail_w, OVERVIEW_HEIGHT), Sense::click_and_drag());
    let painter = ui.painter_at(rect);

    painter.rect_filled(rect, 0.0, pal.panel_bg);
    painter.line_segment(
        [rect.left_bottom(), rect.right_bottom()],
        Stroke::new(1.0, pal.border),
    );

    // Left gutter label.
    let gutter_rect = Rect::from_min_max(rect.min, Pos2::new(rect.min.x + GUTTER_WIDTH, rect.max.y));
    painter.text(
        Pos2::new(gutter_rect.min.x + 12.0, gutter_rect.center().y),
        Align2::LEFT_CENTER,
        "Overview",
        FontId::proportional(11.5),
        pal.text_muted,
    );
    painter.line_segment(
        [gutter_rect.right_top(), gutter_rect.right_bottom()],
        Stroke::new(1.0, pal.border),
    );

    let track_rect = Rect::from_min_max(
        Pos2::new(gutter_rect.max.x + 6.0, rect.min.y + 4.0),
        Pos2::new(rect.max.x - 8.0, rect.max.y - 4.0),
    );
    if track_rect.width() <= 4.0 {
        return;
    }

    painter.rect_filled(track_rect, 3.0, pal.canvas_bg);

    // Downsample alignment columns to pixel buckets across `track_rect.width()`.
    let px_bins = (track_rect.width() as usize).clamp(1, 800);
    let bin_w = track_rect.width() / px_bins as f32;
    let bar_max_h = track_rect.height() * 0.62;
    let bar_bottom = track_rect.min.y + bar_max_h + 1.0;

    for bin in 0..px_bins {
        let c0 = (bin * width_cols) / px_bins;
        let c1 = (((bin + 1) * width_cols) / px_bins).max(c0 + 1).min(width_cols);
        let mut sum_id = 0.0f32;
        for c in c0..c1 {
            sum_id += aln.columns()[c].identity_pct;
        }
        let avg_id = sum_id / (c1 - c0) as f32;
        let h = (avg_id / 100.0) * bar_max_h;
        let bar_col = identity_bar_color(avg_id, pal);
        let x0 = track_rect.min.x + bin as f32 * bin_w;
        painter.rect_filled(
            Rect::from_min_max(
                Pos2::new(x0, bar_bottom - h),
                Pos2::new(x0 + bin_w.max(1.0), bar_bottom),
            ),
            0.0,
            bar_col,
        );
    }

    // Draw miniature annotation bars along the bottom strip of the overview.
    let ann_y0 = bar_bottom + 2.0;
    let ann_y1 = track_rect.max.y - 1.0;
    for ann in aln.annotations() {
        let x0 = track_rect.min.x + (ann.start as f32 / width_cols as f32) * track_rect.width();
        let x1 = (track_rect.min.x + (ann.end as f32 / width_cols as f32) * track_rect.width())
            .max(x0 + 2.5);
        let [r, g, b] = ann.rgb();
        painter.rect_filled(
            Rect::from_min_max(Pos2::new(x0, ann_y0), Pos2::new(x1, ann_y1)),
            1.0,
            Color32::from_rgb(r, g, b),
        );
    }

    // Draw visible viewport indicator box.
    let total_px = (width_cols as f32 * col_w).max(1.0);
    let view_w_px = (avail_w - GUTTER_WIDTH).max(1.0);
    let frac_start = (app.horizontal_scroll_px / total_px).clamp(0.0, 1.0);
    let frac_len = (view_w_px / total_px).clamp(0.02, 1.0);
    let vx0 = track_rect.min.x + frac_start * track_rect.width();
    let vx1 = (vx0 + frac_len * track_rect.width()).min(track_rect.max.x);
    let view_box = Rect::from_min_max(
        Pos2::new(vx0, track_rect.min.y),
        Pos2::new(vx1, track_rect.max.y),
    );
    painter.rect_filled(view_box, 2.0, pal.selection_fill);
    painter.rect_stroke(view_box, 2.0, Stroke::new(1.5, pal.accent));

    // Click or drag in the Overview bar to pan the alignment viewport directly!
    if response.clicked() || response.dragged() {
        if let Some(pos) = response.interact_pointer_pos() {
            if track_rect.contains(Pos2::new(pos.x, track_rect.center().y)) {
                let frac = ((pos.x - track_rect.min.x) / track_rect.width()).clamp(0.0, 1.0);
                let target_px = frac * total_px - view_w_px * 0.5;
                let max_scroll = (total_px - view_w_px).max(0.0);
                app.horizontal_scroll_px = target_px.clamp(0.0, max_scroll);
            }
        }
    }
}

fn show_main_grid(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let Some(doc) = app.docs.get(app.active_doc) else {
        return;
    };

    let width_cols = doc.alignment.width();
    let row_count = doc.alignment.row_count();
    let is_nuc = doc.alignment.alphabet.is_nucleotide();

    // Handle Ctrl+Wheel zoom and horizontal scroll BEFORE taking a long-lived immutable borrow on `doc.alignment`.
    let avail_rect = ui.available_rect_before_wrap();
    let data_viewport_w = (avail_rect.width() - GUTTER_WIDTH).max(60.0);

    if ui.rect_contains_pointer(avail_rect) {
        let (scroll_delta, ctrl_down, shift_down) = ctx.input(|i| {
            (
                i.smooth_scroll_delta,
                i.modifiers.ctrl || i.modifiers.command,
                i.modifiers.shift,
            )
        });
        if ctrl_down && scroll_delta.y.abs() > 0.1 {
            let factor = if scroll_delta.y > 0.0 { 1.1 } else { 0.91 };
            app.display.col_width = (app.display.col_width * factor).clamp(4.0, 26.0);
        } else if scroll_delta.x.abs() > 0.1 {
            app.horizontal_scroll_px -= scroll_delta.x;
        } else if shift_down && scroll_delta.y.abs() > 0.1 {
            app.horizontal_scroll_px -= scroll_delta.y;
        }
    }

    let col_w = app.display.col_width;
    let row_h = app.display.row_height;
    let total_seq_w = (width_cols as f32 * col_w).max(1.0);
    let max_h_scroll = (total_seq_w - data_viewport_w).max(0.0);

    // Honor programmatic scroll requests (e.g., motif jump or click in 3D protein window!).
    if let Some(target_col) = app.scroll_to_column.take() {
        let target_px = target_col as f32 * col_w - data_viewport_w * 0.35;
        app.horizontal_scroll_px = target_px.clamp(0.0, max_h_scroll);
    }
    app.horizontal_scroll_px = app.horizontal_scroll_px.clamp(0.0, max_h_scroll);

    // Now borrow `doc` and `aln` immutably across the painting pass.
    let doc = &app.docs[app.active_doc];
    let aln = &doc.alignment;
    let restriction_sites = doc.restriction_cache.as_deref().unwrap_or(&[]);

    // Pack annotations into non-overlapping lanes.
    let (ann_lanes, ann_lane_map) = if app.display.show_annotations {
        pack_annotation_lanes(aln.annotations(), 1)
    } else {
        (0, Vec::new())
    };

    let show_trans = app.display.show_translation && is_nuc;
    let show_restr = app.display.show_restriction_sites && is_nuc;
    let show_comp = app.display.show_complement && is_nuc;
    let effective_row_h = if show_comp { row_h * 1.85 } else { row_h };

    // Compute vertical layout heights.
    let identity_h = if app.display.show_identity_graph {
        IDENTITY_HEIGHT
    } else {
        0.0
    };
    let consensus_h = if app.display.show_consensus {
        CONSENSUS_HEIGHT
    } else {
        0.0
    };
    let restr_h = if show_restr { RESTRICTION_HEIGHT } else { 0.0 };
    let ann_h = ann_lanes as f32 * ANN_LANE_HEIGHT;
    let trans_h = if show_trans { TRANSLATION_HEIGHT } else { 0.0 };
    let header_tracks_h = identity_h + consensus_h + RULER_HEIGHT + restr_h + ann_h + trans_h;
    let sequences_h = row_count as f32 * effective_row_h;

    // Determine column highlighted via bidirectional sync from the 3D Protein Window!
    let synced_3d_cols: Option<(usize, usize)> = app.viewer3d.alignment_map.as_ref().and_then(|map| {
        let active_res = app.viewer3d.hovered_residue.or(app.viewer3d.selected_residue)?;
        map.residue_to_cols.get(active_res).copied().flatten()
    });

    // Reserve a crisp horizontal scrollbar at the very bottom of the canvas.
    let scrollbar_h = 16.0;
    let scroll_area_h = (avail_rect.height() - scrollbar_h).max(120.0);

    // Deferred actions triggered during painting/hit-testing so we don't hold immutable borrow of `aln`.
    let mut clicked_ann: Option<(u64, usize, usize)> = None;
    let mut double_clicked_ann_3d: Option<(u64, String)> = None;
    let mut new_hover_col: Option<usize> = None;
    let mut new_hover_row: Option<usize> = None;
    let mut new_selected_cols: Option<(usize, usize)> = app.selected_cols;
    let mut clicked_col_sync_3d: Option<usize> = None;

    egui::ScrollArea::vertical()
        .id_salt("alignment_vertical_scroll")
        .max_height(scroll_area_h)
        .show(ui, |ui| {
            let content_h = (header_tracks_h + sequences_h + 16.0).max(scroll_area_h);
            let (canvas_rect, response) = ui.allocate_exact_size(
                Vec2::new(avail_rect.width(), content_h),
                Sense::click_and_drag(),
            );
            let painter = ui.painter_at(canvas_rect);

            let data_left = canvas_rect.min.x + GUTTER_WIDTH;
            let data_rect = Rect::from_min_max(
                Pos2::new(data_left, canvas_rect.min.y),
                canvas_rect.max,
            );
            let data_painter = ui.painter_at(data_rect);

            // Compute visible column window (horizontal virtualization -> O(visible_cols) per frame).
            let scroll_x = app.horizontal_scroll_px;
            let first_col = ((scroll_x / col_w).floor() as usize).min(width_cols);
            let vis_cols = ((data_rect.width() / col_w).ceil() as usize) + 2;
            let last_col = (first_col + vis_cols).min(width_cols);

            let col_to_x = |col: usize| -> f32 { data_left + col as f32 * col_w - scroll_x };

            // Determine hovered column and row under mouse pointer.
            if let Some(ptr) = response.hover_pos() {
                if ptr.x >= data_left && ptr.x <= canvas_rect.max.x && width_cols > 0 {
                    let rel_x = ptr.x - data_left + scroll_x;
                    if rel_x >= 0.0 {
                        let c = (rel_x / col_w) as usize;
                        if c < width_cols {
                            new_hover_col = Some(c);
                        }
                    }
                }
                let seq_top = canvas_rect.min.y + header_tracks_h;
                if ptr.y >= seq_top {
                    let r = ((ptr.y - seq_top) / effective_row_h) as usize;
                    if r < row_count {
                        new_hover_row = Some(r);
                    }
                }
            }

            // Drag or click in sequence/ruler area to select columns & sync with 3D viewer!
            let ann_top_y = canvas_rect.min.y + identity_h + consensus_h + RULER_HEIGHT + restr_h;
            let ann_bottom_y = ann_top_y + ann_h;

            if response.clicked() {
                if let (Some(ptr), Some(col)) = (response.interact_pointer_pos(), new_hover_col) {
                    if ptr.y < ann_top_y || ptr.y >= ann_bottom_y {
                        new_selected_cols = Some((col, col + 1));
                        clicked_col_sync_3d = Some(col);
                    }
                }
            } else if response.dragged() {
                if let (Some(origin), Some(curr)) = (
                    ui.input(|i| i.pointer.press_origin()),
                    response.interact_pointer_pos(),
                ) {
                    if origin.x >= data_left && (origin.y < ann_top_y || origin.y >= ann_bottom_y) {
                        let c0 = (((origin.x - data_left + scroll_x).max(0.0)) / col_w) as usize;
                        let c1 = (((curr.x - data_left + scroll_x).max(0.0)) / col_w) as usize;
                        let lo = c0.min(c1).min(width_cols.saturating_sub(1));
                        let hi = (c0.max(c1) + 1).min(width_cols);
                        new_selected_cols = Some((lo, hi));
                    }
                }
            }

            // Paint background & vertical gutter divider.
            painter.rect_filled(canvas_rect, 0.0, pal.canvas_bg);
            painter.rect_filled(
                Rect::from_min_max(canvas_rect.min, Pos2::new(data_left, canvas_rect.max.y)),
                0.0,
                pal.panel_bg,
            );

            // Track Y cursor descending from top of canvas.
            let mut y = canvas_rect.min.y;

            // --- Track 1: Identity Histogram (Geneious Prime signature track) ---
            if app.display.show_identity_graph {
                draw_gutter_label(&painter, canvas_rect.min.x, y, IDENTITY_HEIGHT, "Identity", pal.text_muted, pal.border);
                for col in first_col..last_col {
                    let id = aln.columns()[col].identity_pct;
                    let bar_h = (id / 100.0) * (IDENTITY_HEIGHT - 6.0);
                    let x0 = col_to_x(col);
                    let color = identity_bar_color(id, pal);
                    data_painter.rect_filled(
                        Rect::from_min_max(
                            Pos2::new(x0 + 0.5, y + IDENTITY_HEIGHT - 3.0 - bar_h),
                            Pos2::new(x0 + col_w - 0.5, y + IDENTITY_HEIGHT - 3.0),
                        ),
                        1.0,
                        color,
                    );
                }
                painter.line_segment(
                    [Pos2::new(canvas_rect.min.x, y + IDENTITY_HEIGHT), Pos2::new(canvas_rect.max.x, y + IDENTITY_HEIGHT)],
                    Stroke::new(1.0, pal.border),
                );
                y += IDENTITY_HEIGHT;
            }

            // --- Track 2: Consensus Sequence ---
            let font_size = (col_w * 0.80).clamp(8.5, 14.5);
            let mono_font = FontId::monospace(font_size);

            if app.display.show_consensus {
                draw_gutter_label(&painter, canvas_rect.min.x, y, CONSENSUS_HEIGHT, "Consensus", pal.text_primary, pal.border);
                for col in first_col..last_col {
                    let c_byte = aln.columns()[col].consensus;
                    let x0 = col_to_x(col);
                    let cell = Rect::from_min_max(Pos2::new(x0, y + 1.0), Pos2::new(x0 + col_w, y + CONSENSUS_HEIGHT - 1.0));
                    let (bg, fg) = residue_colors(
                        c_byte,
                        c_byte,
                        aln.alphabet,
                        app.display.color_scheme,
                        crate::theme::HighlightMode::AllResidues,
                        app.theme,
                    );
                    if let Some(bg_col) = bg {
                        data_painter.rect_filled(cell, 1.0, bg_col);
                    }
                    if col_w >= 8.0 {
                        let ch_str = [c_byte];
                        let s = std::str::from_utf8(&ch_str).unwrap_or("?");
                        data_painter.text(cell.center(), Align2::CENTER_CENTER, s, mono_font.clone(), fg);
                    }
                }
                painter.line_segment(
                    [Pos2::new(canvas_rect.min.x, y + CONSENSUS_HEIGHT), Pos2::new(canvas_rect.max.x, y + CONSENSUS_HEIGHT)],
                    Stroke::new(1.0, pal.border),
                );
                y += CONSENSUS_HEIGHT;
            }

            // --- Track 3: Coordinate Ruler ---
            draw_gutter_label(&painter, canvas_rect.min.x, y, RULER_HEIGHT, "Position", pal.text_muted, pal.border);
            let ruler_step: usize = if col_w >= 11.0 {
                10
            } else if col_w >= 6.5 {
                20
            } else {
                50
            };
            for col in first_col..last_col {
                let pos_1 = col + 1;
                let cx = col_to_x(col) + col_w * 0.5;
                if pos_1 == 1 || pos_1 % ruler_step == 0 {
                    data_painter.line_segment(
                        [Pos2::new(cx, y + RULER_HEIGHT - 6.0), Pos2::new(cx, y + RULER_HEIGHT)],
                        Stroke::new(1.0, pal.text_muted),
                    );
                    data_painter.text(
                        Pos2::new(cx, y + 7.0),
                        Align2::CENTER_CENTER,
                        pos_1.to_string(),
                        FontId::monospace(9.5),
                        pal.text_muted,
                    );
                } else if pos_1 % (ruler_step / 2).max(1) == 0 {
                    data_painter.line_segment(
                        [Pos2::new(cx, y + RULER_HEIGHT - 3.5), Pos2::new(cx, y + RULER_HEIGHT)],
                        Stroke::new(1.0, pal.border),
                    );
                }
            }
            painter.line_segment(
                [Pos2::new(canvas_rect.min.x, y + RULER_HEIGHT), Pos2::new(canvas_rect.max.x, y + RULER_HEIGHT)],
                Stroke::new(1.0, pal.border),
            );
            y += RULER_HEIGHT;

            // --- Track 3b: Optional Restriction Enzyme Cleavage Sites Track ---
            if show_restr {
                draw_gutter_label(
                    &painter,
                    canvas_rect.min.x,
                    y,
                    RESTRICTION_HEIGHT,
                    "Restriction Sites",
                    Color32::from_rgb(192, 68, 128),
                    pal.border,
                );
                let mut last_label_x = f32::NEG_INFINITY;
                for hit in restriction_sites {
                    if hit.cut_col < first_col || hit.cut_col >= last_col {
                        continue;
                    }
                    let cx = col_to_x(hit.cut_col);
                    let site_col = Color32::from_rgb(196, 62, 124);
                    data_painter.line_segment(
                        [Pos2::new(cx, y + RESTRICTION_HEIGHT - 8.0), Pos2::new(cx, y + RESTRICTION_HEIGHT)],
                        Stroke::new(1.8, site_col),
                    );
                    if cx > last_label_x + 36.0 {
                        data_painter.text(
                            Pos2::new(cx + 2.0, y + 7.0),
                            Align2::LEFT_CENTER,
                            hit.enzyme_name,
                            FontId::proportional(10.0),
                            site_col,
                        );
                        last_label_x = cx;
                    }
                }
                painter.line_segment(
                    [Pos2::new(canvas_rect.min.x, y + RESTRICTION_HEIGHT), Pos2::new(canvas_rect.max.x, y + RESTRICTION_HEIGHT)],
                    Stroke::new(1.0, pal.border),
                );
                y += RESTRICTION_HEIGHT;
            }

            // --- Track 4: Packed Annotation Arrow Tracks ---
            if ann_lanes > 0 {
                draw_gutter_label(
                    &painter,
                    canvas_rect.min.x,
                    y,
                    ann_h,
                    "Annotations",
                    pal.accent,
                    pal.border,
                );

                for (ann_idx, ann) in aln.annotations().iter().enumerate() {
                    if ann.end <= first_col || ann.start >= last_col {
                        continue;
                    }
                    let lane = ann_lane_map[ann_idx];
                    let lane_y = y + lane as f32 * ANN_LANE_HEIGHT + 2.0;
                    let bar_h = ANN_LANE_HEIGHT - 4.0;
                    let x0 = col_to_x(ann.start);
                    let x1 = col_to_x(ann.end).max(x0 + 8.0);
                    let ann_rect = Rect::from_min_max(Pos2::new(x0, lane_y), Pos2::new(x1, lane_y + bar_h));

                    let [r, g, b] = ann.rgb();
                    let fill = Color32::from_rgb(r, g, b);
                    let is_selected = app.selected_annotation_id == Some(ann.id);
                    let stroke = if is_selected {
                        Stroke::new(2.0, pal.text_primary)
                    } else {
                        Stroke::new(1.0, Color32::from_black_alpha(60))
                    };

                    // Paint directional arrow polygon (Geneious Prime feature arrow style).
                    let tip = (bar_h * 0.45).min((ann_rect.width() * 0.35).max(0.0));
                    let pts = match ann.strand {
                        Strand::Forward => vec![
                            Pos2::new(ann_rect.min.x, ann_rect.min.y),
                            Pos2::new(ann_rect.max.x - tip, ann_rect.min.y),
                            Pos2::new(ann_rect.max.x, ann_rect.center().y),
                            Pos2::new(ann_rect.max.x - tip, ann_rect.max.y),
                            Pos2::new(ann_rect.min.x, ann_rect.max.y),
                        ],
                        Strand::Reverse => vec![
                            Pos2::new(ann_rect.min.x + tip, ann_rect.min.y),
                            Pos2::new(ann_rect.max.x, ann_rect.min.y),
                            Pos2::new(ann_rect.max.x, ann_rect.max.y),
                            Pos2::new(ann_rect.min.x + tip, ann_rect.max.y),
                            Pos2::new(ann_rect.min.x, ann_rect.center().y),
                        ],
                    };
                    data_painter.add(Shape::convex_polygon(pts, fill, stroke));

                    // Clip label to visible portion of annotation arrow.
                    let vis_x0 = ann_rect.min.x.max(data_left + 4.0);
                    let vis_x1 = ann_rect.max.x.min(canvas_rect.max.x - 4.0);
                    if vis_x1 - vis_x0 > 24.0 {
                        let label_text = if ann.uniprot.is_some() {
                            format!("{} [3D]", ann.name)
                        } else {
                            ann.name.clone()
                        };
                        let text_col = if (u16::from(r) * 299 + u16::from(g) * 587 + u16::from(b) * 114) / 1000 > 155 {
                            Color32::from_rgb(20, 26, 36)
                        } else {
                            Color32::WHITE
                        };
                        let clipped = data_painter.with_clip_rect(Rect::from_min_max(
                            Pos2::new(vis_x0, ann_rect.min.y),
                            Pos2::new(vis_x1, ann_rect.max.y),
                        ));
                        clipped.text(
                            Pos2::new(vis_x0 + 5.0, ann_rect.center().y),
                            Align2::LEFT_CENTER,
                            label_text,
                            FontId::proportional(11.0),
                            text_col,
                        );
                    }

                    // Hit test single-click (select) & double-click (launch AlphaFold 3D window!).
                    if let Some(ptr) = response.interact_pointer_pos() {
                        if ann_rect.contains(ptr) {
                            if response.double_clicked() {
                                clicked_ann = Some((ann.id, ann.start, ann.end));
                                if let Some(up) = &ann.uniprot {
                                    double_clicked_ann_3d = Some((ann.id, up.accession.clone()));
                                }
                            } else if response.clicked() {
                                clicked_ann = Some((ann.id, ann.start, ann.end));
                            }
                        }
                    }
                }

                painter.line_segment(
                    [Pos2::new(canvas_rect.min.x, y + ann_h), Pos2::new(canvas_rect.max.x, y + ann_h)],
                    Stroke::new(1.0, pal.border),
                );
                y += ann_h;
            }

            // --- Track 5: Optional Nucleotide CDS Translation Track ---
            if show_trans {
                let frame_lbl = format!("Translation (+{})", app.display.translation_frame + 1);
                draw_gutter_label(&painter, canvas_rect.min.x, y, TRANSLATION_HEIGHT, &frame_lbl, pal.text_muted, pal.border);
                if let Some(ref_seq) = aln.sequences().first() {
                    let codons = translate_aligned_region(
                        &ref_seq.residues,
                        0,
                        ref_seq.len(),
                        Strand::Forward,
                        app.display.translation_frame,
                    );
                    for codon in codons {
                        if codon.end_col <= first_col || codon.start_col >= last_col {
                            continue;
                        }
                        let x0 = col_to_x(codon.start_col);
                        let x1 = col_to_x(codon.end_col);
                        let c_rect = Rect::from_min_max(
                            Pos2::new(x0 + 1.0, y + 2.0),
                            Pos2::new(x1 - 1.0, y + TRANSLATION_HEIGHT - 2.0),
                        );
                        let (bg, fg) = residue_colors(
                            codon.amino_acid,
                            codon.amino_acid,
                            helix_core::Alphabet::Protein,
                            app.display.color_scheme,
                            crate::theme::HighlightMode::AllResidues,
                            app.theme,
                        );
                        if let Some(bg_c) = bg {
                            data_painter.rect_filled(c_rect, 2.0, bg_c);
                        }
                        if c_rect.width() >= 9.0 {
                            let b = [codon.amino_acid];
                            let s = std::str::from_utf8(&b).unwrap_or("X");
                            data_painter.text(c_rect.center(), Align2::CENTER_CENTER, s, FontId::monospace(10.5), fg);
                        }
                    }
                }
                painter.line_segment(
                    [Pos2::new(canvas_rect.min.x, y + TRANSLATION_HEIGHT), Pos2::new(canvas_rect.max.x, y + TRANSLATION_HEIGHT)],
                    Stroke::new(1.0, pal.border),
                );
                y += TRANSLATION_HEIGHT;
            }

            // --- Track 6: Aligned Sequence Rows (with optional Complement 3'->5' strand) ---
            let seq_start_y = y;
            for (row_idx, seq) in aln.sequences().iter().enumerate() {
                let ry = seq_start_y + row_idx as f32 * effective_row_h;

                // Subtle alternating row shading.
                if row_idx % 2 == 1 {
                    painter.rect_filled(
                        Rect::from_min_max(Pos2::new(canvas_rect.min.x, ry), Pos2::new(canvas_rect.max.x, ry + effective_row_h)),
                        0.0,
                        pal.row_alt_bg,
                    );
                }

                // Gutter sequence label + UniProt tag.
                let gutter_cell = Rect::from_min_max(
                    Pos2::new(canvas_rect.min.x, ry),
                    Pos2::new(data_left, ry + row_h),
                );
                let gutter_clipped = painter.with_clip_rect(gutter_cell);
                let row_label = if seq.uniprot_accession.is_some() {
                    format!("{}  *", seq.name)
                } else {
                    seq.name.clone()
                };
                gutter_clipped.text(
                    Pos2::new(gutter_cell.min.x + 10.0, gutter_cell.center().y),
                    Align2::LEFT_CENTER,
                    row_label,
                    FontId::proportional(11.5),
                    pal.text_primary,
                );
                if show_comp {
                    painter.text(
                        Pos2::new(canvas_rect.min.x + 18.0, ry + row_h + row_h * 0.35),
                        Align2::LEFT_CENTER,
                        "3'->5' Comp",
                        FontId::proportional(9.5),
                        pal.text_muted,
                    );
                }

                // Draw visible columns for this sequence row.
                for col in first_col..last_col {
                    let residue = seq.residues[col];
                    let consensus = aln.columns()[col].consensus;
                    let x0 = col_to_x(col);
                    let cell = Rect::from_min_max(
                        Pos2::new(x0, ry + 1.0),
                        Pos2::new(x0 + col_w, ry + row_h - 1.0),
                    );

                    let (bg, fg) = residue_colors(
                        residue,
                        consensus,
                        aln.alphabet,
                        app.display.color_scheme,
                        app.display.highlight_mode,
                        app.theme,
                    );
                    if let Some(bg_col) = bg {
                        data_painter.rect_filled(cell, 1.0, bg_col);
                    }

                    if col_w >= 7.5 {
                        let agrees_with_consensus =
                            !is_gap(residue) && !is_gap(consensus) && residue == consensus;
                        let draw_byte = if app.display.dots_for_matches
                            && row_idx > 0
                            && agrees_with_consensus
                        {
                            b'.'
                        } else {
                            residue
                        };
                        let buf = [draw_byte];
                        let s = std::str::from_utf8(&buf).unwrap_or("?");
                        data_painter.text(cell.center(), Align2::CENTER_CENTER, s, mono_font.clone(), fg);

                        if show_comp {
                            let comp_byte = [complement_nucleotide(residue)];
                            let comp_s = std::str::from_utf8(&comp_byte).unwrap_or(" ");
                            data_painter.text(
                                Pos2::new(cell.center().x, ry + row_h + row_h * 0.35),
                                Align2::CENTER_CENTER,
                                comp_s,
                                FontId::monospace((font_size * 0.85).max(8.0)),
                                pal.text_muted,
                            );
                        }
                    }
                }
            }

            // --- Overlay: Motif search hits ---
            for (hit_idx, hit) in app.motif_hits.iter().enumerate() {
                if hit.end_col <= first_col || hit.start_col >= last_col {
                    continue;
                }
                let hx0 = col_to_x(hit.start_col);
                let hx1 = col_to_x(hit.end_col);
                let hy0 = seq_start_y + hit.row as f32 * effective_row_h;
                let h_rect = Rect::from_min_max(Pos2::new(hx0, hy0), Pos2::new(hx1, hy0 + row_h));
                let stroke_w = if hit_idx == app.active_motif_hit { 2.4 } else { 1.5 };
                data_painter.rect_stroke(
                    h_rect,
                    2.0,
                    Stroke::new(stroke_w, pal.motif_highlight),
                );
            }

            // --- Overlay: Active Base Cursor Indicator (like `1,179` vertical cursor in screenshot) ---
            if app.cursor_col >= first_col && app.cursor_col < last_col {
                let cx = col_to_x(app.cursor_col);
                let cursor_col_color = if app.display.allow_editing {
                    Color32::from_rgb(224, 96, 36)
                } else {
                    pal.accent
                };
                data_painter.line_segment(
                    [Pos2::new(cx, seq_start_y - 4.0), Pos2::new(cx, seq_start_y + sequences_h)],
                    Stroke::new(2.0, cursor_col_color),
                );
            }

            // --- Overlay: Selected column range ---
            if let Some((sel_s, sel_e)) = new_selected_cols {
                if sel_e > first_col && sel_s < last_col {
                    let sx0 = col_to_x(sel_s);
                    let sx1 = col_to_x(sel_e);
                    let sel_rect = Rect::from_min_max(
                        Pos2::new(sx0, canvas_rect.min.y),
                        Pos2::new(sx1, seq_start_y + sequences_h),
                    );
                    data_painter.rect_filled(sel_rect, 0.0, pal.selection_fill);
                    data_painter.rect_stroke(
                        sel_rect,
                        0.0,
                        Stroke::new(1.2, pal.selection_stroke),
                    );
                }
            }

            // --- Overlay: Bidirectional 3D Protein Sync Highlight Column ---
            if let Some((s3d, e3d)) = synced_3d_cols {
                if e3d > first_col && s3d < last_col {
                    let sx0 = col_to_x(s3d);
                    let sx1 = col_to_x(e3d);
                    let sync_rect = Rect::from_min_max(
                        Pos2::new(sx0, canvas_rect.min.y),
                        Pos2::new(sx1, seq_start_y + sequences_h),
                    );
                    data_painter.rect_stroke(
                        sync_rect,
                        1.0,
                        Stroke::new(2.2, pal.sync_3d_highlight),
                    );
                }
            }

            // --- Overlay: Hover column guide ---
            if let Some(hcol) = new_hover_col {
                let hx0 = col_to_x(hcol);
                let hover_rect = Rect::from_min_max(
                    Pos2::new(hx0, canvas_rect.min.y),
                    Pos2::new(hx0 + col_w, seq_start_y + sequences_h),
                );
                data_painter.rect_filled(hover_rect, 0.0, pal.hover_col_fill);
            }

            // Crisp gutter vertical border.
            painter.line_segment(
                [Pos2::new(data_left, canvas_rect.min.y), Pos2::new(data_left, canvas_rect.max.y)],
                Stroke::new(1.0, pal.border),
            );
        });

    // --- Custom Horizontal Scrollbar at bottom of alignment canvas ---
    let (sb_rect, sb_resp) =
        ui.allocate_exact_size(Vec2::new(avail_rect.width(), scrollbar_h), Sense::click_and_drag());
    let sb_painter = ui.painter_at(sb_rect);
    sb_painter.rect_filled(sb_rect, 0.0, pal.panel_bg);
    sb_painter.line_segment(
        [sb_rect.left_top(), sb_rect.right_top()],
        Stroke::new(1.0, pal.border),
    );

    let track = Rect::from_min_max(
        Pos2::new(sb_rect.min.x + GUTTER_WIDTH + 4.0, sb_rect.min.y + 3.0),
        Pos2::new(sb_rect.max.x - 4.0, sb_rect.max.y - 3.0),
    );
    if track.width() > 16.0 {
        sb_painter.rect_filled(track, 4.0, pal.header_bg);
        let thumb_frac = (data_viewport_w / total_seq_w).clamp(0.06, 1.0);
        let thumb_w = (track.width() * thumb_frac).max(24.0);
        let travel = (track.width() - thumb_w).max(1.0);
        let pos_frac = if max_h_scroll > 0.0 {
            app.horizontal_scroll_px / max_h_scroll
        } else {
            0.0
        };
        let tx0 = track.min.x + pos_frac * travel;
        let thumb = Rect::from_min_max(Pos2::new(tx0, track.min.y), Pos2::new(tx0 + thumb_w, track.max.y));
        let thumb_col = if sb_resp.hovered() || sb_resp.dragged() {
            pal.accent
        } else {
            pal.text_muted
        };
        sb_painter.rect_filled(thumb, 4.0, thumb_col);

        if sb_resp.dragged() || sb_resp.clicked() {
            if let Some(ptr) = sb_resp.interact_pointer_pos() {
                let frac = ((ptr.x - track.min.x - thumb_w * 0.5) / travel).clamp(0.0, 1.0);
                app.horizontal_scroll_px = frac * max_h_scroll;
            }
        }
    }

    // Apply deferred interaction updates cleanly.
    app.hovered_column = new_hover_col;
    app.hovered_row = new_hover_row;
    app.selected_cols = new_selected_cols;

    if let Some((ann_id, start, end)) = clicked_ann {
        app.selected_annotation_id = Some(ann_id);
        app.selected_cols = Some((start, end));
        app.show_inspector_panel = true;
        app.inspector_tab = InspectorTab::Annotations;
        app.rebuild_3d_alignment_map();
    }

    if let Some((ann_id, acc)) = double_clicked_ann_3d {
        app.open_3d_for_accession(ctx, &acc, Some(ann_id));
    }

    // Clicking a column in the alignment updates the active cursor & highlights the mapped residue in 3D!
    if let Some(col) = clicked_col_sync_3d {
        app.cursor_col = col;
        if let Some(r) = new_hover_row {
            app.cursor_row = r;
        }
        if let Some(ref map) = app.viewer3d.alignment_map {
            if let Some(res_idx) = map.residue_for_column(col) {
                app.viewer3d.selected_residue = Some(res_idx);
            }
        }
    }
}

fn draw_gutter_label(
    painter: &egui::Painter,
    left_x: f32,
    top_y: f32,
    height: f32,
    label: &str,
    color: Color32,
    border: Color32,
) {
    painter.text(
        Pos2::new(left_x + 12.0, top_y + height * 0.5),
        Align2::LEFT_CENTER,
        label,
        FontId::proportional(11.5),
        color,
    );
    painter.line_segment(
        [Pos2::new(left_x, top_y + height), Pos2::new(left_x + GUTTER_WIDTH, top_y + height)],
        Stroke::new(0.5, border),
    );
}

fn identity_bar_color(identity_pct: f32, pal: &Palette) -> Color32 {
    if identity_pct >= 99.9 {
        pal.identity_high
    } else if identity_pct >= 30.0 {
        pal.identity_mid
    } else {
        pal.identity_low
    }
}
