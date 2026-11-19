//! Interactive Circular Phage / Plasmid Genome Visualization widget.
//!
//! Renders the active alignment as a concentric circular genome map inspired by Geneious Prime
//! and CGView: curved directional annotation arc-arrows with outer callout labels & `[3D]` badges,
//! coordinate ruler, multi-strain aligned sequence rings highlighting variants/indels, pairwise
//! identity ring, circular GC-content/GC-skew ring, and bidirectional 3D protein sync.

use crate::app::{HelixApp, InspectorTab};
use crate::theme::{residue_colors, HighlightMode, Palette};
use egui::{
    Align2, Color32, Context, FontId, Frame, Pos2, RichText, Sense, Shape, Stroke, Ui, Vec2,
};
use helix_core::sequence::is_gap;
use helix_core::{pack_annotation_lanes, Strand};
use std::f32::consts::{FRAC_PI_2, TAU};

/// Renders the circular phage/genome header bar and interactive radial canvas.
pub fn show_circular_view(app: &mut HelixApp, ctx: &Context, ui: &mut Ui) {
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

    // 1. Top document header + circular controls strip.
    show_circular_header(app, ctx, ui, &pal);

    // 2. Main interactive circular map canvas.
    show_circular_canvas(app, ctx, ui, &pal);
}

fn show_circular_header(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let mut launch_3d_req: Option<(String, Option<u64>)> = None;

    Frame::none()
        .fill(pal.panel_bg)
        .stroke(Stroke::new(1.0, pal.border))
        .inner_margin(egui::Margin::symmetric(12.0, 6.0))
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                if let Some(aln) = app.current_alignment() {
                    Frame::none()
                        .fill(Color32::from_rgb(16, 138, 134))
                        .rounding(egui::Rounding::same(4.0))
                        .inner_margin(egui::Margin::symmetric(6.0, 2.0))
                        .show(ui, |ui| {
                            ui.label(
                                RichText::new("CIRCULAR PHAGE")
                                    .strong()
                                    .size(10.0)
                                    .color(Color32::WHITE),
                            );
                        });

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
                    let mut stats_txt = format!(
                        "{} seqs  |  {} {unit} circular  |  {:.1}% id",
                        aln.row_count(),
                        aln.width(),
                        aln.mean_pairwise_identity_pct()
                    );
                    if aln.alphabet.is_nucleotide() {
                        if let Some(gc) = aln.sequences().first().and_then(|s| s.gc_content_pct()) {
                            stats_txt.push_str(&format!("  |  GC: {gc:.1}%"));
                        }
                    }
                    ui.label(RichText::new(stats_txt).small().color(pal.text_muted));
                }

                ui.separator();
                ui.checkbox(&mut app.circular.show_labels, "Labels");
                ui.checkbox(&mut app.circular.show_sequence_rings, "Strain Rings");
                ui.checkbox(&mut app.circular.show_identity_ring, "Identity");
                if app
                    .current_alignment()
                    .is_some_and(|a| a.alphabet.is_nucleotide())
                {
                    ui.checkbox(&mut app.circular.show_gc_ring, "GC Content");
                }
                if ui.small_button("Reset Circle").clicked() {
                    app.circular.rotation_rad = 0.0;
                    app.circular.zoom = 1.0;
                }

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
                        }
                    }
                });
            });
        });

    if let Some((acc, ann_id)) = launch_3d_req {
        app.open_3d_for_accession(ctx, &acc, ann_id);
    }
}

fn show_circular_canvas(app: &mut HelixApp, ctx: &Context, ui: &mut Ui, pal: &Palette) {
    let avail_rect = ui.available_rect_before_wrap();
    let (rect, response) = ui.allocate_exact_size(avail_rect.size(), Sense::click_and_drag());
    let painter = ui.painter_at(rect);

    painter.rect_filled(rect, 0.0, pal.canvas_bg);

    // Handle drag-to-rotate, scroll-to-zoom, and double-click-empty reset BEFORE borrowing `aln`.
    let center = rect.center();
    if response.dragged() {
        let delta = response.drag_delta();
        if let Some(ptr) = response.interact_pointer_pos() {
            // Compute tangential angular drag around the circle center.
            let prev = ptr - delta;
            let a0 = (prev.y - center.y).atan2(prev.x - center.x);
            let a1 = (ptr.y - center.y).atan2(ptr.x - center.x);
            let mut da = a1 - a0;
            if da > std::f32::consts::PI {
                da -= TAU;
            } else if da < -std::f32::consts::PI {
                da += TAU;
            }
            app.circular.rotation_rad = (app.circular.rotation_rad + da).rem_euclid(TAU);
        }
    }

    if response.hovered() {
        let scroll_y = ctx.input(|i| i.smooth_scroll_delta.y);
        if scroll_y.abs() > 0.1 {
            let factor = (1.0 + scroll_y * 0.0025).clamp(0.8, 1.25);
            app.circular.zoom = (app.circular.zoom * factor).clamp(0.55, 2.4);
        }
    }

    let Some(doc) = app.docs.get(app.active_doc) else {
        return;
    };
    let aln = &doc.alignment;
    let width_cols = aln.width().max(1);

    let min_dim = rect.width().min(rect.height()).max(180.0);
    let base_r = min_dim * 0.31 * app.circular.zoom;
    let rot = app.circular.rotation_rad;

    // Helper mapping a 0-based column coordinate (floating) to angle in radians:
    // Column 0 starts at 12 o'clock (-PI/2) and proceeds clockwise.
    let col_to_angle = |col_f: f32| -> f32 { -FRAC_PI_2 + rot + (col_f / width_cols as f32) * TAU };
    let polar = |r: f32, angle: f32| -> Pos2 {
        Pos2::new(center.x + r * angle.cos(), center.y + r * angle.sin())
    };

    // Pack annotations into non-overlapping radial lanes.
    let (ann_lanes, ann_lane_map) = if app.display.show_annotations {
        pack_annotation_lanes(aln.annotations(), 1)
    } else {
        (0, Vec::new())
    };
    let lane_thick = (base_r * 0.068).clamp(11.0, 20.0);
    let lane_gap = 4.0f32;
    let ann_base_r = base_r * 1.03;

    // Hit-test pointer angle and radius against annotations and genome coordinates.
    let mut hovered_col: Option<usize> = None;
    let mut hovered_ann_idx: Option<usize> = None;

    if let Some(ptr) = response.hover_pos() {
        let dx = ptr.x - center.x;
        let dy = ptr.y - center.y;
        let dist = (dx * dx + dy * dy).sqrt();
        let angle = dy.atan2(dx);
        let rel_angle = (angle - (-FRAC_PI_2 + rot)).rem_euclid(TAU);
        let col = ((rel_angle / TAU) * width_cols as f32).floor() as usize;
        let clamped_col = col.min(width_cols.saturating_sub(1));

        // Check annotation rings first.
        if ann_lanes > 0 {
            for (idx, ann) in aln.annotations().iter().enumerate() {
                let lane = ann_lane_map[idx];
                let r_in = ann_base_r + lane as f32 * (lane_thick + lane_gap);
                let r_out = r_in + lane_thick;
                if dist >= r_in - 2.0 && dist <= r_out + 2.0 && ann.contains_column(clamped_col) {
                    hovered_ann_idx = Some(idx);
                    break;
                }
            }
        }

        // Track genome coordinate whenever hovering within the circular map annulus.
        let max_ring_r = ann_base_r + ann_lanes.max(1) as f32 * (lane_thick + lane_gap) + 8.0;
        if dist >= base_r * 0.25 && dist <= max_ring_r {
            hovered_col = Some(clamped_col);
        }
    }

    // --- Layer 1: Subtle background circular guide grid ---
    let guide_stroke = Stroke::new(1.0, pal.border);
    painter.circle_stroke(center, base_r * 0.96, guide_stroke);

    // --- Layer 2: Coordinate Ruler Ring (ticks + labels) ---
    let ruler_r = base_r * 0.96;
    let major_ticks = 12usize;
    let minor_ticks = 48usize;
    for t in 0..minor_ticks {
        let frac = t as f32 / minor_ticks as f32;
        let a = -FRAC_PI_2 + rot + frac * TAU;
        let is_major = t % (minor_ticks / major_ticks) == 0;
        let tick_len = if is_major { 7.0 } else { 3.5 };
        let stroke = if is_major {
            Stroke::new(1.3, pal.text_muted)
        } else {
            Stroke::new(0.8, pal.border)
        };
        painter.line_segment([polar(ruler_r - tick_len, a), polar(ruler_r, a)], stroke);

        if is_major {
            let coord_1 = ((frac * width_cols as f32).round() as usize).max(1);
            let lbl_pos = polar(ruler_r - 15.0, a);
            painter.text(
                lbl_pos,
                Align2::CENTER_CENTER,
                coord_1.to_string(),
                FontId::monospace(9.5),
                pal.text_muted,
            );
        }
    }

    // --- Layer 3: Concentric Aligned Sequence Rings (Phage Isolate Variant Comparison) ---
    let mut inner_cursor_r = base_r * 0.87;
    if app.circular.show_sequence_rings && aln.row_count() > 0 {
        let max_rows = aln.row_count().min(8);
        let seq_band_total = base_r * 0.24;
        let ring_step = seq_band_total / max_rows as f32;
        let ring_thick = (ring_step * 0.78).clamp(3.0, 12.0);
        let bins = width_cols.clamp(60, 240);

        for (row_idx, seq) in aln.sequences().iter().take(max_rows).enumerate() {
            let r_out = inner_cursor_r - row_idx as f32 * ring_step;
            let r_in = (r_out - ring_thick).max(base_r * 0.28);

            // Base track ring for this isolate.
            let base_track_col = if row_idx == 0 {
                pal.header_bg
            } else if row_idx % 2 == 1 {
                pal.panel_bg
            } else {
                pal.row_alt_bg
            };
            draw_annular_ring(
                &painter,
                center,
                r_in,
                r_out,
                base_track_col,
                Stroke::new(0.5, pal.border),
            );

            // Paint per-bin residue / variant wedges around the ring.
            for b in 0..bins {
                let c0 = (b * width_cols) / bins;
                let c1 = (((b + 1) * width_cols) / bins).max(c0 + 1).min(width_cols);
                let mut mismatch_byte: Option<u8> = None;
                let mut gap_found = false;

                for c in c0..c1 {
                    let res = seq.residues[c];
                    let cons = aln.columns()[c].consensus;
                    if is_gap(res) {
                        gap_found = true;
                    } else if !is_gap(cons) && res != cons {
                        mismatch_byte = Some(res);
                        break;
                    }
                }

                let wedge_col = if let Some(mm) = mismatch_byte {
                    let (bg, _) = residue_colors(
                        mm,
                        b' ',
                        aln.alphabet,
                        app.display.color_scheme,
                        HighlightMode::AllResidues,
                        app.theme,
                    );
                    bg.unwrap_or(pal.identity_low)
                } else if gap_found {
                    pal.text_muted
                } else if row_idx == 0 && app.display.highlight_mode == HighlightMode::AllResidues {
                    // Reference track subtle tint.
                    Color32::from_rgba_unmultiplied(pal.accent.r(), pal.accent.g(), pal.accent.b(), 45)
                } else {
                    continue;
                };

                let a0 = col_to_angle(c0 as f32);
                let a1 = col_to_angle(c1 as f32).max(a0 + 0.018);
                draw_arc_quad(&painter, center, r_in, r_out, a0, a1, wedge_col);
            }
        }
        inner_cursor_r -= seq_band_total + base_r * 0.03;
    }

    // --- Layer 4: Pairwise Conservation / Identity Ring ---
    if app.circular.show_identity_ring {
        let id_r_out = inner_cursor_r;
        let id_band = base_r * 0.12;
        let id_r_in = (id_r_out - id_band).max(base_r * 0.22);
        let bins = width_cols.clamp(60, 220);

        painter.circle_stroke(center, id_r_in, Stroke::new(0.6, pal.border));
        for b in 0..bins {
            let c0 = (b * width_cols) / bins;
            let c1 = (((b + 1) * width_cols) / bins).max(c0 + 1).min(width_cols);
            let avg_id = aln.columns()[c0..c1]
                .iter()
                .map(|c| c.identity_pct)
                .sum::<f32>()
                / (c1 - c0) as f32;
            let bar_r_out = id_r_in + (avg_id / 100.0) * (id_r_out - id_r_in);
            let bar_col = if avg_id >= 99.9 {
                pal.identity_high
            } else if avg_id >= 30.0 {
                pal.identity_mid
            } else {
                pal.identity_low
            };
            let a0 = col_to_angle(c0 as f32);
            let a1 = col_to_angle(c1 as f32);
            draw_arc_quad(&painter, center, id_r_in, bar_r_out, a0, a1, bar_col);
        }
        inner_cursor_r = id_r_in - base_r * 0.03;
    }

    // --- Layer 5: Circular GC-Content & GC-Skew Ring (Nucleotide / Phage Genomes) ---
    if app.circular.show_gc_ring && aln.alphabet.is_nucleotide() {
        if let Some(gc_prof) = aln.circular_gc_profile(144) {
            let gc_band = base_r * 0.13;
            let gc_mid = (inner_cursor_r - gc_band * 0.5).max(base_r * 0.20);
            let gc_half = gc_band * 0.48;

            painter.circle_stroke(center, gc_mid, Stroke::new(0.8, pal.border));
            let n = gc_prof.len();
            for (i, &(gc_frac, _gc_skew)) in gc_prof.iter().enumerate() {
                let a0 = -FRAC_PI_2 + rot + (i as f32 / n as f32) * TAU;
                let a1 = -FRAC_PI_2 + rot + ((i + 1) as f32 / n as f32) * TAU;
                // Plot deviation from 50% GC: GC-rich (>50%) outward in teal, AT-rich (<50%) inward in amber/violet.
                let delta = ((gc_frac - 0.5) * 2.4).clamp(-1.0, 1.0);
                let (r0, r1, col) = if delta >= 0.0 {
                    (
                        gc_mid,
                        gc_mid + delta * gc_half,
                        Color32::from_rgb(38, 166, 154),
                    )
                } else {
                    (
                        gc_mid + delta * gc_half,
                        gc_mid,
                        Color32::from_rgb(171, 104, 212),
                    )
                };
                draw_arc_quad(&painter, center, r0, r1, a0, a1, col);
            }
            inner_cursor_r = (gc_mid - gc_half - base_r * 0.03).max(base_r * 0.18);
        }
    }

    // --- Layer 6: Outer Packed Directional Annotation Arc-Arrows + Callout Labels ---
    let mut clicked_ann: Option<(u64, usize, usize)> = None;
    let mut double_clicked_ann_3d: Option<(u64, String)> = None;

    if ann_lanes > 0 {
        let label_ring_r = ann_base_r + ann_lanes as f32 * (lane_thick + lane_gap) + 18.0;

        for (ann_idx, ann) in aln.annotations().iter().enumerate() {
            let lane = ann_lane_map[ann_idx];
            let r_in = ann_base_r + lane as f32 * (lane_thick + lane_gap);
            let r_out = r_in + lane_thick;

            let a0 = col_to_angle(ann.start as f32);
            let mut a1 = col_to_angle(ann.end as f32);
            // Guarantee even 1-column variant sites are clearly visible as crisp wedges.
            if a1 - a0 < 0.035 {
                let mid = (a0 + a1) * 0.5;
                a1 = mid + 0.0175;
            }
            let a_start = a0.min(a1 - 0.025);

            let [r, g, b] = ann.rgb();
            let fill = Color32::from_rgb(r, g, b);
            let is_selected = app.selected_annotation_id == Some(ann.id);
            let is_hovered = hovered_ann_idx == Some(ann_idx);
            let stroke = if is_selected {
                Stroke::new(2.2, pal.text_primary)
            } else if is_hovered {
                Stroke::new(1.8, pal.accent)
            } else {
                Stroke::new(0.8, Color32::from_black_alpha(70))
            };

            draw_curved_annotation_arrow(
                &painter,
                center,
                r_in,
                r_out,
                a_start,
                a1,
                ann.strand,
                fill,
                stroke,
            );

            // Outer radial callout line + readable label with `[3D]` badge if linked to UniProt.
            if app.circular.show_labels {
                let mid_a = (a_start + a1) * 0.5;
                let anchor = polar(r_out + 2.0, mid_a);
                let elbow_r = label_ring_r + (ann_idx % 2) as f32 * 14.0;
                let elbow = polar(elbow_r, mid_a);
                let on_right = mid_a.cos() >= 0.0;
                let tip = Pos2::new(elbow.x + if on_right { 12.0 } else { -12.0 }, elbow.y);

                let line_col = if is_selected || is_hovered {
                    pal.accent
                } else {
                    pal.border
                };
                painter.line_segment([anchor, elbow], Stroke::new(1.0, line_col));
                painter.line_segment([elbow, tip], Stroke::new(1.0, line_col));

                let lbl = if ann.uniprot.is_some() {
                    format!("{} [3D]", ann.name)
                } else {
                    ann.name.clone()
                };
                let align = if on_right {
                    Align2::LEFT_CENTER
                } else {
                    Align2::RIGHT_CENTER
                };
                let text_pos = Pos2::new(tip.x + if on_right { 4.0 } else { -4.0 }, tip.y);
                let text_col = if is_selected || is_hovered {
                    pal.accent
                } else {
                    pal.text_primary
                };
                painter.text(text_pos, align, lbl, FontId::proportional(11.5), text_col);
            }

            if is_hovered {
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

    // --- Layer 7: Selected Range Wedge, 3D Sync Beam & Hover Crosshair ---
    let outer_highlight_r = ann_base_r + ann_lanes.max(1) as f32 * (lane_thick + lane_gap);
    let inner_highlight_r = (inner_cursor_r * 0.85).max(base_r * 0.22);

    if let Some((sel_s, sel_e)) = app.selected_cols {
        let a0 = col_to_angle(sel_s as f32);
        let a1 = col_to_angle(sel_e as f32).max(a0 + 0.02);
        draw_arc_quad(
            &painter,
            center,
            inner_highlight_r,
            outer_highlight_r,
            a0,
            a1,
            pal.selection_fill,
        );
    }

    // Bidirectional 3D Protein Window sync beam!
    let synced_3d_cols: Option<(usize, usize)> = app.viewer3d.alignment_map.as_ref().and_then(|map| {
        let active_res = app.viewer3d.hovered_residue.or(app.viewer3d.selected_residue)?;
        map.residue_to_cols.get(active_res).copied().flatten()
    });
    if let Some((s3d, e3d)) = synced_3d_cols {
        let mid_a = col_to_angle((s3d + e3d) as f32 * 0.5);
        painter.line_segment(
            [polar(inner_highlight_r, mid_a), polar(outer_highlight_r + 6.0, mid_a)],
            Stroke::new(2.5, pal.sync_3d_highlight),
        );
    }

    if let Some(hcol) = hovered_col {
        let ha = col_to_angle(hcol as f32 + 0.5);
        painter.line_segment(
            [polar(inner_highlight_r, ha), polar(outer_highlight_r, ha)],
            Stroke::new(1.2, pal.accent),
        );
    }

    // --- Layer 8: Center Genome Medallion ---
    let hub_r = (inner_cursor_r * 0.82).clamp(38.0, base_r * 0.34);
    painter.circle_filled(center, hub_r, pal.panel_bg);
    painter.circle_stroke(center, hub_r, Stroke::new(1.2, pal.border));

    let short_title = truncate_label(&aln.name, 22);
    painter.text(
        Pos2::new(center.x, center.y - 16.0),
        Align2::CENTER_CENTER,
        short_title,
        FontId::proportional(12.5),
        pal.text_primary,
    );
    let unit = if aln.alphabet.is_nucleotide() { "bp" } else { "aa" };
    painter.text(
        Pos2::new(center.x, center.y + 2.0),
        Align2::CENTER_CENTER,
        format!("{} {unit}  (Circular)", aln.width()),
        FontId::monospace(11.0),
        pal.accent,
    );
    painter.text(
        Pos2::new(center.x, center.y + 18.0),
        Align2::CENTER_CENTER,
        format!(
            "{} strains  |  {:.0}% ID",
            aln.row_count(),
            aln.mean_pairwise_identity_pct()
        ),
        FontId::proportional(10.0),
        pal.text_muted,
    );

    // Subtle circular legend & interaction hint in bottom corners.
    painter.text(
        Pos2::new(rect.min.x + 12.0, rect.max.y - 10.0),
        Align2::LEFT_BOTTOM,
        "Drag: Rotate Circle  |  Scroll: Zoom  |  Click: Select Coordinate  |  Double-Click [3D] Arc: Open AlphaFold 3D",
        FontId::proportional(11.0),
        pal.text_muted,
    );

    // Floating hover tooltip when hovering an annotation arc.
    if let (Some(ann_idx), Some(ptr)) = (hovered_ann_idx, response.hover_pos()) {
        if let Some(ann) = aln.annotations().get(ann_idx) {
            draw_annotation_tooltip(&painter, ptr, ann, pal);
        }
    }

    // Commit deferred interaction events cleanly after immutable `aln` borrow ends.
    let click_sync_col = if response.clicked() && hovered_ann_idx.is_none() {
        hovered_col
    } else {
        None
    };
    let dbl_click_reset = response.double_clicked() && hovered_ann_idx.is_none();

    app.hovered_column = hovered_col;

    if dbl_click_reset {
        app.circular.rotation_rad = 0.0;
        app.circular.zoom = 1.0;
    }

    if let Some((ann_id, start, end)) = clicked_ann {
        app.selected_annotation_id = Some(ann_id);
        app.selected_cols = Some((start, end));
        app.scroll_to_column = Some(start);
        app.show_inspector_panel = true;
        app.inspector_tab = InspectorTab::Annotations;
        app.rebuild_3d_alignment_map();
    }

    if let Some((ann_id, acc)) = double_clicked_ann_3d {
        app.open_3d_for_accession(ctx, &acc, Some(ann_id));
    }

    if let Some(col) = click_sync_col {
        app.selected_cols = Some((col, col + 1));
        app.scroll_to_column = Some(col);
        if let Some(ref map) = app.viewer3d.alignment_map {
            if let Some(res_idx) = map.residue_for_column(col) {
                app.viewer3d.selected_residue = Some(res_idx);
            }
        }
    }
}

/// Draws a filled 360-degree annular ring using smooth polygonal segments.
fn draw_annular_ring(
    painter: &egui::Painter,
    center: Pos2,
    r_in: f32,
    r_out: f32,
    fill: Color32,
    stroke: Stroke,
) {
    let segs = 96usize;
    for i in 0..segs {
        let a0 = (i as f32 / segs as f32) * TAU;
        let a1 = ((i + 1) as f32 / segs as f32) * TAU;
        draw_arc_quad(painter, center, r_in, r_out, a0, a1, fill);
    }
    painter.circle_stroke(center, r_in, stroke);
    painter.circle_stroke(center, r_out, stroke);
}

/// Draws a single annular arc wedge from angle `a0` to `a1` between radii `r_in` and `r_out`.
fn draw_arc_quad(
    painter: &egui::Painter,
    center: Pos2,
    r_in: f32,
    r_out: f32,
    a0: f32,
    a1: f32,
    fill: Color32,
) {
    if r_out <= r_in + 0.2 || a1 <= a0 {
        return;
    }
    let span = a1 - a0;
    let steps = ((span / 0.14).ceil() as usize).clamp(1, 32);
    for s in 0..steps {
        let t0 = a0 + (s as f32 / steps as f32) * span;
        let t1 = a0 + ((s + 1) as f32 / steps as f32) * span;
        let pts = vec![
            Pos2::new(center.x + r_in * t0.cos(), center.y + r_in * t0.sin()),
            Pos2::new(center.x + r_out * t0.cos(), center.y + r_out * t0.sin()),
            Pos2::new(center.x + r_out * t1.cos(), center.y + r_out * t1.sin()),
            Pos2::new(center.x + r_in * t1.cos(), center.y + r_in * t1.sin()),
        ];
        painter.add(Shape::convex_polygon(pts, fill, Stroke::NONE));
    }
}

/// Draws a curved directional annotation arrow along an annular arc (clockwise for Forward `+`,
/// counter-clockwise for Reverse `-`).
#[allow(clippy::too_many_arguments)]
fn draw_curved_annotation_arrow(
    painter: &egui::Painter,
    center: Pos2,
    r_in: f32,
    r_out: f32,
    a0: f32,
    a1: f32,
    strand: Strand,
    fill: Color32,
    stroke: Stroke,
) {
    let span = (a1 - a0).max(0.02);
    let r_mid = (r_in + r_out) * 0.5;
    let arrow_ang = (span * 0.28).min(0.11);

    let (body_a0, body_a1, tip_base_a, tip_a) = match strand {
        Strand::Forward => (a0, a1 - arrow_ang, a1 - arrow_ang, a1),
        Strand::Reverse => (a0 + arrow_ang, a1, a0 + arrow_ang, a0),
    };

    // Draw smooth curved arrow body segments.
    let body_span = (body_a1 - body_a0).max(0.0);
    let segs = ((body_span / 0.10).ceil() as usize).clamp(1, 48);
    if body_span > 0.002 {
        for s in 0..segs {
            let t0 = body_a0 + (s as f32 / segs as f32) * body_span;
            let t1 = body_a0 + ((s + 1) as f32 / segs as f32) * body_span;
            let pts = vec![
                Pos2::new(center.x + r_in * t0.cos(), center.y + r_in * t0.sin()),
                Pos2::new(center.x + r_out * t0.cos(), center.y + r_out * t0.sin()),
                Pos2::new(center.x + r_out * t1.cos(), center.y + r_out * t1.sin()),
                Pos2::new(center.x + r_in * t1.cos(), center.y + r_in * t1.sin()),
            ];
            painter.add(Shape::convex_polygon(pts, fill, Stroke::NONE));
        }
    }

    // Draw pointed arrowhead triangle at the strand tip.
    let head_pts = vec![
        Pos2::new(center.x + r_in * tip_base_a.cos(), center.y + r_in * tip_base_a.sin()),
        Pos2::new(center.x + r_out * tip_base_a.cos(), center.y + r_out * tip_base_a.sin()),
        Pos2::new(center.x + r_mid * tip_a.cos(), center.y + r_mid * tip_a.sin()),
    ];
    painter.add(Shape::convex_polygon(head_pts, fill, Stroke::NONE));

    // Crisp outline along outer and inner arc edges.
    if stroke.width > 0.0 {
        let mut outline = Vec::with_capacity(segs * 2 + 5);
        for s in 0..=segs {
            let t = body_a0 + (s as f32 / segs as f32) * body_span;
            outline.push(Pos2::new(center.x + r_out * t.cos(), center.y + r_out * t.sin()));
        }
        if strand == Strand::Forward {
            outline.push(Pos2::new(center.x + r_mid * tip_a.cos(), center.y + r_mid * tip_a.sin()));
        }
        for s in (0..=segs).rev() {
            let t = body_a0 + (s as f32 / segs as f32) * body_span;
            outline.push(Pos2::new(center.x + r_in * t.cos(), center.y + r_in * t.sin()));
        }
        if strand == Strand::Reverse {
            outline.push(Pos2::new(center.x + r_mid * tip_a.cos(), center.y + r_mid * tip_a.sin()));
        }
        painter.add(Shape::closed_line(outline, stroke));
    }
}

fn draw_annotation_tooltip(
    painter: &egui::Painter,
    pointer: Pos2,
    ann: &helix_core::Annotation,
    pal: &Palette,
) {
    let mut lines = format!(
        "{} ({}, strand {})\nRange: {}..{} ({} cols)",
        ann.name,
        ann.kind.label(),
        ann.strand.symbol(),
        ann.start + 1,
        ann.end,
        ann.len()
    );
    if let Some(up) = &ann.uniprot {
        lines.push_str(&format!("\nUniProt: {}  (Double-click for AlphaFold 3D)", up.accession));
    }
    if !ann.notes.is_empty() {
        let short_note = truncate_label(&ann.notes, 72);
        lines.push_str(&format!("\n{short_note}"));
    }

    let galley = painter.layout(lines, FontId::proportional(11.5), pal.text_primary, 290.0);
    let pos = pointer + Vec2::new(14.0, 14.0);
    let box_rect = egui::Rect::from_min_size(pos, galley.size() + Vec2::new(16.0, 12.0));
    painter.rect_filled(box_rect, 6.0, pal.panel_bg);
    painter.rect_stroke(box_rect, 6.0, Stroke::new(1.2, pal.accent));
    painter.galley(pos + Vec2::new(8.0, 6.0), galley, pal.text_primary);
}

fn truncate_label(text: &str, max_chars: usize) -> String {
    let mut chars = text.chars();
    let taken: String = chars.by_ref().take(max_chars).collect();
    if chars.next().is_some() {
        format!("{taken}...")
    } else {
        taken
    }
}
