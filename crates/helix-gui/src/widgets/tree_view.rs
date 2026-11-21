//! Interactive Phylogenetic Tree viewer (Phylogram & Cladogram rendering with branch
//! length labels, evolutionary distance scale bar, and Newick export).

use crate::app::HelixApp;
use crate::theme::Palette;
use egui::{Align2, Color32, FontId, Frame, Pos2, RichText, Sense, Stroke, Ui};
use helix_core::phylo::PhyloNode;

/// Renders the phylogenetic tree stored in the currently active document.
pub fn show_tree_view(app: &mut HelixApp, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);
    let Some(doc) = app.current_doc() else {
        return;
    };
    let Some(ref tree) = doc.tree else {
        ui.centered_and_justified(|ui| {
            ui.vertical_centered(|ui| {
                ui.label(
                    RichText::new("No phylogenetic tree attached to this document.")
                        .color(pal.text_muted),
                );
                ui.add_space(6.0);
                if ui.button("Build Phylogenetic Tree from This Alignment...").clicked() {
                    app.tree_modal.open = true;
                }
            });
        });
        return;
    };

    let newick = tree.to_newick();
    let tree_name = tree.name.clone();
    let method_lbl = tree.method.label();
    let model_lbl = tree.model.label();
    let leaf_count = tree.root.leaf_count().max(1);
    let max_depth = tree.root.max_depth().max(0.001);
    let root_clone = tree.root.clone();

    // Header strip with Newick copy & method metadata.
    Frame::none()
        .fill(pal.panel_bg)
        .stroke(Stroke::new(1.0, pal.border))
        .inner_margin(egui::Margin::symmetric(12.0, 6.0))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(&tree_name).strong().size(13.5));
                ui.label(
                    RichText::new(format!(
                        "Method: {method_lbl}  |  Model: {model_lbl}  |  Taxa: {leaf_count}"
                    ))
                    .small()
                    .color(pal.text_muted),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("Copy Newick").clicked() {
                        ui.ctx().copy_text(newick.clone());
                        app.status_message = String::from("Copied Newick tree string to clipboard.");
                    }
                });
            });
        });

    let avail = ui.available_rect_before_wrap();
    let (rect, _) = ui.allocate_exact_size(avail.size(), Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, pal.canvas_bg);

    let pad_left = 36.0;
    let pad_right = 200.0;
    let pad_top = 34.0;
    let pad_bottom = 52.0;

    let draw_w = (rect.width() - pad_left - pad_right).max(80.0);
    let draw_h = (rect.height() - pad_top - pad_bottom).max(80.0);
    let y_step = draw_h / leaf_count.max(2) as f32;

    let mut leaf_cursor = 0usize;
    draw_subtree(
        &painter,
        &root_clone,
        0.0,
        max_depth,
        rect.min.x + pad_left,
        rect.min.y + pad_top,
        draw_w,
        y_step,
        &mut leaf_cursor,
        &pal,
    );

    // Draw evolutionary distance scale bar along the bottom.
    let scale_val = nice_scale_bar(max_depth * 0.25);
    let bar_px = ((scale_val / max_depth) as f32 * draw_w).clamp(32.0, 180.0);
    let bar_y = rect.max.y - 24.0;
    let bar_x0 = rect.min.x + pad_left;
    let bar_x1 = bar_x0 + bar_px;
    painter.line_segment([Pos2::new(bar_x0, bar_y), Pos2::new(bar_x1, bar_y)], Stroke::new(2.0, pal.text_primary));
    painter.line_segment([Pos2::new(bar_x0, bar_y - 4.0), Pos2::new(bar_x0, bar_y + 4.0)], Stroke::new(1.5, pal.text_primary));
    painter.line_segment([Pos2::new(bar_x1, bar_y - 4.0), Pos2::new(bar_x1, bar_y + 4.0)], Stroke::new(1.5, pal.text_primary));
    painter.text(
        Pos2::new((bar_x0 + bar_x1) * 0.5, bar_y + 11.0),
        Align2::CENTER_CENTER,
        format!("{scale_val:.4} subs/site"),
        FontId::monospace(10.0),
        pal.text_muted,
    );
}

#[allow(clippy::too_many_arguments)]
fn draw_subtree(
    painter: &egui::Painter,
    node: &PhyloNode,
    depth_so_far: f64,
    max_depth: f64,
    origin_x: f32,
    origin_y: f32,
    draw_w: f32,
    y_step: f32,
    leaf_cursor: &mut usize,
    pal: &Palette,
) -> (f32, f32) {
    let x = origin_x + ((depth_so_far / max_depth) as f32).clamp(0.0, 1.0) * draw_w;

    if node.is_leaf() {
        let y = origin_y + (*leaf_cursor as f32 + 0.5) * y_step;
        *leaf_cursor += 1;
        painter.circle_filled(Pos2::new(x, y), 3.5, pal.accent);
        painter.text(
            Pos2::new(x + 10.0, y),
            Align2::LEFT_CENTER,
            &node.label,
            FontId::proportional(12.5),
            pal.text_primary,
        );
        return (x, y);
    }

    let mut child_pts = Vec::with_capacity(node.children.len());
    for child in &node.children {
        let child_depth = depth_so_far + child.branch_length.max(0.0);
        let (cx, cy) = draw_subtree(
            painter,
            child,
            child_depth,
            max_depth,
            origin_x,
            origin_y,
            draw_w,
            y_step,
            leaf_cursor,
            pal,
        );
        child_pts.push((cx, cy, child.branch_length));
    }

    let min_y = child_pts.iter().map(|p| p.1).fold(f32::INFINITY, f32::min);
    let max_y = child_pts
        .iter()
        .map(|p| p.1)
        .fold(f32::NEG_INFINITY, f32::max);
    let mid_y = (min_y + max_y) * 0.5;

    let branch_stroke = Stroke::new(1.8, pal.text_primary);
    // Vertical clade connector bar.
    painter.line_segment([Pos2::new(x, min_y), Pos2::new(x, max_y)], branch_stroke);

    // Horizontal child branches + subtle branch length annotations.
    for (cx, cy, blen) in child_pts {
        painter.line_segment([Pos2::new(x, cy), Pos2::new(cx, cy)], branch_stroke);
        if cx - x > 38.0 {
            painter.text(
                Pos2::new((x + cx) * 0.5, cy - 8.0),
                Align2::CENTER_CENTER,
                format!("{blen:.3}"),
                FontId::monospace(9.0),
                pal.text_muted,
            );
        }
    }

    painter.circle_filled(Pos2::new(x, mid_y), 3.0, Color32::from_rgb(16, 138, 134));
    (x, mid_y)
}

fn nice_scale_bar(raw: f64) -> f64 {
    if raw <= 0.0 {
        return 0.01;
    }
    let exp = raw.log10().floor();
    let base = 10.0f64.powf(exp);
    let frac = raw / base;
    let nice = if frac < 1.5 {
        1.0
    } else if frac < 3.5 {
        2.0
    } else if frac < 7.5 {
        5.0
    } else {
        10.0
    };
    nice * base
}
