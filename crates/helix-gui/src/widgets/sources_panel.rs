//! Left Geneious Prime-inspired Sources folder tree and Workspace Document Browser.

use crate::app::HelixApp;
use crate::theme::Palette;
use egui::{Color32, Frame, RichText, Rounding, ScrollArea, Stroke, Ui};
use helix_core::Alphabet;

/// Renders the left Sources / Documents sidebar.
pub fn show_sources_panel(app: &mut HelixApp, ui: &mut Ui) {
    let pal = Palette::for_mode(app.theme);

    ui.horizontal(|ui| {
        ui.label(
            RichText::new("SOURCES & DOCUMENTS")
                .strong()
                .small()
                .color(pal.text_muted),
        );
    });
    ui.add_space(4.0);
    ui.separator();

    // Collect distinct folders preserving insertion order.
    let mut folders: Vec<String> = Vec::new();
    for doc in &app.docs {
        if !folders.iter().any(|f| f == &doc.folder) {
            folders.push(doc.folder.clone());
        }
    }

    let mut switch_to: Option<usize> = None;

    ScrollArea::vertical().show(ui, |ui| {
        for folder in &folders {
            egui::CollapsingHeader::new(
                RichText::new(folder).strong().color(pal.text_primary),
            )
            .default_open(true)
            .show(ui, |ui| {
                for (idx, doc) in app.docs.iter().enumerate() {
                    if &doc.folder != folder {
                        continue;
                    }
                    let selected = idx == app.active_doc;
                    let aln = &doc.alignment;

                    let card_fill = if selected {
                        pal.accent_soft
                    } else {
                        Color32::TRANSPARENT
                    };
                    let card_stroke = if selected {
                        Stroke::new(1.0, pal.accent)
                    } else {
                        Stroke::new(1.0, Color32::TRANSPARENT)
                    };

                    let resp = Frame::none()
                        .fill(card_fill)
                        .stroke(card_stroke)
                        .rounding(Rounding::same(5.0))
                        .inner_margin(egui::Margin::symmetric(8.0, 6.0))
                        .show(ui, |ui| {
                            ui.horizontal(|ui| {
                                // Compact Geneious-style alphabet badge.
                                let (badge_txt, badge_col) = match aln.alphabet {
                                    Alphabet::Protein => ("PRO", Color32::from_rgb(132, 94, 212)),
                                    Alphabet::Dna => ("DNA", Color32::from_rgb(38, 148, 112)),
                                    Alphabet::Rna => ("RNA", Color32::from_rgb(214, 128, 46)),
                                };
                                Frame::none()
                                    .fill(badge_col)
                                    .rounding(Rounding::same(3.0))
                                    .inner_margin(egui::Margin::symmetric(4.0, 1.0))
                                    .show(ui, |ui| {
                                        ui.label(
                                            RichText::new(badge_txt)
                                                .strong()
                                                .size(9.5)
                                                .color(Color32::WHITE),
                                        );
                                    });

                                ui.label(
                                    RichText::new(&aln.name)
                                        .strong()
                                        .size(12.5)
                                        .color(pal.text_primary),
                                );
                            });

                            ui.horizontal(|ui| {
                                let unit = if aln.alphabet.is_nucleotide() {
                                    "bp"
                                } else {
                                    "aa"
                                };
                                ui.label(
                                    RichText::new(format!(
                                        "{} seqs  |  {} {unit}  |  {:.0}% id",
                                        aln.row_count(),
                                        aln.width(),
                                        aln.mean_pairwise_identity_pct()
                                    ))
                                    .size(10.5)
                                    .color(pal.text_muted),
                                );
                            });
                        })
                        .response
                        .interact(egui::Sense::click());

                    if resp.clicked() {
                        switch_to = Some(idx);
                    }
                    ui.add_space(2.0);
                }
            });
        }

        ui.add_space(12.0);
        ui.separator();
        ui.add_space(4.0);

        // Quick Database / 3D shortcuts card.
        ui.label(
            RichText::new("EXTERNAL DATABASES")
                .strong()
                .small()
                .color(pal.text_muted),
        );
        ui.add_space(4.0);
        Frame::none()
            .fill(pal.header_bg)
            .stroke(Stroke::new(1.0, pal.border))
            .rounding(Rounding::same(6.0))
            .inner_margin(egui::Margin::same(8.0))
            .show(ui, |ui| {
                ui.label(
                    RichText::new("UniProtKB + AlphaFold DB")
                        .strong()
                        .size(12.0),
                );
                ui.label(
                    RichText::new(
                        "Click any annotation arrow with a '3D' tag in the alignment view to fetch & inspect its AlphaFold structure in a separate window.",
                    )
                    .size(11.0)
                    .color(pal.text_muted),
                );
                ui.add_space(6.0);
                if ui.button("Open Offline 3D Preview").clicked() {
                    app.open_offline_3d_preview();
                }
            });
    });

    if let Some(idx) = switch_to {
        app.select_document(idx);
    }
}
