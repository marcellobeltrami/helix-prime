//! Molecular biology & analysis modal dialogs: In-Silico Cloning, PCR Primers,
//! Align/Assemble, Phylogenetic Tree Builder, Local BLAST, and Automated Workflows.

use crate::app::{DocTab, HelixApp, LineageEntry, OperationStatus};
use crate::theme::Palette;
use egui::{Color32, Context, Frame, RichText, ScrollArea, Stroke};
use helix_core::aligners::{
    detect_minimap2_binary, run_alignment_algorithm, AlignerAlgorithm, AlignerConfig,
    Minimap2Preset,
};
use helix_core::assembly::{assemble_contigs_denovo, map_reads_to_reference};
use helix_core::blast::{run_blast, BlastProgram};
use helix_core::cloning::{simulate_cloning, CloningMethod, CloningRequest};
use helix_core::features_db::auto_annotate_from_reference_db;
use helix_core::phylo::{build_phylogenetic_tree, DistanceModel, TreeMethod};
use helix_core::primers::{analyze_primer, design_primer_pairs, test_primer_binding};
use helix_core::restriction::COMMON_ENZYMES;
use helix_core::{Alphabet, Alignment, Annotation, FeatureKind, Sequence, Strand};

/// Renders any currently open molecular biology tool dialogs.
pub fn show_all_tool_modals(app: &mut HelixApp, ctx: &Context) {
    show_cloning_modal(app, ctx);
    show_primers_modal(app, ctx);
    show_assemble_modal(app, ctx);
    show_tree_modal(app, ctx);
    show_blast_modal(app, ctx);
    show_workflow_modal(app, ctx);
}

fn show_cloning_modal(app: &mut HelixApp, ctx: &Context) {
    if !app.cloning_modal.open {
        return;
    }
    let pal = Palette::for_mode(app.theme);
    let dna_docs: Vec<(usize, String)> = app
        .docs
        .iter()
        .enumerate()
        .filter(|(_, d)| !d.deleted && d.alignment.alphabet.is_nucleotide() && d.tree.is_none())
        .map(|(i, d)| (i, d.alignment.name.clone()))
        .collect();

    if dna_docs.is_empty() {
        app.cloning_modal.open = false;
        app.status_message = String::from("Cloning requires at least one DNA document.");
        return;
    }

    if !dna_docs.iter().any(|(i, _)| *i == app.cloning_modal.vector_doc_idx) {
        app.cloning_modal.vector_doc_idx = dna_docs[0].0;
    }
    if !dna_docs.iter().any(|(i, _)| *i == app.cloning_modal.insert_doc_idx) {
        app.cloning_modal.insert_doc_idx = dna_docs.last().map_or(dna_docs[0].0, |(i, _)| *i);
    }

    let mut open = app.cloning_modal.open;
    let mut run_cloning = false;
    let mut close_btn = false;

    egui::Window::new("In-Silico Molecular Cloning")
        .open(&mut open)
        .collapsible(false)
        .resizable(true)
        .default_width(470.0)
        .show(ctx, |ui| {
            egui::Grid::new("cloning_grid")
                .num_columns(2)
                .spacing([12.0, 8.0])
                .show(ui, |ui| {
                    ui.label("Cloning Method:");
                    egui::ComboBox::from_id_salt("cloning_method_cb")
                        .selected_text(app.cloning_modal.method.label())
                        .show_ui(ui, |ui| {
                            for m in CloningMethod::ALL {
                                ui.selectable_value(&mut app.cloning_modal.method, m, m.label());
                            }
                        });
                    ui.end_row();

                    ui.label("Backbone Vector:");
                    let vec_lbl = dna_docs
                        .iter()
                        .find(|(i, _)| *i == app.cloning_modal.vector_doc_idx)
                        .map_or("Select...", |(_, n)| n.as_str());
                    egui::ComboBox::from_id_salt("cloning_vec_cb")
                        .selected_text(vec_lbl)
                        .show_ui(ui, |ui| {
                            for (idx, name) in &dna_docs {
                                ui.selectable_value(&mut app.cloning_modal.vector_doc_idx, *idx, name);
                            }
                        });
                    ui.end_row();

                    ui.label("Insert Fragment:");
                    let ins_lbl = dna_docs
                        .iter()
                        .find(|(i, _)| *i == app.cloning_modal.insert_doc_idx)
                        .map_or("Select...", |(_, n)| n.as_str());
                    egui::ComboBox::from_id_salt("cloning_ins_cb")
                        .selected_text(ins_lbl)
                        .show_ui(ui, |ui| {
                            for (idx, name) in &dna_docs {
                                ui.selectable_value(&mut app.cloning_modal.insert_doc_idx, *idx, name);
                            }
                        });
                    ui.end_row();

                    if app.cloning_modal.method == CloningMethod::RestrictionDigest {
                        ui.label("5' Restriction Enzyme:");
                        egui::ComboBox::from_id_salt("cloning_enz1_cb")
                            .selected_text(&app.cloning_modal.enzyme_5prime)
                            .show_ui(ui, |ui| {
                                for e in COMMON_ENZYMES {
                                    ui.selectable_value(
                                        &mut app.cloning_modal.enzyme_5prime,
                                        e.name.to_owned(),
                                        format!("{} ({})", e.name, e.recognition_str()),
                                    );
                                }
                            });
                        ui.end_row();

                        ui.label("3' Enzyme (optional):");
                        egui::ComboBox::from_id_salt("cloning_enz2_cb")
                            .selected_text(if app.cloning_modal.enzyme_3prime.is_empty() {
                                "(Single Digest)"
                            } else {
                                &app.cloning_modal.enzyme_3prime
                            })
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut app.cloning_modal.enzyme_3prime,
                                    String::new(),
                                    "(Single Digest)",
                                );
                                for e in COMMON_ENZYMES {
                                    ui.selectable_value(
                                        &mut app.cloning_modal.enzyme_3prime,
                                        e.name.to_owned(),
                                        format!("{} ({})", e.name, e.recognition_str()),
                                    );
                                }
                            });
                        ui.end_row();
                    }

                    ui.label("Orientation:");
                    ui.checkbox(&mut app.cloning_modal.reverse_insert, "Insert Reverse Complement");
                    ui.end_row();

                    ui.label("Output Construct Name:");
                    ui.text_edit_singleline(&mut app.cloning_modal.construct_name);
                    ui.end_row();
                });

            if let Some(err) = &app.cloning_modal.error {
                ui.add_space(6.0);
                ui.label(RichText::new(err).small().color(pal.identity_low));
            }

            ui.add_space(10.0);
            ui.horizontal(|ui| {
                if ui.button("Simulate Cloning & Create Construct").clicked() {
                    run_cloning = true;
                }
                if ui.button("Close").clicked() {
                    close_btn = true;
                }
            });
        });

    if close_btn {
        open = false;
    }
    if run_cloning {
        let v_idx = app.cloning_modal.vector_doc_idx;
        let i_idx = app.cloning_modal.insert_doc_idx;
        let req = CloningRequest {
            method: app.cloning_modal.method,
            construct_name: app.cloning_modal.construct_name.clone(),
            enzyme_5prime: app.cloning_modal.enzyme_5prime.clone(),
            enzyme_3prime: app.cloning_modal.enzyme_3prime.clone(),
            reverse_insert: app.cloning_modal.reverse_insert,
        };
        let op_id = app.start_operation("In-Silico Cloning", req.method.label());
        let res = simulate_cloning(&app.docs[v_idx].alignment, &app.docs[i_idx].alignment, &req);
        match res {
            Ok(outcome) => {
                let v_id = app.docs[v_idx].id;
                let v_name = app.docs[v_idx].alignment.name.clone();
                let i_id = app.docs[i_idx].id;
                let i_name = app.docs[i_idx].alignment.name.clone();
                let summary = outcome.summary.clone();
                let construct_name = outcome.construct.name.clone();
                app.add_workspace_doc(
                    "Cloning",
                    true,
                    outcome.construct,
                    None,
                    Some(LineageEntry {
                        operation: req.method.label().to_owned(),
                        parent_ids: vec![v_id, i_id],
                        parent_names: vec![v_name, i_name],
                        details: summary.clone(),
                    }),
                );
                app.finish_operation(op_id, OperationStatus::Completed, format!("Created {construct_name}: {summary}"));
                app.status_message = format!("Cloned '{construct_name}' ({summary}) — see Lineage tab!");
                app.cloning_modal.error = None;
                open = false;
            }
            Err(err) => {
                app.finish_operation(op_id, OperationStatus::Failed, err.to_string());
                app.cloning_modal.error = Some(err.to_string());
            }
        }
    }
    app.cloning_modal.open = open;
}

fn show_primers_modal(app: &mut HelixApp, ctx: &Context) {
    if !app.primers_modal.open {
        return;
    }
    let pal = Palette::for_mode(app.theme);
    let mut open = app.primers_modal.open;
    let mut do_design = false;
    let mut do_test = false;
    let mut apply_pair_idx: Option<usize> = None;

    egui::Window::new("PCR Primers — Design & Test")
        .open(&mut open)
        .collapsible(false)
        .resizable(true)
        .default_width(580.0)
        .show(ctx, |ui| {
            ui.label(RichText::new("DESIGN PCR PRIMER PAIRS").strong().small().color(pal.accent));
            ui.horizontal_wrapped(|ui| {
                ui.label("Length:");
                ui.add(egui::DragValue::new(&mut app.primers_modal.config.min_len).range(14..=30));
                ui.label("to");
                ui.add(egui::DragValue::new(&mut app.primers_modal.config.max_len).range(16..=36));
                ui.label("Target Tm (°C):");
                ui.add(egui::DragValue::new(&mut app.primers_modal.config.target_tm).range(45.0..=75.0));
                if ui.button("Design Primers on Active Document").clicked() {
                    do_design = true;
                }
            });

            if !app.primers_modal.designed_pairs.is_empty() {
                ui.add_space(6.0);
                ScrollArea::vertical().max_height(190.0).show(ui, |ui| {
                    for (idx, pair) in app.primers_modal.designed_pairs.iter().enumerate() {
                        Frame::none()
                            .fill(pal.canvas_bg)
                            .stroke(Stroke::new(1.0, pal.border))
                            .rounding(egui::Rounding::same(5.0))
                            .inner_margin(egui::Margin::same(7.0))
                            .show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.label(
                                        RichText::new(format!(
                                            "Pair #{}  |  Product: {} bp  |  ΔTm: {:.1} °C",
                                            idx + 1,
                                            pair.amplicon_len,
                                            pair.tm_diff
                                        ))
                                        .strong(),
                                    );
                                    if ui.small_button("Annotate & Save to Primers Folder").clicked() {
                                        apply_pair_idx = Some(idx);
                                    }
                                });
                                let fwd_str = String::from_utf8_lossy(&pair.forward.sequence);
                                let rev_str = String::from_utf8_lossy(&pair.reverse.sequence);
                                ui.label(
                                    RichText::new(format!(
                                        "Fwd ({}..{}): 5'-{fwd_str}-3'  (Tm {:.1}°C, GC {:.0}%, Score {:.0})",
                                        pair.fwd_start + 1,
                                        pair.fwd_end,
                                        pair.forward.tm_celsius,
                                        pair.forward.gc_pct,
                                        pair.forward.quality_score
                                    ))
                                    .monospace()
                                    .size(10.5),
                                );
                                ui.label(
                                    RichText::new(format!(
                                        "Rev ({}..{}): 5'-{rev_str}-3'  (Tm {:.1}°C, GC {:.0}%, Score {:.0})",
                                        pair.rev_start + 1,
                                        pair.rev_end,
                                        pair.reverse.tm_celsius,
                                        pair.reverse.gc_pct,
                                        pair.reverse.quality_score
                                    ))
                                    .monospace()
                                    .size(10.5),
                                );
                            });
                        ui.add_space(3.0);
                    }
                });
            }

            ui.add_space(8.0);
            ui.separator();
            ui.label(RichText::new("TEST PRIMER BINDING & THERMODYNAMICS").strong().small().color(pal.accent));
            ui.horizontal(|ui| {
                ui.label("5'->3' Oligo:");
                ui.add(egui::TextEdit::singleline(&mut app.primers_modal.test_primer_seq).desired_width(250.0));
                ui.label("Max MM:");
                ui.add(egui::DragValue::new(&mut app.primers_modal.test_max_mismatches).range(0..=4));
                if ui.button("Test Binding").clicked() {
                    do_test = true;
                }
            });
            let metrics = analyze_primer(app.primers_modal.test_primer_seq.as_bytes());
            ui.label(
                RichText::new(format!(
                    "Len: {} bp | Tm: {:.1}°C | GC: {:.1}% | 3' Clamp: {} | Self-dimer run: {} bp | Score: {:.0}/100",
                    metrics.length,
                    metrics.tm_celsius,
                    metrics.gc_pct,
                    if metrics.has_gc_clamp { "Yes" } else { "No" },
                    metrics.max_self_comp_run,
                    metrics.quality_score
                ))
                .small()
                .color(pal.text_muted),
            );
            if !app.primers_modal.binding_hits.is_empty() {
                ui.label(
                    RichText::new(format!(
                        "Found {} binding site(s) on active template:",
                        app.primers_modal.binding_hits.len()
                    ))
                    .small()
                    .strong(),
                );
                for h in app.primers_modal.binding_hits.iter().take(8) {
                    ui.label(
                        RichText::new(format!(
                            "  Pos {}..{} (strand {}, {} mismatch(es))",
                            h.start_col + 1,
                            h.end_col,
                            h.strand.symbol(),
                            h.mismatches
                        ))
                        .monospace()
                        .small(),
                    );
                }
            }
        });

    if do_design {
        if let Some(aln) = app.current_alignment() {
            if aln.alphabet.is_nucleotide() {
                if let Some(seq) = aln.sequences().first() {
                    let (s, e) = app.selected_cols.unwrap_or((0, aln.width()));
                    let (s, e) = if e.saturating_sub(s) < 60 {
                        (0, aln.width())
                    } else {
                        (s, e)
                    };
                    app.primers_modal.designed_pairs =
                        design_primer_pairs(&seq.residues, s, e, app.primers_modal.config, 6);
                }
            }
        }
    }

    if do_test {
        if let Some(aln) = app.current_alignment() {
            if let Some(seq) = aln.sequences().first() {
                app.primers_modal.binding_hits = test_primer_binding(
                    &seq.residues,
                    app.primers_modal.test_primer_seq.as_bytes(),
                    app.primers_modal.test_max_mismatches,
                );
            }
        }
    }

    if let Some(pidx) = apply_pair_idx {
        if let Some(pair) = app.primers_modal.designed_pairs.get(pidx).cloned() {
            let parent_info = app.current_doc().map(|d| (d.id, d.alignment.name.clone()));
            if let Some(doc) = app.current_doc_mut() {
                doc.alignment.add_annotation(
                    Annotation::new(0, format!("Primer_Fwd_{}", pidx + 1), FeatureKind::Primer, pair.fwd_start, pair.fwd_end)
                        .with_strand(Strand::Forward)
                        .with_notes(format!("Tm {:.1}°C, GC {:.0}%", pair.forward.tm_celsius, pair.forward.gc_pct)),
                );
                doc.alignment.add_annotation(
                    Annotation::new(0, format!("Primer_Rev_{}", pidx + 1), FeatureKind::Primer, pair.rev_start, pair.rev_end)
                        .with_strand(Strand::Reverse)
                        .with_notes(format!("Tm {:.1}°C, GC {:.0}%", pair.reverse.tm_celsius, pair.reverse.gc_pct)),
                );
            }
            if let Some((pid, pname)) = parent_info {
                if let Ok(primer_aln) = Alignment::new(
                    format!("{pname}_Primers_{}", pidx + 1),
                    Alphabet::Dna,
                    vec![
                        Sequence::new(format!("{pname}_Fwd"), &pair.forward.sequence)
                            .with_description(format!("Tm {:.1}°C", pair.forward.tm_celsius)),
                        Sequence::new(format!("{pname}_Rev"), &pair.reverse.sequence)
                            .with_description(format!("Tm {:.1}°C", pair.reverse.tm_celsius)),
                    ],
                ) {
                    app.add_workspace_doc(
                        "Primers",
                        false,
                        primer_aln,
                        None,
                        Some(LineageEntry {
                            operation: String::from("Design PCR Primers"),
                            parent_ids: vec![pid],
                            parent_names: vec![pname],
                            details: format!("Expected amplicon {} bp", pair.amplicon_len),
                        }),
                    );
                }
            }
            app.status_message = format!("Annotated PCR primer pair #{} and saved to Primers folder.", pidx + 1);
            open = false;
        }
    }

    app.primers_modal.open = open;
}

fn show_assemble_modal(app: &mut HelixApp, ctx: &Context) {
    if !app.assemble_modal.open {
        return;
    }
    let pal = Palette::for_mode(app.theme);
    let all_docs: Vec<(usize, String)> = app
        .docs
        .iter()
        .enumerate()
        .filter(|(_, d)| !d.deleted && d.tree.is_none())
        .map(|(i, d)| (i, d.alignment.name.clone()))
        .collect();
    if all_docs.is_empty() {
        app.assemble_modal.open = false;
        return;
    }
    if !all_docs.iter().any(|(i, _)| *i == app.assemble_modal.reference_doc_idx) {
        app.assemble_modal.reference_doc_idx = all_docs[0].0;
    }
    if !all_docs.iter().any(|(i, _)| *i == app.assemble_modal.reads_doc_idx) {
        app.assemble_modal.reads_doc_idx = all_docs[0].0;
    }

    let mut open = app.assemble_modal.open;
    let mut run = false;
    let mut browse_minimap2 = false;

    egui::Window::new("Align / Assemble Sequences")
        .open(&mut open)
        .collapsible(false)
        .resizable(true)
        .default_width(540.0)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.radio_value(&mut app.assemble_modal.mode, 0, "Multiple Alignment (MSA / minimap2)");
                ui.radio_value(&mut app.assemble_modal.mode, 1, "De Novo Contig Assembly");
                ui.radio_value(&mut app.assemble_modal.mode, 2, "Map to Reference");
            });
            ui.separator();

            if app.assemble_modal.mode == 0 {
                egui::Grid::new("msa_algo_grid")
                    .num_columns(2)
                    .spacing([12.0, 8.0])
                    .show(ui, |ui| {
                        ui.label("Alignment Algorithm:");
                        egui::ComboBox::from_id_salt("msa_algo_cb")
                            .selected_text(app.assemble_modal.algorithm.label())
                            .width(310.0)
                            .show_ui(ui, |ui| {
                                for alg in AlignerAlgorithm::ALL {
                                    ui.selectable_value(
                                        &mut app.assemble_modal.algorithm,
                                        alg,
                                        alg.label(),
                                    );
                                }
                            });
                        ui.end_row();

                        ui.label("Input Document:");
                        let cur_lbl = all_docs
                            .iter()
                            .find(|(i, _)| *i == app.active_doc)
                            .map_or("(Active Document)", |(_, n)| n.as_str());
                        ui.label(RichText::new(cur_lbl).strong());
                        ui.end_row();

                        if app.assemble_modal.algorithm == AlignerAlgorithm::Minimap2 {
                            ui.label("minimap2 Preset (-x):");
                            egui::ComboBox::from_id_salt("mm2_preset_cb")
                                .selected_text(app.assemble_modal.minimap2_preset.label())
                                .width(270.0)
                                .show_ui(ui, |ui| {
                                    for p in Minimap2Preset::ALL {
                                        ui.selectable_value(
                                            &mut app.assemble_modal.minimap2_preset,
                                            p,
                                            p.label(),
                                        );
                                    }
                                });
                            ui.end_row();

                            ui.label("minimap2 Binary:");
                            ui.horizontal(|ui| {
                                ui.add(
                                    egui::TextEdit::singleline(
                                        &mut app.assemble_modal.minimap2_custom_path,
                                    )
                                    .desired_width(210.0)
                                    .hint_text("Auto (PATH / tools/minimap2.exe)"),
                                );
                                if ui.small_button("Browse...").clicked() {
                                    browse_minimap2 = true;
                                }
                            });
                            ui.end_row();

                            let custom = if app.assemble_modal.minimap2_custom_path.trim().is_empty() {
                                None
                            } else {
                                Some(std::path::Path::new(
                                    app.assemble_modal.minimap2_custom_path.trim(),
                                ))
                            };
                            let detected = detect_minimap2_binary(custom);
                            ui.label("Engine Status:");
                            match detected {
                                Some(p) => {
                                    ui.label(
                                        RichText::new(format!("External binary found: {}", p.display()))
                                            .small()
                                            .color(pal.identity_high),
                                    );
                                }
                                None => {
                                    ui.label(
                                        RichText::new("External minimap2 not on PATH -> uses built-in Rust (w,k) minimizer-chaining engine")
                                            .small()
                                            .color(pal.accent),
                                    );
                                }
                            }
                            ui.end_row();
                        } else {
                            ui.label("Gap Open / Extend:");
                            ui.horizontal(|ui| {
                                ui.add(
                                    egui::DragValue::new(&mut app.assemble_modal.gap_open)
                                        .range(2..=20)
                                        .prefix("Open -"),
                                );
                                ui.add(
                                    egui::DragValue::new(&mut app.assemble_modal.gap_extend)
                                        .range(1..=10)
                                        .prefix("Ext -"),
                                );
                            });
                            ui.end_row();
                        }
                    });

                let desc = match app.assemble_modal.algorithm {
                    AlignerAlgorithm::DeBruijn => "Builds a shared k-mer De Bruijn graph to chain collinear unitig anchors and fills variable spans with DP.",
                    AlignerAlgorithm::ClustalW => "Computes k-mer distances, builds a guide tree, and performs progressive profile-profile alignment.",
                    AlignerAlgorithm::ClustalOmega => "Uses mBed seed embedding to build the guide tree + conservation-weighted HMM-style profile alignment.",
                    AlignerAlgorithm::Muscle => "3-stage pipeline: draft progressive MSA -> improved Kimura guide tree -> iterative Sum-of-Pairs refinement.",
                    AlignerAlgorithm::Mafft => "Physicochemical signal cross-correlation lag guide tree + conservation-weighted iterative refinement.",
                    AlignerAlgorithm::Minimap2 => "Runs external minimap2 (-a -x preset, parsing SAM CIGAR) if installed, or built-in Rust minimizer-chaining.",
                };
                ui.add_space(4.0);
                ui.label(RichText::new(desc).small().color(pal.text_muted));
            } else if app.assemble_modal.mode == 1 {
                ui.horizontal(|ui| {
                    ui.label("Reads Document:");
                    let lbl = all_docs
                        .iter()
                        .find(|(i, _)| *i == app.assemble_modal.reads_doc_idx)
                        .map_or("", |(_, n)| n.as_str());
                    egui::ComboBox::from_id_salt("denovo_reads_cb")
                        .selected_text(lbl)
                        .show_ui(ui, |ui| {
                            for (i, n) in &all_docs {
                                ui.selectable_value(&mut app.assemble_modal.reads_doc_idx, *i, n);
                            }
                        });
                });
                ui.horizontal(|ui| {
                    ui.label("Min Overlap (bp):");
                    ui.add(egui::DragValue::new(&mut app.assemble_modal.min_overlap).range(6..=50));
                });
                ui.horizontal(|ui| {
                    ui.label("Contig Name:");
                    ui.text_edit_singleline(&mut app.assemble_modal.output_name);
                });
            } else {
                ui.horizontal(|ui| {
                    ui.label("Reference Genome/Seq:");
                    let rlbl = all_docs
                        .iter()
                        .find(|(i, _)| *i == app.assemble_modal.reference_doc_idx)
                        .map_or("", |(_, n)| n.as_str());
                    egui::ComboBox::from_id_salt("map_ref_cb")
                        .selected_text(rlbl)
                        .show_ui(ui, |ui| {
                            for (i, n) in &all_docs {
                                ui.selectable_value(&mut app.assemble_modal.reference_doc_idx, *i, n);
                            }
                        });
                });
                ui.horizontal(|ui| {
                    ui.label("Query Reads:");
                    let qlbl = all_docs
                        .iter()
                        .find(|(i, _)| *i == app.assemble_modal.reads_doc_idx)
                        .map_or("", |(_, n)| n.as_str());
                    egui::ComboBox::from_id_salt("map_reads_cb")
                        .selected_text(qlbl)
                        .show_ui(ui, |ui| {
                            for (i, n) in &all_docs {
                                ui.selectable_value(&mut app.assemble_modal.reads_doc_idx, *i, n);
                            }
                        });
                });
            }

            if let Some(note) = &app.assemble_modal.last_note {
                ui.add_space(4.0);
                ui.label(RichText::new(note).small().color(pal.identity_high));
            }
            if let Some(err) = &app.assemble_modal.error {
                ui.label(RichText::new(err).small().color(pal.identity_low));
            }
            ui.add_space(8.0);
            if ui.button("Run Alignment / Assembly").clicked() {
                run = true;
            }
        });

    if browse_minimap2 {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("minimap2 executable", &["exe", ""])
            .pick_file()
        {
            app.assemble_modal.minimap2_custom_path = path.display().to_string();
        }
    }

    if run {
        match app.assemble_modal.mode {
            0 => {
                let alg = app.assemble_modal.algorithm;
                let op_id = app.start_operation(
                    format!("Align ({})", alg.short_name()),
                    alg.label(),
                );
                let custom_path = if app.assemble_modal.minimap2_custom_path.trim().is_empty() {
                    None
                } else {
                    Some(std::path::PathBuf::from(
                        app.assemble_modal.minimap2_custom_path.trim(),
                    ))
                };
                let cfg = AlignerConfig {
                    algorithm: alg,
                    minimap2_path: custom_path,
                    minimap2_preset: app.assemble_modal.minimap2_preset,
                    gap_open: app.assemble_modal.gap_open,
                    gap_extend: app.assemble_modal.gap_extend,
                };
                if let Some(doc) = app.current_doc() {
                    let pid = doc.id;
                    let pname = doc.alignment.name.clone();
                    let alpha = doc.alignment.alphabet;
                    let seqs = doc.alignment.sequences().to_vec();
                    let out_name = format!("{pname} ({})", alg.short_name());
                    match run_alignment_algorithm(out_name, alpha, &seqs, &cfg) {
                        Ok(out) => {
                            let note = out.engine_note.clone();
                            app.add_workspace_doc(
                                "Alignments",
                                false,
                                out.alignment,
                                None,
                                Some(LineageEntry {
                                    operation: format!("Multiple Alignment ({})", alg.short_name()),
                                    parent_ids: vec![pid],
                                    parent_names: vec![pname],
                                    details: note.clone(),
                                }),
                            );
                            app.finish_operation(op_id, OperationStatus::Completed, note.clone());
                            app.assemble_modal.last_note = Some(note.clone());
                            app.assemble_modal.error = None;
                            app.status_message = note;
                            open = false;
                        }
                        Err(err) => {
                            app.finish_operation(op_id, OperationStatus::Failed, err.to_string());
                            app.assemble_modal.error = Some(err.to_string());
                        }
                    }
                }
            }
            1 => {
                let r_idx = app.assemble_modal.reads_doc_idx;
                let op_id = app.start_operation("De Novo Contig Assembly", "Overlap-Layout-Consensus");
                let p_id = app.docs[r_idx].id;
                let p_name = app.docs[r_idx].alignment.name.clone();
                let alpha = app.docs[r_idx].alignment.alphabet;
                let reads = app.docs[r_idx].alignment.sequences().to_vec();
                let out_name = app.assemble_modal.output_name.clone();
                match assemble_contigs_denovo(out_name, &reads, alpha, app.assemble_modal.min_overlap) {
                    Ok(contig) => {
                        let w = contig.width();
                        app.add_workspace_doc(
                            "Contig Assembly",
                            false,
                            contig,
                            None,
                            Some(LineageEntry {
                                operation: String::from("De Novo Contig Assembly (OLC)"),
                                parent_ids: vec![p_id],
                                parent_names: vec![p_name],
                                details: format!("Assembled {} bp consensus contig", w),
                            }),
                        );
                        app.finish_operation(op_id, OperationStatus::Completed, format!("Assembled {w} bp contig"));
                        app.assemble_modal.error = None;
                        open = false;
                    }
                    Err(err) => {
                        app.finish_operation(op_id, OperationStatus::Failed, err.to_string());
                        app.assemble_modal.error = Some(err.to_string());
                    }
                }
            }
            _ => {
                let ref_idx = app.assemble_modal.reference_doc_idx;
                let r_idx = app.assemble_modal.reads_doc_idx;
                let op_id = app.start_operation("Map Reads to Reference", "k-mer seed & placement");
                let ref_id = app.docs[ref_idx].id;
                let ref_name = app.docs[ref_idx].alignment.name.clone();
                let r_id = app.docs[r_idx].id;
                let r_name = app.docs[r_idx].alignment.name.clone();
                let alpha = app.docs[ref_idx].alignment.alphabet;
                let ref_seq = app.docs[ref_idx].alignment.sequences()[0].clone();
                let reads = app.docs[r_idx].alignment.sequences().to_vec();
                match map_reads_to_reference(format!("{r_name} mapped to {ref_name}"), &ref_seq, &reads, alpha) {
                    Ok((mapped_aln, count)) => {
                        app.add_workspace_doc(
                            "Contig Assembly",
                            false,
                            mapped_aln,
                            None,
                            Some(LineageEntry {
                                operation: String::from("Map to Reference"),
                                parent_ids: vec![ref_id, r_id],
                                parent_names: vec![ref_name, r_name],
                                details: format!("Mapped {count} reads to reference"),
                            }),
                        );
                        app.finish_operation(op_id, OperationStatus::Completed, format!("Mapped {count} reads"));
                        app.assemble_modal.error = None;
                        open = false;
                    }
                    Err(err) => {
                        app.finish_operation(op_id, OperationStatus::Failed, err.to_string());
                        app.assemble_modal.error = Some(err.to_string());
                    }
                }
            }
        }
    }
    app.assemble_modal.open = open;
}

fn show_tree_modal(app: &mut HelixApp, ctx: &Context) {
    if !app.tree_modal.open {
        return;
    }
    let pal = Palette::for_mode(app.theme);
    let mut open = app.tree_modal.open;
    let mut build = false;

    egui::Window::new("Build Phylogenetic Tree")
        .open(&mut open)
        .collapsible(false)
        .resizable(false)
        .default_width(360.0)
        .show(ctx, |ui| {
            egui::Grid::new("tree_modal_grid")
                .num_columns(2)
                .spacing([12.0, 8.0])
                .show(ui, |ui| {
                    ui.label("Tree Build Method:");
                    egui::ComboBox::from_id_salt("tree_method_cb")
                        .selected_text(app.tree_modal.method.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut app.tree_modal.method, TreeMethod::NeighborJoining, TreeMethod::NeighborJoining.label());
                            ui.selectable_value(&mut app.tree_modal.method, TreeMethod::Upgma, TreeMethod::Upgma.label());
                        });
                    ui.end_row();

                    ui.label("Genetic Distance Model:");
                    egui::ComboBox::from_id_salt("tree_model_cb")
                        .selected_text(app.tree_modal.model.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut app.tree_modal.model, DistanceModel::JukesCantor, DistanceModel::JukesCantor.label());
                            ui.selectable_value(&mut app.tree_modal.model, DistanceModel::PDistance, DistanceModel::PDistance.label());
                        });
                    ui.end_row();
                });

            if let Some(err) = &app.tree_modal.error {
                ui.label(RichText::new(err).small().color(pal.identity_low));
            }
            ui.add_space(8.0);
            if ui.button("Build Tree").clicked() {
                build = true;
            }
        });

    if build {
        if let Some(doc) = app.current_doc() {
            let pid = doc.id;
            let pname = doc.alignment.name.clone();
            let aln_clone = doc.alignment.clone();
            let method = app.tree_modal.method;
            let model = app.tree_modal.model;
            let op_id = app.start_operation("Build Phylogenetic Tree", format!("{} / {}", method.label(), model.label()));
            match build_phylogenetic_tree(&aln_clone, method, model) {
                Ok(tree) => {
                    let mut tree_aln = aln_clone;
                    tree_aln.name = format!("{pname} ({})", method.label());
                    app.add_workspace_doc(
                        "Tree Documents",
                        false,
                        tree_aln,
                        Some(tree),
                        Some(LineageEntry {
                            operation: format!("Build Tree ({}, {})", method.label(), model.label()),
                            parent_ids: vec![pid],
                            parent_names: vec![pname],
                            details: String::from("Inferred phylogenetic tree from MSA"),
                        }),
                    );
                    app.active_tab = DocTab::TreeView;
                    app.finish_operation(op_id, OperationStatus::Completed, "Tree built successfully");
                    app.tree_modal.error = None;
                    open = false;
                }
                Err(err) => {
                    app.finish_operation(op_id, OperationStatus::Failed, err.to_string());
                    app.tree_modal.error = Some(err.to_string());
                }
            }
        }
    }
    app.tree_modal.open = open;
}

fn show_blast_modal(app: &mut HelixApp, ctx: &Context) {
    if !app.blast_modal.open {
        return;
    }
    let pal = Palette::for_mode(app.theme);
    let mut open = app.blast_modal.open;
    let mut trigger_blast = false;
    let mut jump_to: Option<(usize, usize, usize)> = None;

    egui::Window::new("BLAST Similarity Search (totalBLAST / blastn / blastp / blastx / tblastn / tblastx / MegaBLAST)")
        .open(&mut open)
        .collapsible(false)
        .resizable(true)
        .default_width(720.0)
        .show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Program:").strong());
                egui::ComboBox::from_id_salt("blast_program_cb")
                    .selected_text(app.blast_modal.config.program.label())
                    .width(290.0)
                    .show_ui(ui, |ui| {
                        for prog in BlastProgram::ALL {
                            ui.selectable_value(
                                &mut app.blast_modal.config.program,
                                prog,
                                prog.label(),
                            );
                        }
                    });

                ui.label("Min ID %:");
                ui.add(
                    egui::DragValue::new(&mut app.blast_modal.config.min_identity_pct)
                        .range(20.0..=100.0)
                        .suffix("%"),
                );

                ui.label("Max E:");
                ui.add(
                    egui::DragValue::new(&mut app.blast_modal.config.max_e_value)
                        .range(1e-20..=1000.0)
                        .speed(0.5),
                );
            });
            ui.add_space(4.0);
            ui.label(RichText::new("Query Sequence (DNA / RNA / Protein):").small().strong());
            ui.add(
                egui::TextEdit::multiline(&mut app.blast_modal.query)
                    .font(egui::TextStyle::Monospace)
                    .desired_rows(3)
                    .desired_width(f32::INFINITY),
            );
            ui.horizontal(|ui| {
                if ui.button(format!("Run {}", app.blast_modal.config.program.short_tag())).clicked() {
                    trigger_blast = true;
                }
                ui.label(
                    RichText::new("totalBLAST automatically searches BOTH nucleotide strands AND 6-frame translated protein databases!")
                        .small()
                        .color(pal.text_muted),
                );
            });
            ui.separator();

            if app.blast_modal.hits.is_empty() {
                ui.label(
                    RichText::new("Click Run above to search workspace documents across DNA & 6-frame translated Protein space.")
                        .small()
                        .color(pal.text_muted),
                );
            } else {
                ScrollArea::vertical().max_height(210.0).show(ui, |ui| {
                    egui::Grid::new("blast_hits_grid")
                        .num_columns(7)
                        .striped(true)
                        .spacing([10.0, 4.0])
                        .show(ui, |ui| {
                            for hdr in ["Mode / Frame", "Document", "Subject", "Identity / Pos", "Bit Score", "E-value", "Action"] {
                                ui.label(RichText::new(hdr).strong().small().color(pal.text_muted));
                            }
                            ui.end_row();

                            for (idx, hit) in app.blast_modal.hits.iter().enumerate() {
                                let sel = app.blast_modal.selected_hit == Some(idx);
                                ui.label(
                                    RichText::new(&hit.program_tag)
                                        .monospace()
                                        .small()
                                        .strong()
                                        .color(pal.accent),
                                );
                                if ui.selectable_label(sel, &hit.doc_name).clicked() {
                                    app.blast_modal.selected_hit = Some(idx);
                                }
                                ui.label(RichText::new(&hit.subject_name).small());
                                ui.label(
                                    RichText::new(format!(
                                        "{:.0}% / {:.0}% ({} len)",
                                        hit.identity_pct, hit.positives_pct, hit.align_len
                                    ))
                                    .monospace()
                                    .small()
                                    .color(pal.identity_high),
                                );
                                ui.label(RichText::new(format!("{:.1}", hit.bit_score)).monospace().small());
                                ui.label(RichText::new(format!("{:.2e}", hit.e_value)).monospace().small());
                                if ui.small_button("Open").clicked() {
                                    jump_to = Some((hit.doc_index, hit.subject_start, hit.subject_end));
                                }
                                ui.end_row();
                            }
                        });
                });

                if let Some(sidx) = app.blast_modal.selected_hit {
                    if let Some(hit) = app.blast_modal.hits.get(sidx) {
                        ui.add_space(6.0);
                        Frame::none()
                            .fill(pal.canvas_bg)
                            .stroke(Stroke::new(1.0, pal.border))
                            .inner_margin(egui::Margin::same(6.0))
                            .show(ui, |ui| {
                                ui.label(
                                    RichText::new(format!(
                                        "HSP Alignment Preview [{}] ({}..{} on {}):",
                                        hit.program_tag,
                                        hit.subject_start + 1,
                                        hit.subject_end,
                                        hit.subject_name
                                    ))
                                    .small()
                                    .strong(),
                                );
                                let preview = format!(
                                    "Query {:>4}  {}  {}\n            {}\nSubj  {:>4}  {}  {}",
                                    hit.query_start + 1,
                                    hit.query_snippet,
                                    hit.query_end,
                                    hit.match_midline,
                                    hit.subject_start + 1,
                                    hit.subject_snippet,
                                    hit.subject_end
                                );
                                ui.label(RichText::new(preview).monospace().size(10.5));
                            });
                    }
                }
            }
        });

    if trigger_blast {
        let prog = app.blast_modal.config.program;
        let op_id = app.start_operation(
            format!("BLAST ({})", prog.short_tag()),
            "Searching all local workspace documents",
        );
        let subjects: Vec<(usize, &Alignment)> = app
            .docs
            .iter()
            .enumerate()
            .filter(|(_, d)| !d.deleted && d.tree.is_none())
            .map(|(i, d)| (i, &d.alignment))
            .collect();
        let hits = run_blast(&app.blast_modal.query, &subjects, app.blast_modal.config);
        let n = hits.len();
        app.blast_modal.hits = hits;
        app.blast_modal.selected_hit = if n > 0 { Some(0) } else { None };
        app.finish_operation(
            op_id,
            OperationStatus::Completed,
            format!("{} found {n} hit(s)", prog.short_tag()),
        );
    }

    if let Some((doc_idx, s, e)) = jump_to {
        app.select_document(doc_idx);
        app.selected_cols = Some((s, e));
        app.cursor_col = s;
        app.scroll_to_column = Some(s);
        app.active_tab = DocTab::SequenceView;
        open = false;
    }
    app.blast_modal.open = open;
}

fn show_workflow_modal(app: &mut HelixApp, ctx: &Context) {
    if !app.workflow_modal.open {
        return;
    }
    let pal = Palette::for_mode(app.theme);
    let mut open = app.workflow_modal.open;
    let mut run_wf = false;

    egui::Window::new("Automated Workflows Pipeline")
        .open(&mut open)
        .collapsible(false)
        .resizable(false)
        .default_width(430.0)
        .show(ctx, |ui| {
            ui.label(
                RichText::new("Chain multiple bioinformatics steps to execute automatically on the active document:")
                    .small()
                    .color(pal.text_muted),
            );
            ui.add_space(6.0);
            ui.checkbox(
                &mut app.workflow_modal.step_auto_annotate,
                "1. Auto-Annotate from Reference Plasmid/Phage Features",
            );
            ui.checkbox(
                &mut app.workflow_modal.step_realign,
                "2. Progressive Needleman-Wunsch Multiple Sequence Alignment",
            );
            ui.checkbox(
                &mut app.workflow_modal.step_translate,
                "3. Translate DNA/RNA -> Protein Document (Frame +1)",
            );
            ui.checkbox(
                &mut app.workflow_modal.step_build_tree,
                "4. Build Neighbor-Joining Phylogenetic Tree (for MSAs >= 2 seqs)",
            );
            ui.checkbox(
                &mut app.workflow_modal.step_design_primers,
                "5. Auto-Design Top PCR Primer Pair & Annotate",
            );

            if let Some(report) = &app.workflow_modal.last_report {
                ui.add_space(8.0);
                Frame::none()
                    .fill(pal.canvas_bg)
                    .stroke(Stroke::new(1.0, pal.border))
                    .inner_margin(egui::Margin::same(6.0))
                    .show(ui, |ui| {
                        ui.label(RichText::new(report).small().color(Color32::from_rgb(38, 152, 92)));
                    });
            }

            ui.add_space(10.0);
            if ui.button("Run Workflow Pipeline").clicked() {
                run_wf = true;
            }
        });

    if run_wf {
        let op_id = app.start_operation("Workflow Pipeline", "Executing chained steps...");
        let mut steps_log = Vec::new();
        let start_doc_idx = app.active_doc;

        if app.workflow_modal.step_auto_annotate {
            if let Some(doc) = app.docs.get_mut(start_doc_idx) {
                let added = auto_annotate_from_reference_db(&mut doc.alignment, 90.0);
                steps_log.push(format!("Auto-Annotate: +{added} features"));
            }
        }
        if app.workflow_modal.step_realign {
            if let Some(doc) = app.docs.get_mut(start_doc_idx) {
                if doc.alignment.row_count() >= 2 && doc.alignment.realign_progressive().is_ok() {
                    doc.invalidate_caches();
                    steps_log.push(String::from("Aligned (NW)"));
                }
            }
        }
        if app.workflow_modal.step_design_primers {
            if let Some(doc) = app.docs.get_mut(start_doc_idx) {
                if doc.alignment.alphabet.is_nucleotide() {
                    if let Some(seq) = doc.alignment.sequences().first() {
                        let pairs = design_primer_pairs(
                            &seq.residues,
                            0,
                            doc.alignment.width(),
                            helix_core::PrimerDesignConfig::default(),
                            1,
                        );
                        if let Some(top) = pairs.first() {
                            doc.alignment.add_annotation(
                                Annotation::new(0, "WF_Primer_Fwd", FeatureKind::Primer, top.fwd_start, top.fwd_end)
                                    .with_strand(Strand::Forward),
                            );
                            doc.alignment.add_annotation(
                                Annotation::new(0, "WF_Primer_Rev", FeatureKind::Primer, top.rev_start, top.rev_end)
                                    .with_strand(Strand::Reverse),
                            );
                            steps_log.push(format!("Primers: {} bp amplicon", top.amplicon_len));
                        }
                    }
                }
            }
        }
        if app.workflow_modal.step_translate {
            if let Some(doc) = app.docs.get(start_doc_idx) {
                if doc.alignment.alphabet.is_nucleotide() {
                    let pid = doc.id;
                    let pname = doc.alignment.name.clone();
                    if let Ok(prot) = doc.alignment.translate_to_protein_alignment(0) {
                        app.add_workspace_doc(
                            "Protein Documents",
                            false,
                            prot,
                            None,
                            Some(LineageEntry {
                                operation: String::from("Workflow: Translate"),
                                parent_ids: vec![pid],
                                parent_names: vec![pname],
                                details: String::from("Automated workflow translation"),
                            }),
                        );
                        steps_log.push(String::from("Translated -> Protein Documents"));
                    }
                }
            }
        }
        if app.workflow_modal.step_build_tree {
            if let Some(doc) = app.docs.get(start_doc_idx) {
                if doc.alignment.row_count() >= 2 {
                    let pid = doc.id;
                    let pname = doc.alignment.name.clone();
                    let aln_clone = doc.alignment.clone();
                    if let Ok(tree) =
                        build_phylogenetic_tree(&aln_clone, TreeMethod::NeighborJoining, DistanceModel::JukesCantor)
                    {
                        let mut tree_aln = aln_clone;
                        tree_aln.name = format!("{pname} (Workflow Tree)");
                        app.add_workspace_doc(
                            "Tree Documents",
                            false,
                            tree_aln,
                            Some(tree),
                            Some(LineageEntry {
                                operation: String::from("Workflow: Build Tree (NJ)"),
                                parent_ids: vec![pid],
                                parent_names: vec![pname],
                                details: String::from("Automated workflow phylogenetic tree"),
                            }),
                        );
                        steps_log.push(String::from("Built NJ Tree"));
                    }
                }
            }
        }

        let summary = if steps_log.is_empty() {
            String::from("No applicable steps ran for this document.")
        } else {
            format!("Completed: {}", steps_log.join(" -> "))
        };
        app.finish_operation(op_id, OperationStatus::Completed, summary.clone());
        app.workflow_modal.last_report = Some(summary.clone());
        app.status_message = summary;
    }
    app.workflow_modal.open = open;
}
