//! HelixPrime — Open-Source Modern Alignment Inspection & 3D Protein Visualization Suite.
//!
//! Inspired by Geneious Prime, combining high-FPS multiple sequence alignment inspection
//! with live UniProtKB annotation lookup and a separate native OS 3D Protein Viewer window
//! powered by the AlphaFold Protein Structure Database (EBI) and RCSB PDB.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod renderer3d;
mod sys_mem;
mod theme;
mod widgets;

use app::HelixApp;

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("HelixPrime — Open-Source Alignment Inspector & 3D Protein Suite")
            .with_inner_size([1440.0, 880.0])
            .with_min_inner_size([960.0, 620.0]),
        ..Default::default()
    };

    eframe::run_native(
        "HelixPrime",
        options,
        Box::new(|cc| Ok(Box::new(HelixApp::new(cc)))),
    )
}
