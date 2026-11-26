//! HelixPrime — Open-Source Modern Alignment Inspection & 3D Protein Visualization Suite.
//!
//! Inspired by Geneious Prime, combining high-FPS multiple sequence alignment inspection
//! with live UniProtKB annotation lookup and a separate native OS 3D Protein Viewer window
//! powered by the AlphaFold Protein Structure Database (EBI) and RCSB PDB.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod persistence;
mod renderer3d;
mod sys_mem;
mod theme;
mod widgets;

use app::HelixApp;
use std::sync::Arc;

/// Returns the embedded 64x64 RGBA HelixPrime window icon for both the main desktop window
/// and the separate 3D Protein Viewer OS viewport.
pub(crate) fn app_window_icon() -> Arc<egui::IconData> {
    const ICON_RGBA_64: &[u8] = include_bytes!("../../../assets/helix-prime-64.rgba");
    Arc::new(egui::IconData {
        rgba: ICON_RGBA_64.to_vec(),
        width: 64,
        height: 64,
    })
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("HelixPrime — Open-Source Alignment Inspector & 3D Protein Suite")
            .with_icon(app_window_icon())
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
