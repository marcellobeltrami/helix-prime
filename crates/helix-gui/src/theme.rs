//! Visual theme, Geneious Prime-inspired modern light/dark palettes, and bioinformatics
//! residue color schemes.

use egui::{Color32, Rounding, Stroke, Style, Visuals};
use helix_core::sequence::is_gap;
use helix_core::Alphabet;

/// Application UI theme mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ThemeMode {
    /// Clean modern light slate/white desktop aesthetic inspired by Geneious Prime.
    #[default]
    Light,
    /// Modern low-glare dark laboratory aesthetic.
    Dark,
}

/// Residue color palette selector for the alignment canvas.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ColorScheme {
    /// Balanced biochemical palette inspired by Geneious Prime.
    #[default]
    Geneious,
    /// Classic ClustalX sequence alignment color rules.
    Clustal,
    /// Kyte-Doolittle hydrophobicity / polarity spectrum.
    Polarity,
    /// Zappo physicochemical amino acid classes.
    Zappo,
    /// Monochrome text (minimal visual noise for disagreement highlighting).
    Plain,
}

impl ColorScheme {
    pub const ALL: [Self; 5] = [
        Self::Geneious,
        Self::Clustal,
        Self::Polarity,
        Self::Zappo,
        Self::Plain,
    ];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Geneious => "Geneious Standard",
            Self::Clustal => "ClustalX",
            Self::Polarity => "Polarity / Hydrophobicity",
            Self::Zappo => "Zappo Biochemical",
            Self::Plain => "Monochrome",
        }
    }
}

/// Residue highlighting rule matching Geneious Prime's "Highlighting" control section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HighlightMode {
    /// Color all residues according to [`ColorScheme`].
    #[default]
    AllResidues,
    /// Vividly highlight residues that disagree with consensus; mute matching residues.
    Disagreements,
    /// Vividly highlight residues that agree with consensus; mute mismatches.
    Agreements,
}

impl HighlightMode {
    pub const ALL: [Self; 3] = [Self::AllResidues, Self::Disagreements, Self::Agreements];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::AllResidues => "Color All Residues",
            Self::Disagreements => "Highlight Disagreements",
            Self::Agreements => "Highlight Agreements",
        }
    }
}

/// Semantic UI palette derived from [`ThemeMode`].
#[derive(Debug, Clone, Copy)]
pub struct Palette {
    pub canvas_bg: Color32,
    pub panel_bg: Color32,
    pub header_bg: Color32,
    pub row_alt_bg: Color32,
    pub border: Color32,
    pub text_primary: Color32,
    pub text_muted: Color32,
    pub accent: Color32,
    pub accent_soft: Color32,
    pub selection_fill: Color32,
    pub selection_stroke: Color32,
    pub hover_col_fill: Color32,
    pub sync_3d_highlight: Color32,
    pub motif_highlight: Color32,
    pub identity_high: Color32,
    pub identity_mid: Color32,
    pub identity_low: Color32,
    pub viewport_3d_bg_top: Color32,
    pub viewport_3d_bg_bottom: Color32,
}

impl Palette {
    #[must_use]
    pub const fn for_mode(mode: ThemeMode) -> Self {
        match mode {
            ThemeMode::Light => Self {
                canvas_bg: Color32::from_rgb(255, 255, 255),
                panel_bg: Color32::from_rgb(245, 247, 250),
                header_bg: Color32::from_rgb(236, 240, 246),
                row_alt_bg: Color32::from_rgb(249, 250, 252),
                border: Color32::from_rgb(214, 220, 230),
                text_primary: Color32::from_rgb(28, 36, 48),
                text_muted: Color32::from_rgb(100, 112, 128),
                accent: Color32::from_rgb(26, 115, 184),
                accent_soft: Color32::from_rgba_premultiplied(26, 115, 184, 28),
                selection_fill: Color32::from_rgba_premultiplied(37, 132, 214, 48),
                selection_stroke: Color32::from_rgb(26, 115, 184),
                hover_col_fill: Color32::from_rgba_premultiplied(37, 132, 214, 22),
                sync_3d_highlight: Color32::from_rgb(0, 178, 169),
                motif_highlight: Color32::from_rgb(245, 158, 11),
                identity_high: Color32::from_rgb(46, 160, 92),
                identity_mid: Color32::from_rgb(214, 174, 42),
                identity_low: Color32::from_rgb(212, 76, 76),
                viewport_3d_bg_top: Color32::from_rgb(17, 24, 39),
                viewport_3d_bg_bottom: Color32::from_rgb(30, 41, 59),
            },
            ThemeMode::Dark => Self {
                canvas_bg: Color32::from_rgb(22, 26, 33),
                panel_bg: Color32::from_rgb(28, 33, 42),
                header_bg: Color32::from_rgb(35, 41, 52),
                row_alt_bg: Color32::from_rgb(26, 31, 39),
                border: Color32::from_rgb(52, 60, 74),
                text_primary: Color32::from_rgb(230, 236, 245),
                text_muted: Color32::from_rgb(142, 154, 172),
                accent: Color32::from_rgb(66, 165, 245),
                accent_soft: Color32::from_rgba_premultiplied(66, 165, 245, 36),
                selection_fill: Color32::from_rgba_premultiplied(66, 165, 245, 55),
                selection_stroke: Color32::from_rgb(88, 180, 255),
                hover_col_fill: Color32::from_rgba_premultiplied(255, 255, 255, 16),
                sync_3d_highlight: Color32::from_rgb(45, 212, 191),
                motif_highlight: Color32::from_rgb(251, 191, 36),
                identity_high: Color32::from_rgb(52, 199, 114),
                identity_mid: Color32::from_rgb(234, 190, 56),
                identity_low: Color32::from_rgb(235, 87, 87),
                viewport_3d_bg_top: Color32::from_rgb(13, 18, 28),
                viewport_3d_bg_bottom: Color32::from_rgb(24, 33, 49),
            },
        }
    }
}

/// Applies the modern Geneious Prime-inspired visual styling to an `egui::Context`.
pub fn apply_theme(ctx: &egui::Context, mode: ThemeMode) {
    let pal = Palette::for_mode(mode);
    let mut style = Style::default();
    let mut visuals = match mode {
        ThemeMode::Light => Visuals::light(),
        ThemeMode::Dark => Visuals::dark(),
    };

    visuals.panel_fill = pal.panel_bg;
    visuals.window_fill = pal.panel_bg;
    visuals.extreme_bg_color = pal.canvas_bg;
    visuals.faint_bg_color = pal.row_alt_bg;
    visuals.selection.bg_fill = pal.accent_soft;
    visuals.selection.stroke = Stroke::new(1.0, pal.accent);
    visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, pal.border);
    visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, pal.text_primary);
    visuals.widgets.inactive.rounding = Rounding::same(5.0);
    visuals.widgets.hovered.rounding = Rounding::same(5.0);
    visuals.widgets.active.rounding = Rounding::same(5.0);
    visuals.widgets.noninteractive.rounding = Rounding::same(5.0);
    visuals.window_rounding = Rounding::same(8.0);

    style.visuals = visuals;
    style.spacing.item_spacing = egui::vec2(7.0, 5.0);
    style.spacing.button_padding = egui::vec2(8.0, 4.0);
    ctx.set_style(style);
}

/// Returns `(background_color, foreground_text_color)` for a residue tile in the alignment grid.
#[must_use]
pub fn residue_colors(
    residue: u8,
    consensus: u8,
    alphabet: Alphabet,
    scheme: ColorScheme,
    highlight: HighlightMode,
    mode: ThemeMode,
) -> (Option<Color32>, Color32) {
    let pal = Palette::for_mode(mode);
    if is_gap(residue) {
        return (None, pal.text_muted);
    }

    let agrees = !is_gap(consensus) && residue == consensus;
    let muted_fg = match mode {
        ThemeMode::Light => Color32::from_rgb(130, 140, 154),
        ThemeMode::Dark => Color32::from_rgb(112, 124, 140),
    };

    match highlight {
        HighlightMode::Disagreements if agrees => return (None, muted_fg),
        HighlightMode::Agreements if !agrees => return (None, muted_fg),
        _ => {}
    }

    if scheme == ColorScheme::Plain {
        if highlight == HighlightMode::Disagreements && !agrees {
            let bg = match mode {
                ThemeMode::Light => Color32::from_rgb(254, 215, 215),
                ThemeMode::Dark => Color32::from_rgb(130, 48, 56),
            };
            return (Some(bg), pal.text_primary);
        }
        return (None, pal.text_primary);
    }

    let base_rgb = if alphabet.is_nucleotide() {
        nucleotide_rgb(residue, scheme)
    } else {
        amino_acid_rgb(residue, scheme)
    };

    let Some([r, g, b]) = base_rgb else {
        return (None, pal.text_primary);
    };

    // Soften tile saturation slightly in light mode for crisp readability, or deepen in dark mode.
    let bg = match mode {
        ThemeMode::Light => Color32::from_rgb(
            blend_u8(r, 255, 0.22),
            blend_u8(g, 255, 0.22),
            blend_u8(b, 255, 0.22),
        ),
        ThemeMode::Dark => Color32::from_rgb(
            blend_u8(r, 20, 0.30),
            blend_u8(g, 24, 0.30),
            blend_u8(b, 32, 0.30),
        ),
    };

    let fg = match mode {
        ThemeMode::Light => Color32::from_rgb(20, 26, 36),
        ThemeMode::Dark => Color32::from_rgb(245, 248, 252),
    };

    (Some(bg), fg)
}

const fn nucleotide_rgb(base: u8, scheme: ColorScheme) -> Option<[u8; 3]> {
    match scheme {
        ColorScheme::Clustal => match base {
            b'A' => Some([235, 118, 92]),
            b'C' => Some([98, 168, 232]),
            b'G' => Some([240, 182, 74]),
            b'T' | b'U' => Some([106, 204, 128]),
            _ => Some([180, 186, 196]),
        },
        _ => match base {
            // Classic Geneious Prime nucleotide palette: A = coral/red, C = cobalt blue,
            // G = amber/gold, T/U = emerald green.
            b'A' => Some([238, 114, 114]),
            b'C' => Some([104, 166, 238]),
            b'G' => Some([244, 192, 78]),
            b'T' | b'U' => Some([112, 206, 134]),
            b'N' => Some([192, 198, 206]),
            _ => Some([188, 152, 224]),
        },
    }
}

const fn amino_acid_rgb(aa: u8, scheme: ColorScheme) -> Option<[u8; 3]> {
    match scheme {
        ColorScheme::Geneious => match aa {
            // Hydrophobic aliphatic: warm peach/gold
            b'A' | b'V' | b'L' | b'I' | b'M' => Some([242, 198, 112]),
            // Aromatic: soft amber-coral
            b'F' | b'W' | b'Y' => Some([235, 154, 114]),
            // Positively charged (basic): crisp sky blue
            b'K' | b'R' | b'H' => Some([106, 168, 238]),
            // Negatively charged (acidic): soft crimson/rose
            b'D' | b'E' => Some([238, 116, 132]),
            // Polar uncharged: fresh mint/green
            b'S' | b'T' | b'N' | b'Q' => Some([114, 204, 146]),
            // Cysteine: vivid yellow
            b'C' => Some([245, 224, 96]),
            // Proline & Glycine: structural lavender/teal
            b'P' => Some([186, 146, 232]),
            b'G' => Some([120, 202, 212]),
            b'*' => Some([150, 156, 168]),
            _ => Some([195, 200, 210]),
        },
        ColorScheme::Clustal => match aa {
            b'A' | b'I' | b'L' | b'M' | b'F' | b'W' | b'V' => Some([128, 179, 255]),
            b'K' | b'R' => Some([240, 102, 102]),
            b'E' | b'D' => Some([204, 118, 224]),
            b'N' | b'Q' | b'S' | b'T' => Some([102, 204, 128]),
            b'C' => Some([240, 148, 182]),
            b'G' => Some([242, 168, 86]),
            b'P' => Some([232, 224, 96]),
            b'H' | b'Y' => Some([88, 196, 196]),
            _ => Some([190, 196, 204]),
        },
        ColorScheme::Polarity => match aa {
            // Polar/charged (cool blues) -> Hydrophobic (warm reds/oranges)
            b'D' | b'E' | b'K' | b'R' => Some([78, 148, 236]),
            b'N' | b'Q' | b'H' => Some([118, 184, 240]),
            b'S' | b'T' | b'Y' | b'G' => Some([168, 214, 196]),
            b'A' | b'P' | b'C' => Some([242, 206, 118]),
            b'V' | b'M' | b'L' | b'I' | b'F' | b'W' => Some([236, 124, 94]),
            _ => Some([190, 196, 204]),
        },
        ColorScheme::Zappo => match aa {
            b'I' | b'L' | b'V' | b'A' | b'M' => Some([255, 175, 175]),
            b'F' | b'W' | b'Y' => Some([255, 200, 100]),
            b'K' | b'R' | b'H' => Some([110, 160, 255]),
            b'D' | b'E' => Some([245, 96, 96]),
            b'S' | b'T' | b'N' | b'Q' => Some([110, 224, 110]),
            b'P' | b'G' => Some([224, 128, 224]),
            b'C' => Some([245, 240, 96]),
            _ => Some([195, 200, 210]),
        },
        ColorScheme::Plain => None,
    }
}

fn blend_u8(a: u8, b: u8, t: f32) -> u8 {
    let af = f32::from(a);
    let bf = f32::from(b);
    (af + (bf - af) * t).clamp(0.0, 255.0) as u8
}

/// Official DeepMind / EBI AlphaFold pLDDT confidence color scale.
#[must_use]
pub fn plddt_color(plddt: f32) -> Color32 {
    if plddt >= 90.0 {
        Color32::from_rgb(0, 83, 214) // Very high (>90)
    } else if plddt >= 70.0 {
        Color32::from_rgb(101, 203, 243) // Confident (70–90)
    } else if plddt >= 50.0 {
        Color32::from_rgb(255, 219, 19) // Low (50–70)
    } else {
        Color32::from_rgb(255, 125, 69) // Very low (<50)
    }
}

/// Maps MSA conservation (`0.0..=100.0`) onto a crisp perceptually distinct scientific ramp.
#[must_use]
pub fn conservation_color(identity_pct: Option<f32>) -> Color32 {
    let Some(pct) = identity_pct else {
        return Color32::from_rgb(122, 132, 148); // Unmapped residues: neutral slate
    };
    let t = (pct / 100.0).clamp(0.0, 1.0);
    if t >= 0.85 {
        // High conservation: vivid crimson/magenta (ConSurf-style)
        let u = (t - 0.85) / 0.15;
        Color32::from_rgb(
            blend_u8(230, 195, u),
            blend_u8(92, 30, u),
            blend_u8(92, 95, u),
        )
    } else if t >= 0.50 {
        let u = (t - 0.50) / 0.35;
        Color32::from_rgb(
            blend_u8(238, 230, u),
            blend_u8(204, 92, u),
            blend_u8(84, 92, u),
        )
    } else {
        let u = t / 0.50;
        Color32::from_rgb(
            blend_u8(72, 238, u),
            blend_u8(168, 204, u),
            blend_u8(224, 84, u),
        )
    }
}

/// Rainbow spectral ramp from N-terminus (`0.0`) to C-terminus (`1.0`).
#[must_use]
pub fn rainbow_color(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    Color32::from_rgb((r * 235.0) as u8, (g * 225.0) as u8, (b * 240.0) as u8)
}

/// Standard CPK chemical element colors for atomic / heteroatom visualization.
#[must_use]
pub fn cpk_element_color(element: &str) -> Color32 {
    match element {
        "C" => Color32::from_rgb(168, 178, 192),
        "N" => Color32::from_rgb(76, 134, 245),
        "O" => Color32::from_rgb(236, 78, 78),
        "S" => Color32::from_rgb(238, 206, 62),
        "P" => Color32::from_rgb(245, 148, 52),
        "ZN" => Color32::from_rgb(148, 142, 212),
        "FE" => Color32::from_rgb(224, 118, 62),
        "MG" => Color32::from_rgb(102, 214, 128),
        _ => Color32::from_rgb(196, 140, 216),
    }
}
