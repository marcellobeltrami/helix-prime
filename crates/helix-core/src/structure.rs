//! 3D protein structure data model, PDB/mmCIF parsers, geometric secondary-structure
//! inference, and bidirectional alignment-to-3D-residue coordinate mapping.

use crate::alignment::Alignment;
use crate::annotation::Annotation;
use crate::error::{HelixError, Result};
use crate::sequence::{is_gap, three_to_one_aa, translate_aligned_region};
use std::ops::{Add, Sub};

/// Compact 3D Euclidean vector used for atomic coordinates, normals, and spline control points.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    #[inline]
    #[must_use]
    pub const fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    #[inline]
    #[must_use]
    pub fn scale(self, s: f32) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }

    #[inline]
    #[must_use]
    pub fn dot(self, rhs: Self) -> f32 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    #[inline]
    #[must_use]
    pub fn cross(self, rhs: Self) -> Self {
        Self::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }

    #[inline]
    #[must_use]
    pub fn length_sq(self) -> f32 {
        self.dot(self)
    }

    #[inline]
    #[must_use]
    pub fn length(self) -> f32 {
        self.length_sq().sqrt()
    }

    #[inline]
    #[must_use]
    pub fn distance(self, rhs: Self) -> f32 {
        (self - rhs).length()
    }

    #[inline]
    #[must_use]
    pub fn normalize(self) -> Self {
        let len = self.length();
        if len > 1e-6 {
            self.scale(1.0 / len)
        } else {
            Self::new(0.0, 1.0, 0.0)
        }
    }

    #[inline]
    #[must_use]
    pub fn lerp(self, rhs: Self, t: f32) -> Self {
        Self::new(
            self.x + (rhs.x - self.x) * t,
            self.y + (rhs.y - self.y) * t,
            self.z + (rhs.z - self.z) * t,
        )
    }
}

impl std::ops::Add for Vec3 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

/// Secondary structure state for a protein residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    #[default]
    Coil,
}

impl SecondaryStructure {
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Helix => "Alpha Helix",
            Self::Sheet => "Beta Strand",
            Self::Coil => "Coil / Loop",
        }
    }
}

/// Compact atomic coordinate record ordered by alignment to avoid struct padding bloat.
#[derive(Debug, Clone, Copy)]
pub struct Atom {
    /// 3D Cartesian coordinate in Angstroms.
    pub pos: Vec3,
    /// B-factor (experimental crystallographic factor, or pLDDT 0..=100 for AlphaFold models).
    pub b_factor: f32,
    /// Index into [`ProteinStructure::residues`], or `u32::MAX` for unattached heteroatoms.
    pub residue_idx: u32,
    /// Up to 4-byte ASCII PDB atom name (trimmed, e.g. `b"CA  "`).
    pub name: [u8; 4],
    /// Up to 2-byte ASCII chemical element symbol (e.g. `b"C "`, `b"N "`, `b"O "`, `b"ZN"`).
    pub element: [u8; 2],
    /// `true` for `HETATM` records (ligands, metal ions, heme groups).
    pub is_hetatm: bool,
    /// `true` if this atom is an alpha-carbon (`CA`) or nucleotide `P`/`C4'`.
    pub is_ca: bool,
}

impl Atom {
    /// Returns the element symbol as a trimmed uppercase ASCII string slice.
    #[must_use]
    pub fn element_str(&self) -> &str {
        std::str::from_utf8(&self.element).unwrap_or("C").trim()
    }
}

/// Polymer residue in a 3D protein structure, referencing a contiguous slice of [`Atom`]s.
#[derive(Debug, Clone, Copy)]
pub struct Residue3D {
    /// Index of the first atom belonging to this residue in `ProteinStructure::atoms`.
    pub atom_start: u32,
    /// Index of the alpha-carbon (`CA`) atom if present.
    pub ca_atom: Option<u32>,
    /// Sequence residue number from the PDB/mmCIF file (`auth_seq_id` / `resSeq`).
    pub seq_id: i32,
    /// AlphaFold pLDDT confidence score (`0.0..=100.0`) or experimental B-factor.
    pub plddt: f32,
    /// Chain identifier (e.g. `'A'`).
    pub chain_id: char,
    /// Number of atoms belonging to this residue.
    pub atom_count: u16,
    /// 3-letter residue name (e.g. `b"ALA"`).
    pub name: [u8; 3],
    /// 1-letter IUPAC amino acid code (`A`..`Z`, or `X`).
    pub one_letter: u8,
    /// Assigned or inferred secondary structure element.
    pub secondary_structure: SecondaryStructure,
}

impl Residue3D {
    /// Returns the 3-letter residue name as a UTF-8 string slice.
    #[must_use]
    pub fn three_letter(&self) -> &str {
        std::str::from_utf8(&self.name).unwrap_or("UNK").trim()
    }
}

/// Immutable, cache-friendly 3D protein structure loaded from AlphaFold DB, PDB, or local file.
///
/// # Invariants
/// - `residues` is non-empty and every `Residue3D::atom_start + atom_count <= atoms.len()`.
/// - `radius > 0.0` (clamped to at least `1.0` Angstrom so camera framing never divides by zero).
#[derive(Debug, Clone)]
pub struct ProteinStructure {
    /// Display title (e.g., `"Cellular tumor antigen p53"`).
    pub title: String,
    /// Associated UniProt accession, if known (e.g., `"P04637"`).
    pub accession: Option<String>,
    /// Database entry identifier (e.g., `"AF-P04637-F1"`).
    pub entry_id: String,
    /// Organism scientific name, if available.
    pub organism: Option<String>,
    /// Data source tag (e.g., `"AlphaFold DB v4"`, `"Local PDB"`, `"Offline Preview"`).
    pub source_label: String,
    /// `true` when B-factors represent AlphaFold per-residue pLDDT (`0..=100`).
    pub is_alphafold: bool,
    /// Contiguous atom buffer.
    pub atoms: Box<[Atom]>,
    /// Contiguous polymer residue buffer.
    pub residues: Box<[Residue3D]>,
    /// Centroid of all C-alpha (or fallback all) atoms.
    pub center: Vec3,
    /// Bounding sphere radius in Angstroms around `center`.
    pub radius: f32,
    /// Mean pLDDT confidence score across all residues.
    pub mean_plddt: f32,
}

impl ProteinStructure {
    /// Returns the 3D position of residue `idx`'s C-alpha atom (or first atom fallback).
    #[must_use]
    pub fn residue_ca_pos(&self, idx: usize) -> Option<Vec3> {
        let res = self.residues.get(idx)?;
        if let Some(ca_idx) = res.ca_atom {
            return self.atoms.get(ca_idx as usize).map(|a| a.pos);
        }
        self.atoms.get(res.atom_start as usize).map(|a| a.pos)
    }

    /// Returns the contiguous slice of atoms belonging to residue `idx`.
    #[must_use]
    pub fn residue_atoms(&self, idx: usize) -> &[Atom] {
        let Some(res) = self.residues.get(idx) else {
            return &[];
        };
        let start = res.atom_start as usize;
        let end = start.saturating_add(usize::from(res.atom_count)).min(self.atoms.len());
        &self.atoms[start..end]
    }

    /// Extracts the 1-letter amino acid sequence of all polymer residues in chain order.
    #[must_use]
    pub fn one_letter_sequence(&self) -> Vec<u8> {
        self.residues.iter().map(|r| r.one_letter).collect()
    }
}

/// Secondary-structure range parsed from PDB `HELIX`/`SHEET` or mmCIF `_struct_conf`.
#[derive(Debug, Clone, Copy)]
struct SsRange {
    chain: char,
    start_seq: i32,
    end_seq: i32,
    kind: SecondaryStructure,
}

/// Parses a PDB-format string (including AlphaFold DB `.pdb` files) into a [`ProteinStructure`].
pub fn parse_pdb(
    content: &str,
    entry_id: impl Into<String>,
    accession: Option<String>,
    source_label: impl Into<String>,
    is_alphafold: bool,
) -> Result<ProteinStructure> {
    let mut title = String::new();
    let mut ss_ranges: Vec<SsRange> = Vec::new();
    let mut atoms: Vec<Atom> = Vec::with_capacity(2048);
    let mut residues: Vec<Residue3D> = Vec::with_capacity(256);

    let mut current_res_key: Option<(char, i32, [u8; 3])> = None;

    for line in content.lines() {
        // Stop after the first model in multi-model NMR ensembles.
        if line.starts_with("ENDMDL") && !atoms.is_empty() {
            break;
        }

        if line.starts_with("TITLE") {
            if let Some(frag) = line.get(10..) {
                let t = frag.trim();
                if !t.is_empty() {
                    if !title.is_empty() {
                        title.push(' ');
                    }
                    title.push_str(t);
                }
            }
            continue;
        }

        if line.starts_with("HELIX ") && line.len() >= 37 {
            let chain = line.as_bytes().get(19).copied().map(char::from).unwrap_or('A');
            if let (Some(s), Some(e)) = (
                parse_trimmed_i32(line.get(21..25)),
                parse_trimmed_i32(line.get(33..37)),
            ) {
                ss_ranges.push(SsRange {
                    chain,
                    start_seq: s,
                    end_seq: e,
                    kind: SecondaryStructure::Helix,
                });
            }
            continue;
        }

        if line.starts_with("SHEET ") && line.len() >= 37 {
            let chain = line.as_bytes().get(21).copied().map(char::from).unwrap_or('A');
            if let (Some(s), Some(e)) = (
                parse_trimmed_i32(line.get(22..26)),
                parse_trimmed_i32(line.get(33..37)),
            ) {
                ss_ranges.push(SsRange {
                    chain,
                    start_seq: s,
                    end_seq: e,
                    kind: SecondaryStructure::Sheet,
                });
            }
            continue;
        }

        let is_atom = line.starts_with("ATOM  ");
        let is_hetatm = line.starts_with("HETATM");
        if (!is_atom && !is_hetatm) || line.len() < 54 {
            continue;
        }

        // Skip alternate conformations other than primary (' ' or 'A').
        let alt_loc = line.as_bytes().get(16).copied().unwrap_or(b' ');
        if alt_loc != b' ' && alt_loc != b'A' {
            continue;
        }

        let Some(x) = parse_trimmed_f32(line.get(30..38)) else {
            continue;
        };
        let Some(y) = parse_trimmed_f32(line.get(38..46)) else {
            continue;
        };
        let Some(z) = parse_trimmed_f32(line.get(46..54)) else {
            continue;
        };
        let b_factor = parse_trimmed_f32(line.get(60..66)).unwrap_or(50.0);

        let atom_name_raw = line.get(12..16).unwrap_or("    ").trim();
        let res_name_raw = line.get(17..20).unwrap_or("UNK").trim();
        let chain_id = line
            .as_bytes()
            .get(21)
            .copied()
            .map(char::from)
            .filter(|c| !c.is_ascii_whitespace())
            .unwrap_or('A');
        let seq_id = parse_trimmed_i32(line.get(22..26)).unwrap_or(0);

        // Skip solvent water molecules to keep 3D viewport crisp and uncluttered.
        if is_hetatm && matches!(res_name_raw, "HOH" | "WAT" | "DOD") {
            continue;
        }

        let elem_raw = line
            .get(76..78)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| infer_element_from_atom_name(atom_name_raw));

        // Skip hydrogens in heavy-atom views to keep atom counts fast and readable.
        if elem_raw.eq_ignore_ascii_case("H") {
            continue;
        }

        let name_bytes = copy_ascii_fixed::<4>(atom_name_raw);
        let elem_bytes = copy_ascii_fixed::<2>(elem_raw);
        let res_bytes = copy_ascii_fixed::<3>(res_name_raw);
        let is_ca = !is_hetatm && atom_name_raw == "CA";

        let atom_idx = atoms.len() as u32;
        let residue_idx = if is_hetatm {
            u32::MAX
        } else {
            let key = (chain_id, seq_id, res_bytes);
            if current_res_key != Some(key) {
                current_res_key = Some(key);
                residues.push(Residue3D {
                    atom_start: atom_idx,
                    ca_atom: if is_ca { Some(atom_idx) } else { None },
                    seq_id,
                    plddt: b_factor,
                    chain_id,
                    atom_count: 1,
                    name: res_bytes,
                    one_letter: three_to_one_aa(res_name_raw),
                    secondary_structure: SecondaryStructure::Coil,
                });
            } else if let Some(last) = residues.last_mut() {
                last.atom_count = last.atom_count.saturating_add(1);
                if is_ca {
                    last.ca_atom = Some(atom_idx);
                    last.plddt = b_factor;
                }
            }
            (residues.len().saturating_sub(1)) as u32
        };

        atoms.push(Atom {
            pos: Vec3::new(x, y, z),
            b_factor,
            residue_idx,
            name: name_bytes,
            element: elem_bytes,
            is_hetatm,
            is_ca,
        });
    }

    finalize_structure(
        entry_id.into(),
        title,
        accession,
        source_label.into(),
        is_alphafold,
        atoms,
        residues,
        &ss_ranges,
    )
}

/// Parses an mmCIF (`.cif`) protein structure file into a [`ProteinStructure`].
pub fn parse_mmcif(
    content: &str,
    entry_id: impl Into<String>,
    accession: Option<String>,
    source_label: impl Into<String>,
    is_alphafold: bool,
) -> Result<ProteinStructure> {
    let mut title = String::new();
    let mut ss_ranges: Vec<SsRange> = Vec::new();
    let mut atoms: Vec<Atom> = Vec::with_capacity(2048);
    let mut residues: Vec<Residue3D> = Vec::with_capacity(256);
    let mut current_res_key: Option<(char, i32, [u8; 3])> = None;

    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0usize;

    while i < lines.len() {
        let trimmed = lines[i].trim();
        if trimmed.starts_with("_struct.title") {
            let rest = trimmed.trim_start_matches("_struct.title").trim();
            if !rest.is_empty() {
                title = rest.trim_matches('\'').trim_matches('"').to_owned();
            } else if let Some(next_line) = lines.get(i + 1) {
                title = next_line.trim().trim_matches(';').trim_matches('\'').to_owned();
            }
            i += 1;
            continue;
        }

        if trimmed == "loop_" {
            i += 1;
            let mut cols: Vec<&str> = Vec::new();
            while i < lines.len() && lines[i].trim().starts_with('_') {
                cols.push(lines[i].trim());
                i += 1;
            }

            let is_atom_site = cols.first().is_some_and(|c| c.starts_with("_atom_site."));
            let is_struct_conf = cols.first().is_some_and(|c| c.starts_with("_struct_conf."));
            let is_sheet_range = cols
                .first()
                .is_some_and(|c| c.starts_with("_struct_sheet_range."));

            if is_atom_site {
                let c_group = find_col(&cols, "_atom_site.group_PDB");
                let c_atom = find_col(&cols, "_atom_site.label_atom_id");
                let c_alt = find_col(&cols, "_atom_site.label_alt_id");
                let c_comp = find_col(&cols, "_atom_site.label_comp_id");
                let c_chain = find_col(&cols, "_atom_site.auth_asym_id")
                    .or_else(|| find_col(&cols, "_atom_site.label_asym_id"));
                let c_seq = find_col(&cols, "_atom_site.auth_seq_id")
                    .or_else(|| find_col(&cols, "_atom_site.label_seq_id"));
                let c_x = find_col(&cols, "_atom_site.Cartn_x");
                let c_y = find_col(&cols, "_atom_site.Cartn_y");
                let c_z = find_col(&cols, "_atom_site.Cartn_z");
                let c_bfac = find_col(&cols, "_atom_site.B_iso_or_equiv");
                let c_elem = find_col(&cols, "_atom_site.type_symbol");
                let c_model = find_col(&cols, "_atom_site.pdbx_PDB_model_num");

                while i < lines.len() {
                    let row_line = lines[i].trim();
                    if row_line.is_empty()
                        || row_line.starts_with('#')
                        || row_line.starts_with("loop_")
                        || row_line.starts_with('_')
                    {
                        break;
                    }
                    let toks = tokenize_mmcif_row(row_line);
                    i += 1;
                    if toks.len() < cols.len() {
                        continue;
                    }
                    // Only load model 1.
                    if let Some(cm) = c_model {
                        if toks[cm] != "1" && toks[cm] != "." {
                            continue;
                        }
                    }
                    if let Some(ca) = c_alt {
                        if !matches!(toks[ca], "." | "?" | "A") {
                            continue;
                        }
                    }
                    let group = c_group.map_or("ATOM", |idx| toks[idx]);
                    let is_hetatm = group == "HETATM";
                    let res_name = c_comp.map_or("UNK", |idx| toks[idx]);
                    if is_hetatm && matches!(res_name, "HOH" | "WAT") {
                        continue;
                    }
                    let elem = c_elem.map_or("C", |idx| toks[idx]);
                    if elem.eq_ignore_ascii_case("H") {
                        continue;
                    }

                    let (Some(ix), Some(iy), Some(iz)) = (c_x, c_y, c_z) else {
                        continue;
                    };
                    let (Ok(x), Ok(y), Ok(z)) = (
                        toks[ix].parse::<f32>(),
                        toks[iy].parse::<f32>(),
                        toks[iz].parse::<f32>(),
                    ) else {
                        continue;
                    };
                    let b_factor = c_bfac
                        .and_then(|idx| toks[idx].parse::<f32>().ok())
                        .unwrap_or(50.0);
                    let atom_name = c_atom.map_or("CA", |idx| toks[idx]);
                    let chain_id = c_chain
                        .and_then(|idx| toks[idx].chars().next())
                        .unwrap_or('A');
                    let seq_id = c_seq
                        .and_then(|idx| toks[idx].parse::<i32>().ok())
                        .unwrap_or(0);

                    let name_bytes = copy_ascii_fixed::<4>(atom_name);
                    let elem_bytes = copy_ascii_fixed::<2>(elem);
                    let res_bytes = copy_ascii_fixed::<3>(res_name);
                    let is_ca = !is_hetatm && atom_name == "CA";
                    let atom_idx = atoms.len() as u32;

                    let residue_idx = if is_hetatm {
                        u32::MAX
                    } else {
                        let key = (chain_id, seq_id, res_bytes);
                        if current_res_key != Some(key) {
                            current_res_key = Some(key);
                            residues.push(Residue3D {
                                atom_start: atom_idx,
                                ca_atom: if is_ca { Some(atom_idx) } else { None },
                                seq_id,
                                plddt: b_factor,
                                chain_id,
                                atom_count: 1,
                                name: res_bytes,
                                one_letter: three_to_one_aa(res_name),
                                secondary_structure: SecondaryStructure::Coil,
                            });
                        } else if let Some(last) = residues.last_mut() {
                            last.atom_count = last.atom_count.saturating_add(1);
                            if is_ca {
                                last.ca_atom = Some(atom_idx);
                                last.plddt = b_factor;
                            }
                        }
                        (residues.len().saturating_sub(1)) as u32
                    };

                    atoms.push(Atom {
                        pos: Vec3::new(x, y, z),
                        b_factor,
                        residue_idx,
                        name: name_bytes,
                        element: elem_bytes,
                        is_hetatm,
                        is_ca,
                    });
                }
                continue;
            } else if is_struct_conf {
                let c_type = find_col(&cols, "_struct_conf.conf_type_id");
                let c_chain = find_col(&cols, "_struct_conf.beg_auth_asym_id")
                    .or_else(|| find_col(&cols, "_struct_conf.beg_label_asym_id"));
                let c_beg = find_col(&cols, "_struct_conf.beg_auth_seq_id")
                    .or_else(|| find_col(&cols, "_struct_conf.beg_label_seq_id"));
                let c_end = find_col(&cols, "_struct_conf.end_auth_seq_id")
                    .or_else(|| find_col(&cols, "_struct_conf.end_label_seq_id"));

                while i < lines.len() {
                    let row_line = lines[i].trim();
                    if row_line.is_empty()
                        || row_line.starts_with('#')
                        || row_line.starts_with("loop_")
                        || row_line.starts_with('_')
                    {
                        break;
                    }
                    let toks = tokenize_mmcif_row(row_line);
                    i += 1;
                    if toks.len() < cols.len() {
                        continue;
                    }
                    let conf_type = c_type.map_or("", |idx| toks[idx]);
                    let kind = if conf_type.starts_with("HELX") {
                        SecondaryStructure::Helix
                    } else if conf_type.starts_with("STRN") {
                        SecondaryStructure::Sheet
                    } else {
                        continue;
                    };
                    let chain = c_chain
                        .and_then(|idx| toks[idx].chars().next())
                        .unwrap_or('A');
                    if let (Some(s), Some(e)) = (
                        c_beg.and_then(|idx| toks[idx].parse::<i32>().ok()),
                        c_end.and_then(|idx| toks[idx].parse::<i32>().ok()),
                    ) {
                        ss_ranges.push(SsRange {
                            chain,
                            start_seq: s,
                            end_seq: e,
                            kind,
                        });
                    }
                }
                continue;
            } else if is_sheet_range {
                let c_chain = find_col(&cols, "_struct_sheet_range.beg_auth_asym_id")
                    .or_else(|| find_col(&cols, "_struct_sheet_range.beg_label_asym_id"));
                let c_beg = find_col(&cols, "_struct_sheet_range.beg_auth_seq_id")
                    .or_else(|| find_col(&cols, "_struct_sheet_range.beg_label_seq_id"));
                let c_end = find_col(&cols, "_struct_sheet_range.end_auth_seq_id")
                    .or_else(|| find_col(&cols, "_struct_sheet_range.end_label_seq_id"));

                while i < lines.len() {
                    let row_line = lines[i].trim();
                    if row_line.is_empty()
                        || row_line.starts_with('#')
                        || row_line.starts_with("loop_")
                        || row_line.starts_with('_')
                    {
                        break;
                    }
                    let toks = tokenize_mmcif_row(row_line);
                    i += 1;
                    if toks.len() < cols.len() {
                        continue;
                    }
                    let chain = c_chain
                        .and_then(|idx| toks[idx].chars().next())
                        .unwrap_or('A');
                    if let (Some(s), Some(e)) = (
                        c_beg.and_then(|idx| toks[idx].parse::<i32>().ok()),
                        c_end.and_then(|idx| toks[idx].parse::<i32>().ok()),
                    ) {
                        ss_ranges.push(SsRange {
                            chain,
                            start_seq: s,
                            end_seq: e,
                            kind: SecondaryStructure::Sheet,
                        });
                    }
                }
                continue;
            }
        }
        i += 1;
    }

    finalize_structure(
        entry_id.into(),
        title,
        accession,
        source_label.into(),
        is_alphafold,
        atoms,
        residues,
        &ss_ranges,
    )
}

/// Applies secondary structure records (or geometric C-alpha inference when none are present)
/// and computes bounding sphere + mean pLDDT.
#[allow(clippy::too_many_arguments)]
fn finalize_structure(
    entry_id: String,
    title: String,
    accession: Option<String>,
    source_label: String,
    is_alphafold: bool,
    atoms: Vec<Atom>,
    mut residues: Vec<Residue3D>,
    ss_ranges: &[SsRange],
) -> Result<ProteinStructure> {
    if atoms.is_empty() || residues.is_empty() {
        return Err(HelixError::Structure(
            "structure file contains no polymer ATOM records".to_owned(),
        ));
    }

    if ss_ranges.is_empty() {
        infer_secondary_structure_from_ca(&atoms, &mut residues);
    } else {
        for res in &mut residues {
            for r in ss_ranges {
                if res.chain_id == r.chain && res.seq_id >= r.start_seq && res.seq_id <= r.end_seq {
                    res.secondary_structure = r.kind;
                    break;
                }
            }
        }
    }

    // Compute centroid and bounding radius over C-alpha positions (fallback: all atoms).
    let mut sum = Vec3::ZERO;
    let mut count = 0usize;
    for res in &residues {
        if let Some(ca) = res.ca_atom {
            if let Some(atom) = atoms.get(ca as usize) {
                sum = sum.add(atom.pos);
                count += 1;
            }
        }
    }
    if count == 0 {
        for atom in &atoms {
            sum = sum.add(atom.pos);
            count += 1;
        }
    }
    let center = if count > 0 {
        sum.scale(1.0 / count as f32)
    } else {
        Vec3::ZERO
    };

    let mut max_dist_sq = 1.0f32;
    for atom in &atoms {
        let d_sq = atom.pos.sub(center).length_sq();
        if d_sq > max_dist_sq {
            max_dist_sq = d_sq;
        }
    }
    let radius = max_dist_sq.sqrt().max(5.0);

    let mean_plddt = residues.iter().map(|r| r.plddt).sum::<f32>() / residues.len() as f32;
    let display_title = if title.is_empty() {
        entry_id.clone()
    } else {
        title
    };

    Ok(ProteinStructure {
        title: display_title,
        accession,
        entry_id,
        organism: None,
        source_label,
        is_alphafold,
        atoms: atoms.into_boxed_slice(),
        residues: residues.into_boxed_slice(),
        center,
        radius,
        mean_plddt,
    })
}

/// Fast geometric C-alpha secondary structure assignment (P-SEA inspired) used when PDB files
/// (such as AlphaFold v4 PDB exports) omit explicit `HELIX` and `SHEET` header records.
fn infer_secondary_structure_from_ca(atoms: &[Atom], residues: &mut [Residue3D]) {
    let n = residues.len();
    if n < 5 {
        return;
    }

    let ca_pos = |idx: usize| -> Option<Vec3> {
        let ca = residues.get(idx)?.ca_atom?;
        atoms.get(ca as usize).map(|a| a.pos)
    };

    let mut assigned = vec![SecondaryStructure::Coil; n];

    // Step 1: Identify right-handed alpha-helices via characteristic Cα(i..i+2), Cα(i..i+3),
    // and Cα(i..i+4) distances (~5.4 Å, ~5.1 Å, ~6.2 Å).
    for i in 0..n.saturating_sub(4) {
        if residues[i].chain_id != residues[i + 4].chain_id {
            continue;
        }
        let (Some(p0), Some(p1), Some(p2), Some(p3), Some(p4)) = (
            ca_pos(i),
            ca_pos(i + 1),
            ca_pos(i + 2),
            ca_pos(i + 3),
            ca_pos(i + 4),
        ) else {
            continue;
        };

        // Ensure backbone continuity (peptide Cα-Cα distance ~3.8 Å).
        if p0.distance(p1) > 4.4
            || p1.distance(p2) > 4.4
            || p2.distance(p3) > 4.4
            || p3.distance(p4) > 4.4
        {
            continue;
        }

        let d02 = p0.distance(p2);
        let d03 = p0.distance(p3);
        let d04 = p0.distance(p4);
        let d13 = p1.distance(p3);

        // Right-handed helix chirality check via scalar triple product.
        let v1 = p1.sub(p0);
        let v2 = p2.sub(p1);
        let v3 = p3.sub(p2);
        let chiral = v1.cross(v2).dot(v3);

        if (4.7..=6.2).contains(&d02)
            && (4.5..=6.0).contains(&d03)
            && (5.3..=7.2).contains(&d04)
            && (4.7..=6.2).contains(&d13)
            && chiral > 0.0
        {
            for slot in assigned.iter_mut().take(i + 5).skip(i) {
                *slot = SecondaryStructure::Helix;
            }
        }
    }

    // Step 2: Identify extended beta-strands via extended Cα(i..i+2) (~6.3–7.4 Å) and
    // Cα(i..i+3) (~9.1–11.2 Å) geometry.
    for i in 0..n.saturating_sub(3) {
        if residues[i].chain_id != residues[i + 3].chain_id {
            continue;
        }
        if assigned[i] == SecondaryStructure::Helix {
            continue;
        }
        let (Some(p0), Some(p1), Some(p2), Some(p3)) =
            (ca_pos(i), ca_pos(i + 1), ca_pos(i + 2), ca_pos(i + 3))
        else {
            continue;
        };
        if p0.distance(p1) > 4.4 || p1.distance(p2) > 4.4 || p2.distance(p3) > 4.4 {
            continue;
        }
        let d02 = p0.distance(p2);
        let d13 = p1.distance(p3);
        let d03 = p0.distance(p3);

        if (6.0..=7.5).contains(&d02) && (6.0..=7.5).contains(&d13) && (8.8..=11.4).contains(&d03)
        {
            for slot in assigned.iter_mut().take(i + 4).skip(i) {
                if *slot == SecondaryStructure::Coil {
                    *slot = SecondaryStructure::Sheet;
                }
            }
        }
    }

    for (res, ss) in residues.iter_mut().zip(assigned) {
        res.secondary_structure = ss;
    }
}

/// Bidirectional coordinate map linking a 3D [`ProteinStructure`] to an [`Alignment`] and
/// optional source [`Annotation`].
#[derive(Debug, Clone)]
pub struct StructureAlignmentMap {
    /// For each 3D residue index `r`, the linked alignment column span `(start_col, end_col)`.
    /// (`end_col - start_col == 1` for protein alignments, `>= 3` for translated DNA CDS).
    pub residue_to_cols: Vec<Option<(usize, usize)>>,
    /// Per-3D-residue conservation percentage (`0.0..=100.0`) sampled from the MSA columns.
    pub residue_conservation: Vec<Option<f32>>,
    /// Number of 3D residues successfully mapped onto alignment columns.
    pub mapped_count: usize,
}

impl StructureAlignmentMap {
    /// Builds a coordinate map connecting `structure` to `alignment` via `annotation`.
    ///
    /// - For **protein** alignments: matches the annotated sequence span against the 3D
    ///   structure's 1-letter residue sequence (falling back to UniProt residue offset when
    ///   sequences diverge).
    /// - For **DNA/RNA** alignments: translates the annotated CDS/Gene region (respecting
    ///   strand orientation) and maps each codon's 3-column span to its corresponding 3D amino acid.
    #[must_use]
    pub fn build(
        structure: &ProteinStructure,
        alignment: &Alignment,
        annotation: Option<&Annotation>,
    ) -> Self {
        let res_count = structure.residues.len();
        let mut residue_to_cols = vec![None; res_count];
        let mut residue_conservation = vec![None; res_count];

        let row_idx = annotation
            .and_then(|a| a.target_sequence)
            .unwrap_or(0)
            .min(alignment.row_count().saturating_sub(1));

        let Some(seq) = alignment.sequences().get(row_idx) else {
            return Self {
                residue_to_cols,
                residue_conservation,
                mapped_count: 0,
            };
        };

        let (span_start, span_end) = match annotation {
            Some(ann) => (ann.start, ann.end.min(alignment.width())),
            None => (0, alignment.width()),
        };

        // Extract (amino_acid, col_start, col_end) tuples from the alignment.
        let aln_residues: Vec<(u8, usize, usize)> = if alignment.alphabet.is_nucleotide() {
            let strand = annotation.map(|a| a.strand).unwrap_or_default();
            translate_aligned_region(&seq.residues, span_start, span_end, strand, 0)
                .into_iter()
                .map(|c| (c.amino_acid, c.start_col, c.end_col))
                .collect()
        } else {
            let mut list = Vec::with_capacity(span_end.saturating_sub(span_start));
            for col in span_start..span_end {
                let b = seq.residues[col];
                if !is_gap(b) {
                    list.push((b, col, col + 1));
                }
            }
            list
        };

        if aln_residues.is_empty() {
            return Self {
                residue_to_cols,
                residue_conservation,
                mapped_count: 0,
            };
        }

        // Strategy 1: Slide the alignment amino acid sequence across the 3D structure sequence
        // to find the highest-identity ungapped anchor offset (handles domain sub-fragments cleanly!).
        let struct_seq = structure.one_letter_sequence();
        let best_offset = find_best_sequence_offset(&aln_residues, &struct_seq);

        let uniprot_hint_offset = annotation
            .and_then(|a| a.uniprot.as_ref())
            .map(|u| u.uniprot_start_residue.saturating_sub(1) as isize);

        // Choose the sequence-matched offset if identity >= 45%, otherwise respect explicit UniProt offset.
        let offset: isize = match best_offset {
            Some((off, identity)) if identity >= 0.45 => off,
            _ => uniprot_hint_offset.unwrap_or(0),
        };

        let mut mapped_count = 0usize;
        for (i, &(_aa, col_start, col_end)) in aln_residues.iter().enumerate() {
            let struct_idx = i as isize + offset;
            if struct_idx >= 0 && (struct_idx as usize) < res_count {
                let s_idx = struct_idx as usize;
                residue_to_cols[s_idx] = Some((col_start, col_end));
                // Mean conservation across the column(s) belonging to this residue/codon.
                let col_slice = &alignment.columns()[col_start..col_end.min(alignment.width())];
                if !col_slice.is_empty() {
                    let avg_id = col_slice.iter().map(|c| c.identity_pct).sum::<f32>()
                        / col_slice.len() as f32;
                    residue_conservation[s_idx] = Some(avg_id);
                }
                mapped_count += 1;
            }
        }

        Self {
            residue_to_cols,
            residue_conservation,
            mapped_count,
        }
    }

    /// Finds the 3D residue index mapped to alignment column `col`, if any.
    #[must_use]
    pub fn residue_for_column(&self, col: usize) -> Option<usize> {
        self.residue_to_cols.iter().position(|entry| {
            entry.is_some_and(|(start, end)| col >= start && col < end)
        })
    }
}

/// Finds the relative offset `struct_idx = aln_idx + offset` maximizing residue identity.
fn find_best_sequence_offset(
    aln_residues: &[(u8, usize, usize)],
    struct_seq: &[u8],
) -> Option<(isize, f32)> {
    let a_len = aln_residues.len();
    let s_len = struct_seq.len();
    if a_len == 0 || s_len == 0 {
        return None;
    }

    // Require at least 8 overlapping residues (or the shorter length if < 8).
    let min_overlap = a_len.min(s_len).min(8);
    let min_offset = -(a_len as isize) + min_overlap as isize;
    let max_offset = (s_len as isize) - min_overlap as isize;

    let mut best_score = 0usize;
    let mut best_id = 0.0f32;
    let mut best_off = 0isize;

    for offset in min_offset..=max_offset {
        let a_start = (-offset).max(0) as usize;
        let s_start = offset.max(0) as usize;
        let overlap = (a_len - a_start).min(s_len - s_start);
        if overlap < min_overlap {
            continue;
        }
        let mut matches = 0usize;
        for k in 0..overlap {
            if aln_residues[a_start + k].0 == struct_seq[s_start + k] {
                matches += 1;
            }
        }
        let id = matches as f32 / overlap as f32;
        if matches > best_score || (matches == best_score && id > best_id) {
            best_score = matches;
            best_id = id;
            best_off = offset;
        }
    }

    if best_score == 0 {
        None
    } else {
        Some((best_off, best_id))
    }
}

fn find_col(cols: &[&str], name: &str) -> Option<usize> {
    cols.iter().position(|&c| c == name)
}

/// Splits an mmCIF data row into tokens, respecting single- and double-quoted strings.
fn tokenize_mmcif_row(line: &str) -> Vec<&str> {
    let mut tokens = Vec::with_capacity(24);
    let bytes = line.as_bytes();
    let mut i = 0usize;
    let len = bytes.len();

    while i < len {
        while i < len && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len {
            break;
        }
        if bytes[i] == b'\'' || bytes[i] == b'"' {
            let quote = bytes[i];
            i += 1;
            let start = i;
            while i < len && bytes[i] != quote {
                i += 1;
            }
            tokens.push(&line[start..i]);
            if i < len {
                i += 1;
            }
        } else {
            let start = i;
            while i < len && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            tokens.push(&line[start..i]);
        }
    }
    tokens
}

fn parse_trimmed_f32(slice: Option<&str>) -> Option<f32> {
    slice?.trim().parse::<f32>().ok()
}

fn parse_trimmed_i32(slice: Option<&str>) -> Option<i32> {
    slice?.trim().parse::<i32>().ok()
}

fn copy_ascii_fixed<const N: usize>(s: &str) -> [u8; N] {
    let mut buf = [b' '; N];
    for (slot, b) in buf.iter_mut().zip(s.bytes()) {
        *slot = b.to_ascii_uppercase();
    }
    buf
}

fn infer_element_from_atom_name(atom_name: &str) -> &str {
    let trimmed = atom_name.trim_start_matches(|c: char| c.is_ascii_digit());
    if trimmed.starts_with("ZN") {
        "ZN"
    } else if trimmed.starts_with("FE") {
        "FE"
    } else if trimmed.starts_with("MG") {
        "MG"
    } else if trimmed.starts_with("CA") && trimmed.len() == 2 {
        "C"
    } else {
        trimmed.get(0..1).unwrap_or("C")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{Alphabet, Sequence};

    const MINI_PDB: &str = "\
TITLE     MINI HELIX TEST
HELIX    1   1 ALA A    1  GLY A    4  1                                   4
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 94.20           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 95.10           C
ATOM      3  CA  VAL A   2       3.800   2.100   1.200  1.00 88.50           C
ATOM      4  CA  LEU A   3       6.100   0.400   2.900  1.00 72.30           C
ATOM      5  CA  GLY A   4       8.400   2.500   4.100  1.00 61.00           C
END
";

    #[test]
    fn pdb_parser_extracts_residues_plddt_and_secondary_structure() {
        let st = parse_pdb(MINI_PDB, "TEST-1", Some("P00001".to_owned()), "Unit", true)
            .expect("parse mini pdb");
        assert_eq!(st.residues.len(), 4);
        assert_eq!(st.one_letter_sequence(), b"AVLG");
        assert_eq!(st.residues[0].secondary_structure, SecondaryStructure::Helix);
        assert!((st.residues[0].plddt - 95.1).abs() < 0.01);
        assert!(st.radius >= 5.0);
    }

    #[test]
    fn structure_alignment_map_aligns_dna_cds_codons_to_3d_residues() {
        let st = parse_pdb(MINI_PDB, "TEST-1", None, "Unit", true).expect("parse pdb");
        // DNA coding for A-V-L-G: GCT GTG CTG GGT
        let aln = Alignment::new(
            "DNA",
            Alphabet::Dna,
            vec![Sequence::new("Seq1", b"GCTGTGCTGGGT")],
        )
        .expect("aln");
        let map = StructureAlignmentMap::build(&st, &aln, None);
        assert_eq!(map.mapped_count, 4);
        assert_eq!(map.residue_to_cols[0], Some((0, 3)));
        assert_eq!(map.residue_to_cols[1], Some((3, 6)));
        assert_eq!(map.residue_for_column(4), Some(1));
    }
}
