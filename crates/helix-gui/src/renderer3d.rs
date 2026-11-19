//! Pure-Rust GPU-accelerated 3D protein renderer built on `egui::Mesh`.
//!
//! Renders smooth Catmull-Rom cartoon ribbons (coiled alpha-helices, arrowed beta-strands,
//! tubular loops), backbone tubes, C-alpha traces, ball-and-stick, and spacefill models
//! with Blinn-Phong lighting, depth fog, painter's Z-sorting, and interactive residue picking.
//! Zero `unsafe` code and 100% compatible with native separate OS viewports in `eframe`.

use crate::theme::{conservation_color, cpk_element_color, plddt_color, rainbow_color};
use egui::{Color32, Mesh, Painter, Pos2, Rect, Shape, Stroke};
use helix_core::{ProteinStructure, SecondaryStructure, StructureAlignmentMap, Vec3};
use std::ops::{Add, Sub};

/// 3D molecular representation style.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RenderStyle {
    /// Smooth spline ribbon with broad coiled helices, flat arrowed beta-strands, and loop tubes.
    #[default]
    Cartoon,
    /// Continuous cylindrical tube along the C-alpha spline.
    BackboneTube,
    /// Direct C-alpha bonds with residue spheres.
    CaTrace,
    /// Heavy-atom covalent sticks and spheres (great for inspecting active sites).
    BallAndStick,
    /// Shaded van der Waals space-filling atomic spheres.
    Spacefill,
}

impl RenderStyle {
    pub const ALL: [Self; 5] = [
        Self::Cartoon,
        Self::BackboneTube,
        Self::CaTrace,
        Self::BallAndStick,
        Self::Spacefill,
    ];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Cartoon => "Cartoon Ribbon",
            Self::BackboneTube => "Backbone Tube",
            Self::CaTrace => "C-Alpha Trace",
            Self::BallAndStick => "Ball & Stick",
            Self::Spacefill => "Spacefill (VDW)",
        }
    }
}

/// 3D color mapping strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ColorBy3D {
    /// Official DeepMind/EBI AlphaFold pLDDT confidence coloring.
    #[default]
    AlphaFoldPlddt,
    /// Secondary structure (Alpha Helix, Beta Strand, Coil).
    SecondaryStructure,
    /// Pairwise conservation % projected directly from the active MSA alignment!
    AlignmentConservation,
    /// N-terminus (blue) to C-terminus (red) spectral gradient.
    Rainbow,
    /// Chain color (or CPK element color in atomic modes).
    ChainOrElement,
}

impl ColorBy3D {
    pub const ALL: [Self; 5] = [
        Self::AlphaFoldPlddt,
        Self::SecondaryStructure,
        Self::AlignmentConservation,
        Self::Rainbow,
        Self::ChainOrElement,
    ];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::AlphaFoldPlddt => "AlphaFold pLDDT",
            Self::SecondaryStructure => "Secondary Structure",
            Self::AlignmentConservation => "MSA Conservation",
            Self::Rainbow => "Rainbow (N -> C)",
            Self::ChainOrElement => "Chain / Element",
        }
    }
}

/// Orbital 3D camera with Euler yaw/pitch, screen-space pan, and zoom factor.
#[derive(Debug, Clone)]
pub struct Camera3D {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub pan_x: f32,
    pub pan_y: f32,
    pub auto_rotate: bool,
}

impl Default for Camera3D {
    fn default() -> Self {
        Self {
            yaw: 0.55,
            pitch: 0.25,
            zoom: 1.0,
            pan_x: 0.0,
            pan_y: 0.0,
            auto_rotate: false,
        }
    }
}

impl Camera3D {
    /// Resets orbit, zoom, and pan to the default framing angle.
    pub fn reset(&mut self) {
        self.yaw = 0.55;
        self.pitch = 0.25;
        self.zoom = 1.0;
        self.pan_x = 0.0;
        self.pan_y = 0.0;
    }

    /// Rotates a model-space point (relative to protein centroid) into view space
    /// where `+Z` points toward the viewer.
    #[inline]
    fn to_view_space(&self, p: Vec3) -> Vec3 {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();

        // Yaw around Y axis.
        let x1 = p.x * cy + p.z * sy;
        let z1 = -p.x * sy + p.z * cy;

        // Pitch around X axis.
        let y2 = p.y * cp - z1 * sp;
        let z2 = p.y * sp + z1 * cp;

        Vec3::new(x1, y2, z2)
    }
}

/// Depth-sortable 3D shaded triangle primitive.
#[derive(Clone, Copy)]
struct Tri3D {
    p0: Pos2,
    p1: Pos2,
    p2: Pos2,
    depth: f32,
    color: Color32,
}

/// Reusable scratch buffers for zero-churn 60 FPS rendering and residue picking.
#[derive(Default)]
pub struct Renderer3D {
    tris: Vec<Tri3D>,
    picked_residues: Vec<(usize, Pos2, f32)>,
}

/// Visual highlight parameters synchronized with the main alignment window.
pub struct RenderHighlights {
    /// Currently hovered 3D residue index (either from 3D hover or alignment hover).
    pub hovered_residue: Option<usize>,
    /// Currently selected 3D residue index.
    pub selected_residue: Option<usize>,
    /// Optional residue range `[start, end)` corresponding to the active annotation.
    pub annotation_residue_span: Option<(usize, usize)>,
}

/// Output of a single rendered frame, including the residue under the pointer if any.
pub struct FramePickResult {
    pub hovered_residue: Option<usize>,
}

impl Renderer3D {
    /// Renders `structure` inside `rect` using `painter`, returning the residue index
    /// under `pointer_pos` (if hovering near the 3D backbone).
    #[allow(clippy::too_many_arguments)]
    pub fn render(
        &mut self,
        painter: &Painter,
        rect: Rect,
        structure: &ProteinStructure,
        map: Option<&StructureAlignmentMap>,
        camera: &Camera3D,
        style: RenderStyle,
        color_by: ColorBy3D,
        highlights: &RenderHighlights,
        pointer_pos: Option<Pos2>,
    ) -> FramePickResult {
        self.tris.clear();
        self.picked_residues.clear();

        // Perspective projection scaling from Angstroms to screen pixels.
        let min_dim = rect.width().min(rect.height()).max(64.0);
        let cam_dist = structure.radius * 3.2;
        let focal = min_dim * 0.44 * camera.zoom * (cam_dist / structure.radius.max(1.0));
        let center_px = rect.center() + egui::vec2(camera.pan_x, camera.pan_y);

        let project = |world: Vec3| -> (Pos2, Vec3, f32) {
            let rel = world.sub(structure.center);
            let v = camera.to_view_space(rel);
            let denom = (cam_dist - v.z).max(2.0);
            let scale = focal / denom;
            let sx = center_px.x + v.x * scale;
            let sy = center_px.y - v.y * scale;
            (Pos2::new(sx, sy), v, scale)
        };

        // Record projected C-alpha positions for accurate mouse picking and highlight halos.
        self.picked_residues.reserve(structure.residues.len());
        for (res_idx, _res) in structure.residues.iter().enumerate() {
            if let Some(ca_world) = structure.residue_ca_pos(res_idx) {
                let (pos2, v, _scale) = project(ca_world);
                self.picked_residues.push((res_idx, pos2, v.z));
            }
        }

        // Build 3D geometry according to selected molecular representation.
        match style {
            RenderStyle::Cartoon => {
                self.emit_cartoon_or_tube(structure, map, color_by, highlights, &project, true);
            }
            RenderStyle::BackboneTube => {
                self.emit_cartoon_or_tube(structure, map, color_by, highlights, &project, false);
            }
            RenderStyle::CaTrace => {
                self.emit_ca_trace(structure, map, color_by, highlights, &project);
            }
            RenderStyle::BallAndStick => {
                self.emit_atomic(structure, map, color_by, highlights, &project, false);
            }
            RenderStyle::Spacefill => {
                self.emit_atomic(structure, map, color_by, highlights, &project, true);
            }
        }

        // Always render HETATM ligands / metal ions (e.g., Heme Fe or Zn2+) as golden/CPK spheres.
        for atom in structure.atoms.iter() {
            if atom.is_hetatm {
                let base_color = cpk_element_color(atom.element_str());
                self.emit_shaded_sphere(atom.pos, 1.35, base_color, &project, 10);
            }
        }

        // Render vivid highlight sphere on the selected or hovered residue.
        if let Some(sel_idx) = highlights.selected_residue {
            if let Some(pos) = structure.residue_ca_pos(sel_idx) {
                self.emit_shaded_sphere(
                    pos,
                    1.65,
                    Color32::from_rgb(0, 230, 215),
                    &project,
                    12,
                );
            }
        }
        if let Some(hov_idx) = highlights.hovered_residue {
            if highlights.selected_residue != Some(hov_idx) {
                if let Some(pos) = structure.residue_ca_pos(hov_idx) {
                    self.emit_shaded_sphere(
                        pos,
                        1.45,
                        Color32::from_rgb(255, 235, 80),
                        &project,
                        10,
                    );
                }
            }
        }

        // Painter's algorithm sort: back-to-front by view-space Z (smaller Z = farther from viewer).
        self.tris
            .sort_unstable_by(|a, b| a.depth.total_cmp(&b.depth));

        // Pack sorted triangles into a single GPU-rendered egui::Mesh (batched draw call).
        let mut mesh = Mesh::default();
        mesh.vertices.reserve(self.tris.len() * 3);
        mesh.indices.reserve(self.tris.len() * 3);
        for tri in &self.tris {
            let base = mesh.vertices.len() as u32;
            mesh.vertices.push(egui::epaint::Vertex {
                pos: tri.p0,
                uv: egui::epaint::WHITE_UV,
                color: tri.color,
            });
            mesh.vertices.push(egui::epaint::Vertex {
                pos: tri.p1,
                uv: egui::epaint::WHITE_UV,
                color: tri.color,
            });
            mesh.vertices.push(egui::epaint::Vertex {
                pos: tri.p2,
                uv: egui::epaint::WHITE_UV,
                color: tri.color,
            });
            mesh.indices.push(base);
            mesh.indices.push(base + 1);
            mesh.indices.push(base + 2);
        }
        painter.add(Shape::mesh(mesh));

        // Draw crisp screen-space indicator ring around selected/hovered residue on top of mesh.
        if let Some(sel_idx) = highlights.selected_residue {
            if let Some(&(_, pos2, _)) = self.picked_residues.iter().find(|r| r.0 == sel_idx) {
                painter.circle_stroke(pos2, 12.0, Stroke::new(2.0, Color32::from_rgb(45, 225, 210)));
            }
        }

        // Screen-space residue picking: find front-most C-alpha within 16px of pointer.
        let mut best_pick: Option<(usize, f32)> = None;
        if let Some(ptr) = pointer_pos {
            if rect.contains(ptr) {
                let pick_radius_sq = 16.0f32 * 16.0f32;
                for &(res_idx, screen_pos, z) in &self.picked_residues {
                    let dist_sq = screen_pos.distance_sq(ptr);
                    if dist_sq <= pick_radius_sq {
                        let replace = match best_pick {
                            None => true,
                            Some((_, best_z)) => z > best_z,
                        };
                        if replace {
                            best_pick = Some((res_idx, z));
                        }
                    }
                }
            }
        }

        FramePickResult {
            hovered_residue: best_pick.map(|(idx, _)| idx),
        }
    }

    /// Generates smooth Catmull-Rom cartoon ribbons (helices, arrowed beta-sheets, coil tubes)
    /// or uniform backbone tubes across all polymer chains.
    fn emit_cartoon_or_tube<F>(
        &mut self,
        structure: &ProteinStructure,
        map: Option<&StructureAlignmentMap>,
        color_by: ColorBy3D,
        highlights: &RenderHighlights,
        project: &F,
        cartoon_shapes: bool,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        let res_count = structure.residues.len();
        if res_count < 2 {
            return;
        }

        // Adaptive spline smoothness so both 150-residue and 1,500-residue proteins stay 60 FPS.
        let subdivs: usize = if res_count > 700 { 2 } else { 4 };

        // Partition residues into continuous polymer chain segments (break on chain change or gap > 5.0 Å).
        let mut chain_start = 0usize;
        while chain_start < res_count {
            let mut chain_end = chain_start + 1;
            while chain_end < res_count {
                let prev = &structure.residues[chain_end - 1];
                let curr = &structure.residues[chain_end];
                if prev.chain_id != curr.chain_id {
                    break;
                }
                let broken = match (
                    structure.residue_ca_pos(chain_end - 1),
                    structure.residue_ca_pos(chain_end),
                ) {
                    (Some(a), Some(b)) => a.distance(b) > 5.2,
                    _ => true,
                };
                if broken {
                    break;
                }
                chain_end += 1;
            }

            self.emit_chain_ribbon(
                structure,
                map,
                color_by,
                highlights,
                project,
                chain_start,
                chain_end,
                subdivs,
                cartoon_shapes,
            );
            chain_start = chain_end;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn emit_chain_ribbon<F>(
        &mut self,
        structure: &ProteinStructure,
        map: Option<&StructureAlignmentMap>,
        color_by: ColorBy3D,
        highlights: &RenderHighlights,
        project: &F,
        start: usize,
        end: usize,
        subdivs: usize,
        cartoon_shapes: bool,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        let seg_len = end - start;
        if seg_len < 2 {
            return;
        }

        // Gather C-alpha control points for this continuous chain segment.
        let mut pts: Vec<Vec3> = Vec::with_capacity(seg_len);
        for idx in start..end {
            if let Some(p) = structure.residue_ca_pos(idx) {
                pts.push(p);
            }
        }
        if pts.len() < 2 {
            return;
        }

        // Build spline samples with smooth tangent/normal transport frames (parallel transport
        // avoids ribbon twisting flips at inflection points).
        struct SplineSample {
            pos: Vec3,
            normal: Vec3,
            binormal: Vec3,
            res_idx: usize,
            half_width: f32,
            half_thick: f32,
            color: Color32,
        }

        let total_samples = (pts.len() - 1) * subdivs + 1;
        let mut samples: Vec<SplineSample> = Vec::with_capacity(total_samples);

        let mut prev_normal = Vec3::new(0.0, 1.0, 0.0);

        for i in 0..(pts.len() - 1) {
            let p0 = pts[i.saturating_sub(1)];
            let p1 = pts[i];
            let p2 = pts[i + 1];
            let p3 = pts[(i + 2).min(pts.len() - 1)];

            let res_a = start + i;
            let res_b = start + i + 1;
            let steps = if i == pts.len() - 2 {
                subdivs + 1
            } else {
                subdivs
            };

            for s in 0..steps {
                let t = s as f32 / subdivs as f32;
                let pos = catmull_rom(p0, p1, p2, p3, t);
                let tang = catmull_rom_tangent(p0, p1, p2, p3, t).normalize();

                // Parallel-transport previous normal onto the plane perpendicular to `tang`.
                let proj = prev_normal.sub(tang.scale(prev_normal.dot(tang)));
                let normal = if proj.length_sq() > 1e-4 {
                    proj.normalize()
                } else {
                    let fallback = if tang.y.abs() < 0.9 {
                        Vec3::new(0.0, 1.0, 0.0)
                    } else {
                        Vec3::new(1.0, 0.0, 0.0)
                    };
                    tang.cross(fallback).normalize()
                };
                prev_normal = normal;
                let binormal = tang.cross(normal).normalize();

                let res_idx = if t < 0.5 { res_a } else { res_b };
                let ss = structure.residues[res_idx].secondary_structure;

                // Determine cross-section dimensions: helices = broad ribbon, sheets = flat arrowed
                // ribbon, coils = sleek round tube.
                let (half_width, half_thick) = if !cartoon_shapes {
                    (0.55, 0.55)
                } else {
                    match ss {
                        SecondaryStructure::Helix => (1.18, 0.34),
                        SecondaryStructure::Sheet => {
                            // Flare into a crisp arrowhead on the final residue of a beta-strand.
                            let is_sheet_end = res_idx + 1 >= end
                                || structure.residues[res_idx + 1].secondary_structure
                                    != SecondaryStructure::Sheet;
                            if is_sheet_end {
                                let arrow_t = ((t - 0.1).max(0.0) / 0.9).clamp(0.0, 1.0);
                                let w = 1.85 * (1.0 - arrow_t) + 0.35 * arrow_t;
                                (w, 0.24)
                            } else {
                                (1.05, 0.24)
                            }
                        }
                        SecondaryStructure::Coil => (0.36, 0.36),
                    }
                };

                let color = residue_3d_color(structure, res_idx, map, color_by, highlights);

                samples.push(SplineSample {
                    pos,
                    normal,
                    binormal,
                    res_idx,
                    half_width,
                    half_thick,
                    color,
                });
            }
        }

        // Extrude a 4-sided chamfered ribbon prism between consecutive spline frames.
        const PROFILE: [(f32, f32); 4] = [
            (1.0, 0.0),
            (0.0, 1.0),
            (-1.0, 0.0),
            (0.0, -1.0),
        ];

        for pair in samples.windows(2) {
            let a = &pair[0];
            let b = &pair[1];

            let mut ring_a = [Vec3::ZERO; 4];
            let mut ring_b = [Vec3::ZERO; 4];

            for k in 0..4 {
                let (pu, pv) = PROFILE[k];
                ring_a[k] = a
                    .pos
                    .add(a.normal.scale(pu * a.half_width))
                    .add(a.binormal.scale(pv * a.half_thick));
                ring_b[k] = b
                    .pos
                    .add(b.normal.scale(pu * b.half_width))
                    .add(b.binormal.scale(pv * b.half_thick));
            }

            // Use slightly brighter highlight tint if this segment belongs to the selected residue.
            let seg_color = if highlights.selected_residue == Some(a.res_idx) {
                Color32::from_rgb(0, 224, 210)
            } else {
                a.color
            };

            for k in 0..4 {
                let next = (k + 1) % 4;
                self.emit_lit_quad(
                    ring_a[k],
                    ring_a[next],
                    ring_b[next],
                    ring_b[k],
                    seg_color,
                    project,
                    structure.radius,
                );
            }
        }
    }

    /// Renders direct C-alpha bond cylinders and residue joints.
    fn emit_ca_trace<F>(
        &mut self,
        structure: &ProteinStructure,
        map: Option<&StructureAlignmentMap>,
        color_by: ColorBy3D,
        highlights: &RenderHighlights,
        project: &F,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        for i in 0..structure.residues.len() {
            let color = residue_3d_color(structure, i, map, color_by, highlights);
            if let Some(p1) = structure.residue_ca_pos(i) {
                self.emit_shaded_sphere(p1, 0.62, color, project, 6);
                if i + 1 < structure.residues.len()
                    && structure.residues[i].chain_id == structure.residues[i + 1].chain_id
                {
                    if let Some(p2) = structure.residue_ca_pos(i + 1) {
                        if p1.distance(p2) <= 5.2 {
                            self.emit_cylinder(p1, p2, 0.34, color, project, structure.radius);
                        }
                    }
                }
            }
        }
    }

    /// Renders heavy-atom Ball-and-Stick or van der Waals Spacefill representations.
    fn emit_atomic<F>(
        &mut self,
        structure: &ProteinStructure,
        map: Option<&StructureAlignmentMap>,
        color_by: ColorBy3D,
        highlights: &RenderHighlights,
        project: &F,
        spacefill: bool,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        // Cap detail smoothly for massive >2,500-atom AlphaFold proteins so orbit stays responsive.
        let atom_count = structure.atoms.len();
        let stride = if atom_count > 4500 {
            3
        } else if atom_count > 2200 {
            2
        } else {
            1
        };

        for (res_idx, _res) in structure.residues.iter().enumerate() {
            let res_color = residue_3d_color(structure, res_idx, map, color_by, highlights);
            let res_atoms = structure.residue_atoms(res_idx);

            for (local_i, atom) in res_atoms.iter().enumerate() {
                if stride > 1 && !atom.is_ca && (local_i % stride != 0) {
                    continue;
                }
                let col = if color_by == ColorBy3D::ChainOrElement {
                    cpk_element_color(atom.element_str())
                } else {
                    res_color
                };
                let r = if spacefill {
                    vdw_radius(atom.element_str())
                } else {
                    0.46
                };
                let sides = if atom_count > 2000 { 6 } else { 8 };
                self.emit_shaded_sphere(atom.pos, r, col, project, sides);
            }

            // Covalent intra-residue bonds in Ball-and-Stick mode.
            if !spacefill && stride == 1 {
                for a_i in 0..res_atoms.len() {
                    for b_i in (a_i + 1)..res_atoms.len() {
                        let d_sq = res_atoms[a_i].pos.sub(res_atoms[b_i].pos).length_sq();
                        if d_sq > 0.4 && d_sq <= 3.1 {
                            self.emit_cylinder(
                                res_atoms[a_i].pos,
                                res_atoms[b_i].pos,
                                0.17,
                                res_color,
                                project,
                                structure.radius,
                            );
                        }
                    }
                }
            }
        }
    }

    /// Emits a 4-sided 3D cylinder connecting `a` and `b`.
    fn emit_cylinder<F>(
        &mut self,
        a: Vec3,
        b: Vec3,
        radius: f32,
        color: Color32,
        project: &F,
        depth_scale: f32,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        let axis = b.sub(a).normalize();
        let up = if axis.y.abs() < 0.9 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let u = axis.cross(up).normalize().scale(radius);
        let v = axis.cross(u).normalize().scale(radius);

        let offs = [u, v, u.scale(-1.0), v.scale(-1.0)];
        for k in 0..4 {
            let next = (k + 1) % 4;
            self.emit_lit_quad(
                a.add(offs[k]),
                a.add(offs[next]),
                b.add(offs[next]),
                b.add(offs[k]),
                color,
                project,
                depth_scale,
            );
        }
    }

    /// Emits a shaded 3D billboard sphere approximated by a radial fan with a specular highlight.
    ///
    /// Much faster than full 3D icospheres when drawing thousands of atoms while maintaining
    /// smooth spherical shading and depth-sorted compositing with ribbon triangles.
    fn emit_shaded_sphere<F>(
        &mut self,
        world_pos: Vec3,
        world_radius: f32,
        base_color: Color32,
        project: &F,
        segments: usize,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        let (center_2d, view_pos, px_per_angstrom) = project(world_pos);
        let r_px = (world_radius * px_per_angstrom).clamp(2.0, 90.0);

        // Offset highlight center slightly upper-left to mimic key light direction.
        let light_center = Pos2::new(center_2d.x - r_px * 0.28, center_2d.y - r_px * 0.28);

        let depth_fog = ((view_pos.z + 40.0) / 80.0).clamp(0.42, 1.08);
        let core_color = scale_color(base_color, 1.25 * depth_fog, 38.0);
        let rim_color = scale_color(base_color, 0.55 * depth_fog, 0.0);

        let segs = segments.max(6);
        let step = std::f32::consts::TAU / segs as f32;

        for i in 0..segs {
            let a0 = i as f32 * step;
            let a1 = (i + 1) as f32 * step;
            let p0 = Pos2::new(center_2d.x + a0.cos() * r_px, center_2d.y + a0.sin() * r_px);
            let p1 = Pos2::new(center_2d.x + a1.cos() * r_px, center_2d.y + a1.sin() * r_px);

            // Slight directional rim tint so upper-left facets catch more light than lower-right.
            let mid_ang = (a0 + a1) * 0.5;
            let directional = (-mid_ang.cos() * 0.65 - mid_ang.sin() * 0.65).clamp(-1.0, 1.0);
            let facet_color = if directional > 0.1 {
                lerp_color(rim_color, core_color, 0.45 + directional * 0.45)
            } else {
                lerp_color(rim_color, core_color, 0.25)
            };

            self.tris.push(Tri3D {
                p0: light_center,
                p1: p0,
                p2: p1,
                depth: view_pos.z + world_radius * 0.25,
                color: facet_color,
            });
        }
    }

    /// Projects a 3D quad into two Blinn-Phong shaded triangles with depth cueing.
    #[allow(clippy::too_many_arguments)]
    fn emit_lit_quad<F>(
        &mut self,
        w0: Vec3,
        w1: Vec3,
        w2: Vec3,
        w3: Vec3,
        base_color: Color32,
        project: &F,
        depth_scale: f32,
    ) where
        F: Fn(Vec3) -> (Pos2, Vec3, f32),
    {
        let (p0, v0, _) = project(w0);
        let (p1, v1, _) = project(w1);
        let (p2, v2, _) = project(w2);
        let (p3, v3, _) = project(w3);

        // Compute view-space face normal for Blinn-Phong shading.
        let e1 = v1.sub(v0);
        let e2 = v2.sub(v0);
        let mut normal = e1.cross(e2).normalize();
        // Two-sided lighting so ribbon undersides stay illuminated cleanly.
        if normal.z < 0.0 {
            normal = normal.scale(-1.0);
        }

        // Key light from upper-left-front + subtle rim light.
        let key_dir = Vec3::new(-0.42, 0.58, 0.70).normalize();
        let half_vec = key_dir.add(Vec3::new(0.0, 0.0, 1.0)).normalize();

        let diffuse = normal.dot(key_dir).clamp(0.0, 1.0);
        let spec = normal.dot(half_vec).clamp(0.0, 1.0).powf(18.0) * 0.35;
        let ambient = 0.30;

        // Depth fog gently darkens back-most segments for crisp 3D spatial perception.
        let avg_z = (v0.z + v1.z + v2.z + v3.z) * 0.25;
        let fog = ((avg_z + depth_scale * 1.3) / (depth_scale * 2.4)).clamp(0.48, 1.08);
        let shade = (ambient + diffuse * 0.78) * fog;
        let lit = scale_color(base_color, shade, spec * 180.0 * fog);

        self.tris.push(Tri3D {
            p0,
            p1,
            p2,
            depth: (v0.z + v1.z + v2.z) / 3.0,
            color: lit,
        });
        self.tris.push(Tri3D {
            p0,
            p1: p2,
            p2: p3,
            depth: (v0.z + v2.z + v3.z) / 3.0,
            color: lit,
        });
    }
}

/// Evaluates the color of 3D residue `res_idx` under the active coloring mode.
fn residue_3d_color(
    structure: &ProteinStructure,
    res_idx: usize,
    map: Option<&StructureAlignmentMap>,
    color_by: ColorBy3D,
    highlights: &RenderHighlights,
) -> Color32 {
    let res = &structure.residues[res_idx];
    let base = match color_by {
        ColorBy3D::AlphaFoldPlddt => plddt_color(res.plddt),
        ColorBy3D::SecondaryStructure => match res.secondary_structure {
            SecondaryStructure::Helix => Color32::from_rgb(230, 86, 134),
            SecondaryStructure::Sheet => Color32::from_rgb(240, 186, 46),
            SecondaryStructure::Coil => Color32::from_rgb(92, 168, 224),
        },
        ColorBy3D::AlignmentConservation => {
            let cons = map.and_then(|m| m.residue_conservation.get(res_idx).copied().flatten());
            conservation_color(cons)
        }
        ColorBy3D::Rainbow => {
            let denom = structure.residues.len().saturating_sub(1).max(1) as f32;
            rainbow_color(res_idx as f32 / denom)
        }
        ColorBy3D::ChainOrElement => chain_color(res.chain_id),
    };

    // Gently emphasize residues inside the selected alignment annotation span when viewing a sub-domain.
    if let Some((span_start, span_end)) = highlights.annotation_residue_span {
        if span_end > span_start && span_end - span_start < structure.residues.len() {
            if res_idx >= span_start && res_idx < span_end {
                return base;
            }
            return lerp_color(base, Color32::from_rgb(76, 86, 102), 0.62);
        }
    }

    base
}

fn chain_color(chain: char) -> Color32 {
    const CHAIN_PALETTE: [Color32; 6] = [
        Color32::from_rgb(72, 152, 232),
        Color32::from_rgb(78, 198, 134),
        Color32::from_rgb(236, 148, 72),
        Color32::from_rgb(186, 118, 232),
        Color32::from_rgb(236, 98, 128),
        Color32::from_rgb(72, 202, 204),
    ];
    let idx = (chain as usize).saturating_sub('A' as usize) % CHAIN_PALETTE.len();
    CHAIN_PALETTE[idx]
}

fn vdw_radius(element: &str) -> f32 {
    match element {
        "C" => 1.65,
        "N" => 1.55,
        "O" => 1.50,
        "S" => 1.80,
        "P" => 1.80,
        "ZN" | "FE" | "MG" => 1.40,
        _ => 1.55,
    }
}

/// Evaluates a standard Catmull-Rom cubic spline at `t in [0, 1]`.
fn catmull_rom(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
    let t2 = t * t;
    let t3 = t2 * t;
    let x = 0.5
        * ((2.0 * p1.x)
            + (-p0.x + p2.x) * t
            + (2.0 * p0.x - 5.0 * p1.x + 4.0 * p2.x - p3.x) * t2
            + (-p0.x + 3.0 * p1.x - 3.0 * p2.x + p3.x) * t3);
    let y = 0.5
        * ((2.0 * p1.y)
            + (-p0.y + p2.y) * t
            + (2.0 * p0.y - 5.0 * p1.y + 4.0 * p2.y - p3.y) * t2
            + (-p0.y + 3.0 * p1.y - 3.0 * p2.y + p3.y) * t3);
    let z = 0.5
        * ((2.0 * p1.z)
            + (-p0.z + p2.z) * t
            + (2.0 * p0.z - 5.0 * p1.z + 4.0 * p2.z - p3.z) * t2
            + (-p0.z + 3.0 * p1.z - 3.0 * p2.z + p3.z) * t3);
    Vec3::new(x, y, z)
}

/// First derivative (tangent vector) of the Catmull-Rom spline at `t in [0, 1]`.
fn catmull_rom_tangent(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
    let t2 = t * t;
    let x = 0.5
        * ((-p0.x + p2.x)
            + 2.0 * (2.0 * p0.x - 5.0 * p1.x + 4.0 * p2.x - p3.x) * t
            + 3.0 * (-p0.x + 3.0 * p1.x - 3.0 * p2.x + p3.x) * t2);
    let y = 0.5
        * ((-p0.y + p2.y)
            + 2.0 * (2.0 * p0.y - 5.0 * p1.y + 4.0 * p2.y - p3.y) * t
            + 3.0 * (-p0.y + 3.0 * p1.y - 3.0 * p2.y + p3.y) * t2);
    let z = 0.5
        * ((-p0.z + p2.z)
            + 2.0 * (2.0 * p0.z - 5.0 * p1.z + 4.0 * p2.z - p3.z) * t
            + 3.0 * (-p0.z + 3.0 * p1.z - 3.0 * p2.z + p3.z) * t2);
    Vec3::new(x, y, z)
}

fn scale_color(c: Color32, factor: f32, specular_add: f32) -> Color32 {
    let r = (f32::from(c.r()) * factor + specular_add).clamp(0.0, 255.0) as u8;
    let g = (f32::from(c.g()) * factor + specular_add).clamp(0.0, 255.0) as u8;
    let b = (f32::from(c.b()) * factor + specular_add).clamp(0.0, 255.0) as u8;
    Color32::from_rgb(r, g, b)
}

fn lerp_color(a: Color32, b: Color32, t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let mix = |x: u8, y: u8| -> u8 { (f32::from(x) + (f32::from(y) - f32::from(x)) * t) as u8 };
    Color32::from_rgb(mix(a.r(), b.r()), mix(a.g(), b.g()), mix(a.b(), b.b()))
}
