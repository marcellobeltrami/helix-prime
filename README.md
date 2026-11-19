# 🧬 HelixPrime

An open-source, native Rust alignment inspection and 3D protein visualization suite inspired by **Geneious Prime**, combining high-FPS multiple sequence alignment analysis with live **UniProtKB** annotation lookup and a **separate-window 3D Protein Viewer** connected to the **AlphaFold Protein Structure Database (EBI)** and **RCSB PDB**.

## Features

### Geneious Prime-Inspired Modern Desktop GUI (`helix-gui`)
- **Workspace & Sources Browser**: Collapsible folder tree with document cards displaying alphabet badges (`PRO` / `DNA` / `RNA`), sequence count, alignment width, and mean pairwise identity %.
- **High-Performance Alignment Canvas** (horizontally virtualized for smooth 60 FPS scrolling on large alignments):
  - **Overview Minimap Bar**: Full-length pairwise identity profile, miniature annotation markers, and click/drag viewport navigation box.
  - **Identity % Histogram**: Signature green (`100%`), amber (`30–99%`), and red (`<30%`) conservation bars per column.
  - **Consensus Track & Ruler**: Configurable consensus calling threshold & gap handling with labeled coordinate ticks.
  - **Packed Directional Annotation Tracks**: Non-overlapping interval-packed feature arrows (`Gene`, `CDS`, `Domain`, `Active Site`, `Binding Site`, `Variant`, `Helix`, `Beta Strand`, `Motif`, `Primer`) with `[3D]` badges for UniProt-linked features.
  - **CDS / Reading-Frame Translation Track**: Live gap-aware nucleotide-to-amino-acid translation (`+1`, `+2`, `+3`) for DNA/RNA alignments.
  - **Coloring & Highlighting Schemes**: *Geneious Standard*, *ClustalX*, *Polarity / Hydrophobicity*, *Zappo*, *Monochrome*, plus *Highlight Disagreements*, *Highlight Agreements*, and *Dots (`.`) for consensus matches*.
- **Right-Side Inspector Drawer**:
  - **Display**: Color schemes, highlighting modes, track visibility, zoom sliders, and consensus threshold.
  - **Annotations & UniProtKB Search**: Edit/create annotations, search UniProtKB live by gene/protein name, import curated UniProt domains/active sites directly onto the alignment, and launch AlphaFold 3D models in one click.
  - **Find Motif**: Gap-skipping IUPAC-aware motif search with hit list and Prev/Next jump controls.
  - **Statistics**: Alignment summary, reference GC %, and color-coded **Pairwise Identity Matrix**.
- **Light & Dark Modern Themes**: Toggle **Dark Mode** / **Light Mode** from the top toolbar or Display Inspector (synchronized across both the main alignment window and the separate 3D Protein Viewer window).
- **Phage Visualization (Circular Genome View)**: Toggle **Phage Circular** on the toolbar or Display Inspector to wrap any genome alignment into an interactive concentric circle featuring curved directional annotation arrows with outer callout labels and `[3D]` badges, coordinate ruler, multi-strain variant rings, pairwise identity ring, circular GC-content ring, drag-to-rotate, scroll-to-zoom, and bidirectional 3D protein sync. Includes a built-in **Bacteriophage Genome (Circular)** multi-strain alignment with UniProt/AlphaFold links (`P03713`, `P00720`, `P03034`, `P03700`).

### Separate-Window 3D Protein Viewer (UniProt + AlphaFold DB + RCSB PDB)
- Opens in a **true separate native OS window** (`egui::ViewportId`), allowing side-by-side alignment inspection and 3D structural analysis across monitors.
- **Direct Annotation-to-3D Connection**: Click or double-click any UniProt-linked annotation (protein domain OR DNA CDS!) to fetch its full-atom **AlphaFold DB v4** structure (with local disk caching) and **UniProtKB** metadata/PDB cross-references.
- **Bidirectional Alignment ↔ 3D Sync**:
  - Hovering or clicking a residue in the 3D window displays residue identity, chain, secondary structure, pLDDT confidence, and mapped MSA column + conservation %, and jumps the alignment viewport straight to that column/codon!
  - Clicking a column/codon in the alignment highlights the corresponding residue in the 3D protein window.
- **Pure-Rust 3D Renderer** (zero `unsafe` code):
  - **Styles**: *Cartoon Ribbon* (coiled alpha-helices, arrowed beta-strands, tubular loops), *Backbone Tube*, *C-Alpha Trace*, *Ball & Stick*, and *Spacefill (VDW)*.
  - **Coloring**: *AlphaFold pLDDT Confidence* (official 4-band DeepMind/EBI scale), *Secondary Structure*, *MSA Conservation* (projected directly from your multiple sequence alignment!), *Rainbow (N → C)*, and *Chain / CPK Element*.
  - **Local Structure Support**: Also loads local `.pdb` and `.cif` (mmCIF) files directly, and includes an instant **Offline Preview** fold.

### Core Bioinformatics Library (`helix-core`)
- Compact, cache-friendly Rust data structures with zero `.unwrap()` panics in library code.
- **File I/O**: FASTA (with automatic UniProt `sp|ACC|ENTRY OS=...` header parsing), ClustalW/Omega (`.aln`), GenBank (`.gb`/`.gbk` including `/db_xref="UniProtKB:..."`), GFF3 annotation import, PDB, and mmCIF parsers.
- **Progressive Multiple Sequence Aligner**: Built-in Needleman-Wunsch center-star aligner so unaligned FASTA imports can be aligned directly inside the GUI.

## Workspace Layout

```
genome-visualizer/
├── Cargo.toml
├── sample_data/                    # Sample FASTA, GenBank, and GFF3 files
└── crates/
    ├── helix-core/                 # Pure Rust bioinformatics, I/O, UniProt/AlphaFold clients & 3D model
    └── helix-gui/                  # Native eframe/egui modern desktop GUI + 3D renderer
```

## Building & Running

Requires Rust 1.84+ (`rustup` stable):

```bash
# Run the desktop application (release mode recommended for smoothest 60 FPS 3D interaction)
cargo run --release -p helix-gui

# Run all unit & integration tests
cargo test --workspace

# Run Clippy lints
cargo clippy --workspace --all-targets
```

> **Windows note**: If using the MSVC toolchain (`x86_64-pc-windows-msvc`), ensure the Visual Studio C++ Build Tools + Windows 10/11 SDK are installed. Alternatively, you can build with the self-contained GNU toolchain: `cargo +1.84.0-x86_64-pc-windows-gnu run --release -p helix-gui`.
