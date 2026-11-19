//! `helix-core`: Core bioinformatics data model, file parsers, progressive aligner,
//! UniProtKB/AlphaFold DB REST clients, and 3D protein structure mapping engine.

pub mod alignment;
pub mod annotation;
pub mod database;
pub mod demo;
pub mod error;
pub mod io;
pub mod sequence;
pub mod structure;

pub use alignment::{Alignment, ColumnStats, ConsensusConfig, MotifHit};
pub use annotation::{pack_annotation_lanes, Annotation, FeatureKind, UniProtLink};
pub use database::{
    fetch_alphafold_structure, fetch_rcsb_pdb, fetch_uniprot_entry, search_uniprot, AlphaFoldMeta,
    UniProtEntry, UniProtFeature, UniProtSearchHit,
};
pub use error::{HelixError, Result};
pub use sequence::{Alphabet, Sequence, Strand};
pub use structure::{
    parse_mmcif, parse_pdb, Atom, ProteinStructure, Residue3D, SecondaryStructure,
    StructureAlignmentMap, Vec3,
};
