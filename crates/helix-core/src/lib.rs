//! `helix-core`: Core bioinformatics data model, file parsers, progressive aligner,
//! cloning simulator, primer designer, phylogenetic trees, local BLAST, NCBI & UniProt/AlphaFold
//! REST clients, and 3D protein structure mapping engine.

pub mod aligners;
pub mod alignment;
pub mod annotation;
pub mod assembly;
pub mod blast;
pub mod cloning;
pub mod database;
pub mod demo;
pub mod error;
pub mod features_db;
pub mod io;
pub mod ncbi;
pub mod phylo;
pub mod primers;
pub mod restriction;
pub mod sequence;
pub mod structure;

pub use aligners::{
    detect_minimap2_binary, parse_sam_into_alignment, run_alignment_algorithm, AlignerAlgorithm,
    AlignerConfig, AlignmentRunOutput, Minimap2Preset,
};
pub use alignment::{Alignment, ColumnStats, ConsensusConfig, MotifHit};
pub use annotation::{pack_annotation_lanes, Annotation, FeatureKind, UniProtLink};
pub use assembly::{assemble_contigs_denovo, map_reads_to_reference};
pub use blast::{run_blast, run_local_blast, BlastConfig, BlastHit, BlastProgram};
pub use cloning::{simulate_cloning, CloningMethod, CloningOutcome, CloningRequest};
pub use database::{
    fetch_alphafold_structure, fetch_protein_bundle, fetch_rcsb_pdb, fetch_uniprot_entry,
    search_uniprot, AlphaFoldMeta, FetchedProteinBundle, ProteinStructureSource, UniProtEntry,
    UniProtFeature, UniProtSearchHit,
};
pub use error::{HelixError, Result};
pub use features_db::{
    auto_annotate_from_reference_db, ReferenceFeature, PLASMID_REFERENCE_FEATURES,
};
pub use ncbi::{download_ncbi_record, search_ncbi, NcbiDatabase, NcbiSearchHit};
pub use phylo::{build_phylogenetic_tree, DistanceModel, PhyloNode, PhyloTree, TreeMethod};
pub use primers::{
    analyze_primer, design_primer_pairs, test_primer_binding, DesignedPrimerPair, PrimerBindingHit,
    PrimerDesignConfig, PrimerMetrics,
};
pub use restriction::{
    find_enzyme_by_name, find_restriction_sites, OverhangKind, RestrictionEnzyme,
    RestrictionSiteHit, COMMON_ENZYMES,
};
pub use sequence::{Alphabet, Sequence, Strand};
pub use structure::{
    parse_mmcif, parse_pdb, Atom, ProteinStructure, Residue3D, SecondaryStructure,
    StructureAlignmentMap, Vec3,
};
