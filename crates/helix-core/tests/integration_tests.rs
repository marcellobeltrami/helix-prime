//! End-to-end integration tests for `helix-core` covering file parsers, alignment
//! edge cases, progressive alignment, GFF3 import, mmCIF parsing, geometric secondary
//! structure inference, and alignment-to-3D coordinate mapping.

use helix_core::demo::{build_default_workspace, build_offline_preview_structure};
use helix_core::io::{
    clustal::parse_clustal, fasta::parse_fasta, genbank::parse_genbank,
    gff3::import_gff3_annotations, mafft::{parse_mafft, MafftFlavor},
};
use helix_core::{
    parse_mmcif, parse_pdb, Alignment, Alphabet, Sequence, SecondaryStructure,
    StructureAlignmentMap,
};

#[test]
fn default_workspace_invariants_hold() {
    let docs = build_default_workspace();
    assert_eq!(docs.len(), 12);
    for doc in &docs {
        let aln = &doc.alignment;
        assert!(aln.row_count() >= 1);
        assert!(aln.width() >= 20);
        assert_eq!(aln.columns().len(), aln.width());
        for seq in aln.sequences() {
            assert_eq!(seq.len(), aln.width());
        }
        for ann in aln.annotations() {
            assert!(ann.start < ann.end);
            assert!(ann.end <= aln.width());
        }
    }

    // Verify circular phage genome produces a valid circular GC profile.
    let phage = docs
        .iter()
        .find(|d| d.alignment.name.contains("Bacteriophage"))
        .expect("phage doc")
        .alignment
        .clone();
    let gc = phage
        .circular_gc_profile(72)
        .expect("phage DNA alignment produces circular GC profile");
    assert_eq!(gc.len(), 72);
}

#[test]
fn offline_preview_structure_maps_cleanly_to_hbb_protein_and_dna_cds() {
    let docs = build_default_workspace();
    let preview = build_offline_preview_structure();
    assert_eq!(preview.residues.len(), 147);
    assert!(preview.mean_plddt > 75.0);

    // Map against HBB protein alignment.
    let hbb_prot = &docs
        .iter()
        .find(|d| d.alignment.name.contains("Hemoglobin Beta Orthologs"))
        .expect("hbb prot")
        .alignment;
    let map_prot = StructureAlignmentMap::build(&preview, hbb_prot, hbb_prot.annotations().first());
    assert!(
        map_prot.mapped_count >= 140,
        "expected >= 140 mapped residues, got {}",
        map_prot.mapped_count
    );

    // Map against HBB DNA CDS alignment (444 bp -> 147 aa + stop codon).
    let hbb_dna = &docs
        .iter()
        .find(|d| d.alignment.name.contains("HBB Coding Sequence"))
        .expect("hbb dna")
        .alignment;
    let cds_ann = hbb_dna.annotations().first();
    let map_dna = StructureAlignmentMap::build(&preview, hbb_dna, cds_ann);
    assert!(
        map_dna.mapped_count >= 140,
        "expected >= 140 translated CDS codons mapped, got {}",
        map_dna.mapped_count
    );
    // Residue 6 (0-based idx 6, Glu7 / Sickle-cell HbS position) maps to DNA columns 18..21.
    assert_eq!(map_dna.residue_to_cols[6], Some((18, 21)));
    assert_eq!(map_dna.residue_for_column(19), Some(6));
}

#[test]
fn editing_extract_rc_and_translate_preserve_invariants() {
    let mut aln = Alignment::new(
        "EditTest",
        Alphabet::Dna,
        vec![Sequence::new("S1", b"ATGGTGCACCTG")],
    )
    .expect("aln");
    assert!(aln.set_residue(0, 0, b'g'));
    assert_eq!(aln.sequences()[0].residues[0], b'G');
    assert!(aln.insert_residue_at(0, 0, b'A'));
    assert_eq!(aln.width(), 13);
    assert!(aln.delete_column_at(1));
    assert_eq!(aln.width(), 12);

    let extracted = aln.extract_region("Sub", 0, 9).expect("extract");
    assert_eq!(extracted.width(), 9);

    let rc = extracted.reverse_complement_copy().expect("rc");
    assert_eq!(rc.width(), 9);

    let prot = extracted.translate_to_protein_alignment(0).expect("translate");
    assert_eq!(prot.alphabet, Alphabet::Protein);
    assert_eq!(prot.width(), 3);
}

#[test]
fn mmcif_parser_extracts_atoms_plddt_and_secondary_structure() {
    let cif = "\
data_AF-TEST
_struct.title 'Mini AlphaFold mmCIF Test'
loop_
_struct_conf.conf_type_id
_struct_conf.beg_auth_asym_id
_struct_conf.beg_auth_seq_id
_struct_conf.end_auth_seq_id
HELX_P A 1 3
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.type_symbol
ATOM CA MET A 1 0.0 0.0 0.0 91.5 C
ATOM CA VAL A 2 3.8 0.0 0.0 93.2 C
ATOM CA HIS A 3 7.6 1.0 0.0 88.4 C
ATOM CA LEU A 4 11.4 1.5 0.0 74.0 C
";
    let st = parse_mmcif(cif, "AF-TEST", Some("P99999".to_owned()), "Test", true)
        .expect("parse mmcif");
    assert_eq!(st.title, "Mini AlphaFold mmCIF Test");
    assert_eq!(st.residues.len(), 4);
    assert_eq!(st.one_letter_sequence(), b"MVHL");
    assert_eq!(st.residues[0].secondary_structure, SecondaryStructure::Helix);
    assert_eq!(st.residues[3].secondary_structure, SecondaryStructure::Coil);
    assert!((st.residues[1].plddt - 93.2).abs() < 0.01);
}

#[test]
fn gff3_import_and_error_paths_are_handled_cleanly() {
    let mut aln = Alignment::new(
        "Test",
        Alphabet::Protein,
        vec![Sequence::new("HBB_HUMAN", b"MVHLTPEEKSAVTALWGKVNV")],
    )
    .expect("valid alignment");

    let gff = "##gff-version 3\n\
               HBB_HUMAN\tUniProt\tActive_site\t5\t7\t.\t+\t.\tName=SiteA;Dbxref=UniProtKB:P68871\n";
    let count = import_gff3_annotations(&mut aln, gff).expect("import gff3");
    assert_eq!(count, 1);
    assert_eq!(aln.annotations()[0].name, "SiteA");
    assert_eq!(
        aln.annotations()[0]
            .uniprot
            .as_ref()
            .map(|u| u.accession.as_str()),
        Some("P68871")
    );

    // Malformed inputs must return typed errors rather than panicking.
    assert!(Alignment::new("Empty", Alphabet::Dna, vec![]).is_err());
    assert!(parse_fasta("Empty", "   \n\n").is_err());
    assert!(parse_fasta("NoHeader", "ACGTACGT").is_err());
    assert!(parse_clustal("Bad", "NOT A CLUSTAL FILE").is_err());
    assert!(parse_genbank("Empty", "LOCUS TEST\n//").is_err());
    assert!(parse_pdb("END\n", "EMPTY", None, "Test", false).is_err());
    assert!(import_gff3_annotations(&mut aln, "bad\trow\tonly3").is_err());
}

#[test]
fn sample_data_files_parse_cleanly() {
    let fasta_str = include_str!("../../../sample_data/hemoglobin_orthologs.fasta");
    let fasta_aln = parse_fasta("sample_fasta", fasta_str).expect("parse sample fasta");
    assert_eq!(fasta_aln.row_count(), 4);
    assert_eq!(fasta_aln.width(), 147);

    let gb_str = include_str!("../../../sample_data/hbb_annotated.gb");
    let gb_aln = parse_genbank("sample_gb", gb_str).expect("parse sample genbank");
    assert_eq!(gb_aln.width(), 444);
    assert!(gb_aln.annotations().len() >= 4);

    let gff_str = include_str!("../../../sample_data/hbb_features.gff3");
    let mut fasta_copy = fasta_aln;
    let added = import_gff3_annotations(&mut fasta_copy, gff_str).expect("import sample gff3");
    assert_eq!(added, 5);

    let mafft_str = include_str!("../../../sample_data/phage_capsid_mafft.mafft");
    let mafft_rep = parse_mafft("sample_mafft", mafft_str).expect("parse sample mafft");
    assert_eq!(mafft_rep.flavor, MafftFlavor::ClustalOut);
    assert_eq!(mafft_rep.alignment.row_count(), 4);
    assert_eq!(mafft_rep.alignment.width(), 120);
    assert_eq!(mafft_rep.adjusted_direction_count, 1);
}
