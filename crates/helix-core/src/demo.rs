//! Built-in curated biological alignments and an offline preview 3D protein structure
//! so users can immediately explore Geneious Prime-style inspection and UniProt/AlphaFold 3D viewing.

use crate::alignment::Alignment;
use crate::annotation::{Annotation, FeatureKind, UniProtLink};
use crate::sequence::{Alphabet, Sequence, Strand};
use crate::structure::{Atom, ProteinStructure, Residue3D, SecondaryStructure, Vec3};
use std::ops::{Add, Sub};

/// A workspace document folder entry mimicking Geneious Prime's Sources hierarchy.
#[derive(Debug, Clone)]
pub struct DemoDocument {
    pub folder: &'static str,
    pub alignment: Alignment,
}

/// Constructs the default set of biologically curated alignments loaded on startup.
#[must_use]
pub fn build_default_workspace() -> Vec<DemoDocument> {
    vec![
        DemoDocument {
            folder: "Protein Alignments",
            alignment: build_hbb_protein_alignment(),
        },
        DemoDocument {
            folder: "Protein Alignments",
            alignment: build_tp53_dbd_alignment(),
        },
        DemoDocument {
            folder: "Nucleotide CDS",
            alignment: build_hbb_dna_cds_alignment(),
        },
        DemoDocument {
            folder: "Phage Genomes",
            alignment: build_bacteriophage_circular_alignment(),
        },
    ]
}

/// Real Hemoglobin Subunit Beta (HBB, UniProt `P68871`) ortholog alignment (147 residues).
fn build_hbb_protein_alignment() -> Alignment {
    let seqs = vec![
        Sequence::new(
            "HBB_HUMAN (WT)",
            b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        )
        .with_description("Hemoglobin subunit beta (Wild-type canonical)")
        .with_organism("Homo sapiens")
        .with_uniprot("P68871"),
        Sequence::new(
            "HBB_HUMAN (HbS)",
            b"MVHLTPVEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        )
        .with_description("Sickle-cell anemia variant (p.Glu7Val / mature Glu6Val)")
        .with_organism("Homo sapiens")
        .with_uniprot("P68871"),
        Sequence::new(
            "HBB_PANTR (Chimp)",
            b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        )
        .with_description("Hemoglobin subunit beta")
        .with_organism("Pan troglodytes")
        .with_uniprot("P68873"),
        Sequence::new(
            "HBB_MOUSE",
            b"MVHLTDAEKAAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNAKVKAHGKKVITAFNDGLNHLDSLKGTFASLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGKDFTPAAQAAFQKVVAGVATALAHKYH",
        )
        .with_description("Hemoglobin subunit beta-1")
        .with_organism("Mus musculus")
        .with_uniprot("P02088"),
        Sequence::new(
            "HBB_BOVIN (Cow)",
            b"MLTAEEKAAVTAFWGKVKVDEVGGEALGRLLVVYPWTQRFFESFGDLSTADAVMNNPKVKAHGKKVLDSFSNGMKHLDDLKGTFAALSELHCDKLHVDPENFKLLGNVLVVVLARNFGKEFTPVLQADFQKVVAGVANALAHRYH",
        )
        .with_description("Hemoglobin subunit beta (N-terminal truncation)")
        .with_organism("Bos taurus")
        .with_uniprot("P02070"),
        Sequence::new(
            "HBB_CHICK (Chicken)",
            b"MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMVRAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLVRVVAHALARKYH",
        )
        .with_description("Hemoglobin subunit beta")
        .with_organism("Gallus gallus")
        .with_uniprot("P02112"),
    ];

    // Insert 2 N-terminal alignment gaps in Bovine HBB so homologous globin positions align cleanly.
    let mut aligned_seqs = seqs;
    let mut bov_padded = Vec::with_capacity(147);
    bov_padded.extend_from_slice(b"M--");
    bov_padded.extend_from_slice(&aligned_seqs[4].residues[1..]);
    aligned_seqs[4].residues = bov_padded;

    let mut aln = Alignment::new("Hemoglobin Beta Orthologs (HBB)", Alphabet::Protein, aligned_seqs)
        .unwrap_or_else(|_| unreachable!("hardcoded non-empty sequence list"))
        .with_description(
            "Vertebrate hemoglobin subunit beta orthologs + Sickle-cell HbS variant (UniProt P68871)",
        );

    let hbb_link = |start_res: usize| {
        UniProtLink::new("P68871")
            .with_start_residue(start_res)
            .with_metadata("HBB_HUMAN", "Hemoglobin subunit beta", "Homo sapiens")
    };

    aln.add_annotation(
        Annotation::new(0, "Globin Domain (P68871)", FeatureKind::Domain, 1, 147)
            .with_target_sequence(0)
            .with_notes("Full-length beta-globin fold (8 alpha-helices A–H). Click 'View 3D (AlphaFold)' to inspect in 3D.")
            .with_uniprot(hbb_link(2)),
    );
    aln.add_annotation(
        Annotation::new(0, "HbS Variant (E7V)", FeatureKind::Variant, 6, 7)
            .with_target_sequence(1)
            .with_notes("Sickle-cell disease missense substitution Glu7Val (classic mature chain Glu6Val) creating a hydrophobic polymerization patch.")
            .with_uniprot(hbb_link(7)),
    );
    aln.add_annotation(
        Annotation::new(0, "Distal His64 (O2 gate)", FeatureKind::BindingSite, 63, 64)
            .with_target_sequence(0)
            .with_notes("Distal histidine stabilizing bound dioxygen via hydrogen bonding.")
            .with_uniprot(hbb_link(64)),
    );
    aln.add_annotation(
        Annotation::new(0, "Proximal His93 (Heme Fe)", FeatureKind::ActiveSite, 92, 93)
            .with_target_sequence(0)
            .with_notes("Proximal histidine directly coordinating the heme prosthetic group iron atom.")
            .with_uniprot(hbb_link(93)),
    );
    aln.add_annotation(
        Annotation::new(0, "Heme Pocket Helix F", FeatureKind::Helix, 85, 98)
            .with_target_sequence(0)
            .with_notes("Conserved F-helix cradling the heme prosthetic group.")
            .with_uniprot(hbb_link(86)),
    );
    aln.add_annotation(
        Annotation::new(0, "2,3-BPG Binding", FeatureKind::BindingSite, 82, 83)
            .with_target_sequence(0)
            .with_notes("Lys83 contributes to allosteric 2,3-bisphosphoglycerate binding.")
            .with_uniprot(hbb_link(83)),
    );

    aln
}

/// Real TP53 Tumor Suppressor DNA-Binding Domain (residues 94..292, UniProt `P04637`).
fn build_tp53_dbd_alignment() -> Alignment {
    let human_dbd: &[u8] = b"SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRK";
    // Introduce clinically famous somatic cancer hotspot mutations (R175H, R248Q, R273H) in a tumor track.
    // Residue 94 corresponds to index 0, so residue R => index (R - 94):
    // R175 -> idx 81, R248 -> idx 154, R273 -> idx 179.
    let mut tumor_hotspots = human_dbd.to_vec();
    if let Some(slot) = tumor_hotspots.get_mut(81) {
        *slot = b'H'; // R175H structural zinc-loop mutant
    }
    if let Some(slot) = tumor_hotspots.get_mut(154) {
        *slot = b'Q'; // R248Q DNA minor-groove contact mutant
    }
    if let Some(slot) = tumor_hotspots.get_mut(179) {
        *slot = b'H'; // R273H DNA phosphate-backbone contact mutant
    }

    let mouse_dbd: &[u8] = b"SSSVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFCQLAKTCPVQLWVSATPPAGSRVRAMAIYKKSQHMTEVVRRCPHHERCSDGDGLAPPQHLIRVEGNLYPEYLEDRQTFRHSVVVPYEPPEAGSEYTTIHYKYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRDSFEVRVCACPGRDRRTEEENFRK";
    let zebrafish_dbd: &[u8] = b"SSSVPSQKTYPGNYGFRLGFLQSGTAKSVTCTYSPSLNKLFCQLAKTCPVQLWVSSTPPPGTRVRAMAIYKKSEHVTEVVRRCPHHERCSDGDGLAPPQHLIRVEGNQRANYREDNITLRHSVFVPYEAPQLGAEWTTVLLNYMCNSSCMGGMNRRPILTIITLETQEGQLLGRRSFEVRVCACPGRDRKTEESNFKK";

    let seqs = vec![
        Sequence::new("P53_HUMAN (WT)", human_dbd)
            .with_description("Cellular tumor antigen p53 DNA-binding core domain (aa 94-292)")
            .with_organism("Homo sapiens")
            .with_uniprot("P04637"),
        Sequence::new("P53_TUMOR (Hotspots)", &tumor_hotspots)
            .with_description("Composite somatic hotspot profile (R175H, R248Q, R273H)")
            .with_organism("Homo sapiens")
            .with_uniprot("P04637"),
        Sequence::new("P53_MOUSE", mouse_dbd)
            .with_description("Cellular tumor antigen p53 DNA-binding domain")
            .with_organism("Mus musculus")
            .with_uniprot("P02340"),
        Sequence::new("P53_DANRE (Zebrafish)", zebrafish_dbd)
            .with_description("Cellular tumor antigen p53 DNA-binding domain")
            .with_organism("Danio rerio")
            .with_uniprot("P79734"),
    ];

    let mut aln = Alignment::new("TP53 DNA-Binding Domain (P04637)", Alphabet::Protein, seqs)
        .unwrap_or_else(|_| unreachable!("hardcoded non-empty sequence list"))
        .with_description(
            "Evolutionary conservation & oncogenic hotspots across the p53 DNA-binding domain (aa 94-292)",
        );

    let p53_link = |start_res: usize| {
        UniProtLink::new("P04637")
            .with_start_residue(start_res)
            .with_metadata("P53_HUMAN", "Cellular tumor antigen p53", "Homo sapiens")
    };

    aln.add_annotation(
        Annotation::new(0, "p53 DNA-Binding Domain", FeatureKind::Domain, 0, human_dbd.len())
            .with_target_sequence(0)
            .with_notes("Immunoglobulin-like beta-sandwich scaffold + L2/L3 zinc-stabilized DNA-binding loops (residues 94–292 in UniProt P04637).")
            .with_uniprot(p53_link(94)),
    );
    aln.add_annotation(
        Annotation::new(0, "L2 Loop / Zn2+ (C176, H179)", FeatureKind::BindingSite, 70, 96)
            .with_target_sequence(0)
            .with_notes("Loop L2 coordinates a structural Zn2+ ion via Cys176 and His179.")
            .with_uniprot(p53_link(164)),
    );
    aln.add_annotation(
        Annotation::new(0, "Hotspot R175H", FeatureKind::Variant, 81, 82)
            .with_target_sequence(1)
            .with_notes("Most frequent structural cancer mutation in TP53; destabilizes Zn2+ coordination.")
            .with_uniprot(p53_link(175)),
    );
    aln.add_annotation(
        Annotation::new(0, "L3 DNA-Contact Loop (R248)", FeatureKind::ActiveSite, 143, 158)
            .with_target_sequence(0)
            .with_notes("Loop L3 inserts Arg248 directly into the DNA minor groove.")
            .with_uniprot(p53_link(237)),
    );
    aln.add_annotation(
        Annotation::new(0, "Hotspot R248Q", FeatureKind::Variant, 154, 155)
            .with_target_sequence(1)
            .with_notes("DNA minor-groove contact mutation abolishing sequence-specific transcription.")
            .with_uniprot(p53_link(248)),
    );
    aln.add_annotation(
        Annotation::new(0, "H2 DNA-Recognition Helix (R273)", FeatureKind::Helix, 178, 193)
            .with_target_sequence(0)
            .with_notes("C-terminal H2 helix making major-groove DNA contacts including Arg273 and Cys277.")
            .with_uniprot(p53_link(272)),
    );

    aln
}

/// Real Human HBB coding nucleotide alignment (exons 1–3 concatenated, 444 bp) comparing
/// Wild-Type, Sickle-Cell (`HbS`, c.20A>T), Hemoglobin C (`HbC`, c.19G>A), and Chimpanzee.
fn build_hbb_dna_cds_alignment() -> Alignment {
    let wt_cds: &[u8] = concat!(
        "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC",
        "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG",
        "AGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAG",
        "GTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGAC",
        "AACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT",
        "CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC",
        "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT",
        "GCCCTGGCCCACAAGTATCACTAA"
    )
    .as_bytes();

    // HbS: c.20A>T (GAG -> GTG, Glu -> Val)
    let mut hbs_cds = wt_cds.to_vec();
    hbs_cds[19] = b'T';

    // HbC: c.19G>A (GAG -> AAG, Glu -> Lys)
    let mut hbc_cds = wt_cds.to_vec();
    hbc_cds[18] = b'A';

    // Pan troglodytes (synonymous third-codon-position differences)
    let mut chimp_cds = wt_cds.to_vec();
    chimp_cds[8] = b'C'; // CAT -> CAC (His)
    chimp_cds[143] = b'C';
    chimp_cds[311] = b'T';

    let seqs = vec![
        Sequence::new("NM_000518 (HBB WT)", wt_cds)
            .with_description("Homo sapiens hemoglobin subunit beta (HBB) coding sequence")
            .with_organism("Homo sapiens")
            .with_uniprot("P68871"),
        Sequence::new("HBB_HbS (c.20A>T)", &hbs_cds)
            .with_description("Sickle-cell allele (GAG -> GTG, p.Glu7Val)")
            .with_organism("Homo sapiens")
            .with_uniprot("P68871"),
        Sequence::new("HBB_HbC (c.19G>A)", &hbc_cds)
            .with_description("Hemoglobin C allele (GAG -> AAG, p.Glu7Lys)")
            .with_organism("Homo sapiens")
            .with_uniprot("P68871"),
        Sequence::new("PANTR_HBB (Chimp)", &chimp_cds)
            .with_description("Pan troglodytes HBB coding sequence")
            .with_organism("Pan troglodytes")
            .with_uniprot("P68873"),
    ];

    let mut aln = Alignment::new("HBB Coding Sequence (DNA -> 3D)", Alphabet::Dna, seqs)
        .unwrap_or_else(|_| unreachable!("hardcoded non-empty sequence list"))
        .with_description(
            "Human HBB coding DNA alignment illustrating 6-frame translation and codon-to-3D AlphaFold mapping",
        );

    let hbb_link = UniProtLink::new("P68871")
        .with_start_residue(1)
        .with_metadata("HBB_HUMAN", "Hemoglobin subunit beta", "Homo sapiens");

    aln.add_annotation(
        Annotation::new(0, "HBB CDS (UniProt P68871)", FeatureKind::Cds, 0, wt_cds.len())
            .with_strand(Strand::Forward)
            .with_target_sequence(0)
            .with_notes("Complete coding sequence translated in Frame +1 -> maps every codon directly onto AlphaFold 3D structure AF-P68871-F1.")
            .with_uniprot(hbb_link.clone()),
    );
    aln.add_annotation(
        Annotation::new(0, "Sickle/HbC Codon 7 (GAG)", FeatureKind::Variant, 18, 21)
            .with_strand(Strand::Forward)
            .with_target_sequence(1)
            .with_notes("Clinically critical codon: GAG (Glu) -> GTG (Val, HbS) or AAG (Lys, HbC).")
            .with_uniprot(hbb_link.clone().with_start_residue(7)),
    );
    aln.add_annotation(
        Annotation::new(0, "Proximal His93 Codon", FeatureKind::ActiveSite, 276, 279)
            .with_strand(Strand::Forward)
            .with_target_sequence(0)
            .with_notes("CAC codon encoding the heme-coordinating proximal histidine (His93).")
            .with_uniprot(hbb_link.with_start_residue(93)),
    );

    aln
}

/// Curated modular bacteriophage circular genome alignment comparing wild-type and host-range /
/// lysis-evolved phage isolates, with real UniProt/AlphaFold links (T4 Lysozyme `P00720`,
/// Lambda Capsid `P03713`, Lambda cI Repressor `P03034`, and Lambda Integrase `P03700`).
fn build_bacteriophage_circular_alignment() -> Alignment {
    // Modular 600 bp bacteriophage genome cassette with realistic GC% variation:
    // - AT-rich ori/cos packaging region (0..60)
    // - GC-rich Head/Capsid structural cassette (60..210)
    // - Tail fiber & receptor-binding module (210..330)
    // - Lysis cassette: Holin + Endolysin/Lysozyme (330..450)
    // - Reverse-strand lysogeny switch: cI Repressor & Integrase (450..600)
    let wt_genome: &[u8] = concat!(
        // 0..60: cos/ori AT-rich packaging & replication origin
        "AATTTAAATGGGCGGCGACCTTAATTAATTTATATGCTTAAGCGATCGATAATTTATAAA",
        // 60..210: Major Capsid Protein E module (GC-rich structural core)
        "ATGCGTCAGCTCGGCCCGTTCGTGGCCGAGCTGGACGGCCTGCCGTCCAAGCTGACCGAG",
        "GGCACCAAGCTGATCAAGGCCGACCTGGCCCAGTTCGTGCGCGAGAACGGCTACACCGTG",
        "CTGGGCGTGAACGACCTGACCGGCAAGCTG",
        // 210..330: Tail Fiber / Receptor-Binding Domain (hypervariable host-range locus)
        "ATGACCAACGTCAGCGGCACCTGGCAGTACAACGGCTCCAACACCATCACCGGCAACGTC",
        "AGCGGCTGGCAGACCATCTACAACGCCAACACCAGCGGCATCGTGAACAAGACCTACTGA",
        // 330..450: Lysis Cassette (Holin + Endolysin / Lysozyme T4-like)
        "ATGAATATATTTGAAATGTTACGTATAGATGAAGGTCTTAGACTTAAAATCTATAAAGAC",
        "ACAGAAGGCTATTACACTATTGGCATCGGTCATTTGCTTACAAAAAGTCCATCACTTAAT",
        // 450..600: Lysogeny Control (cI Repressor & Integrase on reverse strand)
        "TTAGCCAAACGTCTCTTCAGGCCACTGACTAGCGATAACTTTCCCCACAACGGAACAACT",
        "CTCATTGCATGGGATCATTGGGTACTGTGGGTTTAGTGGTTGTAAAAACACCTGACCGCT",
        "ATCCCTGATCAGTTTCTTGAAGGTAAACTC"
    )
    .as_bytes();

    // Isolate B: Host-range adapted tail-fiber variant (clustered mutations in 225..315).
    let mut host_range = wt_genome.to_vec();
    for &(idx, base) in &[
        (228, b'T'),
        (235, b'A'),
        (244, b'C'),
        (258, b'G'),
        (267, b'T'),
        (279, b'A'),
        (292, b'C'),
        (304, b'G'),
        (312, b'T'),
        (356, b'G'),
    ] {
        host_range[idx] = base;
    }

    // Isolate C: Thermostable capsid & accelerated-lysis variant.
    let mut fast_lysis = wt_genome.to_vec();
    for &(idx, base) in &[
        (78, b'A'),
        (114, b'T'),
        (165, b'C'),
        (244, b'T'),
        (267, b'A'),
        (345, b'C'),
        (368, b'A'),
        (394, b'G'),
        (418, b'T'),
        (510, b'C'),
    ] {
        fast_lysis[idx] = base;
    }

    // Isolate D: Clear-plaque virulent mutant (mutations in cI repressor + small indel in tail).
    let mut virulent = wt_genome.to_vec();
    for &(idx, base) in &[
        (260, b'-'),
        (261, b'-'),
        (262, b'-'),
        (472, b'A'),
        (486, b'T'),
        (498, b'G'),
        (515, b'A'),
        (532, b'C'),
        (548, b'T'),
    ] {
        virulent[idx] = base;
    }

    let seqs = vec![
        Sequence::new("Phage_WT (Circular)", wt_genome)
            .with_description("Reference bacteriophage circular genome assembly")
            .with_organism("Escherichia phage"),
        Sequence::new("Phage_HR2 (Tail Var)", &host_range)
            .with_description("Host-range evolved isolate with receptor-binding mutations")
            .with_organism("Escherichia phage"),
        Sequence::new("Phage_LY4 (Fast Lysis)", &fast_lysis)
            .with_description("Accelerated endolysin & thermostable capsid isolate")
            .with_organism("Escherichia phage"),
        Sequence::new("Phage_VIR (cI- Mut)", &virulent)
            .with_description("Obligately lytic clear-plaque mutant in cI immunity module")
            .with_organism("Escherichia phage"),
    ];

    let mut aln = Alignment::new("Bacteriophage Genome (Circular)", Alphabet::Dna, seqs)
        .unwrap_or_else(|_| unreachable!("hardcoded non-empty phage sequences"))
        .with_description(
            "Circular bacteriophage genome alignment (Capsid, Tail Fiber, Endolysin, cI Repressor, Integrase) — toggle Phage Circular View to inspect as a ring!",
        );

    aln.add_annotation(
        Annotation::new(0, "cos / ori (Origin)", FeatureKind::Motif, 0, 58)
            .with_strand(Strand::Forward)
            .with_notes("AT-rich cohesive packaging site (cos) and bidirectional replication origin."),
    );
    aln.add_annotation(
        Annotation::new(0, "Major Capsid E", FeatureKind::Cds, 60, 208)
            .with_strand(Strand::Forward)
            .with_notes("Icosahedral head major capsid subunit (Lambda gpE). Double-click to inspect AlphaFold 3D structure.")
            .with_uniprot(
                UniProtLink::new("P03713")
                    .with_start_residue(1)
                    .with_metadata("CAPSD_LAMBD", "Major capsid protein", "Escherichia phage lambda"),
            ),
    );
    aln.add_annotation(
        Annotation::new(0, "Tail Fiber (RBP)", FeatureKind::Gene, 210, 328)
            .with_strand(Strand::Forward)
            .with_notes("Hypervariable receptor-binding tail fiber determining bacterial host range."),
    );
    aln.add_annotation(
        Annotation::new(0, "Endolysin (T4 Lysozyme)", FeatureKind::Cds, 330, 448)
            .with_strand(Strand::Forward)
            .with_notes("Peptidoglycan-degrading endolysin / T4 lysozyme (UniProt P00720). Double-click for AlphaFold 3D model.")
            .with_uniprot(
                UniProtLink::new("P00720")
                    .with_start_residue(1)
                    .with_metadata("LYS_BPT4", "Endolysin (Lysozyme)", "Escherichia phage T4"),
            ),
    );
    aln.add_annotation(
        Annotation::new(0, "cI Repressor", FeatureKind::Cds, 452, 535)
            .with_strand(Strand::Reverse)
            .with_notes("Lysogenic switch repressor (UniProt P03034); mutated in obligately lytic clear-plaque isolates.")
            .with_uniprot(
                UniProtLink::new("P03034")
                    .with_start_residue(1)
                    .with_metadata("RPC1_LAMBD", "Repressor protein cI", "Escherichia phage lambda"),
            ),
    );
    aln.add_annotation(
        Annotation::new(0, "Integrase (int)", FeatureKind::Cds, 538, 598)
            .with_strand(Strand::Reverse)
            .with_notes("Tyrosine site-specific recombinase catalyzing prophage integration (UniProt P03700).")
            .with_uniprot(
                UniProtLink::new("P03700")
                    .with_start_residue(1)
                    .with_metadata("INT_LAMBD", "Integrase", "Escherichia phage lambda"),
            ),
    );

    aln
}

/// Builds a realistic procedural 3D helical-bundle + beta-hairpin protein fold matching
/// the HBB canonical sequence (`P68871`) for immediate offline 3D inspection before fetching
/// the live AlphaFold DB model over the network.
#[must_use]
pub fn build_offline_preview_structure() -> ProteinStructure {
    let seq: &[u8] = b"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
    let mut atoms: Vec<Atom> = Vec::with_capacity(seq.len() * 4);
    let mut residues: Vec<Residue3D> = Vec::with_capacity(seq.len());

    // Construct an 8-segment globin-inspired compact 3D fold (helices A-H connected by turns)
    // with realistic ~3.8 Å C-alpha spacing and backbone N/CA/C/O atoms.
    struct Segment {
        len: usize,
        ss: SecondaryStructure,
        dir: Vec3,
        plddt_base: f32,
    }

    let segments = [
        Segment { len: 4,  ss: SecondaryStructure::Coil,  dir: Vec3::new(1.0, 0.2, 0.0),   plddt_base: 64.0 },
        Segment { len: 16, ss: SecondaryStructure::Helix, dir: Vec3::new(0.85, 0.35, 0.35), plddt_base: 93.0 },
        Segment { len: 4,  ss: SecondaryStructure::Coil,  dir: Vec3::new(0.0, 0.8, 0.5),   plddt_base: 78.0 },
        Segment { len: 18, ss: SecondaryStructure::Helix, dir: Vec3::new(-0.75, 0.25, 0.6), plddt_base: 95.0 },
        Segment { len: 5,  ss: SecondaryStructure::Sheet, dir: Vec3::new(-0.2, -0.9, 0.2), plddt_base: 86.0 },
        Segment { len: 5,  ss: SecondaryStructure::Coil,  dir: Vec3::new(0.5, -0.4, -0.7), plddt_base: 72.0 },
        Segment { len: 17, ss: SecondaryStructure::Helix, dir: Vec3::new(0.65, -0.55, -0.5), plddt_base: 94.0 },
        Segment { len: 5,  ss: SecondaryStructure::Coil,  dir: Vec3::new(0.2, 0.3, -0.9),  plddt_base: 76.0 },
        Segment { len: 20, ss: SecondaryStructure::Helix, dir: Vec3::new(-0.4, 0.8, -0.45), plddt_base: 96.0 },
        Segment { len: 5,  ss: SecondaryStructure::Sheet, dir: Vec3::new(0.8, 0.1, -0.5),  plddt_base: 88.0 },
        Segment { len: 4,  ss: SecondaryStructure::Coil,  dir: Vec3::new(0.5, -0.7, 0.5),  plddt_base: 74.0 },
        Segment { len: 20, ss: SecondaryStructure::Helix, dir: Vec3::new(-0.7, -0.6, 0.35), plddt_base: 95.0 },
        Segment { len: 5,  ss: SecondaryStructure::Coil,  dir: Vec3::new(0.4, 0.2, 0.9),   plddt_base: 81.0 },
        Segment { len: 19, ss: SecondaryStructure::Helix, dir: Vec3::new(0.8, 0.45, -0.35), plddt_base: 92.0 },
    ];

    let mut cursor = Vec3::ZERO;
    let mut seq_idx = 0usize;

    for seg in &segments {
        let axis = seg.dir.normalize();
        let ref_up = if axis.y.abs() < 0.9 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let u = axis.cross(ref_up).normalize();
        let v = axis.cross(u).normalize();

        for local in 0..seg.len {
            if seq_idx >= seq.len() {
                break;
            }
            let aa = seq[seq_idx];
            let ca = match seg.ss {
                SecondaryStructure::Helix => {
                    // Canonical alpha-helix: 3.6 residues/turn (100 deg/res), 1.5 Å rise, 2.3 Å radius.
                    let angle = (local as f32) * 1.745_329;
                    let radial = u.scale(angle.cos() * 2.3).add(v.scale(angle.sin() * 2.3));
                    cursor = cursor.add(axis.scale(1.52));
                    cursor.add(radial)
                }
                SecondaryStructure::Sheet => {
                    // Pleated beta-strand: 3.3 Å step with alternating +-0.8 Å pleat.
                    let pleat = if local % 2 == 0 { 0.85 } else { -0.85 };
                    cursor = cursor.add(axis.scale(3.3));
                    cursor.add(u.scale(pleat))
                }
                SecondaryStructure::Coil => {
                    let wave = ((local as f32) * 1.1).sin() * 1.1;
                    cursor = cursor.add(axis.scale(3.1)).add(v.scale(wave * 0.4));
                    cursor
                }
            };

            let plddt = (seg.plddt_base + ((seq_idx as f32 * 0.7).sin() * 3.5)).clamp(45.0, 99.0);
            let atom_start = atoms.len() as u32;
            let res_idx = residues.len() as u32;

            // Add N, CA, C, O heavy backbone atoms around CA.
            let n_pos = ca.sub(axis.scale(0.9)).add(u.scale(0.55));
            let c_pos = ca.add(axis.scale(0.95)).add(v.scale(0.5));
            let o_pos = c_pos.add(u.scale(-1.1));

            let three = crate::sequence::one_to_three_aa(aa);
            let mut res_name = [b' '; 3];
            for (slot, b) in res_name.iter_mut().zip(three.bytes()) {
                *slot = b;
            }

            atoms.push(Atom {
                pos: n_pos,
                b_factor: plddt,
                residue_idx: res_idx,
                name: *b"N   ",
                element: *b"N ",
                is_hetatm: false,
                is_ca: false,
            });
            let ca_atom_idx = atoms.len() as u32;
            atoms.push(Atom {
                pos: ca,
                b_factor: plddt,
                residue_idx: res_idx,
                name: *b"CA  ",
                element: *b"C ",
                is_hetatm: false,
                is_ca: true,
            });
            atoms.push(Atom {
                pos: c_pos,
                b_factor: plddt,
                residue_idx: res_idx,
                name: *b"C   ",
                element: *b"C ",
                is_hetatm: false,
                is_ca: false,
            });
            atoms.push(Atom {
                pos: o_pos,
                b_factor: plddt,
                residue_idx: res_idx,
                name: *b"O   ",
                element: *b"O ",
                is_hetatm: false,
                is_ca: false,
            });

            residues.push(Residue3D {
                atom_start,
                ca_atom: Some(ca_atom_idx),
                seq_id: (seq_idx + 1) as i32,
                plddt,
                chain_id: 'A',
                atom_count: 4,
                name: res_name,
                one_letter: aa,
                secondary_structure: seg.ss,
            });

            seq_idx += 1;
        }
    }

    // Recenter around origin for smooth orbital camera.
    let mut sum = Vec3::ZERO;
    for a in &atoms {
        sum = sum.add(a.pos);
    }
    let center = sum.scale(1.0 / atoms.len().max(1) as f32);
    let mut max_r_sq = 25.0f32;
    for a in &atoms {
        let d_sq = a.pos.sub(center).length_sq();
        if d_sq > max_r_sq {
            max_r_sq = d_sq;
        }
    }
    let mean_plddt = residues.iter().map(|r| r.plddt).sum::<f32>() / residues.len().max(1) as f32;

    ProteinStructure {
        title: "Hemoglobin Subunit Beta (Offline Preview — Click 'Fetch AlphaFold DB' for live model)".to_owned(),
        accession: Some("P68871".to_owned()),
        entry_id: "PREVIEW-P68871".to_owned(),
        organism: Some("Homo sapiens".to_owned()),
        source_label: "Offline Schematic Preview".to_owned(),
        is_alphafold: true,
        atoms: atoms.into_boxed_slice(),
        residues: residues.into_boxed_slice(),
        center,
        radius: max_r_sq.sqrt(),
        mean_plddt,
    }
}
