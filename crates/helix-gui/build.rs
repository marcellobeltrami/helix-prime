//! Build script that embeds BOTH:
//! 1. Multi-resolution application icon (`RT_ICON` Type 3 + `RT_GROUP_ICON` Type 14, ID 1)
//!    from `../../assets/helix-prime.ico`, so Windows Explorer, shortcuts, and Taskbar show
//!    the crisp HelixPrime DNA squircle icon.
//! 2. `helix-prime.manifest` (`RT_MANIFEST` Type 24, ID 1, `<requestedExecutionLevel level="asInvoker"/>`)
//!    so Windows never prompts for administrator / UAC elevation when launching `helix-prime.exe`.
//!
//! Generates a minimal, spec-compliant x86_64 COFF `.obj` file in pure Rust so it works
//! identically on both `x86_64-pc-windows-gnu` and `x86_64-pc-windows-msvc` with zero
//! external resource compiler (`windres` / `rc.exe`) dependencies.

use std::env;
use std::fs;
use std::path::PathBuf;

const RT_ICON: u32 = 3;
const RT_GROUP_ICON: u32 = 14;
const RT_MANIFEST: u32 = 24;
const LANG_EN_US: u32 = 0x0409;

struct ResourceLeaf {
    res_type: u32,
    res_id: u32,
    lang: u32,
    codepage: u32,
    data: Vec<u8>,
}

fn main() {
    println!("cargo:rerun-if-changed=helix-prime.manifest");
    println!("cargo:rerun-if-changed=../../assets/helix-prime.ico");

    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    if target_os != "windows" || target_arch != "x86_64" {
        return;
    }

    let manifest_bytes = fs::read("helix-prime.manifest")
        .expect("helix-prime.manifest must exist alongside build.rs");
    let ico_bytes = fs::read("../../assets/helix-prime.ico")
        .expect("assets/helix-prime.ico must exist for Windows icon embedding");

    let mut leaves = parse_ico_to_resources(&ico_bytes);
    leaves.push(ResourceLeaf {
        res_type: RT_MANIFEST,
        res_id: 1,
        lang: LANG_EN_US,
        codepage: 65001, // UTF-8
        data: manifest_bytes,
    });

    let coff_obj = build_x64_rsrc_coff(&leaves);

    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR is set by Cargo"));
    let obj_path = out_dir.join("helix_prime_resources.o");
    fs::write(&obj_path, coff_obj).expect("failed writing COFF resource object");

    println!("cargo:rustc-link-arg-bin=helix-prime={}", obj_path.display());
}

/// Parses a standard Windows `.ico` file into `RT_ICON` leaves (IDs `1..=N`) plus a single
/// `RT_GROUP_ICON` leaf (ID `1`) referencing them.
fn parse_ico_to_resources(ico: &[u8]) -> Vec<ResourceLeaf> {
    assert!(ico.len() >= 6, "ICO file too small");
    let reserved = u16::from_le_bytes([ico[0], ico[1]]);
    let ico_type = u16::from_le_bytes([ico[2], ico[3]]);
    let count = u16::from_le_bytes([ico[4], ico[5]]) as usize;
    assert_eq!(reserved, 0, "invalid ICO reserved header");
    assert_eq!(ico_type, 1, "expected ICO type 1");
    assert!(count > 0, "ICO must contain at least one icon image");
    assert!(ico.len() >= 6 + count * 16, "truncated ICO directory");

    let mut leaves = Vec::with_capacity(count + 1);
    let mut group_data = Vec::with_capacity(6 + count * 14);
    group_data.extend_from_slice(&ico[0..6]); // GRPICONDIR header (reserved=0, type=1, count)

    for i in 0..count {
        let off = 6 + i * 16;
        let bytes_in_res = u32::from_le_bytes([ico[off + 8], ico[off + 9], ico[off + 10], ico[off + 11]]) as usize;
        let img_offset = u32::from_le_bytes([ico[off + 12], ico[off + 13], ico[off + 14], ico[off + 15]]) as usize;
        let img_end = img_offset
            .checked_add(bytes_in_res)
            .expect("ICO entry overflow");
        assert!(img_end <= ico.len(), "ICO image slice out of bounds");

        let icon_id = (i + 1) as u16;
        // GRPICONDIRENTRY = first 12 bytes of ICONDIRENTRY + u16 resource ID.
        group_data.extend_from_slice(&ico[off..off + 12]);
        group_data.extend_from_slice(&icon_id.to_le_bytes());

        leaves.push(ResourceLeaf {
            res_type: RT_ICON,
            res_id: u32::from(icon_id),
            lang: LANG_EN_US,
            codepage: 0,
            data: ico[img_offset..img_end].to_vec(),
        });
    }

    leaves.push(ResourceLeaf {
        res_type: RT_GROUP_ICON,
        res_id: 1,
        lang: LANG_EN_US,
        codepage: 0,
        data: group_data,
    });

    leaves
}

/// Constructs a valid x86_64 PE/COFF object file containing a `.rsrc` section with an
/// arbitrary 3-level resource tree (`Type -> ID -> Language -> DataEntry`), properly sorted
/// by ascending numeric ID at every level as required by the Windows PE specification.
fn build_x64_rsrc_coff(leaves: &[ResourceLeaf]) -> Vec<u8> {
    // Group leaves into sorted `(res_type, Vec<(res_id, &ResourceLeaf)>)`.
    let mut types: Vec<u32> = leaves.iter().map(|l| l.res_type).collect();
    types.sort_unstable();
    types.dedup();

    let mut groups: Vec<(u32, Vec<&ResourceLeaf>)> = Vec::with_capacity(types.len());
    for t in types {
        let mut items: Vec<&ResourceLeaf> = leaves.iter().filter(|l| l.res_type == t).collect();
        items.sort_by_key(|l| (l.res_id, l.lang));
        groups.push((t, items));
    }

    // Compute exact byte offsets for all 4 sections of `.rsrc`:
    //   [Level 1 dir] -> [Level 2 dirs] -> [Level 3 dirs] -> [Data entries] -> [Payloads]
    let l1_size = 16 + groups.len() * 8;
    let l2_start = l1_size;
    let mut l2_offsets = Vec::with_capacity(groups.len());
    let mut cursor = l2_start;
    for (_, items) in &groups {
        l2_offsets.push(cursor);
        cursor += 16 + items.len() * 8;
    }

    let total_leaves: usize = groups.iter().map(|(_, items)| items.len()).sum();
    let l3_start = cursor;
    let data_entries_start = l3_start + total_leaves * (16 + 8);
    let payloads_start = data_entries_start + total_leaves * 16;

    let mut rsrc = Vec::with_capacity(payloads_start + 65_536);

    // Level 1 directory (Resource Types)
    push_rsrc_dir(&mut rsrc, groups.len() as u16);
    for (g_idx, (res_type, _)) in groups.iter().enumerate() {
        push_u32(&mut rsrc, *res_type);
        push_u32(&mut rsrc, 0x8000_0000 | (l2_offsets[g_idx] as u32));
    }
    debug_assert_eq!(rsrc.len(), l2_start);

    // Level 2 directories (Resource IDs per Type)
    let mut leaf_cursor = 0usize;
    for (_, items) in &groups {
        push_rsrc_dir(&mut rsrc, items.len() as u16);
        for item in items {
            let l3_off = l3_start + leaf_cursor * 24;
            push_u32(&mut rsrc, item.res_id);
            push_u32(&mut rsrc, 0x8000_0000 | (l3_off as u32));
            leaf_cursor += 1;
        }
    }
    debug_assert_eq!(rsrc.len(), l3_start);

    // Level 3 directories (Language per ID)
    let mut leaf_cursor = 0usize;
    for (_, items) in &groups {
        for item in items {
            let de_off = data_entries_start + leaf_cursor * 16;
            push_rsrc_dir(&mut rsrc, 1);
            push_u32(&mut rsrc, item.lang);
            push_u32(&mut rsrc, de_off as u32);
            leaf_cursor += 1;
        }
    }
    debug_assert_eq!(rsrc.len(), data_entries_start);

    // Precompute 4-byte-aligned payload offsets and emit IMAGE_RESOURCE_DATA_ENTRY table.
    let mut payload_cursor = payloads_start;
    let mut reloc_offsets: Vec<u32> = Vec::with_capacity(total_leaves);
    for (_, items) in &groups {
        for item in items {
            let entry_off = rsrc.len() as u32;
            reloc_offsets.push(entry_off); // OffsetToData field is at +0 of DATA_ENTRY
            push_u32(&mut rsrc, payload_cursor as u32); // Relocated via IMAGE_REL_AMD64_ADDR32NB
            push_u32(&mut rsrc, item.data.len() as u32);
            push_u32(&mut rsrc, item.codepage);
            push_u32(&mut rsrc, 0);

            payload_cursor += item.data.len();
            while payload_cursor % 4 != 0 {
                payload_cursor += 1;
            }
        }
    }
    debug_assert_eq!(rsrc.len(), payloads_start);

    // Emit payloads aligned to 4-byte boundaries.
    for (_, items) in &groups {
        for item in items {
            rsrc.extend_from_slice(&item.data);
            while rsrc.len() % 4 != 0 {
                rsrc.push(0);
            }
        }
    }
    debug_assert_eq!(rsrc.len(), payload_cursor);

    // --- Assemble full x86_64 COFF object ---
    let ptr_to_raw_data: u32 = 20 + 40; // COFF header (20) + 1 section header (40)
    let ptr_to_relocs: u32 = ptr_to_raw_data + rsrc.len() as u32;
    let num_relocs = reloc_offsets.len() as u16;
    let ptr_to_symtab: u32 = ptr_to_relocs + u32::from(num_relocs) * 10;

    let mut obj = Vec::with_capacity(ptr_to_symtab as usize + 18 + 4);

    // COFF File Header (20 bytes)
    push_u16(&mut obj, 0x8664); // IMAGE_FILE_MACHINE_AMD64
    push_u16(&mut obj, 1); // NumberOfSections
    push_u32(&mut obj, 0); // TimeDateStamp
    push_u32(&mut obj, ptr_to_symtab);
    push_u32(&mut obj, 1); // NumberOfSymbols
    push_u16(&mut obj, 0); // SizeOfOptionalHeader
    push_u16(&mut obj, 0); // Characteristics

    // Section Header `.rsrc` (40 bytes)
    obj.extend_from_slice(b".rsrc\0\0\0");
    push_u32(&mut obj, 0); // VirtualSize
    push_u32(&mut obj, 0); // VirtualAddress
    push_u32(&mut obj, rsrc.len() as u32);
    push_u32(&mut obj, ptr_to_raw_data);
    push_u32(&mut obj, ptr_to_relocs);
    push_u32(&mut obj, 0); // PointerToLinenumbers
    push_u16(&mut obj, num_relocs);
    push_u16(&mut obj, 0); // NumberOfLinenumbers
    push_u32(&mut obj, 0x4000_0040); // IMAGE_SCN_CNT_INITIALIZED_DATA | IMAGE_SCN_MEM_READ

    // Section Raw Data
    obj.extend_from_slice(&rsrc);

    // Relocation Table (10 bytes per IMAGE_RESOURCE_DATA_ENTRY)
    for reloc_va in reloc_offsets {
        push_u32(&mut obj, reloc_va);
        push_u32(&mut obj, 0); // SymbolTableIndex (symbol 0 = `.rsrc`)
        push_u16(&mut obj, 0x0003); // IMAGE_REL_AMD64_ADDR32NB
    }

    // Symbol Table (18 bytes): Symbol 0 = `.rsrc` static section symbol
    obj.extend_from_slice(b".rsrc\0\0\0");
    push_u32(&mut obj, 0); // Value
    push_u16(&mut obj, 1); // SectionNumber (1-based: section 1 = `.rsrc`)
    push_u16(&mut obj, 0); // Type
    obj.push(3); // StorageClass = IMAGE_SYM_CLASS_STATIC
    obj.push(0); // NumberOfAuxSymbols

    // Empty COFF String Table (4-byte length header = 4)
    push_u32(&mut obj, 4);

    obj
}

fn push_rsrc_dir(buf: &mut Vec<u8>, id_entries: u16) {
    push_u32(buf, 0); // Characteristics
    push_u32(buf, 0); // TimeDateStamp
    push_u16(buf, 0); // MajorVersion
    push_u16(buf, 0); // MinorVersion
    push_u16(buf, 0); // NumberOfNamedEntries
    push_u16(buf, id_entries);
}

fn push_u16(buf: &mut Vec<u8>, v: u16) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn push_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}
