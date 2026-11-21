//! Lightweight process & system RAM monitor for the bottom-left status bar
//! (matching Geneious Prime's `XXX / YYY MB Memory` indicator).

/// Snapshot of process resident memory and total physical RAM in megabytes.
#[derive(Debug, Clone, Copy, Default)]
pub struct MemorySnapshot {
    pub used_mb: u64,
    pub total_mb: u64,
}

/// Samples current process working-set memory and total system physical memory in MiB.
#[must_use]
pub fn sample_memory_mb() -> MemorySnapshot {
    #[cfg(windows)]
    {
        windows_mem::sample()
    }
    #[cfg(not(windows))]
    {
        MemorySnapshot {
            used_mb: 0,
            total_mb: 0,
        }
    }
}

#[cfg(windows)]
mod windows_mem {
    use super::MemorySnapshot;

    #[repr(C)]
    struct MemoryStatusEx {
        dw_length: u32,
        dw_memory_load: u32,
        ull_total_phys: u64,
        ull_avail_phys: u64,
        ull_total_page_file: u64,
        ull_avail_page_file: u64,
        ull_total_virtual: u64,
        ull_avail_virtual: u64,
        ull_avail_extended_virtual: u64,
    }

    #[repr(C)]
    struct ProcessMemoryCounters {
        cb: u32,
        page_fault_count: u32,
        peak_working_set_size: usize,
        working_set_size: usize,
        quota_peak_paged_pool_usage: usize,
        quota_paged_pool_usage: usize,
        quota_peak_non_paged_pool_usage: usize,
        quota_non_paged_pool_usage: usize,
        pagefile_usage: usize,
        peak_pagefile_usage: usize,
    }

    extern "system" {
        fn GlobalMemoryStatusEx(lp_buffer: *mut MemoryStatusEx) -> i32;
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(
            process: isize,
            ppsmem_counters: *mut ProcessMemoryCounters,
            cb: u32,
        ) -> i32;
    }

    pub(super) fn sample() -> MemorySnapshot {
        let mut status = MemoryStatusEx {
            dw_length: std::mem::size_of::<MemoryStatusEx>() as u32,
            dw_memory_load: 0,
            ull_total_phys: 0,
            ull_avail_phys: 0,
            ull_total_page_file: 0,
            ull_avail_page_file: 0,
            ull_total_virtual: 0,
            ull_avail_virtual: 0,
            ull_avail_extended_virtual: 0,
        };

        let mut pmc = ProcessMemoryCounters {
            cb: std::mem::size_of::<ProcessMemoryCounters>() as u32,
            page_fault_count: 0,
            peak_working_set_size: 0,
            working_set_size: 0,
            quota_peak_paged_pool_usage: 0,
            quota_paged_pool_usage: 0,
            quota_peak_non_paged_pool_usage: 0,
            quota_non_paged_pool_usage: 0,
            pagefile_usage: 0,
            peak_pagefile_usage: 0,
        };

        // SAFETY: Both `status` and `pmc` are properly aligned, stack-allocated, initialized
        // structs with exact byte size fields (`dw_length` and `cb`) required by Win32.
        let (used_bytes, total_bytes) = unsafe {
            let total = if GlobalMemoryStatusEx(&mut status) != 0 {
                status.ull_total_phys
            } else {
                0
            };
            let proc_handle = GetCurrentProcess();
            let used = if K32GetProcessMemoryInfo(proc_handle, &mut pmc, pmc.cb) != 0 {
                pmc.working_set_size as u64
            } else {
                0
            };
            (used, total)
        };

        const MIB: u64 = 1024 * 1024;
        MemorySnapshot {
            used_mb: used_bytes / MIB,
            total_mb: total_bytes / MIB,
        }
    }
}
