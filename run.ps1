# HelixPrime Windows Build & Launch Helper
# Uses the bundled Rust 1.84 GNU toolchain and a local temp target dir (avoiding OneDrive lock quirks).

param(
    [switch]$Release,
    [switch]$TestOnly
)

$ErrorActionPreference = "Stop"
$toolchainBin = "$env:USERPROFILE\.rustup\toolchains\1.84.0-x86_64-pc-windows-gnu\lib\rustlib\x86_64-pc-windows-gnu\bin\self-contained"
$targetDir = Join-Path $env:LOCALAPPDATA "helix-prime-target"
$extraLibDir = Join-Path $targetDir "gnu-extra"

if (Test-Path $toolchainBin) {
    $cleanPath = ($env:PATH -split ';' | Where-Object { $_ -notlike '*MinGW*' }) -join ';'
    $env:PATH = "$toolchainBin;$cleanPath"

    # Provide libshlwapi.a alias from winapi-x86_64-pc-windows-gnu if present in Cargo cache.
    if (-not (Test-Path (Join-Path $extraLibDir "libshlwapi.a"))) {
        $src = Get-ChildItem "$env:USERPROFILE\.cargo\registry\src\*\winapi-x86_64-pc-windows-gnu*\lib\libwinapi_shlwapi.a" -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($src) {
            New-Item -ItemType Directory -Force $extraLibDir | Out-Null
            Copy-Item $src.FullName (Join-Path $extraLibDir "libshlwapi.a")
        }
    }

    if ($TestOnly) {
        cargo +1.84.0-x86_64-pc-windows-gnu test --workspace --target-dir $targetDir
        exit $LASTEXITCODE
    }

    $argsList = @("+1.84.0-x86_64-pc-windows-gnu", "rustc", "-p", "helix-gui", "--target-dir", $targetDir)
    if ($Release) { $argsList += "--release" }
    $argsList += @("--", "-L", $extraLibDir)

    & cargo @argsList
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    $profile = if ($Release) { "release" } else { "debug" }
    $exe = Join-Path $targetDir "$profile\helix-prime.exe"
    Write-Host "Launching $exe ..." -ForegroundColor Cyan
    & $exe
} else {
    if ($TestOnly) {
        cargo test --workspace --target-dir $targetDir
    } elseif ($Release) {
        cargo run --release -p helix-gui --target-dir $targetDir
    } else {
        cargo run -p helix-gui --target-dir $targetDir
    }
}
