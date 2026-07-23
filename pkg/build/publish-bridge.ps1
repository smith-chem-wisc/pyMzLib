<#
.SYNOPSIS
    Publishes the mzLib bridge as a self-contained single-file executable and stages it inside
    the Python package, ready for a wheel build.

.DESCRIPTION
    This is the step that makes decision D2 true: after it runs, the Python package carries its
    own .NET runtime and a consumer needs nothing installed.

    One runtime identifier produces one platform's payload; a released wheel is per-platform,
    so this script is run once per target in CI (win-x64, linux-x64, osx-arm64, ...).

.PARAMETER Runtime
    The .NET runtime identifier to publish for. Defaults to this machine's.

.PARAMETER Configuration
    Build configuration. Release unless you are debugging the bridge itself.

.EXAMPLE
    .\publish-bridge.ps1
    Publishes for the current platform and stages it for a local wheel build.

.EXAMPLE
    .\publish-bridge.ps1 -Runtime linux-x64
    Cross-publishes the Linux payload (no Linux machine required).
#>
[CmdletBinding()]
param(
    [string]$Runtime,
    [ValidateSet('Release', 'Debug')]
    [string]$Configuration = 'Release'
)

$ErrorActionPreference = 'Stop'

# Forward slashes throughout: this script also runs on the Linux and macOS CI runners under pwsh,
# and .NET accepts them on Windows too.
$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$project = Join-Path $repoRoot 'pkg/bridge/MzLibBridge.csproj'
$stageRoot = Join-Path $repoRoot 'pkg/python/src/pymzlib/_dotnet'

if (-not $Runtime) {
    $arch = if ([Environment]::Is64BitOperatingSystem) { 'x64' } else { 'x86' }
    $Runtime = switch ($true) {
        $IsLinux   { "linux-$arch"; break }
        $IsMacOS   { "osx-$arch"; break }
        default    { "win-$arch" }
    }
}

if (-not (Test-Path (Join-Path $repoRoot 'code/mzLib/mzLib/mzLib.sln'))) {
    throw "The mzLib source is missing at code/mzLib. Recreate the worktree from the pin in code/PINNED.md (or, in CI, clone mzLib at that commit) before publishing."
}

$stageDir = Join-Path $stageRoot $Runtime
Write-Host "Publishing mzlib-bridge for $Runtime -> $stageDir"

if (Test-Path $stageDir) { Remove-Item -Recurse -Force $stageDir }
New-Item -ItemType Directory -Force -Path $stageDir | Out-Null

# The solution platform stays x64 to match the mzLib projects, but the codegen target must follow
# the runtime identifier: publishing osx-arm64 while PlatformTarget is x64 fails with NETSDK1032.
$platformTarget = if ($Runtime -match 'arm64$') { 'arm64' } else { 'x64' }

dotnet publish $project `
    -c $Configuration `
    -r $Runtime `
    --self-contained true `
    -p:Platform=x64 `
    -p:PlatformTarget=$platformTarget `
    -p:PublishSingleFile=true `
    -p:IncludeNativeLibrariesForSelfExtract=true `
    -p:EnableCompressionInSingleFile=true `
    -p:DebugType=none `
    -o $stageDir `
    --nologo
if ($LASTEXITCODE -ne 0) { throw "dotnet publish failed with exit code $LASTEXITCODE" }

# Only the executable belongs in the wheel; publish also emits pdb/config files we do not ship.
Get-ChildItem $stageDir -File |
    Where-Object { $_.Extension -in '.pdb', '.json', '.xml' } |
    Remove-Item -Force

$payload = Get-ChildItem $stageDir -File
$sizeMb = [math]::Round(($payload | Measure-Object Length -Sum).Sum / 1MB, 1)
Write-Host "Staged $($payload.Count) file(s), $sizeMb MB total:"
$payload | ForEach-Object { Write-Host "  $($_.Name)  $([math]::Round($_.Length / 1MB, 1)) MB" }

# The wheel's whole promise is that this runs on a machine with no .NET. A payload that cannot
# even report its own version here will certainly fail there.
$exeName = if ($Runtime.StartsWith('win')) { 'mzlib-bridge.exe' } else { 'mzlib-bridge' }
$exePath = Join-Path $stageDir $exeName
if (-not (Test-Path $exePath)) { throw "Expected $exeName in $stageDir but it is not there." }

$isNative = ($Runtime.StartsWith('win') -and $IsWindows -ne $false) -or
            ($Runtime.StartsWith('linux') -and $IsLinux) -or
            ($Runtime.StartsWith('osx') -and $IsMacOS)
if ($isNative) {
    $probe = & $exePath version | ConvertFrom-Json
    if (-not $probe.ok) { throw "Staged bridge failed its version probe." }
    Write-Host "Probe OK: bridge $($probe.data.bridge), protocol $($probe.data.protocol), runtime $($probe.data.runtime)"
} else {
    Write-Host "Cross-published for $Runtime; skipping the version probe (cannot run it here)."
}
