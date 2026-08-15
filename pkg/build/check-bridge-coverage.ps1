<#
.SYNOPSIS
    Asserts that the bridge's own C# coverage meets a threshold.

.DESCRIPTION
    `dotnet test --collect:"XPlat Code Coverage"` measures every assembly the test host loads,
    which here means the whole of mzLib — tens of thousands of lines with their own test suite in
    their own repository. The resulting headline number is around 0.6%, which is not a low score
    so much as a meaningless one, and coverlet's Include filter is not honoured by the collector
    version in use.

    Rather than gate on a number that answers no question, this reads the Cobertura report and
    computes the rate over the bridge's own classes only. Codecov does the equivalent for its
    display, since it drops coverage for files that are not in the repository and mzLib's sources
    are not (they live in the gitignored code/mzLib checkout).

.PARAMETER ResultsPath
    Directory to search for coverage.cobertura.xml. Defaults to the test project's TestResults.

.PARAMETER Configuration
    The build configuration the coverage report came from. `dotnet test` builds Debug unless told
    otherwise, so this defaults to Debug and the freshness check is judged against the Debug test
    assembly. Pass -Configuration Release if the report came from a Release test run. Scoping the
    freshness check to one configuration is what stops a stray build of the OTHER configuration —
    newer on disk than a perfectly valid report — from being misread as staleness.

.PARAMETER Threshold
    Minimum line rate, as a percentage. Matches the Python side's fail_under.

.EXAMPLE
    .\check-bridge-coverage.ps1 -Threshold 90
#>
[CmdletBinding()]
param(
    [string]$ResultsPath,
    # Non-empty, so a caller that forwards an unset variable (-Configuration '' / $null) fails fast
    # rather than collapsing "bin/$Configuration" back to "bin/" and silently reinstating the
    # unscoped, globally-newest search this fix removes.
    [ValidateNotNullOrEmpty()]
    [string]$Configuration = 'Debug',
    # Lower than the Python side's 90 on purpose. The three verb handlers construct a
    # PrideArchiveClient internally and talk to EBI, so they cannot be unit-tested without
    # injecting an HttpClient. That refactor has since landed (see Program.PrideClientFactory), so
    # they are covered end-to-end by the Python live tests and not by C# units, and the honest
    # number is ~88%. Raise this to 90 once the handlers are injectable; do not lower it.
    [double]$Threshold = 85
)

$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
if (-not $ResultsPath) {
    $ResultsPath = Join-Path $repoRoot 'pkg/bridge.tests/TestResults'
}

$report = Get-ChildItem $ResultsPath -Recurse -Filter 'coverage.cobertura.xml' -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1

if (-not $report) {
    throw "No coverage.cobertura.xml under '$ResultsPath'. Run dotnet test with --collect:`"XPlat Code Coverage`" first."
}

# Refuse a report that predates the code it claims to measure.
#
# `dotnet test` without --collect leaves whatever report was there last time, and this script would
# happily gate on it: "newest file on disk" is not the same as "from this run". That failure is
# silent and convincing — it prints a healthy percentage with a brand-new class missing from the
# listing entirely, which is exactly how untested code gets through a gate that appears to be green.
#
# Filtered to the report's configuration, NOT the globally-newest DLL under bin: the folder holds
# every configuration, TFM and platform, so a Release build performed after a Debug test run leaves a
# Release assembly newer than a perfectly valid Debug report, and comparing against that globally
# newest DLL reports STALE for a report that describes the current code exactly. The freshness
# question is only meaningful against the assembly of the SAME configuration the report came from.
#
# Matched by a path SEGMENT equal to $Configuration rather than a bin/<Configuration> subtree, because
# a platform-qualified build nests it one level deeper: CI runs `-c Release -p:Platform=x64`, which
# emits bin/x64/Release/net10.0/, and a plain build emits bin/Release/net10.0/. Both carry a `Release`
# segment; `bin/Debug/...` does not, so the other configuration is still excluded.
$configSegment = "[\\/]$([regex]::Escape($Configuration))[\\/]"
$assembly = Get-ChildItem (Join-Path $repoRoot 'pkg/bridge.tests/bin') -Recurse -Filter 'MzLibBridge.Tests.dll' -ErrorAction SilentlyContinue |
    Where-Object { $_.FullName -match $configSegment } |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1

if (-not $assembly) {
    # No assembly means the freshness check cannot run at all. Silently skipping it would restore
    # exactly the hole this guard closes -- a green gate measured against an unknown-age report --
    # so say so rather than pass quietly.
    throw @"
Cannot verify coverage freshness: no '$Configuration' MzLibBridge.Tests.dll found under pkg/bridge.tests/bin.
Build the test project (configuration '$Configuration') first, then re-run with collection enabled.
If the report came from a different configuration, pass -Configuration <name> to match it.
"@
}

if ($report.LastWriteTime -lt $assembly.LastWriteTime) {
    throw @"
Coverage report is STALE: it predates the '$Configuration' test assembly, so it cannot describe the current code.
  report:   $($report.FullName)
            written $($report.LastWriteTime)
  assembly: $($assembly.FullName)
            written $($assembly.LastWriteTime)
Re-run with collection enabled:
  dotnet test pkg/bridge.tests/MzLibBridge.Tests.csproj --filter "TestCategory!=ExternalService" ``
      --collect:"XPlat Code Coverage" --settings pkg/bridge.tests/coverlet.runsettings
If a stray build of another configuration is newer than a valid report, pass -Configuration to match
the report's, or remove the stray output.
"@
}

# The freshness check trusts -Configuration to name the report's origin, but the report itself is
# chosen configuration-blind: `dotnet test` writes every run into the same TestResults regardless of
# configuration. So a report produced by a DIFFERENT configuration's run is being judged fresh against
# THIS configuration's (possibly much older) assembly, and the OK could be describing code the report
# never measured. Timestamps alone cannot tell a stray build apart from the report's true origin, so
# rather than pass silently, warn when ANOTHER configuration's assembly (one whose path carries no
# `$Configuration` segment) is newer than the report — the one observable sign the pairing may have
# drifted. Non-fatal: CI passes -Configuration matching its build so this never fires there, and a
# local developer's cross-configuration builds are their call; failing would reinstate the false STALE.
$otherNewer = Get-ChildItem (Join-Path $repoRoot 'pkg/bridge.tests/bin') -Recurse -Filter 'MzLibBridge.Tests.dll' -ErrorAction SilentlyContinue |
    Where-Object { $_.FullName -notmatch $configSegment -and $_.LastWriteTime -gt $report.LastWriteTime } |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1
if ($otherNewer) {
    Write-Warning @"
Coverage freshness was judged against the '$Configuration' assembly, but a newer test assembly exists
in another build output:
  $($otherNewer.FullName)
  written $($otherNewer.LastWriteTime)
If the report actually came from that build, this result does not describe it. Re-run this check with
-Configuration matching the report's, or re-run the tests with collection so the pairing is fresh.
"@
}

[xml]$coverage = Get-Content $report.FullName

# Only the bridge's own classes, and only ones a person wrote. The compiler emits a nested class
# per async method and per lambda closure (`Program/<PrideFilesAsync>d__5`, `Program/<>c__...`);
# counting those measures the C# compiler, not our tests, and double-counts the same source lines.
$ours = @($coverage.coverage.packages.package.classes.class |
    Where-Object { $_.name -like 'MzLibBridge.*' -and $_.name -notmatch '/<' })
if ($ours.Count -eq 0) {
    throw "The report contains no MzLibBridge classes. Did the test run actually execute the bridge?"
}

# Line rate weighted by statement count, not a mean of per-class rates: a one-line class must not
# count as much as a two-hundred-line one.
$total = 0; $covered = 0
foreach ($class in $ours) {
    foreach ($line in @($class.lines.line)) {
        if ($null -eq $line) { continue }
        $total++
        if ([int]$line.hits -gt 0) { $covered++ }
    }
}

$rate = if ($total -gt 0) { 100.0 * $covered / $total } else { 0 }

Write-Host ''
Write-Host 'Bridge coverage (MzLibBridge.* only)' -ForegroundColor Cyan
foreach ($class in ($ours | Sort-Object name)) {
    $classRate = 100.0 * [double]$class.'line-rate'
    Write-Host ('  {0,-45} {1,6:N1}%' -f $class.name, $classRate)
}
Write-Host ('  {0,-45} {1,6:N1}%  ({2}/{3} lines)' -f 'TOTAL', $rate, $covered, $total) -ForegroundColor Cyan
Write-Host ''

if ($rate -lt $Threshold) {
    Write-Host ("FAIL: bridge coverage {0:N1}% is below the {1:N0}% threshold." -f $rate, $Threshold) -ForegroundColor Red
    exit 1
}

Write-Host ("OK: bridge coverage {0:N1}% meets the {1:N0}% threshold." -f $rate, $Threshold) -ForegroundColor Green
