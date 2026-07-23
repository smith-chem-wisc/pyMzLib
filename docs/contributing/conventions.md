# Conventions

The patterns every pyMzLib capability follows. A new tranche (a verb + its Python surface) should
match these rather than invent its own; they are why the PRIDE, peptidoform, and FlashLFQ tranches
feel like one library. Read this before [adding a capability](adding-a-capability.md).

## 1. One JSON envelope per call

Every bridge verb writes exactly one envelope to **stdout** and nothing else:

```json
{"ok": true,  "data": { … }}
{"ok": false, "error": {"type": "…", "message": "…"}}
```

Diagnostics — progress, warnings, anything a human might want — go to **stderr**. stdout must stay
parseable without stripping log lines. When a verb calls mzLib code that prints its own progress
(FlashLFQ does), redirect the console for the duration of that call so the envelope is never
corrupted; do not rely on a `Silent` flag alone.

## 2. Error classification lives in the bridge, not in Python

The bridge labels a failure by **kind**, so every consumer — Python today, a Rust binding tomorrow —
inherits the distinction:

- **Availability** (`ServiceUnavailable`): a timeout, a dropped socket, HTTP 408/429/5xx. "Try later."
- **Correctness** (the .NET exception type name): a 404, a malformed request, a parse failure. "This
  is a bug."
- **Usage** (`usage`, exit code 2): a malformed or missing argument. Raised before any work.

Python maps these to `ServiceUnavailableError`, `BridgeError`, and `UsageError`. Never launder a
correctness failure as an outage — a `NullReferenceException` reported as `ServiceUnavailable` is a
real bug that every test suite would then skip.

## 3. Large or variadic input on stdin, not argv

argv has a ~32 KB ceiling and a real experiment blows past it. A selection or a list of files goes
on **stdin**, one record per line; single scalar inputs (an accession, a PSM-file path) go on argv.
PRIDE download takes file names this way; FlashLFQ takes its spectra runs this way
(`path[⇥condition[⇥biorep[⇥techrep[⇥fraction]]]]`).

## 4. Injectable seams for testing

The parts worth testing are downstream of the network or the filesystem, so the boundary is a
replaceable factory: `Program.PrideClientFactory`, `Peptidoform.UniProtXmlSource`. A test swaps in a
stub over a fake `HttpMessageHandler` or a fixture file and exercises the whole verb offline. When a
new verb reaches out to a service or reads a file it cannot bundle, give it the same seam.

## 5. Follow mzLib / FlashLFQ names — do not invent

A value on the wire or in the Python API should mean exactly what it means in the mzLib source, the
MetaMorpheus output columns, and the papers. Use the snake_case of the mzLib name and stop there:
`match_between_runs` (not "mbr_enabled"), `ppm_tolerance` (not "tolerance"), `protein_groups`,
`detection_types`, `FlashLfqResults`, `ProteinGroup`. Renaming forces every reader to hold a
translation table in their head — and makes a support question ("what's `ppm_tolerance`?") answerable
straight from FlashLFQ's own docs.

!!! note "Known future rename — not yet actioned"
    The `peptidoform` module corresponds to mzLib's `PeptideWithSetModifications`. Aligning the name
    is on the table for a later pass (it is a public rename, so it is deliberately deferred, not done
    silently here). "Peptidoform" is also the HUPO-PSI term for the same concept, so this one is a
    judgement call to make explicitly when the time comes, not a slip to fix.

## 6. Disclose the traps; do not hide them

The recurring value of this library is refusing to let a caller be silently wrong. Where mzLib
applies an invisible rule, **surface it** rather than swallow it: the peptidoform tranche reports the
modification census (annotated vs applied), the silent isoform cap, and the fixed-charge convention;
the FlashLFQ tranche reports MBR peak counts and turns a median-polish NaN protein intensity into an
explicit `None` instead of a wrong number. A short answer and a truncated answer must not look alike.

## 7. Zero third-party Python runtime dependencies (D2)

`pymzlib` declares no runtime dependency, so it can never join a dependency conflict in someone's
environment. Do not add one without raising it explicitly — the .NET payload ships *inside* the
wheel instead. The bridge follows the same rule: take no NuGet package it does not need, because
every dependency is weight in the published wheel (and, as FlashLFQ showed, can drag native readers
worth pruning).

## 8. The bridge is a translation layer, not a second implementation

Compose mzLib's own methods; do not re-implement its logic in C#. FlashLFQ identifications come from
mzLib's `MakeIdentifications`, not a hand-rolled psmtsv parser; PRIDE filtering reuses
`WhereCategory`/`WhereExtension`. The bridge's job is to shape arguments in and results out — every
line of domain logic it contains is a line that can drift from mzLib.

## 9. The wire contract is language-neutral (D6)

Nothing in the envelope may assume a Python caller. All transport knowledge stays inside
`pkg/python/src/pymzlib/_bridge.py`, so a parallel Rust binding stays cheap.

## 10. Tests: offline fixtures by default, live canaries quarantined

The default test path is offline — a recorded fixture or a stub, so the suite passes whether or not
EBI or UniProt is up. Live calls go in a separate `ExternalService` category that **skips** (not
fails) on an outage, following mzLib's own convention. Coverage is gated; a new public surface needs
coverage, and the boundary (argument assembly, error mapping) is the half that matters.
