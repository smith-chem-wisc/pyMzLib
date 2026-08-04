# Changelog

Notable changes to pyMzLib. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow semantic versioning judged on the **Python** API — a change to the internal JSON
envelope is not a breaking change unless Python callers can see it.

## [Unreleased]

### Added
- Readers: `pymzlib.readers.identify()`, `read_results()`, and `formats()` — identify any of the 29
  result-file types mzLib recognises (returning the projections each supports), and read the three
  quantifiable formats into a uniform record view, each with per-format caveats about what its
  numbers do and do not mean.
- **Exhaustive readers coverage** — all 29 file types are now readable, up from 3.
  `pymzlib.readers.read_records()` reads *any* format mzLib recognises into that format's own
  fields (so TopPIC, Crux, MSFragger's peptide/protein tables and the FlashDeconv formats become
  reachable for the first time), naming every field it could not project rather than dropping it
  silently. Alongside it, the three remaining cross-format views: `read_features()`
  (`ms1_features`), `read_matches()` (`spectral_match`) and `read_spectra()` (scan headers, with
  peaks opt-in). Each reports per-format caveats — that `_ms1.feature` rows are expanded across
  charge states while Dinosaur's are not, that TopFD changed its retention-time unit mid-version so
  the unit is honestly `"unknown"`, that Casanovo's `is_decoy` is `None` because de novo sequencing
  has no decoys, and that nothing from a typed view is FDR-filtered.

### Fixed
- A failure inside mzLib's parallel spectra readers now reports its real cause
  (`MzLibException: Reading profile mode mzmls not supported`) instead of the wrapper
  (`AggregateException: One or more errors occurred.`), and a usage failure raised inside one still
  exits 2 rather than being reclassified as a fault.
- Non-finite numbers (an unbounded mzML scan window reports infinity) cross the wire as `null`
  instead of failing the whole read with a JSON serialization error.
- PRIDE Archive support: `pymzlib.pride.list_files()`, `download()`, `total_size_bytes()`,
  and the `PrideFile` type.
- PRIDE complete file listing: `pymzlib.pride.list_ftp_files()` and `approximate_total_size_bytes()`,
  with the `PrideFtpFile` type — the authoritative file list read by walking a project's FTP
  directory tree (mzLib #1121), for the projects where PRIDE's REST manifest is incomplete (for
  PXD000001 it omits five of 13 files, including the two largest). Sizes are PRIDE's rounded index
  sizes, so the
  total is an estimate over the whole project, the opposite trade-off from `total_size_bytes()`
  (exact-ish but over an incomplete manifest) (#12).
- Peptidoforms: `pymzlib.peptidoform.fragments()` — digest an annotated UniProt protein, apply its
  modifications, and fragment every peptide, disclosing the rules it applied (modification census,
  the silent isoform cap, the fixed-charge convention).
- FlashLFQ label-free quantification: `pymzlib.flashlfq.quantify()` — quantify a search's peptides
  across mzML runs with match-between-runs, returning typed `FlashLfqResults` / `Peptide` /
  `ProteinGroup` / `Peak`. Match-between-runs transfers are exposed via `result.peaks`.
- A self-contained .NET payload bundled in the wheel, so no .NET installation is required.
- Documentation site, including the reasoning behind each design decision.

### Fixed
- `Peptide.intensity()` now returns `0.0` (never `None`) when the wire value is `null`, matching the
  documented "0.0 when missing, never None" invariant (#7).
- Documentation corrections back-ported from the mzLibRust bake-off: the glycation-exclusion
  rationale, ETD's spurious y-ion over-count, `max_threads` as a correctness (not only performance)
  knob, PRIDE's decompressed-size / incomplete-manifest reporting, and the trypsin vs `trypsin|P`
  peptide-count figure.

[Unreleased]: https://github.com/smith-chem-wisc/pyMzLib/commits/main