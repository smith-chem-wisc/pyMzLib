# Changelog

Notable changes to pyMzLib. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow semantic versioning judged on the **Python** API — a change to the internal JSON
envelope is not a breaking change unless Python callers can see it.

## [Unreleased]

### Added
- **DIA-NN and SDRF are readable**, following the pin to mzLib 1.0.585. `DiaNnReport` is the
  fourth format offering the `quantifiable` view, so DIA data can now feed
  `pymzlib.flashlfq.quantify()`; `Sdrf` (HUPO-PSI experimental design, `.sdrf.tsv`) reads through
  `read_records()` like any other format. This takes the supported count from 29 to 31 (mzLib
  #1120, #1138).

  DIA-NN retention times cross as `'minutes'` rather than `'unknown'`: DIA-NN writes minutes and
  mzLib converts nothing, which its own reader states. Note that mzLib dispatches this format on
  the file's **header**, not its name — a renamed DIA-NN report still reads, and a `report.tsv`
  that is not one still does not.
- **The raw bridge is now a release asset.** Every `v*` tag attaches `mzlib-bridge-<rid>.tar.gz`
  for all four platforms alongside the wheels, plus a `SHA256SUMS` covering both. A consumer with
  no reason to install a Python package — mzLibRust, a shell script, a container build — can
  unpack one and set `MZLIB_BRIDGE`. tar rather than zip so the executable bit survives extraction.
  `SHA256SUMS` also makes re-pinning mechanical for bindings that record these digests by hand
  (#31).
- Readers: `pymzlib.readers.identify()`, `read_results()`, and `formats()` — identify any of the 31
  result-file types mzLib recognises (returning the projections each supports), and read the four
  quantifiable formats into a uniform record view, each with per-format caveats about what its
  numbers do and do not mean.
- **Exhaustive readers coverage** — all 31 file types are now readable, up from 3.
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
- A download that dies part-way through now raises `ServiceUnavailableError` instead of a plain
  `BridgeError`. A request that fails outright carries a status code; one that fails *after* the
  response has begun does not — the server already said 200 — so it surfaced as a bare
  `IOException` (`Received an unexpected EOF or 0 bytes from the transport stream`) and escaped the
  availability classification entirely. A retry loop written around `ServiceUnavailableError`, as
  the PRIDE guide recommends, did not catch the one failure most worth retrying. Disk failures
  during a download are deliberately *not* reclassified: a full disk is still reported as itself
  (#30).
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
- Median-polish protein quantification: `pymzlib.flashlfq.median_polish()` — roll a
  `QuantifiedPeptides.tsv` up to protein intensities with FlashLFQ's own median-polish algorithm,
  without re-running peak-finding. Returns a list of `ProteinGroup`; takes an optional experimental
  design (condition/replicate grouping) and a `use_shared_peptides` toggle, and can write a
  `QuantifiedProteins.tsv`.
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