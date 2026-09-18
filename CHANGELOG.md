# Changelog

Notable changes to pyMzLib. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow semantic versioning judged on the **Python** API — a change to the internal JSON
envelope is not a breaking change unless Python callers can see it.

## [Unreleased]

### Fixed
- **Reading a Thermo `.raw` file works.** In 0.1.0 every `.raw` read failed with
  `BridgeError: Method invocation failed on Method[ThreadedFileFactory]`, on every platform. Thermo's
  RawFileReader looks for its own assemblies on disk and cannot run from inside the single-file
  bridge executable, which is where they were packed. They now ship as two loose DLLs beside it. No
  test caught this because the bridge's C# suite runs mzLib from loose assemblies; a new test reads a
  small `.raw` through the packaged bridge on every platform CI tests the wheel on.

  mzLibR and mzLibRust run the same bridge, so they get the fix too. The `mzlib-bridge-<rid>.tar.gz`
  asset now carries the two DLLs; extracting the whole archive, as the release docs already say to,
  picks them up.

## [0.1.0] - 2026-09-18

The first release on PyPI: `pip install mzlib`.


### Changed
- **The PyPI package is now `pip install mzlib`, not `pip install pymzlib`.** The import name is
  unchanged — it is still `import pymzlib`, and the package directory is still `src/pymzlib`. A
  distribution name and an import name are independent in Python (`pip install scikit-learn` gives
  you `import sklearn`), so no user code changes.

  `pymzlib` was verified unclaimed on PyPI on 2026-07-23 and claimed by an unrelated uploader on
  2026-07-29, before this project had published anything. A name nobody has uploaded to is not
  reserved. `mzlib` was claimed on 2026-08-16 and is ours.

  The wheels attached to a release are renamed to match (`mzlib-<version>-py3-none-<platform>.whl`).
  This does not reach mzLibR or mzLibRust: both take the `mzlib-bridge-<rid>.tar.gz` asset rather
  than unzipping a wheel, and the import package they look for inside a wheel is still `pymzlib/`.

- **The wheels are ~60 MB, down from up to 166 MB, and now fit on PyPI.** Nearly all of the old size
  was libtorch, a machine-learning library mzLib pulls in for a retention-time predictor that no
  pyMzLib function calls. It is no longer shipped; bridge output is byte-identical without it.
  Linux 165.9 → 56.4 MiB, Windows 129.4 → 59.3, macOS Apple Silicon 101.8 → 59.9. The .NET runtime
  is still inside every wheel, so there is still nothing to install. The same saving reaches
  mzLibRust and mzLibR through the `mzlib-bridge-<rid>.tar.gz` assets.

### Added
- **Find PRIDE projects by keyword**: `pymzlib.pride.search()`. Every other function in that module
  takes an accession you already have; this is the one that produces them, so you can go from a
  subject to a dataset without leaving Python. Paging is handled and no accession is repeated
  (mzLib #1187).

  **A hit is not a project's metadata.** PRIDE serves search from a separate Elasticsearch
  projection in which every controlled-vocabulary field is flattened to a display string —
  instruments come back as `["Q Exactive"]` rather than terms with accessions, contacts as display
  names, publications as one pre-formatted citation. That is PRIDE's wire, not a simplification
  chosen here, so `PrideProjectSearchResult` is its own type with string collections. What it adds
  over the metadata endpoint is `highlights`: which fields matched, and with what.

  Two honesty notes carried through rather than smoothed over. Dates are `datetime.date`, not
  `datetime` as on `PrideFile`, because this endpoint sends a bare calendar date with no time or
  offset. And a zero or an empty list means *not reported* — PRIDE omits nothing as null, and
  several fields are sparse — so `download_count == 0` does not mean nobody downloaded it.

- **SDRF-Proteomics experimental design**: `pymzlib.sdrf.read()` reads one `.sdrf.tsv`, and
  `pymzlib.sdrf.pool()` merges several into one analysis table with a `comment[source document]`
  column recording provenance. This is the first module here that answers *what was searched*
  rather than *what the search found*, which is what makes results from two experiments
  comparable at all.

  **It is row-major where every other reader is columnar**, and that is forced by the format: SDRF
  column names are data rather than a schema and they repeat - 649 files in the curated corpus
  carry `comment[modification parameters]` more than once, up to eight times in one file - so a
  name-keyed table would keep one occurrence and silently drop the rest. `columns` is a list that
  may contain duplicates, `rows` is a list of cell lists, and `value()` / `all()` reach them by
  position. Rows are ragged and stay ragged, cells are raw strings that are never interpreted, and
  a reserved word (`"not available"`) is reported as the real value it is, distinct from `None` for
  a column the document does not have.

  `pool()` takes a `{path: label}` mapping. Passing a plain list falls back to mzLib's
  `containing-folder/file-stem`, which depends on where the files sit, so the result carries a
  caveat saying it is not reproducible elsewhere. A partially-labelled set is refused outright.

  Validation is not offered yet. mzLib's `SdrfValidator` and `SdrfDriftLint` became public in
  mzLib #1207, which mzLib 1.0.589 includes, so the bridge can call them. They will be exposed
  from mzLib rather than reimplemented, because reimplementing a specification's rules once per
  binding is how three copies drift apart.

- **DIA-NN and SDRF are readable**, following the pin to mzLib 1.0.585. `DiaNnReport` is the
  fourth format offering the `quantifiable` view, so DIA data can now feed
  `pymzlib.flashlfq.quantify()`; `Sdrf` (HUPO-PSI experimental design, `.sdrf.tsv`) became a
  recognised type. This takes the supported count from 29 to 31 (mzLib #1120, #1138).

  *Corrected since:* SDRF is **not** usefully readable through `read_records()`, which
  semicolon-joins each row into one unsplittable string. Use `pymzlib.sdrf.read()`.

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
- **A PRIDE project's file manifest could be silently truncated.** `pymzlib.pride.list_files()`
  stopped as soon as it held as many files as the server's `total_records` header claimed, so a
  server understating that count had the tail of its manifest dropped with no error - and
  `file_count` and the total size were computed from what arrived, so the answer looked like a
  smaller project rather than a failure. Fixed upstream in mzLib #1173 and pinned by a regression
  test here.
- **Negative-mode MGF precursors reported a positive charge.** An MGF `CHARGE` line carries its
  sign as a trailing character (`CHARGE=2-`), which mzLib's reader dropped, so `read_spectra()`
  published charge `+2` and `Positive` polarity for a negative-mode scan. A neutral mass computed
  from those was wrong by two proton masses and looked entirely ordinary. Fixed upstream in mzLib
  #1164, pinned by a regression test here.
- **`quantify(match_between_runs=True)` was not reproducible.** FlashLFQ built its PEP training
  rows in a nondeterministic order, so the same inputs could give different MBR results between
  runs. Fixed upstream in mzLib #1155, which carries its own determinism test.

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

[Unreleased]: https://github.com/smith-chem-wisc/pyMzLib/compare/v0.1.0...main
[0.1.0]: https://github.com/smith-chem-wisc/pyMzLib/releases/tag/v0.1.0