# Changelog

Notable changes to pyMzLib. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow semantic versioning judged on the **Python** API — a change to the internal JSON
envelope is not a breaking change unless Python callers can see it.

## [Unreleased]

### Fixed
- **`peptidoform.fragments()`: `peptides_at_isoform_cap` counts the raw digest, as it always said
  it did.** It was recomputed on the de-duplicated list, so a locus that mzLib truncated at the cap,
  and that the mzLib#1108 de-duplication then dropped below it, was reported as untruncated. The
  correct count had been computed and left unused.
- **`pride.list_ftp_files()` raises `ProjectNotFoundError` only for a missing project.** It re-mapped
  every mzLib error to "no such project, check for a typo", including a cyclic FTP listing that
  exceeded the depth limit. That now stays a `BridgeError`.
- **`pride.download_files()` accepts a selection whose first line carries a UTF-8 BOM.** The bridge
  kept the BOM, so the first name matched nothing and failed as "not in project", after every
  other selected file had downloaded.
- **Stale FlashLFQ warnings.** The `max_threads` docstring still said results change with
  threads. mzLib#1155 fixed that cause (mzLib#1111), and 1.0.592 includes it. The `median_polish`
  docstring and the FlashLFQ guide still said `QuantifiedProteins.tsv` labels disagree until a
  re-pin. mzLib#1129 is in the pin, so they agree; only the column order can differ.

## [0.2.0] - 2026-09-24

Built from mzLib 1.0.592. Adds `pymzlib.proteins`, SDRF validation, and reading many files in one call.

### Added
- **Reference fact tables rendered from the bridge's per-verb specs.** Parameters (with units and
  ranges), result fields (with units and what null means), errors, caveats, the wrapped mzLib code
  and the same verb's spelling in mzLibRust and mzLibR now come from one spec per wire verb,
  vendored in `docs/specs/`. They are rendered into `docs/reference/`, starting with `readers
  formats`, `read-records` and `read-spectra`. CI fails when a table is stale, or when a docstring
  lacks a spec param or field, or its unit.
- **The docstring examples run in CI.** Each runs against a replay bridge that answers from
  fixtures recorded from the real bridge. Before this change none ran. Now 54 example lines run;
  20 are still skipped, each with its reason (the network, a download, pandas, or no recording
  yet). Running them corrected two examples: the PXD000001 total size was shown for RAW files but
  computed over all files, and the census text had changed.
- **A weekly external link check** (`links.yml`, lychee).
- **Read many files in one call** (bridge `design/BULK.md`). Every reader has a `_many` twin -
  `read_spectra_many()`, `read_records_many()`, `read_results_many()`, `read_features_many()`,
  `read_matches_many()`, `identify_many()`, and the three below - that takes a list of paths and
  returns one long table (`ReadBatch`) whose first columns are `source_index` and `source_path`,
  with each file's facts in `.files`. One bridge process reads the whole list, `threads` files at a
  time; the table is byte-identical at any `threads`, which defaults to 1 because each file in
  flight is held whole in memory. `on_error="skip"` records an unreadable file in its `FileReport`
  and reads the rest; `out=` streams the table to disk one file at a time. pyMzLib deliberately
  has no thread pool of its own.
- **`read_protein_groups()`, `read_quantified_peptides()` and `read_occupancy()`** (mzLib #1347):
  MetaMorpheus's `AllQuantifiedProteinGroups.tsv` as one row per group per sample group (spectral
  count, intensity), FlashLFQ's `QuantifiedPeptides.tsv` as one row per peptide per sample
  (intensity, detection type), and the protein-group table's PTM site occupancy as one row per
  site. Long rather than wide, so the columns are the same in every experiment.
- **`absent_fields` on every reader result**: the columns this file's format has no source for.
  Each is `None` in every row, whatever default mzLib filled in: `mbr_score` on a current FlashLFQ
  peaks table (#1345), `q_value` on an MSPathFinder targets file (where mzLib reads a "perfect" 0),
  apex `intensity` on a FLASHDeconv feature file. Every result also now carries `reader`,
  `failed_fields` and `excluded_fields`, and each `excluded_fields` entry names the `verb` that
  does carry the field.
- **`read_spectra()` reports the run**: `source` has the instrument model and its accession, the
  serial number, and the acquisition start time with whether it is UTC (mzLib #1349).
- **mzIdentML confidence in `read_matches()`**: `q_value`, `rank` and `pass_threshold` columns
  (#1306), `scores=True` for each engine's scores as long rows, and `skipped_count`/`skipped` for
  the items mzLib could not represent (#1313).
- **`bridge_version()["verbs"]`** lists every command the bridge dispatches, generated from its
  dispatch table at build time. A function newer than the bridge in use now fails with a
  `UsageError` naming the pyMzLib release it needs, instead of "Unknown command".
- **SDRF: validate, lint, assess, samples and ages** (mzLib 1.0.592). Five functions, each one mzLib
  type projected, and none re-deriving a rule:
  `sdrf.validate()` (`SdrfValidator` — structural findings with severity, rule, line and column),
  `sdrf.lint()` (`SdrfDriftLint` — concepts several files wrote differently, one row per spelling),
  `sdrf.assess()` (`SdrfSampleInformativeness`, #1325 — `Informative`, `Partial` or `Skeleton`, with
  the per-column counts behind it), `sdrf.samples()` (`SdrfSampleBlock`, #1334 — one row per sample
  x column, a column the sample's rows disagree about withheld and named) and `sdrf.parse_ages()`
  (`SdrfAge`, #1326/#1333 — ages in **years** with a precision, a bare `63` refused rather than
  assumed). `samples()` carries the parsed age on each `characteristics[age]` row. A lower bound
  (`>=90Y`) has no upper end: `max_years` is `None` and `precision` is `LowerBound`; a refused cell
  names why in `refusal`. See the [SDRF guide](https://smith-chem-wisc.github.io/pyMzLib/guides/sdrf/).
- **`validate_many()`, `assess_many()`, `samples_many()`** read a whole corpus in one bridge call,
  with `threads` (default 1, `-1` = every core; the result is identical at any value) and
  `on_error="fail"|"skip"`. A skipped file keeps its entry, with a `FileError`.
- **`sdrf.lint()` ignores `comment[searched data file]`** (mzLib #1335), which `read()` and `pool()`
  carry verbatim like any column: file names that differ between studies are data, not drift.
- **`pymzlib.proteins`, a new module for protein databases** (mzLib 1.0.592). Guide:
  [Protein databases](docs/guides/proteins.md).
  - `read()` loads UniProt XML or FASTA (decoys off) and returns one row per protein: organism,
    NCBI taxonomy id (XML, or FASTA `OX=`), gene names, length, monoisotopic mass (Da, unmodified
    sequence). On request, GO terms (mzLib #1336: aspect, term, ECO evidence codes) and Ensembl
    transcript links as long tables. Filter by accession; misses are listed.
  - `resolve_genes()` resolves proteins to stable Ensembl gene ids with mzLib's
    `EnsemblGeneResolver` (#1338), against a GTF you supply, with an outcome per protein and the
    sha256 of every input on every row. Optional Ensembl xref for a second opinion.
  - `classify_peptides()` classifies peptides as `Unique`, `SharedWithinGene`,
    `SharedAcrossGenes` or `NotInDatabase` with mzLib's `PeptideUniquenessClassifier` (#1348),
    I and L treated as one residue.
  - All three read one database or many in one call (`threads`, default 1, same answer at any
    value), with `contaminants=` marked. A FASTA's lack of GO and Ensembl data is reported in
    `absent_fields`, not left as an empty table.

### Fixed
- **Docstring examples render as code in the API reference.** They sat under `Example:`, which
  griffe reads as an admonition, so the site showed every `>>>` line as nested blockquotes. They are
  now under `Examples:`, and a test keeps them there.

### Changed
- **The bridge is built from mzLib 1.0.592** (was the 8931f219 commit, via 1.0.591), so mzLib's
  fixes to MGF, mzML and mzIdentML reading and writing, semi-specific digestion and RNA databases
  arrive with it.
- **Four more formats: 36, up from 32** (mzLib 1.0.592). `read_records()` and `identify()` now take
  mzIdentML (`.mzid`, and `.mzid.gz` read without unpacking it, mzLib #1313), MetaMorpheus's
  `AllQuantifiedProteinGroups.tsv` and its `AllQuantifiedPeptides.tsv`, which FlashLFQ also writes
  (#1347). Each was "file type not supported" before. mzIdentML also offers the `spectral_match`
  view, so `read_matches()` reads it; the two quantification tables have no uniform view. Their
  per-sample values (`sample_groups`, `samples`) and mzIdentML's engine `scores` are dictionaries,
  which `read_records()` names in `excluded_fields`; the functions above project them.
- **`read_matches()` on mzIdentML carries its own caveats**: `is_decoy` is `None` because the
  file's `isDecoy` attribute defaults to false when a writer omits it; every identification item
  is a row, not only accepted ones; items mzLib cannot represent (crosslinks, unresolvable
  modifications, substitutions) are skipped, and `skipped` lists each with its reason.
- **`excluded_fields` names a read-only dictionary as a dictionary.** It said "a list of composite
  values" for any `IReadOnlyDictionary`, which described the new readers' per-sample fields wrongly.
- **`read_records()` on an older MetaMorpheus `.psmtsv`/`.osmtsv` now fills `pro_forma`** (mzLib
  #1346). Files from MetaMorpheus 1.1.11 and earlier have no ProForma column, so it was always
  `None`; mzLib now converts `full_sequence`. Ambiguous (`|`-joined) or unconvertible rows are still
  `None`. It costs time: a 271,551-row, 295 MB file read in 15.6 s, up from 12.9 s (+21%).
- **`read_records()` on a FlashLFQ `QuantifiedPeaks.tsv`: `mbr_score` is `None` where it was `0`**
  (mzLib #1345), for every peak without a score. Current FlashLFQ and MetaMorpheus 1.1.11 peaks
  tables, which have no `MBR Score` column, read now instead of failing with
  `HeaderValidationException`, and seven columns join the record: `organism`, `peak_fwhm`,
  `peak_fwhm_status`, `pip_q_value`, `pip_pep`, `decoy_peptide` and `random_rt`.
- **The psmtsv caveat's line citations moved** to `SpectrumMatchFromTsv.cs:119` and `:194`, where
  mzLib #1346 left the lines they describe. The claims themselves are unchanged.
- **`read_records()` reads Pytheas match output**, a 32nd format (mzLib #1277). It has no uniform
  view. A `.txt` file is read as Pytheas when a `#theoretical_digest` line appears in its first five
  lines, and as Crux otherwise.
- **MGF: `one_based_precursor_scan_number` is no longer always null.** mzLib #1227 reads it from
  `PRECURSORSCAN`, an mzLib extension header. It is still null on MGF files mzLib did not write. The
  caveat now says so. It also said scan numbers came from the `TITLE` line, which was wrong: they
  come from `SCANS`, or from file order.

## [0.1.1] - 2026-09-18

Fixes Thermo `.raw` reading, which failed on every platform in 0.1.0.

### Changed
- **Every file-path argument accepts a `pathlib.Path` as well as a string.** `pride.download`
  returns `Path` objects, but the readers, `sdrf` and `flashlfq` refused them with "A file path is
  required", so its result could not go straight into `read_spectra` without wrapping it in `str()`.
  Any `os.PathLike` now works for input files, `out=` and `output_directory=`, and in FlashLFQ's
  `spectra` list. `bytes` paths are still refused.

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

[Unreleased]: https://github.com/smith-chem-wisc/pyMzLib/compare/v0.2.0...main
[0.2.0]: https://github.com/smith-chem-wisc/pyMzLib/releases/tag/v0.2.0
[0.1.1]: https://github.com/smith-chem-wisc/pyMzLib/releases/tag/v0.1.1
[0.1.0]: https://github.com/smith-chem-wisc/pyMzLib/releases/tag/v0.1.0