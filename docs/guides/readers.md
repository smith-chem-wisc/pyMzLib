# Readers

**Spectra files are read here, not just search output.** [`read_spectra()`](#read_spectra-scans-and-peaks)
reads **mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign: scan headers always,
peaks on request. Alongside those, mzLib maintains a parser for the output of a dozen search and
deconvolution tools (MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo,
FlashDeconv, Dinosaur, DIA-NN, FlashLFQ, Pytheas, and any mzIdentML writer). pyMzLib lets you point
at a file, ask what it is, and read it.

| You want to | Call | mzLib does it with |
|---|---|---|
| Know what a file is, and which views it offers | [`identify()`](#start-with-views-not-with-the-file-type) | `SupportedFileTypes.ParseFileType` |
| Read **everything** a file has, under mzLib's names | [`read_records()`](#read_records-any-format-its-own-fields) | the format's own `ResultFile<T>` |
| Compare search results across tools, or feed FlashLFQ | [`read_results()`](#read_results-the-quantifiable-view) | `IQuantifiableRecord` |
| Read deconvolved MS1 features | [`read_features()`](#read_features-deconvolved-ms1-features) | `ISingleChargeMs1Feature` |
| Compare identifications from MsPathFinderT, Casanovo or mzIdentML | [`read_matches()`](#read_matches-identifications) | `ISpectralMatch` |
| Read scans and peaks | [`read_spectra()`](#read_spectra-scans-and-peaks) | `MsDataFile` |
| Read per-sample protein groups, peptides or PTM occupancy | [`read_protein_groups()`, `read_quantified_peptides()`, `read_occupancy()`](#quantification-tables-protein-groups-peptides-and-occupancy) | `ProteinGroupFromTsvFile`, `QuantifiedPeptideFile` |
| Read many files into one table | [`*_many()`](#many-files-at-once) | the same readers, in one bridge process |

Every `>>>` example on this page runs in CI against output recorded from the real bridge, on the
small files mzLib ships with its own tests. Blocks titled **Not run** say why.

```pycon
>>> import pymzlib
>>> info = pymzlib.readers.identify("FraggerPsm_FragPipev21.1_psm.tsv")
>>> info.file_type, info.views
('MsFraggerPsm', ['quantifiable'])
>>> t = pymzlib.readers.read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv")
>>> t.record_type, t.record_count
('ToppicPrsm', 4)

```

**Every format mzLib recognises is readable.** What differs between them is not whether you can
read them but what the columns mean, which is the subject of this page.

## Ways to read, and how to choose

There is one universal function, four cross-format views, and three functions for mzLib's
quantification tables. The choice is a real one:

| function | reads | columns | use it when |
|---|---|---|---|
| [`read_records()`](#read_records-any-format-its-own-fields) | **every format** | **this format's own fields**, under mzLib's names | you want *everything* a file has |
| [`read_results()`](#read_results-the-quantifiable-view) | the `quantifiable` view | uniform: sequence, RT, charge, mass, proteins | you are feeding [FlashLFQ](flashlfq.md) or comparing search results |
| [`read_features()`](#read_features-deconvolved-ms1-features) | the `ms1_features` view | uniform: m/z, charge, RT range, intensity | you are working with deconvolved MS1 features |
| [`read_matches()`](#read_matches-identifications) | the `spectral_match` view | uniform: scan, sequences, accession, mods | you are comparing identifications from MsPathFinderT, Casanovo or mzIdentML |
| [`read_spectra()`](#read_spectra-scans-and-peaks) | the `spectra` view | uniform: scan headers, peaks on request | the file is spectra rather than results |
| [`read_protein_groups()`](#quantification-tables-protein-groups-peptides-and-occupancy) | MetaMorpheus protein-group tables | long: one row per protein group per sample group | you want MetaMorpheus's per-sample protein (or RNA transcript) intensities and spectral counts |
| [`read_quantified_peptides()`](#quantification-tables-protein-groups-peptides-and-occupancy) | FlashLFQ peptide tables | long: one row per peptide per sample | you want FlashLFQ's per-sample peptide (or oligonucleotide) intensities |
| [`read_occupancy()`](#quantification-tables-protein-groups-peptides-and-occupancy) | MetaMorpheus protein-group tables | long: one row per modified site | you want PTM site occupancy from a MetaMorpheus protein-group (or transcript-group) table |

**Every one of them has a `_many` twin** (`read_spectra_many()`, `read_records_many()`,
`identify_many()` and so on) that reads a list of files into one table in one call. See
[Many files at once](#many-files-at-once).

The rule of thumb:

> **A typed view when you need numbers that mean the same thing across files.
> `read_records()` when you need everything one file has.**

Both matter. A `.psmtsv` read through `read_results()` gives you a handful of columns you can
safely compare against an MSFragger file. The same `.psmtsv` through `read_records()` gives you
dozens more, including the q-values, PEP and scores the uniform view does not carry, but under
MetaMorpheus's names, which no other format has.

## Start with `views`, not with the file type

It would be convenient if mzLib read every format into one uniform table. **It does not.** The
formats fall into disjoint families, and many belong to no family at all. `formats()` says which,
straight from mzLib:

```pycon
>>> from collections import Counter
>>> Counter(view for f in pymzlib.readers.formats() for view in (f.views or ["(none)"]))
Counter({'(none)': 19, 'spectra': 7, 'spectral_match': 6, 'quantifiable': 4, 'ms1_features': 2})

```

| view | what it means | which formats |
|---|---|---|
| `quantifiable` | a cross-format record view: sequence, retention time, charge, mass, protein groups. What [`flashlfq.quantify()`](flashlfq.md) accepts. | MetaMorpheus `.psmtsv`/`.osmtsv`, MSFragger `psm.tsv`, DIA-NN `report.tsv` |
| `ms1_features` | deconvolved MS1 features | TopFD `_ms1.feature`, Dinosaur |
| `spectral_match` | records are identifications, but share no *file*-level interface | MsPathFinderT's three files, Casanovo, mzIdentML `.mzid`/`.mzid.gz` |
| `spectra` | the file is spectra, not results | `.raw`, `.mzML`, `.mgf`, both `.d` types, both msalign types |
| *(none)* | mzLib parses it into a format-specific shape with nothing in common | TopPIC's four tables, Crux, Pytheas, MSFragger peptide/protein, FlashDeconv, MetaMorpheus protein groups and peptides and their RNA counterparts, and more |

`views == []` is a real and common answer, not an error: it is the commonest answer. It means
"mzLib reads this, but there is no uniform projection of it", and `read_records()` is exactly the
function for that case:

```python title="Not run: a pattern for your own files"
info = pymzlib.readers.identify(path)

if info.is_quantifiable:
    table = pymzlib.readers.read_results(path)     # comparable columns
else:
    table = pymzlib.readers.read_records(path)     # this format's own columns
```

## `read_records()`: any format, its own fields

This is the exhaustive verb. If `identify()` succeeds on a path, this reads it:

```pycon
>>> t.file_type, t.views
('ToppicPrsm', [])
>>> t.column_names[:4]
['file_name_without_extension', 'file_path', 'prsm_id', 'spectrum_id']
>>> "e_value" in t.column_names
True

```

Column names are mzLib's own property names converted to `snake_case`, so they are
**cross-referenceable against the mzLib source**: a column called `e_value` is `ToppicPrsm.EValue`,
and `record_type` tells you which class to look in. Acronyms survive the conversion intact:
`EValue` becomes `e_value`, `MIScore` becomes `mi_score`, `FixedPTMs` becomes `fixed_ptms`.

```python title="Not run: needs pandas, which pyMzLib does not depend on"
import pandas as pd
frame = pd.DataFrame(t.columns)
frame[frame.e_value < 1e-10][["base_sequence", "protein_accession", "e_value"]]
```

### Nothing is silently dropped

Some fields cannot become a column. A nested object or a dictionary has no faithful column shape,
and flattening one would mean publishing a schema mzLib does not have. Those fields are **named,
with the reason**, rather than quietly omitted, because a column that simply vanished cannot be
told apart from a field the format does not have.

Some exclusions carry the numbers you came for. The per-sample values of mzLib's quantification
tables are dictionaries keyed by sample, so `read_records()` names them and does not project them,
and each entry names the verb that does:

```pycon
>>> g = pymzlib.readers.read_records("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv", limit=2)
>>> [(e["field"], e["verb"]) for e in g.excluded_fields]
[('sample_groups', 'readers read-protein-groups')]

```

| file type | excluded field | what it holds | read it with |
|---|---|---|---|
| `MetaMorpheusQuantifiedProteinGroups` | `sample_groups` | per-sample intensity, spectral count and modification occupancy | [`read_protein_groups()`](#quantification-tables-protein-groups-peptides-and-occupancy), and [`read_occupancy()`](#site-occupancy) for the occupancy |
| `FlashLFQQuantifiedPeptide` | `samples` | per-run intensity, detection type and retention time | [`read_quantified_peptides()`](#quantification-tables-protein-groups-peptides-and-occupancy) |
| `MzIdentML`, `MzIdentMLGz` | `scores` | the search engine's own scores, e.g. `MS-GF:SpecEValue` | [`read_matches(scores=True)`](#engine-scores-as-long-rows) |

Their scalar fields (protein group name, gene, organism, q-value, sequence) are columns as usual.

`failed_fields` is the other half. Several mzLib properties are *computed* and assume a
UniProt-style FASTA header: Crux's and MsPathFinderT's `accession` are both
`protein_id.split("|")[1]`, so on a database with plain headers they throw. Those cells arrive
`None` rather than taking the whole file down, and the field is named, so a failure never looks like
missing data. mzLib's Crux test file was searched against UniProt headers, so nothing failed:

```pycon
>>> crux = pymzlib.readers.read_records("crux.txt", limit=3)
>>> crux.failed_fields
[]

```

On a non-UniProt database the same read reports `['accession: IndexOutOfRangeException']`.

!!! tip "`-1` is not treated as missing here"
    `read_results()` maps mzLib's documented `-1` "absent" sentinel to `None` for two specific
    interface fields. `read_records()` deliberately does **not** generalise that, because in a
    format's own columns `-1` is often a real measurement: a mass difference, a delta, a log ratio,
    TopPIC's `feature_score`. Nulling those would destroy data. Non-finite values (`NaN`, infinity)
    still cross as `None`, since JSON cannot carry them.

### Four ways a field can have no value

A `None` in a table can mean four different things, and each is named, so you never have to guess
which:

| where it is named | what it means | what you see in the table |
|---|---|---|
| `absent_fields` | the function defines the field, but **this file's format has no column for it** | `None` in every row |
| `failed_fields` | the field exists, but **reading it threw** on some rows | `None` in those rows |
| `excluded_fields` | the field has **no column shape** (a dictionary, a nested object), or would cross only lossily (an SDRF row's `header` and `cells`; use [`sdrf.read()`](sdrf.md)) | not a column at all |
| none of them | the value is **genuinely missing for that row**: the precursor of an MS1 scan, a blank cell | `None` in that row |

`absent_fields` is the one that protects numbers. When a file lacks an optional column, mzLib does
not say so: it fills in its default, and for a number that default is usually zero. pyMzLib reads
mzLib's own declaration of which columns are optional, checks it against the file's header, and
reports and blanks every column the file does not have:

```pycon
>>> peaks = pymzlib.readers.read_records("FlashLFQ_MzLib1.0.591_QuantifiedPeaks.tsv", limit=2)
>>> peaks.absent_fields                       # the table has no MBR Score column (mzLib#1345)
['mbr_score']
>>> m = pymzlib.readers.read_matches("MsPathFinderT_WithMods_IcTda.tsv")
>>> m.absent_fields
['q_value', 'rank', 'pass_threshold']
>>> m.columns["q_value"][:3]                  # mzLib would have said 0.0: "perfect"
[None, None, None]

```

The MsPathFinderT case is why this exists. A file without a `QValue` column (its `_IcTarget.tsv` is
written *before* target-decoy analysis) still reads, and mzLib types that property as a plain
`double`, so it reads **0**, a perfect q-value, for every match. Filtering on `q_value <= 0.01`
would keep everything. A column in `absent_fields` is always `None`, whatever mzLib filled in.

Every function reports all three lists, so a check can be written once:

```pycon
>>> def trustworthy(result, column):
...     return column not in result.absent_fields and not any(
...         f.startswith(column + ":") for f in result.failed_fields)
>>> trustworthy(m, "q_value"), trustworthy(m, "sequence")
(False, True)

```

Values that changed between mzLib releases, for the same file, are on
[Upgrading: what changed in your results](../upgrading.md).

## `read_results()`: the quantifiable view

The formats offering `quantifiable` (MetaMorpheus `.psmtsv` and `.osmtsv`, MSFragger `psm.tsv`,
DIA-NN `report.tsv`) read into one fixed shape that is safe to compare between files:

```pycon
>>> r = pymzlib.readers.read_results("FraggerPsm_FragPipev21.1_psm.tsv", limit=2)
>>> sorted(r.columns)[:5]
['base_sequence', 'charge_state', 'file_name', 'full_sequence', 'gene_name']
>>> r.retention_time_unit
'minutes'

```

Data comes back **columnar**: one list per field, rather than one object per record. pyMzLib has
no third-party dependencies, so it can never hand you a DataFrame; a map of lists is the one shape
that becomes one in a single call (`pd.DataFrame(r.columns)`). If you would rather loop, `r.records`
gives the same data as one dict per row.

### Nothing is ever silently short

There is **no default row limit** on any of the reading functions. A result file can carry a
million rows, and a library whose default answer is "here's some of it" eventually puts a truncated
table in a paper. Ask for a limit and you are told when it bites:

```pycon
>>> r.returned_count, r.record_count, r.truncated
(2, 5, True)

```

`rows_not_read` is the other half of that promise. mzLib drops a malformed row silently: it
collects a warning per unreadable line and the reader discards the list, so a half-corrupt file
reads "successfully" with fewer rows than it contains. pyMzLib counts the difference and reports
it. This file is whole:

```pycon
>>> r.rows_not_read
0

```

### Large files: write, don't page

```python title="Not run: writes a file; the recordings here were made without out="
r = pymzlib.readers.read_records("huge_prsm.tsv", out="records.tsv")
r.output.path, r.output.row_count          # ('C:/.../records.tsv', 842130)
```

The table goes to disk and the envelope carries only a summary. It is **tab-separated**, because
these fields contain commas (MSFragger's mapped proteins are a comma-separated list inside a single
field) and because every mzLib reader and writer uses tabs. Read it with
`pandas.read_csv(path, sep="\t")`, or `csv.reader(f, delimiter="\t")` with no dependencies at all.

`out=` must end in `.tsv`, in any case. Any other extension, or none, is a `UsageError` raised
before the file is read, and no extension is ever added for you: `out="records"` is refused, not
written to `records.tsv`. The same rule holds for every `read_*()` function and its `_many` twin.

A read too large to come back as one answer is a `UsageError` naming its record count; add
`limit=` or write it with `out=`. The answer is one JSON document held in one .NET string, which
holds about 1.07 billion characters however much memory the machine has.

!!! warning "`offset` is a window, not a cursor"
    mzLib's readers look lazy and are not: every one of them reads the whole file into a list. So
    `offset` does not resume where you left off; it re-reads and re-parses the entire file and then
    skips. Paging a large file is quadratic. Use `out=` instead.

## `read_features()`: deconvolved MS1 features

TopFD/FLASHDeconv `_ms1.feature` and Dinosaur `.feature.tsv` offer the `ms1_features` view. The
columns are `mz`, `charge`, `retention_time_start`, `retention_time_end`, `intensity` and
`number_of_isotopes`.

```pycon
>>> f = pymzlib.readers.read_features("Ms1Feature_TopFDv1.6.2_ms1.feature", limit=5)
>>> f.record_count, f.retention_time_unit
(25, 'unknown')
>>> f.absent_fields
['number_of_isotopes']

```

!!! warning "One row is not one line of the file, for `_ms1.feature`"
    An `_ms1.feature` row is a deconvolved **neutral mass spanning a charge range**, and mzLib
    expands it into one single-charge feature per charge in `[ChargeStateMin, ChargeStateMax]`. A
    hundred-feature file can read as a thousand rows. Dinosaur is one-for-one. Either way
    `read_records()` gives you the file's own rows.

`intensity` is the **apex** intensity, not the sum over the feature. Both formats carry a summed
intensity column too, and `read_records()` has it.

!!! warning "`intensity` is `None` for every FLASHDeconv `_ms1.feature`"
    mzLib takes the per-charge intensity from `Apex_intensity`, an *optional* column that the
    FLASHDeconv/OpenMS `_ms1.feature` layout does not have, and substitutes **zero** when it is
    absent. A whole column of zeros cannot be told apart from real measurements of nothing, so
    pyMzLib crosses those as `None`, names `intensity` in `absent_fields`, and says why in
    `caveats`. TopFD files, which do write the column, are unaffected. `read_records()` has the
    file's own summed `intensity` either way.

### `retention_time_unit` is `'unknown'` for `_ms1.feature`, and that is the honest answer

TopFD wrote retention times in **seconds** through v1.6.2 and in **minutes** from v1.7.0, *within the
same file type*, with nothing in the file to tell you which. mzLib normalises neither, and its own
deconvolution code falls back on a heuristic (divide everything by 60 if the largest end time
exceeds 500). pyMzLib will not launder a guess into a stated fact:

```pycon
>>> f.retention_time_start_in_minutes                      # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
pymzlib._bridge.UsageError: Cannot convert retention_time_start for 'Ms1Feature' (...): mzLib gives no basis to say what unit it is in. ...

```

Dinosaur reports `'minutes'` and converts without complaint.

## `read_matches()`: identifications

MsPathFinderT's targets, decoys and combined results, Casanovo's `.mztab`, and mzIdentML (`.mzid`,
and `.mzid.gz` read without unpacking it) offer the `spectral_match` view. These are the
identification formats that share no *file*-level interface, so `read_results()` cannot reach them.

```pycon
>>> m.file_type, m.record_count
('MsPathFinderTAllResults', 5)
>>> sorted(m.columns)[:4]
['accession', 'base_sequence', 'file_name_without_extension', 'full_sequence']

```

Beside the identity fields, the view carries the three confidence fields a format may record:
`q_value` (mzIdentML's PSM-level q-value, mzLib#1306; MSPathFinder's `QValue`), and mzIdentML's
`rank` and `pass_threshold`. A format without one names it in `absent_fields`, as above, and the
column is `None` throughout.

!!! danger "Nothing here is FDR-filtered"
    Every row the file holds is a row here. Filter on `q_value`, where `absent_fields` does not name
    it, and, for mzIdentML, on `rank == 1` and `pass_threshold`, before you count or report
    anything. The engines' other scores are in `read_records()`, or as long rows with
    [`scores=True`](#engine-scores-as-long-rows).

Three `is_decoy` traps, all reported in `caveats`:

- **MsPathFinderT** infers decoys from the protein *name*: mzLib reports a decoy when `ProteinName`
  starts with `XXX`. A database whose decoys carry a different prefix reads **entirely as targets**.
- **Casanovo** is de novo and writes no target/decoy label at all. mzLib's record leaves the field
  at its default `False` and never assigns it, so `False` would mean *unknown*. pyMzLib crosses it
  as `None` instead, the same rule `read_results()` already applies to MSFragger.
- **mzIdentML**'s `isDecoy` attribute is optional and defaults to false, so a writer that omits it
  reads as a target. `read_matches()` crosses `None`; `read_records()` carries mzLib's own boolean
  for when you know your writer sets it.

```pycon
>>> mz = pymzlib.readers.read_matches("PXD078927_msgf_1_1_0.mzid", limit=3)
>>> mz.absent_fields, mz.columns["is_decoy"]
(['is_decoy'], [None, None, None])

```

Casanovo also numbers scans by mzTab **index**, not by the instrument's scan number; when Casanovo
was run on an MGF the two are unrelated, so do not join on it.

mzIdentML has three more things to know, also in `caveats`:

- **Every identification item is a row**, not only the ones the submitter accepted. Lower-ranked
  candidates and items that failed the threshold are included; filter on `rank` and
  `pass_threshold`.
- **Some items are skipped, not read**: crosslinks, modifications without a resolvable UNIMOD
  accession, substitutions, and two modifications on one residue (mzLib#1313). They are reported:
  `skipped_count` says how many, and `skipped` lists each with its reason, so `record_count +
  skipped_count` is the number of items in the file. This file skipped none:

  ```pycon
  >>> mz.record_count, mz.skipped_count
  (12, 0)

  ```
- **Scan numbers come from the nativeID.** `scan=N` gives N, but `index=N` (peak-list input) is a
  zero-based position and gives N + 1, which is not an instrument scan number.

### Engine scores as long rows

mzIdentML carries each search engine's own scores (`MS-GF:SpecEValue`, `Mascot:score`,
`Scaffold:Peptide Probability`), and no two engines share a name. There is no fixed set of columns
to put them in, so `scores=True` makes the table **long** instead: one row per match *and* score,
with `match_index` (the match's position in its file), `score_name` and `score_value`.

```pycon
>>> s = pymzlib.readers.read_matches("PXD078927_msgf_1_1_0.mzid", scores=True, limit=1)
>>> s.returned_count, s.row_count                 # one match, its scores as rows
(1, 7)
>>> "MS-GF:SpecEValue" in s.columns["score_name"]
True

```

`returned_count` still counts matches, the unit `limit` and `offset` count in, and `row_count`
counts rows. A format with no engine scores keeps one row per match and names `score_name` and
`score_value` in `absent_fields`. To compare scores, pivot on `match_index`:

```python title="Not run: needs pandas, which pyMzLib does not depend on"
frame = pd.DataFrame(s.columns)
scores = frame.pivot(index="match_index", columns="score_name", values="score_value")
scores["MS-GF:SpecEValue"].lt(1e-10).sum()
```

## `read_spectra()`: scans and peaks

The `spectra` view's retention times **are** in minutes for every format: mzLib's spectra readers
convert at the boundary, unlike its result-file readers.

```pycon
>>> sp = pymzlib.readers.read_spectra("sliced_ethcd.mzML", ms_order=2, limit=2)
>>> sp.scan_count, sp.record_count, sp.returned_count       # file total, filtered, returned
(6, 5, 2)
>>> sp.retention_time_unit
'minutes'

```

`ms_order` filters **before** the offset/limit window, so `ms_order=2, limit=10` means the first ten
MS2 scans rather than the MS2 scans among the first ten. `scan_count` always reports the file's real
total, so a filter that matched nothing can never look like an empty file.

### Peaks are opt-in

A scan header is tens of bytes; its peak list is thousands, and a mid-size mzML holds tens of
thousands of scans. Returning peaks by default would make the ordinary "what is in this file?" call
serialise hundreds of megabytes.

```pycon
>>> one = pymzlib.readers.read_spectra("sliced_ethcd.mzML", peaks=True, ms_order=1, limit=1)
>>> len(one.columns["mz"][0]), len(one.columns["intensity"][0])     # one list per scan
(484, 484)
>>> one.columns["peak_count"][0]
484

```

Without `peaks=True`, `peak_count` still tells you how many peaks each scan has.

### What the file says about the run

Every spectra read reports what the file records about where and when it was acquired, from mzLib's
`SourceFile` (mzLib#1349):

```pycon
>>> src = sp.source
>>> src.instrument_model, src.instrument_model_accession, src.instrument_serial_number
('Orbitrap Fusion', 'MS:1002416', 'FSN10189')
>>> src.acquisition_start_time, src.acquisition_start_time_is_utc
('2021-03-16T17:09:07Z', True)
>>> src.acquired_at
datetime.datetime(2021, 3, 16, 17, 9, 7, tzinfo=datetime.timezone.utc)

```

These are what a batch or instrument confound is built from: which unit of a model, and in what
order the runs were acquired. Three things to know:

- **Match instruments on the accession, not the name.** mzML carries both the model's name and its
  PSI-MS accession; Thermo `.raw` records only the name, so the accession is `None` there rather
  than looked up.
- **A time is UTC only when it says so.** `acquisition_start_time` ends in `Z`, and
  `acquisition_start_time_is_utc` is true, only when the file fixed the instant, as this mzML does.
  A Thermo `.raw` records the acquisition computer's local clock with no offset; `acquired_at` then
  returns a *naive* `datetime` rather than inventing a time zone. A `.raw` and ProteoWizard's mzML
  of it can disagree by the site's UTC offset, because ProteoWizard converts using the *converting*
  machine's time zone and writes `Z`.
- **MGF and msalign record none of it**: every field is `None`.

!!! info "Two of the spectra formats need Windows"
    Bruker `.d` and timsTOF `.d` are read through vendor native libraries (`baf2sql`, `timsdata`)
    and are **Windows-x64 only**. Thermo `.raw` uses managed vendor assemblies and works
    everywhere. msalign files hold **deconvolved neutral masses**, not raw m/z: do not
    re-deconvolve them.

## Quantification tables: protein groups, peptides and occupancy

MetaMorpheus and FlashLFQ write their quantification as wide tables, one column per sample:
`Intensity_QE-002106_GM1_a-calib`, `SpectralCount_QE-002106_GM1_a-calib`, and so on. mzLib reads
them (mzLib#1347), and pyMzLib returns them **long**: one row per record per sample, with the
sample in a column.

| function | reads | one row per | per-sample columns |
|---|---|---|---|
| `read_protein_groups()` | MetaMorpheus `AllQuantifiedProteinGroups.tsv`, `AllProteinGroups.tsv`, `<file>_ProteinGroups.tsv`; RNA `AllQuantifiedTranscriptGroups.tsv` | protein group × sample group | `spectral_count`, `intensity` |
| `read_quantified_peptides()` | FlashLFQ `QuantifiedPeptides.tsv`, MetaMorpheus `AllQuantifiedPeptides.tsv`; RNA `AllQuantifiedOligos.tsv` | peptide × sample | `intensity`, `detection_type`, `retention_time` (IsoTracker only) |
| `read_occupancy()` | any protein-group or transcript-group table `read_protein_groups()` reads | group × sample group × basis × modified site | `fraction`, `numerator`, `denominator` |

**Why long.** A wide table's column names would be made up from your sample labels, so they would
differ in every experiment and could not be documented, checked or shared across files, and two
searches with different samples could not be read into one table. A long table has the same
columns every time; it is what pandas, polars and every plotting library want, and the wide matrix
is one `pivot` away.

```pycon
>>> pg = pymzlib.readers.read_protein_groups("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv", limit=1)
>>> pg.record_count, pg.returned_count, pg.row_count       # groups, groups returned, rows
(6, 1, 18)

```

`record_count` and `returned_count` count **groups** (or peptides), the unit `limit` and `offset`
count in, and `row_count` counts rows. The group's other fields (coverage, masses, member counts)
are in `read_records()` on the same file; join on `protein_group_name`.

```python title="Not run: needs pandas, which pyMzLib does not depend on"
import pandas as pd
df = pd.DataFrame(pymzlib.readers.read_protein_groups("AllQuantifiedProteinGroups.tsv").columns)

# The table is UNFILTERED, as MetaMorpheus writes it: filter before anything else.
confident = df[(df.q_value <= 0.01) & (df.decoy_contaminant_target == "T")]
matrix = confident.pivot(index="protein_group_name", columns="sample_label", values="intensity")
```

Four things these tables will not do for you:

- **Filter.** Decoys, contaminants and groups above 1% FDR are all rows; see above.
- **Tell a blank from a zero, unless you look.** A protein-group intensity that was not measured is
  a blank cell and arrives `None`. FlashLFQ's peptide table is different: it writes a literal **0**
  for a peptide it did not quantify, and mzLib keeps it. `detection_type` (`MSMS`, `MBR`,
  `NotDetected`, ...) is what separates a measurement from a missing value; filter on it before a
  mean or a log.
- **Know your design.** `sample_label` is the header label verbatim. Condition, replicate and
  channel cannot be recovered from it: map it yourself, from an [SDRF](sdrf.md) or a sample sheet.
- **Name a leading protein.** `protein_group_name` lists the members sorted by accession, so the
  first is only the one that sorts first.

### RNA searches, and searches without quantification

The same three functions read two more kinds of table:

- **An RNA search** writes `AllQuantifiedTranscriptGroups.tsv` and `AllQuantifiedOligos.tsv`
  (mzLib#1388). They are their own file types, but mzLib reads them with subclasses of the
  protein-group and peptide readers, so the columns keep their protein names:
  `protein_group_name` holds the transcript group, `sequence` the oligonucleotide, and
  `read_occupancy()` returns RNA modifications such as `2'-O-methyluridine on U`.
- **A search with label-free quantification off** writes `AllProteinGroups.tsv`, and each file's
  `<file>_ProteinGroups.tsv` (mzLib#1365). MetaMorpheus writes `Intensity_` columns only when it
  quantified, so these rows have `spectral_count` and no `intensity`: the field is named in
  `absent_fields`, which means *no basis*, not zero.

```pycon
>>> rna = pymzlib.readers.read_protein_groups("MetaMorpheus_RNA_AllQuantifiedTranscriptGroups.tsv")
>>> rna.record_count, rna.row_count
(3, 111)
>>> fluc = next(r for r in rna.records
...             if r["protein_group_name"] == "FLuc" and r["sample_label"] == "1:1_1")
>>> fluc["spectral_count"], round(fluc["intensity"])
(355, 89077471)

```

```python title="Not run: mzLib's tests ship no unquantified protein-group table to record from"
g = pymzlib.readers.read_protein_groups("AllProteinGroups.tsv")
"intensity" in g.absent_fields                      # True: this search was not quantified
```

### Site occupancy

MetaMorpheus writes two occupancy cells per group per sample group, one from PSM counts and one
from intensities, each a list of modified sites encoded as text. `read_occupancy()` returns one row
per site, with `basis` saying which cell it came from:

```pycon
>>> o = pymzlib.readers.read_occupancy("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
>>> o.record_count, o.row_count
(6, 95)
>>> sorted(set(o.columns["basis"]))
['count', 'intensity']

```

The two bases round differently, so trust the ratio of the counts on `count` rows and `fraction`
on `intensity` rows. `position` counts 0 as the protein N-terminus and residues from 1. A cell the
writer cut short keeps its complete sites with `cell_is_truncated` set, and `truncated_cell_count`
counts every cut cell, including those with no complete site left to show, so a short table is
never mistaken for a complete one. A cell that is not an occupancy cell at all is refused, named in
`failed_fields`, and gives no rows.

## Many files at once

Every reading function has a `_many` twin that takes a list of paths and returns **one** long
table: `read_spectra_many()`, `read_records_many()`, `read_results_many()`,
`read_features_many()`, `read_matches_many()`, `read_protein_groups_many()`,
`read_quantified_peptides_many()`, `read_occupancy_many()` and `identify_many()`.

```pycon
>>> batch = pymzlib.readers.read_records_many(
...     ["PXD078927_msgf_1_1_0.mzid", "PXD078927_msgf_1_1_0.mzid.gz"]
... )
>>> batch.read_count, batch.record_count
(2, 24)
>>> sorted(set(batch.columns["source_index"]))
[0, 1]

```

**When to use it.** Whenever you would otherwise write a loop over files. The list goes to one
bridge process, so the .NET start-up (about 120 ms) and the JIT warm-up are paid once rather than
once per file, and mzLib reads `threads` files at a time. For two hundred small files that is the
difference between most of a minute of start-up and none. pyMzLib has no thread or process pool of
its own, deliberately: the bridge is the one place that can count every thread in use.

**The table.** Its first two columns are `source_index` (the file's position in your list) and
`source_path`. Rows are grouped by file in list order, and within a file in the file's own order.
`batch.files[i]` is a `FileReport` with everything a single-file read would have told you about
file `i`: its `file_type`, `record_count`, `caveats`, `absent_fields`, `source`, `skipped`, and so
on. `batch.record_count` and `batch.row_count` are summed over the files.

**`threads`.** The result is **byte-identical at any `threads`**, which is tested, so it only trades
memory for speed and never changes an answer. It defaults to 1 because every mzLib reader holds a
whole file in memory while it works, so `threads=8` can mean eight whole files at once. Some
guidance:

- Many small files (search results, feature tables, MGFs): `threads` near your core count, or
  `-1` for one per core.
- Large spectra files: mzLib's mzML and Thermo readers already parallelise *inside* one file, so
  extra files in flight add memory faster than speed. Start at 2 and watch memory.
- Hundreds of files, or peaks: add `out=`; see below.

**`on_error`.** By default the first file that cannot be read stops the batch and raises, with a
message starting `Input <i> ('<path>')`. With `on_error="skip"` it is recorded instead, and the
rest are read:

```pycon
>>> mixed = pymzlib.readers.read_spectra_many(
...     ["sliced_ethcd.mzML", "no-such-run.mzML", "withZeros.mgf"], on_error="skip"
... )
>>> mixed.read_count, mixed.record_count
(2, 8)
>>> [(f.path.split("/")[-1].split("\\")[-1], f.error.kind) for f in mixed.failed_files]
[('no-such-run.mzML', 'usage')]
>>> [f.record_count for f in mixed.files]
[6, None, 2]

```

A failed file keeps its `files` entry, with `record_count` `None`, not zero, because nothing was
counted, and contributes no rows.

**`out=`: hundreds of files in bounded memory.** With `out`, the table is written to a
tab-separated file **one file at a time, in order**, so memory holds at most `threads` files
however long the list is. The result carries `output` and the per-file facts, and `columns` is
`None`. A batch that stops on an error deletes its partial table rather than leave one that looks
complete.

```python title="Not run: writes a file"
batch = pymzlib.readers.read_records_many(mzid_paths, out="all_matches.tsv", threads=4)
table = pd.read_csv(batch.output.path, sep="\t")
```

Three rules, each an error rather than a surprise:

- **`read_records_many()` needs one record type.** Its columns are the format's own, so a list
  mixing, say, mzIdentML and TopPIC files is refused before anything is parsed, naming the groups.
  The typed views mix formats freely.
- **No `limit` or `offset`.** They window one file; read that file on its own to window it.
- **A path may appear once.** A repeated path is almost always a mistake in how the list was built,
  and reading it twice would double its rows in every total.

## The numbers do not mean the same thing across formats

This is the trap most likely to produce a wrong result, so every read reports it in `caveats`, per
format, citing the mzLib source each one comes from:

```pycon
>>> for c in r.caveats[:2]:
...     print("-", c[:90])                     # doctest: +ELLIPSIS
- ...
- ...

```

`retention_time_unit` gives you the same fact as a **value**, so you can convert in code instead of
hard-coding a table, and `retention_time_in_minutes` converts, or raises rather than guess when the
unit is `'unknown'`:

```pycon
>>> r.retention_time_unit
'minutes'
>>> len(r.retention_time_in_minutes)
2

```

The reason it differs at all is that mzLib's **result-file** readers largely pass each tool's
columns through without normalising them, while its **spectra** readers convert. Where that has
been fixed, it was fixed upstream in mzLib rather than papered over here: MSFragger wrote seconds
until [mzLib#1116](https://github.com/smith-chem-wisc/mzLib/pull/1116) made the reader divide by 60,
and this library's caveat and unit changed with it. Today every quantifiable format reports
`'minutes'`, `read_spectra()` is always minutes, and the one unresolved case is TopFD
`_ms1.feature`, which is `'unknown'`.

**So: identifying a file is safe. Comparing a raw field across formats needs a look at
`retention_time_unit` and `caveats` first.**

!!! info "Why the caveats cite line numbers"
    Each one names the mzLib source it came from, and a test asserts that the cited line still
    mentions what the caveat claims. That is not decoration: citations have gone stale against the
    pinned mzLib before, when an upstream change inserted lines above them. A caveat that reads
    authoritatively and is wrong is worse than no caveat, so the anchoring is checked mechanically.

## Errors

`readers` never returns a sentinel for "unknown". mzLib has no such concept, so a file is
dispatchable or it is a `UsageError`:

| situation | what happens |
|---|---|
| Path does not exist | `UsageError` naming the path |
| Extension mzLib does not recognise | `UsageError` pointing at `formats()` |
| Recognised, but lacks the view a verb needs | `UsageError` naming the views it *does* have, **and pointing at `read_records()`** |
| An option given without a value | `UsageError`, never a silent default |
| A list given to a single-file function | `UsageError` pointing at its `_many` twin |
| A `_many` call's file cannot be read (default `on_error="fail"`) | that file's own error, its message starting `Input <i> ('<path>')` |
| A quantification function on a bridge that predates it (`PYMZLIB_BRIDGE` pointing at an old build) | `UsageError` naming the pyMzLib release it needs, before any process starts |
| `out=` equal to the input path | `UsageError`: a read must not overwrite what it is reading |
| `out=` not ending in `.tsv` (any case), including no extension | `UsageError` before the file is read; no extension is added |
| An answer too large for one JSON document | `UsageError` naming the record count: use `limit`/`offset` or `out=` |

A list given to a single-file function is caught before the bridge starts:

```pycon
>>> pymzlib.readers.read_records(["crux.txt", "crux2.txt"])      # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
pymzlib._bridge.UsageError: ...read_records_many...

```

```python title="Not run: the replay bridge answers only successful calls, and this one fails in mzLib"
try:
    pymzlib.readers.read_features("run_prsm.tsv")
except pymzlib.UsageError as e:
    print(e)
    # 'ToppicPrsm' files do not offer the ms1_features view, so read-features cannot read them:
    # it has no cross-format view at all. Every file type can be read with read-records, which
    # returns that format's own fields.
```

Failures inside mzLib's parallel readers are unwrapped before they reach you, so you see
`MzLibException: Reading profile mode mzmls not supported` rather than
`AggregateException: One or more errors occurred.`

[Every error each verb can raise](../errors.md) is listed on one page.

## Every supported format

`formats()` returns mzLib's own table at run time, so it reflects your installed version rather
than this page's age. Every row is readable with `read_records()`; the views say which typed
functions also apply. This is the list the installed mzLib gives:

```pycon
>>> for f in pymzlib.readers.formats():
...     print(f"{f.file_type:40} {f.extension:28} {', '.join(f.views) or '(none)'}")
Ms1Feature                               _ms1.feature                 ms1_features
Ms2Feature                               _ms2.feature                 (none)
TopFDMzrt                                .mzrt.csv                    (none)
Ms1Tsv_FlashDeconv                       _ms1.tsv                     (none)
Tsv_FlashDeconv                          .tsv                         (none)
Tsv_Dinosaur                             .feature.tsv                 ms1_features
ThermoRaw                                .raw                         spectra
MzML                                     .mzML                        spectra
Mgf                                      .mgf                         spectra
Ms1Align                                 _ms1.msalign                 spectra
Ms2Align                                 _ms2.msalign                 spectra
psmtsv                                   .psmtsv                      quantifiable
osmtsv                                   .osmtsv                      quantifiable
ToppicPrsm                               _prsm.tsv                    (none)
ToppicPrsmSingle                         _prsm_single.tsv             (none)
ToppicProteoform                         _proteoform.tsv              (none)
ToppicProteoformSingle                   _proteoform_single.tsv       (none)
MsFraggerPsm                             psm.tsv                      quantifiable
MsFraggerPeptide                         peptide.tsv                  (none)
MsFraggerProtein                         protein.tsv                  (none)
FlashLFQQuantifiedPeak                   Peaks.tsv                    (none)
MsPathFinderTTargets                     _IcTarget.tsv                spectral_match
MsPathFinderTDecoys                      _IcDecoy.tsv                 spectral_match
MsPathFinderTAllResults                  _IcTDA.tsv                   spectral_match
CruxResult                               .txt                         (none)
ExperimentAnnotation                     experiment_annotation.tsv    (none)
BrukerD                                  .d                           spectra
BrukerTimsTof                            .d                           spectra
CasanovoMzTab                            .mztab                       spectral_match
DiaNnReport                              report.tsv                   quantifiable
Sdrf                                     .sdrf.tsv                    (none)
PytheasResult                            .txt                         (none)
MzIdentML                                .mzid                        spectral_match
MzIdentMLGz                              .mzid.gz                     spectral_match
MetaMorpheusQuantifiedProteinGroups      QuantifiedProteinGroups.tsv  (none)
FlashLFQQuantifiedPeptide                QuantifiedPeptides.tsv       (none)
MetaMorpheusQuantifiedTranscriptGroups   QuantifiedTranscriptGroups.tsv (none)
FlashLFQQuantifiedOligo                  QuantifiedOligos.tsv         (none)

```

Note that **extensions are not unique**: both Bruker types are `.d` (told apart by what the
directory contains), and several formats share `.tsv`, disambiguated by filename suffix and
sometimes by reading the first line. Renaming a file changes how it parses. That is not
hypothetical: mzLib's own Dinosaur test fixture is named `.features.tsv` and cannot be dispatched
until it is renamed to `.feature.tsv`.

`DiaNnReport` is one of two rows where the extension is **not** how dispatch works. mzLib reports
its extension as `report.tsv`, the conventional DIA-NN name, but matches on the header instead: a
file is a DIA-NN report if its first line carries `File.Name`, `Precursor.Id` and
`Stripped.Sequence`. That is deliberate upstream: whoever ran the search routinely renames the
report, and `File.Name` is what separates the long-format report from the `pr_matrix` reports
DIA-NN writes beside it, which carry the other two columns but one column per run. So a renamed
DIA-NN report still reads, and a `report.tsv` that is not one still will not.

`PytheasResult` is the other. Pytheas (RNA oligonucleotide search) writes its match output with no
fixed name, so mzLib looks for a `#theoretical_digest` line in the first five lines of a `.txt`
file. Only a `.txt` without one is read as a `CruxResult`. The records are Pytheas's own match
lines, one per candidate. Charges are negative, and `molecule_location` reads `decoy` on decoy
matches.

The MetaMorpheus quantification tables dispatch on a filename suffix:

- **Protein groups:** any name ending `ProteinGroups.tsv`. That covers
  `AllQuantifiedProteinGroups.tsv` and, since mzLib#1365, the two names MetaMorpheus uses when it
  did not quantify: `AllProteinGroups.tsv` and each file's `<file>_ProteinGroups.tsv`. Those used to
  fail with `Tsv file type not supported`.
- **Peptides:** any name ending `QuantifiedPeptides.tsv`, which is also what FlashLFQ writes.
- **RNA** (mzLib#1388): any name ending `TranscriptGroups.tsv`, and `QuantifiedOligos.tsv`. They are
  their own file types, but mzLib reads them with subclasses of the protein-group and peptide
  readers, so the same three functions read them and the columns keep their protein names:
  `protein_group_name` holds a transcript group, `sequence` an oligonucleotide.

Their per-sample values are read by [their own functions](#quantification-tables-protein-groups-peptides-and-occupancy).

## What is not covered

- **Confidence in `read_results()`.** The quantifiable view exposes no q-value, PEP or score,
  because the mzLib interface it projects carries none; `read_records()` has those columns.
  `read_matches()` carries `q_value` (and mzIdentML's `rank` and `pass_threshold`) where the format
  records one. **Nothing from a typed view is FDR-filtered.**
- **Peptide-level occupancy.** mzLib parses occupancy cells for protein groups only; its
  peptide-table reader exposes none, so `read_occupancy()` reads protein-group tables.
- **Format conversion.** mzLib can write most formats, but the psmtsv family throws
  `NotImplementedException`, so a general read-A-write-B is not offered.
- **The ion-mobility axis.** timsTOF data is read with its mobility dimension collapsed into scans;
  a 1/K0 value is not reported.
- **Bruker off Windows.** Both `.d` types need vendor native libraries that exist only for
  Windows-x64. This is a vendor constraint, not a pyMzLib one.
- **Live objects across calls.** Each call is independent; there is no handle you can hold onto and
  re-query.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the tool
that wrote the file you read:

--8<-- "docs/reference/_generated/cite.readers.md"
