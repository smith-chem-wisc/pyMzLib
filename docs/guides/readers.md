# Readers

**Spectra files are read here, not just search output.** [`read_spectra()`](#read_spectra-scans-and-peaks) reads
**mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign — scan headers always, peaks on
request.

mzLib recognises **36 file types** in all: those instrument and deconvolution formats, plus the
output of a dozen search tools — MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux,
Casanovo, FlashDeconv, Dinosaur, DIA-NN, FlashLFQ, Pytheas, and any mzIdentML writer — and
maintains a parser for each. pyMzLib lets you point at a file, ask what it is, and read it.

```python
import pymzlib

info = pymzlib.readers.identify("psm.tsv")
print(info.file_type, info.views)      # MsFraggerPsm ['quantifiable']

table = pymzlib.readers.read_records("toppic_prsm.tsv")
print(table.record_type, len(table.column_names))    # ToppicPrsm 36
```

**All 36 formats are readable.** What differs between them is not whether you can read them but
what the columns mean — which is the whole subject of this page.

## Ways to read, and how to choose

There is one universal function, four cross-format views, and three functions for mzLib's
quantification tables. The choice is a real one, so it is worth stating plainly before anything
else:

| function | reads | columns | use it when |
|---|---|---|---|
| [`read_records()`](#read_records-any-format-its-own-fields) | **all 36** | **this format's own fields**, under mzLib's names | you want *everything* a file has |
| [`read_results()`](#read_results-the-quantifiable-view) | 4 | uniform: sequence, RT, charge, mass, proteins | you are feeding [FlashLFQ](flashlfq.md) or comparing search results |
| [`read_features()`](#read_features-deconvolved-ms1-features) | 2 | uniform: m/z, charge, RT range, intensity | you are working with deconvolved MS1 features |
| [`read_matches()`](#read_matches-identifications) | 6 | uniform: scan, sequences, accession, mods | you are comparing identifications from MsPathFinderT, Casanovo or mzIdentML |
| [`read_spectra()`](#read_spectra-scans-and-peaks) | 7 | uniform: scan headers, peaks on request | the file is spectra rather than results |
| [`read_protein_groups()`](#quantification-tables-protein-groups-peptides-and-occupancy) | 1 | long: one row per protein group per sample group | you want MetaMorpheus's per-sample protein intensities and spectral counts |
| [`read_quantified_peptides()`](#quantification-tables-protein-groups-peptides-and-occupancy) | 1 | long: one row per peptide per sample | you want FlashLFQ's per-sample peptide intensities |
| [`read_occupancy()`](#quantification-tables-protein-groups-peptides-and-occupancy) | 1 | long: one row per modified site | you want PTM site occupancy from a MetaMorpheus protein-group table |

**Every one of them has a `_many` twin** - `read_spectra_many()`, `read_records_many()`,
`identify_many()` and so on - that reads a list of files into one table in one call. See
[Many files at once](#many-files-at-once).

The rule of thumb:

> **A typed view when you need numbers that mean the same thing across files.
> `read_records()` when you need everything one file has.**

Both matter. A `.psmtsv` read through `read_results()` gives you 10 columns you can safely compare
against an MSFragger file. The same `.psmtsv` through `read_records()` gives you **73**, including
the q-values, PEP and scores the uniform view does not carry — but those column names are
MetaMorpheus's, and no other format has them.

## Start with `views`, not with the file type

It would be convenient if mzLib read all 36 formats into one uniform table. **It does not.** The
formats fall into disjoint families, and seventeen belong to no family at all:

| view | what it means | which formats |
|---|---|---|
| `quantifiable` | a cross-format record view — sequence, retention time, charge, mass, protein groups. What [`flashlfq.quantify()`](flashlfq.md) accepts. | **4**: MetaMorpheus `.psmtsv`/`.osmtsv`, MSFragger `psm.tsv`, DIA-NN `report.tsv` |
| `ms1_features` | deconvolved MS1 features | **2**: TopFD `_ms1.feature`, Dinosaur |
| `spectral_match` | records are identifications, but share no *file*-level interface | **6**: MsPathFinderT ×3, Casanovo, mzIdentML `.mzid`/`.mzid.gz` |
| `spectra` | the file is spectra, not results | **7**: `.raw`, `.mzML`, `.mgf`, `.d` ×2, msalign ×2 |
| *(none)* | mzLib parses it into a format-specific shape with nothing in common | **17**: TopPIC ×4, Crux, Pytheas, MSFragger peptide/protein, FlashDeconv, MetaMorpheus protein groups and peptides, and more |

`views == []` is a real and common answer, not an error — it is the commonest answer, in fact. It
means "mzLib reads this, but there is no uniform projection of it", and `read_records()` is exactly
the function for that case.

```python
info = pymzlib.readers.identify(path)

if info.is_quantifiable:
    table = pymzlib.readers.read_results(path)     # comparable columns
else:
    table = pymzlib.readers.read_records(path)     # this format's own columns
    print(f"{info.file_type}: {len(table.column_names)} native columns")
```

## `read_records()`: any format, its own fields

This is the exhaustive verb. If `identify()` succeeds on a path, this reads it.

```python
t = pymzlib.readers.read_records("toppic_prsm.tsv")

t.file_type        # 'ToppicPrsm'
t.record_type      # 'ToppicPrsm'  - the mzLib class the columns came from
t.views            # []            - no uniform view at all, and still readable
t.column_names     # ['file_name_without_extension', 'prsm_id', ..., 'e_value', ...]
```

Column names are mzLib's own property names converted to `snake_case`, which means they are
**cross-referenceable against the mzLib source**: a column called `e_value` is `ToppicPrsm.EValue`,
and `record_type` tells you which class to look in. Acronyms survive the conversion intact —
`EValue` → `e_value`, `MIScore` → `mi_score`, `FixedPTMs` → `fixed_ptms`.

```python
import pandas as pd
frame = pd.DataFrame(t.columns)
frame[frame.e_value < 1e-10][["base_sequence", "protein_accession", "e_value"]]
```

### Nothing is silently dropped

Some fields cannot become a column. A nested object or a dictionary has no faithful column shape,
and flattening one would mean publishing a schema mzLib does not have. Those fields are **named,
with the reason**, rather than quietly omitted — because a column that simply vanished is
indistinguishable from a field the format does not have:

```python
for field in t.excluded_fields:
    print(field["field"], "—", field["reason"])
```

```
alternative_identifications — a list of composite values has no faithful column shape
```

Some of these exclusions carry the numbers you came for. The per-sample values of mzLib's newer
tables are dictionaries keyed by sample, so `read_records()` names them here and does not project
them:

| file type | excluded field | what it holds | read it with |
|---|---|---|---|
| `MetaMorpheusQuantifiedProteinGroups` | `sample_groups` | per-sample intensity, spectral count and modification occupancy | [`read_protein_groups()`](#quantification-tables-protein-groups-peptides-and-occupancy), and [`read_occupancy()`](#site-occupancy) for the occupancy |
| `FlashLFQQuantifiedPeptide` | `samples` | per-run intensity, detection type and retention time | [`read_quantified_peptides()`](#quantification-tables-protein-groups-peptides-and-occupancy) |
| `MzIdentML`, `MzIdentMLGz` | `scores` | the search engine's own scores, e.g. `MS-GF:SpecEValue` | [`read_matches(scores=True)`](#engine-scores-as-long-rows) |

Their scalar fields (protein group name, gene, organism, q-value, sequence) are columns as usual.
Each `excluded_fields` entry carries a `verb` naming the bridge command that does project it, so
the pointer is data, not only this table:

```python
{e["field"]: e["verb"] for e in t.excluded_fields}
# {'sample_groups': 'readers read-protein-groups'}
```

`failed_fields` is the other half. Several mzLib properties are *computed* and assume a
UniProt-style FASTA header — Crux's and MsPathFinderT's `accession` are both
`protein_id.split("|")[1]` — so on a database with plain headers they raise. Those cells arrive
`None` rather than taking the whole file down with them, but a failure must not look like missing
data:

```python
t = pymzlib.readers.read_records("crux.txt")
t.failed_fields          # ['accession: IndexOutOfRangeException']  (on a non-UniProt database)
```

!!! tip "`-1` is not treated as missing here"
    `read_results()` maps mzLib's documented `-1` "absent" sentinel to `None` for two specific
    interface fields. `read_records()` deliberately does **not** generalise that, because in a
    format's own columns `-1` is frequently a real measurement — a mass difference, a delta, a log
    ratio, TopPIC's `feature_score`. Nulling those would destroy data. Non-finite values (`NaN`,
    infinity) still cross as `None`, since JSON cannot carry them at all.

### Four ways a field can have no value

A `None` in a table can mean four different things, and each is named, so you never have to guess
which:

| where it is named | what it means | what you see in the table |
|---|---|---|
| `absent_fields` | the function defines the field, but **this file's format has no column for it** | `None` in every row |
| `failed_fields` | the field exists, but **reading it threw** on some rows | `None` in those rows |
| `excluded_fields` | the field has **no column shape** (a dictionary, a nested object) | not a column at all |
| none of them | the value is **genuinely missing for that row** - the precursor of an MS1 scan, a blank cell | `None` in that row |

`absent_fields` is the one that protects numbers. When a file lacks an optional column, mzLib does
not say so: it fills in its default, and for a number that default is usually zero. pyMzLib reads
mzLib's own declaration of which columns are optional, checks it against the file's header, and
reports - and blanks - every column the file does not have:

```python
peaks = pymzlib.readers.read_records("QuantifiedPeaks.tsv")          # a current FlashLFQ table
peaks.absent_fields                                                  # ['mbr_score']  (mzLib #1345)

targets = pymzlib.readers.read_matches("run_IcTarget.tsv")           # MSPathFinder, pre-FDR
targets.absent_fields                  # ['q_value', 'rank', 'pass_threshold']
targets.columns["q_value"][:3]         # [None, None, None] - mzLib would have said 0.0: "perfect"
```

The MSPathFinder case is the reason this exists. Its `_IcTarget.tsv` is written *before*
target-decoy analysis and has no `QValue` column, and mzLib types that property as a plain
`double`, so it reads **0** - a perfect q-value - for every match. Filtering on `q_value <= 0.01`
would keep everything. A column in `absent_fields` is always `None`, whatever mzLib filled in.

Every function reports all three lists, so a check can be written once:

```python
def trustworthy(result, column):
    return column not in result.absent_fields and not any(
        f.startswith(column + ":") for f in result.failed_fields)
```

### Values that changed with mzLib 1.0.592

Two formats read differently from the same file than they did before mzLib 1.0.592:

- **MetaMorpheus `.psmtsv`/`.osmtsv`: `pro_forma` is filled on older files** (mzLib #1346).
  MetaMorpheus 1.1.11 and earlier wrote no `ProForma` column, so the column was always `None`.
  mzLib now converts `full_sequence` to ProForma 2.0 when the file has no value, naming a
  modification by its UNIMOD accession where it has one and by name otherwise. A row whose
  `full_sequence` is ambiguous (`|`-joined), or which the converter cannot parse, is still `None`.
  The conversion costs time: reading a 271,551-row, 295 MB `AllPSMs.psmtsv` went from 12.9 s to
  15.6 s (+21%), and a 2,796-row file from 0.35 s to 0.81 s, most of that the one-time load of
  the modification databases. `read_results()` does not carry `pro_forma` and is unaffected.
- **FlashLFQ `QuantifiedPeaks.tsv`: `mbr_score` is `None` where it was `0`** (mzLib #1345). A peak
  with no MBR score (every MSMS peak, and every peak in a file without the column) used to read
  `0.0`, a number that says "scored, and scored nothing". Current FlashLFQ and MetaMorpheus 1.1.11
  peaks tables, which have no `MBR Score` column, used to fail with `HeaderValidationException`;
  they now read, and seven columns join the record: `organism`, `peak_fwhm`, `peak_fwhm_status`,
  `pip_q_value`, `pip_pep`, `decoy_peptide` and `random_rt`.

## `read_results()`: the quantifiable view

The four types offering `quantifiable` — MetaMorpheus `.psmtsv` and `.osmtsv`, MSFragger
`psm.tsv`, DIA-NN `report.tsv` — read into a fixed 10-column shape that is safe to compare between files.

```python
r = pymzlib.readers.read_results("AllPSMs.psmtsv")
frame = pd.DataFrame(r.columns)            # or pl.DataFrame(r.columns)
```

Data comes back **columnar** — one array per field, rather than one object per record. pyMzLib has
no third-party dependencies, so it can never hand you a DataFrame; a map of arrays is the one shape
that becomes one in a single call. If you would rather loop, `r.records` gives the same data as one
dict per row.

### Nothing is ever silently short

There is **no default row limit**, on any of the reading functions. A result file can carry a million
rows, and a library whose default answer is "here's some of it" is a library that eventually puts a
truncated table in a paper. Ask for a limit and you are told when it bites:

```python
r = pymzlib.readers.read_results("AllPSMs.psmtsv", limit=100)
r.returned_count, r.record_count, r.truncated      # (100, 84213, True)
```

`rows_not_read` is the other half of that promise. mzLib drops a malformed row silently — it
collects a warning per unreadable line and the reader discards the list — so a half-corrupt file
reads "successfully" with fewer rows than it contains. pyMzLib counts the difference and reports it:

```python
if r.rows_not_read:
    print(f"warning: {r.rows_not_read} rows in the file did not parse")
```

### Large files: write, don't page

```python
r = pymzlib.readers.read_records("huge_prsm.tsv", out="records.tsv")
r.output.path, r.output.row_count          # ('C:/.../records.tsv', 842130)
```

The table goes to disk and the envelope carries only a summary. It is **tab-separated**, because
these fields contain commas — MSFragger's mapped proteins are a comma-separated list inside a
single field — and because every mzLib reader and writer uses tabs. Read it with
`pandas.read_csv(path, sep="\t")`, or `csv.reader(f, delimiter="\t")` with no dependencies at all.

!!! warning "`offset` is a window, not a cursor"
    mzLib's readers look lazy and are not — every one of them materializes the whole file into a
    list. So `offset` does not resume where you left off; it re-reads and re-parses the entire file
    and then skips. Paging a large file is quadratic. Use `out=` instead.

## `read_features()`: deconvolved MS1 features

Two formats offer the `ms1_features` view: TopFD/FLASHDeconv `_ms1.feature` and Dinosaur
`.feature.tsv`. The columns are `mz`, `charge`, `retention_time_start`, `retention_time_end`,
`intensity`, `number_of_isotopes`.

```python
f = pymzlib.readers.read_features("sample_ms1.feature")
f.record_count, f.retention_time_unit       # (25, 'unknown')
```

!!! warning "One row is not one line of the file — for `_ms1.feature`"
    An `_ms1.feature` row is a deconvolved **neutral mass spanning a charge range**, and mzLib
    expands it into one single-charge feature per charge in `[ChargeStateMin, ChargeStateMax]`. A
    hundred-feature file can read as a thousand rows. Dinosaur is one-for-one. Either way
    `read_records()` gives you the file's own rows.

`intensity` is the **apex** intensity, not the sum over the feature — both formats carry a summed
intensity column too, and `read_records()` has it.

!!! warning "`intensity` is `None` for every FLASHDeconv `_ms1.feature`"
    mzLib takes the per-charge intensity from `Apex_intensity`, which is an *optional* column that
    the FLASHDeconv/OpenMS `_ms1.feature` layout does not have — and substitutes **zero** when it is
    absent. A whole column of zeros is indistinguishable from real measurements of nothing, so
    pyMzLib crosses those as `None`, names `intensity` in `absent_fields`, and says why in
    `caveats`. TopFD files, which do write the column, are unaffected. `read_records()` has the
    file's own summed `intensity` either way.

### `retention_time_unit` is `'unknown'` for `_ms1.feature`, and that is the honest answer

TopFD wrote retention times in **seconds** through v1.6.2 and in **minutes** from v1.7.0 — *within
the same file type*, with nothing in the file to tell you which. mzLib does not normalise either,
and its own deconvolution code resorts to a heuristic (divide everything by 60 if the largest end
time exceeds 500). pyMzLib will not launder a guess into a stated fact:

```python
f.retention_time_start_in_minutes
# UsageError: Cannot convert retention time for 'Ms1Feature': mzLib gives no basis to say what
# unit it is in. TopFD changed from seconds to minutes at v1.7.0 without changing the file type...
```

Dinosaur reports `'minutes'` and converts without complaint.

## `read_matches()`: identifications

Six formats offer the `spectral_match` view: MsPathFinderT's targets, decoys and combined results,
Casanovo's `.mztab`, and mzIdentML (`.mzid`, and `.mzid.gz` read without unpacking it). These are
the identification formats that share no *file*-level interface, so `read_results()` cannot reach
them.

```python
m = pymzlib.readers.read_matches("results_IcTda.tsv")
m.columns["modifications"][0]        # '12:Oxidation on M'
```

Beside the identity fields, the view carries the three confidence fields a format may record:
`q_value` (mzIdentML's PSM-level q-value, mzLib #1306; MSPathFinder's `QValue`), and mzIdentML's
`rank` and `pass_threshold`. A format without one names it in `absent_fields`, and the column is
`None` throughout.

!!! danger "Nothing here is FDR-filtered"
    Every row the file holds is a row here. Filter on `q_value` - where `absent_fields` does not
    name it - and, for mzIdentML, on `rank == 1` and `pass_threshold`, before you count or report
    anything. The engines' other scores are in `read_records()`, or as long rows with
    [`scores=True`](#engine-scores-as-long-rows).

Three `is_decoy` traps, all reported in `caveats`:

- **MsPathFinderT** infers decoys from the protein *name* — mzLib reports a decoy when
  `ProteinName` starts with `XXX`. A database whose decoys carry a different prefix reads
  **entirely as targets**.
- **Casanovo** is de novo and writes no target/decoy label at all. mzLib's record leaves the field
  at its default `False` and never assigns it, so `False` would mean *unknown*. pyMzLib crosses it
  as `None` instead — the same rule `read_results()` already applies to MSFragger.
- **mzIdentML**'s `isDecoy` attribute is optional and defaults to false, so a writer that omits it
  reads as a target. `read_matches()` crosses `None`; `read_records()` carries mzLib's own boolean
  for when you know your writer sets it.

Casanovo also numbers scans by mzTab **index**, not by the instrument's scan number; when Casanovo
was run on an MGF the two are unrelated, so do not join on it.

mzIdentML has three more things to know, also in `caveats`:

- **Every identification item is a row**, not only the ones the submitter accepted. Lower-ranked
  candidates and items that failed the threshold are included; filter on `rank` and
  `pass_threshold`.
- **Some items are skipped, not read**: crosslinks, modifications without a resolvable UNIMOD
  accession, substitutions, and two modifications on one residue (mzLib #1313). They are reported:
  `skipped_count` says how many, and `skipped` lists each with its reason, so `record_count +
  skipped_count` is the number of items in the file.

  ```python
  m = pymzlib.readers.read_matches("xlink_search.mzid")
  m.record_count, m.skipped_count      # (0, 16) - a crosslink search: nothing is a linear match
  m.skipped[0].reason                  # 'crosslink identification'
  ```
- **Scan numbers come from the nativeID.** `scan=N` gives N, but `index=N` (peak-list input) is a
  zero-based position and gives N + 1, which is not an instrument scan number.

### Engine scores as long rows

mzIdentML carries each search engine's own scores - `MS-GF:SpecEValue`, `Mascot:score`,
`Scaffold:Peptide Probability` - and no two engines share a name. There is no fixed set of columns
to put them in, so `scores=True` makes the table **long** instead: one row per match *and* score,
with `match_index` (the match's position in its file), `score_name` and `score_value`.

```python
s = pymzlib.readers.read_matches("search.mzid.gz", scores=True)
s.returned_count, s.row_count     # (12, 84): 12 matches, 7 scores each
frame = pd.DataFrame(s.columns)
scores = frame.pivot(index="match_index", columns="score_name", values="score_value")
scores["MS-GF:SpecEValue"].lt(1e-10).sum()
```

`returned_count` still counts matches - the unit `limit` and `offset` count in - and `row_count`
counts rows. A format with no engine scores keeps one row per match and names `score_name` and
`score_value` in `absent_fields`.

## `read_spectra()`: scans and peaks

Seven formats offer the `spectra` view. Retention times here **are** in minutes for every one of
them — mzLib's spectra readers convert at the boundary, unlike its result-file readers.

```python
s = pymzlib.readers.read_spectra("run.mzML", ms_order=2, limit=5)
s.scan_count, s.record_count           # (14238, 11902)  - file total, then filtered
s.columns["selected_ion_mz"][:3]       # [447.7391, 551.2903, 638.8215]
```

`ms_order` filters **before** the offset/limit window, so `ms_order=2, limit=10` means the first ten
MS2 scans rather than the MS2 scans among the first ten. `scan_count` always reports the file's real
total, so a filter that matched nothing can never look like an empty file.

### Peaks are opt-in

A scan header is tens of bytes; its peak list is thousands, and a mid-size mzML holds tens of
thousands of scans. Returning peaks by default would make the ordinary "what is in this file?" call
serialise hundreds of megabytes.

```python
s = pymzlib.readers.read_spectra("run.mzML", peaks=True, ms_order=1, limit=1)
mz  = s.columns["mz"][0]           # one array per scan
ion = s.columns["intensity"][0]
```

Without `peaks=True`, `peak_count` still tells you how many peaks each scan has.

### What the file says about the run

Every spectra read reports what the file records about where and when it was acquired, from mzLib's
`SourceFile` (mzLib #1349):

```python
s = pymzlib.readers.read_spectra("run.raw", limit=0)
s.source.instrument_model               # 'Orbitrap Fusion Lumos'
s.source.instrument_serial_number       # 'EXRFSN20410'
s.source.acquisition_start_time         # '2023-10-25T10:40:10.188556'
s.source.acquisition_start_time_is_utc  # False - the instrument PC's clock
s.source.acquired_at                    # datetime(2023, 10, 25, 10, 40, 10, 188556), naive
```

These are what a batch or instrument confound is built from - which unit of a model, and in what
order the runs were acquired. Three things to know:

- **Match instruments on the accession, not the name.** mzML carries both the model's name and its
  PSI-MS accession (`instrument_model_accession`, e.g. `MS:1002416`); Thermo `.raw` records only
  the name, so the accession is `None` there rather than looked up.
- **A time is UTC only when it says so.** `acquisition_start_time` ends in `Z`, and
  `acquisition_start_time_is_utc` is true, only when the file fixed the instant. A Thermo `.raw`
  records the acquisition computer's local clock with no offset; `acquired_at` then returns a
  *naive* `datetime` rather than inventing a time zone. A `.raw` and ProteoWizard's mzML of it can
  disagree by the site's UTC offset, because ProteoWizard converts using the *converting*
  machine's time zone and writes `Z`.
- **MGF and msalign record none of it**: every field is `None`.

!!! info "Two of the seven need Windows"
    Bruker `.d` and timsTOF `.d` are read through vendor native libraries (`baf2sql`, `timsdata`)
    and are **Windows-x64 only**. Thermo `.raw` uses managed vendor assemblies and works
    everywhere. msalign files hold **deconvolved neutral masses**, not raw m/z — do not
    re-deconvolve them.

## Quantification tables: protein groups, peptides and occupancy

MetaMorpheus and FlashLFQ write their quantification as wide tables, one column per sample:
`Intensity_QE-002106_GM1_a-calib`, `SpectralCount_QE-002106_GM1_a-calib`, ... mzLib 1.0.592 reads
them (mzLib #1347), and pyMzLib returns them **long** - one row per record per sample, with the
sample in a column:

| function | reads | one row per | per-sample columns |
|---|---|---|---|
| `read_protein_groups()` | MetaMorpheus `AllQuantifiedProteinGroups.tsv` | protein group x sample group | `spectral_count`, `intensity` |
| `read_quantified_peptides()` | FlashLFQ `QuantifiedPeptides.tsv`, MetaMorpheus `AllQuantifiedPeptides.tsv` | peptide x sample | `intensity`, `detection_type`, `retention_time` (IsoTracker only) |
| `read_occupancy()` | MetaMorpheus `AllQuantifiedProteinGroups.tsv` | group x sample group x basis x modified site | `fraction`, `numerator`, `denominator` |

**Why long.** A wide table's column names would be made up from your sample labels, so they would
differ in every experiment and could not be documented, checked or shared across files - and two
searches with different samples could not be read into one table. A long table has the same
columns every time, it is what pandas, polars and every plotting library want, and the wide matrix
is one `pivot` away:

```python
import pandas as pd

g = pymzlib.readers.read_protein_groups("AllQuantifiedProteinGroups.tsv")
g.record_count, g.returned_count, g.row_count      # (6, 6, 108): 6 groups x 18 sample groups
df = pd.DataFrame(g.columns)

# The table is UNFILTERED, as MetaMorpheus writes it: filter before anything else.
confident = df[(df.q_value <= 0.01) & (df.decoy_contaminant_target == "T")]

matrix = confident.pivot(index="protein_group_name", columns="sample_label", values="intensity")
```

`record_count` and `returned_count` count **groups** (or peptides) - the unit `limit` and `offset`
count in - and `row_count` counts rows. The group's other fields - coverage, masses, member
counts - are in `read_records()` on the same file; join on `protein_group_name`.

Four things these tables will not do for you:

- **Filter.** Decoys, contaminants and groups above 1% FDR are all rows; see above.
- **Tell a blank from a zero - unless you look.** A protein-group intensity that was not measured
  is a blank cell and arrives `None`. FlashLFQ's peptide table is different: it writes a literal
  **0** for a peptide it did not quantify, and mzLib keeps it. `detection_type` (`MSMS`, `MBR`,
  `NotDetected`, ...) is what separates a measurement from a missing value; filter on it before a
  mean or a log.
- **Know your design.** `sample_label` is the header label verbatim. Condition, replicate and
  channel cannot be recovered from it - map it yourself, from an [SDRF](sdrf.md) or a sample sheet.
- **Name a leading protein.** `protein_group_name` lists the members sorted by accession, so the
  first is only the one that sorts first.

### Site occupancy

MetaMorpheus writes two occupancy cells per group per sample group - one from PSM counts, one from
intensities - each a list of modified sites encoded as text. `read_occupancy()` returns one row per
site, with `basis` saying which cell it came from:

```python
o = pymzlib.readers.read_occupancy("AllQuantifiedProteinGroups.tsv")
sites = pd.DataFrame(o.columns)

by_count = sites[sites.basis == "count"].assign(occupancy=lambda d: d.numerator / d.denominator)
by_intensity = sites[sites.basis == "intensity"]            # use `fraction` here
```

The two bases round differently, so trust the ratio of the counts on `count` rows and `fraction`
on `intensity` rows. `position` counts 0 as the protein N-terminus and residues from 1. A cell the
writer cut short keeps its complete sites with `cell_is_truncated` set, and `truncated_cell_count`
counts every cut cell - including those with no complete site left to show - so a short table is
never mistaken for a complete one. A cell that is not an occupancy cell at all is refused, named
in `failed_fields`, and gives no rows.

## Many files at once

Every reading function has a `_many` twin that takes a list of paths and returns **one** long
table: `read_spectra_many()`, `read_records_many()`, `read_results_many()`,
`read_features_many()`, `read_matches_many()`, `read_protein_groups_many()`,
`read_quantified_peptides_many()`, `read_occupancy_many()`, and `identify_many()`.

```python
from pathlib import Path
import pandas as pd

runs = sorted(Path("raw").glob("*.raw"))
batch = pymzlib.readers.read_spectra_many(runs, ms_order=1, threads=4)

scans = pd.DataFrame(batch.columns)          # source_index, source_path, then read_spectra's columns
runs_meta = pd.DataFrame([
    {"source_index": i, "instrument": f.source.instrument_serial_number, "started": f.source.acquired_at}
    for i, f in enumerate(batch.files)
])
scans = scans.merge(runs_meta, on="source_index")
```

**When to use it.** Whenever you would otherwise write a loop over files. The list goes to one
bridge process, so the .NET start-up (about 120 ms) and the JIT warm-up are paid once rather than
once per file, and mzLib reads `threads` files at a time. For two hundred small files that is the
difference between most of a minute of start-up and none. pyMzLib has no thread or process pool of
its own, deliberately: the bridge is the one place that can count every thread in use.

**The table.** Its first two columns are `source_index` - the file's position in your list - and
`source_path`. Rows are grouped by file in list order, and within a file in the file's own order.
`batch.files[i]` is a `FileReport` with everything a single-file read would have told you about
file `i`: its `file_type`, `record_count`, `caveats`, `absent_fields`, `source`, `skipped`, and so
on. `batch.record_count` and `batch.row_count` are summed over the files.

**`threads`.** The result is **byte-identical at any `threads`** - that is tested - so it only
trades memory for speed, never changes an answer. It defaults to 1 because every mzLib reader
holds a whole file in memory while it works, so `threads=8` can mean eight whole files at once.
Some guidance:

- Many small files (search results, feature tables, MGFs): `threads` near your core count, or
  `-1` for one per core.
- Large spectra files: mzLib's mzML and Thermo readers already parallelise *inside* one file, so
  extra files in flight add memory faster than speed. Start at 2 and watch memory.
- Hundreds of files, or peaks: add `out=` - see below.

**`on_error`.** By default the first file that cannot be read stops the batch and raises, with a
message starting `Input <i> ('<path>')`. With `on_error="skip"` it is recorded instead, and the
rest are read:

```python
batch = pymzlib.readers.read_matches_many(paths, on_error="skip")
for f in batch.failed_files:
    print(f.path, f.error.kind, f.error.message)   # kind: 'usage' (missing, not this view) or 'correctness'
```

A failed file keeps its `files` entry - with `record_count` `None`, not zero, because nothing was
counted - and contributes no rows.

**`out=`: hundreds of files in bounded memory.** With `out`, the table is written to a
tab-separated file **one file at a time, in order**, so memory holds at most `threads` files
however long the list is. The result carries `output` and the per-file facts, and `columns` is
`None`. A batch that stops on an error deletes its partial table rather than leave one that looks
complete.

```python
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

This is the trap most likely to produce a wrong result, so every read reports it in `caveats` —
per format, citing the mzLib source each one comes from:

```python
for c in pymzlib.readers.read_results("psm.tsv").caveats:
    print("-", c)
```

```
- is_decoy is null for this format: MSFragger's psm.tsv carries no target/decoy column, so mzLib
  cannot report decoy status (MsFraggerPsm.cs:231). Null means 'unknown', not 'target'.
- monoisotopic_mass is the THEORETICAL peptide mass (MsFraggerPsm.cs:233, CalculatedPeptideMass),
  not the observed precursor mass. ...
- file_name is the full 'Spectrum File' path including its .pep.xml extension ...
```

`retention_time_unit` gives you the same fact as a **value**, so you can convert programmatically
instead of hard-coding a table:

```python
r = pymzlib.readers.read_results("psm.tsv")
r.retention_time_unit            # 'minutes'
r.retention_time_in_minutes      # converted; raises rather than guess if the unit is 'unknown'
```

The reason it differs at all is that mzLib's **result-file** readers largely pass each tool's
columns through without normalising them, while its **spectra** readers convert. Where that has
been fixed, it was fixed upstream in mzLib rather than papered over here: MSFragger wrote seconds
until [mzLib #1116](https://github.com/smith-chem-wisc/mzLib/pull/1116) made the reader divide by
60, and this library's caveat and unit changed with it. Today all four quantifiable formats report
`'minutes'`, `read_spectra()` is always minutes, and the one genuinely unresolved case is TopFD
`_ms1.feature`, which is `'unknown'`.

**So: identifying a file is safe. Comparing a raw field across formats needs a look at
`retention_time_unit` and `caveats` first.**

!!! info "Why the caveats cite line numbers"
    Each one names the mzLib source it came from, and a test asserts that the cited line still
    mentions what the caveat claims. That is not decoration: two citations had already gone stale
    against the pinned mzLib because #1116 inserted fourteen lines above them. A caveat that reads
    authoritatively and is wrong is worse than no caveat, so the anchoring is checked mechanically.

## Errors

`readers` never returns a sentinel for "unknown" — mzLib has no such concept, so a file is
dispatchable or it is a `UsageError`:

| situation | what happens |
|---|---|
| Path does not exist | `UsageError` naming the path |
| Extension mzLib does not recognise | `UsageError` pointing at `formats()` |
| Recognised, but lacks the view a verb needs | `UsageError` naming the views it *does* have, **and pointing at `read_records()`** |
| An option given without a value | `UsageError` — never a silent default |
| A list given to a single-file function | `UsageError` pointing at its `_many` twin |
| A `_many` call's file cannot be read (default `on_error="fail"`) | that file's own error, its message starting `Input <i> ('<path>')` |
| A quantification function on a bridge that predates it (`PYMZLIB_BRIDGE` pointing at an old build) | `UsageError` naming the pyMzLib release it needs, before any process starts |
| `out=` equal to the input path | `UsageError` — a read must not overwrite what it is reading |

```python
try:
    pymzlib.readers.read_features("run_prsm.tsv")
except pymzlib.UsageError as e:
    print(e)
    # 'ToppicPrsm' files do not offer the ms1_features view, so read-features cannot read them —
    # it has no cross-format view at all. Every file type can be read with read-records, which
    # returns that format's own fields.
```

Failures inside mzLib's parallel readers are unwrapped before they reach you, so you see
`MzLibException: Reading profile mode mzmls not supported` rather than
`AggregateException: One or more errors occurred.`

## Every supported format

Generated from mzLib itself — `pymzlib.readers.formats()` returns this same table at runtime, so it
reflects your installed version rather than this page's age. Every row is readable with
`read_records()`; the `views` column says which typed functions also apply.

| file type | extension | views |
|---|---|---|
| `Ms1Feature` | `_ms1.feature` | `ms1_features` |
| `Ms2Feature` | `_ms2.feature` | (none) |
| `TopFDMzrt` | `.mzrt.csv` | (none) |
| `Ms1Tsv_FlashDeconv` | `_ms1.tsv` | (none) |
| `Tsv_FlashDeconv` | `.tsv` | (none) |
| `Tsv_Dinosaur` | `.feature.tsv` | `ms1_features` |
| `ThermoRaw` | `.raw` | `spectra` |
| `MzML` | `.mzML` | `spectra` |
| `Mgf` | `.mgf` | `spectra` |
| `Ms1Align` | `_ms1.msalign` | `spectra` |
| `Ms2Align` | `_ms2.msalign` | `spectra` |
| `psmtsv` | `.psmtsv` | `quantifiable` |
| `osmtsv` | `.osmtsv` | `quantifiable` |
| `ToppicPrsm` | `_prsm.tsv` | (none) |
| `ToppicPrsmSingle` | `_prsm_single.tsv` | (none) |
| `ToppicProteoform` | `_proteoform.tsv` | (none) |
| `ToppicProteoformSingle` | `_proteoform_single.tsv` | (none) |
| `MsFraggerPsm` | `psm.tsv` | `quantifiable` |
| `MsFraggerPeptide` | `peptide.tsv` | (none) |
| `MsFraggerProtein` | `protein.tsv` | (none) |
| `FlashLFQQuantifiedPeak` | `Peaks.tsv` | (none) |
| `MsPathFinderTTargets` | `_IcTarget.tsv` | `spectral_match` |
| `MsPathFinderTDecoys` | `_IcDecoy.tsv` | `spectral_match` |
| `MsPathFinderTAllResults` | `_IcTDA.tsv` | `spectral_match` |
| `CruxResult` | `.txt` | (none) |
| `ExperimentAnnotation` | `experiment_annotation.tsv` | (none) |
| `BrukerD` | `.d` | `spectra` |
| `BrukerTimsTof` | `.d` | `spectra` |
| `CasanovoMzTab` | `.mztab` | `spectral_match` |
| `DiaNnReport` | `report.tsv` | `quantifiable` |
| `Sdrf` | `.sdrf.tsv` | (none) |
| `PytheasResult` | `.txt` | (none) |
| `MzIdentML` | `.mzid` | `spectral_match` |
| `MzIdentMLGz` | `.mzid.gz` | `spectral_match` |
| `MetaMorpheusQuantifiedProteinGroups` | `QuantifiedProteinGroups.tsv` | (none) |
| `FlashLFQQuantifiedPeptide` | `QuantifiedPeptides.tsv` | (none) |

Note that **extensions are not unique**: both Bruker types are `.d` (told apart by what the
directory contains), and several formats share `.tsv`, disambiguated by filename suffix and
sometimes by reading the first line. Renaming a file changes how it parses — which is not
hypothetical: mzLib's own Dinosaur test fixture is named `.features.tsv` and cannot be dispatched
until it is renamed to `.feature.tsv`.

`DiaNnReport` is one of two rows where the extension column is **not** how dispatch works. mzLib
reports its extension as `report.tsv`, the conventional DIA-NN name, but matches on the header
instead — a file is a DIA-NN report if its first line carries `File.Name`, `Precursor.Id` and
`Stripped.Sequence`. That is deliberate upstream: whoever ran the search routinely renames the
report, and `File.Name` is what separates the long-format report from the `pr_matrix` reports
DIA-NN writes beside it, which carry the other two columns but one column per run. So a renamed
DIA-NN report still reads, and a `report.tsv` that is not one still will not.

`PytheasResult` is the other. Pytheas (RNA oligonucleotide search) writes its match output with no
fixed name, so mzLib looks for a `#theoretical_digest` line in the first five lines of a `.txt`
file. Only a `.txt` without one is read as a `CruxResult`. The records are
Pytheas's own match lines, one per candidate. Charges are negative, and `molecule_location` reads
`decoy` on decoy matches.

The two MetaMorpheus quantification tables dispatch on a filename suffix:
`AllQuantifiedProteinGroups.tsv` (or any name ending `QuantifiedProteinGroups.tsv`) and
`AllQuantifiedPeptides.tsv` (any name ending `QuantifiedPeptides.tsv`, which is also what FlashLFQ
writes). Their per-sample values are read by [their own functions](#quantification-tables-protein-groups-peptides-and-occupancy).

## What is not covered

- **Confidence in `read_results()`.** The quantifiable view exposes no q-value, PEP or score,
  because the mzLib interface it projects carries none; `read_records()` has those columns.
  `read_matches()` carries `q_value` (and mzIdentML's `rank` and `pass_threshold`) where the
  format records one. **Nothing from a typed view is FDR-filtered.**
- **Peptide-level occupancy.** mzLib 1.0.592 parses occupancy cells for protein groups only; its
  peptide-table reader exposes none, so `read_occupancy()` reads protein-group tables.
- **Format conversion.** mzLib can write most formats, but the psmtsv family throws
  `NotImplementedException`, so a general read-A-write-B is not offered.
- **The ion-mobility axis.** timsTOF data is read with its mobility dimension collapsed into scans;
  a 1/K0 value is not reported.
- **Bruker off Windows.** Both `.d` types need vendor native libraries that exist only for
  Windows-x64. This is a vendor constraint, not a pyMzLib one.
- **Live objects across calls.** Each call is independent; there is no handle you can hold onto and
  re-query.
