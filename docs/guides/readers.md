# Readers

mzLib recognises **31 file types** written by a dozen different search and deconvolution tools —
MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur,
FlashLFQ — and maintains a parser for each. pyMzLib lets you point at a file, ask what it is, and
read it.

```python
import pymzlib

info = pymzlib.readers.identify("psm.tsv")
print(info.file_type, info.views)      # MsFraggerPsm ['quantifiable']

table = pymzlib.readers.read_records("toppic_prsm.tsv")
print(table.record_type, len(table.column_names))    # ToppicPrsm 36
```

**All 31 formats are readable.** What differs between them is not whether you can read them but
what the columns mean — which is the whole subject of this page.

## Five ways to read, and how to choose

There is one universal function and four cross-format views. The choice is a real one, so it is
worth stating plainly before anything else:

| function | reads | columns | use it when |
|---|---|---|---|
| [`read_records()`](#read_records-any-format-its-own-fields) | **all 31** | **this format's own fields**, under mzLib's names | you want *everything* a file has |
| [`read_results()`](#read_results-the-quantifiable-view) | 3 | uniform: sequence, RT, charge, mass, proteins | you are feeding [FlashLFQ](flashlfq.md) or comparing search results |
| [`read_features()`](#read_features-deconvolved-ms1-features) | 2 | uniform: m/z, charge, RT range, intensity | you are working with deconvolved MS1 features |
| [`read_matches()`](#read_matches-identifications) | 4 | uniform: scan, sequences, accession, mods | you are comparing identifications from MsPathFinderT or Casanovo |
| [`read_spectra()`](#read_spectra-scans-and-peaks) | 7 | uniform: scan headers, peaks on request | the file is spectra rather than results |

The rule of thumb:

> **A typed view when you need numbers that mean the same thing across files.
> `read_records()` when you need everything one file has.**

Both matter. A `.psmtsv` read through `read_results()` gives you 10 columns you can safely compare
against an MSFragger file. The same `.psmtsv` through `read_records()` gives you **73**, including
the q-values, PEP and scores the uniform view does not carry — but those column names are
MetaMorpheus's, and no other format has them.

## Start with `views`, not with the file type

It would be convenient if mzLib read all 31 formats into one uniform table. **It does not.** The
formats fall into disjoint families, and thirteen belong to no family at all:

| view | what it means | which formats |
|---|---|---|
| `quantifiable` | a cross-format record view — sequence, retention time, charge, mass, protein groups. What [`flashlfq.quantify()`](flashlfq.md) accepts. | **3**: MetaMorpheus `.psmtsv`/`.osmtsv`, MSFragger `psm.tsv` |
| `ms1_features` | deconvolved MS1 features | **2**: TopFD `_ms1.feature`, Dinosaur |
| `spectral_match` | records are identifications, but share no *file*-level interface | **4**: MsPathFinderT ×3, Casanovo |
| `spectra` | the file is spectra, not results | **7**: `.raw`, `.mzML`, `.mgf`, `.d` ×2, msalign ×2 |
| *(none)* | mzLib parses it into a format-specific shape with nothing in common | **13**: TopPIC ×4, Crux, MSFragger peptide/protein, FlashDeconv, and more |

`views == []` is a real and common answer, not an error — it is the majority answer, in fact. It
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

There is **no default row limit**, on any of the five functions. A result file can carry a million
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
    pyMzLib crosses those as `None` and says so in `caveats`. TopFD files, which do write the
    column, are unaffected. `read_records()` has the file's own summed `intensity` either way.

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

Four formats offer the `spectral_match` view: MsPathFinderT's targets, decoys and combined results,
and Casanovo's `.mztab`. These are the identification formats that share no *file*-level interface,
so `read_results()` cannot reach them.

```python
m = pymzlib.readers.read_matches("results_IcTda.tsv")
m.columns["modifications"][0]        # '12:Oxidation on M'
```

!!! danger "Nothing here is FDR-filtered — and there is no confidence column to filter on"
    mzLib's `ISpectralMatch` carries identity fields only. Every one of these formats records an
    E-value or q-value somewhere; `read_records()` will give you those columns. Filter before you
    report.

Two `is_decoy` traps, both reported in `caveats`:

- **MsPathFinderT** infers decoys from the protein *name* — mzLib reports a decoy when
  `ProteinName` starts with `XXX`. A database whose decoys carry a different prefix reads
  **entirely as targets**.
- **Casanovo** is de novo and writes no target/decoy label at all. mzLib's record leaves the field
  at its default `False` and never assigns it, so `False` would mean *unknown*. pyMzLib crosses it
  as `None` instead — the same rule `read_results()` already applies to MSFragger.

Casanovo also numbers scans by mzTab **index**, not by the instrument's scan number; when Casanovo
was run on an MGF the two are unrelated, so do not join on it.

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

!!! info "Two of the seven need Windows"
    Bruker `.d` and timsTOF `.d` are read through vendor native libraries (`baf2sql`, `timsdata`)
    and are **Windows-x64 only**. Thermo `.raw` uses managed vendor assemblies and works
    everywhere. msalign files hold **deconvolved neutral masses**, not raw m/z — do not
    re-deconvolve them.

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

Note that **extensions are not unique**: both Bruker types are `.d` (told apart by what the
directory contains), and several formats share `.tsv`, disambiguated by filename suffix and
sometimes by reading the first line. Renaming a file changes how it parses — which is not
hypothetical: mzLib's own Dinosaur test fixture is named `.features.tsv` and cannot be dispatched
until it is renamed to `.feature.tsv`.

`DiaNnReport` is the one row where the extension column is **not** how dispatch works. mzLib
reports its extension as `report.tsv`, the conventional DIA-NN name, but matches on the header
instead — a file is a DIA-NN report if its first line carries `File.Name`, `Precursor.Id` and
`Stripped.Sequence`. That is deliberate upstream: whoever ran the search routinely renames the
report, and `File.Name` is what separates the long-format report from the `pr_matrix` reports
DIA-NN writes beside it, which carry the other two columns but one column per run. So a renamed
DIA-NN report still reads, and a `report.tsv` that is not one still will not.

## What is not covered

- **Confidence in the typed views.** `read_results()` and `read_matches()` expose no q-value, PEP or
  score, because the mzLib interfaces they project do not carry one. `read_records()` does — those
  columns exist in every one of these formats. **Nothing from a typed view is FDR-filtered.**
- **Format conversion.** mzLib can write most formats, but the psmtsv family throws
  `NotImplementedException`, so a general read-A-write-B is not offered.
- **The ion-mobility axis.** timsTOF data is read with its mobility dimension collapsed into scans;
  a 1/K0 value is not reported.
- **Bruker off Windows.** Both `.d` types need vendor native libraries that exist only for
  Windows-x64. This is a vendor constraint, not a pyMzLib one.
- **Live objects across calls.** Each call is independent; there is no handle you can hold onto and
  re-query.
