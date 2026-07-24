# Readers

mzLib recognises **29 file types** written by a dozen different search and deconvolution tools —
MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur,
FlashLFQ — and maintains a parser for each. pyMzLib lets you point at a file and ask what it is,
then read it.

```python
import pymzlib

info = pymzlib.readers.identify("psm.tsv")
print(info.file_type, info.views)      # MsFraggerPsm ['quantifiable']

results = pymzlib.readers.read_results("AllPSMs.psmtsv")
print(results.record_count)            # 8
```

## Start with `views`, not with the file type

It would be convenient if mzLib read all 29 formats into one uniform table. **It does not**, and
assuming otherwise is the mistake this module exists to prevent. The formats fall into disjoint
families, and many belong to no family at all:

| view | what it means | which formats |
|---|---|---|
| `quantifiable` | a cross-format record view — sequence, retention time, charge, mass, protein groups. The input [`flashlfq.quantify()`](flashlfq.md) accepts, and the only view `read_results()` can read. | **3 of 29** |
| `ms1_features` | deconvolved MS1 features | TopFD, Dinosaur |
| `spectra` | the file is spectra, not results | `.raw`, `.mzML`, `.mgf`, `.d`, msalign |
| `spectral_match` | records are identifications, but share no *file*-level interface | MsPathFinderT, Casanovo |
| *(none)* | mzLib parses it into a format-specific shape with nothing in common | TopPIC, Crux, and more |

An empty list is a real and common answer, not an error. So the useful first question is not "what
format is this" but "what can I do with it":

```python
info = pymzlib.readers.identify(path)

if info.is_quantifiable:
    results = pymzlib.readers.read_results(path)
else:
    print(f"{info.file_type} supports {info.views or 'no uniform view'}")
```

## Reading records

`read_results()` returns **columnar** data — one array per field, rather than one object per
record. pyMzLib has no third-party dependencies, so it can never hand you a DataFrame; a map of
arrays is the one shape that becomes one in a single call:

```python
import pandas as pd

r = pymzlib.readers.read_results("AllPSMs.psmtsv")
frame = pd.DataFrame(r.columns)            # or pl.DataFrame(r.columns)
```

If you would rather loop, `r.records` gives the same data as one dict per row.

### Nothing is ever silently short

There is **no default row limit**. A result file can carry a million rows, and a library whose
default answer is "here's some of it" is a library that eventually puts a truncated table in a
paper. Ask for a limit and you are told when it bites:

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
r = pymzlib.readers.read_results("AllPSMs.psmtsv", out="records.tsv")
r.output.path, r.output.row_count          # ('C:/.../records.tsv', 84213)
```

The table goes to disk and the envelope carries only a summary. It is **tab-separated**, because
these fields contain commas — MSFragger's mapped proteins are a comma-separated list inside a
single field — and because every mzLib reader and writer uses tabs. Read it with
`pandas.read_csv(path, sep="\t")`, or `csv.reader(f, delimiter="\t")` with no dependencies at all.

!!! warning "`offset` is a window, not a cursor"
    mzLib's readers look lazy and are not — every one of them materializes the whole file into a
    list. So `offset` does not resume where you left off; it re-reads and re-parses the entire file
    and then skips. Paging a large file is quadratic. Use `out=` instead.

## The numbers do not mean the same thing across formats

This is the trap most likely to produce a wrong result, so `read_results()` reports it in
`caveats` — per format, citing the mzLib source each one comes from:

```python
for c in pymzlib.readers.read_results("psm.tsv").caveats:
    print("-", c)
```

```
- retention_time is in SECONDS for this format, not minutes: MSFragger's Retention column is
  passed through unconverted (MsFraggerPsm.cs:48). ...
- is_decoy is always false: mzLib does not read MSFragger decoys (MsFraggerPsm.cs:217). False
  means 'unknown', not 'target'.
- monoisotopic_mass is the THEORETICAL peptide mass (CalculatedPeptideMass), not the observed one ...
- file_name is the full 'Spectrum File' path including its .pep.xml extension ...
```

The reason is that mzLib's result-file readers pass each tool's columns through **without
normalising them**. Retention time is minutes for MetaMorpheus, seconds for MSFragger and TopPIC —
and TopFD changed from seconds to minutes between v1.6.2 and v1.7.0, *within the same file type*.
(mzLib's spectra readers do convert; its result readers do not.)

**So: identifying a file is safe. Comparing a raw field across formats is not.** Read the caveats
before you join two tables or plot them together.

!!! danger "Do not quantify an MSFragger `psm.tsv`"
    Because of the retention-time mismatch above, passing an MSFragger file to
    [`flashlfq.quantify()`](flashlfq.md) returns near-zero intensities — FlashLFQ reads the seconds
    as minutes and looks for each peptide about sixty times too early in the gradient. `identify()`
    will still report the file as `quantifiable`, which is honest about mzLib's *interface* and is
    not a claim about the numbers. This is an mzLib defect affecting every caller and is reported
    upstream; until it is fixed, quantify MetaMorpheus output only.

## Errors

`readers` never returns a sentinel for "unknown" — mzLib has no such concept, so a file is
dispatchable or it is a `UsageError`:

| situation | what happens |
|---|---|
| Path does not exist | `UsageError` naming the path |
| Extension mzLib does not recognise | `UsageError` pointing at `formats()` |
| Recognised, but no `quantifiable` view | `UsageError` from `read_results()` naming the views it *does* have |

```python
try:
    pymzlib.readers.read_results("run_prsm.tsv")
except pymzlib.UsageError as e:
    print(e)
    # Cannot read '...' into the uniform record view. 'ToppicPrsm' files have no cross-format
    # record view at all — mzLib parses them into a format-specific shape only. ...
```

## Every supported format

Generated from mzLib itself — `pymzlib.readers.formats()` returns this same table at runtime, so it
reflects your installed version rather than this page's age.

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

Note that **extensions are not unique**: both Bruker types are `.d` (told apart by what the
directory contains), and several formats share `.tsv`, disambiguated by filename suffix and
sometimes by reading the first line. Renaming a file changes how it parses.

## What is not covered yet

- **Per-format rich records.** The uniform view is the common denominator; the format-specific
  columns each reader parses (matched ions, scores, q-values) are not exposed.
- **Format conversion.** mzLib can write most formats, but the psmtsv family throws
  `NotImplementedException`, so a general read-A-write-B is not offered.
- **Reading spectra as data.** `.raw`/`.mzML` identify as `spectra`, but their scans are
  fine-grained numeric data and need a different transport than one JSON envelope per call.
