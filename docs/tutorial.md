# Your first hour with pyMzLib

A complete walkthrough, from an empty folder to a table you can plot. Everything here uses **public
data**, runs on any machine, and takes about ten minutes of typing and fifty of waiting for
downloads if you do the optional parts.

**No proteomics background is assumed.** Terms are defined where they first appear, in a box like
this one:

!!! abstract "Jargon: *proteomics*"
    Measuring the proteins in a biological sample. A mass spectrometer weighs molecules very
    precisely; software then works out which proteins those weights came from. The files in this
    tutorial are the inputs and outputs of that process.

---

## 0. Install

```bash
pip install pymzlib
```

That is genuinely all. There is no .NET to install and no other Python package to reconcile —
pyMzLib carries its own machinery and declares **zero** third-party dependencies, so it cannot
break anything already in your environment.

Check it worked:

```python
import pymzlib
print(pymzlib.bridge_version())
```

```
{'bridge': '0.1.0', 'protocol': 1, 'runtime': '8.0.20'}
```

If that printed, you are done installing. (`pandas` shows up later and is *optional* — every
example works without it.)

---

## 1. Find some data

Public proteomics data lives in the [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/), run by
the European Bioinformatics Institute. Every dataset has an accession like `PXD000001`.

!!! abstract "Jargon: *accession*"
    A permanent ID for a dataset, like a DOI for a paper. `PXD000001` is the very first PRIDE
    submission and is deliberately small, which makes it the standard thing to try first.

```python
files = pymzlib.pride.list_files("PXD000001")
for f in files:
    print(f"{f.file_name:55} {f.size_mb:8.1f} MB   {f.category}")
```

```
PRIDE_Exp_Complete_Ac_22134.xml.gz                        60.4 MB   OTHER
PRIDE_Exp_mzData_Ac_22134.xml.gz                          52.8 MB   PEAK
F063721.dat-mztab.txt                                      2.5 MB   SEARCH
...
```

### The first real lesson

Now ask the same question a different way:

```python
rest = pymzlib.pride.list_files("PXD000001")       # PRIDE's REST API
ftp  = pymzlib.pride.list_ftp_files("PXD000001")   # walking the FTP directory
print(len(rest), len(ftp))
```

```
8 13
```

**The two disagree, and the shorter one is the API's.** PRIDE's REST manifest omits five of this
project's thirteen files, including the two largest. Neither number is a bug in pyMzLib; the point
is that pyMzLib shows you both rather than picking one and letting you assume it is complete.

That is the single idea running through this whole library, so it is worth stating once:

> **Where the underlying data is incomplete, ambiguous, or means different things in different
> files, pyMzLib tells you — rather than papering over it with a plausible number.**

You will meet it again three more times before the end of this page.

---

## 2. Look at a file without downloading it

Say someone hands you a file and you do not know what it is.

```python
info = pymzlib.readers.identify("psm.tsv")
print(info.file_type, info.views)
```

```
MsFraggerPsm ['quantifiable']
```

mzLib recognises **29 file formats** from a dozen different tools. `identify()` works out which one
this is without parsing the contents, so it is instant even on a million-row file.

!!! abstract "Jargon: *PSM*"
    Peptide-Spectrum Match. One row saying "this measured spectrum came from this peptide". A
    search engine (MSFragger, MetaMorpheus, TopPIC…) produces millions of them, and a `.psmtsv` or
    `psm.tsv` file is the table of results.

`views` is the interesting part. It says what **shapes** this file can be read into:

```python
for fmt in pymzlib.readers.formats():
    print(f"{fmt.file_type:26} {str(fmt.views) or '(none)'}")
```

Most formats print `[]`. **That is not an error** — 13 of the 29 have no shape in common with any
other format, because the tools that wrote them simply do not agree on what a result row is.

---

## 3. Read it — two different ways, on purpose

### Everything the file has

```python
t = pymzlib.readers.read_records("toppic_prsm.tsv")
print(t.record_type, len(t.column_names))
```

```
ToppicPrsm 36
```

`read_records()` reads **any** of the 29 formats, giving you that format's own fields under its own
names. Thirty-six columns for TopPIC, twenty-three for Crux, five for an experiment annotation.

```python
import pandas as pd                     # optional, but nice here
frame = pd.DataFrame(t.columns)
frame[frame.e_value < 1e-10][["base_sequence", "protein_accession", "e_value"]].head()
```

### Only what is comparable between files

```python
r = pymzlib.readers.read_results("AllPSMs.psmtsv")
print(r.column_names)
```

```
['file_name', 'base_sequence', 'full_sequence', 'retention_time', 'charge_state',
 'monoisotopic_mass', 'is_decoy', 'protein_accessions', 'gene_name', 'organism']
```

Ten columns, and they mean the same thing for every format that offers them — so you can stack two
different search engines' results in one table.

**The rule of thumb:**

> A typed view (`read_results`, `read_features`, `read_matches`, `read_spectra`) when you need
> numbers that mean the same thing across files. `read_records` when you need everything one file
> has.

### The second lesson

```python
for c in r.caveats:
    print("-", c)
```

Every read comes back with `caveats`: what these particular numbers cannot be trusted to mean, each
citing the line of mzLib source it came from. For an MSFragger file you will be told that
`is_decoy` is `None` rather than `False` —

!!! abstract "Jargon: *decoy*"
    A deliberately wrong protein sequence added to the search, used to estimate how many "matches"
    are just luck. Filtering them out is how false-discovery rate is controlled.

— because MSFragger's `psm.tsv` has no target/decoy column at all. `False` would have been a
*fabricated* answer that you could group by. `None` means *unknown*, and that is the honest one.

---

## 4. Spectra, and the third lesson

```python
s = pymzlib.readers.read_spectra("run.mzML", ms_order=1, limit=5)
print(s.scan_count, s.record_count)
```

```
14238 2336
```

!!! abstract "Jargon: *MS1 and MS2*"
    An MS1 scan weighs whole peptides as they come off the column. The instrument then picks
    interesting ones, smashes them, and weighs the pieces — that is an MS2 scan. `ms_order=1` and
    `ms_order=2` select between them.

Two numbers, deliberately. `scan_count` is the file's real total; `record_count` is what survived
your filter. A library that reported only the second would let a filter that matched **nothing**
look exactly like an empty file.

Peaks — the actual m/z and intensity arrays — are **opt-in**:

```python
s = pymzlib.readers.read_spectra("run.mzML", ms_order=1, limit=1, peaks=True)
mz, intensity = s.columns["mz"][0], s.columns["intensity"][0]
print(len(mz), "peaks in scan", s.columns["one_based_scan_number"][0])
```

They are off by default because a scan header is tens of bytes and its peak list is thousands; a
mid-size file would otherwise hand you hundreds of megabytes for the question "what is in here?".

---

## 5. Predict something

You can also ask a machine-learning model what a peptide *should* look like.

```python
r = pymzlib.prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK", "ELVISLIVESK"])
print(r.columns["retention_time"])
print(r.retention_time_unit)
```

```
[5.5165, 129.7723]
indexed_retention_time
```

!!! abstract "Jargon: *retention time*"
    Peptides are separated by a chromatography column before they reach the mass spectrometer.
    Retention time is when a given peptide comes off — usually stated in minutes.

### The fourth lesson, and the one most likely to bite you

Look at that unit again. It says `indexed_retention_time`, **not** `minutes`.

Most of these models are trained on **iRT** — a dimensionless scale anchored to standard peptides,
not a clock. A value of 129.77 looks exactly like a plausible 130-minute gradient, which is why
this is the trap: nothing about the number itself reveals the mistake. You have to fit the
iRT-to-minutes line on peptides you have measured before you can compare it to a real run.

One model in the family, `Chronologer_RT`, does return minutes. `retention_time_unit` is how you
tell — which is why it is a value you can branch on rather than a sentence in a docstring.

### Check before you send

```python
for m in pymzlib.prediction.models("fragment_intensity")[:3]:
    print(m.model, m.max_peptide_length, m.collision_energy.requirement)
```

```
AlphaPeptDeep_ms2_generic 500 any_value_required
Altimeter_2024_intensities 40 one_of
Prosit_2020_intensity_CID 30 not_applicable
```

`max_peptide_length` is **30** for most models, which excludes a lot of real peptides — they come
back as `None` with a warning rather than an error, so check `failed_row_count`. And
`collision_energy.requirement` is a three-way answer, because "this model has no such input" and
"this model requires one, any value" are genuinely different things that mzLib's raw field
conflates.

!!! warning "Koina is someone else's GPU"
    Predictions run on [Koina](https://koina.wilhelmlab.org/), a free, public, community-run
    inference server. pyMzLib keeps mzLib's polite throttling defaults and does not raise them.
    If you have a big job, use `max_batches` and `throttle_ms` deliberately — and be kind.

---

## What you now know

Four ideas, and they are the same idea:

| you saw | the library said |
|---|---|
| PRIDE's REST manifest is incomplete | here is the FTP walk too, and the size estimate is labelled `approximate_` |
| 13 formats have no shared shape | `views == []`, and `read_records` for when that is your file |
| MSFragger records no decoys | `is_decoy` is `None`, not a fabricated `False` |
| most RT models return iRT | `retention_time_unit` is a value you can branch on |

**When the data is incomplete or ambiguous, you get told.** Everything else in these docs is detail.

## Where next

<div class="grid cards" markdown>

- :material-file-search: **[Readers guide](guides/readers.md)**
  All 29 formats, all five ways to read them, and every caveat.

- :material-brain: **[Prediction guide](guides/prediction.md)**
  The 37 models, their constraints, and what the numbers do and don't mean.

- :material-chart-line: **[FlashLFQ guide](guides/flashlfq.md)**
  Quantify peptides across runs, with match-between-runs.

- :material-help-circle: **[FAQ](faq.md)**
  Short answers, including "why is it slow the first time".

</div>
