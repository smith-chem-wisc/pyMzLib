# SDRF experimental design

Every other reader in pyMzLib answers *what did the search find*. **SDRF-Proteomics answers what
was searched** — which sample, which organism part, which replicate, which instrument settings —
and that is the half you need before results from two experiments can be compared at all.

```python
import pymzlib

doc = pymzlib.sdrf.read("PXD000070.sdrf.tsv")
print(doc.row_count, len(doc.columns))            # 6 31
print(doc.value("characteristics[organism part]")[0])   # human erythrocytes
```

Pool several experiments into one table, giving each a name you choose:

```python
pooled = pymzlib.sdrf.pool({
    "PXD000070.sdrf.tsv": "malaria",
    "PXD026824.sdrf.tsv": "colon",
})
print(pooled.document_count, pooled.row_count)    # 2 24
print(set(pooled.source_documents()))             # {'malaria', 'colon'}
```

## What you can ask of an SDRF

Reading is the first of seven questions. Each is answered by one mzLib type, projected once in the
bridge so that Python, Rust and R get the same answer:

| You want to know | Call | mzLib answers with |
|---|---|---|
| What does the file say, cell for cell? | [`read()`](#this-module-is-row-major-and-every-other-reader-is-columnar) | `SdrfDocument` |
| Several experiments in one table | [`pool()`](#pooling) | `SdrfCollection.Merge` |
| Is this file well-formed? | [`validate()`](#validate-a-deposit) | `SdrfValidator` |
| Does it describe its samples at all? | [`assess()`](#decide-whether-an-sdrf-is-worth-using) | `SdrfSampleInformativeness` |
| Do these files write the same thing the same way? | [`lint()`](#lint-a-pooled-set) | `SdrfDriftLint` |
| What does it say about each sample? | [`samples()`](#one-row-per-sample) | `SdrfSampleBlock` |
| How old was each sample, in years? | [`parse_ages()`](#parse-ages-safely) | `SdrfAge` |

`validate`, `assess` and `samples` also come as [`*_many`](#many-files-in-one-call), which reads a
whole corpus in one call.

**The three judging functions are blind in different places, which is why there are three.** A
file whose sample columns are all `"not available"` passes `validate()` — reserved words are the
specification's *correct* way to say nothing — and passes `lint()`, because a reserved word cannot
drift. Only `assess()` sees that it answers no biological question. Twenty files can each pass
`validate()` and still be unpoolable; only `lint()` looks across files. Every result carries
`caveats` naming the blind spot it leaves.

## This module is row-major, and every other reader is columnar

[`read_records()`](readers.md) and its siblings hand back `columns`, a name-to-values dict. This
one hands back `columns` (a **list of names**) and `rows` (a list of cell lists). That is not a
style choice — it is forced by the format:

**SDRF column names are data, not a schema, and they repeat.** In the 1,236-file curated corpus,
649 files carry `comment[modification parameters]` more than once — up to eight times in a single
file — 258 repeat `comment[sdrf template]`, and one file repeats an *empty* name 23 times. A dict
keyed by name would keep one occurrence of each and silently drop the rest.

So position is what links a name to a cell:

```python
doc.value("comment[modification parameters]")   # the FIRST one, per row
doc.all("comment[modification parameters]")     # ALL EIGHT, per row
doc.indexes_of("comment[modification parameters]")   # [19, 20, 21, ...]
```

`records` gives you the convenient dict-per-row shape, but it is **lossy when a name repeats** —
later positions overwrite earlier ones. `has_repeated_columns` tells you whether that applies
before you trust it.

## Three things to know before you index anything

### Rows are ragged

`len(row)` may be less than `len(columns)`. Real files are like this: `PXD059974` in mzLib's own
fixtures has a 46-column header with **17 of its 22 rows carrying only 42 cells**. mzLib preserves
that rather than padding, so the file round-trips byte for byte, and pyMzLib preserves it in turn.

`value()` returns `None` for a position a row does not reach. `ragged_row_count` counts the short
rows you were given, and the `caveats` list names the count for the **whole document** even when
you paged:

```python
r = pymzlib.sdrf.read("PXD059974.sdrf.tsv", limit=1)
print(r.caveats)
# ['...', '17 of 22 rows are SHORT: they carry fewer cells than column_names has entries (46)...']
```

### Cells are raw strings, never interpreted

The SDRF key=value grammar arrives exactly as written:

```python
doc.value("comment[modification parameters]")[0]
# 'NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed'
```

It is not decoded, and that is deliberate rather than unfinished. A cell containing `=` and `;`
cannot be told apart from a CV term by its shape: `comment[file uri]` cells in the corpus carry
pre-signed download URLs whose query strings contain `Signature=`, `Expires=` and `Id=` — 5,750
occurrences each. Anything that decoded "cells that look like CV terms" would decode those.

!!! note "Why this matters more than it sounds"

    `readers.read_records()` will also open a `.sdrf.tsv`, because mzLib dispatches the type. It
    gives you two columns, `header` and `cells`, each **semicolon-joined** — and since SDRF's own
    grammar is semicolon-delimited, that string cannot be split back into cells. Use
    `pymzlib.sdrf.read()` for SDRF; the joined form is not recoverable.

### A reserved word is a real value

`"not available"` and `"not applicable"` mean *the experiment stated an absence*. That is the
specification's own mechanism and a curator chose to write it. It is **not** the same as a column
the document does not have, which comes back as `None`. Do not collapse the two:

```python
doc.value("characteristics[disease]")[0]      # 'not applicable'  - stated
doc.value("characteristics[nonesuch]")[0]     # None              - no such column
```

## Pooling

`pool()` merges documents into one analysis table. Columns are the union of every document's,
ordered by SDRF's own block structure — sample metadata, then data-file metadata, then factor
values — and a name that repeats is carried at the highest multiplicity any single document used,
so nothing is dropped. A cell a document did not have is filled with `"not available"`.

### Give your documents names

```python
pooled = pymzlib.sdrf.pool({"a.sdrf.tsv": "malaria", "b.sdrf.tsv": "colon"})
```

Pass a plain list instead and mzLib falls back to `containing-folder/file-stem`, which **depends
on where the files happen to sit** — so the same two documents pooled from a different directory
produce a different table. `caveats` says so when that happened. Label every document or none; a
partial set is refused, because a provenance column that means a chosen name on some rows and a
local path on others cannot be undone downstream.

### The pooled table is for analysis, not deposit

`source name` + `assay name` + `comment[label]` is unique within one document, but two experiments
may both have a `"Sample 1"`. A pooled table will therefore usually violate SDRF's uniqueness
rule. Use `comment[source document]` as part of any key:

```python
import pandas as pd
df = pd.DataFrame(pooled.records)
df.groupby(["comment[source document]", "characteristics[organism part]"]).size()
```

Writing it back out gives a well-formed SDRF file, and `limit`/`offset` never shorten it — they
shape what comes back over the wire, not what lands on disk:

```python
pymzlib.sdrf.pool(docs, out="merged.sdrf.tsv", limit=5)   # 5 rows returned, all 24 written
```

## Validate a deposit

Before you upload an SDRF to PRIDE, or trust one you downloaded, check its structure:

```python
v = pymzlib.sdrf.validate("my_study.sdrf.tsv")
print(v.is_valid, v.error_count, v.warning_count)      # False 7 2
for m in v.errors:
    print(m.line_number, m.rule, m.column_name)
# None RequiredColumn technology type
# None RequiredColumn comment[instrument]
# ...
```

This is mzLib's `SdrfValidator`, checking **structure** against SDRF-Proteomics v1.1.0: the
thirteen required columns, the three recommended ones, column order and casing, malformed names
(`characteristics [organism]` with a space), rows whose width differs from the header, replicate
and fraction columns that are not positive integers, reserved-word casing, and the one hard row
rule — `source name` + `assay name` + `comment[label]` must be unique.

**Errors and warnings mean different things, and `is_valid` only looks at errors.** An error means
the file cannot be reliably consumed: a ragged row reads every later cell under the wrong column, a
missing `comment[data file]` leaves nothing to join on. A warning is worth fixing and never worth
refusing. mzLib calibrated every severity against the 1,236 curated files of
[bigbio/sdrf-annotated-datasets](https://github.com/bigbio/sdrf-annotated-datasets): a rule that
fired on most curated files was judged a wrong rule, not 1,236 wrong files.

**Where a finding points.** `line_number` is the line in the file, counting the header as line 1,
so you can jump straight to it in an editor. `row_index` is the same place as a 0-based index into
`read(...).rows`. **Both are `None` for a finding about the whole document** — a missing column, an
empty header — never `0`, which would name a real row. `column_name` is `None` when the finding is
not about one column.

The findings are a table, so the usual summaries are one line of pandas:

```python
import pandas as pd
pd.DataFrame(v.columns).groupby(["severity", "rule"]).size()
# Error    RequiredColumn       7
# Warning  RecommendedColumn    2
```

!!! warning "What a clean result does not mean"
    **Controlled-vocabulary terms are not resolved.** `AC=MS:1003028;NT=Q Exactive` is well-formed
    and passes, even though MS:1003028 is the Orbitrap Exploris 480. And a file that says nothing
    — every sample cell `not available` — validates cleanly. [`assess()`](#decide-whether-an-sdrf-is-worth-using)
    is the check for that.

## Decide whether an SDRF is worth using

A file generated from a list of data files — reserved words in every sample column, one replicate
number, no factor — is valid SDRF. For grouping results by biology it is the same as having no SDRF
at all. `assess()` is the gate:

```python
a = pymzlib.sdrf.assess("PXD000070.sdrf.tsv")
print(a.verdict)                  # Partial
print(a.factor_value_varies, a.sample_is_described, a.biological_replicate_varies)
# False True True
```

mzLib's `SdrfSampleInformativeness` asks three questions of the **sample** half of the file:

| Check | Passes when |
|---|---|
| `factor_value_varies` | some `factor value[...]` column holds two or more different real answers — the file says what the study varied |
| `sample_is_described` | some `characteristics[...]` column **other than organism and biological replicate** holds a real answer. Organism is left out because a search can fill it without a human ever describing the samples |
| `biological_replicate_varies` | `characteristics[biological replicate]` holds two or more different real answers |

| Verdict | Means |
|---|---|
| `Informative` | all three pass: the file can drive an experimental design |
| `Partial` | some pass. **Often legitimate** — a single-condition study has nothing to vary — so whether it is good enough is your call; the three booleans say which check failed |
| `Skeleton` | none pass. On the curated corpus mzLib measured 300 Informative, 742 Partial and 194 Skeleton |

"A real answer" means neither empty nor a reserved word, and different answers are compared
ignoring case and surrounding space, so `liver` and `Liver` are one value. The evidence is in
`columns` — one row per column each check read, with `rows`, `filled`, `absent`,
`distinct_values` and `fill_rate` (a fraction, 0 to 1):

```python
pd.DataFrame(pymzlib.sdrf.assess("skeleton.sdrf.tsv").columns)
#                     role                            column_name  rows  filled  absent  distinct_values  fill_rate
# 0  sample_characteristic               characteristics[disease]     2       0       2                0          0
# 1  sample_characteristic         characteristics[organism part]     2       0       2                0          0
# 2   biological_replicate  characteristics[biological replicate]     2       2       0                1          1
```

To gate a corpus, assess it in one call and keep what passes:

```python
batch = pymzlib.sdrf.assess_many(paths, threads=-1, on_error="skip")
print(batch.verdict_counts)       # {'informative': 1, 'partial': 0, 'skeleton': 1}
usable = batch.paths_with("Informative", "Partial")
```

`paths_with()` refuses a verdict it does not know — `"informative"` in lower case would otherwise
select nothing, silently.

## Lint a pooled set

Validity is a property of one file; comparability is a property of the relationship between files.
Before you [`pool()`](#pooling) experiments and group by a column, check that they wrote it the same
way:

```python
drift = pymzlib.sdrf.lint({"cohort.sdrf.tsv": "cohort", "partner.sdrf.tsv": "partner"})
pd.DataFrame(drift.columns)[["finding_index", "kind", "column_name", "value", "documents"]]
#    finding_index                   kind                column_name                       value  documents
# 0              0  AccessionNameConflict        comment[instrument]                Exploris 480  [partner]
# 1              0  AccessionNameConflict        comment[instrument]       Orbitrap Exploris 480   [cohort]
# 2              1   MixedTermAndFreeText             comment[label]  controlled vocabulary term   [cohort]
# 3              1   MixedTermAndFreeText             comment[label]                   free text  [partner]
# 4              2      ColumnNameVariant                       None        Characteristics[sex]  [partner]
# 5              2      ColumnNameVariant                       None        characteristics[sex]   [cohort]
# 6              3       ValueCaseVariant  characteristics[organism]                Homo sapiens   [cohort]
# 7              3       ValueCaseVariant  characteristics[organism]                homo sapiens  [partner]
```

The table is **one row per finding × variant**: finding 3 is one concept written two ways, so it is
two rows sharing a `finding_index`, the majority spelling first (`variant_rank` 0). `documents` lists
the labels of the files that used each spelling — which is why you should
[label them](#give-your-documents-names). `drift.findings()` gives the same rows grouped.

| `kind` | What was inconsistent | Why it matters |
|---|---|---|
| `NameAccessionConflict` | one name, several accessions | **the serious one**: a join on the accession splits the concept, a join on the name merges two |
| `AccessionNameConflict` | one accession, several names | harmless to a join on the accession; the earliest sign two files were annotated by different hands |
| `MixedTermAndFreeText` | a CV term in some files, free text in others | the two never compare equal, so half the corpus drops out of any query |
| `ColumnNameVariant` | column names differing only by case or spacing | they are different columns to every consumer. `column_name` is `None`: the finding is about a name |
| `ValueCaseVariant` | free-text values differing only by case or spacing | equality grouping treats `Homo sapiens` and `homo sapiens` as two populations |

!!! note "The majority is a description, not advice"
    `variant_rank` 0 is what most of *these* files did. In the curated corpus
    `characteristics[organism]` is free text in 963 documents and a CV term in 275, so following the
    majority would throw the accession away. Where variants are merely spellings of one thing, match
    the majority; where they differ in kind, keep the richer one.

**Columns whose values are unique per row are never compared** — `source name`, `assay name`,
`comment[data file]`, `comment[searched data file]` (since mzLib 1.0.592, #1335), `comment[file uri]`
and `comment[source document]` — so two studies naming different raw files is data, not drift.
Reserved words and empty cells are skipped too, so a set of files that say nothing lints clean.

## One row per sample

An SDRF has one row per *data file*; a sample measured in eight fractions is eight rows.
`samples()` turns that back into samples — mzLib's `SdrfSampleBlock`, keyed on `source name`:

```python
s = pymzlib.sdrf.samples("cohort.sdrf.tsv")
print(s.sample_count, s.row_count)     # 6 12
print(s.conflicts())                   # [('S6', 'characteristics[disease]')]
```

The table is long — one row per sample × column × position — so repeated columns and ragged files
cost nothing:

| column | meaning |
|---|---|
| `source_name` | the sample, as its first row spelled it. Matched ignoring case and surrounding space, as mzLib keys it |
| `sample_row_count` | how many rows (runs, fractions) carry the sample |
| `column_name` | the header's own spelling — `characteristics[Age]` stays itself |
| `column_kind` | `source_name`, `characteristic` or `factor_value` |
| `position` | 0, 1, … for the first, second, … occurrence of a repeated column |
| `status` | `agreed`, or `conflicting` |
| `value` | the cell, verbatim — reserved words included |
| `age_*` | on `characteristics[age]` rows only: [the parsed age](#parse-ages-safely) |

**A conflicting column is withheld, not guessed.** When a sample's rows disagree — S6's first
fraction says `COVID-19` and its second says `normal` — copying either through would be a "first
row wins" rule, which is how a sample silently acquires the wrong disease. mzLib withholds the
column and names it; it appears here once, with `status == "conflicting"` and `value` `None`. That
is not the same as a column the file lacks, which has no row at all.

Across files a source name is not unique — two studies may both have `Sample 1` — so key samples on
`(source_index, source_name)` in a [`samples_many()`](#many-files-in-one-call) result.

## Parse ages safely

`characteristics[age]` is free text in practice. mzLib's `SdrfAge` reads it into **years** and,
more importantly, refuses what it cannot read without guessing:

```python
ages = pymzlib.sdrf.parse_ages(["58Y", "40Y-85Y", ">=90Y", "63", "6-8 weeks"])
pd.DataFrame(ages.columns)
#         cell  years  min_years  max_years   precision  follows_specification  refusal
# 0        58Y  58.00      58.00      58.00       Exact                   True     None
# 1    40Y-85Y  62.50      40.00      85.00       Range                   True     None
# 2      >=90Y  90.00      90.00       None  LowerBound                   True     None
# 3         63   None       None       None        None                   None  no_unit
# 4  6-8 weeks   0.13       0.11       0.15       Range                  False     None
```

**Units are years**: a month is 1/12 year, a week 7/365.25, a day 1/365.25, an hour
1/(365.25 × 24). `years` is the single figure for an age axis — the age, a range's **midpoint**, or
a bound — so for a wide cohort range it is a poor stand-in for any one sample. Filter on `precision`
before you plot:

| `precision` | Example | `years` | `min_years` | `max_years` |
|---|---|---|---|---|
| `Exact` | `58Y`, `30Y6M`, and `40Y-40Y` (a degenerate range, since mzLib #1333) | the age | the age | the age |
| `Range` | `40Y-85Y`, `6-8 weeks` | the midpoint | the low end | the high end |
| `LowerBound` | `>=90Y`, the usual de-identification cap | the bound | the bound | **`None` = unbounded** |
| `UpperBound` | `<1Y` | the bound | `0` | the bound |

**`max_years` is `None` for two different reasons, and `precision` tells them apart.** For a lower
bound, mzLib's maximum is +infinity, which JSON cannot carry; for a refused cell there is no age at
all. `follows_specification` is `False` for unambiguous words (`3 year`, `6-8 weeks`) — read, but
not the specification's `nYnMnD` grammar.

**A refused cell has no age. Never fill it with a guess.** The reason is in `refusal`:

| `refusal` | The cell | What to do |
|---|---|---|
| `no_unit` | a bare number, `63` — 11% of real age cells | ask the curator: 63 years and 63 days are both plausible in one study |
| `reserved_word` | `not available`, `anonymized`, … in any case | the experiment stated there is no age; respect it |
| `empty` | blank, or `None` in your list | the cell is missing |
| `unreadable` | free text, `about forty` | ask the curator |

### Worked example: ages across a corpus

`samples_many()` parses every `characteristics[age]` cell for you, on the same row as the sample it
belongs to. From there it is ordinary pandas:

```python
import pandas as pd
import pymzlib

batch = pymzlib.sdrf.samples_many(paths, threads=-1, on_error="skip")
df = pd.DataFrame(batch.columns)

ages = df[(df.column_name.str.lower() == "characteristics[age]") & (df.status == "agreed")]
point = ages[ages.age_precision == "Exact"]            # one age per sample: safe to plot
capped = ages[ages.age_precision == "LowerBound"]      # ">=90Y": at least 90, no upper bound
refused = ages[ages.age_refusal.notna()]               # no age - report, do not impute

print(point.groupby("source_path").age_years.describe())
print(refused.groupby("age_refusal").size())
```

!!! warning "pandas turns `None` into `NaN`"
    In a float column pandas stores `None` as `NaN`, so `age_max_years` cannot tell you *why* it is
    missing once it is in a DataFrame. Always decide with `age_precision` and `age_refusal`, which
    are strings and keep their meaning; never with `age_max_years.isna()`.

For ages that did not come from `samples()` — a column of `read()`, or PRIDE project metadata —
pass the cells to `parse_ages()` directly. It takes any number of cells in one bridge call, keeps
blank and `None` cells in place, and returns row `i` for cell `i`:

```python
doc = pymzlib.sdrf.read("PXD026824.sdrf.tsv")
ages = pymzlib.sdrf.parse_ages(doc.value("characteristics[age]"))
```

## Many files in one call

`validate_many()`, `assess_many()` and `samples_many()` read a list of files in **one** bridge
process:

```python
batch = pymzlib.sdrf.validate_many(paths, threads=-1, on_error="skip")
print(batch.file_count, batch.read_count, batch.valid_count)
for f in batch.files:
    if f.error:
        print(f.path, f.error.kind, f.error.message)
```

- **Why not a loop?** Each bridge call costs about 120 ms before mzLib does anything; a
  1,236-file corpus spends two and a half minutes just starting processes. One call spends it once,
  and the parallelism happens inside .NET, where no Python thread pool is needed.
- **`threads`** is how many files are read at once: `1` by default, because each worker holds a
  whole document in memory; `-1` for every core. **The result is the same at any value** — rows are
  in input order, then the file's own order — which is tested byte for byte.
- **`on_error="fail"`** (the default) raises on the first unreadable file *in input order*, so the
  same bad batch fails with the same message at any thread count. **`"skip"`** records the failure
  on that file's entry in `files` — a `FileError` with `kind` `usage` (missing) or `correctness`
  (mzLib could not read it) — and carries on. Every input keeps its entry either way.
- The table gains `source_index` (the file's 0-based position in your list) and `source_path` as its
  first two columns.

`lint()` and `pool()` take many files by nature and use the `{path: label}` form instead.

## `comment[searched data file]`

mzLib 1.0.592's SDRF builder (#1327) can write `comment[searched data file]` beside
`comment[data file]`: the file a search actually read (an `.mzML`) next to the one the instrument
acquired (a `.raw`). pyMzLib has no builder, but reads and pools such files like any other — the
column is carried verbatim by `read()` and `pool()`. When it is present, **join search results on it
rather than on `comment[data file]`**: it names the file whose scans the identifications refer to.
`lint()` treats it as unique per row, so differing file names are never reported as drift.

## What this does not do

- **Build an SDRF.** mzLib's `SdrfBuilder` (and its 1.0.592 additions — the searched data file
  #1327, bare `comment[label]` #1307, atomic writing #1324) takes mzLib-native types — `Tolerance`,
  `DigestionAgent`, `Modification`, `DissociationType` — that have no representation on a
  language-neutral wire yet. Designing one is a contract change, not a projection, so it is not
  offered. `pool(out=...)` writes an SDRF, and gets #1324's atomic write for free.
- **Resolve controlled-vocabulary terms.** Neither `validate()` nor `lint()` checks that an
  accession is real or names what it claims.
- **Audit quantification design.** mzLib's `SdrfQuantAuditor` (is this TMT experiment channel-level
  or kit-only, which files are on disk) is public but not projected yet: it is a separate question
  with its own eleven-fact report, and needs its own verb.
- **Repair anything.** Every function here reports; none rewrites a cell. Fixing a curated file is
  a curation decision, and a reader that made it silently would make it irreversibly.

## References

- Dai C., Füllgrabe A., Pfeuffer J., *et al.* A proteomics sample metadata representation for
  multiomics integration and big data analysis. *Nature Communications* **12**, 5854 (2021).
  [doi:10.1038/s41467-021-26111-3](https://doi.org/10.1038/s41467-021-26111-3) — the SDRF-Proteomics
  format this module reads.
- The specification itself:
  [bigbio/proteomics-sample-metadata](https://github.com/bigbio/proteomics-sample-metadata). mzLib
  validates against v1.1.0.
- The curated corpus mzLib's severities and verdict counts were calibrated on:
  [bigbio/sdrf-annotated-datasets](https://github.com/bigbio/sdrf-annotated-datasets).
