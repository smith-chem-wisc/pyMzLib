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

## What this does not do

**It does not validate.** mzLib models SDRF's structural rules in `SdrfValidator` and its
vocabulary-drift rules in `SdrfDriftLint`, but both are `internal` to mzLib's Readers assembly as
of the pinned commit, so the bridge cannot reach them. Reimplementing them here — and then again
in the Rust and R bindings — is exactly the per-binding repair pyMzLib exists to avoid, and three
copies of a specification's rules is how they drift apart. The fix belongs upstream. Until it
lands, this module reads, pools and reports honestly, and makes no claim about whether a document
is *correct*.

Building an SDRF from a search (`SdrfBuilder`) and measuring how much a corpus actually says
(`SdrfCoverage`) are both reachable and not yet exposed — see
[Adding a capability](../contributing/adding-a-capability.md).
