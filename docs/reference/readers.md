# readers: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibreaders). For choosing between the readers, see the
[Readers guide](../guides/readers.md).

Parameters are listed under their Python keyword. Where it differs from the wire option, the wire
spelling is given too. A unit of "—" means the value has none: a count of things, a flag or a
name.

**Performance, for every verb here:** each call starts one bridge process, about 120 ms before
mzLib does any work, and mzLib parses the whole file on every call whatever `limit` and `offset`
say. For a large file, make one call with `out` instead of paging.

**Many files at once, for every verb here but `formats`:** the wire verb also takes
`--paths-stdin`, and pyMzLib spells that form `<function>_many` - `read_spectra_many`,
`read_records_many`, `identify_many` and so on. One process reads the whole list, `threads` files
at a time, into one long table whose first columns are `source_index` and `source_path`, with one
`files[]` entry per input. The output is byte-identical at any thread count. Each verb's
parameters and fields for that form are in its spec's `bulk` section, and the
[readers guide](../guides/readers.md#many-files-at-once) walks through it.

## `formats`

--8<-- "docs/reference/_generated/readers.formats.md"

## `read_records`

--8<-- "docs/reference/_generated/readers.read-records.md"

## `read_spectra`

--8<-- "docs/reference/_generated/readers.read-spectra.md"

## `identify`

--8<-- "docs/reference/_generated/readers.identify.md"

## `read_results`

--8<-- "docs/reference/_generated/readers.read-results.md"

## `read_features`

--8<-- "docs/reference/_generated/readers.read-features.md"

## `read_matches`

--8<-- "docs/reference/_generated/readers.read-matches.md"

## `read_protein_groups`

--8<-- "docs/reference/_generated/readers.read-protein-groups.md"

## `read_quantified_peptides`

--8<-- "docs/reference/_generated/readers.read-quantified-peptides.md"

## `read_occupancy`

--8<-- "docs/reference/_generated/readers.read-occupancy.md"
