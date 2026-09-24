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

## `formats`

--8<-- "docs/reference/_generated/readers.formats.md"

## `read_records`

--8<-- "docs/reference/_generated/readers.read-records.md"

## `read_spectra`

--8<-- "docs/reference/_generated/readers.read-spectra.md"
