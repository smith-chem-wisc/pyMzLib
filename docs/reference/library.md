# library: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzliblibrary). For worked tasks, see the
[spectral-library guide](../guides/library.md).

Parameters are listed under their Python keyword. Where it differs from the wire option, the wire
spelling is given too. A unit of "—" means the value has none: a count of things, a flag or a
name.

## `read`

`read` reads one local `.msl` file and makes no network call. It runs through the same table
engine as the [readers](readers.md) verbs, so `limit`, `offset` and `out` behave as they do there,
with one difference: they count **precursors**, and with `fragments=True` each precursor gives one
row per fragment ion. Python's `read_many()` is the wire's `--paths-stdin` form, one process and one
long table for a list of libraries; its parameters and fields are in the spec's `bulk` section, and
its per-library facts arrive as `LibraryFileReport` objects. `index_only_threshold` is sent only
when given, so mzLib's own default applies otherwise.

**Performance:** each call starts one bridge process, about 120 ms before mzLib does any work, and
opens the library again whatever `limit` and `offset` say. A library under the index-only threshold
is loaded whole, fragments included, even when only precursors are returned. Index-only without
fragments reads no fragment block at all; index-only with fragments reads each block twice. For a
whole library with fragments, use `out`.

--8<-- "docs/reference/_generated/library.read.md"
