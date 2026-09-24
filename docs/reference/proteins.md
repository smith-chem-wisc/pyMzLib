# proteins: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibproteins). For worked tasks, see the
[Protein databases guide](../guides/proteins.md).

**Python spells the database inputs differently from the wire.** Each function takes one database
or a list as its first argument, and contaminant databases as `contaminants=`. pyMzLib chooses
`--path`, `--contaminant` or `--paths-stdin` from those, so the tables' `path`, `contaminant` and
`paths_stdin` rows describe the wire, not Python keywords. `accessions_stdin` is Python's
`accessions=` list. `genes resolve` is a wire verb of the `genes` module, and Python projects it as
`pymzlib.proteins.resolve_genes`.

**Multi-table output is JSON only.** `proteins read` has no `out` option: its three tables would
need three files. **Performance:** each call starts one bridge process, about 120 ms before mzLib
does any work, and each database is read whole. Pass all your databases in one call rather than
looping, and raise `threads` only if memory allows one whole database per thread.

## `read`

--8<-- "docs/reference/_generated/proteins.read.md"

## `resolve_genes`

--8<-- "docs/reference/_generated/genes.resolve.md"

## `classify_peptides`

--8<-- "docs/reference/_generated/proteins.classify-peptides.md"
