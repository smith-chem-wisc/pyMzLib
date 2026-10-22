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

**`annotate_go` takes the table and the database positionally** (`groups`, `database`), as the wire's
`--groups` and `--database`; every other option is a keyword with the wire's name in snake case.
Its rows are windowed by `limit`/`offset`, while `out=` always receives the whole table, so a large
run is `out=` plus `limit=0`.

**Multi-table output is JSON only.** `proteins read` has no `out` option: its three tables would
need three files. **Performance:** each call starts one bridge process, about 120 ms before mzLib
does any work, and each database is read whole. Pass all your databases in one call rather than
looping, and raise `threads` only if memory allows one whole database per thread.

## `read`

--8<-- "docs/reference/_generated/proteins.read.md"

## `resolve_genes`

--8<-- "docs/reference/_generated/genes.resolve.md"

## `build_orthology`

**`genes orthology` is a wire verb of the `genes` module, projected as
`pymzlib.proteins.build_orthology`.** Both of its input lists travel on stdin: the species lines
first, then a line holding only `--`, then one homology dump per line. Python renders the species
lines from two dicts, `gtfs=` (a species' GTF) and `gene_sets=` (a compact gene table instead), so
the table's `species` row is those two keywords, and `dumps` is Python's `dumps=` list. The wire's
`restrict_gene_trees` flag is Python's `restrict_gene_trees=True`, and `pair_status` is a
`(from_species, to_species)` tuple, sent as `FROM,TO`. `out` is a **directory**, not a `.tsv`.

The relationships themselves are in the Parquet files `out=` writes, and mzLib ships no reader for
them; pyMzLib reads none either, because it has no third-party dependency. The
[Protein databases guide](../guides/proteins.md#orthologs-across-species-an-ensembl-compara-snapshot)
shows how to query them with DuckDB or pyarrow. **Performance:** every input is read whole into
memory. For human, mouse and rat at release 116, mzLib #1381 reports about 20 s to build and 40 s
more to write 13 files (9.0 MB), and the vertebrate gene-tree file holds 4.2 million genes unless
`restrict_gene_trees=True` keeps only yours. No mzLib call here takes a thread count.

--8<-- "docs/reference/_generated/genes.orthology.md"

## `classify_peptides`

--8<-- "docs/reference/_generated/proteins.classify-peptides.md"

## `annotate_go`

--8<-- "docs/reference/_generated/proteins.annotate-go.md"

## `update_go`

--8<-- "docs/reference/_generated/proteins.update-go.md"

## `align` and `align_many`

**The pairs always travel on stdin, one per line, source and target separated by a tab.**
`align_many(pairs)` sends your list as those lines; `align(source, target)` sends one line and
returns that pair's `Alignment`. Every option you leave as `None` is not sent, so mzLib's
`PairwiseAligner` default applies, and the result's `gap_open`, `gap_extend`, `end_gaps` and
`max_cells` say what was used. The wire's `no_residue_map` flag is Python's `residue_map=False`.
**Performance:** time and memory grow with the product of the two lengths (one byte per table cell),
and the pairs of one call are aligned one after another, so send a batch as one `align_many` call
rather than looping over `align`.

--8<-- "docs/reference/_generated/proteins.align.md"

## `parse_accessions`

**Plural in Python, singular on the wire,** as `sdrf parse-age` is: each stdin line is one
accession, and Python's list becomes those lines exactly as given, never stripped or split.

--8<-- "docs/reference/_generated/proteins.parse-accession.md"
