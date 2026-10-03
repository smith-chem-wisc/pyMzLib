# peptidoform: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibpeptidoform). For worked examples, see the
[Peptidoforms guide](../guides/peptidoforms.md).

**Python groups five census fields into one object.** `annotated_modification_sites`,
`annotated_modifications_loaded`, `uniprot_annotated_features`, `unresolved_modifications` and
`uniprot_features_by_type` arrive as `Digest.modification_census` (`sites`, `applied`,
`annotated`, `unresolved`, `by_type`), which can also `explain()` itself in a sentence. The wire's
`no-modifications` flag is `modifications=False`, and `max-mods` is `max_modifications=`.

## `fragments`

--8<-- "docs/reference/_generated/peptidoform.fragments.md"

## `convert`

`convert` reads no file and makes no network call. Python's `source=`, `target=` and `mode=` are
the wire's `--from`, `--to` and `--mode`, and the sequences travel on stdin, one per line. Python
always sends all three options, with the wire's own defaults. `SequenceConversions.sequences` and
`not_converted` regroup the table's rows as `ConvertedSequence` objects, and `outputs` is the
`output` column.

--8<-- "docs/reference/_generated/peptidoform.convert.md"
