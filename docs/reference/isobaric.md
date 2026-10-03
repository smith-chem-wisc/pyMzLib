# isobaric: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibisobaric). For worked tasks, see the
[isobaric kits guide](../guides/isobaric.md).

`kits` takes no file and makes no network call. Python's `kits("TMT18")` is the wire's
`--kit TMT18`; with no argument every kit is listed. The result's `kits` field crosses as a list of
`{kit, channel_count}` summaries, and Python's `IsobaricKits.kits` regroups the table under them.

## `kits`

--8<-- "docs/reference/_generated/isobaric.kits.md"
