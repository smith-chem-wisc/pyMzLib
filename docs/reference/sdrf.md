# sdrf: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibsdrf). For which verb answers which question - and a worked
example of each - see the [SDRF guide](../guides/sdrf.md).

`validate`, `assess` and `samples` take one path or many (`--paths-stdin`). Python spells the two
forms as two functions, `validate()` and `validate_many()`, so the bulk-only options (`threads`,
`on_error`) appear only on the `_many` form; the tables list the wire options for both. A unit of
"—" means the value has none: a count of things, a flag or a name.

## `read`

--8<-- "docs/reference/_generated/sdrf.read.md"

## `pool`

--8<-- "docs/reference/_generated/sdrf.pool.md"

## `validate` and `validate_many`

--8<-- "docs/reference/_generated/sdrf.validate.md"

## `lint`

--8<-- "docs/reference/_generated/sdrf.lint.md"

## `assess` and `assess_many`

--8<-- "docs/reference/_generated/sdrf.assess.md"

## `samples` and `samples_many`

--8<-- "docs/reference/_generated/sdrf.samples.md"

## `parse_ages`

--8<-- "docs/reference/_generated/sdrf.parse-age.md"

## `design`

**Python spells two options differently from the wire.** `condition_columns=` is a list; pyMzLib
joins it with tabs for `condition-columns`, because a tab is the one character an SDRF column name
cannot contain. `searched_files=` is a list; pyMzLib sends it on stdin and sets
`searched-files-stdin`. The result's design coordinates are 0-based on the wire and in Python;
`ExperimentalDesign.tsv`, which `out=` writes, is 1-based.

--8<-- "docs/reference/_generated/sdrf.design.md"

## `isobaric_design`

**The numbers are 1-based, on the wire and in Python.** `fraction`, `technical_replicate`,
`biological_replicate` and `plex_id` are the SDRF's own numbers, as mzLib's `SdrfIsobaricDesign`
keeps them and `TmtDesign.txt` holds them. `design` above is 0-based. Neither the bridge nor
pyMzLib converts either one, so never compare the two tables without adding or removing 1.

`kit=`, `plex_column=`, `plex_file_name_pattern=` and `single_plex=` are sent only when given.
`condition_columns=` and `searched_files=` cross exactly as for `design`. `out=` names a **folder**,
not a file: mzLib writes MetaMorpheus's `TmtDesign.txt` inside it.

--8<-- "docs/reference/_generated/sdrf.isobaric-design.md"

## `restrict`

**Python spells stdin as `searched_files=`.** A list sends one searched file per line. A
`{searched: acquired}` mapping adds `<TAB>acquired` to each line whose value is not `None`.
`RestrictedSdrf.columns` is the wire's `column_names`, as for `read`.

--8<-- "docs/reference/_generated/sdrf.restrict.md"
