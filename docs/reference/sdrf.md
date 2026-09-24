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
