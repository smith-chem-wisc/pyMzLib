# pride: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR publish, because all three render
them from one spec per wire verb. The parameter, result, error and caveat tables are generated;
edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibpride). For which function answers which question, with
worked examples, see the [PRIDE Archive guide](../guides/pride.md).

**Python returns the rows, not the envelope.** `list_files()`, `list_ftp_files()`, `search()` and
`download()` each return a list: of files, of hits, of written paths. The wire's count fields are
that list's `len()`, and its size totals are `total_size_bytes(files)` and
`approximate_total_size_bytes(files)`. Python also spells three download options differently:
`destination=` for `dest`, `extensions=` (a list) for `ext`, and `overwrite=False` for
`no-overwrite`. Selecting exact files by name is its own function, `download_files()`.

## `list_files`

--8<-- "docs/reference/_generated/pride.files.md"

## `list_ftp_files`

--8<-- "docs/reference/_generated/pride.ftp-files.md"

## `search`

--8<-- "docs/reference/_generated/pride.search.md"

## `download` and `download_files`

--8<-- "docs/reference/_generated/pride.download.md"
