# ptmqtl: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR will publish, because all three
render them from one spec per wire verb. The parameter, result, error and caveat tables are
generated; edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibptmqtl). For worked tasks, see the
[PTM occupancy and PTM-QTL guide](../guides/ptmqtl.md).

**Python spells lists as lists.** `modifications` and `categories` are tab-separated on the wire,
because a modification name holds spaces and commas but never a tab; pyMzLib takes a list and joins
it. The flags `keep-ceiling` and `remove-run-level` are `keep_ceiling=True` and
`remove_run_level=True`; `False`, the default, sends nothing, so mzLib's own default applies. Every
other option left as `None` is not sent either, and the result reports the value mzLib used.

**Performance:** each call starts one bridge process, about 120 ms before mzLib does any work. None
of these verbs takes `--threads`: mzLib's occupancy and pair calls take no degree of parallelism,
and `trait-effects` fits its sites with mzLib's mixed model at mzLib's own setting, whose results do
not depend on it. **Every number is mzLib's**; the bridge reads the tables and shapes the result,
and pyMzLib does no arithmetic.

## `occupancy`

`occupancy()` sends `--samples`, `--modifications`, `--categories` and `--out` only when they are
set. With `out=`, `SiteOccupancy.columns` is `None` and `SiteOccupancy.written` says where the
table went; that file is what `trait-effects` reads as `--occupancy`.

--8<-- "docs/reference/_generated/ptmqtl.occupancy.md"

## `pairs`

`pairs()` always sends `--type`: it has no default, in Python or on the wire, because mzLib has
two methods and no preference between them.

--8<-- "docs/reference/_generated/ptmqtl.pairs.md"

## `trait-effects`

`TraitEffects.status_counts` is keyed by `pymzlib.ptmqtl.TRAIT_STATUSES`, and `min_runs`,
`min_replicates`, `logit_epsilon` and `exclude_ceiling` report the options mzLib used, yours or
its defaults.

--8<-- "docs/reference/_generated/ptmqtl.trait-effects.md"
