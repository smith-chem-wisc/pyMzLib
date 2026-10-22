# stats: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR will publish, because all three
render them from one spec per wire verb. The parameter, result, error and caveat tables are
generated; edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibstats). For worked tasks, see the
[Differential abundance guide](../guides/stats.md).

**Python spells stdin as a list argument.** Each verb reads its list input on stdin, and pyMzLib
writes it for you: `stats fit`'s `stdin` is `coefficients=`, `stats adjust`'s is `p_values=` (with
`None` for an untested entry), `stats meta`'s is `studies=` and `stats combine-p`'s is `tests=`
(with `None` for a blank field). The `stdin` rows below describe the wire. `stats correlate` reads
no stdin: its three tables are files, as `stats fit`'s are.

**Performance:** each call starts one bridge process, about 120 ms before mzLib does any work.
`stats fit` then fits every feature in one pass: a 10,000-feature by 24-sample table takes about
0.6 s on one thread and about 0.4 s with `threads=-1`. Test every coefficient you need in one call,
since the model is fitted once and moderated per coefficient. `stats correlate` on the same
10,000 by 24 table takes about 0.4 s (0.5 s with one covariate), and `stats combine-p` on 50,000
tests about 0.8 s, most of it writing and reading the text; neither takes threads, because mzLib's
calls take no degree of parallelism. **Every number is mzLib's**; the bridge parses the tables and
shapes the result, and pyMzLib does no arithmetic.

## `fit`

--8<-- "docs/reference/_generated/stats.fit.md"

## `adjust`

--8<-- "docs/reference/_generated/stats.adjust.md"

## `meta`

--8<-- "docs/reference/_generated/stats.meta.md"

## `combine-p`

`combine_p()` always sends `--method`: it has no default, in Python or on the wire, because mzLib
has none. Python's tuples and mappings become the stdin lines; a `None` weight or effect is a
blank field. `CombinedPValues.p_values` is the `p_value` column keyed by feature.

--8<-- "docs/reference/_generated/stats.combine-p.md"

## `correlate`

`correlate()` sends `--given` only when `given=` is set. `Correlations.method_counts` is keyed by
`pymzlib.stats.SPEARMAN_METHODS`.

--8<-- "docs/reference/_generated/stats.correlate.md"
