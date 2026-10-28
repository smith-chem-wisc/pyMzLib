# stats: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR will publish, because all three
render them from one spec per wire verb. The parameter, result, error and caveat tables are
generated; edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibstats). For worked tasks, see the
[Differential abundance guide](../guides/stats.md).

**Python spells stdin as a list argument.** Each verb reads its list input on stdin, and pyMzLib
writes it for you: `stats fit`'s `stdin` is `coefficients=`, `stats adjust`'s is `p_values=` (with
`None` for an untested entry), `stats meta`'s is `studies=`, `stats combine-p`'s is `tests=`
(with `None` for a blank field), and `stats fit-logistic`'s, `stats fit-mixed`'s and
`stats fit-latent-mixed`'s are `coefficients=`, as `stats fit`'s. The `stdin` rows below describe
the wire. `stats correlate`, `stats entrapment-fdp` and `stats fit-link` read no stdin: their
tables are files, as `stats fit`'s are.

**Performance:** each call starts one bridge process, about 120 ms before mzLib does any work.
`stats fit` then fits every feature in one pass: a 10,000-feature by 24-sample table takes about
0.6 s on one thread and about 0.4 s with `threads=-1`. Test every coefficient and contrast you need in one call,
since the model is fitted once and moderated per coefficient and per contrast. Three contrasts
add about 0.7 s on that table. `estimator="marginal_likelihood"` costs about 4.5 s per prior fit
there, and mzLib fits the prior once per tested coefficient but once for all the contrasts
together: with it, test several coefficients as one-coefficient contrasts (mzLib gives identical
numbers) rather than on `coefficients`. `stats correlate` on the same
10,000 by 24 table takes about 0.4 s (0.5 s with one covariate), and `stats combine-p` on 50,000
tests about 0.8 s, most of it writing and reading the text; neither takes threads, because mzLib's
calls take no degree of parallelism. `stats fit-logistic` on a 10,000 by 24 table of 0s and 1s
takes about 1.7 s on one thread and 0.5 s with `threads=-1`; `stats fit-mixed` on a 10,000 by 24
table in six groups about 3.9 s and 2.5 s, since each feature is a one-dimensional likelihood
search. `stats fit-latent-mixed` integrates every counted sample numerically, so it is much
slower: on mzLib's reference input, with every sample counted, about 0.8 s per feature on one
thread when sigma² is estimated and 1.6 s when it is given, divided across cores by `threads=-1`; a
feature with nothing counted costs about as much as in `fit_mixed`. `stats fit-link` on 4,000 cells
takes under a second. `stats entrapment-fdp` on 100,000 targets and 100,000 entrapments at mzLib's default
thresholds takes about 0.5 s; it takes no threads, because mzLib's call takes no degree of
parallelism. **Every number is mzLib's**; the bridge parses the tables and shapes the result, and
pyMzLib does no arithmetic.

## `fit`

`fit()` sends `--contrasts`, `--confidence` and `--estimator` only when `contrasts=`,
`confidence=` and `estimator=` are set, so leaving them out is mzLib's default, and
`ModeratedContrasts.confidence_level` and `VariancePrior.estimator` report what mzLib used either
way. `coefficients` may be empty when `contrasts=` is given. The `contrasts` object is
`ModeratedFit.contrasts`, a `ModeratedContrasts` (`None` without `contrasts=`), and
`ModeratedContrasts.rows()` gives one contrast's block; `pymzlib.stats.PRIOR_ESTIMATORS` lists the
estimator names. **The three options need the pyMzLib 0.6.0 bridge.** An older bridge lists
`stats fit` but ignores them, so `fit()` checks its answer and raises `UsageError` rather than
return a fit without what you asked for.

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

## `fit-logistic`

`fit_logistic()` sends `--max-iterations` and `--tolerance` only when they are set, so leaving them
out is mzLib's default, and `LogisticModelFit.max_iterations` and `.tolerance` are then `None`.
`LogisticModelFit.status_counts` is keyed by `pymzlib.stats.LOGISTIC_STATUSES`, and
`LogisticModelFit.rows()` gives one coefficient's block.

--8<-- "docs/reference/_generated/stats.fit-logistic.md"

## `fit-mixed`

`fit_mixed()` sends `--estimator` only when `estimator=` is set; `MixedModelFit.estimator` reports
what mzLib used either way. `MixedModelFit.status_counts` is keyed by
`pymzlib.stats.MIXED_STATUSES`, and `MixedModelFit.rows()` gives one coefficient's block.

--8<-- "docs/reference/_generated/stats.fit-mixed.md"

## `fit-latent-mixed`

`fit_latent_mixed()` always sends `--link`: it has no default, in Python or on the wire, because
mzLib has none. `link="identity"` sends `identity`; a tuple sends only the numbers it holds, so a
depth slope or overdispersion left out is mzLib's default; a `DetectionLink`, or `fit_link()`'s
result, sends all four. `--responses`, `--counts`, `--depth`, `--residual-variance` and
`--estimator` are sent only when set, and `LatentMixedModelFit.estimator` and `.link` report what
mzLib used either way. `LatentMixedModelFit.status_counts` is keyed by
`pymzlib.stats.LATENT_MIXED_STATUSES`, and `LatentMixedModelFit.rows()` gives one coefficient's
block. It is a verb of its own, not an option of `fit-mixed`, because mzLib fits it with a
different entry point, optimiser and set of statuses: you choose the model by choosing the call.

--8<-- "docs/reference/_generated/stats.fit-latent-mixed.md"

## `fit-link`

`fit_link()` sends `--with-depth`, `--with-overdispersion` and `--detection-only` (`true` or
`false`) only when they are set, so leaving them out is mzLib's default; `LinkFit` reports the
values used either way. `LinkFit.link` is a `DetectionLink`, or `None` when mzLib's search did not
converge.

--8<-- "docs/reference/_generated/stats.fit-link.md"

## `entrapment-fdp`

`entrapment_fdp()` always sends `--ratio`: it has no default, in Python or on the wire, because
mzLib has none. It sends `--thresholds` (comma-separated) only when `thresholds=` is set.
`EntrapmentFdpEstimates.at()` gives one threshold's row.

--8<-- "docs/reference/_generated/stats.entrapment-fdp.md"
