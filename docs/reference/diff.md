# diff: wire verbs

The facts on this page are the same facts mzLibRust and mzLibR will publish, because all three
render them from one spec per wire verb. The parameter, result, error and caveat tables are
generated; edit the spec, not this page. [How the reference facts are generated](../contributing/reference-facts.md)
explains the mechanism. For the Python signatures and result classes, see the
[API reference](../reference.md#pymzlibdiff). For worked tasks, see the
[differential design and observations guide](../guides/diff.md).

**Python spells stdin as keyword arguments.** Both verbs read a list on stdin, and pyMzLib writes
it for you. `diff design`'s model is `factors=` (`FactorSpec`s, or bare column names),
`covariates=` (`CovariateSpec`s), `interactions=` (pairs), `strata=` and `batch=`, one stdin line
per term. `diff observations`' runs are `runs=`, each an `ObservationRun` or a
`(file, sample[, fraction[, technical_replicate]])` tuple. A `None` field is not sent, so mzLib's
own default applies. The `stdin` rows below describe the wire.

**The pivot files.** `diff design --out DIR` writes, per stratum, `design_<i>.tsv` (read by
`stats fit --design` as it is) and `contrasts_<i>.tsv` (the weights mzLib's
`EmpiricalBayes.ModerateContrasts` takes, on the design file's columns). `diff observations
--responses-out` writes the peptide-by-sample table `stats fit --responses` reads. The verbs
compose through these files; no binding joins them.

**Performance:** each call starts one bridge process, about 120 ms before mzLib does any work.
`diff design` builds one small matrix per stratum: about 0.2 s in all for the 27-run MALAT1
sheet, and about 0.25 s for a 2,000-sample sheet in four strata. `diff observations` reads the
whole peptide table once, so its time follows the table: a synthetic 20,000-peptide by 50-run
table takes about 3 s with `out=`, and about 4.3 s returning its million rows through JSON, before
Python parses them. Write a large table with `out=`. Neither verb takes threads, because no mzLib
call they make takes a degree of parallelism. **Every number is mzLib's**; the bridge parses the tables and shapes the result, and
pyMzLib does no arithmetic.

## `design`

`design()` sends `--out` only when `out=` is set. `Design.stratum(key)` looks a stratum up by its
key; `Design.records` is the contrast table row by row.

--8<-- "docs/reference/_generated/diff.design.md"

## `observations`

`observations()` always sends `--basis`: it has no default, in Python or on the wire, because
mzLib has none. `--normalize` is sent only for `normalize=True`, and `--min-peptides` only when
`min_peptides=` is set. The `written` field arrives as a `pymzlib.readers.WrittenTable`.

--8<-- "docs/reference/_generated/diff.observations.md"
