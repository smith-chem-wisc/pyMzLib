# Differential design and observations

**The problem.** Before any differential test there are two tables to get right, and both are
easy to get quietly wrong. The **design**: which coefficient means what, which level is the
baseline, whether two things you want to separate can be separated by these samples at all. And
the **observations**: what each peptide measured in each biological sample, with fractions summed
and technical replicates averaged correctly, and a missing value kept as missing rather than
written as 0.

`pymzlib.diff` hands both to mzLib. `design` turns the sample sheet you already keep into mzLib's
design matrix, which equals R's `model.matrix`, and says which comparisons cannot be estimated and
why. `observations` turns a FlashLFQ or MetaMorpheus peptide table into every peptide in every
sample, as a log2 intensity or `None` with the reason. Both write files that
[`stats.fit`](stats.md) reads, so the analysis needs no glue code of yours.

Every `>>>` example on this page is executed in CI, against recordings of the real bridge, on the
files pyMzLib ships in `pkg/python/tests/fixtures/diff/`.

## Which function do I want?

| You want to know | Call | mzLib computes it with | Method |
|---|---|---|---|
| What is my design matrix, and which comparisons can these samples answer? | [`diff.design`][pymzlib.diff.design] | `AnalysisDesign.ByStratum`, `.DefaultContrasts`, `.Estimability`, `.FittedWeights` | treatment coding as R's `model.matrix` (Wilkinson and Rogers 1973) |
| What did each peptide measure in each sample, and if nothing, why? | [`diff.observations`][pymzlib.diff.observations] | `FlashLfqObservations.FromPeptideTable`, `ObservationTable` | fractions summed, technical replicates averaged in log2 |
| How far apart do the samples sit on the peptides they share? | [`diff.observations`][pymzlib.diff.observations] | `SharedPeptideNormalization.Measure` / `.Apply` | median of ratios to a geometric-mean reference (as Anders and Huber 2010) |

Then test with [`stats.fit`](stats.md): it reads the design file and the contrasts file `design`
writes, and the responses file `observations` writes.

## 1. From a sample sheet to a design

The [differential-abundance guide](stats.md#2-a-real-experiment-with-a-known-answer) analyses a
MALAT1 dilution series - 27 runs at 500, 250 or 125 ng - with a hand-made design of 0/1 columns.
Here is the same experiment as a scientist writes it down: one row per run, and its load.

```pycon
>>> with open("malat_dilution_samples.tsv") as fh:
...     header, *runs = [line.rstrip("\n").split("\t") for line in fh]
>>> header
['run', 'malat_load']
>>> from collections import Counter
>>> Counter(load for _, load in runs)
Counter({'250ng': 12, '125ng': 8, '500ng': 7})

```

Say which column is a factor and which level is the baseline. Everything else is mzLib's:

```pycon
>>> design = pymzlib.diff.design("malat_dilution_samples.tsv",
...                              factors=[pymzlib.diff.FactorSpec("malat_load", "500ng")],
...                              out="malat_design")
>>> stratum = design.stratum("all")
>>> stratum.column_names
['intercept', 'malat_load=125ng', 'malat_load=250ng']
>>> stratum.rank, stratum.aliased_terms, stratum.warnings
(3, [], [])
>>> for row in design.records:
...     print(row["contrast_id"], row["label"], row["estimability"], row["fitted_weights"])
c1 malat_load=125ng vs malat_load=500ng estimable [0, 1, 0]
c2 malat_load=250ng vs malat_load=500ng estimable [0, 0, 1]

```

What came back:

- **The columns** are named `factor=level`. The reference comes first and has no column; the
  other levels follow in ordinal (character) order, so `125ng` before `250ng`. That is R's
  `model.matrix` with the same level order.
- **The default contrasts** are each level against the reference. Each is one design column here,
  so `fitted_weights` has a single 1.
- **`out=`** wrote `design_1.tsv`, the design table [`stats.fit`](stats.md) reads, and
  `contrasts_1.tsv`, each estimable contrast's weights on the same columns, which
  `stats.fit(contrasts=...)` reads.

Fit it. Each contrast against the reference is its design column, tested by name:

```pycon
>>> fit = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", stratum.design_file,
...                         ["malat_load=250ng", "malat_load=125ng"])
>>> fit.status_counts
{'fitted': 212, 'too_few_observations': 214, 'rank_deficient': 66}
>>> import statistics
>>> def median_change(coefficient):
...     changes = [r["estimate"] for r in fit.rows(coefficient)
...                if r["status"] == "fitted" and r["feature"].startswith("MALAT1:")]
...     return round(statistics.median(changes), 2)
>>> median_change("malat_load=250ng"), median_change("malat_load=125ng")
(-0.89, -1.01)

```

The same numbers as the hand-made design, because it is the same model: the bridge's tests hold
the two fits equal, feature by feature.

**Or test the contrasts file itself**, and get a confidence interval with each answer. Name the
coefficients too, and one call returns both tables:

```pycon
>>> both = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", stratum.design_file,
...                          ["malat_load=250ng", "malat_load=125ng"], contrasts=stratum.contrasts_file)
>>> labels = dict(zip(design.columns["contrast_id"], design.columns["label"]))
>>> for cid in both.contrasts.tested_contrasts:
...     print(cid, labels[cid], both.contrasts.weights[cid])
c1 malat_load=125ng vs malat_load=500ng [0, 1, 0]
c2 malat_load=250ng vs malat_load=500ng [0, 0, 1]
>>> all(a["t"] == b["t"] and a["bh_adjusted"] == b["bh_adjusted"]
...     for a, b in zip(both.rows("malat_load=250ng"), both.contrasts.rows("c2")))
True
>>> row = next(r for r in both.contrasts.rows("c2")
...            if r["status"] == "fitted" and r["feature"].startswith("MALAT1:"))
>>> row["feature"], row["observed"], round(row["estimate"], 2), round(row["confidence_low"], 2), round(row["confidence_high"], 2)
('MALAT1:AAAAAAAAACUCG', 4, -0.75, -3.68, 2.17)

```

A contrast of one column is that column's test, number for number, so `c2` repeats the
`malat_load=250ng` rows and adds the interval. This oligo was measured in 4 runs, and its 95%
interval runs from -3.7 to +2.2 log2 units: it cannot tell halving from no change. The contrasts
file earns its keep for comparisons that are no single column, such as 250 ng against 125 ng.
Declare no reference level and `design` writes every pair; the
[differential-abundance guide](stats.md#3-contrasts-any-comparison-with-a-confidence-interval)
runs that chain on this series, with the known answers.

## 2. Covariates, units and batch

A covariate is a number in a unit. Say the unit, and how the slope should read: here age in
years, as a slope per decade, centred at 50 years. Batch is a fixed block. The sheet is mzLib's own
test design, and its R output ships beside it:

```pycon
>>> design = pymzlib.diff.design(
...     "design_factor_batch_covariate_samples.tsv",
...     factors=[pymzlib.diff.FactorSpec("condition", "ctrl")],
...     covariates=[pymzlib.diff.CovariateSpec("age", "years", scale=10, scaled_unit="decade", centre=50)],
...     batch="batch")
>>> stratum = design.strata[0]
>>> stratum.column_names
['intercept', 'condition=treated', 'batch=b2', 'age (per decade)']
>>> slope = design.records[1]
>>> slope["label"], slope["covariate_unit"], slope["covariate_scale"]
('age (per decade)', 'decade', 'years/10, centred at 50 years')
>>> import csv
>>> with open("design_factor_batch_covariate_matrix.tsv", newline="") as fh:
...     r_header, *r_rows = list(csv.reader(fh, delimiter="\t"))
>>> r_header
['(Intercept)', 'conditiontreated', 'batchb2', 'I((age - 50)/10)']
>>> all(abs(ours - float(theirs)) < 1e-12
...     for row, r_row in zip(stratum.matrix, r_rows) for ours, theirs in zip(row, r_row))
True

```

Every value equals R's `model.matrix(~ condition + batch + I((age - 50)/10))`. The bridge's tests
hold all fourteen of mzLib's reference designs to R the same way.

**Leave a field out and mzLib chooses.** `CovariateSpec("age", "years")` sends no scale, scaled unit
or centre, so the slope is per year and the centre is the stratum's mean; `covariate_centres`
reports the centre used. Batch enters the model only when you name its column, and drops out with
a warning in a stratum that has only one batch.

## 3. When the samples cannot answer

Here every sample of condition A ran in batch 1 and every sample of B in batch 2:

```pycon
>>> design = pymzlib.diff.design("design_confounded_samples.tsv",
...                              factors=[pymzlib.diff.FactorSpec("condition", "A")], batch="batch")
>>> stratum = design.strata[0]
>>> stratum.column_names, stratum.rank, stratum.aliased_terms
(['intercept', 'condition=B', 'batch=2'], 2, ['condition', 'batch'])
>>> design.records[0]["label"], design.records[0]["estimability"], design.records[0]["fitted_weights"]
('condition=B vs condition=A', 'not_estimable:confounded', None)

```

**No analysis of these samples can separate the condition from the batch**, and mzLib says so
instead of returning a number. The verdicts are:

- `estimable`;
- `not_estimable:confounded`: the comparison lies along a direction two or more terms share, as
  here;
- `not_estimable:rank_deficient`: it needs a column no sample informs, such as an interaction cell
  with no samples.

The design file still holds a design that can be fitted: `fitted_column_names` leaves out the later
of two aliased columns, because a fit refuses a design without full rank. A contrast that is not
estimable has no fitted weights, so it never reaches the contrasts file.

## 4. Strata: one design per tissue

When samples are different kinds of thing - liver and muscle - a comparison is made within each
kind. Name the stratum column, and each stratum gets its own design:

```pycon
>>> design = pymzlib.diff.design("stratified_samples.tsv",
...                              factors=[pymzlib.diff.FactorSpec("condition", "A")],
...                              strata=["organism part"])
>>> [(s.stratum, s.sample_ids) for s in design.strata]
[('organism part=kidney', ['s7', 's8']), ('organism part=liver', ['s1', 's3', 's5']), ('organism part=muscle', ['s2', 's4', 's6'])]
>>> design.stratum("organism part=kidney").not_run
[NotRunContrast(label='condition=B vs condition=A', reason='not_in_stratum', id='c1')]
>>> list(zip(design.columns["stratum"], design.columns["contrast_id"], design.columns["label"]))
[('organism part=liver', 'c1', 'condition=B vs condition=A'), ('organism part=muscle', 'c1', 'condition=B vs condition=A')]

```

Kidney has only condition A, so it cannot run `c1`, and `not_run` says so under that id. **A
contrast id is numbered once over the whole analysis** (mzLib's `AnalysisDesign.ByStratum`), so `c1`
is the same comparison in liver and in muscle, and results from each stratum's fit stack by id. A
stratum that lacks one side of a comparison skips its id rather than renumbering. With `out=`, stratum *i* writes `design_i.tsv`, numbered in the order of
`strata`, because a stratum key can hold characters a file name cannot.

## 5. From a peptide table to observations

`observations` reads the peptide table FlashLFQ (`QuantifiedPeptides.tsv`) or MetaMorpheus
(`AllQuantifiedPeptides.tsv`) wrote, and needs the runs: which biological sample each file belongs
to and, if it applies, its fraction and technical replicate. **Grouping comes only from what you
give here**; mzLib never reads it out of a file name. This MetaMorpheus table is mzLib's test data:
three peptides in 18 runs, each run its own sample:

```pycon
>>> with open("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv") as fh:
...     columns = fh.readline().rstrip("\n").split("\t")
>>> files = [c[len("Intensity_"):].removesuffix("-calib") for c in columns if c.startswith("Intensity_")]
>>> files[:2]
['QE-002106_GM1_a', 'QE-002107_GM1_b']
>>> table = pymzlib.diff.observations(
...     "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", [(f, f) for f in files], basis="mbr_kept",
...     protein_groups="MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
>>> table.peptide_count, table.sample_count, table.row_count
(3, 18, 54)
>>> table.state_counts
{'quantified': 8, 'mbr_transferred': 1, 'ambiguous_peak': 0, 'identified_not_quantified': 0, 'not_detected': 45}

```

MetaMorpheus quantifies calibrated files, so its columns end in `-calib`; give the file name without
it, and mzLib matches either. Each row is one peptide in one sample:

```pycon
>>> row = next(r for r in table.records if r["state"] == "mbr_transferred")
>>> row["base_sequence"], row["sample"], round(row["log2_intensity"], 2), row["unique_protein_group"]
('CACASHVAK', 'QE-002118_GM6_a', 17.14, 'Q16643')
>>> next(r for r in table.records if r["log2_intensity"] is None)["state"]
'not_detected'

```

**`basis` decides whether match-between-runs values count**, and mzLib has no default: it is your
scientific choice. Under `"msms_only"` the same cell has no value, and still says why:

```pycon
>>> only = pymzlib.diff.observations(
...     "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", [(f, f) for f in files], basis="msms_only",
...     protein_groups="MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
>>> same = next(r for r in only.records if r["sample"] == "QE-002118_GM6_a" and r["base_sequence"] == "CACASHVAK")
>>> same["log2_intensity"], same["state"]
(None, 'mbr_transferred')

```

The protein-group table adds each group's q-value and flags. It is how a contaminant is recognised,
and mzLib warns about groups the peptides name that it does not describe:

```pycon
>>> [(g.name, g.genes) for g in table.protein_groups if g.is_contaminant or g.is_decoy]
[('P02769', 'ALB'), ('DECOY_P62750', 'RPL23A')]
>>> table.warnings
['2 protein group(s) named by peptides are not in the protein-group table, so their q-value and decoy/contaminant status are unknown: Q16643, Q9UI26.']

```

What each `state` means:

- `quantified`: a value from MS/MS-identified peaks only.
- `mbr_transferred`: the value includes a match-between-runs peak. Under `"msms_only"` it is left
  out, and `log2_intensity` is `None`.
- `ambiguous_peak`: a peak was found but withheld as ambiguous. No value.
- `identified_not_quantified`: identified by MS/MS, but no peak could be quantified. No value.
- `not_detected`: no value.

A sample whose runs disagree takes the most informative state among them. With fractions,
mzLib sums a technical replicate's fractions on the linear scale, then averages the sample's
replicates in log2; give each run as `(file, sample, fraction, technical_replicate)`, or an
`ObservationRun`.

## 6. Normalization, and the assumption it makes

`observations` always measures how far each sample sits from the others, and subtracts it only
when you ask. mzLib's rule: over the peptides with a value in every sample (decoys, contaminants
and entrapments excluded), each sample's shift is the median of its log2 value minus the peptide's
mean. Two runs of FlashLFQ's E. coli spike-in vignette, which their names say were spiked at 1x and
2x:

```pycon
>>> runs = [("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"),
...         ("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x")]
>>> spike = pymzlib.diff.observations("EcoliSpikeIn_QuantifiedPeptides.tsv", runs, basis="msms_only")
>>> spike.normalization.setting, spike.normalization.reference_set_size
('none', 86)
>>> spike.normalization.warnings
['86 reference peptide(s) have a value in every sample and 86 in at least half, under the floor of 100, so no sample is shifted (none).']

```

**Below mzLib's floor of 100 shared peptides, nothing is shifted**, because a median of a few
peptides can move a whole sample by itself. You can lower the floor with `min_peptides`. Here is
what that does on this table:

```pycon
>>> shifted = pymzlib.diff.observations("EcoliSpikeIn_QuantifiedPeptides.tsv", runs, basis="msms_only",
...                                     min_peptides=50, normalize=True,
...                                     responses_out="spike_in_responses.tsv")
>>> shifted.normalization.setting
'shared_peptide_median'
>>> {sample: round(shift, 2) for sample, shift in shifted.normalization.shifts.items()}
{'spike_1x': -0.47, 'spike_2x': 0.47}
>>> pairs = [(a["log2_intensity"], b["log2_intensity"])
...          for a, b in zip(spike.records[0::2], spike.records[1::2])
...          if a["log2_intensity"] is not None and b["log2_intensity"] is not None]
>>> len(pairs), round(statistics.median(b - a for a, b in pairs), 2)
(86, 0.94)

```

**The normalization removed the spike-in.** This slice of the vignette holds almost only E. coli
peptides (241 of 243), so nearly every peptide was doubled. The 86 measured in both runs sit 0.94
log2 units higher in the 2x run, close to the doubling the names promise, and the shifts of -0.47
and +0.47 take all of it away. Median normalization assumes most peptides do not change between
samples. When most of them do - a table of the spiked organism alone, a loading series, a
pull-down - that assumption is false, and normalizing deletes the result. Read the shifts before
you subtract them.

`responses_out=` also wrote the table's own peptide-by-sample shape, which
[`stats.fit`](stats.md) reads as its responses:

```pycon
>>> shifted.responses_written.row_count, shifted.responses_written.sample_count
(243, 2)

```

Its rows are peptides by full sequence, its columns the samples, and an empty cell means no value.
For more than one stratum, narrow `observations` to one stratum's samples with `samples=` (a
design file from section 4 will do): the shifts are measured over every sample of the table, so
measure each stratum on its own.

## What these functions will not do

- **No protein roll-up and no row choice.** mzLib does not yet turn the observation table into a
  protein-by-sample table, choose which rows a test uses (decoys, contaminants and q-value cuts
  included), or assemble fits into its differential result table. pyMzLib does not do it either.
  Every peptide is a row; filter with `protein_groups`, `unique_protein_group` and
  `Observations.protein_groups`.
- **Only the default contrasts.** `design` lists each level against its reference (every pair when
  there is none) and each covariate's slope. mzLib's custom contrasts have no wire form yet.
- **The peptide table does not record FlashLFQ's design**, so mzLib cannot check that your runs
  group the files as FlashLFQ did. Give the runs the design FlashLFQ quantified with.
- **A table written in a comma-decimal locale does not read correctly.** FlashLFQ writes in the
  machine's culture, and mzLib reads with a decimal point.
- **pyMzLib does no arithmetic.** Every number on this page came from mzLib; pyMzLib parses the
  bridge's answer into the classes documented in the [API reference](../reference.md#pymzlibdiff).

## See also

- [diff wire verbs](../reference/diff.md): every parameter, field, unit, error and caveat, as all
  three bindings publish them.
- [Differential abundance](stats.md): the moderated test these files feed.
- [FlashLFQ](flashlfq.md) and [Readers](readers.md): where peptide tables come from.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
methods it reproduces:

--8<-- "docs/reference/_generated/cite.diff.md"
