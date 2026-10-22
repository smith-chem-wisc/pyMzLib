# Differential abundance

**The problem.** Most proteomics differential-abundance analyses end in limma: a linear model per
protein, with each protein's noisy variance shrunk toward what all the proteins together say
(empirical-Bayes moderation, Smyth 2004). From Python that has meant installing R and driving it
through `rpy2`, or reimplementing the method and hoping it agrees.

`pymzlib.stats` is that method computed by mzLib, in C#, with nothing to install. It agrees with
limma to better than one part in 10^8, and the first section of this page lets you check that
yourself rather than take it on trust.

Every `>>>` example on this page is executed in CI, against recordings of the real bridge, on the
files pyMzLib ships in `pkg/python/tests/fixtures/stats/`.

## Which function do I want?

| You want to know | Call | mzLib computes it with | Method |
|---|---|---|---|
| Which features change with a coefficient, and how surely? | [`stats.fit`][pymzlib.stats.fit] | `LinearModel.Fit`, then `EmpiricalBayes.Moderate` | limma `lmFit` + `eBayes(legacy = TRUE)` |
| Which features differ between two levels (any contrast), and within what interval? | [`stats.fit`][pymzlib.stats.fit] with `contrasts=` | `EmpiricalBayes.ModerateContrasts` | limma `contrasts.fit` + `eBayes` + `topTable(confint = 0.95)` |
| Which of my p-values survive a false discovery rate? | [`stats.adjust`][pymzlib.stats.adjust] | `MultipleTesting.BenjaminiHochberg` | Benjamini and Hochberg 1995 |
| What is one feature's effect, pooled across studies? | [`stats.meta`][pymzlib.stats.meta] | `RandomEffectsMeta.Pool` | DerSimonian and Laird 1986, as metafor `rma(method = "DL")` |
| What do several studies' p-values for a feature say together? | [`stats.combine_p`][pymzlib.stats.combine_p] | `PValueCombination.Fisher` / `.Stouffer`, with `.OneSided` for signed tests | Fisher 1925; Stouffer 1949; as R's `pchisq` / `qnorm` |
| Which features rise or fall with a per-sample trait (age, dose, load)? | [`stats.correlate`][pymzlib.stats.correlate] | `SpearmanCorrelation.Correlate` / `.PartialCorrelate` | Spearman 1904, as R's `cor.test`; partial as ppcor `pcor.test` |
| Does a feature's detection (or any 0/1 outcome) follow a coefficient? | [`stats.fit_logistic`][pymzlib.stats.fit_logistic] | `LogisticRegression.Fit` | logistic regression by IRLS, as R's `glm(family = binomial)` |
| Which features change, when the samples come in groups (days, donors, studies)? | [`stats.fit_mixed`][pymzlib.stats.fit_mixed] | `MixedModel.Fit` | random-intercept mixed model, as `nlme::lme` with one random intercept per group |
| Does my search really control its FDR? | [`stats.entrapment_fdp`][pymzlib.stats.entrapment_fdp] | `EntrapmentFdp.Sweep` | entrapment FDP estimates, Wen et al. 2025, as FDRBench |

`fit` already adjusts its own p-values, so you need `adjust` only for p-values from somewhere else,
and `meta` only when you have one estimate per study and want them combined. `combine_p` is for
when you have only p-values (and the direction of each effect), and `correlate` for a trait measured
on a scale rather than in groups. `fit_logistic` models a yes/no outcome, such as whether a protein
was detected; `fit_mixed` is `fit` for samples that share a baseline within groups; and
`entrapment_fdp` checks a search's identifications rather than its quantities.

## 1. Check it against limma

mzLib's test suite carries a frozen limma run: 400 features, 12 samples with missing values, a
design of intercept, age in decades and sex, and limma 3.68.5's `eBayes` output for the age
coefficient. pyMzLib ships the input and limma's output. Fit the input:

```pycon
>>> import csv
>>> fit = pymzlib.stats.fit("limma_reference_responses.tsv", "limma_reference_design.tsv",
...                         ["age_decades"])
>>> fit.feature_count, fit.sample_count, fit.coefficient_names
(400, 12, ['Intercept', 'age_decades', 'sex'])

```

and compare every feature with what limma said:

```pycon
>>> with open("limma_ebayes_notrend.tsv", newline="") as fh:
...     limma = list(csv.DictReader(fh, delimiter="\t"))
>>> def worst(ours, column):
...     return max(abs(a - float(row[column])) / abs(float(row[column]))
...                for a, row in zip(ours, limma))
>>> for ours, theirs in [("t", "t_age"), ("p_value", "p_age"), ("bh_adjusted", "bh_age"),
...                      ("posterior_variance", "s2_post"), ("df_total", "df_total")]:
...     print(f"{ours:>18}  worst relative difference from limma < 1e-8: {worst(fit.columns[ours], theirs) < 1e-8}")
                 t  worst relative difference from limma < 1e-8: True
           p_value  worst relative difference from limma < 1e-8: True
       bh_adjusted  worst relative difference from limma < 1e-8: True
posterior_variance  worst relative difference from limma < 1e-8: True
          df_total  worst relative difference from limma < 1e-8: True

```

That is the whole claim, run. mzLib's own tests hold the same agreement for the trended prior
(`trend=True`) and for metafor's DerSimonian-Laird; the bridge's tests repeat both through the
wire, so a column mix-up between mzLib and Python would fail them.

!!! note "Which limma: the *legacy* estimator"
    This matches `eBayes(legacy = TRUE)`. When features have different residual degrees of freedom,
    limma 3.61 and later default to a different prior estimator. Omitting missing values per feature
    is exactly what makes residual df differ, so on most label-free data **default** limma would
    report a slightly different prior and different moderated statistics. The fit tells you when
    your input is such a case:

    ```pycon
    >>> fit.residual_df_differ
    True

    ```

    The reference above is such a case too, and it matches because limma was run with
    `legacy = TRUE`. For such data mzLib also has an estimator of its own, built for unequal
    degrees of freedom: `estimator="marginal_likelihood"`, at the end of
    [section 2](#unequal-degrees-of-freedom-the-marginal-likelihood-prior).

## 2. A real experiment with a known answer

mzLib's RNA test data (mzLib #1388) contains a dilution series, which makes it a rare thing: real
mass-spectrometry data where the right answer is known in advance. Across 27 runs, the lncRNA
**MALAT1** was loaded at 500, 250 or 125 ng, while the luciferase transcript **FLuc** stayed at
500 ng and a 20-mer spike at 50 ng. Halving MALAT1 should read as a log2 fold change of -1, and
quartering it as -2. FLuc and the 20-mer should not change.

The two input files were derived from mzLib's `MetaMorpheus_RNA_AllQuantifiedOligos.tsv` by a
rule that is written once, in the bridge's tests (`MalatDilutionFixture_IsLog2OfMzLibsOligoTable`),
and that fails CI if the committed files ever stop matching it:

1. keep the 27 runs of the series, dropping runs their acquirer labelled `dontuse`;
2. take log2 of each oligo's FlashLFQ intensity, leaving a cell blank where FlashLFQ detected
   nothing;
3. in each run, subtract the median log2 intensity of the oligos that map to FLuc alone. Raw run
   medians in this table span about 7 log2 units, so without this the experiment measures how
   much was injected, not how much MALAT1 there was.

The responses file has one row per oligo, named `<transcript>:<sequence>`, and one column per run.
The design says which runs had less MALAT1. It has an intercept, so each coefficient is the change
from the 500 ng runs:

```pycon
>>> with open("malat_dilution_design.tsv") as fh:
...     header, *runs = [line.rstrip("\n").split("\t") for line in fh]
>>> header
['run', 'intercept', 'malat_250ng', 'malat_125ng']
>>> from collections import Counter
>>> Counter(tuple(run[2:]) for run in runs)
Counter({('1', '0'): 12, ('0', '1'): 8, ('0', '0'): 7})

```

Test both coefficients in one call. The model is fitted once and moderated per coefficient:

```pycon
>>> fit = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
...                         ["malat_250ng", "malat_125ng"])
>>> fit.status_counts
{'fitted': 212, 'too_few_observations': 214, 'rank_deficient': 66}

```

**Most oligos could not be fitted, and the fit says why rather than dropping them.** An oligo seen
in two or three runs has no residual variance to estimate (`too_few_observations`). One seen only in
the 500 ng runs cannot have a 250 ng effect (`rank_deficient`). Both are rows in the table with
`None` statistics, and neither enters the Benjamini-Hochberg family, so they cannot dilute or
inflate it.

Now the question the experiment was designed to answer. Take the median log2 fold change per
transcript:

```pycon
>>> import statistics
>>> def median_change(coefficient, transcript):
...     changes = [r["estimate"] for r in fit.rows(coefficient)
...                if r["status"] == "fitted" and r["feature"].split(":")[0] == transcript]
...     return len(changes), round(statistics.median(changes), 2)
>>> for transcript in ["MALAT1", "FLuc"]:
...     print(transcript, median_change("malat_250ng", transcript), median_change("malat_125ng", transcript))
MALAT1 (116, -0.89) (116, -1.01)
FLuc (91, -0.22) (91, -0.05)

```

**Half the MALAT1 reads as -0.89 against a known -1. A quarter reads as -1.01 against a known -2.**
The first is close; the second is compressed by about half. This table alone cannot say why - the
measurement may compress large changes, or the quarter-load runs may differ in some other way - but
the statistics did not create it: the estimates are ordinary least squares, and moderation changes
only their standard errors. FLuc sits near 0 largely by construction, since step
3 normalised every run to it.

Which individual oligos pass a 5% false discovery rate?

```pycon
>>> for coefficient in fit.tested_coefficients:
...     hits = [r for r in fit.rows(coefficient) if r["bh_adjusted"] is not None and r["bh_adjusted"] < 0.05]
...     print(coefficient, Counter(r["feature"].split(":")[0] for r in hits))
malat_250ng Counter({'MALAT1': 2, 'FLuc': 1})
malat_125ng Counter({'MALAT1': 4, '20mer2': 2})

```

**Few, and not all of them true.** With 27 runs and most oligos missing from many of them, the
moderated test calls very little - that is what moderation is for. And two of the quarter-load
calls are 20-mer oligos, which were loaded at a constant 50 ng. A 5% false discovery rate keeps the
expected share of wrong calls at or below one in twenty, *on average*; it does not promise none; a known-answer experiment is
how you see it happen. Look at the data behind any call you intend to build on.

### Unequal degrees of freedom: the marginal-likelihood prior

This experiment is the case the note in section 1 warns about: oligos are missing from different
runs, so they have different residual degrees of freedom, and the default prior is limma's legacy
estimator:

```pycon
>>> fit.residual_df_differ, fit.prior.estimator, round(fit.prior.df, 1)
(True, 'moments_legacy', 19.6)

```

mzLib has a second estimator for exactly this case, `"marginal_likelihood"` (mzLib #1368). It fits
the prior by maximum likelihood with **each oligo's own degrees of freedom**, so an oligo seen in
five runs counts as five runs' worth of evidence about the variance, not as an average one. Ask for
it by name:

```pycon
>>> ml = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
...                        ["malat_250ng", "malat_125ng"], estimator="marginal_likelihood")
>>> ml.prior.estimator, round(ml.prior.df, 1)
('marginal_likelihood', 11.0)
>>> for coefficient in ml.tested_coefficients:
...     hits = [r for r in ml.rows(coefficient) if r["bh_adjusted"] is not None and r["bh_adjusted"] < 0.05]
...     print(coefficient, Counter(r["feature"].split(":")[0] for r in hits))
malat_250ng Counter({'MALAT1': 1})
malat_125ng Counter({'MALAT1': 4, '20mer2': 2})

```

**The prior df falls from 19.6 to 11.0**: once each oligo's own degrees of freedom count, the data
support less shrinkage toward the common variance. The estimates are the same least-squares numbers;
only the standard errors move. At half load the FLuc call is gone, which is right, since FLuc did not
change, and so is one MALAT1 call. The two 20-mer calls at quarter load stay.

What to know before you choose it:

- **It is mzLib's own estimator, not limma's.** Default limma (3.61 and later) handles unequal
  degrees of freedom with a different method (`fitFDistUnequalDF1`). mzLib's tests hold this one
  close to it, not equal, so expect small differences from default limma.
- **mzLib keeps `"moments_legacy"` the default**, so that `fit` reproduces `eBayes(legacy = TRUE)`
  exactly. Leaving `estimator` out always means mzLib's default, and `fit.prior.estimator` says
  which one you got.
- **A prior df of about 0.05 is a floor, not a finding.** mzLib searches d0 from e^-3 upward, so a
  reported value there means the likelihood wanted less still.
- One prior serves the whole call: every coefficient and every contrast is moderated with it.

## 3. Contrasts: any comparison, with a confidence interval

A design's coefficients answer one kind of question: each level against the baseline. The
questions a scientist asks are often other comparisons: 500 ng against 250 ng when the baseline is
125 ng, or the mean of two groups against a third. Each is a **contrast**, a weighted sum of the
coefficients, and `fit(..., contrasts=...)` tests it with the same moderated t, plus a confidence
interval. That is limma's `contrasts.fit` followed by `eBayes` and `topTable(confint = 0.95)`.

You rarely write the weights yourself. [`diff.design`](diff.md) builds the design from your sample
sheet and writes the contrasts beside it. With no reference level, every pair of levels is a
contrast. Here is the MALAT1 dilution series again, as a sheet of runs and loads:

```pycon
>>> design = pymzlib.diff.design("malat_dilution_samples.tsv",
...                              factors=[pymzlib.diff.FactorSpec("malat_load")], out="malat_pairs")
>>> stratum = design.strata[0]
>>> stratum.column_names
['intercept', 'malat_load=250ng', 'malat_load=500ng']
>>> stratum.warnings
["No reference level was declared for factor 'malat_load', so every pair of its levels is compared (GR-14). Levels: 125ng, 250ng, 500ng."]
>>> labels = dict(zip(design.columns["contrast_id"], design.columns["label"]))
>>> labels
{'c1': 'malat_load=250ng vs malat_load=125ng', 'c2': 'malat_load=500ng vs malat_load=125ng', 'c3': 'malat_load=500ng vs malat_load=250ng'}

```

mzLib took the first level, 125 ng, as the coding baseline, and wrote `design_1.tsv` and
`contrasts_1.tsv`. Fit the design and test every contrast in the file. No coefficient is named, so
the coefficient table is empty and every answer is in `contrasts`:

```pycon
>>> pairs = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", stratum.design_file,
...                           contrasts=stratum.contrasts_file)
>>> pairs.contrasts.tested_contrasts, pairs.contrasts.confidence_level
(['c1', 'c2', 'c3'], 0.95)
>>> pairs.contrasts.weights
{'c1': [0, 1, 0], 'c2': [0, 0, 1], 'c3': [0, -1, 1]}
>>> pairs.row_count, pairs.contrasts.row_count
(0, 1476)

```

The weights are on `pairs.coefficient_names`. `c1` and `c2` are single coefficients. **`c3` is
not**: it is the 500 ng coefficient minus the 250 ng one, a comparison no column of this design
holds. The answers are known in advance, +1, +2 and +1 in log2, so look at the 116 MALAT1 oligos:

```pycon
>>> known = {"c1": 1, "c2": 2, "c3": 1}
>>> for cid in pairs.contrasts.tested_contrasts:
...     malat = [r for r in pairs.contrasts.rows(cid)
...              if r["status"] == "fitted" and r["feature"].startswith("MALAT1:")]
...     covers = sum(r["confidence_low"] <= known[cid] <= r["confidence_high"] for r in malat)
...     clear_of_zero = sum(not r["confidence_low"] <= 0 <= r["confidence_high"] for r in malat)
...     median = statistics.median(r["estimate"] for r in malat)
...     print(f"{labels[cid]:<37} known {known[cid]:+d}  median {median:+.2f}  "
...           f"interval holds it {covers}/{len(malat)}  excludes 0 {clear_of_zero}/{len(malat)}")
malat_load=250ng vs malat_load=125ng  known +1  median +0.18  interval holds it 103/116  excludes 0 7/116
malat_load=500ng vs malat_load=125ng  known +2  median +1.01  interval holds it 105/116  excludes 0 23/116
malat_load=500ng vs malat_load=250ng  known +1  median +0.89  interval holds it 109/116  excludes 0 16/116

```

What the intervals say:

- **Each oligo's 95% interval is wide.** It holds the known change for 103 to 109 of the 116 oligos,
  and also holds 0 for most of them: one oligo, measured in a handful of runs, cannot tell a
  doubling from no change. The medians across oligos are much more precise than any one interval.
- **The compression from section 2 is here too, and the contrasts locate it.** 500 ng against
  250 ng reads +0.89 against a known +1, but 250 ng against 125 ng reads only +0.18 against +1. The
  shortfall is between the two lower loads.
- **The intervals are per feature and not adjusted.** Of 116 correct 95% intervals, about 6 would
  miss by chance alone; 7 to 13 miss here. For a decision across features use `bh_adjusted`, which
  each contrast has over its own fitted features.

Each row of `pairs.contrasts` carries what a coefficient row does, plus `confidence_low` and
`confidence_high`, in response units per contrast unit: a log2 fold change here. As a table:
`pandas.DataFrame(pairs.contrasts.columns)`. One contrast's rows: `pairs.contrasts.rows("c3")`.

**Your own contrasts.** A contrasts file is a header of a contrast id column and then exactly the
design's coefficient columns, in any order, then one row per contrast: its id and one weight per
coefficient. To compare the mean of 250 ng and 125 ng against 500 ng on the design above:

```text
contrast           intercept  malat_load=250ng  malat_load=500ng
low_mean_vs_500    0          0.5               -1
```

(125 ng is the baseline, so its coefficient is 0 and needs no weight of its own.) A column the
design lacks, a design column the file lacks, a weight that is not a number, or a row of all zeros
is refused.

The traps:

- **Contrast ids are numbered within each stratum**, so `c1` can mean different comparisons in two
  strata. Read the labels from `diff.design`, as above.
- **With missing values the standard error is exact per feature**, computed on the samples that
  feature was observed in. limma's `contrasts.fit` approximates it there, so on data with missing
  values these numbers differ from limma's by design. With complete data they agree; the bridge's
  tests hold them to limma 3.68.5 to 1e-8 on mzLib's reference data, with and without missing
  values (against the exact per-feature value for the latter).
- **A contrast of one coefficient is that coefficient's test**, number for number. Asking for both
  in one call costs nothing: the model is fitted once and the prior is shared.
- `confidence=0.9` (a fraction) changes the level; leave it out for mzLib's default, and
  `pairs.contrasts.confidence_level` reports the level used.

## 4. Your own data

`fit` reads two tab-separated files. A **responses** table, one row per feature and one column per
sample, holding the values exactly as they should be modelled:

```text
protein   ctrl_1  ctrl_2  ctrl_3  drug_1  drug_2  drug_3
P02768    24.1    24.3    23.9    25.0    25.2            <- blank: not measured
P68871    20.8    21.0    20.7    20.9    21.1    20.8
```

and a **design**, one row per sample, matched to the responses' columns by name:

```text
sample  intercept  drug
ctrl_1  1          0
ctrl_2  1          0
...
drug_3  1          1
```

Three things decide whether the answer means what you think:

- **Log-transform first.** The model is additive in whatever you pass. Pass log2 intensities and
  `estimate` is a log2 fold change.
- **A 0 is an observation.** mzLib, and many quantification tables, write 0 for "not measured".
  Blank those cells (or write `NA`), or they are fitted as a real zero. The fit counts them in
  `zero_count` and adds a caveat when there are any.
- **The design is used as written.** Include the intercept column yourself. A comparison that is
  not one column is a contrast: see [section 3](#3-contrasts-any-comparison-with-a-confidence-interval).
  [`diff.design`](diff.md) writes both files from a sample sheet.

From pandas, `df.to_csv("responses.tsv", sep="\t")` writes the responses file, with the index as
the feature id.

**Noisier at low abundance?** Pass `trend=True` to let the prior variance follow each feature's
average response, as limma's `trend = TRUE` does. `spline_basis` sets the flexibility of that
curve, and `fit.prior.spline_basis_count` reports what mzLib used.

## 5. Adjusting p-values from anywhere else

`adjust` is Benjamini-Hochberg on any list of p-values. Mark a feature that was not tested with
`None` (or NaN) rather than 1: it then stays out of the family, which is the honest m.

```pycon
>>> adjusted = pymzlib.stats.adjust([0.0002, 0.004, 0.019, None, 0.031, 0.2, float("nan"), 0.74])
>>> adjusted.tested_count
6
>>> [None if v is None else round(v, 4) for v in adjusted.bh_adjusted]
[0.0012, 0.012, 0.038, None, 0.0465, 0.24, None, 0.74]

```

## 6. Pooling across studies

`meta` combines one estimate per study into a random-effects estimate per feature: the shape of a
question like "across every deposit that measured it, does this protein change with age?". Each row
also reports the two things a reader of a pooled result should ask about first: how many studies
point the same way, and how far the answer moves if any one study is dropped.

Here are metafor's own reference cases, which mzLib reproduces to 1e-8:

```pycon
>>> with open("metafor_dl_inputs.tsv", newline="") as fh:
...     studies = [(r["case"], float(r["yi"]), float(r["sei"]))
...                for r in csv.DictReader(fh, delimiter="\t")]
>>> pooled = pymzlib.stats.meta(studies)
>>> for r in pooled.records:
...     print(f'{r["feature"]:>13} k={r["studies"]}  estimate {r["estimate"]:+.3f} '
...           f'[{r["confidence_low"]:+.3f}, {r["confidence_high"]:+.3f}]  I2={r["i_squared"]:.2f}  '
...           f'agree {r["direction_agree"]}/{r["studies"]}  leave-one-out {r["leave_one_out_max_delta"]:.3f}')
  homogeneous k=5  estimate +0.310 [+0.220, +0.400]  I2=0.00  agree 5/5  leave-one-out 0.014
heterogeneous k=6  estimate +0.251 [+0.001, +0.502]  I2=0.82  agree 5/6  leave-one-out 0.095
  two_studies k=2  estimate -0.200 [-0.680, +0.280]  I2=0.48  agree 1/2  leave-one-out 0.300
 random_eight k=8  estimate +0.107 [-0.135, +0.348]  I2=0.81  agree 6/8  leave-one-out 0.075

```

The `heterogeneous` case is the one to learn from: its interval barely excludes zero, most of the
variation (I² = 0.82) is between studies rather than within them, and one study disagrees in sign.
A pooled estimate like that summarises studies that do not agree.

## 7. Combining p-values across studies

Sometimes all you have from each study is a p-value: a supplementary table, a reanalysis that used
another test, or effect sizes on scales too different to pool with `meta`. `combine_p` turns several
independent p-values for the same feature into one, by one of two rules mzLib implements:

| `method=` | Statistic | What drives it |
|---|---|---|
| `"fisher"` | X = -2 sum ln p, against chi-square on 2k degrees of freedom | the smallest p-values: one very strong study is enough |
| `"stouffer"` | Z = sum w z / sqrt(sum w²), against the standard normal | every study, each by its weight (for example the square root of its sample size) |

There is no default: mzLib prefers neither, so you choose. mzLib's own reference cases show the
difference. In `mixed`, one p-value of 0.001 sits among four that disagree:

```pycon
>>> with open("combination_inputs.tsv", newline="") as fh:
...     rows = list(csv.DictReader(fh, delimiter="\t"))
>>> [(float(r["p"]), float(r["weight"])) for r in rows if r["case"] == "mixed"]
[(0.001, 1.5), (0.999, 2.5), (0.5, 3.5), (0.05, 4.5), (0.95, 5.5)]
>>> fisher = pymzlib.stats.combine_p([(r["case"], float(r["p"])) for r in rows], method="fisher")
>>> stouffer = pymzlib.stats.combine_p([(r["case"], float(r["p"]), float(r["weight"])) for r in rows],
...                                    method="stouffer")
>>> round(fisher.p_values["mixed"], 3), round(stouffer.p_values["mixed"], 3)
(0.019, 0.713)

```

Fisher says the five studies together are significant, because one of them is. Weighted Stouffer
says they are not, because the heaviest study (p = 0.95) points the other way. Both agree with R's
`pchisq` and `qnorm` to 1e-8 (the [`combine_p`][pymzlib.stats.combine_p] example checks it).

### The trap: two-sided p-values have no direction

Both rules combine **one-sided** p-values. A two-sided p-value - which is what `fit` reports - says
how surprising an effect is, not which way it went. Combine those and two studies that disagree
count as two that agree. Here are three proteins, each tested in three cohorts of 16, 25 and 9
people (illustrative numbers, written here, not measured):

```pycon
>>> cohorts = [
...     # (protein, two-sided p, weight = sqrt(cohort size), log2 fold change)
...     ("up_everywhere", 0.003, 4, 0.8), ("up_everywhere", 0.02, 5, 0.5), ("up_everywhere", 0.3, 3, 0.3),
...     ("flips_sign", 0.001, 4, 1.2), ("flips_sign", 0.004, 5, -0.9), ("flips_sign", 0.4, 3, 0.2),
...     ("no_effect", 0.6, 4, -0.1), ("no_effect", 0.45, 5, 0.1), ("no_effect", 0.85, 3, -0.05),
... ]
>>> ignored = pymzlib.stats.combine_p([c[:3] for c in cohorts], method="stouffer")
>>> signed = pymzlib.stats.combine_p(cohorts, method="stouffer")
>>> for protein in signed.p_values:
...     print(f"{protein:>13}  direction ignored {ignored.p_values[protein]:<8.2g}  signed {signed.p_values[protein]:.2g}")
up_everywhere  direction ignored 0.00062   signed 8.4e-05
   flips_sign  direction ignored 9.5e-05   signed 0.43
    no_effect  direction ignored 0.69      signed 0.44

```

**Ignoring direction, `flips_sign` is the strongest result of the three.** It rose sharply in one
cohort and fell sharply in another. Give each test its effect (the fourth value; only its sign is
used) and mzLib first turns each two-sided p into the one-sided p-value for "effect > 0"
(`PValueCombination.OneSided`: p/2 for a positive effect, 1 - p/2 for a negative one). Opposite
effects then cancel, and `flips_sign` comes out at 0.43. The result says so:

```pycon
>>> signed.signed, signed.weighted, signed.method
(True, True, 'stouffer')

```

What the columns mean:

- `studies` is how many tests were combined, and `omitted` how many were not: a `None` p-value (a
  test that could not be run) or, signed, an effect of 0 or `None` (no direction).
- `statistic` is X (Fisher) or Z (Stouffer). `None` means it was infinite: a p-value of exactly 0
  makes `p_value` 0, and in Stouffer a p-value of exactly 1 makes it 1.
- `p_value` is the **upper tail**: it tests one direction. Signed, that is effect > 0; a Stouffer
  `p_value` near 1 means a consistently *negative* effect. mzLib computes no two-sided combined
  p-value.

Before you trust a combined p-value:

- **The studies must be independent.** Two analyses of the same samples are one study.
- **It is not adjusted across features.** Send `p_values` to `adjust` for that.
- **It is not an effect size.** If every study reports an estimate and a standard error on one
  scale, `meta` pools those and tells you how much the studies disagree; prefer it.

## 8. Correlating features with a trait

When the samples differ by a number rather than a group - age, dose, the amount loaded - the
question is which features rise or fall with it. `correlate` answers it with Spearman's rank
correlation (Spearman 1904): for each feature, mzLib ranks the feature's values and the trait over
the samples where both are known, correlates the ranks, and tests rho = 0, two-sided. Ranks make it
indifferent to a log transform and to one extreme value, and it finds any relationship that only
ever rises (or only ever falls), not only a straight line.

It reads `fit`'s responses table and a **trait** table: one row per sample, the sample name, then
one value column, whose header names the trait. The dilution series from section 2 has one ready,
derived from its design by a rule held in the bridge's tests:

```pycon
>>> with open("malat_dilution_load_ng.tsv") as fh:
...     header, *runs = [line.rstrip("\n").split("\t") for line in fh]
>>> header, Counter(load for _, load in runs)
(['run', 'malat_ng'], Counter({'250': 12, '125': 8, '500': 7}))
>>> result = pymzlib.stats.correlate("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_load_ng.tsv")
>>> result.feature_count, result.sample_count, result.method_counts
(492, 27, {'exact': 4, 'asymptotic': 308, 'not_estimable': 180})

```

Per transcript, with the unadjusted p-value below 0.01:

```pycon
>>> for transcript in ["MALAT1", "FLuc", "20mer2"]:
...     rows = [r for r in result.records
...             if r["rho"] is not None and r["feature"].split(":")[0] == transcript]
...     rho = statistics.median(r["rho"] for r in rows)
...     print(f"{transcript:>6}  {len(rows):>3} oligos  median rho {rho:+.2f}  "
...           f"p < 0.01: {sum(r['p_value'] < 0.01 for r in rows)}")
MALAT1  194 oligos  median rho +0.43  p < 0.01: 9
  FLuc  111 oligos  median rho +0.00  p < 0.01: 2
20mer2    5 oligos  median rho +0.49  p < 0.01: 2

```

MALAT1's oligos rise with its load, and FLuc's, the spike every run was normalised to, do not. The
20-mer was loaded at a constant 50 ng, yet its five oligos lean positive too, like the 20-mer
calls `fit` made in section 2. A known-answer experiment is how you see that: the measurement, not
the statistics, carries it.

What each row says:

- `rho` runs from -1 to 1; ties take their average rank.
- `p_value` is two-sided against rho = 0. Its `method` says how it was computed: `"exact"` (the
  full permutation distribution, for n of 9 or fewer with no ties, as R's `cor.test`),
  `"asymptotic"` (t on n - 2 degrees of freedom, as `cor.test(exact = FALSE)`), or
  `"not_estimable"`, where `rho` and `p_value` are `None`.
- `n` is how many samples the correlation used: those with both a response and a trait.
- `has_ties` says whether either was tied. Here the trait has only three values, so any oligo seen
  in more than three runs has ties, and its p-value is asymptotic.

**Small n cannot be significant.** The four exact rows each used three runs:

```pycon
>>> [(r["n"], r["rho"], round(r["p_value"], 3)) for r in result.records if r["method"] == "exact"]
[(3, -1, 0.333), (3, -0.5, 1), (3, 0.5, 1), (3, 0.5, 1)]

```

Three samples have six orderings, so even a perfect rank correlation has p = 1/3. And the 180
`not_estimable` oligos were seen in fewer than three runs, or only in runs of one load, where the
trait is constant.

### Partial correlation: holding a nuisance variable fixed

If something else also varies by sample - batch, run order, sex - and could produce the
correlation by itself, give it in a **given** table of the same shape, with one column per
covariate. mzLib then correlates what is left of the feature's ranks and the trait's ranks after
each is regressed on the covariates' ranks (`SpearmanCorrelation.PartialCorrelate`, as
`ppcor::pcor.test(method = "spearman")`). On the limma reference input from section 1, age given
sex:

```pycon
>>> plain = pymzlib.stats.correlate("limma_reference_responses.tsv", "limma_reference_age.tsv")
>>> partial = pymzlib.stats.correlate("limma_reference_responses.tsv", "limma_reference_age.tsv",
...                                   given="limma_reference_sex.tsv")
>>> partial.partial, partial.covariate_names
(True, ['sex'])
>>> plain.method_counts, partial.method_counts
({'exact': 22, 'asymptotic': 378, 'not_estimable': 0}, {'exact': 0, 'asymptotic': 400, 'not_estimable': 0})
>>> for p, q in zip(plain.records[:5], partial.records[:5]):
...     print(f'{p["feature"]}  n={p["n"]}  rho {p["rho"]:+.2f} -> {q["rho"]:+.2f}   p {p["p_value"]:.2g} -> {q["p_value"]:.2g}')
g1  n=12  rho +0.38 -> +0.37   p 0.23 -> 0.26
g2  n=11  rho +0.93 -> +0.97   p 4e-05 -> 4.2e-06
g3  n=11  rho +0.65 -> +0.64   p 0.029 -> 0.046
g4  n=11  rho +0.12 -> +0.12   p 0.73 -> 0.74
g5  n=11  rho +0.43 -> +0.37   p 0.19 -> 0.29

```

The partial test is always asymptotic, on n - 2 - k degrees of freedom for k covariates, so the 22
features that had exact p-values lose them. mzLib drops a covariate that the others already carry
and counts only the ones kept; it does not report that count.

Before you trust a correlation:

- **Partial out only nuisance variables.** A covariate that follows the trait takes the trait's
  signal with it: correlate dose given a batch that ran one dose per day, and nothing is left.
- **An asymptotic p-value of exactly 0 means rho is exactly -1 or 1.** It is the t-test's limit,
  not a measured probability. Eleven of the 400 reference features have one:

  ```pycon
  >>> [(r["feature"], r["n"], r["rho"], r["p_value"]) for r in plain.records if r["p_value"] == 0][:3]
  [('g6', 12, -1, 0), ('g34', 12, 1, 0), ('g50', 11, 1, 0)]

  ```

- **For 10 to 1289 samples without ties, R's default `cor.test` uses an Edgeworth series** that
  mzLib does not, so p-values there differ slightly from R's default output (not from
  `cor.test(exact = FALSE)`).
- **It is not adjusted across features.** Send the column to `adjust`; a `None` (not estimable)
  stays out of the family:

  ```python title="Not run: no recording holds the adjustment of these p-values"
  adjusted = pymzlib.stats.adjust(result.columns["p_value"])
  ```

## 9. Is a feature's detection related to the condition?

Missing values in label-free data are rarely random. A protein near the detection limit is missed
more often where there is less of it, and an analysis that ignores the missing cells is then
biased. Before you choose how to handle them, ask the data: **does each feature's chance of being
detected depend on the condition?** That is a logistic regression per feature, on a table of 1
(measured) and 0 (not measured). `fit_logistic` runs one per feature with mzLib's
`LogisticRegression.Fit`, which is R's `glm(family = binomial)` per feature, and tests each
coefficient with a Wald z.

The dilution series from section 2 has a detection table ready. It is derived from the responses
file by one rule, held in the bridge's tests: 1 where the responses file has a value, 0 where it is
blank. The design is section 2's:

```pycon
>>> with open("malat_dilution_detected.tsv") as fh:
...     header, *oligos = [line.rstrip("\n").split("\t") for line in fh]
>>> sorted({cell for oligo in oligos for cell in oligo[1:]})
['0', '1']
>>> detected = pymzlib.stats.fit_logistic("malat_dilution_detected.tsv", "malat_dilution_design.tsv",
...                                       ["malat_250ng", "malat_125ng"])
>>> detected.status_counts
{'fitted': 191, 'too_few_observations': 0, 'rank_deficient': 0, 'separated': 301, 'not_converged': 0}

```

Every oligo has an outcome in every run, because 0 is an outcome here, not a missing value. So none
is `too_few_observations`. But most are **separated**, and that word needs reading before the
rest. Per transcript, at a quarter of the MALAT1 load:

```pycon
>>> for transcript in ["MALAT1", "FLuc", "20mer2"]:
...     rows = [r for r in detected.rows("malat_125ng") if r["feature"].split(":")[0] == transcript]
...     fitted = [r for r in rows if r["status"] == "fitted"]
...     line = f"{transcript:>6}  {len(rows):>3} oligos  {len(fitted):>3} fitted"
...     if fitted:
...         change = statistics.median(r["estimate"] for r in fitted)
...         line += f"  median log-odds {change:+.2f}  p < 0.05: {sum(r['p_value'] < 0.05 for r in fitted)}"
...     print(line)
MALAT1  345 oligos  113 fitted  median log-odds -0.18  p < 0.05: 4
  FLuc  138 oligos   77 fitted  median log-odds +0.41  p < 0.05: 3
20mer2    6 oligos    0 fitted

```

`estimate` is a change in the **log-odds** of detection per unit of the coefficient. For a 0/1
column it is the log of the odds ratio against the 500 ng runs. Quartering MALAT1 moves its
oligos' odds of detection by a median factor of exp(-0.18), about 0.84, and only 4 of 113 oligos
pass an unadjusted p of 0.05. In this table, MALAT1's missing values follow its load only weakly:
they are mostly about how close each oligo sits to the detection limit in every run. FLuc's odds
rise a little as MALAT1 falls, plausibly from less competition for the same signal.

What each row says:

- `z` is the Wald statistic, `estimate / standard_error`. `p_value` is two-sided and **not
  adjusted across features**: send one coefficient's column to `adjust`.
- `deviance` is -2 times the log-likelihood at the estimate, the same in every coefficient's block.
- `events` counts the runs with a 1, and `observed` the runs with a 0 or a 1.
- `iterations` is how many reweighted least-squares steps mzLib took. It is 0 when no fit was
  attempted.

**A separated feature has no estimate, and that is not the same as no effect.** When the design
splits a feature's 0s from its 1s, the likelihood keeps rising as a coefficient runs to infinity,
so no maximum-likelihood estimate exists. R's `glm` prints a huge coefficient and a warning there;
mzLib reports `separated` and `None`. With three load groups, an oligo seen in a handful of runs is
almost always missing from every run of one group, and so is separated:

```pycon
>>> separated = [r for r in detected.rows("malat_125ng") if r["status"] == "separated"]
>>> Counter(r["events"] for r in separated).most_common(4)
[(0, 64), (2, 55), (1, 54), (3, 37)]

```

Sixty-four were never detected at all, and most of the rest in one to three runs. Here separation
means "too sparse to say". In a cleaner experiment the same status can mark the strongest result
in the table: a protein detected in every treated sample and in no control. Read `events` and
`observed` before you discard one. Firth's penalized likelihood (R's `logistf`), which gives a
finite estimate there, is not implemented.

**Check it against R.** pyMzLib ships mzLib's own reference input (simulated binary features) and
R's `glm` output for it:

```pycon
>>> def relative_worst(pairs):
...     return max(abs(ours - float(theirs)) / abs(float(theirs)) for ours, theirs in pairs)
>>> reference = pymzlib.stats.fit_logistic("glm_reference_responses.tsv", "glm_reference_design.tsv",
...                                        ["Intercept", "x", "grp"])
>>> with open("glm_reference_fit.tsv", newline="") as fh:
...     glm = list(csv.DictReader(fh, delimiter="\t"))
>>> rows = list(zip(reference.rows("grp"), glm))
>>> relative_worst((ours["estimate"], r["b2"]) for ours, r in rows if r["status"] == "fitted") < 1e-8
True
>>> [ours["status"] for ours, r in rows if r["status"] == "separated"]
['separated', 'separated', 'separated', 'separated']

```

The estimates agree with R to better than one part in 10^8, and the standard errors and p-values
to 1e-6, because `summary.glm` takes its weights from the step before the estimate. Where R's own
fitted probabilities came within 1e-10 of 0 or 1, mzLib reports `separated` too.

Before you trust a logistic fit:

- **Wald tests are weakest exactly where they matter.** Near separation, and in small samples, the
  standard error grows faster than the estimate, and z shrinks. No likelihood-ratio test is
  offered.
- **A 0 is an outcome.** Blank only a cell you genuinely do not know.
- **`max_iterations` and `tolerance` are mzLib's when you leave them out,** and the result then
  echoes `None` for them. mzLib's default tolerance is tighter than `glm`'s default `epsilon` of
  1e-8, so default `glm` output can differ in the late digits.

## 10. When samples come in groups

Runs acquired on different days, replicates from one donor, studies pooled into one analysis: each
group can sit at its own baseline. `fit` treats every sample as independent, so a day-to-day shift
either inflates its residual variance or, when the days are unbalanced, leaks into the effect. A
**mixed model** gives every group its own baseline, a random intercept with variance tau², and
estimates each coefficient *within* groups. `fit_mixed` runs one per feature with mzLib's
`MixedModel.Fit`, which is `nlme::lme(y ~ ..., random = ~ 1 | group)` per feature.

The dilution series' runs were acquired on five days, and each run's name begins with its date. A
**groups** table says which group each sample is in: the sample name, then one label column whose
header names the grouping. This one is derived from the run names by a rule held in the bridge's
tests:

```pycon
>>> with open("malat_dilution_day.tsv") as fh:
...     header, *runs = [line.rstrip("\n").split("\t") for line in fh]
>>> header, Counter(day for _, day in runs)
(['run', 'day'], Counter({'2026_09_24': 8, '2026_09_18': 6, '2026_09_27': 6, '2026_09_26': 4, '2026_09_25': 3}))
>>> mixed = pymzlib.stats.fit_mixed("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
...                                 "malat_dilution_day.tsv", ["malat_250ng", "malat_125ng"])
>>> mixed.estimator, mixed.group_count
('reml', 5)
>>> mixed.status_counts
{'fitted': 210, 'too_few_observations': 214, 'rank_deficient': 66, 'too_few_groups': 2, 'variance_ratio_at_limit': 0, 'not_converged': 0}

```

The design is section 2's, unchanged. Do not add a column per day: the random intercept replaces
it. `estimator` is `"reml"` because none was given and that is mzLib's default. Two oligos were
seen on only one day, or on no day twice, so their day variance cannot be told from their residual
variance: `too_few_groups`.

Compare the same MALAT1 oligos under both models:

```pycon
>>> linear = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
...                            ["malat_250ng", "malat_125ng"])
>>> def malat(result, coefficient):
...     return {r["feature"]: r for r in result.rows(coefficient)
...             if r["status"] == "fitted" and r["feature"].startswith("MALAT1:")}
>>> malat(mixed, "malat_250ng").keys() == malat(linear, "malat_250ng").keys()
True
>>> for coefficient in ["malat_250ng", "malat_125ng"]:
...     for name, result in [("fit", linear), ("fit_mixed", mixed)]:
...         rows = malat(result, coefficient).values()
...         estimate = statistics.median(r["estimate"] for r in rows)
...         error = statistics.median(r["standard_error"] for r in rows)
...         print(f"{coefficient}  {name:>9}  {len(rows)} oligos  median estimate {estimate:+.2f}  "
...               f"median standard error {error:.2f}")
malat_250ng        fit  116 oligos  median estimate -0.89  median standard error 1.00
malat_250ng  fit_mixed  116 oligos  median estimate -1.05  median standard error 0.71
malat_125ng        fit  116 oligos  median estimate -1.01  median standard error 1.09
malat_125ng  fit_mixed  116 oligos  median estimate -1.09  median standard error 0.70

```

**Halving MALAT1 now reads as -1.05 against a known -1, with smaller standard errors.** Part of
what section 2 counted as noise was the days. The quarter load still reads about -1 against a
known -2, so that compression is not a day effect. How much the days matter is in the variance
columns, which are the same in every coefficient's block:

```pycon
>>> by_day = malat(mixed, "malat_250ng").values()
>>> round(statistics.median(r["group_variance"] for r in by_day), 2)
0.32
>>> round(statistics.median(r["residual_variance"] for r in by_day), 2)
0.55

```

A typical MALAT1 oligo's day-to-day variance is more than half its run-to-run variance within a
day.

What each row says:

- `estimate` is the effect within groups, in response units per coefficient unit: a log2 fold
  change here. `standard_error` and `t` are as in `fit`.
- `df` is the t test's denominator degrees of freedom, by nlme's **containment** rule: observed runs
  minus observed days minus the coefficients that vary within a day. It is small for an oligo seen
  in few runs:

  ```pycon
  >>> Counter(r["df"] for r in by_day).most_common(4)
  [(1, 26), (2, 18), (None, 12), (4, 11)]

  ```

  A `None` df means the coefficient was constant within every day the oligo was seen on, so it can
  only be tested between days, and too few days are left to test it. `t` is still reported, and
  `p_value` is `None`.
- `residual_variance` is sigma², within groups, and `group_variance` is tau², between them, both in
  squared response units. A `group_variance` of 0 is the boundary: no day effect at all, and the fit
  is ordinary least squares. nlme cannot reach 0 and prints a tiny positive number there instead.
- `log_likelihood` is restricted under `"reml"`, so compare it only between models with the same
  fixed effects. Pass `estimator="ml"` to compare designs.

**Check it against nlme.** mzLib's reference input has 29 simulated samples in seven unequal
groups, with a group-level covariate `wp` (constant within each group) beside `age`, which varies
within groups:

```pycon
>>> reference = pymzlib.stats.fit_mixed("lmm_reference_responses.tsv", "lmm_reference_design.tsv",
...                                     "lmm_reference_groups.tsv", ["Intercept", "age", "wp"])
>>> with open("lmm_reference_fit.tsv", newline="") as fh:
...     nlme = [r for r in csv.DictReader(fh, delimiter="\t") if r["method"] == "REML"]
>>> relative_worst((ours["estimate"], r["b1"]) for ours, r in zip(reference.rows("age"), nlme)) < 1e-5
True
>>> all(ours["df"] == float(r["df2"]) for ours, r in zip(reference.rows("wp"), nlme))
True

```

Estimates, standard errors, variances and log-likelihoods agree with nlme to 1e-5 relative, which
is nlme's own optimizer precision, and every df is identical.

Before you trust a mixed model:

- **The df are nlme's containment df,** not Satterthwaite's or Kenward-Roger's, so `lmerTest` and
  `pbkrtest` report different p-values.
- **One random intercept, nothing more.** No random slopes, and no nested or crossed groupings.
- **A feature the fixed effects fit exactly is reported fitted,** with a residual variance near
  1e-16 and a t in the millions: mzLib's guard catches only an exactly zero residual. A constant
  feature is the usual case. Check `residual_variance` before you believe a tiny p-value.
- **Not adjusted across features.** As everywhere on this page, send `p_value` to `adjust`.

## 11. Does a search really control its FDR?

A search engine reports q-values, and a 1% threshold should mean that about 1% of the discoveries
are false. **Entrapment** checks that claim. Peptides that cannot be in the sample, from another
species or made from the real ones, are added to the search database, so every entrapment the
search reports is a known false discovery. `entrapment_fdp` turns those counts into estimates of
the false discovery proportion (FDP) at each q-value threshold, with mzLib's
`EntrapmentFdp.Sweep`, following Wen et al. 2025 (the method FDRBench implements).

It reads two tables, one row per unit you are estimating for (precursor, peptide or protein group),
decoys removed, each with its best q-value:

- **targets**: the columns `score` and `q_value`, and no others;
- **entrapments**: `score` and `q_value`, then `partner_score` and `partner_q_value`, the score and
  q-value of the one target it was made from when the search identified that target, both blank
  otherwise. An optional `foreign` column marks entrapments with no partner by construction.

pyMzLib ships mzLib's hand-built test search: seven targets and four entrapments, one entrapment
for each case the paired estimate tells apart.

```pycon
>>> with open("entrapment_reference_entrapments.tsv", newline="") as fh:
...     for row in csv.DictReader(fh, delimiter="\t"):
...         print(dict(row))
{'score': '85', 'q_value': '0.003', 'partner_score': '100', 'partner_q_value': '0.001'}
{'score': '75', 'q_value': '0.006', 'partner_score': '70', 'partner_q_value': '0.008'}
{'score': '65', 'q_value': '0.01', 'partner_score': '', 'partner_q_value': ''}
{'score': '55', 'q_value': '0.03', 'partner_score': '20', 'partner_q_value': '0.5'}
>>> fdp = pymzlib.stats.entrapment_fdp("entrapment_reference_targets.tsv",
...                                    "entrapment_reference_entrapments.tsv", ratio=1.0)
>>> fdp.thresholds_given, fdp.threshold_count, fdp.thresholds[:3], fdp.thresholds[-1]
(False, 100, [0.001, 0.002, 0.003], 0.1)
>>> for q in [0.003, 0.006, 0.01, 0.05]:
...     r = fdp.at(q)
...     print(f"q {q:<5}  targets {r['target_count']}  entrapments {r['entrapment_count']}  "
...           f"lower {r['lower_bound']:.3f}  combined {r['combined']:.3f}  paired {r['paired']:.3f}")
q 0.003  targets 2  entrapments 1  lower 0.333  combined 0.667  paired 0.333
q 0.006  targets 3  entrapments 2  lower 0.400  combined 0.800  paired 0.600
q 0.01   targets 4  entrapments 3  lower 0.429  combined 0.857  paired 0.857
q 0.05   targets 6  entrapments 4  lower 0.400  combined 0.800  paired 0.800

```

Without `thresholds`, mzLib sweeps its default thresholds and the result echoes them. The three
estimates answer different questions:

- **`lower_bound`** is N_E / (N_T + N_E). Above the FDR you claimed, the search **fails** to
  control it. Below, it proves nothing.
- **`combined`** is N_E (1 + 1/r) / (N_T + N_E), with r the entrapment-to-target database size
  ratio. It is an upper bound: at or below the claimed FDR, it is **evidence of control**.
- **`paired`** is a tighter upper bound for 1:1 pairing. It asks of each discovered entrapment
  whether its own target was discovered too, and ranked ahead of it. At q = 0.003 it equals the
  lower bound, because the one entrapment's partner was found first.

None is capped at 1. This toy search fails badly: at 1%, three of its seven discoveries are
entrapments.

**`paired` needs a ratio of exactly 1.** With any other ratio, or when every entrapment is
foreign, mzLib does not report it:

```pycon
>>> doubled = pymzlib.stats.entrapment_fdp("entrapment_reference_targets.tsv",
...                                        "entrapment_reference_entrapments.tsv", ratio=2.0,
...                                        thresholds=[0.01, 0.05])
>>> [(r["q_value_threshold"], round(r["combined"], 3), r["paired"]) for r in doubled.records]
[(0.01, 0.643, None), (0.05, 0.6, None)]

```

For 1:1 pairing pass 1 itself, not a ratio computed from the database sizes after the generator's
exclusions: 0.9998 is not 1, and `paired` disappears.

Before you trust an entrapment estimate:

- **One row per unit, at its best q-value, decoys removed.** mzLib does not deduplicate. Estimate
  at the level your FDR claim is made at.
- **The partner columns are copies.** mzLib does not check them against the targets file, so build
  both files from the same table.
- **Exclusions go in opposite directions.** An entrapment that is really a target peptide comes
  out of the entrapments file. A target that cannot be traced to one target stays in the targets
  file, but is nobody's partner (mzLib's `EntrapmentFdp.ClassifyExclusion`).
- **The "sample" estimator N_E / (r N_T) is deliberately absent:** it is invalid in both
  directions. The k-fold paired estimator, for k entrapments per target, is not implemented.

## What these functions will not do

- **No imputation.** A missing value is left out of its feature's fit. If your pipeline imputes,
  do it before writing the responses file, and remember the imputed values then count as data.
- **Not limma's current default** when residual df differ (see the note in section 1):
  `estimator="marginal_likelihood"` is mzLib's own estimator for that case, close to it but not
  equal. Not `robust = TRUE`, the B-statistic or observation weights.
- **`bh_adjusted` is not a q-value from target-decoy search.** It controls the false discovery rate
  among the features you tested, which is a different question from identification FDR.
- **`meta` uses a normal reference** (no Hartung-Knapp adjustment), so its intervals do not widen
  for two or three studies. Treat such pooled estimates with care.
- **`combine_p` gives no two-sided combined p-value and no effect size**, and `correlate` does not
  report how many covariates mzLib kept. Neither adjusts across features: that is `adjust`.
- **`fit_logistic` has no Firth correction and no likelihood-ratio test**, so a separated feature
  gets no estimate. **`fit_mixed` has one random intercept** and nlme's containment df: no random
  slopes, no Satterthwaite or Kenward-Roger. Neither adjusts across features.
- **`entrapment_fdp` estimates; it does not build the entrapment database** or decide which
  generator exclusions apply. It has no "sample" estimator and no k-fold paired estimator.
- **pyMzLib does no arithmetic.** Every number on this page came from mzLib's
  `StatisticalModels`; pyMzLib parses the bridge's answer into the classes documented in the
  [API reference](../reference.md#pymzlibstats).

## See also

- [stats wire verbs](../reference/stats.md): every parameter, field, unit, error and caveat, as all
  three bindings publish them.
- [Differential design and observations](diff.md): the design and contrasts files `fit` reads,
  from a sample sheet.
- [FlashLFQ](flashlfq.md) and [Readers](readers.md): where intensity tables come from.
- [Protein databases](proteins.md): which proteins mzLib marks as entrapment (`is_entrapment`),
  for building the `entrapment_fdp` tables.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
methods it reproduces:

--8<-- "docs/reference/_generated/cite.stats.md"
