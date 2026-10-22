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
| Which of my p-values survive a false discovery rate? | [`stats.adjust`][pymzlib.stats.adjust] | `MultipleTesting.BenjaminiHochberg` | Benjamini and Hochberg 1995 |
| What is one feature's effect, pooled across studies? | [`stats.meta`][pymzlib.stats.meta] | `RandomEffectsMeta.Pool` | DerSimonian and Laird 1986, as metafor `rma(method = "DL")` |
| What do several studies' p-values for a feature say together? | [`stats.combine_p`][pymzlib.stats.combine_p] | `PValueCombination.Fisher` / `.Stouffer`, with `.OneSided` for signed tests | Fisher 1925; Stouffer 1949; as R's `pchisq` / `qnorm` |
| Which features rise or fall with a per-sample trait (age, dose, load)? | [`stats.correlate`][pymzlib.stats.correlate] | `SpearmanCorrelation.Correlate` / `.PartialCorrelate` | Spearman 1904, as R's `cor.test`; partial as ppcor `pcor.test` |

`fit` already adjusts its own p-values, so you need `adjust` only for p-values from somewhere else,
and `meta` only when you have one estimate per study and want them combined. `combine_p` is for
when you have only p-values (and the direction of each effect), and `correlate` for a trait measured
on a scale rather than in groups.

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
    `legacy = TRUE`.

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

## 3. Your own data

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
- **The design is used as written.** Include the intercept column yourself; write each comparison
  you want as its own column (contrasts are not implemented).

From pandas, `df.to_csv("responses.tsv", sep="\t")` writes the responses file, with the index as
the feature id.

**Noisier at low abundance?** Pass `trend=True` to let the prior variance follow each feature's
average response, as limma's `trend = TRUE` does. `spline_basis` sets the flexibility of that
curve, and `fit.prior.spline_basis_count` reports what mzLib used.

## 4. Adjusting p-values from anywhere else

`adjust` is Benjamini-Hochberg on any list of p-values. Mark a feature that was not tested with
`None` (or NaN) rather than 1: it then stays out of the family, which is the honest m.

```pycon
>>> adjusted = pymzlib.stats.adjust([0.0002, 0.004, 0.019, None, 0.031, 0.2, float("nan"), 0.74])
>>> adjusted.tested_count
6
>>> [None if v is None else round(v, 4) for v in adjusted.bh_adjusted]
[0.0012, 0.012, 0.038, None, 0.0465, 0.24, None, 0.74]

```

## 5. Pooling across studies

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

## 6. Combining p-values across studies

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

## 7. Correlating features with a trait

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

## What these functions will not do

- **No imputation.** A missing value is left out of its feature's fit. If your pipeline imputes,
  do it before writing the responses file, and remember the imputed values then count as data.
- **Not limma's current default** when residual df differ (see the note in section 1); not
  `robust = TRUE`, contrasts, the B-statistic or observation weights.
- **`bh_adjusted` is not a q-value from target-decoy search.** It controls the false discovery rate
  among the features you tested, which is a different question from identification FDR.
- **`meta` uses a normal reference** (no Hartung-Knapp adjustment), so its intervals do not widen
  for two or three studies. Treat such pooled estimates with care.
- **`combine_p` gives no two-sided combined p-value and no effect size**, and `correlate` does not
  report how many covariates mzLib kept. Neither adjusts across features: that is `adjust`.
- **pyMzLib does no arithmetic.** Every number on this page came from mzLib's
  `StatisticalModels`; pyMzLib parses the bridge's answer into the classes documented in the
  [API reference](../reference.md#pymzlibstats).

## See also

- [stats wire verbs](../reference/stats.md): every parameter, field, unit, error and caveat, as all
  three bindings publish them.
- [FlashLFQ](flashlfq.md) and [Readers](readers.md): where intensity tables come from.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
methods it reproduces:

--8<-- "docs/reference/_generated/cite.stats.md"
