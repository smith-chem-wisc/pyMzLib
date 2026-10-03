# Differential abundance

**The problem.** Most proteomics differential-abundance analyses end in limma: a linear model per
protein, with each protein's noisy variance shrunk toward what all the proteins together say
(empirical-Bayes moderation, Smyth 2004). From Python that has meant installing R and driving it
through `rpy2`, or reimplementing the method and hoping it agrees.

`pymzlib.stats` is that method computed by mzLib, in C#, with nothing to install. It agrees with
limma to better than one part in 10^8, and the first section of this page lets you check that
yourself rather than take it on trust.

Every code block on this page is executed in CI, against recordings of the real bridge, on the
files pyMzLib ships in `pkg/python/tests/fixtures/stats/`.

## Which function do I want?

| You want to know | Call | mzLib computes it with | Method |
|---|---|---|---|
| Which features change with a coefficient, and how surely? | [`stats.fit`][pymzlib.stats.fit] | `LinearModel.Fit`, then `EmpiricalBayes.Moderate` | limma `lmFit` + `eBayes(legacy = TRUE)` |
| Which of my p-values survive a false discovery rate? | [`stats.adjust`][pymzlib.stats.adjust] | `MultipleTesting.BenjaminiHochberg` | Benjamini and Hochberg 1995 |
| What is one feature's effect, pooled across studies? | [`stats.meta`][pymzlib.stats.meta] | `RandomEffectsMeta.Pool` | DerSimonian and Laird 1986, as metafor `rma(method = "DL")` |

`fit` already adjusts its own p-values, so you need `adjust` only for p-values from somewhere else,
and `meta` only when you have one estimate per study and want them combined.

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

## What these functions will not do

- **No imputation.** A missing value is left out of its feature's fit. If your pipeline imputes,
  do it before writing the responses file, and remember the imputed values then count as data.
- **Not limma's current default** when residual df differ (see the note in section 1); not
  `robust = TRUE`, contrasts, the B-statistic or observation weights.
- **`bh_adjusted` is not a q-value from target-decoy search.** It controls the false discovery rate
  among the features you tested, which is a different question from identification FDR.
- **`meta` uses a normal reference** (no Hartung-Knapp adjustment), so its intervals do not widen
  for two or three studies. Treat such pooled estimates with care.
- **pyMzLib does no arithmetic.** Every number on this page came from mzLib's
  `StatisticalModels`; pyMzLib parses the bridge's answer into the classes documented in the
  [API reference](../reference.md#pymzlibstats).

## See also

- [stats wire verbs](../reference/stats.md): every parameter, field, unit, error and caveat, as all
  three bindings publish them.
- [FlashLFQ](flashlfq.md) and [Readers](readers.md): where intensity tables come from.

## References

- Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential expression
  in microarray experiments. *Stat Appl Genet Mol Biol* 3:3.
  [doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027)
- Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and
  microarray studies. *Nucleic Acids Res* 43:e47.
  [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful
  approach to multiple testing. *J R Stat Soc B* 57:289.
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. *Control Clin Trials* 7:177.
  [doi:10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456(86)90046-2)
- Viechtbauer W (2010). Conducting meta-analyses in R with the metafor package. *J Stat Softw*
  36(3). [doi:10.18637/jss.v036.i03](https://doi.org/10.18637/jss.v036.i03)
