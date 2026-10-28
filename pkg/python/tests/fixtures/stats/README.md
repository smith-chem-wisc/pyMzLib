# Inputs for the `stats` examples and tests

Small tables the `pymzlib.stats` docstrings, the [differential-abundance guide](../../../../../docs/guides/stats.md)
and the tests read. None is hand-written. Every file is held to its source by a C# test in
`pkg/bridge.tests/StatisticsTests.cs`, which fails if a committed file stops being exactly what
its rule produces. To rewrite them from a fresh mzLib checkout:

```text
PYMZLIB_REGENERATE_STATS_FIXTURES=1 dotnet test pkg/bridge.tests --filter "FullyQualifiedName~StatisticsTests"
```

| File | What it is | Held by |
|---|---|---|
| `malat_dilution_log2_vs_fluc.tsv` | One row per RNA oligo, one column per run, from mzLib's `Test/FileReadingTests/ExternalFileTypes/MetaMorpheus_RNA_AllQuantifiedOligos.tsv` (mzLib #1388). The 27 runs of the MALAT1 500/250/125 ng series at constant 500 ng FLuc and 50 ng 20mer-2 (147 min gradient), runs labelled `dontuse` dropped. Cells are log2(intensity) minus that run's median log2 intensity of the FLuc-only oligos; blank where FlashLFQ detected nothing. Feature ids are `<protein groups>:<full sequence>`. | `MalatDilutionFixture_IsLog2OfMzLibsOligoTable` |
| `malat_dilution_design.tsv` | One row per run: `intercept`, and 0/1 indicators `malat_250ng`, `malat_125ng`. The 500 ng runs are the baseline. | same |
| `limma_reference_responses.tsv` | mzLib's `Test/StatisticalModels/ReferenceData/limma_responses.tsv` (400 features, 12 samples, NaN = missing), with feature ids `g1`.. added. Cells verbatim. | `ReferenceFixtures_AreMzLibsLimmaAndMetaforData` |
| `limma_reference_design.tsv` | mzLib's `limma_design.tsv` (Intercept, age_decades, sex), with the sample names added. | same |
| `limma_ebayes_notrend.tsv` | limma 3.68.5's `eBayes(legacy = TRUE)` output for `age_decades` on that input, verbatim. | same |
| `metafor_dl_inputs.tsv`, `metafor_dl_results.tsv` | metafor 5.2.1's DerSimonian-Laird reference cases and results, verbatim. Constructed cases, not proteomics data. | same |
| `limma_metafor_PROVENANCE.txt` | mzLib's record of the R, limma and metafor versions and the seed that produced the reference files. | same |
| `combination_inputs.tsv`, `combination_results.tsv` | mzLib's p-value combination reference cases (case, p, weight) and R 4.6.1's Fisher (`pchisq`), Stouffer and weighted Stouffer (`qnorm`) results for them, verbatim. Constructed p-values, not proteomics data. | `CombinationFixtures_AreMzLibsReferenceCases` |
| `combination_PROVENANCE.txt` | mzLib's `PROVENANCE_glm_lmm.txt`, verbatim: the R version, the nlme version and seed behind the combination, glm and lme references, and the note that the combination files were regenerated on 2026-09-28. | same |
| `malat_dilution_load_ng.tsv` | One row per run of the dilution series: `malat_ng`, the MALAT1 load in ng (250 where `malat_250ng` is 1, 125 where `malat_125ng` is 1, else 500), from `malat_dilution_design.tsv`. The trait `stats.correlate` reads. | `CorrelateFixtures_AreColumnsOfTheCommittedDesigns` |
| `limma_reference_age.tsv`, `limma_reference_sex.tsv` | The `age_decades` and `sex` columns of `limma_reference_design.tsv`, verbatim, as one-column per-sample tables: the trait and the covariate of the partial-correlation example. | same |
| `glm_reference_responses.tsv`, `glm_reference_design.tsv` | mzLib's `glm_responses.tsv` (80 simulated binary features, 40 samples, NaN = missing) and `glm_design.tsv` (Intercept, x, grp), with feature ids `f1`.. and the sample names `V1`.. added. Cells verbatim. | `LogisticFixtures_AreMzLibsGlmReference_AndTheDilutionSeriesDetection` |
| `glm_reference_fit.tsv` | R 4.6.1's `glm(family = binomial, epsilon = 1e-14)` output for that input (status, coefficients, standard errors, p of `x`, deviance), verbatim. | same |
| `malat_dilution_detected.tsv` | The dilution series' detection table: 1 where `malat_dilution_log2_vs_fluc.tsv` has a value (FlashLFQ measured the oligo in that run), 0 where it is blank. What `stats.fit_logistic` reads. | same |
| `lmm_reference_responses.tsv`, `lmm_reference_design.tsv`, `lmm_reference_groups.tsv` | mzLib's `lmm_responses.tsv` (120 simulated features, 29 samples), `lmm_design.tsv` (Intercept, age, wp) and `lmm_groups.tsv` (seven groups of unequal size), with feature ids `f1`.. and the sample names `V1`.. added. Cells verbatim. | `MixedFixtures_AreMzLibsLmeReference_AndTheDilutionSeriesDays` |
| `lmm_reference_fit.tsv` | nlme 3.1.169's `lme` output for that input, REML and ML (coefficients, standard errors, containment df, p-values, variances, log-likelihood), verbatim. | same |
| `malat_dilution_day.tsv` | One row per run of the dilution series: `day`, the acquisition date its name begins with (`2026_09_18`, ...). The groups `stats.fit_mixed` reads. | same |
| `entrapment_reference_targets.tsv`, `entrapment_reference_entrapments.tsv` | mzLib's hand-built entrapment search (`Test/StatisticalModels/EntrapmentFdpTests.cs`: seven targets, four entrapments), read out of the test source, with each entrapment's partner written as its target's score and q-value. Constructed numbers, not a real search. | `EntrapmentFixtures_AreMzLibsReferenceSearch` |
| `entrapment_reference_r1.tsv` | The entrapment project's reference script (`fdp_analysis.py` at 92f3108) output for that search at r = 1, as mzLib's test carries it, to five decimals. | same |
| `lnmm_reference_design.tsv`, `lnmm_reference_groups.tsv` | mzLib's `lnmm_design.tsv` (Intercept, a 0/1 `trait` constant within each animal, a covariate `x`) and `lnmm_groups.tsv` (eight animals of 2 to 4 runs), with the sample names `V1`..`V23` of the counts' header added. Cells verbatim. | `LatentMixedFixtures_AreMzLibsLnmmReference` |
| `lnmm_scatter_counts.tsv` | mzLib's `lnmm_scatter_successes.tsv` and `lnmm_scatter_trials.tsv` (20 simulated features, every sample counted with 20 to 40 trials, with run scatter), joined as `successes/trials` per cell, blank where both are NaN; feature ids `f1`.. added. | same |
| `lnmm_scatter_fit.tsv` | glmmTMB 1.1.15.2's REML and ML fits of that input (an observation-level random intercept for the run scatter), verbatim. | same |
| `lnmm_link_counts.tsv`, `lnmm_link_depth.tsv` | mzLib's `lnmm_counted_link_*.tsv` (16 simulated features, every sample counted, through the link -0.4, 0.6, 0.3) joined as `successes/trials`, and `lnmm_counted_depth.tsv` with the sample names added. | same |
| `lnmm_link_residual_variance.tsv` | sigma² = 2e-7 for each of those features: mzLib's `TinyResidualVariance`, which its glmer comparison gives because glmer fits no run scatter. | same |
| `lnmm_counted_fit.tsv` | lme4 2.0.6's `glmer` fits of the counted inputs (25-point adaptive quadrature; ML, and REML by a Laplace step), verbatim. | same |
| `lnmm_PROVENANCE.txt` | mzLib's `PROVENANCE_lnmm.txt`, verbatim: the R, lme4, glmmTMB and TMB versions and the seed behind every `lnmm_` reference file. | same |
| `lnmm_simulated_responses.tsv`, `lnmm_simulated_counts.tsv`, `lnmm_simulated_design.tsv`, `lnmm_simulated_groups.tsv` | mzLib's measured-and-counted simulation (`LogitNormalMixedModelTests.Simulated`, seeds 1, 2 and 3, one feature each): 16 animals of 4 runs, a true trait effect of 1.0 on the latent scale, every other run measured and the rest counted, 5 trials each, through the link (0.4, 0.6). Values written round-trip. | same |
| `link_cells.tsv` | mzLib's `FitLink` test cells (`LogitNormalMixedModelTests.LinkCells`, the link -0.5, 0.5, 0.3, 4,000 cells, seed 5), values written round-trip. | same |

The recorded bridge payloads that answer these calls in the examples are one folder up:
`stats_fit_malat.json`, `stats_fit_malat_marginal.json` (`estimator="marginal_likelihood"`),
`stats_fit_malat_contrasts.json` (every pair of loads, from `diff design`), `stats_fit_limma.json`,
`stats_adjust.json`, `stats_meta_metafor.json`,
`stats_combine_p_*.json`, `stats_correlate_*.json`, `stats_fit_logistic_*.json`,
`stats_fit_mixed_*.json`, `stats_fit_latent_mixed_*.json`, `stats_fit_link_simulated.json` and
`stats_entrapment_fdp_*.json`.
