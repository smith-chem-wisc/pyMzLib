# Inputs for the `diff` examples and tests

Small tables the `pymzlib.diff` docstrings, the [differential design and observations guide](../../../../../docs/guides/diff.md)
and the tests read. None is hand-made data. Every file is held to its source by
`PythonInputs_AreMzLibsFilesOrTheirRule` in `pkg/bridge.tests/DifferentialTests.cs`, which fails if a
committed file stops being exactly what it says it is.

| File | What it is |
|---|---|
| `design_factor_batch_covariate_samples.tsv`, `design_confounded_samples.tsv` | mzLib's `Test/Quantification/Differential/ReferenceData/` sample sheets, verbatim. Constructed designs, not proteomics data. |
| `design_factor_batch_covariate_matrix.tsv` | Base R 4.6.1's `model.matrix` for that sheet (mzLib's `make_design_fixtures.R`), verbatim. |
| `stratified_samples.tsv` | The six samples of mzLib's `AnalysisDesignTests.SamplesAreSplitByStratum`, plus two kidney samples of condition A only, so one stratum cannot run the contrast. |
| `malat_dilution_samples.tsv` | One row per run of the MALAT1 dilution series in `../stats/`: `malat_load` is `250ng` where `malat_250ng` is 1 in `malat_dilution_design.tsv`, `125ng` where `malat_125ng` is 1, else `500ng`. |
| `MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv`, `MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv` | mzLib's `Test/FileReadingTests/ExternalFileTypes/` MetaMorpheus 1.1.11 tables, verbatim: three peptides in 18 calibrated runs, and six protein groups (a contaminant and a decoy among them). |
| `EcoliSpikeIn_QuantifiedPeptides.tsv` | `../flashlfq_vignette_peptides.tsv`, verbatim: the slice of FlashLFQ's E. coli spike-in vignette (two runs, 1x and 2x) that the median-polish parity test reads, under a name mzLib recognises as a peptide table. |

The recorded bridge payloads that answer the examples are one folder up: `diff_design_*.json`,
`diff_observations_*.json`, and `stats_fit_malat_diff_design.json` (stats fit on the design file
`diff design` wrote for the dilution series).
