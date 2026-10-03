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

The recorded bridge payloads that answer these calls in the examples are one folder up:
`stats_fit_malat.json`, `stats_fit_limma.json`, `stats_adjust.json` and `stats_meta_metafor.json`.
