# Inputs for the `ptmqtl` examples and tests

Small tables the `pymzlib.ptmqtl` docstrings, the [PTM occupancy and PTM-QTL guide](../../../../../docs/guides/ptmqtl.md)
and the tests read. None is hand-written. They are mzLib's own hand-computed PtmQtl test cases
(`mzLib/Test/Quantification/PtmQtlTests.cs` at the pin), each on its own protein so that no case
changes another, rendered by `PtmQtlCases` in `pkg/bridge.tests/PtmQtlTests.cs`. The test
`Fixtures_AreMzLibsHandComputedCases` fails if a committed file stops being exactly what those
cases render. To rewrite them:

```text
PYMZLIB_REGENERATE_PTMQTL_FIXTURES=1 dotnet test pkg/bridge.tests --filter "FullyQualifiedName~PtmQtlTests"
```

| File | What it is |
|---|---|
| `cases.tsv` | The single-run cases, one protein each, with mzLib's run names: `RATIO` (OccupancyIsModifiedOverCoveringIntensity, 30 / 200), `CEILING` (SeenOnlyModifiedIsOccupancyOneWithTheCeilingFlag), `STATES` (TheFourStatesAreKeptApart), `SHARED_POSITION` (TwoModificationsAtOnePositionShareTheDenominator), `EXCLUDED`, `PEPTIDE_NTERM`, `NTERM_AFTER_MET`, `NTERM_AFTER_K` and `UNREGISTERED` (ExclusionsAndFilterFollowDefOccPsms and PeptideNTerminalModificationsAndNonMetStartsAreNotProteinNTermini), and `TWO_ON_ONE` (PhysicalPairsCarryCoOccupancy, 0.25 over two runs). |
| `fractions.tsv`, `fractions_samples.tsv` | FractionsOfOneSampleAreSummedBeforeTheRatio (`SUMMED`) and AFormUnquantifiedInOneFractionAndCoveredInAnotherIsAFloorForTheSample (`FLOOR_ACROSS`), with mzLib's map f1, f2 -> A and g1 -> B. |
| `batch.tsv` | RunLevelRemovalTakesOutABatchShift: twelve runs in two batches, `PAIR_S` and `PAIR_T` nearly unrelated within a batch, and four `LEVEL_` sites at the batch's centre. Occupancy in percent of a covering intensity of 100. Each protein has its own peptide, since mzLib calls sites whose covering peptides share a sequence overlapping. |
| `trait.tsv`, `trait_runs.tsv` | mzLib's TraitDesign (TraitEffectIsRecoveredOnTheLogitScale): six replicates at doses 0.25, 1 and 2, two injections each, logit(occupancy) = -2 + 0.5 dose + a replicate offset, over a covering intensity of 1000 (`TRAIT_SITE`). `CENSORED_SITE` is TraitFitUsesOnlyQuantifiedNonCeilingCells: replicate R5 a floor in one injection and a ceiling in the other. The runs table adds the injection covariate of TraitFitReportsCovariateEffects. |
| `trait_occupancy.tsv` | What `ptmqtl occupancy --out` writes from `trait.tsv`, held byte for byte by `TraitOccupancyFixture_IsWhatOccupancyOutWrites`: the input of the `trait_effects` examples. |
| `bounds.tsv` | BoundsEnterTypeAAsTiedRanks: `BOUND_S` and `BOUND_T` with a value in six of ten runs and, in the other four, a covered-zero, a floor, a ceiling, or a count-only cell seen only modified (a bound) or in both forms (no rank). Occupancy in percent of a covering intensity of 100. |
| `depth.tsv`, `depth_psms.tsv` | ARunCovariateIsPartialledOut: twelve runs in two batches; `DEPTH_S` and `DEPTH_T` also rise with each run's target PSMs (`depth_psms.tsv`, the run covariate), and four `DEPTH_LEVEL_` sites sit at the batch's centre. Occupancy in millionths of a covering intensity of 1e6, so every intensity is an exact integer. |

The recorded bridge payloads that answer these calls are one folder up: `ptmqtl_occupancy_*.json`,
`ptmqtl_pairs_*.json` and `ptmqtl_trait_effects*.json`.
