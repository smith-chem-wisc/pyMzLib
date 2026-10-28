# PTM site occupancy and PTM-QTL

**The problem.** A modification site is rarely all-or-nothing. In one run 15% of a serine is
phosphorylated, in the next 40%. That fraction, the site's **occupancy**, is what changes with
age, dose or genotype, and what two sites can share. Measuring it from label-free data is easy to
get wrong: a modified form that was identified but not quantified is not a zero, a site seen only
in its modified form is not truly 100%, and two sites read off one peptide are correlated by
construction.

`pymzlib.ptmqtl` is mzLib's `Quantification.PtmQtl` engine, which keeps those cases apart. It
computes occupancy per run from a plain table, finds pairs of modifications, and fits the effect
of a per-run trait on each site.

Every `>>>` example on this page is executed in CI, against recordings of the real bridge, on the
files pyMzLib ships in `pkg/python/tests/fixtures/ptmqtl/`. Those files are mzLib's own
hand-computed test cases, so each number below is one mzLib's tests assert.

## Which function do I want?

| You want to know | Call | mzLib computes it with |
|---|---|---|
| What fraction of each site carries its modification, run by run? | [`ptmqtl.occupancy`][pymzlib.ptmqtl.occupancy] | `SiteOccupancyCalculator.Calculate` |
| The same for a sample measured as several fractions? | [`ptmqtl.occupancy`][pymzlib.ptmqtl.occupancy] with `samples=` | `SiteOccupancyCalculator.Calculate` with a run-to-sample map |
| Which two modifications sit on one molecule? | [`ptmqtl.pairs`][pymzlib.ptmqtl.pairs] with `type="physical"` | `PtmPairEngine.Physical` |
| Which two sites rise and fall together across runs? | [`ptmqtl.pairs`][pymzlib.ptmqtl.pairs] with `type="covarying"` | `PtmPairEngine.CoVarying` (Spearman, BH) |
| The same, counting runs where a site was seen in only one form? | [`ptmqtl.pairs`][pymzlib.ptmqtl.pairs] with `include_bounds=True` | `PtmPairEngine.CoVarying(includeBounds)`, `PairTestStatus` |
| The same, given each run's depth? | [`ptmqtl.pairs`][pymzlib.ptmqtl.pairs] with `run_covariate=` | `PtmPairEngine.CoVarying(runCovariate)` |
| How many peptidoforms covered a site in a run, and how many carried the mark? | [`ptmqtl.occupancy`][pymzlib.ptmqtl.occupancy], columns `covering_peptidoforms`, `modified_peptidoforms` | `SiteRunOccupancy.CoveringPeptidoforms`, `.ModifiedPeptidoforms` |
| Does a site's occupancy change with a per-run trait? | [`ptmqtl.trait_effects`][pymzlib.ptmqtl.trait_effects] | `SiteTraitEffects.Fit` (a mixed model per site) |
| The same, with the effect in log2 odds? | [`ptmqtl.trait_effects`][pymzlib.ptmqtl.trait_effects] with `scale="log2_odds"` | `SiteTraitOptions.Scale`, `OccupancyScale` |

All three read the same input: one table of peptidoforms and their intensities. `occupancy` can
write its result to a file, and `trait_effects` reads that file.

## 1. The input: one row per peptidoform, per run, per protein

Every function reads a tab-separated table with a header row. It is mzLib's own
`PeptidoformObservation`, as columns:

| column | what it holds |
|---|---|
| `run` | the run (raw file) the peptidoform was identified in |
| `full_sequence` | the peptidoform in MetaMorpheus notation: `PEPS[Common Biological:Phosphorylation on S]K` |
| `protein_accession` | one protein it maps to |
| `start_residue`, `end_residue` | the peptide's first and last residue in that protein, counting from 1 |
| `intensity` | its MS/MS-quantified intensity in that run; **blank** when it was identified but not quantified |
| `previous_residue` | optional: the protein residue before the peptide, one letter |

Four rules decide whether the numbers mean anything:

- **One row per (run, peptidoform, protein), never one per PSM.** Two rows for one peptidoform in
  one run would count its intensity twice; mzLib refuses them.
- **MS/MS identifications only.** A match-between-runs transfer is not an identification. Leave it
  out.
- **Blank, not 0, for "identified but not quantified".** FlashLFQ writes 0 for a peptide it did not
  quantify. Blank those cells.
- **A shared peptide is one row per protein.** Each protein then gets its full intensity.
  MetaMorpheus instead splits a shared peptide across protein groups, so its
  `IntensityOccupancy_` columns can differ.

**mzLib does not build this table for you yet.** It has no reader that turns a MetaMorpheus or
FlashLFQ result into these rows, because that step makes choices: which PSMs count, how a
shared peptide is placed. pyMzLib will not make them on mzLib's behalf. You can assemble the
table from two files pyMzLib already reads: the PSMs (for runs, proteins and positions) and the
FlashLFQ peptide table (for intensities). One way, with your own choices marked:

```python title="Not run: needs your MetaMorpheus search output"
import csv
import pymzlib

psms = pymzlib.readers.read_records("AllPSMs.psmtsv")
peptides = pymzlib.readers.read_quantified_peptides("AllQuantifiedPeptides.tsv")

# Intensity per (run, full sequence), MS/MS-detected only. FlashLFQ's 0 means "not quantified".
intensity = {
    (r["sample_label"], r["sequence"]): r["intensity"]
    for r in peptides.records
    if r["detection_type"] == "MSMS" and r["intensity"]
}

def per_protein(cell, count):
    """MetaMorpheus writes one value when every protein shares it, else one per protein."""
    parts = cell.split("|")
    return parts * count if len(parts) == 1 else parts

rows = set()
for p in psms.records:
    # Your choice of PSMs: targets, q-value below 1%, one unambiguous peptidoform.
    if p["decoy_contam_target"] != "T" or not p["q_value"] < 0.01 or p["ambiguity_level"] != "1":
        continue
    accessions = p["protein_accession"].split("|")
    spans = per_protein(p["start_and_end_residues_in_parent_sequence"], len(accessions))  # "[69 to 84]"
    previous = per_protein(p["previous_residue"], len(accessions))
    for accession, span, before in zip(accessions, spans, previous):
        start, end = (int(x) for x in span.strip("[]").split(" to "))
        rows.add((p["file_name_without_extension"], p["full_sequence"], accession, start, end,
                  before if len(before) == 1 else ""))

with open("observations.tsv", "w", newline="") as fh:
    out = csv.writer(fh, delimiter="\t")
    out.writerow(["run", "full_sequence", "protein_accession", "start_residue", "end_residue",
                  "intensity", "previous_residue"])
    for run, sequence, accession, start, end, previous in sorted(rows):
        value = intensity.get((run, sequence))
        out.writerow([run, sequence, accession, start, end, "" if value is None else value,
                      previous])
```

The set keeps one row per (run, peptidoform, protein) however many PSMs the peptidoform has. A
shared peptide becomes one row per protein, each with its own span. Check that the run names of
the two files match: FlashLFQ labels a calibrated file `<run>-calib`.

## 2. Occupancy, run by run

The cases file holds several of mzLib's hand-computed cases, each on its own protein, named for
what it tests.

```pycon
>>> occ = pymzlib.ptmqtl.occupancy("cases.tsv")
>>> occ.observation_count, occ.site_count
(22, 9)
>>> for r in occ.records:
...     if r["protein_accession"] == "RATIO":
...         print(r["site_key"], r["position"], r["fraction"], r["unmodified_quantified"])
RATIO:S13:Common Biological:Phosphorylation on S 13 0.15 True

```

`RATIO` has 30 units of the phosphopeptide, 70 of the same peptide unmodified, and 100 of a
missed-cleavage peptide that also covers the serine: 30 / 200 = 0.15. Every quantified form that
covers the position counts in the denominator, whatever its own sequence.

**Each cell has a state.** Only a `"quantified"` cell has a fraction:

```pycon
>>> for r in occ.records:
...     if r["protein_accession"] == "STATES":
...         print(f'{r["run"]:<12}{r["state"]:<14}{r["fraction"]}')
countOnly   count_only    None
floor       floor         None
notDetected not_detected  None
quant       quantified    0.1

```

| state | what was seen | read it as |
|---|---|---|
| `quantified` | the modified form was quantified | a measured fraction in (0, 1] |
| `floor` | the modified form was identified but not quantified, while forms covering the site were | low, **not zero** |
| `count_only` | the modified form was identified, but nothing covering the site was quantified | present, amount unknown |
| `not_detected` | peptidoforms covering the site were identified, and none carried the modification | a **bound**: below what this run could detect, **not zero** |

**Each cell also counts peptidoforms.** `covering_peptidoforms` is how many distinct peptidoforms
covering the position were identified in the run, quantified or not. `modified_peptidoforms` is how
many of them carried the modification there. They count identifications, not intensity, so they
are the "k of n" of a cell that has no fraction:

```pycon
>>> for r in occ.records:
...     if r["protein_accession"] == "STATES":
...         print(f'{r["run"]:<12}{r["state"]:<14}{r["modified_peptidoforms"]} of {r["covering_peptidoforms"]}')
countOnly   count_only    1 of 1
floor       floor         1 of 2
notDetected not_detected  0 of 1
quant       quantified    1 of 2

```

The `count_only` cell is 1 of 1: the site was seen only modified. Had the unmodified form been
identified too, it would read 1 of 2. That difference decides how
[`include_bounds`](#counting-one-form-runs-as-bounds) ranks the cell.

**A 1 can be a ceiling.** When a site is seen only in its modified form, every quantified form
covering it carries the modification, so the fraction is 1. The unmodified form may be there below
detection. `unmodified_quantified` is `False` then:

```pycon
>>> [(r["fraction"], r["unmodified_quantified"]) for r in occ.records
...  if r["protein_accession"] == "CEILING"]
[(1, False)]

```

**Two modifications at one position share a denominator.** On `SHARED_POSITION`, 20 units carry
phosphorylation, 30 HexNAc and 50 neither. They are two sites, at 0.2 and 0.3.

**Positions count from 1, in the protein.** A protein N-terminal modification sits at the
residue it is on: residue 2 when the initiator methionine was removed.

```pycon
>>> sorted({r["site_key"] for r in occ.records if r["protein_accession"].startswith(("NTERM", "UNREG"))})
['NTERM_AFTER_MET:S2:UniProt:N-acetylserine on S', 'UNREGISTERED:C1:Renamed:Unregistered on C']

```

!!! warning "Do not join this table with `readers.read_occupancy` on position"
    MetaMorpheus's occupancy file, which `readers.read_occupancy` reads, puts the N-terminus at
    position 0. Here residue 1 is position 1.

**What mzLib leaves out.** The cases file also holds rows that make no site, and so no row here:

- `EXCLUDED`: oxidation and carbamidomethylation. mzLib always skips the `Common Variable` and
  `Common Fixed` categories.
- `PEPTIDE_NTERM`: a peptide N-terminal modification (ammonia loss) is not a protein site.
- `NTERM_AFTER_K`: N-acetylation on a peptide starting at residue 2 whose `previous_residue` is
  K. Residue 2 is the protein's N-terminus only after a methionine is removed.

C-terminal modifications are never sites, because the table does not carry the protein's length.
A modification name that mzLib's registry does not know, like `UNREGISTERED`'s, is treated as
protein N-terminal. To keep only some modifications, pass `modifications=` (exact names) or
`categories=` (for example `["Common Biological"]`).

## 3. A sample measured as several fractions

Fractions of one sample are summed before the ratio, as FlashLFQ sums them. Give the run-to-sample
map as a file:

```pycon
>>> occ = pymzlib.ptmqtl.occupancy("fractions.tsv", samples="fractions_samples.tsv")
>>> occ.run_count, occ.sample_count, occ.run_names
(3, 2, ['A', 'B'])
>>> for r in occ.records:
...     print(r["protein_accession"], r["run"], r["state"], r["fraction"])
FLOOR_ACROSS A floor None
SUMMED A quantified 0.2
SUMMED B not_detected None

```

`SUMMED` reads 0.3 in fraction 1 and 0.1 in fraction 2. Sample A is (30 + 10) / (100 + 100) = 0.2,
not the mean of 0.3 and 0.1. On `FLOOR_ACROSS` the phosphopeptide was identified, unquantified, in
one fraction, and the site was covered in the other: a floor for the sample.

**Map only fractions together.** Technical replicates are separate measurements of one sample;
give each its own name, and let [`trait_effects`](#5-does-occupancy-change-with-a-trait) treat
them as replicates of one biological sample.

## 4. Pairs of modifications

### On one molecule

A peptidoform carrying two modifications is direct evidence that both sit on one molecule:

```pycon
>>> found = pymzlib.ptmqtl.pairs("cases.tsv", type="physical")
>>> r = found.records[0]
>>> r["site_a"], r["site_b"]
('TWO_ON_ONE:K14:Common Biological:Acetylation on K', 'TWO_ON_ONE:S13:Common Biological:Phosphorylation on S')
>>> r["statistic"], r["n"], r["statistic_n"], r["p_value"]
(0.25, 2, 1, None)

```

`n` is the runs where the doubly modified form was identified (2). `statistic` is the median
**co-occupancy** over the runs where it was quantified (1): 25 of the 100 units covering both
positions carried both marks. There is no test for physical pairs, so `p_value` and `q` are `None`.

### Across runs

Two sites that rise and fall together across runs are co-varying. mzLib correlates their
occupancies with Spearman's rho, on the runs where both were quantified, and adjusts the p-values
with Benjamini-Hochberg within each family:

```pycon
>>> plain = pymzlib.ptmqtl.pairs("batch.tsv", type="covarying")
>>> plain.covarying_run_count, plain.min_quantified_fraction, plain.exclude_ceiling
(12, 0.7, None)
>>> def show(result, a="PAIR_S", b="PAIR_T"):
...     r = next(r for r in result.records if r["site_a"].startswith(a) and r["site_b"].startswith(b))
...     print(f'rho {r["statistic"]:.2f}  p {r["p_value"]:.3f}  q {r["q"]:.3f}  {r["fdr_family"]}  n {r["n"]}')
>>> show(plain)
rho 0.76  p 0.004  q 0.004  A:inter  n 12

```

`0.7` is mzLib's default, reported because you did not set it: a site enters pairs only if it is
quantified in at least 70% of the runs. `exclude_ceiling` is `None` because you did not choose, so
mzLib decided inside `PtmPairEngine.CoVarying`; by default it leaves ceilings out.

**The trap: a batch correlates everything.** These twelve runs come in two batches, and every site
reads 0.2 higher in the second. That alone gives `PAIR_S` and `PAIR_T` a rho of 0.76. Within a
batch they are nearly unrelated. `remove_run_level=True` correlates the two sites given each run's
**level**, the median occupancy of the other sites in that run:

```pycon
>>> given = pymzlib.ptmqtl.pairs("batch.tsv", type="covarying", remove_run_level=True)
>>> show(given)
rho 0.03  p 0.934  q 0.934  A:inter  n 12
>>> given.run_level_group, given.min_runs_for_run_level, given.min_sites_for_run_level
('none', 10, 3)

```

Given the batch, the pair is unrelated. A pair is tested this way only when it shares at least
`min_runs_for_run_level` runs that have a level, and a run has one only with at least
`min_sites_for_run_level` other sites there. The four `LEVEL_` sites in this file follow only the
batch, so their pairs with each other stay at rho 1, and their pairs with `PAIR_S` and `PAIR_T` are
`not_estimable` once the level, which they set, is removed.

!!! warning "Remove the run level for a batch or loading shift only"
    If the level follows your design - one condition reads higher everywhere - removing it removes
    the effect you are studying, and the test becomes a correlation within conditions.

Sites of different chemistry can have different batch shifts. `run_level_group="category"` or
`"modification"` gives each group its own level; for any other grouping, such as all
phosphorylations together, pass a two-column file (`modification`, `group`) as
`run_level_groups=`.

**What to read in the table.**

- **`overlapping` pairs have no `q`.** When one peptide covers both sites, their occupancies share
  intensities and are correlated by construction, so mzLib reports them but leaves them out of
  the adjustment.
- **`fdr_family`** is the type (`P` physical, `A` co-varying) and whether the two sites are on one
  protein (`intra`) or two (`inter`). `q` is adjusted within it. It is not a target-decoy q-value.
- **`keep_ceiling=True`** keeps runs where a site was seen only modified. Two unrelated sites that
  both read 1 in low-load runs then correlate through detection alone; that is why mzLib leaves
  them out by default.
- **`test_status`** says whether a co-varying pair was tested, and if not, why: `not_estimable`
  (too few complete runs, or a constant site), `too_few_shared_runs` (a partial correlation on
  fewer than `min_runs_for_run_level` runs), or `both_sites_mostly_bounds` (next section). An
  untested pair has no `statistic`, `p_value` or `q`, and is not in its family's adjustment.

### Counting one-form runs as bounds

A site seen in only one form in a run has no fraction, but it is not missing either. A
`not_detected` cell says the occupancy was below what the run could see; a ceiling says it was at
the top. By default mzLib leaves such runs out of the correlation. A site that is often below
detection then never reaches the 70% guard, and its pairs are never tested.

`include_bounds=True` lets those runs in as **tied ranks**, never as invented values. A
`not_detected` or `floor` cell ties below every fraction of its site. A ceiling, or a `count_only`
cell whose covering peptidoforms all carry the mark (k of n with k = n), ties above. A `count_only`
cell where both forms were identified has no rank and stays out.

The example is mzLib's `BoundsEnterTypeAAsTiedRanks`: two sites with a value in six of ten runs,
and a bound in the other four.

```pycon
>>> plain = pymzlib.ptmqtl.pairs("bounds.tsv", type="covarying")
>>> plain.row_count, plain.include_bounds
(0, False)
>>> bounded = pymzlib.ptmqtl.pairs("bounds.tsv", type="covarying", include_bounds=True)
>>> r = bounded.records[0]
>>> r["n"], r["test_status"], r["bounds_included"], f'{r["statistic"]:.3f}'
(9, 'tested', True, '0.004')
>>> bounded.min_untied_share, bounded.exclude_ceiling
(0.5, None)

```

Six of ten is below the guard, so the plain call tests nothing. With bounds, `BOUND_S` has a value
or a bound in all ten runs and `BOUND_T` in nine, and the pair is tested on the nine they share.
The tenth is `BOUND_T`'s `count_only` cell with both forms identified (1 of 2): it has no rank. The
two sites turn out unrelated (rho 0.004).

**Too many ties make the test unreliable.** When both sites are mostly one bound, Spearman's
t test calls far too many pairs: in mzLib's simulation, 2 to 5 times the nominal rate at
p < 0.001 with no relation at all. So a pair is tested only if one of its sites has at least
`min_untied_share` (half) of the shared runs outside its largest tie. Otherwise its `test_status`
is `both_sites_mostly_bounds`, and it has no p-value.

`exclude_ceiling` stays `None`: with bounds, mzLib treats a ceiling as a high bound.
`keep_ceiling=True` sends `False` and gives the same answer.

### Given each run's depth

Shallow runs miss low-abundance forms everywhere at once. With bounds in, that makes covered-zeros
co-occur across unrelated sites, and any site that rises with depth correlates with every other.
`run_covariate=` takes a per-run value, such as the run's target PSMs, and correlates the two sites
**given** it: a partial Spearman correlation, as for the run level.

The file has a header row whose second field names the covariate, then one row per run:

```tsv
run	target_psms
r00	9000
r01	12000
```

The example is mzLib's `ARunCovariateIsPartialledOut`: twelve runs in two batches, and two sites
that also rise a little with each run's depth.

```pycon
>>> def depth_pair(result):
...     r = next(r for r in result.records
...              if r["site_a"].startswith("DEPTH_S") and r["site_b"].startswith("DEPTH_T"))
...     return f'rho {r["statistic"]:.2f}  p {r["p_value"]:.3f}  n {r["n"]}'
>>> depth_pair(pymzlib.ptmqtl.pairs("depth.tsv", type="covarying"))
'rho 0.77  p 0.003  n 12'
>>> given = pymzlib.ptmqtl.pairs("depth.tsv", type="covarying", run_covariate="depth_psms.tsv")
>>> depth_pair(given)
'rho 0.76  p 0.006  n 12'
>>> given.run_covariate_name, given.run_covariate_runs_missing, given.min_runs_for_run_level
('target_psms', [], 10)
>>> both = pymzlib.ptmqtl.pairs("depth.tsv", type="covarying", run_covariate="depth_psms.tsv",
...                             remove_run_level=True)
>>> depth_pair(both)
'rho 0.08  p 0.832  n 12'

```

Depth alone barely moves rho: in these runs most of the shared rise is the batch. Given the batch's
run level as well, the pair is unrelated. A covariate removes only what it explains.

- **Any increasing transform gives the same answer.** The covariate is ranked, so target PSMs and
  log2 target PSMs are the same input.
- **Leave a run out for "unknown".** It is then left out of every pair, and listed in
  `run_covariate_runs_missing`. A blank or `NA` value is not "unknown": mzLib refuses it, and
  pyMzLib raises `UsageError`.
- **The same floor of shared runs applies.** A pair partial on the covariate is tested only on at
  least `min_runs_for_run_level` shared runs; below that it is `too_few_shared_runs`.

## 5. Does occupancy change with a trait?

`trait_effects` fits, for each site, logit(occupancy) against a per-run trait, with a random
intercept for each biological replicate. Repeated injections of one animal or donor are then not
counted as independent evidence.

The example is mzLib's `TraitDesign`: six replicates at three doses, two injections each, built so
that logit(occupancy) = -2 + 0.5 x dose + a replicate offset. First compute occupancy, and write it
to a file:

```pycon
>>> occ = pymzlib.ptmqtl.occupancy("trait.tsv", out="trait_occupancy.tsv")
>>> occ.written.row_count, occ.state_counts
(24, {'quantified': 23, 'floor': 1, 'count_only': 0, 'not_detected': 0})

```

The runs table gives each run its replicate, then the trait, then any covariates, each column
named by its header:

```tsv
run	replicate	dose	injection
R0_0	R0	0.25	1
R0_1	R0	0.25	0
R1_0	R1	0.25	1
```

and so on for twelve runs. Then fit:

```pycon
>>> fit = pymzlib.ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv")
>>> fit.trait_name, fit.covariate_names
('dose', ['injection'])
>>> for r in fit.records:
...     print(f'{r["protein_accession"]:<14}{r["status"]:<8}effect {r["effect"]:.2f} '
...           f'(SE {r["standard_error"]:.2f}, df {r["df"]})  q {r["q"]:.4f}  '
...           f'{r["runs"]} runs, {r["replicates"]} replicates')
CENSORED_SITE fitted  effect 0.54 (SE 0.07, df 3)  q 0.0036  10 runs, 5 replicates
TRAIT_SITE    fitted  effect 0.50 (SE 0.06, df 4)  q 0.0016  12 runs, 6 replicates

```

`TRAIT_SITE` recovers the 0.5 it was built with. **`effect` is in logit (natural-log odds) units
per unit of the trait**, so `exp(0.5)` = 1.65: each dose unit multiplies the odds of the modified
form by about 1.65. The degrees of freedom are between replicates (six replicates, less the intercept and the
dose), because the dose never varies within a replicate.

`CENSORED_SITE` is the same site with replicate R5's two runs replaced by a floor and a ceiling.
Both leave the fit - only quantified, non-ceiling cells are responses - so it rests on five
replicates and three degrees of freedom.

**Support.** A site needs `min_runs` runs and `min_replicates` replicates with a response, mzLib's
defaults 6 and 4. Ask for six replicates, and `CENSORED_SITE` is below support:

```pycon
>>> strict = pymzlib.ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv",
...                                       min_replicates=6)
>>> for r in strict.records:
...     print(r["protein_accession"], r["status"], r["effect"] and round(r["effect"], 2))
CENSORED_SITE below_support None
TRAIT_SITE fitted 0.5

```

**In log2 odds.** Natural-log odds are mzLib's default. `scale="log2_odds"` reports the same fit in
log2 odds, where an effect of +1 doubles the odds of the modified form per trait unit:

```pycon
>>> log2 = pymzlib.ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv", scale="log2_odds")
>>> fit.scale, log2.scale
('natural_log_odds', 'log2_odds')
>>> for natural, base2 in zip(fit.records, log2.records):
...     print(f'{natural["protein_accession"]:<14}{natural["effect"]:.3f} -> {base2["effect"]:.3f}  '
...           f'p {natural["p_value"]:.5f} -> {base2["p_value"]:.5f}')
CENSORED_SITE 0.543 -> 0.783  p 0.00361 -> 0.00361
TRAIT_SITE    0.501 -> 0.723  p 0.00082 -> 0.00082

```

`TRAIT_SITE`'s 0.50 natural-log units are 0.72 log2 units: 2^0.72 = 1.65, the same odds ratio per
dose unit. mzLib does not refit. It divides the effects and standard errors by ln 2 and the two
variances by (ln 2)², so `df`, `p_value` and `q` are unchanged. Pick the scale your other results
use; the science does not change.

**What it does not tell you.**

- **Only the modified share.** Floors, count-only and undetected cells are left out, never entered
  as 0. Whether the modification is *detected* more often with the trait is a separate question
  this does not answer.
- **`q` covers one call.** It is Benjamini-Hochberg over the sites fitted in this call, not across
  datasets.
- **A trait that never varies, or a covariate that repeats it, cannot be fitted.** mzLib refuses
  the design, and pyMzLib raises `UsageError` with mzLib's reason.

## What these functions will not do

- **Read MetaMorpheus or FlashLFQ output directly.** Build the table as in
  [section 1](#1-the-input-one-row-per-peptidoform-per-run-per-protein); mzLib has no adapter
  yet.
- **Combine datasets.** Pairs and trait effects are per dataset. mzLib can pool pairs across
  datasets (`GlobalPairEngine`), which pyMzLib does not call yet.
- **Place a C-terminal modification**, since the table carries no protein length.
- **Do arithmetic.** Every number on this page came from mzLib's `Quantification.PtmQtl`;
  pyMzLib parses the bridge's answer into the classes in the
  [API reference](../reference.md#pymzlibptmqtl).

## See also

- [ptmqtl wire verbs](../reference/ptmqtl.md): every parameter, field, unit, error and caveat, as
  all three bindings publish them.
- [Readers](readers.md): the PSM and peptide tables the input is built from, and MetaMorpheus's
  own occupancy file.
- [Differential abundance](stats.md): Spearman correlation and Benjamini-Hochberg on their own.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
methods it uses:

--8<-- "docs/reference/_generated/cite.ptmqtl.md"

Further sources, in mzLib. They cite no paper for these rules, so the code is the reference:

- Bounds, the run covariate and the per-cell peptidoform counts: mzLib [#1468][1468]
  (`PtmPairEngine.CoVarying`, `Quantification/PtmQtl/PtmPairs.cs:238`; `PairTestStatus`, `:19`;
  `MinUntiedShare` and the simulation behind it, `:128`; `SiteRunOccupancy.CoveringPeptidoforms`,
  `Quantification/PtmQtl/SiteOccupancy.cs:99`).
- The log2-odds scale: mzLib [#1464][1464] (`OccupancyScale`,
  `Quantification/PtmQtl/SiteTraitEffects.cs:34`; `SiteTraitOptions.Scale`, `:66`).

[1468]: https://github.com/smith-chem-wisc/mzLib/pull/1468
[1464]: https://github.com/smith-chem-wisc/mzLib/pull/1464
