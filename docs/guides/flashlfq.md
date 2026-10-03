# FlashLFQ

[FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) is the Smith lab's label-free
quantification engine, the one MetaMorpheus uses. Give it a search result and the mzML runs it came
from, and it measures how much of each peptide and protein is in each run. The whole pipeline is
mzLib's: mzLib's readers read the result file, mzLib's converter makes FlashLFQ identifications, and
`FlashLfqEngine` quantifies them. **MetaMorpheus is not involved**, and nothing beyond pyMzLib is
installed or run.

| You want to | Call | mzLib does it with |
|---|---|---|
| Quantify peptides and proteins across runs | [`quantify()`](#what-you-get-back) | `FlashLfqEngine.Run` |
| Fill in peptides missing from a run | [`quantify(match_between_runs=True)`](#match-between-runs) | FlashLFQ's match-between-runs |
| Group runs into conditions and replicates | [`quantify(spectra=[{...}])`](#experimental-design-conditions-and-replicates) | `SpectraFileInfo` |
| Re-roll proteins from a peptide table, without the mzML | [`median_polish()`](#re-quantifying-proteins-from-a-peptide-table-median-polish) | `FlashLfqResults.CalculateProteinResultsMedianPolish` |

Every name is FlashLFQ's or mzLib's, unchanged: `match_between_runs`, `ppm_tolerance`,
`protein_groups`, `detection_types`, `FlashLfqResults`, `ProteinGroup`. What you read in the
FlashLFQ paper or in MetaMorpheus's output columns means the same thing here.

Every `>>>` example below runs in CI against a FlashLFQ run recorded through the real bridge, on
mzLib's own test data: two K562 runs and their MetaMorpheus search, with match-between-runs on.
Blocks titled **Not run** say why.

```pycon
>>> import pymzlib
>>> runs = ["20100614_Velos1_TaGe_SA_K562_3.mzML", "20100614_Velos1_TaGe_SA_K562_4.mzML"]
>>> result = pymzlib.flashlfq.quantify(
...     "AllPSMs.psmtsv",            # a MetaMorpheus search result
...     runs,
...     match_between_runs=True,
...     max_threads=1,
... )
>>> result.identification_count, result.peptide_count, result.protein_count
(594, 354, 943)

```

**Any quantifiable search result works**, not only MetaMorpheus's: an MSFragger `psm.tsv` or a
DIA-NN `report.tsv` too. Since mzLib#1116, mzLib converts MSFragger's retention times from seconds
to minutes as it reads them, so FlashLFQ searches the right window.

## What you get back

`quantify()` returns a [`FlashLfqResults`][pymzlib.flashlfq.FlashLfqResults]. Every identification
in the file was given to FlashLFQ, unfiltered; it is FlashLFQ that decides what it can quantify.

A peptide carries its intensity in each run, keyed by the run's base file name:

```pycon
>>> run3, run4 = (f.file_name for f in result.spectra_files)
>>> run3
'20100614_Velos1_TaGe_SA_K562_3'
>>> p = result.peptides[0]
>>> p.sequence                    # the modified sequence: the identity FlashLFQ quantifies
'AHQLVMEGYNWC[Common Fixed:Carbamidomethyl on C]HDR'
>>> p.base_sequence, p.protein_groups
('AHQLVMEGYNWCHDR', 'H0YC23;P62714;P67775')
>>> round(p.intensity(run3)), p.detection_type(run3)
(1930193, 'MSMS')

```

Proteins are the same shape:

```pycon
>>> g = next(g for g in result.proteins if g.intensity(run3))
>>> g.protein_group, g.gene_name
('P50990', 'primary:CCT8')
>>> round(g.intensity(run3)), round(g.intensity(run4))
(13493447, 3022424)

```

Most protein groups read 0 in both runs, and that is a setting, not a failure. With
`use_shared_peptides_for_protein_quant` off (the default), only a group's own peptides count, and
most groups here have none:

```pycon
>>> unique = {p.protein_groups for p in result.peptides if ";" not in p.protein_groups}
>>> zero = [g for g in result.proteins if g.intensity(run3) == 0 and g.intensity(run4) == 0]
>>> len(zero), sum(g.protein_group not in unique for g in zero)
(847, 842)

```

## Match between runs

This is the reason to reach for FlashLFQ, and it is off by default because it makes an inference you
should opt into. With `match_between_runs=True`, a peptide identified in one run but **missing**
from another is still quantified in the second, by transferring the identification across the
aligned retention-time axis. The transferred peaks are counted per run:

```pycon
>>> [(f.file_name[-1], f.peak_count, f.mbr_peak_count) for f in result.spectra_files]
[('3', 340, 62), ('4', 307, 78)]
>>> result.mbr_peak_count            # total transferred peaks: the number to report
140

```

So 140 peaks were quantified in a run where their peptide was never identified: values a run-by-run
analysis would have left missing.

!!! warning "For MBR, read `result.peaks`, not the peptide table"
    The per-peptide roll-up (`result.peptides`) mirrors FlashLFQ's `QuantifiedPeptides.tsv`, and it
    **does not carry most MBR transfers**. Of the 140 here, few appear as `"MBR"` at the peptide
    level, and none of run 3's do: they read `"NotDetected"` with intensity `0.0`.

    ```pycon
    >>> sum(p.detection_type(run3) == "MBR" for p in result.peptides)
    0
    >>> sum(p.detection_type(run4) == "MBR" for p in result.peptides)
    21

    ```

    So a peptide × run matrix built from `p.intensity(run)` silently drops most of the values MBR
    filled in. Build it from the **peaks**, which carry every transfer:

    ```pycon
    >>> peak = result.mbr_peaks[0]
    >>> peak.file_name == run3, peak.sequence, peak.detection_type
    (True, 'RVHVTQEDFEMAVAK', 'MBR')
    >>> import collections
    >>> matrix = collections.defaultdict(dict)
    >>> for peak in result.peaks:
    ...     matrix[peak.sequence][peak.file_name] = peak.intensity
    >>> round(matrix["RVHVTQEDFEMAVAK"][run3])
    502271

    ```

    A peptide with several peaks in one run keeps the last one in this loop; FlashLFQ's own
    roll-up keeps one of them too, not their sum.

MBR has its own tolerances, defaulted the way FlashLFQ defaults them:

```python title="Not run: the defaults, spelled out; the recorded run above used them"
pymzlib.flashlfq.quantify(
    "AllPSMs.psmtsv", runs,
    match_between_runs=True,
    mbr_ppm_tolerance=10.0,          # mass window for a transfer, in ppm
    mbr_q_value_threshold=0.05,      # the FDR control on transfers: keep it
)
```

`mbr_q_value_threshold` is not a setting to loosen casually. It is the **false-discovery control**
that makes a transfer trustworthy. A transfer is "some peak at about the right mass and retention
time", and at a wide-open threshold most such matches are co-eluting noise. FlashLFQ holds
transfers to this q-value with a decoy model, and leaving it at `0.05` is what separates a rescued
peptide from a fabricated one.

!!! note "Threads and reproducibility"
    Before mzLib#1155, match-between-runs with more than one thread could give different numbers
    for identical inputs, because FlashLFQ assembled its PEP training rows in thread-completion
    order. That is fixed, but pyMzLib has not re-measured its own reproduction since, so
    `max_threads=1`, as above, is the conservative choice when a result must reproduce exactly.

## Experimental design: conditions and replicates

The simplest `spectra` is a list of paths. Each run becomes its own biological replicate with no
condition, as MetaMorpheus does when you give it no design file:

```pycon
>>> [(f.condition, f.biological_replicate) for f in result.spectra_files]
[('', 0), ('', 1)]

```

To group runs into conditions and replicates, pass a mapping per run, using FlashLFQ's
`SpectraFileInfo` field names:

```python title="Not run: a two-condition design needs four runs; mzLib's test data has two"
result = pymzlib.flashlfq.quantify(
    "AllPSMs.psmtsv",
    [
        {"path": "control_1.mzML", "condition": "control", "biological_replicate": 1},
        {"path": "control_2.mzML", "condition": "control", "biological_replicate": 2},
        {"path": "treated_1.mzML", "condition": "treated", "biological_replicate": 1},
        {"path": "treated_2.mzML", "condition": "treated", "biological_replicate": 2},
    ],
    match_between_runs=True,
    normalize=True,          # normalize intensities across runs
)
```

You can mix paths and mappings, and supply any of `condition`, `biological_replicate`,
`technical_replicate` and `fraction`. Anything omitted takes FlashLFQ's default. If your design is
in an SDRF, [`pymzlib.sdrf.design()`](sdrf.md#turn-an-sdrf-into-a-quantification-design) turns it
into exactly these mappings.

!!! warning "Match-between-runs needs a *complete* design"
    MBR transfers across the experimental design, and it assumes the design is **complete and
    balanced**: every condition and biological replicate carries the **same set of fractions**, with
    no missing replicates and no missing fractions. A gap breaks the complementarity MBR relies on
    and makes the transfers unreliable. If you fractionated, make sure every sample has every
    fraction before turning MBR on.

!!! warning "Runs are matched to identifications by base file name"
    FlashLFQ links each identification to its run by **base file name**: the "File Name" column in
    the PSM file must match an mzML you pass. If the PSM file names a run you did not provide, the
    call fails before any spectra are read and says which. Base file names must also be unique
    across `spectra`.

## Writing the FlashLFQ TSVs

Point `output_directory` at a folder and FlashLFQ also writes its standard tables there:

```python title="Not run: writes FlashLFQ's tables to disk"
result = pymzlib.flashlfq.quantify(
    "AllPSMs.psmtsv", runs, match_between_runs=True,
    output_directory="flashlfq_out",
)
# flashlfq_out/QuantifiedPeaks.tsv, QuantifiedPeptides.tsv, QuantifiedProteins.tsv
```

With `bayesian_protein_quant=True` it also writes `BayesianProteinQuant.tsv`. That file is the only
place the Bayesian results appear: they are not on the returned object.

## Re-quantifying proteins from a peptide table: median polish

Protein quantification is the *second half* of `quantify()`: once peptides are measured, FlashLFQ
rolls them up to proteins with **median polish**. `median_polish()` runs that half on its own,
starting from a `QuantifiedPeptides.tsv` FlashLFQ already wrote, with no mzML and no peak-finding.
Here it reads the table the run above wrote:

```pycon
>>> proteins = pymzlib.flashlfq.median_polish("K562_QuantifiedPeptides.tsv")
>>> len(proteins)
943

```

The roll-up is mzLib's own `CalculateProteinResultsMedianPolish`, the method `quantify()` runs, so
the numbers match. Every protein group agrees with the full run:

```pycon
>>> full = {g.protein_group: g.intensities for g in result.proteins}
>>> all(g.intensities == full[g.protein_group] for g in proteins)
True

```

Reach for it to **re-roll proteins without re-quantifying peptides**: try another experimental
design, or toggle shared peptides, in seconds rather than re-reading every run.

### The design is how you group replicates

Median polish compares each peptide across the samples grouped by **condition and biological
replicate**, so the design tells it which columns are replicates of which sample. Pass one mapping
per run, keyed by the run's `Intensity_<name>` column:

```python title="Not run: a two-condition design needs four runs; mzLib's test data has two"
proteins = pymzlib.flashlfq.median_polish(
    "flashlfq_out/QuantifiedPeptides.tsv",
    design=[
        {"file_name": "control_1", "condition": "control", "biological_replicate": 0},
        {"file_name": "control_2", "condition": "control", "biological_replicate": 1},
        {"file_name": "treated_1", "condition": "treated", "biological_replicate": 0},
        {"file_name": "treated_2", "condition": "treated", "biological_replicate": 1},
    ],
    use_shared_peptides=True,             # let shared peptides contribute
)
proteins[0].intensity("control_1")       # keyed by "condition_biorep"
```

With **no** `design`, as above, each `Intensity_` column becomes its own biological replicate with
a blank condition, which is what FlashLFQ assumes when it writes the file with no design, and the
intensities are keyed by run base name. When a design *is* given it must name every run in the
table and only runs in the table: a name that matches no column, or a column with no design line,
is rejected rather than guessed.

The returned objects are ordinary `ProteinGroup`s, so an intensity is **`None`** where median polish
could not resolve a number (the degenerate-matrix case below) and **`0.0`** where the protein was not
measured in that sample. That includes a group whose only peptides are shared, when
`use_shared_peptides` is off. Pass `output_directory=...` to also write a `QuantifiedProteins.tsv`;
the list is the primary result and the file is a convenience.

!!! note "The written file's columns can come in a different order"
    The `QuantifiedProteins.tsv` sample labels match the keys of the returned objects. Before
    [mzLib#1129](https://github.com/smith-chem-wisc/mzLib/pull/1129), FlashLFQ inverted the labelling
    rule for unfractionated data
    ([mzLib#1128](https://github.com/smith-chem-wisc/mzLib/issues/1128)); the bridge includes the fix.
    The column **order** can still differ: the file lists conditions in the order they first appear,
    and the returned objects sort them. Look columns up by name, not by position.

## Two limits worth knowing

Both are reported rather than hidden: a wrong number that looks right is worse than an error.

### mzML only, for now

FlashLFQ can read Thermo `.raw` and Bruker data, but pyMzLib's quantification accepts **mzML only**
for now. A path that is not mzML is rejected before any work starts. Convert `.raw` or `.d` files
to mzML first, for example with ThermoRawFileParser or `msconvert`.

### A protein intensity can be `None`, but a peptide intensity never is

FlashLFQ's median polish marks a protein `NaN` when its peptide matrix is **degenerate**: too few
peptides per run to resolve, such as one peptide in each of two runs on an anti-diagonal, or several
runs reporting the same intensity. mzLib's own tests document it, and it protects you from a
fabricated-looking number. `NaN` is not valid JSON, so pyMzLib returns it as `None`, meaning *"could
not be quantified"*. A handful of groups in this run are like that:

```pycon
>>> sorted(g.protein_group for g in result.proteins if g.intensity(run3) is None)[:3]
['Q13033', 'Q15149', 'Q15393']

```

!!! danger "`None` is proteins only. A missing peptide reads `0.0`, not `None`."
    `Peptide.intensity(run)` returns **`0.0`** where the peptide was not quantified, never `None`, so
    the `is None` filter above finds nothing for peptides. Treat a peptide `0.0` as *missing*, not
    as a measured zero: do not log-transform it, and remember that an MBR-transferred peptide can
    read `0.0` here even though it *was* quantified, with its value in `result.peaks`. The
    `detection_type` tells you which case you are in: `"MSMS"`, `"MBR"`,
    `"MSMSIdentifiedButNotQuantified"` (identified here, but no usable peak),
    `"MSMSAmbiguousPeakfinding"` (more than one peptide fits the peak) or `"NotDetected"`.

    ```pycon
    >>> sorted({p.detection_type(run3) for p in result.peptides})
    ['MSMS', 'MSMSAmbiguousPeakfinding', 'MSMSIdentifiedButNotQuantified', 'NotDetected']

    ```

## Parameters

All defaulted as FlashLFQ defaults them, and named as FlashLFQ names them. The
[quant reference](../reference/quant.md) has each one's unit, range and wire spelling.

| Parameter | Meaning |
|---|---|
| `normalize` | Normalize intensities across runs. |
| `ppm_tolerance` | Mass tolerance for peak-finding, in ppm. |
| `isotope_ppm_tolerance` | Mass tolerance for isotope-envelope matching, in ppm. |
| `integrate` | Integrate peaks rather than take the apex. FlashLFQ recommends leaving it off. |
| `match_between_runs` | Transfer identifications to quantify peptides missing from a run. |
| `mbr_ppm_tolerance`, `mbr_q_value_threshold` | Tolerance and confidence cutoff for MBR transfers. |
| `use_shared_peptides_for_protein_quant` | Let peptides shared between groups contribute to protein quant. |
| `bayesian_protein_quant` | Also run FlashLFQ's Bayesian protein fold-change engine; results go to `output_directory` only. |
| `use_pep_q_value` | Carry each identification's PEP q-value as its q-value. It filters nothing. |
| `max_threads` | Worker threads; `-1` lets FlashLFQ choose. |
| `output_directory` | Also write the FlashLFQ TSVs here. |
| `timeout` | Seconds to allow; `None` waits indefinitely. |

## Errors you might hit

| Exception | Means |
|---|---|
| `UsageError` | A run is not mzML, an mzML is missing, the PSM file names a run you did not provide, or an argument is malformed. Raised before any spectra are read. |
| `BridgeError` | FlashLFQ itself failed while quantifying. |

A malformed argument is refused before the bridge starts:

```pycon
>>> pymzlib.flashlfq.quantify("AllPSMs.psmtsv", [])
Traceback (most recent call last):
    ...
pymzlib._bridge.UsageError: At least one spectra file is required.

```

[Every error each verb can raise](../errors.md) is listed on one page.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and FlashLFQ:

--8<-- "docs/reference/_generated/cite.flashlfq.md"
