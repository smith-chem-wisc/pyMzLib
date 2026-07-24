# FlashLFQ

[FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) is the Smith lab's fast label-free
quantification engine — the same one MetaMorpheus uses. Give it a search result and the mzML runs
it came from, and it measures how much of each peptide and protein is in each run. pyMzLib exposes
it in one call.

```python
import pymzlib

result = pymzlib.flashlfq.quantify(
    psms="AllPSMs.psmtsv",                       # a MetaMorpheus search result
    spectra=["run_3.mzML", "run_4.mzML"],
    match_between_runs=True,
)

print(result.peptide_count, result.protein_count)   # 354 943
```

The whole pipeline is mzLib's own: the result file is read by mzLib's `Readers`, turned into
FlashLFQ identifications by mzLib's converter, and quantified by `FlashLfqEngine`. **MetaMorpheus is
not involved** — nothing beyond pyMzLib is installed or run.

Every name here is FlashLFQ's or mzLib's, unchanged: `match_between_runs`, `ppm_tolerance`,
`protein_groups`, `detection_types`, `FlashLfqResults`, `ProteinGroup`. What you read in the
[FlashLFQ paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00608) or the MetaMorpheus output
columns means the same thing here.

!!! danger "Quantify MetaMorpheus output only — an MSFragger `psm.tsv` gives wrong numbers"
    mzLib will accept an MSFragger `psm.tsv` as `psms`, and `pymzlib.readers.identify()` will call
    it `quantifiable`, because it genuinely implements the interface FlashLFQ consumes. The
    **numbers that come back are wrong, and wrong silently.**

    MSFragger records retention time in **seconds**; MetaMorpheus records **minutes**. mzLib's
    result-file readers pass the column straight through without converting it, and FlashLFQ then
    treats the value as minutes and searches a ±2-minute window around it. A peptide eluting at 60
    minutes is looked for at 1 minute and is simply not found — so quantification collapses toward
    zero instead of failing, which is the worst way for a number to be wrong.

    This is a defect in mzLib affecting every caller, not something pyMzLib introduces, and it is
    reported upstream. Until it is fixed, use `.psmtsv` / `.osmtsv`.

## What you get back

`quantify()` returns a [`FlashLfqResults`](../reference.md):

```python
result.identification_count   # 594 — PSMs read from the result file
result.peptide_count          # 354 — quantified peptides
result.protein_count          # 943 — quantified protein groups
result.spectra_files          # one SpectraFileInfo per run
result.peptides               # one Peptide per modified sequence
result.proteins               # one ProteinGroup per protein group
```

A peptide carries its intensity in each run, keyed by the run's base file name:

```python
p = result.peptides[0]
p.sequence            # 'AHQLVMEGYNWC[Common Fixed:Carbamidomethyl on C]HDR' — the modified sequence
p.base_sequence       # 'AHQLVMEGYNWCHDR'
p.protein_groups      # 'H0YC23;P62714;P67775'
p.intensity("run_3")  # 1930193.13  — this peptide's intensity in run_3
p.detection_type("run_3")   # 'MSMS'  — identified and quantified here
```

Proteins are the same shape:

```python
g = result.proteins[0]
g.protein_group       # 'H0YC23'
g.gene_name           # 'primary:PPP2CB'
g.intensity("run_3")  # protein-level intensity in run_3
```

## Match between runs

This is the reason to reach for FlashLFQ, and it is off by default because it makes an inference you
should opt into. With `match_between_runs=True`, a peptide identified in one run but **missing** from
another is still quantified in the second, by transferring the identification across the aligned
retention-time axis. Those transferred peaks are counted per run:

```python
for f in result.spectra_files:
    print(f.file_name, f.peak_count, f.mbr_peak_count)
# run_3   340   62
# run_4   307   78

result.mbr_peak_count   # 140 — total transferred peaks; the number to report
```

So on this pair, 140 peaks were transferred into a run where the peptide was never identified — the
values a run-by-run analysis would have left missing.

!!! warning "For MBR, read `result.peaks` — not the peptide table"
    The per-peptide roll-up (`result.peptides`) mirrors FlashLFQ's `QuantifiedPeptides.tsv`, and it
    **does not fully carry MBR**: of the 140 transfers here, only 52 appear as `"MBR"` at the peptide
    level, and *none* of run_3's 62 transfers do — they read as `"NotDetected"` with intensity `0.0`.
    So if you build a peptide × run matrix from `p.intensity(run)`, you silently drop most of the
    values MBR just filled in. Use the **peaks**, which represent every transfer:

    ```python
    for peak in result.mbr_peaks:                 # exactly the transferred peaks
        print(peak.file_name, peak.sequence, peak.intensity)

    # Build an MBR-inclusive peptide × run matrix from the peaks, not the peptides:
    import collections
    matrix = collections.defaultdict(dict)
    for peak in result.peaks:
        matrix[peak.sequence][peak.file_name] = peak.intensity
    ```

    A transferred peak reports `peak.detection_type == "MBR"`; a directly identified one, `"MSMS"`.

MBR has its own tolerances, defaulted the way FlashLFQ defaults them:

```python
pymzlib.flashlfq.quantify(
    psms="AllPSMs.psmtsv", spectra=[...],
    match_between_runs=True,
    mbr_ppm_tolerance=10.0,          # mass window for a transfer
    mbr_q_value_threshold=0.05,      # the FDR control on transfers — keep it
)
```

`mbr_q_value_threshold` is not a knob to loosen casually: it is the **false-discovery control** that
makes a transfer trustworthy. A transfer is just "some peak at about the right mass and retention
time," and at a wide-open threshold most such matches are coeluting noise — the reason match-between-
runs is hard, and the reason a purpose-built engine is worth using rather than matching peaks by
hand. FlashLFQ holds transfers to this q-value using a decoy model; leaving it at the default `0.05`
is what separates a rescued peptide from a fabricated one.

## Experimental design: conditions and replicates

The simplest `spectra` is a list of paths — each run becomes its own biological replicate with no
condition, exactly as MetaMorpheus defaults when you give it no design file. To group runs into
conditions and replicates, pass a mapping per run instead, using FlashLFQ's `SpectraFileInfo` field
names:

```python
result = pymzlib.flashlfq.quantify(
    psms="AllPSMs.psmtsv",
    spectra=[
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
`technical_replicate`, `fraction`; anything omitted takes FlashLFQ's default.

!!! warning "Match-between-runs needs a *complete* design"
    MBR transfers across the experimental design, and it assumes the design is **complete and
    balanced**: every condition and biological replicate must carry the **same set of fractions** —
    no missing replicates, no missing fractions. A gap (one condition missing a fraction, one
    replicate absent) breaks the complementarity MBR relies on and makes the transfers unreliable.
    If you fractionated, make sure every sample has every fraction before turning MBR on.

!!! warning "Runs are matched to identifications by base file name"
    FlashLFQ links each identification to its run by **base file name** — the "File Name" column in
    the PSM file must match an mzML you pass. If the PSM file names a run you did not provide, the
    call fails up front and tells you which. Base file names must also be unique across `spectra`.

## Writing the FlashLFQ TSVs

Point `output_directory` at a folder and FlashLFQ writes its standard tables there, in addition to
returning the result:

```python
result = pymzlib.flashlfq.quantify(
    psms="AllPSMs.psmtsv", spectra=[...], match_between_runs=True,
    output_directory="flashlfq_out",
)
# flashlfq_out/QuantifiedPeaks.tsv, QuantifiedPeptides.tsv, QuantifiedProteins.tsv
```

## Two limits worth knowing

Both are surfaced rather than hidden — a wrong number that looks right is worse than an error.

### mzML only, for now

FlashLFQ can read Thermo `.raw` and Bruker data, but pyMzLib's quant surface accepts **mzML only**
at present. A non-mzML path is rejected before any work starts. Convert `.raw`/`.d` to mzML first
(e.g. with ThermoRawFileParser or `msconvert`).

### A protein intensity can be `None` — but a peptide intensity never is

FlashLFQ's median-polish protein quant marks a protein `NaN` when its peptide matrix is **degenerate**
— too few peptides per run to resolve (e.g. one peptide in each of two runs, on an anti-diagonal), or
several runs reporting the same intensity. It is a real artifact of the algorithm, documented in
mzLib's own tests, that protects you from a fabricated-looking number. `NaN` is not valid JSON, so
pyMzLib returns it as `None`: *"could not be quantified"*, not a silently wrong value.

```python
unquantified = [g.protein_group for g in result.proteins if g.intensity("run_3") is None]
```

!!! danger "`None` is proteins-only. A peptide that's missing reads `0.0`, not `None`."
    `Peptide.intensity(run)` returns **`0.0`** where the peptide was not quantified — never `None`.
    So the `is None` filter above works for proteins but finds nothing for peptides. Treat a peptide
    `0.0` as *missing*, not as a measured zero: don't log-transform it, and remember (per the MBR
    warning above) that an MBR-transferred peptide can read `0.0` here even though it *was*
    quantified — its value is in `result.peaks`. The full `detection_type` vocabulary tells you
    which case you're in: `"MSMS"`, `"MBR"`, `"MSMSIdentifiedButNotQuantified"` (identified here but
    no usable peak), `"MSMSAmbiguousPeakfinding"` (more than one peptide fits the peak), or
    `"NotDetected"`.

## Parameters

All defaulted as FlashLFQ defaults them; all named as FlashLFQ names them.

| Parameter | Meaning |
|---|---|
| `normalize` | Normalize intensities across runs. |
| `ppm_tolerance` | Mass tolerance for peak-finding (ppm). |
| `isotope_ppm_tolerance` | Mass tolerance for isotope-envelope matching (ppm). |
| `integrate` | Integrate peaks rather than take the apex — FlashLFQ recommends leaving off. |
| `match_between_runs` | Transfer identifications to quantify peptides missing from a run. |
| `mbr_ppm_tolerance` / `mbr_q_value_threshold` | Tolerance and confidence cutoff for MBR transfers. |
| `use_shared_peptides_for_protein_quant` | Let peptides shared between groups contribute to protein quant. |
| `bayesian_protein_quant` | Run FlashLFQ's Bayesian protein fold-change engine. |
| `use_pep_q_value` | Filter identifications on PEP q-value rather than q-value. |
| `max_threads` | Worker threads; `-1` lets FlashLFQ choose. |
| `output_directory` | Also write the FlashLFQ TSVs here. |
| `timeout` | Seconds to allow; `None` waits indefinitely. |

## Errors you might hit

| Exception | Means |
|---|---|
| `UsageError` | A run isn't mzML, an mzML is missing, the PSM file names a run you didn't provide, or an argument is malformed. Raised before any spectra are read. |
| `BridgeError` | FlashLFQ itself failed while quantifying. |

```python
try:
    result = pymzlib.flashlfq.quantify(psms="AllPSMs.psmtsv", spectra=[...], match_between_runs=True)
except pymzlib.UsageError as e:
    print(f"check the inputs: {e}")
except pymzlib.BridgeError as e:
    print(f"{e.error_type}: {e}")
```
