# Prediction

mzLib ships clients for **37 published models** on the [Koina](https://koina.wilhelmlab.org/)
inference server, across five families. pyMzLib calls them.

```python
import pymzlib

r = pymzlib.prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK", "ELVISLIVESK"])
r.columns["retention_time"]      # [5.5165, 129.7723]
r.retention_time_unit            # 'indexed_retention_time'  ← not minutes
```

| family | function | predicts |
|---|---|---|
| `retention_time` | [`retention_time()`](#retention-time) | when a peptide elutes — 10 models |
| `fragment_intensity` | [`fragments()`](#fragments) | MS2 m/z and relative intensity — 21 models |
| `collisional_cross_section` | [`ccs()`](#ccs-and-detectability) | ion-mobility CCS in Å² — 2 models |
| `detectability` | [`detectability()`](#ccs-and-detectability) | flyability, as 4 class probabilities — 1 model |
| `crosslink_intensity` | [`crosslink_fragments()`](#crosslinks) | MS2 for a crosslinked pair — 3 models |

## Start with `models()`

It is enumerated from mzLib rather than transcribed, so it reflects your installed version. Each
entry carries the constraints that decide whether a peptide can be sent **at all** — and those
constraints are more restrictive than most people expect.

```python
for m in pymzlib.prediction.models("retention_time"):
    print(m.model, m.max_peptide_length, m.allowed_unimod_ids)
```

```
AlphaPeptDeep_rt_generic 500 [...]
Chronologer_RT           30  [35, 4]
Deeplc_hela_hf           60  [...]
Prosit_2019_irt          30  [35, 4]
...
```

!!! warning "`max_peptide_length` is 30 for most models"
    That excludes a substantial fraction of real tryptic peptides, and every one of them comes back
    as a `None` with a warning rather than an error. Check `failed_row_count` before you average
    anything.

`allowed_unimod_ids` is the other quiet filter: **an empty list means the model accepts no
modifications at all**, and a peptide carrying one it was not trained on is refused rather than
silently stripped.

### Constraints are a tri-state, not a list

```python
hcd = [m for m in pymzlib.prediction.models() if m.model == "Prosit_2020_intensity_HCD"][0]
hcd.collision_energy.requirement      # 'any_value_required'  — you MUST send one
cid = [m for m in pymzlib.prediction.models() if m.model == "Prosit_2020_intensity_CID"][0]
cid.collision_energy.requirement      # 'not_applicable'      — this model is fixed at NCE 35
```

| `requirement` | meaning |
|---|---|
| `"not_applicable"` | the model has no such input — do not send it |
| `"any_value_required"` | you must send one; any value is accepted |
| `"one_of"` | you must send one of `constraint.values` |

!!! danger "This is where the raw mzLib field reads backwards"
    mzLib encodes all three states in one nullable set: `null` means *not applicable* and an
    **empty** set means *required, any value*. Read straight, that makes `Prosit_2020_intensity_CID`
    look like it accepts any collision energy — it accepts none — and `Prosit_2020_intensity_HCD`
    look like it accepts none, when in fact it requires one. pyMzLib translates the state into a
    word, which is the whole reason the wire contract carries availability rather than raw fields.

## Retention time

```python
r = pymzlib.prediction.retention_time("Prosit_2019_irt", peptides)
```

!!! danger "Most of these models do not return minutes"
    They return **iRT** — indexed retention time, a dimensionless scale anchored to standard
    peptides. Only `Chronologer_RT` returns absolute minutes.
    `Predictions.retention_time_unit` tells you which, as a value:

    ```python
    r.retention_time_unit    # 'indexed_retention_time' or 'minutes'
    ```

    Plotting an iRT against a gradient without first fitting the iRT-to-minutes line on shared
    peptides is the commonest way to misread these numbers, and nothing about the values themselves
    reveals the mistake — an iRT of 130 looks like a perfectly plausible 130-minute gradient.

## Fragments

```python
f = pymzlib.prediction.fragments(
    "Prosit_2020_intensity_HCD", ["PEPTIDEK", "ELVISLIVESK"],
    precursor_charge=2, collision_energy=28,
)
[len(a) for a in f.columns["fragment_mz"]]     # [28, 40]
```

!!! warning "The arrays are ragged, and that is the correct answer"
    Koina returns a fixed-width grid with `-1` marking ions that cannot exist for a given peptide,
    and mzLib drops those. So each row's three fragment arrays are as long as **that peptide's**
    possible ions — 28 for `PEPTIDEK`, 40 for `ELVISLIVESK`, from a model whose nominal count is
    174. Index them per row; treating them as a rectangle gives a wrong answer rather than an
    exception.

`fragment_annotations`, `fragment_mz` and `fragment_intensity` are index-aligned within a row.
Intensities are **relative**, on Koina's own 0–1 scale, and are not comparable with a measured
intensity or between models.

Per-peptide parameters override the shared default:

```python
pymzlib.prediction.fragments(
    "Prosit_2020_intensity_HCD",
    ["PEPTIDEK", {"sequence": "ELVISLIVESK", "precursor_charge": 3}],
    precursor_charge=2, collision_energy=28,
)
```

## CCS and detectability

```python
c = pymzlib.prediction.ccs("IM2Deep", ["PEPTIDEK"], precursor_charge=2)
c.columns["collisional_cross_section"]     # [327.53]
c.collisional_cross_section_unit           # 'square_angstroms'
```

!!! info "CCS is in Å², never 1/K0"
    Converting to the reduced mobility a timsTOF actually reports needs drift-gas temperature and
    pressure, which mzLib does not carry — so no conversion is offered rather than a guessed one.

```python
d = pymzlib.prediction.detectability("pfly_2024_fine_tuned", peptides)
d.columns["high_detectability"]        # a probability, per peptide
```

The four detectability columns are a **distribution over classes** and sum to 1. They are not an
expected intensity and not a probability of detection.

## Crosslinks

```python
x = pymzlib.prediction.crosslink_fragments(
    "Prosit_2023_intensity_XL_CMS2",
    [{"alpha_sequence": "PEPK[UNIMOD:1896]TIDEK", "beta_sequence": "ELVISK"}],
    precursor_charge=3, collision_energy=30,
)
```

!!! danger "This family speaks a different sequence language"
    The other four functions accept mzLib's `FullSequence` notation — `PEPTIDEK[Common
    Variable:Oxidation on M]` — and convert it. The crosslink models **reject it** and require raw
    UNIMOD brackets: `K[UNIMOD:1896]`. That is mzLib's constraint, not a choice made here, and it
    is repeated in `Predictions.caveats`. Same column name, two input languages.

## A failed row is still a row

Too long, an unsupported modification, a missing required collision energy — the value comes back
`None` and `warning` says why, rather than the row vanishing:

```python
f = pymzlib.prediction.fragments("Prosit_2020_intensity_HCD", ["PEPTIDEK"], precursor_charge=2)
f.failed_row_count        # 1
f.warnings                # [(0, 'Input is missing required parameter CollisionEnergy for this model.')]
```

That matters because predictions must line up with the peptides you sent. A library that dropped
unpredictable rows would silently shift every subsequent index.

## Koina is someone else's GPU

The service is **public, shared, community-run, free and unauthenticated**. pyMzLib keeps mzLib's
throttling defaults exactly as they are:

```python
pymzlib.prediction.retention_time(model, peptides, max_batches=200, throttle_ms=250)
```

Both knobs exist for a genuinely large job. They are not raised by default, because a binding that
maximised throughput out of the box would be spending capacity nobody here pays for.

There is also **no cancellation**. mzLib's public `Predict` is synchronous and threads no
cancellation token, so a large batch runs to completion or times out; `timeout=` bounds the wait
but does not stop the work already in flight.

## Predictions are opinions

Nothing here has been matched against a spectrum. No output is FDR-anything. That is obvious
stated plainly and remarkably easy to forget once the numbers are sitting in a DataFrame beside
measured ones — which is why every response carries `caveats` and why this page ends on it.

## What is not covered

- **The local Chronologer predictor.** mzLib also ships Chronologer as a local TorchSharp network:
  x64-only, extracting hundreds of megabytes of weights to a shared temp path, and racing any
  concurrent process doing the same. pyMzLib publishes an arm64 macOS wheel, so exposing it would
  either break that wheel or ship a function that fails on it. The same model is reachable as
  `Chronologer_RT` over Koina — but note the two report **different units**: absolute retention
  time over the network, % acetonitrile locally.
- **The local analytic predictors** (SSRCalc3, CZE). Pure arithmetic, no network, and a reasonable
  next tranche; they return a hydrophobicity index and a migration time respectively, neither of
  which is a retention time, so they need their own units on the wire.
- **Spectral-library output.** mzLib can turn predicted intensities into a `.msl` library, but its
  API reads state left behind by the previous `Predict` call, which a stateless wire contract
  cannot express without fusing the two into one verb. Worth doing; not done here.
