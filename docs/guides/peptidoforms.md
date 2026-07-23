# Peptidoforms

A mass spectrometrist looking at a protein usually wants one thing: *what peptides come off it, and
what fragments would I see for each?* pyMzLib answers that in a single call. Give it a UniProt
accession and it fetches the annotated entry, applies the modifications UniProt records, digests
the protein, and fragments every peptide.

```python
import pymzlib

digest = pymzlib.peptidoform.fragments("P02768")   # human serum albumin
print(len(digest.peptides))                         # 303
```

Everything below the surface is the same machinery MetaMorpheus uses in C#: the same digestion, the
same modification handling, the same fragment masses.

## The defaults are opinions, not placeholders

The one-argument call above already made several decisions. They are the choices this lab makes when
it has no reason to choose otherwise, and every one of them is a parameter you can change:

| Default | Meaning |
|---|---|
| `protease="trypsin|P"` | Tryptic with the Keil rule — cleave after K/R **except** before proline. |
| `dissociation="ETD"` | c and z• fragment ions. |
| `modifications=True` | Apply UniProt's annotated modifications. |
| `missed_cleavages=2` | Up to two missed cleavage sites per peptide. |
| `min_length=7` | **Discards peptides shorter than 7 residues** — see the trap below. |
| `max_modifications=2` | Consider up to two modifications per peptide. |
| `max_isoforms=1024` | Cap on modification isoforms per peptide — see the trap below. |
| `terminus="Both"` | Fragment from both ends. |

## Reading the result

`fragments()` returns a [`Digest`](../reference.md). Its `peptides` are
[`Peptide`](../reference.md) objects, each carrying its own [`Fragment`](../reference.md) ions:

```python
p = digest.peptides[0]

p.base_sequence      # 'DAHKSEVAHR'   — the bare sequence
p.full_sequence      # the sequence with modifications written inline, as mzLib renders them
p.monoisotopic_mass  # neutral monoisotopic mass, modifications included
p.missed_cleavages   # how many cleavage sites it spans
p.is_modified        # True if it carries at least one modification
p.modifications      # each applied modification, with its residue position

for f in p.fragments:
    f.product_type   # 'c' or 'zDot' for ETD
    f.fragment_number  # c3 is the third residue from the N-terminus
    f.neutral_mass   # monoisotopic NEUTRAL mass — not an m/z
```

Some convenience views:

```python
digest.modified_peptides   # only the peptides carrying a modification
digest.fragment_count      # total fragment ions across every peptide
```

## Modifications: what was annotated, and what could actually be used

This is the part a hand workflow gets wrong most often, so it is the part pyMzLib refuses to hide.

A UniProt annotation is only usable if it has both a **defined target** and a **defined mass**. A
glycosylation site annotated as *"N-linked (GlcNAc…) asparagine"* has neither a fixed composition
nor a fixed mass, so there is nothing to add to a peptide — and a peptide carrying an ambiguous-mass
PTM is not identifiable by mass anyway. Those annotations are correctly excluded. The problem is
that without being told, you cannot tell an exhaustive answer from a filtered one.

Every `Digest` carries a [`ModificationCensus`](../reference.md) that says exactly what happened:

```python
c = digest.modification_census
print(c.explain())
# 14 of 38 annotated modifications were applied, across 14 residue positions.
# Excluded by type: 24 × glycosylation site — these have no defined chemical composition, so no
# mass can be assigned …

c.annotated    # 38 — features UniProt lists
c.applied      # 14 — modifications actually placed
c.excluded     # 24 — annotated but unusable
c.sites        # distinct residue positions carrying a modification
c.unresolved   # names UniProt annotated but that are absent from its own ptmlist, so dropped
```

!!! warning "`sites` is not a modification count"
    A histone lists several alternatives at one residue — K9me1, K9me2, K9me3, K9ac are four
    modifications at **one** site. `sites` is always the smaller number.

`unresolved` catches the quietest failure of all. On histone H3.1, seven N6-lactoyllysine sites are
annotated by UniProt but absent from UniProt's *own* modification list, so no mass can be assigned
and they are dropped — while the per-type summary still says "modified residue … loaded". They land
in `unresolved` rather than vanishing.

## m/z, and why trimethyllysine needs care

`Peptide.mz(charge)` converts a peptide's mass to an m/z, and it handles two conventions that are
invisible in the answer if you get them wrong:

```python
p.mz(2)   # m/z at charge 2+
```

**It uses the proton mass (1.007276), not the hydrogen atom (1.007825).** The difference is
0.55 mDa — about 1.1 ppm at m/z 500, which on an Orbitrap is a match versus a miss. Libraries differ
on this and rarely say which they used.

**It does not double-count fixed charges.** Some modifications leave a residue permanently charged:
trimethylation of a lysine ε-amine gives a quaternary ammonium, and UniProt records the delta as
43.054227 (C₃H₇ *minus an electron*) rather than the neutral 43.054775. So the peptide's mass already
carries that charge, and `mz()` adds only `charge − fixed_charges` protons:

```python
p.fixed_charges   # 1 for a singly-trimethylated peptide
p.mz(2)           # correct — adds one proton, not two
```

Adding a full complement of protons would put a 2+ trimethylated peptide half a Thomson too high —
on the single most important histone modification there is. A peptide with a fixed charge is
observable at that charge with no protonation at all, so `charge` may not be **below**
`fixed_charges`.

!!! note "Fragments deliberately have no `mz()`"
    Converting a fragment correctly needs the fixed charge *within that fragment's span* — a c or z
    ion carries only the charged modifications on the residues it actually contains, not the whole
    peptide's `fixed_charges`. That per-fragment accounting isn't provided yet, so `Fragment` exposes
    `neutral_mass` only. For an unmodified or neutrally-modified peptide,
    `(neutral_mass + z * 1.007276) / z` is correct.

## Two defaults that bite silently

Both of these produce a *short* answer that looks exactly like a *complete* one — the failure mode
this whole feature exists to prevent, so they are worth knowing before you trust a count.

### `min_length=7` drops short peptides

The default digest discards every peptide shorter than seven residues. On a histone that is roughly a
third of the tryptic peptides — including short, heavily modified ones you may specifically care
about. If you asked for "every peptide", you did not get every peptide:

```python
# Keep everything, including single-residue peptides
digest = pymzlib.peptidoform.fragments("P68431", min_length=1)
```

### The isoform cap truncates silently

Modification isoforms are enumerated combinatorially, and it grows fast: histone H3.1 yields 49 bare
tryptic peptides, 2,563 at two modifications, and 7,040 at three. `max_isoforms=1024` caps the
isoforms considered *per peptide position*, and when it binds it discards the excess **without
error** — on H3.1 at four modifications, about 30% (13,700 down to 9,536).

`Digest` tells you when this happened, so a truncated answer is visible rather than merely short:

```python
if digest.truncated:
    print(f"{digest.peptides_at_cap} peptides hit the isoform cap — the list is incomplete")
    digest = pymzlib.peptidoform.fragments("P68431", max_isoforms=20000)
```

Always check `digest.truncated` before treating a peptidoform list as exhaustive.

## Protease naming: read this if you come from MaxQuant or Mascot

The default `protease="trypsin|P"` is the **reverse** of what the same-looking name means elsewhere,
and the difference is real — on serum albumin the two options differ by 37 peptides out of about 200.

| pyMzLib / mzLib | Behaviour | The MaxQuant / Mascot name for it |
|---|---|---|
| `"trypsin|P"` | Cleave after K/R **except** before proline (the classic Keil rule) | plain `Trypsin` |
| `"trypsin"` | Cleave after K/R, **including** before proline | `Trypsin/P` |

`"trypsin|P"` is the default here because it is what a mass spectrometrist usually means by "trypsin".
If you were expecting the MaxQuant `Trypsin/P` semantics (ignore the proline rule), you want mzLib's
plain `"trypsin"`.

## Tuning the digest

Every default is reachable. A few common adjustments:

```python
# The bare sequence, no modifications — useful as a control
pymzlib.peptidoform.fragments("P02768", modifications=False)

# HCD/CID b and y ions instead of ETD c/z•
pymzlib.peptidoform.fragments("P02768", dissociation="HCD")

# A tighter length window and no missed cleavages
pymzlib.peptidoform.fragments("P02768", missed_cleavages=0, min_length=6, max_length=40)

# Only N-terminal fragments
pymzlib.peptidoform.fragments("P02768", terminus="N")
```

`max_modifications` is the knob that most affects run time, because isoforms are combinatorial: three
modifications per peptide is dramatically more work than two.

## Errors you might hit

| Exception | Means |
|---|---|
| `UsageError` | The accession, protease, dissociation type, or terminus isn't recognised, or a numeric argument is negative. Raised before any network call. |
| `ServiceUnavailableError` | UniProt was down, rate-limiting, or timed out. **Not your bug** — retry later. |
| `BridgeError` | UniProt answered but something about the request was rejected. |

`ServiceUnavailableError` is a subclass of `BridgeError`, so catching `BridgeError` catches both;
separate them when you want to retry an outage but report the rest:

```python
try:
    digest = pymzlib.peptidoform.fragments("P02768")
except pymzlib.ServiceUnavailableError:
    ...   # UniProt's problem; worth a retry
except pymzlib.BridgeError as e:
    print(f"{e.error_type}: {e}")
```
