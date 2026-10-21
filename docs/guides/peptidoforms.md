# Peptidoforms

A mass spectrometrist looking at a protein usually wants one thing: *which peptides come off it, and
which fragments would I see for each?* pyMzLib answers that in one call. Give it a UniProt
accession and it fetches the annotated entry, applies the modifications UniProt records, digests
the protein, and fragments every peptide. The digestion, the modification handling and the fragment
masses are the ones MetaMorpheus uses. pyMzLib can also rewrite a MetaMorpheus full sequence in
Unimod or ProForma notation, so it can be handed to another tool (see
[Converting full sequences](#converting-full-sequences-to-unimod-or-proforma)).

| You want to know | Call | mzLib does it with |
|---|---|---|
| Which peptides and fragments a protein gives | [`fragments()`](#reading-the-result) | `Protein.Digest`, `PeptideWithSetModifications.Fragment` |
| Which UniProt modifications were used, and which were not | [`Digest.modification_census`](#modifications-what-was-annotated-and-what-could-be-used) | the UniProt XML reader and its PTM list |
| A peptide's m/z at a charge | [`Peptide.mz()`](#mz-and-why-trimethyllysine-needs-care) | `ClassExtensions.ToMz`, with fixed charges accounted for |
| Whether the list is complete | [`Digest.truncated`](#the-isoform-cap-truncates-silently) | `DigestionParams.MaxModificationIsoforms` |
| A MetaMorpheus full sequence in Unimod accessions | [`convert()`](#converting-full-sequences-to-unimod-or-proforma) | `SequenceConversionService.Convert`, Unimod serializer |
| The same sequence in ProForma | [`convert(target="ProForma")`](#proforma-gives-the-same-accessions) | `SequenceConversionService.Convert`, ProForma serializer |
| Which sequences mzLib could not convert, and why | [`SequenceConversions.not_converted`](#when-mzlib-cannot-convert-a-sequence) | `ConversionWarnings` |
| Which notations mzLib can read and write | [`SequenceConversions.source_formats` and `.target_formats`](#converting-full-sequences-to-unimod-or-proforma) | `SequenceConversionService.AvailableSourceFormats` and `AvailableTargetFormats` |

Every `>>>` example on this page runs in CI against output recorded from the real bridge, so the
numbers you see are the numbers pyMzLib gives. Blocks titled **Not run** say why they are not.

```pycon
>>> import pymzlib
>>> digest = pymzlib.peptidoform.fragments("P02768")   # human serum albumin
>>> digest.full_name, digest.organism, digest.sequence_length
('Albumin', 'Homo sapiens', 609)
>>> len(digest.peptides), len(digest.modified_peptides)
(310, 108)

```

## The defaults are mzLib's

The one-argument call leaves the digestion to mzLib. Each option you do not pass (the protease,
missed cleavages, the length window, the modification and isoform caps, and the terminus) takes
the default of mzLib's `DigestionParams`, and pyMzLib sends nothing for it. The
[`fragments` reference](../reference/peptidoform.md#fragments) gives each value at the pinned
mzLib. Two choices are not digestion defaults:

- **`modifications=True`** applies UniProt's annotated modifications.
- **`dissociation`** decides the fragment series. mzLib names no default dissociation type, so HCD
  (b and y ions) is used until it does
  ([smith-chem-wisc/mzLib#1488](https://github.com/smith-chem-wisc/mzLib/issues/1488)).

Two of mzLib's defaults drop peptides without saying so ([below](#two-defaults-that-bite-silently)),
and its default protease may not be the trypsin you expect
([below](#protease-naming-read-this-if-you-come-from-maxquant-or-mascot)).

The result echoes the settings it ran with, so a saved `Digest` says how it was made:

```pycon
>>> digest.protease, digest.dissociation, digest.max_modifications, digest.max_isoforms
('trypsin', 'HCD', 2, 1024)

```

## Reading the result

`fragments()` returns a [`Digest`][pymzlib.peptidoform.Digest]. Its `peptides` are
[`Peptide`][pymzlib.peptidoform.Peptide] objects in digestion order, each carrying its own
[`Fragment`][pymzlib.peptidoform.Fragment] ions. Here is the first modified one:

```pycon
>>> p = digest.modified_peptides[0]
>>> p.base_sequence
'TCVADESAENCDK'
>>> p.full_sequence                     # modifications inline, as mzLib writes them
'TCVADES[UniProt:Phosphoserine on S]AENCDK'
>>> round(p.monoisotopic_mass, 4)       # neutral monoisotopic mass, modifications included
1463.4946
>>> p.missed_cleavages, p.is_modified
(0, True)

```

Each fragment carries a **neutral** mass, not an m/z:

```pycon
>>> f = p.fragments[0]
>>> f.product_type, f.fragment_number, round(f.neutral_mass, 4)
('b', 1, 101.0477)

```

`fragment_number` counts residues from the fragment's own terminus: b1 is the first residue from
the N-terminus, y1 the first from the C-terminus.

## Modifications: what was annotated, and what could be used

This is the part a hand workflow gets wrong most often, so pyMzLib refuses to hide it.

mzLib loads only UniProt's `modified residue` and `lipid moiety-binding region` annotations. It
drops every other feature type, such as glycosylation sites, **on the feature type alone**, before
any mass lookup. That is usually right: a glycation or glycosylation annotation describes a labile,
heterogeneous adduct, and giving it one exact mass would describe a species you cannot observe. The
problem is that without being told, you cannot tell an exhaustive answer from a filtered one
(smith-chem-wisc/mzLib#1112).

Every `Digest` carries a [`ModificationCensus`][pymzlib.peptidoform.ModificationCensus] that says
what happened:

```pycon
>>> c = digest.modification_census
>>> c.annotated, c.applied, c.excluded, c.sites
(38, 14, 24, 14)
>>> c.by_type
[{'type': 'glycosylation site', 'count': 24, 'loaded': False}, {'type': 'modified residue', 'count': 14, 'loaded': True}]
>>> print(c.explain())                                      # doctest: +ELLIPSIS
14 of 38 annotated modifications were applied, across 14 residue positions. Excluded by type: 24 × glycosylation site ...

```

!!! warning "`sites` is not a modification count"
    A histone lists several alternatives at one residue: K9me1, K9me2, K9me3 and K9ac are four
    modifications at **one** site. `sites` is the number of residue positions, so on a histone it
    is the smaller number.

`unresolved` catches the quietest failure. It lists modification names UniProt annotates that are
absent from UniProt's *own* modification list, so no mass can be assigned and mzLib drops them. On
histone H3.1, seven N6-lactoyllysine sites land here instead of vanishing. Albumin has none:

```pycon
>>> c.unresolved
[]

```

## m/z, and why trimethyllysine needs care

`Peptide.mz(charge)` converts a peptide's neutral mass to an m/z, and it handles two conventions
that are invisible in the answer if you get them wrong:

```pycon
>>> round(p.mz(2), 4)
732.7546

```

**It adds the proton mass (1.007276 Da), not the hydrogen atom (1.007825 Da).** The difference is
0.55 mDa, about 1.1 ppm at m/z 500. On an Orbitrap that is the difference between a match and a
miss. Libraries differ on this and rarely say which they used.

**It does not double-count fixed charges.** Some modifications leave a residue permanently charged.
Trimethylation of a lysine ε-amine gives a quaternary ammonium, and UniProt records its delta as
43.054227 Da (C₃H₇ *minus an electron*) rather than the neutral 43.054775 Da. The peptide's mass
already carries that charge, so `mz()` adds only `charge − fixed_charges` protons. No albumin
peptide carries a fixed charge:

```pycon
>>> {q.fixed_charges for q in digest.peptides}
{0}

```

On a trimethylated histone peptide, `fixed_charges` is 1, and `mz(2)` adds one proton, not two.
Adding a full complement would put a 2+ trimethylated peptide half a Thomson too high, on one of the
most studied histone modifications there is. A peptide with a fixed charge is observable at that
charge with no protonation at all, so `charge` may not be **below** `fixed_charges`.

!!! note "Fragments deliberately have no `mz()`"
    Converting a fragment correctly needs the fixed charge *within that fragment's span*: a b, y, c or
    z ion carries only the charged modifications on the residues it contains, not the whole peptide's
    `fixed_charges`. That per-fragment accounting is not provided yet, so `Fragment` has
    `neutral_mass` only. For an unmodified or neutrally modified peptide,
    `(neutral_mass + z × 1.007276) / z` is correct.

## Two defaults that bite silently

Both produce a *short* answer that looks exactly like a *complete* one, so they are worth knowing
before you trust a count.

### The minimum length drops short peptides

mzLib's default minimum length discards every shorter peptide. On histone H3.1 the default digest
keeps 49 of the 96 peptides that `min_length=1` gives, and the ones it drops include short, heavily
modified peptides you may care about most. If you asked for "every peptide", pass `min_length=1`:

```python title="Not run: a full histone digest at min_length=1 is a multi-megabyte recording"
digest = pymzlib.peptidoform.fragments("P68431", min_length=1)   # histone H3.1, every peptide
```

### The isoform cap truncates silently

Modification isoforms are enumerated combinatorially, and the count grows fast: histone H3.1's 49
peptides give 2,419 peptidoforms at two modifications and 6,580 at three. `max_isoforms` caps the
isoforms considered *per peptide position*, and when the cap binds mzLib discards the excess
**without an error**. On H3.1 at four modifications, mzLib's default cap keeps 8,644 of 12,808
peptidoforms, losing about a third. (These histone counts were measured on 2026-10-21. UniProt
revises its annotations, so a later digest can differ.)

`Digest` says when this happened. Albumin at the defaults is complete:

```pycon
>>> digest.truncated, digest.peptides_at_cap
(False, 0)

```

When `truncated` is true, raise the cap and digest again:

```python title="Not run: a histone digest at four modifications is a multi-megabyte recording"
digest = pymzlib.peptidoform.fragments("P68431", max_modifications=4)
if digest.truncated:
    print(f"{digest.peptides_at_cap} peptide positions hit the isoform cap")
    digest = pymzlib.peptidoform.fragments("P68431", max_modifications=4, max_isoforms=20000)
```

Check `digest.truncated` before you treat a peptidoform list as exhaustive.

## Protease naming: read this if you come from MaxQuant or Mascot

Protease names, and what each one cuts, are mzLib's: they come from its `proteases.tsv`. Two of
them look like names you know from MaxQuant and Mascot and mean the **reverse**:

| mzLib name | What it cuts | The MaxQuant and Mascot name for it |
|---|---|---|
| `"trypsin"` (mzLib's default) | After K/R, **including** before proline | `Trypsin/P` |
| `"trypsin|P"` | After K/R **except** before proline (the Keil rule) | plain `Trypsin` |

In mzLib's `"trypsin|P"`, the `P` names a residue that *prevents* the cut. In MaxQuant's
`Trypsin/P`, it names a rule that is *ignored*.

The difference is real. On serum albumin, without modifications, `"trypsin"` gives 202 distinct
peptides and `"trypsin|P"` gives 195. The counts differ by 7, but the lists differ by 37: 22
peptides are only in the first and 15 only in the second, because every K/R before a proline moves
a cleavage site. If you want the Keil rule (MaxQuant's plain `Trypsin`), pass
`protease="trypsin|P"`.

## Tuning the digest

Every option can be set. The common adjustments:

```python title="Not run: each of these digests is its own megabyte-scale recording"
# The bare sequence, no modifications: a clean control
pymzlib.peptidoform.fragments("P02768", modifications=False)

# ETD c and z• ions instead of HCD b and y
pymzlib.peptidoform.fragments("P02768", dissociation="ETD")

# The Keil rule: no cleavage before proline
pymzlib.peptidoform.fragments("P02768", protease="trypsin|P")

# A tighter length window and no missed cleavages
pymzlib.peptidoform.fragments("P02768", missed_cleavages=0, min_length=6, max_length=40)

# Only N-terminal fragments
pymzlib.peptidoform.fragments("P02768", terminus="N")
```

`modifications=False` is a clean control. It keeps UniProt's signal-peptide and propeptide
boundaries, which mzLib digests at, so the bare digest has the same distinct backbones and only
the modified variants go away.

`max_modifications` is the setting that most affects run time, because isoforms are combinatorial:
three modifications per peptide is much more work than two.

## Converting full sequences to Unimod or ProForma

A MetaMorpheus result names each modification the way mzLib's databases do:
`[UniProt:N-acetylserine on S]`, `[Common Variable:Oxidation on M]`. Most other tools want a Unimod
accession instead, such as `[UNIMOD:1]`. [`convert()`][pymzlib.peptidoform.convert] hands each
sequence to mzLib's `SequenceConversionService` and returns one row per input, in order. pyMzLib
maps no modification itself; every output is mzLib's.

This example uses `BottomUpExample.psmtsv`, a MetaMorpheus search result that ships with mzLib's
tests. Read its full sequences, then convert them:

```pycon
>>> psms = pymzlib.readers.read_results("BottomUpExample.psmtsv")
>>> unimod = pymzlib.peptidoform.convert(psms.columns["full_sequence"])
>>> unimod.source_format, unimod.target_format, unimod.mode
('mzLib', 'Unimod', 'ReturnNull')
>>> for row in unimod.sequences[:5]:
...     print(row.status, row.output)
converted YPIEH[UNIMOD:34]GIVTNWDDMEK
converted VAPEEHPVLLTEAPLNPK
converted AYHEQLSVAEITNAC[UNIMOD:4]FEPANQMVK
converted VEDMAELTC[UNIMOD:4]LNEASVLHNLK
converted YPIEH[UNIMOD:34]GIVTNWDDM[UNIMOD:35]EK
>>> unimod.converted_count == unimod.record_count
True

```

UniProt's tele-methylhistidine became `UNIMOD:34` (Methyl). Carbamidomethyl became `UNIMOD:4`, and
oxidation became `UNIMOD:35`. Unmodified sequences pass through unchanged. `unimod.columns` is
already the shape `pandas.DataFrame` wants, if you want a table instead of a loop.

Every result lists the notations mzLib has registered, so you never need to guess a name:

```pycon
>>> unimod.source_formats
['MassShift', 'Modomics', 'ProForma', 'mzLib']
>>> unimod.target_formats
['Chronologer', 'Essential', 'MassShift', 'ProForma', 'Unimod', 'mzLib']

```

### ProForma gives the same accessions

Ask for ProForma and this file gives the same accessions. The UniProt modification is `UNIMOD:34`
here too:

```pycon
>>> proforma = pymzlib.peptidoform.convert(psms.columns["full_sequence"], target="ProForma")
>>> proforma.outputs[0]
'YPIEH[UNIMOD:34]GIVTNWDDMEK'
>>> proforma.outputs == unimod.outputs
True

```

That is new with [mzLib#1403](https://github.com/smith-chem-wisc/mzLib/pull/1403).
Before it, mzLib's ProForma serializer looked modifications up only in MetaMorpheus's own list,
which has no UniProt entries, so it wrote `[UniProt:Tele-methylhistidine on H]` back by name. It
now falls back to every list mzLib loads, as the Unimod serializer does. pyMzLib never patched
around the gap, so the fix reached it, mzLibRust, mzLibR and MetaMorpheus at once. The `pro_forma`
column that [`readers.read_records()`](readers.md) computes for a `.psmtsv` written without one
comes from the same serializer, so it gained the accessions too.

The two targets still differ on a modification that no list mzLib loads knows. ProForma writes it
back under its own name and the row says `converted`; Unimod fails the row, or drops the
modification, depending on `mode`. In ProForma output, treat any bracket that is not a `UNIMOD:`
term as unresolved.

### When mzLib cannot convert a sequence

What happens to a modification the target cannot write depends on `mode`, which is mzLib's
`SequenceConversionHandlingMode`:

| `mode` | The row | `output` |
|---|---|---|
| `"ReturnNull"` (default) | `failed` | `None` |
| `"RemoveIncompatibleElements"` | `converted_with_warnings` | the sequence without that modification |
| `"UsePrimarySequence"` | `converted_with_warnings` | the sequence without that modification |
| `"ThrowException"` | no rows: the call raises `UsageError` naming the first such sequence | none |

`status` is mzLib's own verdict:

- `converted`: mzLib returned an output and recorded nothing against it.
- `converted_with_warnings`: mzLib returned an output but noted something, such as a dropped
  modification or a skipped character.
- `failed`: mzLib returned nothing.

These four sequences include one modification that has no Unimod accession:

```pycon
>>> result = pymzlib.peptidoform.convert([
...     "[UniProt:N-acetylserine on S]SEQK",
...     "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
...     "PEPM[Common Variable:Oxidation on M]K",
...     "PEPK[Made Up:Not a modification on K]R",
... ])
>>> result.outputs
['[UNIMOD:1]SEQK', 'PEPK[UNIMOD:36]R', 'PEPM[UNIMOD:35]K', None]
>>> (bad,) = result.not_converted
>>> bad.input, bad.status, bad.failure_reason
('PEPK[Made Up:Not a modification on K]R', 'failed', None)
>>> bad.incompatible_items
['Made Up:Not a modification on K @3(K)']

```

`failure_reason` is `None` here. Under `ReturnNull`, mzLib's Unimod serializer records the
incompatible modification but no reason code, so read `incompatible_items` to see what failed.

**Split ambiguous sequences first.** When MetaMorpheus cannot tell candidates apart, it joins their
full sequences with `|`. mzLib's parser does not refuse this. It skips each `|` with a warning and
joins the candidates into one sequence. The row is `converted_with_warnings`, but the output is not
a real peptide. Split on `|` before you convert.

The full list of parameters, fields and caveats is in the
[`peptidoform convert` reference](../reference/peptidoform.md#convert).

## Errors you might hit

| Exception | Means |
|---|---|
| `UsageError` | `fragments()`: the accession, protease, dissociation type or terminus is not recognised, or a number is negative; raised before any network call. `convert()`: no sequences, a blank one, a format name mzLib has not registered (the message lists them), an unknown `mode`, or, under `mode="ThrowException"`, a sequence mzLib could not convert. |
| `ServiceUnavailableError` | UniProt was down, rate-limiting, or timed out. **Not your bug**: retry later. |
| `BridgeError` | UniProt answered, but something about the request was rejected. |

A malformed accession is refused before anything is fetched:

```pycon
>>> pymzlib.peptidoform.fragments("not-an-accession")      # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
pymzlib._bridge.UsageError: 'not-an-accession' is not a valid UniProtKB accession. ...

```

`ServiceUnavailableError` is a subclass of `BridgeError`, so catching `BridgeError` catches both.
Separate them when you want to retry an outage but report everything else:

```python title="Not run: shows how to handle a UniProt outage, which a recording cannot produce"
try:
    digest = pymzlib.peptidoform.fragments("P02768")
except pymzlib.ServiceUnavailableError:
    ...   # UniProt's problem; worth a retry
except pymzlib.BridgeError as e:
    print(f"{e.error_type}: {e}")
```

[Every error each verb can raise](../errors.md) is listed on one page.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
sources the digest reads:

--8<-- "docs/reference/_generated/cite.peptidoforms.md"
