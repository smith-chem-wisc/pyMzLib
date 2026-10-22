# Spectral libraries

A spectral library says, for each **precursor** (a peptide or proteoform at one charge), where a
search should expect to find it: its m/z, its retention time, its ion mobility, and the fragment
ions its MS2 spectrum should hold, with their relative intensities. MetaMorpheus and mzLib write
libraries in mzLib's binary `.msl` format, predicted (Prosit, through Koina) or built from a search.
pyMzLib reads them with mzLib's own reader, so a library reads in Python exactly as MetaMorpheus
reads it.

| You want to | Call | mzLib does it with |
|---|---|---|
| See every precursor in a library | [`read()`](#every-precursor) | `MslFileTypeHandler.Open`, `MslLibrary.GetEntry` |
| Know whether a retention time is iRT or minutes | [`retention_time_unit`](#retention-time-irt-or-minutes-row-by-row) | `MslLibraryEntry.RtIsCalibrated` |
| Count targets, decoys and entrapment precursors | [`decoy_contaminant_target`](#targets-decoys-and-entrapment) | `DecoyContaminantTargetLabel.For` |
| See a precursor's fragment ions | [`read(fragments=True)`](#fragment-ions) | `MslFragmentIon` |
| Read a library too large for memory | [`read(index_only_threshold=...)`](#large-libraries-index-only) | `MslLibrary.LoadIndexOnly` |
| Read several libraries into one table | [`read_many()`](#many-libraries-at-once) | the same reader, in one bridge process |

Every `>>>` example on this page runs in CI against output recorded from the real bridge. Two small
libraries are used. `myPrositLib.msl` is mzLib's own test library `myPrositLib.msp` (five
Prosit-predicted precursors) converted to `.msl` by mzLib. `synthetic_library.msl` holds six
precursors written to cover every field: a target, its decoy, an entrapment target whose retention
time is in minutes, an entrapment decoy, a predicted precursor with nothing measured, and a
proteoform. Both are in the pyMzLib repository under `pkg/python/tests/fixtures/library/`.

## Every precursor

```pycon
>>> import pymzlib
>>> library = pymzlib.library.read("synthetic_library.msl")
>>> library.record_count, library.target_count, library.decoy_count, library.proteoform_count
(6, 4, 2, 1)
>>> for row in library.records:
...     print(row["precursor_index"], row["full_sequence"], row["charge_state"],
...           round(row["precursor_mz"], 4), row["decoy_contaminant_target"])
0 PEPTIDEK 2 465.74 T
1 KEDITPEP 2 465.7401 D
2 PEPTLDEK 2 458.7462 ET
3 KEDLTPEP 2 458.7463 ED
4 ELVISLIVESK 3 406.242 T
5 MAHVRSGKHQIPAEALQ 2 952.493 T

```

One row per precursor, **in the library's stored order**. `precursor_index` is the position in the
file and the key mzLib's `MslLibrary.GetEntry` takes. It is not m/z order: mzLib's own
`GetAllEntries` enumerates by ascending m/z, so sort on `precursor_mz` if that is the order you want.

`library.columns` maps each column name to its values, the shape `pandas.DataFrame(library.columns)`
takes. `library.records` is the same data as one dict per row.

### What each column is

The names are mzLib's `MslLibraryEntry` properties. Where mzLib stores a value it documents as
"not available", the column holds `None`, and the table says which:

| column | what it is | `None` means |
|---|---|---|
| `precursor_index` | 0-based position in the library | never `None` |
| `full_sequence`, `base_sequence` | mzLib bracket notation, e.g. `PEPTM[Common Variable:Oxidation on M]IDE`; the bare residues | never `None` |
| `precursor_mz` | precursor m/z, (M + zH)/z | never `None` |
| `charge_state` | precursor charge | never `None` |
| `retention_time`, `retention_time_unit` | the time, and whether it is `"iRT"` or `"minutes"` | never `None` |
| `ion_mobility` | 1/K0, in V·s/cm² | mzLib stored 0.0, "not available" |
| `is_decoy`, `is_entrapment`, `decoy_contaminant_target` | target/decoy/entrapment flags, and mzLib's label `T`, `D`, `ET` or `ED` | never `None` |
| `protein_accession`, `protein_name`, `gene_name` | the source protein | mzLib stored an empty string: no protein |
| `is_proteotypic` | the library writer's proteotypic flag | never `None` |
| `q_value` | the library's q-value for the precursor | mzLib stored NaN: no q-value, as in a predicted library |
| `elution_group_id` | shared by every precursor of one base sequence | never `None` |
| `molecule_type` | `Peptide`, `Proteoform`, `Oligonucleotide` or `Glycopeptide` | never `None` |
| `source_type` | `Predicted`, `Empirical`, or `EmpiricalRefined` (predicted intensities, measured time or mobility) | never `None` |
| `dissociation_type` | e.g. `HCD`; `Unknown` when the library does not say | never `None` |
| `nce` | normalized collision energy, a whole number | mzLib stored 0: unknown |

The predicted precursor shows the "not available" values, and the empirical target shows them
filled in:

```pycon
>>> for name in ("ion_mobility", "q_value", "nce", "protein_accession", "source_type"):
...     print(f"{name:18} {library.columns[name][0]!r:22} {library.columns[name][4]!r}")
ion_mobility       0.949999988079071      None
q_value            0.001                  None
nce                28                     None
protein_accession  'P02768'               None
source_type        'Empirical'            'Predicted'

```

A whole column can have no source at all. `myPrositLib.msl` has no protein table and no ion
mobility, so those columns are named in `absent_fields` and are `None` in every row:

```pycon
>>> prosit = pymzlib.library.read("myPrositLib.msl")
>>> prosit.absent_fields
['ion_mobility', 'protein_accession', 'protein_name', 'gene_name']
>>> prosit.columns["nce"], prosit.columns["dissociation_type"][0]
([None, None, None, None, None], 'Unknown')

```

The `.msp` it was converted from records `Collision_energy=30.0`, yet `nce` is `None`. That is not
a reading error. mzLib converts a text library through `LibrarySpectrum`, which has no field for
NCE, dissociation type, q-value, ion mobility or protein, so `MslLibraryEntry.FromLibrarySpectrum`
writes mzLib's "not available" value for each. A converted library has none of them.

## Retention time: iRT or minutes, row by row

A library stores each precursor's retention time either as an **indexed retention time** (iRT, a
unitless scale anchored on reference peptides) or as **minutes** in one particular run. Which one
is a flag on each precursor (`MslLibraryEntry.RtIsCalibrated`), and one library can hold both. So
the unit is a column, not a property of the file:

```pycon
>>> for row in library.records[:4]:
...     print(row["full_sequence"], round(row["retention_time"], 2), row["retention_time_unit"])
PEPTIDEK 35.2 iRT
KEDITPEP 35.2 iRT
PEPTLDEK 21.4 minutes
KEDLTPEP 33 iRT

```

**Never plot or compare `retention_time` without `retention_time_unit`.** pyMzLib converts nothing:
iRT becomes minutes only through a calibration against a run, which is not in the library. Keep
the rows of one unit before you use the numbers:

```pycon
>>> irt = [r for r in library.records if r["retention_time_unit"] == "iRT"]
>>> len(irt), len(library.records)
(5, 6)

```

iRT can be negative; Prosit predicts `-15.77` for one of the `myPrositLib.msl` precursors. A
`retention_time` of exactly 0 deserves suspicion: mzLib's conversion from a text library writes 0
when the spectrum had no retention time, and documents no sentinel for it, so 0 is passed through
as a value.

## Targets, decoys and entrapment

**The library is unfiltered.** Decoys and entrapment precursors are rows like any other.
`decoy_contaminant_target` is mzLib's own label, written by `DecoyContaminantTargetLabel.For`: `T`
target, `D` decoy, `ET` entrapment target, `ED` entrapment decoy. The `.msl` format has no
contaminant flag, so `C` never appears.

```pycon
>>> from collections import Counter
>>> Counter(library.columns["decoy_contaminant_target"])
Counter({'T': 3, 'D': 1, 'ET': 1, 'ED': 1})
>>> targets = [r["full_sequence"] for r in library.records if not r["is_decoy"]]
>>> targets
['PEPTIDEK', 'PEPTLDEK', 'ELVISLIVESK', 'MAHVRSGKHQIPAEALQ']

```

Entrapment is a target for the search: it competes with the real targets, so `target_count`
includes the entrapment targets and `decoy_count` includes the entrapment decoys. To count only
the real targets, keep `decoy_contaminant_target == "T"`.

## Fragment ions

With `fragments=True` the table is **one row per fragment ion**, the precursor's columns repeated
on each row, then the fragment's own. `limit` and `offset` still count precursors:

```pycon
>>> ions = pymzlib.library.read("synthetic_library.msl", fragments=True, limit=1)
>>> ions.returned_count, ions.row_count
(1, 6)
>>> for row in ions.records:
...     print(f'{row["annotation"]:11} {row["fragment_mz"]:9.4f} {row["fragment_intensity"]:5.2f} '
...           f'{row["residue_position"]!s:5} {row["neutral_loss"]}')
y1           147.1128  0.40 8     0
b2           227.1026  0.30 2     0
y2           276.1554  1.00 7     0
bIy[2-5]     326.1710  0.10 2     0
y3-18.0106   373.2082  0.20 6     -18.010565
y4-43.0058   461.2242  0.15 None  -43.005814

```

What the numbers are:

- **`fragment_intensity` is relative within its precursor, 0 to 1.** mzLib's writer scales each
  precursor's most intense fragment to 1, so intensities are comparable within a precursor and
  not across precursors.
- **`neutral_loss` is in Da, negative for a loss, and 0 means no loss.** It is 0, not `None`,
  because 0 Da is the value of "no loss", so `fragment_mz` and `neutral_loss` stay usable together.
- **`annotation` is mzLib's library rendering** (`MslFragmentIon.Annotation`): the loss signed to
  four decimals and, for a charge above 1, the charge appended (`y4^2+`). It differs from
  `Product.Annotation` elsewhere in mzLib, so do not match the two strings against each other.
- **An internal ion** (`bIy[2-5]`) has `is_internal` true, its C-terminal type in
  `secondary_product_type` and its end residue in `secondary_fragment_number`. Both are `None` on
  terminal ions, where mzLib stores 0 and null.
- **`residue_position` is `None` for a fragment with a custom neutral loss**, like the
  `-43.0058` above. The `.msl` format reuses that field to index its table of custom losses, and
  mzLib reads it back as 0. The loss itself survives. Named losses (water, ammonia, phosphate)
  keep their position.
- **`fragment_count`** on each row is how many fragment ions the precursor has. A precursor with
  none is still one row, with `fragment_count` 0 and every fragment cell `None`.

Fragments are stored, and returned, in m/z order: mzLib's writer sorts them. The table is a
spectrum to plot as is:

```python title="Not run: needs matplotlib, which pyMzLib does not depend on"
import matplotlib.pyplot as plt

ions = pymzlib.library.read("my_library.msl", fragments=True, offset=120, limit=1)
plt.vlines(ions.columns["fragment_mz"], 0, ions.columns["fragment_intensity"])
for mz, height, label in zip(ions.columns["fragment_mz"], ions.columns["fragment_intensity"],
                             ions.columns["annotation"]):
    plt.annotate(label, (mz, height))
plt.xlabel("m/z")
plt.ylabel("relative intensity")
plt.show()
```

A library has tens of fragment ions per precursor. For a whole library, write the table to disk
with `out=` (a `.tsv`) and read that, rather than returning it as one document.

## Large libraries: index-only

mzLib opens a library one of two ways. **Loaded whole**, every precursor and every fragment is in
memory. **Index-only**, the precursors are in memory and each precursor's fragment block is read
from disk when it is asked for. mzLib chooses by file size (`MslFileTypeHandler.Open`): a library
at or above `MslFileTypeHandler.DefaultIndexOnlyThresholdBytes`, 1 GiB in the mzLib this pyMzLib is
built from, opens index-only. pyMzLib sends nothing unless you pass `index_only_threshold`, a size
in bytes; `0` opens every library index-only:

```pycon
>>> small = pymzlib.library.read("synthetic_library.msl", index_only_threshold=0, limit=2)
>>> small.is_index_only, small.index_only_threshold
(True, 0)
>>> small.columns["full_sequence"]
['PEPTIDEK', 'KEDITPEP']

```

**The rows are the same in both modes**; only memory and speed differ. Two things to know:

- A **compressed** library (`is_compressed`) always loads whole: its fragment offsets point into
  the decompressed data, not into the file, so mzLib decompresses it all and `is_index_only` stays
  `False` whatever the threshold.
- Index-only **with fragments** reads every fragment block from disk twice: once to count the rows,
  once to write them. mzLib keeps each precursor's stored fragment count to itself, so there is no
  cheaper way to size the table. Without fragments, index-only reads no fragment block at all,
  which is the fast way to list the precursors of a very large library.

## Many libraries at once

`read_many()` reads a list of libraries into one long table, in one bridge process. Its first two
columns say where each row came from (`source_index`, `source_path`), and `files` holds one report
per library:

```pycon
>>> batch = pymzlib.library.read_many(["synthetic_library.msl", "myPrositLib.msl"])
>>> batch.record_count, [f.record_count for f in batch.files]
(11, [6, 5])
>>> batch.columns["source_index"]
[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
>>> [f.absent_fields for f in batch.files]
[[], ['ion_mobility', 'protein_accession', 'protein_name', 'gene_name']]

```

`threads` reads that many libraries at once. The answer is identical at any value; the default is
1 because a library loaded whole is held in memory while its rows are written. `on_error="skip"`
records an unreadable library in its report's `error` and reads the rest.

## Traps

- **Stored order, not m/z order.** Rows follow `precursor_index`. Sort on `precursor_mz` yourself.
- **The extra digits are the float.** `precursor_mz`, `retention_time` and `ion_mobility` are
  stored as 32-bit floats and widened by mzLib, so 465.74 reads back as `465.739990234375`. Round
  them for display; they carry no more precision than about seven significant digits.
- **`nce` is a whole number.** The format stores NCE × 10 and mzLib reads it back with integer
  division, so an NCE of 27.5 reads as 27.
- **Converted libraries are thin.** A library converted from `.msp` has no NCE, dissociation type,
  q-value, ion mobility or protein (see [above](#what-each-column-is)).
- **A broken file is an error, not a short table.** mzLib checks the whole file (magic number,
  version, footer, precursor count and a CRC-32 checksum) before it returns a single precursor. A
  damaged library raises `BridgeError` with `error_type` `InvalidDataException`, and a file that
  is not a library raises it with `FormatException`. A file whose name does not end in `.msl` is a
  `UsageError` before anything is opened.

## What this does not do

- **Read a text library** (`.msp`, pDeep, ms2pip). mzLib can turn one into a library in memory
  (`MslConverter.FromMspFile`), but no wire verb projects it yet.
- **Write a library.** mzLib writes them (`MslLibrary.Save`); no wire verb projects that yet.
- **Query by m/z or retention-time window.** mzLib's `QueryMzWindow` and `QueryWindow` are what a
  search engine calls. To find a precursor here, read the table and filter it.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)):

--8<-- "docs/reference/_generated/cite.library.md"

Further sources:

- mzLib [#1383][1383] added the entrapment flag and made `rt_is_calibrated` work, [#1421][1421]
  added `GetEntry` without fragments, and [#1420][1420] made libraries over 2 GB open.

[1383]: https://github.com/smith-chem-wisc/mzLib/pull/1383
[1420]: https://github.com/smith-chem-wisc/mzLib/pull/1420
[1421]: https://github.com/smith-chem-wisc/mzLib/pull/1421
