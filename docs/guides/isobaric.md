# Isobaric kits

An isobaric experiment - TMT, TMTpro, iTRAQ, DiLeu - is read out of the low-m/z **reporter ions**:
one ion per channel, sometimes only 6 mDa apart. Mislabel one channel, or shift its m/z, and every
sample after it is assigned to the wrong condition. Nothing about the result looks wrong.

So pyMzLib keeps no table of its own. `pymzlib.isobaric.kits()` returns **mzLib's**, the one
MetaMorpheus quantifies TMT searches with:

| You want | Call | mzLib |
|---|---|---|
| Every kit mzLib knows, with every channel | [`kits()`](#every-kit) | `IsobaricMassTag`, `IsobaricMassTagType` ([#1375][1375]) |
| One kit, by its name or a MetaMorpheus modification name | [`kits("TMT18")`](#one-kit-tmtpro-18-plex) | `IsobaricMassTag.TryGetTagType` |
| A channel's reporter-ion m/z and matching window | [`.channels`](#what-the-numbers-are) | `ReporterIonMzs`, `ReporterIonMzRanges` |

Every example on this page is executed in CI against output recorded from the real bridge.

## Every kit

```pycon
>>> import pymzlib
>>> catalogue = pymzlib.isobaric.kits()
>>> for kit in catalogue.kits:
...     print(f"{kit.name:8} {kit.channel_count:3}  {kit.labels[0]} ... {kit.labels[-1]}")
TMT6       6  126 ... 131
TMT10     10  126 ... 131N
TMT11     11  126 ... 131C
TMT16     16  126 ... 134N
TMT18     18  126 ... 135N
iTRAQ4     4  114 ... 117
iTRAQ8     8  113 ... 121
diLeu4     4  115 ... 118
diLeu12   12  115a ... 118d

```

`catalogue.columns` is the same data as one long table, one row per (kit, channel), ready for
`pandas.DataFrame(catalogue.columns)`.

**`TMT16` and `TMT18` are TMTpro.** TMT16 is the lowest sixteen channels of the 18-plex set and has
no modification entry of its own, so its m/z are TMT18's first sixteen. **iTRAQ 8-plex has no 120
channel**: the phenylalanine immonium ion sits at m/z 120.081, so the eighth reagent is 121.

```pycon
>>> catalogue.kit_named("iTRAQ8").labels
['113', '114', '115', '116', '117', '118', '119', '121']

```

## One kit: TMTpro 18-plex

Ask for a kit the way MetaMorpheus names it. The match is on the **whole name**, case-insensitive,
optionally followed by `" on <motif>"`: `"TMT18"`, `"tmt10"`, `"TMT6-plex"` and
`"iTRAQ-4plex on K"` all work. A substring never does, so `"TMT10plex"` is refused with a list of
the kits that exist, rather than quietly resolved to whichever kit's name it happens to contain.

```pycon
>>> tmtpro = pymzlib.isobaric.kits("TMT18")
>>> kit = tmtpro.kits[0]
>>> kit.name, kit.channel_count
('TMT18', 18)
>>> for channel in kit.channels[:4]:
...     print(f"{channel.index:2}  {channel.label:5} {channel.reporter_ion_mz:.5f}")
 0  126   126.12773
 1  127N  127.12476
 2  127C  127.13108
 3  128N  128.12812

```

127N and 127C are 6.3 mDa apart. That is why the next section matters.

## What the numbers are

```pycon
>>> last = kit.channels[-1]
>>> last.label, round(last.reporter_ion_mz, 5), round(last.mz_min, 5), round(last.mz_max, 5)
('135N', 135.1516, 135.1486, 135.1546)
>>> tmtpro.absolute_tolerance
0.003

```

- **`reporter_ion_mz` is theoretical, at charge 1.** mzLib does not type m/z values in: each is a
  `DI HCD` diagnostic-ion line of the kit's `Multiplex Label` modification in its embedded
  `TMT.txt`, plus one proton, sorted ascending and paired with the labels by position. It is not a
  calibrated or observed value.
- **`mz_min` and `mz_max` are mzLib's matching window**, `reporter_ion_mz` ± `absolute_tolerance`
  (0.003 Da). When mzLib reads reporter intensities out of a spectrum, it takes the most intense
  peak inside each window. mzLib takes the 0.003 Da figure from Li *et al.* (2021), cited below.
- **`index` is 0-based** and ascends with m/z; `label` at that index names the channel whose
  reporter ion it is.

## What this does not do

- **Read reporter intensities out of your spectra.** mzLib can
  (`IsobaricMassTag.GetReporterIonIntensities`), but no wire verb projects it yet. Read the spectra
  with [`pymzlib.readers.read_spectra`](readers.md) and match against `mz_min`/`mz_max` only for a
  quick look; for a quantification, use MetaMorpheus.
- **Build an isobaric design.** [`pymzlib.sdrf.design()`](sdrf.md#turn-an-sdrf-into-a-quantification-design)
  is label-free only.

## References

- Li J., Cai Z., Bomgarden R.D., *et al.* TMTpro-18plex: The Expanded and Complete Set of TMTpro
  Reagents for Sample Multiplexing. *Journal of Proteome Research* **20**, 2964-2972 (2021).
  [doi:10.1021/acs.jproteome.1c00168](https://doi.org/10.1021/acs.jproteome.1c00168) - the 18-plex
  channel set, and the source of the 0.003 Da matching tolerance mzLib uses.
- mzLib [#1375][1375] moved `IsobaricMassTag` into `Omics.Modifications`, named every channel, and
  derived each m/z from `TMT.txt`.

[1375]: https://github.com/smith-chem-wisc/mzLib/pull/1375
