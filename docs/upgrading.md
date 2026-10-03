# Upgrading: what changed in your results

pyMzLib bundles one mzLib build, and every pyMzLib release moves it forward. Most of what a new
mzLib brings is new capability. This page is the rest: **changes that give you a different answer
from the same file and the same call**. Read the section for the release you are moving to before
you compare numbers across the upgrade.

Each entry names the mzLib pull request it comes from, so you can read the change itself. The
[changelog](https://github.com/smith-chem-wisc/pyMzLib/blob/main/CHANGELOG.md) has everything else.
Which mzLib a bridge was built from is in `pymzlib.bridge_version()["mzlib"]`.

## Next release (built from mzLib 1.0.593)

### Files that used to be refused now read

- **MetaMorpheus protein-group tables written without quantification**
  ([mzLib#1365](https://github.com/smith-chem-wisc/mzLib/pull/1365)). `AllProteinGroups.tsv` and
  each file's `<file>_ProteinGroups.tsv` raised `Tsv file type not supported`. They now read as
  `MetaMorpheusQuantifiedProteinGroups`, with spectral counts and **no intensities**: `intensity`
  is in `absent_fields`, meaning no basis, not zero.
- **RNA quantification tables** ([mzLib#1388](https://github.com/smith-chem-wisc/mzLib/pull/1388)).
  `AllQuantifiedTranscriptGroups.tsv` and `AllQuantifiedOligos.tsv` are new file types, so
  `formats()` lists more entries than before. Code that asserts a count of formats needs updating.

### Same file, different value

- **Four modifications write their Unimod accession in `pro_forma`**
  ([mzLib#1328](https://github.com/smith-chem-wisc/mzLib/pull/1328)). In `read_records()` on a
  MetaMorpheus `.psmtsv` or `.osmtsv`, `GG (Ubiquitination Site)` is now `[UNIMOD:121]`, both
  `Myristoylation` entries `[UNIMOD:45]` and `EQIGG` `[UNIMOD:846]`, where they were written by name.
  **Masses are unchanged.** A string comparison of `pro_forma` across the upgrade will see these
  peptides as different; a pipeline keyed on Unimod will now see ubiquitination sites it silently
  dropped before.

### Same failure, different error

- **A read fault on an existing `.mzid`** ([mzLib#1362](https://github.com/smith-chem-wisc/mzLib/pull/1362))
  is now a `BridgeError` with `error_type` `MzLibException`, naming the file, where it was an
  `IOException`. A missing file is still a `UsageError`. A plain `.mzid` is now streamed from disk,
  so reading a large one needs far less memory.
- **PRIDE errors name the file and the host, never the URL**
  ([mzLib#1350](https://github.com/smith-chem-wisc/mzLib/pull/1350)), so a reviewer token in a
  query string cannot reach a log. Every transport failure is still a `ServiceUnavailableError`.
  Code that parsed a URL out of an error message will find none.

## pyMzLib 0.2.0 (built from mzLib 1.0.592)

### Files that used to be refused now read

- **mzIdentML** (`.mzid`, and `.mzid.gz` without unpacking it,
  [mzLib#1313](https://github.com/smith-chem-wisc/mzLib/pull/1313)), and MetaMorpheus's
  `AllQuantifiedProteinGroups.tsv` and `AllQuantifiedPeptides.tsv`
  ([mzLib#1347](https://github.com/smith-chem-wisc/mzLib/pull/1347)).
- **Pytheas match output** ([mzLib#1277](https://github.com/smith-chem-wisc/mzLib/pull/1277)). A
  `.txt` whose first five lines carry `#theoretical_digest` reads as Pytheas; any other `.txt` is
  still Crux.
- **Current FlashLFQ and MetaMorpheus `QuantifiedPeaks.tsv`**
  ([mzLib#1345](https://github.com/smith-chem-wisc/mzLib/pull/1345)), which have no `MBR Score`
  column, read instead of failing with `HeaderValidationException`.

### Same file, different value

- **`pro_forma` is filled on older MetaMorpheus files**
  ([mzLib#1346](https://github.com/smith-chem-wisc/mzLib/pull/1346)). MetaMorpheus 1.1.11 and
  earlier wrote no ProForma column, so `read_records()` gave `None`; mzLib now converts
  `full_sequence`. Ambiguous (`|`-joined) or unconvertible rows are still `None`. Reading is slower:
  a 295 MB `AllPSMs.psmtsv` took 21% longer.
- **`mbr_score` is `None` where it was `0`** in a FlashLFQ `QuantifiedPeaks.tsv`
  ([mzLib#1345](https://github.com/smith-chem-wisc/mzLib/pull/1345)), for every peak without a
  score. A `0` used to say "scored, and scored nothing". Seven columns also join the record:
  `organism`, `peak_fwhm`, `peak_fwhm_status`, `pip_q_value`, `pip_pep`, `decoy_peptide` and
  `random_rt`.
- **MGF `one_based_precursor_scan_number`** is no longer always `None`
  ([mzLib#1227](https://github.com/smith-chem-wisc/mzLib/pull/1227)): mzLib reads it from the
  `PRECURSORSCAN` header it writes. It is still `None` on MGF files mzLib did not write.
- **FlashLFQ with match-between-runs is reproducible across thread counts**
  ([mzLib#1155](https://github.com/smith-chem-wisc/mzLib/pull/1155)). Before, more than one thread
  could give different peptide and protein numbers for identical inputs.

## Earlier releases

- **MSFragger retention times are minutes** from pyMzLib 0.1.0
  ([mzLib#1116](https://github.com/smith-chem-wisc/mzLib/pull/1116)). Before, an MSFragger
  `psm.tsv` passed seconds through, and FlashLFQ quantified it in the wrong window.
