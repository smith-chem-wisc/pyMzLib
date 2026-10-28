# pyMzLib

**mzLib for Python** — mass spectrometry and proteomics from the
[mzLib](https://github.com/smith-chem-wisc/mzLib) C# library, with no .NET to install.

```bash
pip install mzlib
```

Install it as `mzlib`, import it as `pymzlib` — the same split as `pip install scikit-learn` /
`import sklearn`. The project and this repository are pyMzLib.

```python
import pymzlib

# What's in a PRIDE Archive project?
files = pymzlib.pride.list_files("PXD000001")
print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")

# Download just its Thermo .raw file (220 MB), then read the first five MS2 scans, peaks included.
# The same call reads mzML, Bruker .d, timsTOF .d, MGF and msalign.
raw_files = pymzlib.pride.download("PXD000001", "downloads", category="RAW")
scans = pymzlib.readers.read_spectra(raw_files[0], ms_order=2, limit=5, peaks=True)
print(scans.scan_count, scans.columns["selected_ion_mz"])

# Digest an annotated UniProt protein and fragment its peptides
digest = pymzlib.peptidoform.fragments("P02768")
print(digest.modification_census.explain())
```

That's the whole installation. No .NET runtime, no configuration, and **no third-party Python
dependencies** — so pyMzLib cannot conflict with anything already in your environment.

📖 **[Documentation](https://smith-chem-wisc.github.io/pyMzLib/)** ·
🐛 **[Issues](https://github.com/smith-chem-wisc/pyMzLib/issues)** ·
🧪 **[mzLib](https://github.com/smith-chem-wisc/mzLib)**

## What's covered

Coverage is deliberately partial and grows by demand, the way
[pyOpenMS](https://pyopenms.readthedocs.io/) grew.

| Area | Status |
|---|---|
| [PRIDE Archive](https://smith-chem-wisc.github.io/pyMzLib/guides/pride/) — search for projects by keyword, list a project's files, filtered download, the MD5 and exact size to check each download against | ✅ |
| [Peptidoforms](https://smith-chem-wisc.github.io/pyMzLib/guides/peptidoforms/) — digest an annotated protein, apply its modifications, fragment every peptide; convert full sequences to Unimod or ProForma | ✅ |
| [Quantification](https://smith-chem-wisc.github.io/pyMzLib/guides/flashlfq/) — FlashLFQ label-free quant with match-between-runs, and median-polish protein roll-up | ✅ |
| [Readers](https://smith-chem-wisc.github.io/pyMzLib/guides/readers/) — read spectra from **mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign; identify and read every file type mzLib knows, search results included; many files in one call; MetaMorpheus protein groups, FlashLFQ peptides and PTM site occupancy as long tables, RNA transcript groups and oligos included | ✅ |
| [SDRF experimental design](https://smith-chem-wisc.github.io/pyMzLib/guides/sdrf/) — read an SDRF-Proteomics file, pool several experiments into one analysis table; validate, lint and assess them; one row per sample with ages in years; the label-free design FlashLFQ and MetaMorpheus take, and the isobaric (TMT) channel design MetaMorpheus takes, or every reason either was refused; restrict an SDRF to one search's files | ✅ |
| [Protein databases](https://smith-chem-wisc.github.io/pyMzLib/guides/proteins/) — organism, taxon, GO terms and Ensembl genes per accession from a UniProt XML or FASTA; resolve proteins to Ensembl genes against a pinned release; peptide uniqueness with I = L; Gene Ontology on MetaMorpheus protein groups, every member kept, against a pinned GO release; align two sequences and map residues between isoforms or orthologs; read accessions into entry, isoform and variants; Ensembl Compara orthologs across species, as a Parquet snapshot; map a residue across species through Compara's gene-tree alignment (`proteins`) | ✅ |
| [Isobaric kits](https://smith-chem-wisc.github.io/pyMzLib/guides/isobaric/) — every TMT, TMTpro, iTRAQ and DiLeu channel with its reporter-ion m/z, matching window and SDRF label, as MetaMorpheus quantifies them (`isobaric`) | ✅ |
| [Spectral libraries](https://smith-chem-wisc.github.io/pyMzLib/guides/library/) — read an mzLib `.msl` spectral library: every precursor with its m/z, iRT or minutes, ion mobility and target/decoy/entrapment label, and its fragment ions on request; libraries too large for memory index-only; calibrate a run's minutes onto the library's iRT and back (`library`) | ✅ |
| [Differential abundance](https://smith-chem-wisc.github.io/pyMzLib/guides/stats/) — limma's moderated t-test (checked against limma to 1e-8), for coefficients and for contrasts with confidence intervals, Benjamini-Hochberg, random-effects meta-analysis and p-value combination (Fisher, Stouffer, signed) across studies, Spearman correlation with a trait (partial too), logistic regression of detection, mixed models with a random intercept per group, and entrapment checks of a search's FDR, with no R (`stats`) | ✅ |
| [Differential design and observations](https://smith-chem-wisc.github.io/pyMzLib/guides/diff/) — a sample sheet to mzLib's design matrix per stratum (equal to R's `model.matrix`), with each default contrast's estimability and the design and contrasts files `stats.fit` reads; and every peptide in every sample of a FlashLFQ or MetaMorpheus table as log2 intensity or no value with the reason, with mzLib's shared-peptide normalization (`diff`) | ✅ |
| [PTM occupancy and PTM-QTL](https://smith-chem-wisc.github.io/pyMzLib/guides/ptmqtl/) — site occupancy per run or per sample, every cell with its censoring state and peptidoform counts; PTM pairs on one molecule or co-varying across runs, with a batch or each run's depth removed and one-form runs ranked as bounds; the effect of a per-run trait on a site, replicates kept apart, in natural-log or log2 odds (`ptmqtl`) | ✅ |
| Everything else in mzLib | not yet — [request it](https://github.com/smith-chem-wisc/pyMzLib/issues/new/choose) |

If there's something in mzLib you want from Python, opening an issue is genuinely the fastest
path. Requests are how we decide what to cover next.

## How it works

mzLib runs as a self-contained executable bundled inside the package; pyMzLib starts it, speaks
JSON to it, and hands you ordinary Python objects. The mechanism is an implementation detail and
may change — the public API will not.

Worth knowing: calls carry a small fixed startup cost, so pyMzLib is built for coarse-grained
operations rather than tight per-spectrum loops. The
[architecture page](https://smith-chem-wisc.github.io/pyMzLib/design/architecture/) explains why,
and the [design decisions](https://smith-chem-wisc.github.io/pyMzLib/design/decisions/) record
what was traded for what.

## Contributing

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) and the
[capability recipe](https://smith-chem-wisc.github.io/pyMzLib/contributing/adding-a-capability/),
which walks through exposing a new piece of mzLib end to end.

By participating you agree to the [Code of Conduct](CODE_OF_CONDUCT.md).

## License

[LGPL-3.0-or-later](LICENSE), matching mzLib.

pyMzLib distributes mzLib in compiled form inside its wheels. Under the LGPL you may modify the
mzLib portion and rebuild: the complete source of this project is here, the mzLib commit each
release was built from is recorded in [`code/PINNED.md`](code/PINNED.md) and reported at runtime
by `pymzlib.bridge_version()`, and the build is fully scripted in
[`pkg/build/`](pkg/build/). See
[building from source](https://smith-chem-wisc.github.io/pyMzLib/contributing/building/).

## Citing

pyMzLib does not yet have its own paper. Please cite mzLib and MetaMorpheus for the underlying
science; see [CITATION.cff](CITATION.cff).

---

Developed in the [Smith lab](https://smith.chem.wisc.edu/), University of Wisconsin–Madison.
