# pyMzLib

**mzLib for Python.** [mzLib](https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry
and proteomics library written in C#, developed in the Smith lab at UW–Madison. pyMzLib makes it
callable from Python.

```bash
pip install mzlib
```

Install it as `mzlib`, import it as `pymzlib` — the same split as `pip install scikit-learn` /
`import sklearn`.

That is the whole installation. There is no .NET to install, no runtime to configure, and no
third-party Python package to reconcile with the rest of your environment — pyMzLib declares
**zero runtime dependencies** and carries everything it needs inside the wheel.

## Quick start

```python
import pymzlib

# What's in a PRIDE Archive project?
files = pymzlib.pride.list_files("PXD000001")
print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")

for f in files:
    print(f"{f.category:8s} {f.size_mb:9.1f} MB  {f.file_name}")

# Pull down just the raw file (220 MB).
raw_files = pymzlib.pride.download("PXD000001", "downloads", category="RAW")

# Read it: mzML, Thermo .raw, Bruker .d, timsTOF .d, MGF and msalign all work the same way.
# Scan headers always; peaks only when you ask, because they are thousands of times larger.
scans = pymzlib.readers.read_spectra(raw_files[0], ms_order=2, limit=5, peaks=True)
print(scans.scan_count, scans.columns["selected_ion_mz"])
```

Downloads stream to a temporary name and are moved into place only when complete, so an
interrupted transfer never leaves you with a truncated file. Pass `overwrite=False` to skip
files you already have — a cheap resume for a large project.

## What's covered

Coverage is deliberately partial and grows by demand, the same way pyOpenMS grew.

| Area | Status |
|---|---|
| PRIDE Archive — list project files, filtered download, the MD5 and exact size to check each download against | ✅ |
| Peptidoforms — digest an annotated protein, apply its modifications, fragment every peptide; convert full sequences to Unimod or ProForma | ✅ |
| Quantification — FlashLFQ label-free quant with match-between-runs, and median-polish protein roll-up | ✅ |
| Readers — read spectra from **mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign; identify and read every file type mzLib knows, search results included; many files in one call; MetaMorpheus protein groups, FlashLFQ peptides and PTM site occupancy as long tables, RNA transcript groups and oligos included | ✅ |
| SDRF experimental design — read an SDRF-Proteomics file, pool several experiments into one analysis table; validate, lint and assess them; one row per sample with ages in years; the label-free design FlashLFQ and MetaMorpheus take, and the isobaric (TMT) channel design MetaMorpheus takes, or every reason either was refused; restrict an SDRF to one search's files | ✅ |
| Protein databases (`proteins`) — organism, taxon, GO terms and Ensembl genes per accession; resolve proteins to Ensembl genes against a pinned release; peptide uniqueness with I = L; Gene Ontology on protein groups, every member kept; align sequences and map residues; read accessions into entry, isoform and variants; Ensembl Compara orthologs across species | ✅ |
| Isobaric kits (`isobaric`) — every TMT, TMTpro, iTRAQ and DiLeu channel with its reporter-ion m/z and matching window, as MetaMorpheus quantifies them | ✅ |
| Spectral libraries (`library`) — read an mzLib `.msl` spectral library: every precursor with its m/z, iRT or minutes, ion mobility and target/decoy/entrapment label, and its fragment ions on request | ✅ |
| Differential abundance (`stats`) — limma's moderated t-test (checked against limma to 1e-8), Benjamini-Hochberg, random-effects meta-analysis and p-value combination (Fisher, Stouffer, signed) across studies, Spearman correlation with a trait (partial too), logistic regression of detection, mixed models with a random intercept per group, and entrapment checks of a search's FDR, with no R | ✅ |
| Differential abundance (`stats`) — limma's moderated t-test (checked against limma to 1e-8), Benjamini-Hochberg, random-effects meta-analysis and p-value combination (Fisher, Stouffer, signed) across studies, and Spearman correlation with a trait (partial too), with no R | ✅ |
| Differential design and observations (`diff`) — a sample sheet to mzLib's design matrix per stratum (equal to R's `model.matrix`), with each default contrast's estimability and the design file `stats.fit` reads; and every peptide in every sample of a FlashLFQ or MetaMorpheus table as log2 intensity or no value with the reason, with mzLib's shared-peptide normalization | ✅ |
| PTM occupancy and PTM-QTL (`ptmqtl`) — site occupancy per run or per sample, every cell with its censoring state; PTM pairs on one molecule or co-varying across runs, with a batch removed; the effect of a per-run trait on a site, replicates kept apart | ✅ |
| Everything else in mzLib | not yet — [tell us what you need](https://github.com/smith-chem-wisc/pyMzLib/issues) |

## How it works, and why you probably don't care

mzLib runs as a self-contained executable bundled inside this package; pyMzLib starts it,
speaks JSON to it, and hands you ordinary Python objects. The mechanism is an implementation
detail and may change — the public API will not.

The consequence worth knowing: calls carry a small fixed startup cost (tens of milliseconds),
so pyMzLib is built for coarse-grained operations, not for calling in a tight per-spectrum loop.

## Errors

Everything mzLib can fail at arrives as a typed Python exception:

```python
try:
    pymzlib.pride.list_files("PXD000001")
except pymzlib.UsageError:      # you passed something invalid
    ...
except pymzlib.BridgeError as e:  # mzLib/PRIDE failed; e.error_type says how
    ...
```

## Development

```powershell
# Build the .NET payload and stage it into the package
.\..\build\publish-bridge.ps1

# Set up a dev environment and run the fast tests
python -m venv .venv
.\.venv\Scripts\python -m pip install -e ".[dev]"
.\.venv\Scripts\python -m pytest -m "not network"
```

`pytest -m network` additionally hits the live PRIDE Archive.

## License

LGPL-3.0-or-later, matching mzLib.
