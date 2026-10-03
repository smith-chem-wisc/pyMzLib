# pyMzLib

**mzLib for Python.**
[mzLib](https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry and proteomics library
written in C#, developed in the [Smith lab](https://smith.chem.wisc.edu/) at UW–Madison. pyMzLib
makes it callable from Python.

```bash
pip install mzlib
```

Install it as `mzlib`, import it as `pymzlib` — the same split as `pip install scikit-learn` /
`import sklearn`. The project is pyMzLib.

That is the entire installation. No .NET to install, no runtime to configure, no version to
reconcile — and **no third-party Python dependencies**, so pyMzLib cannot conflict with anything
already in your environment.

```pycon
>>> import pymzlib
>>> files = pymzlib.pride.list_files("PXD000001")
>>> print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")
8 files, 0.51 GB

```

```python title="Not run: downloads from PRIDE"
pymzlib.pride.download("PXD000001", "downloads", category="RAW")
```

---

## Why this exists

A great deal of computational proteomics happens in Python. mzLib holds a decade of carefully
tested proteomics machinery — chemistry, spectra, digestion, deconvolution, repository access —
and until now none of it was reachable from a Python prompt. The result is that Python users
reimplement things mzLib already does correctly.

The hard part was never the calling convention. It was making the result **frictionless**: a
wrapper that asks you to install a .NET runtime, match a Python version, or resolve a dependency
conflict has failed, however complete its API. So that constraint drove the design, and the
measure of success is the two commands above.

## What's covered

Coverage is deliberately partial and grows by demand — the same way
[pyOpenMS](https://pyopenms.readthedocs.io/) grew.

| Area | Status |
|---|---|
| [PRIDE Archive](guides/pride.md) — search for projects by keyword, list a project's files, filtered download | :material-check: available |
| [Peptidoforms](guides/peptidoforms.md) — digest an annotated protein, apply its modifications, fragment every peptide; convert full sequences to Unimod or ProForma | :material-check: available |
| [FlashLFQ](guides/flashlfq.md) — label-free quantification across mzML runs, with match-between-runs | :material-check: available |
| [Readers](guides/readers.md) — read spectra from **mzML**, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF and msalign; identify any format mzLib knows and read it, search results included; many files in one call; MetaMorpheus protein groups, FlashLFQ peptides and PTM site occupancy as long tables, RNA transcript groups and oligos included | :material-check: available |
| [SDRF experimental design](guides/sdrf.md) — read an SDRF-Proteomics file, pool several into one analysis table; validate, lint and assess them; one row per sample with ages in years; the label-free design FlashLFQ and MetaMorpheus take, or every reason it was refused | :material-check: available |
| [Protein databases](guides/proteins.md) — organism, taxon, GO terms and Ensembl genes per accession from a UniProt XML or FASTA; resolve proteins to Ensembl genes against a pinned release; peptide uniqueness with I = L; Gene Ontology on MetaMorpheus protein groups, every member kept, against a pinned GO release | :material-check: available |
| [Isobaric kits](guides/isobaric.md) — every TMT, TMTpro, iTRAQ and DiLeu channel with its reporter-ion m/z and matching window | :material-check: available |
| [Differential abundance](guides/stats.md) — limma's moderated t-test (checked against limma to 1e-8), Benjamini-Hochberg, and random-effects meta-analysis across studies, with no R | :material-check: available |
| Everything else in mzLib | not yet — [tell us what you need](https://github.com/smith-chem-wisc/pyMzLib/issues) |

If there is something in mzLib you want from Python, opening an issue is genuinely the fastest
path. The [extension recipe](contributing/adding-a-capability.md) is short, and requests are how
we decide what to cover next.

## Where to go next

<div class="grid cards" markdown>

- :material-rocket-launch: **[Getting started](getting-started.md)**
  Install it and run something, written for people who don't live in Python.

- :material-database-search: **[PRIDE guide](guides/pride.md)**
  The first covered area, end to end.

- :material-book-open-variant: **[API reference](reference.md)**
  Generated from the source, so it cannot drift.

- :material-cog: **[How it works](design/architecture.md)**
  What's actually happening when you call it, and why it was built this way.

</div>

## Citing

pyMzLib has no paper of its own yet, and no release DOI has been minted. Releases are on
[PyPI](https://pypi.org/project/mzlib/#history) and
[GitHub](https://github.com/smith-chem-wisc/pyMzLib/releases), and
[`CITATION.cff`](https://github.com/smith-chem-wisc/pyMzLib/blob/main/CITATION.cff) gives the
software citation; GitHub's "Cite this repository" button renders it.

pyMzLib is an interface to mzLib, not a reimplementation of it, so cite the science where it
lives:

- **mzLib**, the library every result comes from: cite it as software,
  [smith-chem-wisc/mzLib](https://github.com/smith-chem-wisc/mzLib), with the version
  `pymzlib.bridge_version()["mzlib"]` reports.
- **MetaMorpheus**, whose search, digestion and modification handling mzLib carries: Solntsev S.K.,
  Shortreed M.R., Frey B.L., Smith L.M. Enhanced Global Post-translational Modification Discovery
  with MetaMorpheus. *J. Proteome Res.* **17**, 1844–1851 (2018).
  [doi:10.1021/acs.jproteome.7b00873](https://doi.org/10.1021/acs.jproteome.7b00873)
- **The method or resource behind each function**: every guide ends with a *Cite* section listing
  the papers its functions rest on (FlashLFQ, PRIDE, UniProt, SDRF-Proteomics, the Gene Ontology,
  Ensembl, TMTpro), taken from the same specs the reference pages are built from.
