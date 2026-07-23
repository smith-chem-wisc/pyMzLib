# pyMzLib

**mzLib for Python.**
[mzLib](https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry and proteomics library
written in C#, developed in the [Smith lab](https://smith.chem.wisc.edu/) at UW–Madison. pyMzLib
makes it callable from Python.

```bash
pip install pymzlib
```

!!! warning "Not on PyPI yet"
    pyMzLib has not had its first release. Until it does, install a wheel from the
    [Actions artifacts](https://github.com/smith-chem-wisc/pyMzLib/actions/workflows/wheels.yml)
    or [build from source](contributing/building.md). The command above is what release day looks
    like, not what today looks like.

That is the entire installation. No .NET to install, no runtime to configure, no version to
reconcile — and **no third-party Python dependencies**, so pyMzLib cannot conflict with anything
already in your environment.

```python
import pymzlib

files = pymzlib.pride.list_files("PXD000001")
print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")

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
| [PRIDE Archive](guides/pride.md) — list a project's files, filtered download | :material-check: available |
| [Peptidoforms](guides/peptidoforms.md) — digest an annotated protein, apply its modifications, fragment every peptide | :material-check: available |
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

pyMzLib does not yet have a paper. Cite mzLib and MetaMorpheus for the underlying science; a
release DOI will appear here once the first version ships.
