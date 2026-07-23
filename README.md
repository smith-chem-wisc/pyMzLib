# pyMzLib

**mzLib for Python** — mass spectrometry and proteomics from the
[mzLib](https://github.com/smith-chem-wisc/mzLib) C# library, with no .NET to install.

```bash
pip install pymzlib
```

```python
import pymzlib

files = pymzlib.pride.list_files("PXD000001")
print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")

pymzlib.pride.download("PXD000001", "downloads", category="RAW")

# Digest an annotated UniProt protein and fragment its peptides
digest = pymzlib.peptidoform.fragments("P02768")
print(digest.modification_census.explain())
```

> **Preview release available — not on PyPI yet.** Install the wheel for your OS from the
> [latest release](https://github.com/smith-chem-wisc/pyMzLib/releases/latest). For example, on
> Windows:
>
> ```bash
> pip install https://github.com/smith-chem-wisc/pyMzLib/releases/download/v0.1.0.dev0/pymzlib-0.1.0.dev0-py3-none-win_amd64.whl
> ```
>
> The release page lists the Linux and macOS (Intel / Apple Silicon) wheels too. `pip install
> pymzlib` is what release day will look like.

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
| [PRIDE Archive](https://smith-chem-wisc.github.io/pyMzLib/guides/pride/) — list a project's files, filtered download | ✅ |
| [Peptidoforms](https://smith-chem-wisc.github.io/pyMzLib/guides/peptidoforms/) — digest an annotated protein, apply its modifications, fragment every peptide | ✅ |
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
