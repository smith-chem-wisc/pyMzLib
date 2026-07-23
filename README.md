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
```

> **Not on PyPI yet.** pyMzLib has not had its first release. Until then, grab a wheel from the
> [Actions artifacts](https://github.com/smith-chem-wisc/pyMzLib/actions/workflows/wheels.yml) or
> build from source. The command above is what release day looks like.

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
