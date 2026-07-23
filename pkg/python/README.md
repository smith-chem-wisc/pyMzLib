# pyMzLib

**mzLib for Python.** [mzLib](https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry
and proteomics library written in C#, developed in the Smith lab at UW–Madison. pyMzLib makes it
callable from Python.

```bash
pip install pymzlib
```

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

# Pull down just the raw files.
paths = pymzlib.pride.download("PXD000001", "downloads", category="RAW")
```

Downloads stream to a temporary name and are moved into place only when complete, so an
interrupted transfer never leaves you with a truncated file. Pass `overwrite=False` to skip
files you already have — a cheap resume for a large project.

## What's covered

Coverage is deliberately partial and grows by demand, the same way pyOpenMS grew.

| Area | Status |
|---|---|
| PRIDE Archive — list project files, filtered download | ✅ |
| Everything else in mzLib | not yet — [tell us what you need](https://github.com/smith-chem-wisc/mzLib/issues) |

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
