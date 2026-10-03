# Getting started

This page assumes nothing. If you already work in Python daily, you only need the first section.

## Install

```bash
pip install mzlib
```

Then check it:

```pycon
>>> import pymzlib
>>> v = pymzlib.bridge_version()
>>> v["protocol"]
1
>>> "pride files" in v["verbs"]
True

```

If that works, everything works, including the .NET runtime bundled inside the package. `v["runtime"]`
is that runtime's version and `v["mzlib"]` the mzLib commit the bridge was built from.

Every `>>>` example on this page runs in CI against output recorded from the real bridge. Blocks
titled **Not run** say why.

## Your first real call

```pycon
>>> files = pymzlib.pride.list_files("PXD000001")
>>> for f in files[:4]:
...     print(f"{f.category:8s} {f.size_mb:9.1f} MB  {f.file_name}")
OTHER          0.5 MB  PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz
PEAK          16.4 MB  PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz
PEAK         243.0 MB  TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML
OTHER          1.7 MB  erwinia_carotovora.fasta

```

Continue with the [PRIDE guide](guides/pride.md).

---

## If Python itself is new to you

Nothing below is required to use pyMzLib. It's here because the questions are real and the
answers are not obvious from the outside.

### "Which Python do I have?"

Probably several, and that is the normal state of affairs. Python has no equivalent of a global
assembly cache: each installation carries its own separate set of packages, so "is X installed?"
is only answerable relative to *which* Python you mean.

pyMzLib is built so that this matters as little as possible. It works on Python 3.9 and newer,
declares no dependencies, and therefore cannot be the reason an install fails. Whichever Python
you have, it will very likely just work.

### "Will installing this break something else?"

It cannot. Python allows exactly one version of any given package per environment, so packages
that pull in `numpy` or `pandas` can genuinely deadlock each other's requirements. pyMzLib
requires **nothing**, so it has nothing to disagree with. This was a design constraint, not an
accident.

### The one habit worth adopting: a virtual environment

A virtual environment is a private folder holding one Python and its own packages, so that
projects can't disturb each other. It ships with Python — nothing to install:

=== "Windows"

    ```powershell
    python -m venv .venv
    .\.venv\Scripts\Activate.ps1
    pip install mzlib
    ```

=== "macOS / Linux"

    ```bash
    python -m venv .venv
    source .venv/bin/activate
    pip install mzlib
    ```

Once activated, `python` and `pip` mean the ones inside `.venv`. Delete the folder and every
trace is gone. If a Python environment ever confuses you, deleting `.venv` and recreating it is a
legitimate and complete fix — nothing is registered anywhere else.

!!! tip "If you want the modern version of this"
    [`uv`](https://docs.astral.sh/uv/) is a single executable that manages Python interpreters
    *and* environments, so you never install Python at all: `uv venv && uv pip install mzlib`.
    It is dramatically faster and removes most remaining ways to get confused.

### Reading the results

`list_files()` returns a plain Python list of [`PrideFile`](reference.md) objects. Ordinary list
operations work on it:

```pycon
>>> raw = [f for f in files if f.category == "RAW"]            # filter
>>> big = [f for f in files if f.size_mb > 100]                # filter
>>> total = sum(f.file_size_bytes for f in files)              # aggregate
>>> by_size = sorted(files, key=lambda f: -f.file_size_bytes)  # sort
>>> len(raw), len(big), round(total / 1e9, 2), by_size[0].extension
(1, 2, 0.51, '.mzxml')

```

If you use `pandas`, a manifest converts in one line — pyMzLib doesn't require pandas, but it
doesn't get in the way either. Use `as_dict()` rather than `vars()`: `size_mb`, `extension` and
`downloadable` are computed properties, so `vars()` silently leaves them out.

```python title="Not run: needs pandas, which pyMzLib does not depend on"
import pandas as pd
df = pd.DataFrame([f.as_dict() for f in files])
```

### When something goes wrong

pyMzLib raises typed exceptions, so you can tell a mistake in your call from a failure out in the
world:

```pycon
>>> try:
...     files = pymzlib.pride.list_files("PXD000001")
... except pymzlib.UsageError:
...     pass   # something about the call itself was invalid
... except pymzlib.BridgeError as e:
...     pass   # PRIDE or mzLib failed; e.error_type says which way

```

A mistake in the call is caught before anything runs:

```pycon
>>> pymzlib.pride.list_files("")                       # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
pymzlib._bridge.UsageError: ...

```

See the [FAQ](faq.md) for specific failure messages.
