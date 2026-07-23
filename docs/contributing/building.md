# Building from source

You need the [.NET SDK 8.0 or newer](https://dotnet.microsoft.com/download) and any Python 3.9+.
Nothing else.

## Layout

```
pyMzLib/
  pkg/bridge/          the .NET console app over mzLib
  pkg/python/          the pymzlib package, its tests and fixtures
  pkg/build/           publish-bridge.ps1 — stages the .NET payload into the package
  code/mzLib/          pinned mzLib checkout (not committed; see code/PINNED.md)
  docs/                this site
```

## One-time setup

mzLib itself is not vendored. Check out the pinned commit — the exact SHA is in
[`code/PINNED.md`](https://github.com/smith-chem-wisc/pyMzLib/blob/main/code/PINNED.md), and it
matters: the bridge and the Python package must be built from the same mzLib source.

=== "If you already have an mzLib clone"

    ```powershell
    git -C <your-mzLib-clone> worktree add --detach code/mzLib <PINNED_SHA>
    ```

    A worktree is cheaper than a second clone — it shares the object store.

=== "Fresh"

    ```bash
    git clone --filter=blob:none --no-checkout https://github.com/smith-chem-wisc/mzLib.git code/mzLib
    git -C code/mzLib checkout --detach <PINNED_SHA>
    ```

## Build the bridge

```powershell
./pkg/build/publish-bridge.ps1
```

PowerShell, but cross-platform — `pwsh` runs on Linux and macOS, and that's how CI invokes it.
This publishes a self-contained single-file executable and stages it into
`pkg/python/src/pymzlib/_dotnet/<rid>/`, where the Python package expects it. It ends by running
the executable's `version` verb, because a payload that can't report its own version here will
certainly fail on a user's machine.

Cross-build for another platform — no hardware required:

```powershell
./pkg/build/publish-bridge.ps1 -Runtime linux-x64
./pkg/build/publish-bridge.ps1 -Runtime osx-arm64
```

!!! tip "Iterating on the C# side"
    A full self-contained publish takes a minute. While developing, `dotnet build` the bridge and
    point `PYMZLIB_BRIDGE` at the resulting executable — the Python package will prefer it:

    ```powershell
    $env:PYMZLIB_BRIDGE = "pkg/bridge/bin/x64/Release/net8.0/win-x64/mzlib-bridge.exe"
    ```

## Set up Python and run the tests

```bash
cd pkg/python
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev]"    # Windows
# .venv/bin/python -m pip install -e ".[dev]"      # macOS / Linux
```

```bash
.venv/Scripts/python -m pytest -m "not network"    # fast: ~0.1 s, no network
.venv/Scripts/python -m pytest -m network          # live, hits EBI
.venv/Scripts/python -m pytest                     # everything
```

The offline tests run against a recorded PRIDE manifest in `tests/fixtures/`, so they work on a
plane and can't be broken by an EBI outage. The network tests are the ones that would catch PRIDE
or mzLib changing under us.

!!! warning "Write fixtures without a byte-order mark"
    On Windows PowerShell, `Out-File -Encoding utf8` writes a BOM, and Python's `json` module
    rejects it. Use `[System.IO.File]::WriteAllText($path, $json, (New-Object System.Text.UTF8Encoding($false)))`.

## Build a wheel

```bash
cd pkg/python
python -m pip install build wheel
python -m build --wheel --outdir dist
python -m wheel tags --python-tag py3 --abi-tag none --platform-tag win_amd64 --remove dist/*.whl
```

The retag is required. hatchling sees no compiled Python code and reasonably infers
`py3-none-any`, but the wheel carries a platform binary and must refuse to install on the wrong
operating system. This is also why there's one wheel per OS rather than one per OS × Python
version.

Verify it like a stranger would — install the built artifact into a clean environment, not the
source tree:

```bash
python -m venv /tmp/check
/tmp/check/bin/python -m pip install dist/*.whl
cd /tmp && /tmp/check/bin/python -c "import pymzlib; print(pymzlib.bridge_version())"
```

## Build the docs

```bash
pip install mkdocs-material mkdocstrings[python]
mkdocs serve      # http://127.0.0.1:8000, live-reloading
```

## Code style

**C#** — match mzLib's conventions: XML doc comments on public members documenting the failure
contract, explicit argument validation, `ConfigureAwait(false)` on library awaits.

**Python** — `snake_case`, Google-style docstrings (the API reference is generated from them),
`ruff` for linting, and `from __future__ import annotations` in every module so modern type
syntax stays legal on Python 3.9.

**The rule that isn't negotiable:** no third-party runtime dependency in `pymzlib`. Optional
extras are fine; base requirements are not. See
[D2](../design/decisions.md#d2-it-must-install-with-one-command-and-nothing-else).
