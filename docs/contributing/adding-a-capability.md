# Adding a capability

How to expose another piece of mzLib. The recipe is short by design — if adding a capability
were hard, coverage would never grow.

Worked end to end below with a hypothetical `pymzlib.chemistry.formula_mass()`.

## Before you write anything

**Check that it belongs here.** Two questions:

1. *Is it coarse-grained?* Each call costs tens of milliseconds of process startup and must
   serialize its data as JSON. Fetching a project manifest: perfect. Computing an isotopic
   envelope inside a loop over 100,000 spectra: wrong — that needs a different transport, and
   forcing it through this one produces something slow enough to discredit the whole package.
2. *Does mzLib already do it?* Don't reimplement logic on the bridge side. The bridge is a
   translation layer. If mzLib's version needs a change, change mzLib.

**Check the dependency cost.** Referencing a new mzLib project may drag substantial weight into
the payload. Not a blocker ([D8](../design/decisions.md#d8-payload-size-is-not-a-design-constraint)),
but worth knowing before rather than after.

## 1 · Add the verb to the bridge

In `pkg/bridge/Program.cs`, route it in `DispatchAsync` and write the handler:

```csharp
"chemistry formula-mass" => ChemistryFormulaMass(arguments),
```

The handler below needs `using Chemistry;` at the top of `Program.cs`. Before reaching for a new
`ProjectReference`, check whether the type is already reachable — `Chemistry` arrives transitively
through `UsefulProteomicsDatabases`, so this particular example needs no new reference at all.

```csharp
/// <summary>
/// <c>chemistry formula-mass --formula H2O</c> — the monoisotopic and average mass of a
/// chemical formula.
/// </summary>
private static object ChemistryFormulaMass(Arguments arguments)
{
    string formula = arguments.Required("formula");
    var chemicalFormula = ChemicalFormula.ParseFormula(formula);

    return new
    {
        formula,
        monoisotopic_mass = chemicalFormula.MonoisotopicMass,
        average_mass = chemicalFormula.AverageMass,
        atom_count = chemicalFormula.AtomCount,
    };
}
```

Rules that matter:

- **Return an anonymous object, not an mzLib type.** mzLib's shapes carry things a JSON consumer
  can't use and would couple the wire format to mzLib's internals.
- **`snake_case` keys.** `JsonNamingPolicy.SnakeCaseLower` handles this; don't fight it.
- **Never catch and swallow.** `Main` converts any exception into a structured error envelope.
  Letting it propagate gives the caller the .NET exception type, which is the useful part.
- **Throw `UsageException` for bad input**, so it exits `2` and surfaces as `UsageError` in Python
  rather than looking like a runtime failure.
- **Name it for the contract, not the caller.** Someone writing a Rust binding will read this —
  see [D6](../design/decisions.md#d6-the-wire-contract-stays-language-neutral).

Check it directly, before any Python exists. The executable is not on `PATH` — run the one
`publish-bridge.ps1` staged into the package:

```bash
./pkg/python/src/pymzlib/_dotnet/win-x64/mzlib-bridge.exe chemistry formula-mass --formula H2O
{"ok":true,"data":{"formula":"H2O","monoisotopic_mass":18.01056468403,"average_mass":18.015349999999998,"atom_count":3},"error":null}
```

Note `"error": null` — the envelope always carries all three keys, so a consumer in a statically
typed language can deserialize one fixed shape whatever the outcome.

!!! warning "Re-stage the bridge, or Python keeps using the old one"
    `bridge_path()` prefers the executable staged under `src/pymzlib/_dotnet/`, so after changing
    C# you must either re-run `publish-bridge.ps1` or point `PYMZLIB_BRIDGE` at your
    `dotnet build` output. Otherwise your new verb reports `Unknown command 'chemistry
    formula-mass'` while you stare at the registration you just wrote. The build page frames
    `PYMZLIB_BRIDGE` as a speed optimisation; here it is load-bearing for correctness.

## 2 · Add the Python surface

A new area gets its own module (`pkg/python/src/pymzlib/chemistry.py`); an addition to an
existing area goes in that module.

```python
"""Chemical formula and mass calculations, backed by mzLib's Chemistry library."""

from __future__ import annotations          # required — keeps modern syntax legal on 3.9

from dataclasses import dataclass

from . import _bridge

__all__ = ["FormulaMass", "formula_mass"]


@dataclass(frozen=True)
class FormulaMass:
    """The masses of a chemical formula.

    Attributes:
        formula: The formula as given, e.g. ``"H2O"``.
        monoisotopic_mass: Mass using the most abundant isotope of each element, in daltons.
        average_mass: Mass weighted by natural isotopic abundance, in daltons.
        atom_count: The number of atoms — 3 for water, not 2. Name a wire key for what mzLib
            actually returns: the envelope is a durable, language-neutral contract, so a
            misnamed key outlives the mistake and every copy of this example repeats it.
    """

    formula: str
    monoisotopic_mass: float
    average_mass: float
    atom_count: int


def formula_mass(formula: str, timeout: float | None = 60) -> FormulaMass:
    """Return the monoisotopic and average mass of a chemical formula.

    Args:
        formula: A chemical formula in mzLib's notation, e.g. ``"H2O"``, ``"C6H12O6"``.
        timeout: Seconds to allow.

    Returns:
        The computed masses.

    Raises:
        UsageError: the formula is blank or cannot be parsed.

    Example:
        >>> formula_mass("H2O").monoisotopic_mass      # doctest: +SKIP
        18.01056468403
    """
    if not formula or not formula.strip():
        raise _bridge.UsageError("A chemical formula is required, e.g. 'H2O'.")

    data = _bridge.invoke("chemistry", "formula-mass", "--formula", formula.strip(), timeout=timeout)
    return FormulaMass(
        formula=data["formula"],
        monoisotopic_mass=float(data["monoisotopic_mass"]),
        average_mass=float(data["average_mass"]),
        atom_count=int(data["atom_count"]),
    )
```

Conventions:

- **Validate in Python too, where it's cheap.** Catching a blank argument without starting a
  process gives a better error and a faster one. Not duplicated logic — a fast path in front of
  the authoritative check.
- **Return dataclasses, not dicts.** Attributes autocomplete, typos raise, and the API reference
  documents itself.
- **Google-style docstrings** — the [API reference](../reference.md) is generated from them.
- **Never leak the transport.** No `subprocess`, no exit codes, no raw JSON above `_bridge.py`.
- **No new runtime dependency.** Not negotiable — see
  [D2](../design/decisions.md#d2-it-must-install-with-one-command-and-nothing-else).

Export it in `__init__.py` — the import **and** `__all__`, or `from pymzlib import *` misses it:

```python
from . import chemistry, pride

__all__ = ["pride", "chemistry", ...]
```

## 3 · Test it, both ways

**Offline** (`pytest -m "not network"`) — the Python layer against a recorded payload. These are
the tests that run constantly, so make them cover the behavior you actually care about:

```python
from __future__ import annotations

import pytest

import pymzlib
from pymzlib import _bridge, chemistry


def test_formula_mass_parses(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "formula": "H2O", "monoisotopic_mass": 18.01056468403,
        "average_mass": 18.015349999999998, "atom_count": 3,
    })
    assert chemistry.formula_mass("H2O").monoisotopic_mass == pytest.approx(18.01056, abs=1e-5)


def test_blank_formula_rejected_without_starting_a_process():
    with pytest.raises(pymzlib.UsageError):
        chemistry.formula_mass("   ")
```

**Live** (`@pytest.mark.network`, or `@pytest.mark.slow` if it moves real data) — the real bridge,
which is what would catch mzLib changing underneath. Put these in a `_live.py` module and route
each through `external_service()` so a PRIDE outage skips rather than fails.

Record a fixture from the real bridge rather than hand-writing one; hand-written fixtures encode
what you *think* the format is. (The stub payload above is inline only to keep the example short.)

## 3b · Test the C# half too — CI gates it separately

**This is the step most likely to be skipped, and skipping it fails the build.** The bridge has its
own suite and its own coverage gate at 85%. New C# lines with no C# test will drag the whole
project under it, and nothing in the Python suite will warn you.

Copy the shape from `pkg/bridge.tests/VerbHandlerTests.cs` — a `StubHandler` for anything that
would otherwise reach the network, plus direct tests of whatever the verb computes:

```bash
dotnet test pkg/bridge.tests/MzLibBridge.Tests.csproj --filter "TestCategory!=ExternalService"
./pkg/build/check-bridge-coverage.ps1
```

Both must pass before you open the pull request. The raw coverage number `dotnet test` prints is
meaningless here — it measures all of mzLib and reads about 0.6% — which is exactly why the script
exists. See the testing section of [CONTRIBUTING.md](https://github.com/smith-chem-wisc/pyMzLib/blob/main/CONTRIBUTING.md)
for the full picture, including the skip-versus-fail convention for anything touching a live service.

## 4 · Document it

- A guide page under `docs/guides/` if it's a new area, following
  [the PRIDE guide](../guides/pride.md): what it does, the common cases, the errors, one worked
  example that solves a real problem.
- Add it to the nav in `mkdocs.yml` and to the coverage table on the [home page](../index.md).
- Add `::: pymzlib.<module>` to `docs/reference.md`. The reference is generated from your
  docstrings, but only for modules it is told about — a new module is silently absent otherwise.

## Checklist

- [ ] Coarse-grained enough that per-call startup is irrelevant
- [ ] No logic reimplemented that mzLib already has
- [ ] Bridge verb returns an anonymous object with `snake_case` keys
- [ ] Errors propagate; bad input throws `UsageException`
- [ ] Verb readable as a contract by a non-Python caller
- [ ] Python returns a dataclass, validates cheaply, hides the transport
- [ ] **No new runtime dependency**
- [ ] Offline Python tests with a recorded fixture; live test for the real path
- [ ] **C# tests for the new verb, and `check-bridge-coverage.ps1` still passes** (separate gate)
- [ ] Guide page, nav entry, coverage table, and `reference.md` entry
