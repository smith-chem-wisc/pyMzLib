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
        element_count = chemicalFormula.AtomCount,
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

Check it directly, before any Python exists:

```bash
mzlib-bridge chemistry formula-mass --formula H2O
{"ok":true,"data":{"formula":"H2O","monoisotopic_mass":18.0105646863,…}}
```

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
    """

    formula: str
    monoisotopic_mass: float
    average_mass: float
    element_count: int


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
        18.0105646863
    """
    if not formula or not formula.strip():
        raise _bridge.UsageError("A chemical formula is required, e.g. 'H2O'.")

    data = _bridge.invoke("chemistry", "formula-mass", "--formula", formula.strip(), timeout=timeout)
    return FormulaMass(
        formula=data["formula"],
        monoisotopic_mass=float(data["monoisotopic_mass"]),
        average_mass=float(data["average_mass"]),
        element_count=int(data["element_count"]),
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

Export it in `__init__.py`:

```python
from . import chemistry, pride
```

## 3 · Test it, both ways

**Offline** (`pytest -m "not network"`) — the Python layer against a recorded payload. These are
the tests that run constantly, so make them cover the behavior you actually care about:

```python
def test_formula_mass_parses(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "formula": "H2O", "monoisotopic_mass": 18.0105646863,
        "average_mass": 18.01528, "element_count": 3,
    })
    assert chemistry.formula_mass("H2O").monoisotopic_mass == pytest.approx(18.01056, abs=1e-5)


def test_blank_formula_rejected_without_starting_a_process():
    with pytest.raises(pymzlib.UsageError):
        chemistry.formula_mass("   ")
```

**Live** (`@pytest.mark.network`, or `@pytest.mark.slow` if it moves real data) — the real bridge,
which is what would catch mzLib changing underneath.

Record a fixture from the real bridge rather than hand-writing one; hand-written fixtures encode
what you *think* the format is.

## 4 · Document it

- A guide page under `docs/guides/` if it's a new area, following
  [the PRIDE guide](../guides/pride.md): what it does, the common cases, the errors, one worked
  example that solves a real problem.
- Add it to the nav in `mkdocs.yml` and to the coverage table on the [home page](../index.md).
- The API reference picks up your docstrings automatically — no work needed, provided they're
  written.

## Checklist

- [ ] Coarse-grained enough that per-call startup is irrelevant
- [ ] No logic reimplemented that mzLib already has
- [ ] Bridge verb returns an anonymous object with `snake_case` keys
- [ ] Errors propagate; bad input throws `UsageException`
- [ ] Verb readable as a contract by a non-Python caller
- [ ] Python returns a dataclass, validates cheaply, hides the transport
- [ ] **No new runtime dependency**
- [ ] Offline tests with a recorded fixture; live test for the real path
- [ ] Guide page, nav entry, coverage table
