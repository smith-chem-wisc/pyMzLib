"""Shared reading of the vendored verb specs in docs/specs/.

Used by scripts/render_spec_docs.py (the reference fact tables), pkg/python/tests/test_spec_docs.py
(the docstring lint) and pkg/python/src/conftest.py (the doctest replay bridge), so the three agree
on what a spec says and how its wire names become Python names.

Needs PyYAML, which the docs toolchain already brings (MkDocs depends on it). pyMzLib itself still
has no runtime dependency: nothing under pkg/python/src/pymzlib imports this.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
SPECS = ROOT / "docs" / "specs"
GENERATED = ROOT / "docs" / "reference" / "_generated"
FIXTURES = ROOT / "pkg" / "python" / "tests" / "fixtures"

#: How the spec's language-neutral error kinds surface in Python. The C# bridge reports ``usage``
#: and ``ServiceUnavailable`` by name; anything else is a correctness failure carrying mzLib's own
#: exception type in ``BridgeError.error_type``. See pymzlib._bridge.invoke.
PY_EXCEPTION = {
    "usage": "UsageError",
    "service_unavailable": "ServiceUnavailableError",
    "correctness": "BridgeError",
}

#: Every place pyMzLib deliberately projects a spec param or result field under a different name,
#: or not as a docstring-visible attribute at all. **Nothing else may be skipped by the lint**, and
#: the lint fails on an entry here that names no param or field of the spec, so this list cannot
#: go stale silently.
#:
#: Key: the spec's ``verb``. Value: ``{"param.<wire name>" | "field.<wire name>": {...}}`` with
#: ``python`` - the Python name, or None when it is not a named attribute - and ``why``.
PYTHON_DEVIATIONS: dict[str, dict[str, dict[str, str | None]]] = {
    "readers formats": {
        # formats() returns list[Format] rather than an envelope object: a list is what a caller
        # iterates, and its length already is the count.
        "field.format_count": {
            "python": None,
            "why": "formats() returns a list; this is len(formats())",
        },
        "field.formats": {
            "python": None,
            "why": "formats() returns this list itself, as Format objects",
        },
    },
    # SDRF cannot use the readers' ``columns`` name-to-values map (its names repeat), so a document
    # carries its header as ``columns`` - a list - and ``column_names`` is that list.
    "sdrf read": {
        "field.column_names": {"python": "columns", "why": "the header list is SdrfDocument.columns"},
    },
    "sdrf pool": {
        "param.stdin": {"python": "documents", "why": "the stdin lines are rendered from documents"},
        "field.column_names": {"python": "columns", "why": "the header list is PooledSdrf.columns"},
    },
    "sdrf lint": {
        "param.stdin": {"python": "documents", "why": "the stdin lines are rendered from documents"},
    },
    "sdrf parse-age": {
        "param.stdin": {"python": "cells", "why": "one stdin line per element of cells"},
    },
    # One document and many are two functions (validate / validate_many), a cross-binding decision:
    # the bulk options exist only on the _many form, which takes the list the wire reads on stdin.
    **{
        f"sdrf {verb}": {
            "param.paths-stdin": {"python": None, "why": f"{verb}_many(paths) is the bulk form"},
            "param.threads": {"python": None, "why": f"only {verb}_many() takes threads"},
            "param.on-error": {"python": None, "why": f"only {verb}_many() takes on_error"},
        }
        for verb in ("validate", "assess", "samples")
    },
}


def load_specs() -> list[dict[str, Any]]:
    """Every vendored spec, sorted by file name, each with ``_file`` set to its file name."""
    specs = []
    for path in sorted(SPECS.glob("*.yaml")):
        spec = yaml.safe_load(path.read_text(encoding="utf-8"))
        spec["_file"] = path.name
        specs.append(spec)
    return specs


def source_commit() -> str:
    """The bridge commit the vendored specs came from, per docs/specs/SOURCE."""
    source = SPECS / "SOURCE"
    if not source.exists():
        return "unknown"
    for line in source.read_text(encoding="utf-8").splitlines():
        if line.startswith("commit:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def py_param_name(wire: str) -> str:
    """The Python keyword for a wire option: ``ms-order`` -> ``ms_order``."""
    return wire.replace("-", "_")


def fragment_name(spec: dict[str, Any]) -> str:
    """``readers.read-spectra`` for the spec ``readers.read-spectra.yaml``."""
    return spec["_file"][: -len(".yaml")]


def py_binding(spec: dict[str, Any]) -> str | None:
    return ((spec.get("bindings") or {}).get("py") or {}).get("name")


def shipped_in_python(spec: dict[str, Any]) -> bool:
    """Whether the spec says pyMzLib has shipped this verb (``since.pymzlib`` is set)."""
    return bool((spec.get("since") or {}).get("pymzlib"))


def columns(spec: dict[str, Any]) -> list[dict[str, Any]]:
    cols = (spec.get("result") or {}).get("columns")
    return cols if isinstance(cols, list) else []


def find_fixture(name: str) -> Path | None:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    hits = sorted((ROOT / "pkg").rglob(name))
    return hits[0] if hits else None
