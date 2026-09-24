"""Lint the Python docstrings against the per-verb specs in docs/specs/.

The specs (vendored from the bridge by scripts/sync_specs.py) own the facts about each wire verb:
its parameters, their units, its result fields and theirs. The docstrings own the Python idiom and
the prose. This test holds the prose to the facts:

* every spec param is a keyword of the Python function and is described under ``Args:``;
* every spec result field is described under ``Attributes:`` of the class the function returns
  (or, when the function returns a list of row objects, every per-row field is described on the
  row class);
* a param or field whose spec ``unit`` is not null mentions that unit in its description, so
  "Skip this many" cannot ship without saying "scans".

A binding may deliberately project a field under another name, or not as an attribute at all
(``readers.formats()`` returns a list, so ``format_count`` is ``len()``). Those are declared, with a
reason, in ``PYTHON_DEVIATIONS`` in scripts/spec_facts.py, and nowhere else. A declaration that no
longer names a spec param or field fails here too, so the list cannot rot.

It also checks that the rendered fact tables are current and that every shipped verb's table is
included in a reference page. See docs/contributing/reference-facts.md.

Needs PyYAML and griffe (the docs toolchain: ``pip install -e ".[dev]"`` brings both). The wheel
test matrix installs only pytest, so there the module skips; docs.yml runs it for real.
"""

from __future__ import annotations

import importlib
import inspect
import logging
import re
import sys
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
griffe = pytest.importorskip("griffe")

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))

import render_spec_docs  # noqa: E402
import spec_facts  # noqa: E402

logging.getLogger("griffe").setLevel(logging.ERROR)

SPECS = spec_facts.load_specs()

#: The keys a param or a result field may carry (bridge design/verbs/README.md). Anything else is
#: nearly always YAML's flow-mapping trap: an unquoted comma inside ``{... doc: a, b}`` ends the
#: doc at the comma and turns the rest into a stray key, silently truncating the text every binding
#: renders.
FIELD_KEYS = {
    "name",
    "wire",
    "type",
    "required",
    "unit",
    "range",
    "default",
    "doc",
    "nullable",
    "null_means",
    "present_when",
}

#: Spec defects already reported to the bridge and not yet fixed there: (spec file, param or field
#: name). Vendored specs are never hand-edited, so the fix is a re-sync. Delete an entry when the
#: bridge fixes it; the test fails on an entry that no longer matches a defect.
KNOWN_SPEC_DEFECTS = {
    # Unquoted commas truncate these docs; reported to bridge 2026-09-23.
    ("readers.read-spectra.yaml", "ms_order"),  # envelope "The ms-order filter, echoed."
    ("readers.read-spectra.yaml", "intensity"),  # "Peak intensities, parallel to mz."
    ("readers.read-records.yaml", "columns"),  # "..., in column_names order."
}

UNIT_ALIASES = {
    "min": ["minute"],
    "s": ["second"],
    "ms": ["millisecond"],
    "da": ["dalton"],
    "v": ["volt"],
}


def unit_mentioned(unit: str, text: str) -> bool:
    """Whether ``text`` names ``unit``: the unit itself, its singular, or a spelled-out alias."""
    text = text.lower()
    unit = unit.lower().strip()
    head = unit.split(" (")[0].strip()
    candidates = {unit, head}
    if head.endswith("s") and len(head) > 3:
        candidates.add(head[:-1])
    candidates.update(UNIT_ALIASES.get(head, []))
    for c in candidates:
        if len(c) <= 2 or not c.isalpha():
            if re.search(r"(?<![a-z])" + re.escape(c) + r"(?![a-z])", text):
                return True
        elif re.search(r"(?<![a-z])" + re.escape(c), text):
            return True
    return False


def sections(obj) -> dict:
    doc = inspect.getdoc(obj) or ""
    parsed = griffe.Docstring(doc, parser="google").parse()
    out: dict = {}
    for section in parsed:
        if isinstance(section.value, list):
            out.setdefault(section.kind.value, {}).update(
                {
                    item.name: item.description or ""
                    for item in section.value
                    if hasattr(item, "name")
                }
            )
    return out


def resolve(dotted: str):
    module_name, _, attr = dotted.rpartition(".")
    module = importlib.import_module(module_name)
    return module, getattr(module, attr)


def returned(module, func):
    """(class, is_list) for the function's return annotation, resolved in its module."""
    annotation = inspect.signature(func).return_annotation
    if not isinstance(annotation, str):
        annotation = getattr(annotation, "__name__", str(annotation))
    match = re.fullmatch(r"\s*(?:list|List)\[\s*(\w+)\s*\]\s*", annotation)
    if match:
        return getattr(module, match.group(1), None), True
    return getattr(module, annotation.strip(), None), False


def deviation(spec, key):
    return spec_facts.PYTHON_DEVIATIONS.get(spec["verb"], {}).get(key)


def spec_id(spec):
    return spec["_file"]


# --------------------------------------------------------------------------------------------
# The vendored specs themselves
# --------------------------------------------------------------------------------------------


def test_specs_are_vendored_with_their_source():
    assert SPECS, "docs/specs/ holds no specs; run scripts/sync_specs.py"
    assert (spec_facts.SPECS / "SOURCE").exists(), (
        "docs/specs/SOURCE is missing; re-run scripts/sync_specs.py"
    )


def _items(spec):
    result = spec.get("result") or {}
    yield from spec.get("params") or []
    yield from result.get("envelope_fields") or []
    yield from spec_facts.columns(spec)


def test_no_spec_field_was_truncated_by_yaml():
    found = set()
    for spec in SPECS:
        for item in _items(spec):
            if set(item) - FIELD_KEYS:
                found.add((spec["_file"], item.get("wire") or item.get("name")))
    unreported = found - KNOWN_SPEC_DEFECTS
    assert not unreported, (
        f"Stray keys (an unquoted comma inside a flow mapping truncated the doc): {sorted(unreported)}. "
        "Quote the doc in the bridge's spec and re-sync; never edit docs/specs/ by hand."
    )
    fixed = KNOWN_SPEC_DEFECTS - found
    assert not fixed, f"Fixed in the bridge; delete from KNOWN_SPEC_DEFECTS: {sorted(fixed)}"


def test_every_declared_deviation_names_a_real_param_or_field():
    by_verb = {s["verb"]: s for s in SPECS}
    for verb, entries in spec_facts.PYTHON_DEVIATIONS.items():
        assert verb in by_verb, f"PYTHON_DEVIATIONS names '{verb}', which has no vendored spec"
        spec = by_verb[verb]
        names = {"param." + p["name"] for p in spec.get("params") or []}
        names |= {
            "field." + f["wire"] for f in (spec.get("result") or {}).get("envelope_fields") or []
        }
        names |= {"field." + f["wire"] for f in spec_facts.columns(spec)}
        for key, entry in entries.items():
            assert key in names, (
                f"PYTHON_DEVIATIONS['{verb}']['{key}'] names nothing in {spec['_file']}"
            )
            assert (entry.get("why") or "").strip(), (
                f"PYTHON_DEVIATIONS['{verb}']['{key}'] gives no reason"
            )


# --------------------------------------------------------------------------------------------
# The docstrings against the specs
# --------------------------------------------------------------------------------------------

PY_SPECS = [s for s in SPECS if spec_facts.py_binding(s)]


def _binding(spec):
    try:
        return resolve(spec_facts.py_binding(spec))
    except (ImportError, AttributeError):
        if spec_facts.shipped_in_python(spec):
            pytest.fail(
                f"{spec['_file']} says pyMzLib ships {spec_facts.py_binding(spec)} "
                f"(since.pymzlib {spec['since']['pymzlib']}), but it does not import"
            )
        pytest.skip(
            f"{spec_facts.py_binding(spec)} not written yet, and the spec does not claim it ships"
        )


@pytest.mark.parametrize("spec", PY_SPECS, ids=spec_id)
def test_every_param_is_documented_with_its_unit(spec):
    _, func = _binding(spec)
    signature = inspect.signature(func)
    args = sections(func).get("parameters", {})
    problems = []
    for p in spec.get("params") or []:
        declared = deviation(spec, "param." + p["name"])
        if declared and declared.get("python") is None:
            continue
        name = declared["python"] if declared else spec_facts.py_param_name(p["name"])
        if name not in signature.parameters:
            problems.append(
                f"'{name}' (wire --{p['name']}) is not a parameter of {func.__name__}()"
            )
            continue
        if name not in args:
            problems.append(f"'{name}' is missing from the Args: section")
        elif p.get("unit") and not unit_mentioned(p["unit"], args[name]):
            problems.append(
                f"'{name}' is documented without its unit '{p['unit']}': {args[name]!r}"
            )
    assert not problems, f"{spec['_file']} vs {spec_facts.py_binding(spec)}:\n  " + "\n  ".join(
        problems
    )


@pytest.mark.parametrize("spec", PY_SPECS, ids=spec_id)
def test_every_result_field_is_documented_with_its_unit(spec):
    module, func = _binding(spec)
    cls, is_list = returned(module, func)
    assert cls is not None, (
        f"cannot resolve the return annotation of {func.__name__}(); "
        "annotate it with the result class"
    )
    attributes = sections(cls).get("attributes", {})
    result = spec.get("result") or {}
    # A function returning a list of row objects projects the rows as the result; the envelope is
    # the list itself, so each envelope field must be a declared deviation.
    expected = spec_facts.columns(spec) if is_list else result.get("envelope_fields") or []
    problems = []
    if is_list:
        for f in result.get("envelope_fields") or []:
            if not deviation(spec, "field." + f["wire"]):
                problems.append(
                    f"{func.__name__}() returns a list, so envelope field '{f['wire']}' has "
                    "nowhere to be documented; declare it in PYTHON_DEVIATIONS"
                )
    for f in expected:
        declared = deviation(spec, "field." + f["wire"])
        if declared and declared.get("python") is None:
            continue
        name = declared["python"] if declared else f["wire"]
        if name not in attributes:
            problems.append(f"'{name}' is missing from {cls.__name__}'s Attributes: section")
        elif f.get("unit") and not unit_mentioned(f["unit"], attributes[name]):
            problems.append(
                f"{cls.__name__}.{name} is documented without its unit '{f['unit']}': "
                f"{attributes[name]!r}"
            )
    assert not problems, f"{spec['_file']} vs {cls.__name__}:\n  " + "\n  ".join(problems)


# --------------------------------------------------------------------------------------------
# The rendered tables
# --------------------------------------------------------------------------------------------


def test_rendered_fact_tables_are_current():
    problems = render_spec_docs.check()
    assert not problems, (
        "Stale reference fact tables; run python scripts/render_spec_docs.py:\n  "
        + "\n  ".join(problems)
    )


@pytest.mark.parametrize(
    "spec", [s for s in PY_SPECS if spec_facts.shipped_in_python(s)], ids=spec_id
)
def test_every_shipped_verb_has_its_table_on_a_reference_page(spec):
    include = f'--8<-- "docs/reference/_generated/{spec_facts.fragment_name(spec)}.md"'
    pages = sorted((ROOT / "docs" / "reference").glob("*.md"))
    assert any(include in page.read_text(encoding="utf-8") for page in pages), (
        f"No page under docs/reference/ includes {include}. Add it to docs/reference/"
        f"{spec['module']}.md (and that page to the nav in mkdocs.yml)."
    )


def test_examples_sections_are_named_so_they_render_as_code():
    # Google style's section is "Examples:". griffe reads a singular "Example:" as an admonition,
    # and the >>> lines inside it then render as nested blockquotes instead of code.
    source = ROOT / "pkg" / "python" / "src" / "pymzlib"
    offenders = [
        f"{path.name}:{n}"
        for path in sorted(source.glob("*.py"))
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if line.strip() == "Example:"
    ]
    assert not offenders, f"Rename 'Example:' to 'Examples:' at {offenders}"
