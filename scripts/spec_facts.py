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
    "sdrf design": {
        "param.searched-files-stdin": {
            "python": "searched_files",
            "why": "the flag is set, and the stdin lines rendered, from the searched_files list",
        },
        "param.condition-columns": {
            "python": "condition_columns",
            "why": "a list in Python; pyMzLib joins it with tabs for the wire",
        },
    },
    "isobaric kits": {
        "field.kits": {
            "python": "kit_summaries",
            "why": "the wire's {kit, channel_count} list; IsobaricKits.kits regroups the table under it",
        },
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
    "proteins read": {
        # The protein verbs take one database or a list as their first argument, and contaminant
        # databases as a separate list, so the wire's --path / --paths-stdin / --contaminant choice
        # is made by pymzlib.proteins from the arguments' shape rather than spelled by the caller.
        "param.path": {
            "python": "databases",
            "why": "one database or a list; the list travels as --paths-stdin",
        },
        "param.contaminant": {
            "python": "contaminants",
            "why": "contaminant databases are their own list, tagged on stdin or --contaminant",
        },
        "param.paths-stdin": {
            "python": None,
            "why": "chosen automatically when more than one database is given",
        },
        "param.accessions-stdin": {
            "python": "accessions",
            "why": "the accession list itself; its presence sets the flag and fills stdin",
        },
    },
    "genes resolve": {
        # The protein verbs take one database or a list as their first argument, and contaminant
        # databases as a separate list, so the wire's --path / --paths-stdin / --contaminant choice
        # is made by pymzlib.proteins from the arguments' shape rather than spelled by the caller.
        "param.path": {
            "python": "databases",
            "why": "one database or a list; the list travels as --paths-stdin",
        },
        "param.contaminant": {
            "python": "contaminants",
            "why": "contaminant databases are their own list, tagged on stdin or --contaminant",
        },
        "param.paths-stdin": {
            "python": None,
            "why": "chosen automatically when more than one database is given",
        },
    },
    "proteins classify-peptides": {
        # The protein verbs take one database or a list as their first argument, and contaminant
        # databases as a separate list, so the wire's --path / --paths-stdin / --contaminant choice
        # is made by pymzlib.proteins from the arguments' shape rather than spelled by the caller.
        "param.path": {
            "python": "databases",
            "why": "one database or a list; the list travels as --paths-stdin",
        },
        "param.contaminant": {
            "python": "contaminants",
            "why": "contaminant databases are their own list, tagged on stdin or --contaminant",
        },
        "param.paths-stdin": {
            "python": None,
            "why": "chosen automatically when more than one database is given",
        },
        "param.on-error": {
            "python": None,
            "why": "the only accepted value is fail, the default; there is no skip to choose",
        },
    },
    "stats fit": {
        "param.stdin": {"python": "coefficients", "why": "one stdin line per coefficient to test"},
    },
    "stats adjust": {
        "param.stdin": {"python": "p_values", "why": "one stdin line per p-value; None is a blank line"},
    },
    "stats meta": {
        "param.stdin": {"python": "studies", "why": "one stdin line per (feature, estimate, standard_error)"},
    },
    "peptidoform fragments": {
        # The wire's negative flag reads as a positive keyword in Python, and the wire's census
        # fields arrive grouped as one object, Digest.modification_census.
        "param.no-modifications": {
            "python": "modifications",
            "why": "modifications=False sends --no-modifications",
        },
        "param.max-mods": {
            "python": "max_modifications",
            "why": "spelled out; the wire's abbreviation is not worth keeping in Python",
        },
        "field.annotated_modification_sites": {
            "python": "modification_census",
            "why": "ModificationCensus.sites",
        },
        "field.annotated_modifications_loaded": {
            "python": "modification_census",
            "why": "ModificationCensus.applied",
        },
        "field.uniprot_annotated_features": {
            "python": "modification_census",
            "why": "ModificationCensus.annotated",
        },
        "field.unresolved_modifications": {
            "python": "modification_census",
            "why": "ModificationCensus.unresolved",
        },
        "field.uniprot_features_by_type": {
            "python": "modification_census",
            "why": "ModificationCensus.by_type",
        },
        "field.max_modification_isoforms": {
            "python": "max_isoforms",
            "why": "named after the max_isoforms argument it echoes",
        },
        "field.peptides_at_isoform_cap": {
            "python": "peptides_at_cap",
            "why": "shortened; Digest.truncated is the yes/no reading of it",
        },
        "field.peptide_count": {"python": None, "why": "len(Digest.peptides)"},
    },
    # The PRIDE functions return the rows themselves (a list of files, of paths, of hits), so each
    # envelope field is the list, its length, or an argument echoed back.
    "pride download": {
        "param.dest": {"python": "destination", "why": "spelled out"},
        "param.ext": {
            "python": "extensions",
            "why": "a list in Python; pyMzLib joins it with commas for the wire",
        },
        "param.no-overwrite": {
            "python": "overwrite",
            "why": "overwrite=False sends --no-overwrite",
        },
        "param.names-from-stdin": {
            "python": None,
            "why": "download_files(files) is the select-by-name form",
        },
        "param.stdin": {"python": None, "why": "download_files() renders it from its files list"},
        "field.accession": {"python": None, "why": "the accession argument, echoed"},
        "field.destination_directory": {
            "python": None,
            "why": "the destination argument; each returned Path is inside it",
        },
        "field.downloaded_count": {"python": None, "why": "len() of the returned list"},
        "field.paths": {"python": None, "why": "download() returns this list itself, as Paths"},
    },
    "pride files": {
        "field.accession": {"python": None, "why": "PrideFile.project_accession on every row"},
        "field.file_count": {"python": None, "why": "len(list_files(...))"},
        "field.total_size_bytes": {
            "python": None,
            "why": "pymzlib.pride.total_size_bytes(files), mzLib's TotalSizeBytes",
        },
        "field.files": {"python": None, "why": "list_files() returns this list itself"},
    },
    "pride ftp-files": {
        "field.accession": {"python": None, "why": "PrideFtpFile.project_accession on every row"},
        "field.file_count": {"python": None, "why": "len(list_ftp_files(...))"},
        "field.approximate_total_size_bytes": {
            "python": None,
            "why": "pymzlib.pride.approximate_total_size_bytes(files)",
        },
        "field.files": {"python": None, "why": "list_ftp_files() returns this list itself"},
    },
    "pride search": {
        "field.keyword": {"python": None, "why": "the keyword argument, echoed"},
        "field.result_count": {"python": None, "why": "len(search(...))"},
        "field.results": {"python": None, "why": "search() returns this list itself"},
    },
    # FlashLFQ's keywords are mzLib's FlashLfqEngine parameter names, so a Python call reads like
    # the C# one; the wire's short flags abbreviate the same names.
    "quant flashlfq": {
        "param.stdin": {"python": "spectra", "why": "the stdin lines are rendered from spectra"},
        "param.ppm": {"python": "ppm_tolerance", "why": "mzLib's PpmTolerance"},
        "param.isotope-ppm": {
            "python": "isotope_ppm_tolerance",
            "why": "mzLib's IsotopePpmTolerance",
        },
        "param.mbr": {"python": "match_between_runs", "why": "mzLib's MatchBetweenRuns"},
        "param.mbr-ppm": {"python": "mbr_ppm_tolerance", "why": "mzLib's MbrPpmTolerance"},
        "param.mbr-q": {"python": "mbr_q_value_threshold", "why": "mzLib's MbrQValueThreshold"},
        "param.shared-peptides": {
            "python": "use_shared_peptides_for_protein_quant",
            "why": "mzLib's UseSharedPeptidesForProteinQuant",
        },
        "param.bayesian": {
            "python": "bayesian_protein_quant",
            "why": "mzLib's BayesianProteinQuant",
        },
        "param.use-pep-q": {"python": "use_pep_q_value", "why": "MakeIdentifications usePepQValue"},
        "param.threads": {"python": "max_threads", "why": "mzLib's MaxThreads"},
        "param.out": {"python": "output_directory", "why": "spelled out"},
    },
    "quant median-polish": {
        "param.stdin": {"python": "design", "why": "the stdin lines are rendered from design"},
        "param.shared-peptides": {
            "python": "use_shared_peptides",
            "why": "mzLib's UseSharedPeptidesForProteinQuant, shortened",
        },
        "param.out": {"python": "output_directory", "why": "spelled out"},
        **{
            f"field.{name}": {
                "python": None,
                "why": "not projected: median_polish() returns only the protein list",
            }
            for name in ("peptides_file", "parameters", "samples", "peptide_count", "output_directory")
        },
        "field.protein_count": {"python": None, "why": "len(median_polish(...))"},
        "field.proteins": {"python": None, "why": "median_polish() returns this list itself"},
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
