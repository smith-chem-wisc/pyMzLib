"""Tests for the median-polish protein-quantification entry point.

Offline throughout, in the same spirit as ``test_flashlfq``: the Python layer renders the optional
experimental design onto stdin, assembles the argument list, and parses the returned proteins into
typed objects. A recorded payload stands in for the bridge; the algorithm itself is mzLib's and is
covered by the bridge's own median-polish tests.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, flashlfq

FIXTURE = Path(__file__).parent / "fixtures" / "median_polish_small.json"


@pytest.fixture()
def recorded(monkeypatch):
    """Serve a recorded median-polish payload instead of running the bridge."""
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: payload)
    return payload


@pytest.fixture()
def captured(monkeypatch):
    """Capture the args and stdin a call would send, without running anything."""
    seen = {}
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))

    def fake_invoke(*args, stdin=None, timeout=None, **kwargs):
        seen["args"] = list(args)
        seen["stdin"] = stdin
        return payload

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# --------------------------------------------------------------------------- parsing


def test_returns_a_list_of_protein_groups(recorded):
    proteins = flashlfq.median_polish("QuantifiedPeptides.tsv")
    assert isinstance(proteins, list)
    assert len(proteins) == 3
    assert all(isinstance(p, flashlfq.ProteinGroup) for p in proteins)
    p1 = proteins[0]
    assert p1.protein_group == "P1"
    assert p1.gene_name == "GENE1"
    assert p1.intensity("control_1") == 3005.6
    assert p1.intensity("treated_1") == 6011.3


def test_unresolvable_protein_intensity_arrives_as_none(recorded):
    # P3 could not be resolved by median polish (degenerate matrix): its intensities are None, the
    # documented "could not be quantified" signal, distinct from a 0.0 "not measured".
    proteins = flashlfq.median_polish("QuantifiedPeptides.tsv")
    p3 = next(p for p in proteins if p.protein_group == "P3")
    assert p3.intensity("control_1") is None
    assert p3.intensity("treated_1") is None


def test_a_sample_never_quantified_reads_as_zero(recorded):
    # A sample label the protein does not carry defaults to 0.0, not a KeyError (ProteinGroup's
    # documented "missing is 0.0" behaviour).
    proteins = flashlfq.median_polish("QuantifiedPeptides.tsv")
    assert proteins[0].intensity("never_seen") == 0.0


# --------------------------------------------------------------------------- argument assembly


def test_minimal_call_passes_only_peptides(captured):
    flashlfq.median_polish("QuantifiedPeptides.tsv")
    assert captured["args"] == ["quant", "median-polish", "--peptides", "QuantifiedPeptides.tsv"]
    # No design given -> no stdin at all, so the bridge applies its each-run-its-own-replicate default.
    assert captured["stdin"] is None


def test_shared_peptides_and_out_become_flags(captured):
    flashlfq.median_polish(
        "QuantifiedPeptides.tsv", use_shared_peptides=True, output_directory="proteins"
    )
    args = captured["args"]
    assert "--shared-peptides" in args
    assert args[args.index("--out") + 1] == "proteins"


def test_shared_peptides_flag_absent_when_false(captured):
    flashlfq.median_polish("QuantifiedPeptides.tsv")
    assert "--shared-peptides" not in captured["args"]


def test_design_is_rendered_onto_stdin(captured):
    flashlfq.median_polish(
        "QuantifiedPeptides.tsv",
        design=[
            {"file_name": "run_1", "condition": "control", "biological_replicate": 0},
            {"file_name": "run_2", "condition": "treated", "biological_replicate": 0},
        ],
    )
    assert captured["stdin"] == "run_1\tcontrol\t0\nrun_2\ttreated\t0\n"


# --------------------------------------------------------------------------- stdin rendering


def test_design_stdin_bare_file_names_are_bare():
    assert flashlfq._design_stdin([{"file_name": "run_1"}, {"file_name": "run_2"}]) == "run_1\nrun_2\n"


def test_design_stdin_renders_all_fields():
    line = flashlfq._design_stdin([{"file_name": "run_1", "condition": "c",
                                    "biological_replicate": 2, "technical_replicate": 1, "fraction": 3}])
    assert line == "run_1\tc\t2\t1\t3\n"


def test_design_stdin_trims_trailing_empty_fields():
    assert flashlfq._design_stdin([{"file_name": "run_1", "condition": "c"}]) == "run_1\tc\n"


def test_design_stdin_keeps_empty_middle_field():
    line = flashlfq._design_stdin([{"file_name": "run_1", "fraction": 3}])
    assert line == "run_1\t\t\t\t3\n"


# --------------------------------------------------------------------------- validation


def test_blank_peptides_path_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("   ")


def test_design_must_be_a_list_not_a_string():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("QuantifiedPeptides.tsv", design="run_1")


def test_design_entry_without_file_name_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("QuantifiedPeptides.tsv", design=[{"condition": "c"}])


def test_design_entry_must_be_a_mapping():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("QuantifiedPeptides.tsv", design=["run_1"])


def test_tab_in_file_name_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("QuantifiedPeptides.tsv", design=[{"file_name": "a\tb"}])


def test_negative_replicate_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("QuantifiedPeptides.tsv", design=[{"file_name": "run_1", "biological_replicate": -1}])


def test_empty_output_directory_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.median_polish("QuantifiedPeptides.tsv", output_directory="   ")
