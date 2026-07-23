"""Tests for the FlashLFQ quantification workflow.

Offline throughout: the Python layer's job is to render the spectra design onto the wire, assemble
the argument list, and parse the results back into typed objects. None of that needs the bridge or
real mzML — a recorded payload stands in for the engine, and the argument assembly is asserted
directly. The engine itself is covered by mzLib's own FlashLFQ tests.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, flashlfq

FIXTURE = Path(__file__).parent / "fixtures" / "flashlfq_small.json"


@pytest.fixture()
def recorded(monkeypatch):
    """Serve a recorded quantification payload instead of running the bridge."""
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


def test_results_parse_into_typed_objects(recorded):
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    assert isinstance(result, flashlfq.FlashLfqResults)
    assert result.identification_count == 4
    assert result.peptide_count == 2
    assert result.protein_count == 2
    assert all(isinstance(p, flashlfq.Peptide) for p in result.peptides)
    assert all(isinstance(g, flashlfq.ProteinGroup) for g in result.proteins)
    assert all(isinstance(f, flashlfq.SpectraFileInfo) for f in result.spectra_files)


def test_peptide_intensity_and_detection_lookup(recorded):
    peptide = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"]).peptides[0]
    assert peptide.sequence == "PEPTIDEK"
    assert peptide.protein_groups == "P12345"
    assert peptide.intensity("run_3") == 1000.0
    assert peptide.detection_type("run_4") == "MBR"
    # A run that was never provided reads as not detected, not a KeyError.
    assert peptide.detection_type("never_seen") == "NotDetected"
    # Missing peptide intensity is 0.0, never None (None is proteins-only).
    assert peptide.intensity("never_seen") == 0.0


def test_nan_protein_intensity_arrives_as_none(recorded):
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    unquantified = next(g for g in result.proteins if g.protein_group == "P67890")
    assert unquantified.intensity("run_3") is None
    assert unquantified.intensity("run_4") is None


def test_mbr_peak_count_sums_across_runs(recorded):
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    assert result.mbr_peak_count == 1  # 1 in run_3, 0 in run_4


def test_spectra_file_fields_parse(recorded):
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    run_4 = result.spectra_files[1]
    assert run_4.file_name == "run_4"
    assert run_4.condition == "treated"
    assert run_4.biological_replicate == 1
    assert run_4.peak_count == 2


# --------------------------------------------------------------------------- peaks (the MBR surface)


def test_peaks_parse_into_typed_objects(recorded):
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    assert all(isinstance(pk, flashlfq.Peak) for pk in result.peaks)
    assert len(result.peaks) == 3
    peak = result.peaks[0]
    assert peak.file_name == "run_3"
    assert peak.sequence == "PEPTIDEK"
    assert peak.intensity == 1000.0
    assert peak.detection_type == "MSMS"
    assert peak.retention_time == 30.1
    assert peak.is_mbr is False


def test_mbr_peaks_are_the_transferred_ones(recorded):
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    mbr = result.mbr_peaks
    assert len(mbr) == 1
    assert mbr[0].is_mbr is True
    assert mbr[0].file_name == "run_4"
    assert mbr[0].sequence == "PEPTIDEK"


def test_peaks_carry_mbr_the_peptide_table_drops(recorded):
    # The finding the bake-off surfaced: an MBR transfer visible in the peaks need not appear as
    # "MBR" in the peptide roll-up. The peaks are the authoritative surface.
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    assert any(pk.is_mbr for pk in result.peaks)


def test_mbr_rescued_peptide_count_is_distinct_sequences(recorded):
    # One MBR peak in the fixture (PEPTIDEK in run_4) -> exactly one distinct rescued peptide.
    result = flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML", "run_4.mzML"])
    assert result.mbr_rescued_peptide_count == 1


# --------------------------------------------------------------------------- argument assembly


def test_quantify_passes_expected_args_and_stdin(captured):
    flashlfq.quantify(
        "AllPSMs.psmtsv",
        ["run_3.mzML", {"path": "run_4.mzML", "condition": "treated", "biological_replicate": 1}],
        match_between_runs=True,
        normalize=True,
        ppm_tolerance=8.0,
    )
    args = captured["args"]
    assert args[:4] == ["quant", "flashlfq", "--psms", "AllPSMs.psmtsv"]
    assert "--mbr" in args
    assert "--normalize" in args
    # Numeric options carry their value as the next token, invariant-formatted.
    assert args[args.index("--ppm") + 1] == "8.0"
    # Bare path stays bare; a design mapping renders only the fields it set.
    assert captured["stdin"] == "run_3.mzML\nrun_4.mzML\ttreated\t1\n"


def test_flags_absent_when_false(captured):
    flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML"])
    args = captured["args"]
    assert "--mbr" not in args
    assert "--normalize" not in args
    assert "--bayesian" not in args


def test_output_directory_becomes_out_flag(captured):
    flashlfq.quantify("AllPSMs.psmtsv", ["run_3.mzML"], output_directory="results")
    args = captured["args"]
    assert args[args.index("--out") + 1] == "results"


# --------------------------------------------------------------------------- stdin rendering


def test_spectra_stdin_bare_paths_are_bare():
    assert flashlfq._spectra_stdin(["a.mzML", "b.mzML"]) == "a.mzML\nb.mzML\n"


def test_spectra_stdin_renders_design_fields():
    line = flashlfq._spectra_stdin([{"path": "a.mzML", "condition": "c", "biological_replicate": 2,
                                     "technical_replicate": 1, "fraction": 3}])
    assert line == "a.mzML\tc\t2\t1\t3\n"


def test_spectra_stdin_trims_trailing_empty_fields():
    # condition set, replicate/fraction omitted -> only path and condition on the wire.
    assert flashlfq._spectra_stdin([{"path": "a.mzML", "condition": "c"}]) == "a.mzML\tc\n"


def test_spectra_stdin_keeps_empty_middle_field():
    # fraction set but condition/replicate omitted -> empties preserved as positional placeholders.
    line = flashlfq._spectra_stdin([{"path": "a.mzML", "fraction": 3}])
    assert line == "a.mzML\t\t\t\t3\n"


# --------------------------------------------------------------------------- validation


def test_empty_spectra_raises_before_invoke():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", [])


def test_blank_psms_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("   ", ["a.mzML"])


def test_spectra_must_be_a_list_not_a_string():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", "a.mzML")


def test_spectra_mapping_without_path_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", [{"condition": "c"}])


def test_tab_in_path_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", ["a\tb.mzML"])


def test_negative_replicate_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", [{"path": "a.mzML", "biological_replicate": -1}])


def test_non_finite_tolerance_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", ["a.mzML"], ppm_tolerance=math.nan)


def test_boolean_max_threads_raises():
    with pytest.raises(pymzlib.UsageError):
        flashlfq.quantify("AllPSMs.psmtsv", ["a.mzML"], max_threads=True)
