"""Tests for the mzLib 1.0.592 readers batch: the ``_many`` functions, the per-file facts, the four
ways a value can be missing, the spectra ``source`` block, mzIdentML's confidence fields, and the
three quantification views.

Offline. Every payload is recorded from the real bridge (``tests/fixtures``), and the Python layer's
jobs are pinned exactly: assemble the arguments and the stdin list, refuse a bad call before a
process starts, parse the payload into typed objects, and never loop over files itself.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, readers

FIXTURES = Path(__file__).parent / "fixtures"
VERSION = json.loads((FIXTURES / "bridge_version.json").read_text(encoding="utf-8"))


def load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture(autouse=True)
def _fresh_verb_cache():
    # require_verb caches the bridge's verb list per process; each test starts without it.
    _bridge._VERBS_SEEN.clear()
    yield
    _bridge._VERBS_SEEN.clear()


@pytest.fixture()
def bridge(monkeypatch):
    """A fake bridge: answers ``version`` from the recorded fixture and everything else from a
    payload the test chooses, and records every call it receives."""
    state = {"calls": [], "payload": None, "version": VERSION}

    def fake_invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": args, "stdin": stdin, "timeout": timeout})
        if args == ("version",):
            return state["version"]
        return state["payload"]

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return state


def reads(state) -> list[dict]:
    return [call for call in state["calls"] if call["args"] != ("version",)]


# ---- the single-file shape gains the per-file facts ---------------------------------------------


def test_read_results_names_msfraggers_decoy_flag_absent(bridge):
    bridge["payload"] = load("readers_results_fragger.json")

    result = readers.read_results("psm.tsv", limit=2)

    assert result.reader == "MsFraggerPsmFile"
    assert result.absent_fields == ["is_decoy"]
    assert set(result.columns["is_decoy"]) == {None}
    assert result.failed_fields == [] and result.excluded_fields == [] and result.error is None


def test_read_records_names_a_missing_mbr_score_column_absent(bridge):
    # mzLib #1345: current FlashLFQ peaks tables have no MBR Score column.
    bridge["payload"] = load("readers_records_flashlfq_peaks.json")

    result = readers.read_records("QuantifiedPeaks.tsv", limit=2)

    assert result.absent_fields == ["mbr_score"]
    assert set(result.columns["mbr_score"]) == {None}
    assert result.retention_time_unit is None
    assert result.skipped_count is None and result.skipped is None


def test_read_records_points_each_excluded_dictionary_at_its_verb(bridge):
    bridge["payload"] = load("readers_records_mm_protein_groups.json")

    result = readers.read_records("AllQuantifiedProteinGroups.tsv", limit=2)

    verbs = {entry["field"]: entry["verb"] for entry in result.excluded_fields}
    assert verbs["sample_groups"] == "readers read-protein-groups"


def test_read_spectra_reports_the_instrument_and_when_acquisition_started(bridge):
    bridge["payload"] = load("readers_spectra_mzml.json")

    result = readers.read_spectra("sliced_ethcd.mzML", limit=3)

    source = result.source
    assert source.instrument_model == "Orbitrap Fusion"
    assert source.instrument_model_accession == "MS:1002416"
    assert source.instrument_serial_number == "FSN10189"
    assert source.acquisition_start_time == "2021-03-16T17:09:07Z"
    assert source.acquisition_start_time_is_utc is True
    moment = source.acquired_at
    assert (moment.year, moment.hour, moment.utcoffset().total_seconds()) == (2021, 17, 0)


def test_a_local_acquisition_time_stays_naive():
    # Thermo records the instrument PC's clock with no offset; inventing one would be a fabrication.
    source = readers.SpectraSource(acquisition_start_time="2023-10-25T10:40:10.188556",
                                   acquisition_start_time_is_utc=False)

    assert source.acquired_at.tzinfo is None
    assert readers.SpectraSource().acquired_at is None


def test_read_matches_carries_mzidentml_confidence_and_its_skip_list(bridge):
    bridge["payload"] = load("readers_matches_mzid.json")

    result = readers.read_matches("run.mzid", limit=3)

    assert result.absent_fields == ["is_decoy"]
    assert all(isinstance(q, (int, float)) for q in result.columns["q_value"])
    assert result.columns["rank"] == [1, 2, 3], "every candidate rank is a row"
    assert result.skipped_count == 0 and result.skipped == []
    assert result.scores_included is False
    assert result.row_count == result.returned_count == 3


def test_read_matches_scores_are_long_rows_and_count_separately(bridge):
    bridge["payload"] = load("readers_matches_mzid_scores.json")

    result = readers.read_matches("run.mzid", limit=1, scores=True)

    assert reads(bridge)[0]["args"][-1] == "--scores"
    assert result.scores_included is True
    assert result.returned_count == 1, "returned_count counts matches, the unit limit counts in"
    assert result.row_count == len(result.columns["score_name"]) > 1
    assert set(result.columns["match_index"]) == {0}
    assert "MS-GF:SpecEValue" in result.columns["score_name"]


def test_read_matches_without_scores_sends_no_flag(bridge):
    bridge["payload"] = load("readers_matches_mzid.json")

    readers.read_matches("run.mzid")

    assert "--scores" not in reads(bridge)[0]["args"]


def test_a_casanovo_file_has_no_skip_list_rather_than_an_empty_one(bridge):
    bridge["payload"] = load("readers_matches_casanovo.json")

    result = readers.read_matches("run.mztab")

    assert result.skipped_count is None
    assert {"is_decoy", "q_value", "rank", "pass_threshold"} <= set(result.absent_fields)


def test_a_skipped_mzidentml_item_is_typed():
    item = readers.SkippedMatch._from_wire(
        {"spectrum_identification_item_id": "SII_1", "spectrum_id": "scan=1,scan=2", "reason": "crosslink identification"})

    assert item.reason == "crosslink identification"


# ---- the #1347 quantification views -------------------------------------------------------------


def test_read_protein_groups_is_long_with_one_row_per_sample_group(bridge):
    bridge["payload"] = load("readers_protein_groups.json")

    result = readers.read_protein_groups("AllQuantifiedProteinGroups.tsv", limit=1)

    assert reads(bridge)[0]["args"] == (
        "readers", "read-protein-groups", "--path", "AllQuantifiedProteinGroups.tsv", "--limit", "1")
    assert isinstance(result, readers.ProteinGroupRecords)
    assert result.returned_count == 1
    assert result.row_count == len(result.sample_labels) == 18
    assert result.columns["sample_label"] == result.sample_labels
    assert result.excluded_fields[0]["verb"] == "readers read-occupancy"
    assert any("UNFILTERED" in caveat for caveat in result.caveats)


def test_read_quantified_peptides_keeps_flashlfqs_zero_and_says_why(bridge):
    bridge["payload"] = load("readers_quantified_peptides.json")

    result = readers.read_quantified_peptides("AllQuantifiedPeptides.tsv", limit=1)

    zero_rows = [row for row in result.records if row["intensity"] == 0]
    assert zero_rows and all(row["detection_type"] != "MSMS" for row in zero_rows)
    assert result.absent_fields == ["peak_order", "retention_time"]
    assert set(result.columns["retention_time"]) == {None}
    assert result.retention_time_unit == "minutes"


def test_read_occupancy_rows_carry_their_basis_and_site(bridge):
    bridge["payload"] = load("readers_occupancy.json")

    result = readers.read_occupancy("AllQuantifiedProteinGroups.tsv")

    assert result.truncated_cell_count == 0
    assert set(result.columns["basis"]) == {"count", "intensity"}
    assert result.row_count == len(result.columns["basis"])
    names = set(result.columns["modification"])
    assert "N6,N6,N6-trimethyllysine on K" in names, "a comma inside a name survives the parse"


@pytest.mark.parametrize("function, verb", [
    (readers.read_protein_groups, "readers read-protein-groups"),
    (readers.read_quantified_peptides, "readers read-quantified-peptides"),
    (readers.read_occupancy, "readers read-occupancy"),
])
def test_a_quant_view_on_an_old_bridge_names_the_release_it_needs(bridge, function, verb):
    # A bridge from before the verbs list existed: every new verb is refused without spawning it.
    bridge["version"] = {"bridge": "1.0.0.0", "protocol": 1, "runtime": "10.0", "mzlib": None}

    with pytest.raises(pymzlib.UsageError, match=r"needs the bridge from pyMzLib 0\.2\.0 or later"):
        function("AllQuantifiedProteinGroups.tsv")

    assert reads(bridge) == [], "the verb itself must not be called"
    assert verb not in str(bridge["calls"][0]["args"])


def test_the_verb_check_costs_one_version_call_per_process(bridge):
    bridge["payload"] = load("readers_protein_groups.json")

    readers.read_protein_groups("a.tsv")
    readers.read_protein_groups("b.tsv")
    readers.read_occupancy("c.tsv")

    assert [call["args"] for call in bridge["calls"]].count(("version",)) == 1


def test_a_bridge_that_lists_other_verbs_but_not_this_one_is_refused(bridge):
    bridge["version"] = dict(VERSION, verbs=["version", "readers read-records"])

    with pytest.raises(pymzlib.UsageError, match="read-occupancy"):
        readers.read_occupancy("x.tsv")


def test_bridge_version_exposes_the_verbs(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: VERSION)

    info = pymzlib.bridge_version()

    assert "readers read-occupancy" in info["verbs"]
    assert "readers read-spectra" in info["verbs"]


# ---- the _many functions ------------------------------------------------------------------------


def test_many_sends_one_process_the_whole_list_on_stdin(bridge):
    bridge["payload"] = load("readers_many_spectra.json")
    paths = ["a.mzML", Path("b.mzML"), "c.mgf"]

    readers.read_spectra_many(paths, threads=2, on_error="skip", ms_order=2, peaks=True)

    calls = reads(bridge)
    assert len(calls) == 1, "N files is ONE bridge call, never a Python-side loop"
    assert calls[0]["args"] == (
        "readers", "read-spectra", "--paths-stdin", "--threads", "2", "--on-error", "skip",
        "--ms-order", "2", "--peaks")
    assert calls[0]["stdin"] == "a.mzML\nb.mzML\nc.mgf\n"


def test_many_defaults_are_one_thread_and_fail(bridge):
    bridge["payload"] = load("readers_many_records_mzid.json")

    readers.read_records_many(["a.mzid", "b.mzid.gz"], out="batch.tsv")

    assert reads(bridge)[0]["args"] == (
        "readers", "read-records", "--paths-stdin", "--threads", "1", "--on-error", "fail",
        "--out", "batch.tsv")


def test_many_parses_one_long_table_and_the_per_file_facts(bridge):
    bridge["payload"] = load("readers_many_spectra.json")

    batch = readers.read_spectra_many(["run.mzML", "missing.mzML", "run.mgf"], on_error="skip")

    assert isinstance(batch, readers.ReadBatch)
    assert batch.verb == "readers read-spectra"
    assert (batch.file_count, batch.read_count, batch.failed_count) == (3, 2, 1)
    assert batch.column_names[:2] == ["source_index", "source_path"]
    assert sorted(set(batch.columns["source_index"])) == [0, 2], "the unread file has no rows"
    assert batch.row_count == len(batch.columns["source_index"])
    assert [f.ok for f in batch.files] == [True, False, True]
    failed = batch.failed_files[0]
    assert failed.error.kind == "usage" and "File not found" in failed.error.message
    assert failed.record_count is None, "an unread file has no count, not a count of zero"
    assert batch.files[0].source.instrument_model == "Orbitrap Fusion"
    assert batch.files[2].source.instrument_model is None
    assert batch.files[0].scan_count == batch.files[0].record_count


def test_many_in_minutes_converts_each_files_rows_by_its_own_unit(bridge):
    bridge["payload"] = load("readers_many_spectra.json")

    batch = readers.read_spectra_many(["run.mzML", "missing.mzML", "run.mgf"], on_error="skip")

    assert batch.in_minutes("retention_time") == batch.columns["retention_time"]


def test_in_minutes_refuses_an_unknown_unit_rather_than_guessing():
    batch = readers.ReadBatch(
        verb="readers read-features", file_count=1, read_count=1, failed_count=0, on_error="fail",
        record_count=1, returned_count=1,
        columns={"source_index": [0], "retention_time_start": [12.0]},
        files=[readers.FileReport(path="x_ms1.feature", file_type="Ms1Feature", retention_time_unit="unknown")])

    with pytest.raises(pymzlib.UsageError, match="no basis to say"):
        batch.in_minutes("retention_time_start")


def test_many_records_table_is_the_shared_record_types(bridge):
    bridge["payload"] = load("readers_many_records_mzid.json")

    batch = readers.read_records_many(["a.mzid", "a.mzid.gz"], threads=2)

    assert {f.record_type for f in batch.files} == {"MzIdentMLRecord"}
    assert batch.record_count == sum(f.record_count for f in batch.files)
    assert [f.skipped_count for f in batch.files] == [0, 0]


def test_identify_many_answers_each_path_in_order(bridge):
    bridge["payload"] = load("readers_many_identify.json")

    batch = readers.identify_many(["a.mzid", "gone.mzML", "groups.tsv"], on_error="skip")

    assert isinstance(batch, readers.IdentifyBatch)
    assert [f.file_type for f in batch.files] == ["MzIdentML", "", "MetaMorpheusQuantifiedProteinGroups"]
    assert batch.failed_files[0].error.kind == "usage"
    assert reads(bridge)[0]["args"] == (
        "readers", "identify", "--paths-stdin", "--threads", "1", "--on-error", "skip")


@pytest.mark.parametrize("function", [
    readers.read_protein_groups_many, readers.read_quantified_peptides_many, readers.read_occupancy_many,
])
def test_the_quant_views_many_forms_check_the_bridge_too(bridge, function):
    bridge["version"] = {"bridge": "1.0.0.0", "protocol": 1}

    with pytest.raises(pymzlib.UsageError, match="0.2.0"):
        function(["a.tsv", "b.tsv"])


@pytest.mark.parametrize("paths, message", [
    ("one.mzML", "takes a LIST of paths"),
    (Path("one.mzML"), "takes a LIST of paths"),
    ([], "empty"),
    (["a.mzML", "  "], "Path 1 of the list is blank"),
    (["a.mzML", 3], "Path 1 of the list is blank"),
    (["a\nb.mzML"], "line break"),
    (b"a.mzML", "takes a LIST of paths"),
])
def test_many_refuses_a_bad_list_before_starting_a_process(bridge, paths, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        readers.read_spectra_many(paths)

    assert bridge["calls"] == []


@pytest.mark.parametrize("threads", [0, -2, True, 1.5, "4"])
def test_many_refuses_a_thread_count_that_means_nothing(bridge, threads):
    with pytest.raises(pymzlib.UsageError, match="threads must be"):
        readers.read_records_many(["a.mzid"], threads=threads)

    assert bridge["calls"] == []


def test_many_accepts_one_per_core(bridge):
    bridge["payload"] = load("readers_many_records_mzid.json")

    readers.read_records_many(["a.mzid"], threads=-1)

    assert reads(bridge)[0]["args"][3:5] == ("--threads", "-1")


def test_many_refuses_an_unknown_on_error(bridge):
    with pytest.raises(pymzlib.UsageError, match="on_error must be"):
        readers.read_matches_many(["a.mzid"], on_error="ignore")


def test_many_accepts_a_generator_of_paths(bridge):
    bridge["payload"] = load("readers_many_records_mzid.json")

    readers.read_records_many(p for p in ["a.mzid", "b.mzid"])

    assert reads(bridge)[0]["stdin"] == "a.mzid\nb.mzid\n"


def test_every_many_function_has_a_single_file_twin_and_the_same_defaults():
    import inspect

    for name in readers.__all__:
        if not name.endswith("_many"):
            continue
        single = getattr(readers, name[: -len("_many")])
        many = getattr(readers, name)
        params = inspect.signature(many).parameters
        assert params["threads"].default == 1, name
        assert params["on_error"].default == "fail", name
        assert "limit" not in params and "offset" not in params, f"{name}: a batch has no window"
        assert callable(single)


def test_the_single_file_functions_do_not_take_a_list():
    # The spelling decision: a list goes to <verb>_many, never to the single-path function.
    with pytest.raises(pymzlib.UsageError, match=r"read_records_many\(\)"):
        readers.read_records(["a.mzid", "b.mzid"])
