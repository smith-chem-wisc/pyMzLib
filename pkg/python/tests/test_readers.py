"""Tests for result-file identification.

Offline throughout. The Python layer's job here is small and worth pinning exactly: assemble the
argument list, parse the wire payload into typed objects, and answer "can this file be quantified"
without the caller having to know that the answer is a string in a list. The dispatch itself is
mzLib's and is covered by the bridge's C# suite.

The ``readers formats`` fixture is recorded from the real bridge rather than hand-written, so it
reflects what mzLib actually dispatches.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, readers

FORMATS_FIXTURE = Path(__file__).parent / "fixtures" / "readers_formats.json"


@pytest.fixture()
def recorded_formats(monkeypatch):
    """Serve the recorded ``readers formats`` payload instead of running the bridge."""
    payload = json.loads(FORMATS_FIXTURE.read_text(encoding="utf-8"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: payload)
    return payload


@pytest.fixture()
def captured(monkeypatch):
    """Capture the args a call would send, without running anything."""
    seen: dict = {}

    def fake_invoke(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return {"path": "C:/x/AllPSMs.psmtsv", "file_type": "psmtsv", "extension": ".psmtsv",
                "reader": "PsmFromTsvFile", "views": ["quantifiable"]}

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# ---- identify ------------------------------------------------------------------------------


def test_identify_parses_the_payload_into_a_typed_result(captured):
    info = readers.identify("AllPSMs.psmtsv")

    assert isinstance(info, readers.FileInfo)
    assert info.file_type == "psmtsv"
    assert info.reader == "PsmFromTsvFile"
    assert info.views == ["quantifiable"]


def test_identify_sends_the_expected_verb_and_option(captured):
    readers.identify("  AllPSMs.psmtsv  ")

    assert captured["args"] == ("readers", "identify", "--path", "AllPSMs.psmtsv")


def test_is_quantifiable_is_true_only_for_the_quantifiable_view(captured):
    assert readers.identify("AllPSMs.psmtsv").is_quantifiable is True


def test_a_file_with_no_uniform_view_is_not_quantifiable(monkeypatch):
    # TopPIC is read perfectly well by mzLib and has no cross-format view. A caller must be able to
    # find that out here rather than discover it as a failure inside flashlfq.quantify().
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "path": "C:/x/run_prsm.tsv", "file_type": "ToppicPrsm", "extension": "_prsm.tsv",
        "reader": "ToppicSearchResultFile", "views": [],
    })

    info = readers.identify("run_prsm.tsv")

    assert info.views == []
    assert info.is_quantifiable is False


def test_a_spectra_file_is_not_quantifiable(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "path": "C:/x/run.mzML", "file_type": "MzML", "extension": ".mzML",
        "reader": "MsDataFileToResultFileAdapter", "views": ["spectra"],
    })

    assert readers.identify("run.mzML").is_quantifiable is False


def test_missing_wire_fields_do_not_crash_the_parser(monkeypatch):
    # The wire contract is versioned but a consumer should degrade rather than raise on a payload
    # that predates a field.
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"file_type": "psmtsv"})

    info = readers.identify("x.psmtsv")

    assert info.views == []
    assert info.is_quantifiable is False
    assert info.extension is None


@pytest.mark.parametrize("bad", ["", "   ", None, 42])
def test_a_blank_or_non_string_path_is_rejected_without_starting_a_process(bad):
    with pytest.raises(pymzlib.UsageError):
        readers.identify(bad)


# ---- formats -------------------------------------------------------------------------------


def test_formats_parses_every_entry(recorded_formats):
    formats = readers.formats()

    assert len(formats) == recorded_formats["format_count"]
    assert all(isinstance(f, readers.Format) for f in formats)
    assert all(f.file_type for f in formats)


def test_exactly_three_formats_offer_the_quantifiable_view(recorded_formats):
    # The headline fact about this module, pinned on the Python side too: mzLib reads 29 formats
    # and only these three can feed flashlfq.quantify(). Documented in the module docstring, so if
    # mzLib widens the set this test fails and the docs get corrected rather than quietly lying.
    quantifiable = [f.file_type for f in readers.formats() if f.is_quantifiable]

    assert sorted(quantifiable) == sorted(["psmtsv", "osmtsv", "MsFraggerPsm"])


def test_most_formats_have_no_uniform_view(recorded_formats):
    viewless = [f for f in readers.formats() if not f.views]

    assert len(viewless) > 10, "an empty view list is the common case, not an edge case"


def test_extensions_are_not_unique_across_formats(recorded_formats):
    # Documented in Format.extension and worth pinning: both Bruker types are ".d", so extension is
    # not a key. Anyone tempted to build a dict keyed on it should find this test instead of a bug.
    extensions = [f.extension for f in readers.formats()]

    assert len(extensions) != len(set(extensions))


def test_formats_sends_no_options(monkeypatch):
    seen = {}
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: seen.update(args=a) or {"formats": []})

    readers.formats()

    assert seen["args"] == ("readers", "formats")


# ---- module wiring -------------------------------------------------------------------------


def test_readers_is_exported_from_the_package():
    # A module missing from __all__ is invisible to `from pymzlib import *` and to the docs build.
    assert "readers" in pymzlib.__all__
    assert pymzlib.readers is readers


# ---- read_results --------------------------------------------------------------------------


READ_PAYLOAD = {
    "path": "C:/x/AllPSMs.psmtsv",
    "file_type": "psmtsv",
    "record_count": 3,
    "returned_count": 2,
    "offset": 0,
    "truncated": True,
    "rows_not_read": 0,
    "caveats": ["only the FIRST candidate of an ambiguous identification is kept"],
    "column_names": ["base_sequence", "retention_time", "is_decoy"],
    "columns": {
        "base_sequence": ["PEPTIDE", "SEQUENCE"],
        "retention_time": [12.5, None],
        "is_decoy": [False, True],
    },
    "output": None,
}


@pytest.fixture()
def recorded_read(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: READ_PAYLOAD)
    return READ_PAYLOAD


def test_read_results_parses_the_columnar_payload(recorded_read):
    result = readers.read_results("AllPSMs.psmtsv")

    assert result.record_count == 3
    assert result.returned_count == 2
    assert result.columns["base_sequence"] == ["PEPTIDE", "SEQUENCE"]
    assert result.column_names == ["base_sequence", "retention_time", "is_decoy"]


def test_truncation_is_visible(recorded_read):
    # A short answer and a complete one must not look alike.
    assert readers.read_results("AllPSMs.psmtsv").truncated is True


def test_caveats_are_carried_through(recorded_read):
    assert readers.read_results("AllPSMs.psmtsv").caveats


def test_records_gives_the_same_data_row_wise(recorded_read):
    rows = readers.read_results("AllPSMs.psmtsv").records

    assert len(rows) == 2, "one dict per RETURNED record, not per record in the file"
    assert rows[0] == {"base_sequence": "PEPTIDE", "retention_time": 12.5, "is_decoy": False}
    assert rows[1]["retention_time"] is None, "mzLib's -1 sentinel must arrive as None, not -1.0"


def test_records_is_empty_when_the_table_was_written_to_disk(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "path": "C:/x/AllPSMs.psmtsv", "file_type": "psmtsv", "record_count": 8,
        "returned_count": 0, "offset": 0, "truncated": False, "columns": None,
        "output": {"path": "C:/out/records.tsv", "format": "tsv", "row_count": 8},
    })

    result = readers.read_results("AllPSMs.psmtsv", out="records.tsv")

    assert result.columns is None
    assert result.records == []
    assert isinstance(result.output, readers.WrittenTable)
    assert result.output.format == "tsv", "tab-separated: these fields contain commas"
    assert result.output.row_count == 8


def test_limit_and_offset_and_out_are_sent(monkeypatch):
    seen = {}
    monkeypatch.setattr(_bridge, "invoke",
                        lambda *a, **k: seen.update(args=a) or dict(READ_PAYLOAD))

    readers.read_results("AllPSMs.psmtsv", limit=5, offset=10, out=" records.tsv ")

    assert seen["args"] == (
        "readers", "read-results", "--path", "AllPSMs.psmtsv",
        "--limit", "5", "--offset", "10", "--out", "records.tsv",
    )


def test_defaults_send_no_limit_or_offset(monkeypatch):
    # There is no default row cap: the ordinary call must ask for the whole file.
    seen = {}
    monkeypatch.setattr(_bridge, "invoke",
                        lambda *a, **k: seen.update(args=a) or dict(READ_PAYLOAD))

    readers.read_results("AllPSMs.psmtsv")

    assert seen["args"] == ("readers", "read-results", "--path", "AllPSMs.psmtsv")


def test_zero_limit_is_sent_rather_than_treated_as_absent(monkeypatch):
    seen = {}
    monkeypatch.setattr(_bridge, "invoke",
                        lambda *a, **k: seen.update(args=a) or dict(READ_PAYLOAD))

    readers.read_results("AllPSMs.psmtsv", limit=0)

    assert "--limit" in seen["args"], "limit=0 is a request for no rows, not a missing argument"


@pytest.mark.parametrize("bad", [-1, 1.5, True, "5"])
def test_a_bad_limit_is_rejected_without_starting_a_process(bad):
    with pytest.raises(pymzlib.UsageError):
        readers.read_results("AllPSMs.psmtsv", limit=bad)


@pytest.mark.parametrize("bad", [-1, 2.5, True, "3"])
def test_a_bad_offset_is_rejected_without_starting_a_process(bad):
    with pytest.raises(pymzlib.UsageError):
        readers.read_results("AllPSMs.psmtsv", offset=bad)


@pytest.mark.parametrize("bad", ["", "   ", 7])
def test_a_bad_out_path_is_rejected_without_starting_a_process(bad):
    with pytest.raises(pymzlib.UsageError):
        readers.read_results("AllPSMs.psmtsv", out=bad)


def test_a_blank_read_path_is_rejected():
    with pytest.raises(pymzlib.UsageError):
        readers.read_results("  ")


# ---- bake-off round 1 fixes ------------------------------------------------------------------


def test_retention_time_unit_is_a_value_not_prose(monkeypatch):
    # Round-1 finding: with the unit stated only in a prose caveat, a reader had to grep a sentence
    # for the word "SECONDS" and hard-code a units table by hand.
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: dict(READ_PAYLOAD, **{
        "file_type": "MsFraggerPsm", "retention_time_unit": "seconds",
        "columns": {"retention_time": [60.0, 120.0, None]},
        "column_names": ["retention_time"], "returned_count": 3,
    }))

    result = readers.read_results("psm.tsv")

    assert result.retention_time_unit == "seconds"
    assert result.retention_time_in_minutes == [1.0, 2.0, None]


def test_retention_time_in_minutes_passes_minutes_through(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: dict(READ_PAYLOAD, **{
        "retention_time_unit": "minutes",
        "columns": {"retention_time": [12.5]}, "column_names": ["retention_time"],
        "returned_count": 1,
    }))

    assert readers.read_results("AllPSMs.psmtsv").retention_time_in_minutes == [12.5]


def test_retention_time_in_minutes_refuses_to_guess_an_unknown_unit(monkeypatch):
    # Guessing is the exact failure this module exists to prevent: a silently unconverted axis.
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: dict(READ_PAYLOAD, **{
        "retention_time_unit": "unknown",
        "columns": {"retention_time": [1.0]}, "column_names": ["retention_time"],
        "returned_count": 1,
    }))

    with pytest.raises(pymzlib.UsageError):
        readers.read_results("x.psmtsv").retention_time_in_minutes


def test_missing_unit_on_the_wire_defaults_to_unknown_not_minutes(recorded_read):
    assert readers.read_results("AllPSMs.psmtsv").retention_time_unit == "unknown"


def test_dir_shows_the_api_and_not_the_modules_imports():
    # Round-1 finding: dir() listed Any, annotations, dataclass and field beside the functions -
    # four false leads out of eleven names, on a package whose discovery story is dir()/help().
    names = dir(readers)

    assert "identify" in names and "read_results" in names and "formats" in names
    for leaked in ("Any", "annotations", "dataclass", "field", "_bridge"):
        assert leaked not in names


def test_docstrings_are_ascii_so_help_is_readable_on_a_windows_console():
    # Round-1 finding: every em-dash and arrow rendered as "?" in a cp1252 console, making the
    # primary discovery surface look broken on exactly the machines this audience uses.
    import pymzlib.readers as module

    targets = [module, module.identify, module.read_results, module.formats,
               module.ResultRecords, module.FileInfo, module.Format, module.WrittenTable]
    for target in targets:
        text = target.__doc__ or ""
        assert text.isascii(), f"{getattr(target, '__name__', target)} docstring has non-ASCII"
