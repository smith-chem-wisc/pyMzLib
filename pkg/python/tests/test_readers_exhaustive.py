"""Tests for the four verbs that make readers coverage exhaustive.

``read_results`` reaches three of mzLib's twenty-nine file types. These four reach the rest:
``read_records`` reads any of them into that format's own fields, and ``read_features``,
``read_matches`` and ``read_spectra`` project the three typed views that ``read_results`` is not.

Offline throughout. The payloads are **recorded from the real bridge** against real mzLib fixtures
rather than hand-written, so a wire-shape change shows up here as a parse failure rather than as a
fixture that agrees with a Python file and with nothing else. What is under test is only the
Python layer's job: assemble the arguments, parse the payload into typed objects, and refuse a bad
argument before spawning anything. The reading itself is mzLib's, and the bridge's C# suite proves
all twenty-nine types are reachable.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymzlib import _bridge, readers

FIXTURES = Path(__file__).parent / "fixtures"


def payload(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def recorded(monkeypatch):
    """Serve a recorded payload, and capture the arguments the call would have sent."""
    seen: dict = {}

    def serve(name: str):
        data = payload(name)

        def fake_invoke(*args, **kwargs):
            seen["args"] = args
            seen["kwargs"] = kwargs
            return data

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        return seen

    return serve


@pytest.fixture()
def captured(monkeypatch):
    """Capture the arguments a call would send, returning an empty payload."""
    seen: dict = {}

    def fake_invoke(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return {}

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# ---- read_records --------------------------------------------------------------------------


def test_read_records_parses_a_format_with_no_uniform_view(recorded):
    recorded("readers_records_toppic.json")

    result = readers.read_records("prsm.tsv")

    assert isinstance(result, readers.NativeRecords)
    assert result.file_type == "ToppicPrsm"
    assert result.record_type == "ToppicPrsm"
    # The whole point: TopPIC has no cross-format view at all, and is still readable.
    assert result.views == []
    assert "e_value" in result.column_names


def test_read_records_column_names_are_this_formats_own_not_a_uniform_set(recorded):
    recorded("readers_records_crux.json")
    crux = readers.read_records("crux.txt")

    recorded("readers_records_toppic.json")
    toppic = readers.read_records("prsm.tsv")

    # Two formats, two entirely different column sets. A caller who assumes otherwise is the
    # failure mode this verb's documentation exists to prevent.
    assert set(crux.column_names) != set(toppic.column_names)
    assert "x_corr_score" in crux.column_names
    assert "x_corr_score" not in toppic.column_names


def test_read_records_names_fields_that_could_not_become_columns(recorded):
    recorded("readers_records_toppic.json")

    result = readers.read_records("prsm.tsv")

    names = [field["field"] for field in result.excluded_fields]
    assert "alternative_identifications" in names
    # Every exclusion carries its reason: a column that simply vanished would be
    # indistinguishable from a field the format does not have.
    assert all(entry.get("reason") for entry in result.excluded_fields)


def test_read_records_rows_view_matches_the_columns(recorded):
    recorded("readers_records_crux.json")

    result = readers.read_records("crux.txt")
    rows = result.records

    assert len(rows) == result.returned_count
    assert rows[0]["x_corr_score"] == result.columns["x_corr_score"][0]
    # Row order follows column_names, so a table built from rows and one built from columns agree.
    assert list(rows[0]) == [n for n in result.column_names if n in result.columns]


def test_read_records_sends_the_expected_verb_and_options(captured):
    readers.read_records("  a.tsv  ", limit=5, offset=2, out=" out.tsv ")

    assert captured["args"] == (
        "readers", "read-records", "--path", "a.tsv",
        "--limit", "5", "--offset", "2", "--out", "out.tsv",
    )


# ---- read_features -------------------------------------------------------------------------


def test_read_features_parses_the_ms1_feature_view(recorded):
    recorded("readers_features_topfd.json")

    result = readers.read_features("sample_ms1.feature")

    assert isinstance(result, readers.FeatureRecords)
    assert result.column_names == [
        "mz", "charge", "retention_time_start", "retention_time_end",
        "intensity", "number_of_isotopes",
    ]


def test_read_features_reports_an_unknown_retention_time_unit_for_ms1_feature(recorded):
    recorded("readers_features_topfd.json")

    result = readers.read_features("sample_ms1.feature")

    # Not a gap in this library: TopFD wrote seconds through v1.6.2 and minutes from v1.7.0
    # without changing the file type, so no honest answer exists.
    assert result.retention_time_unit == "unknown"


def test_converting_an_unknown_retention_time_raises_rather_than_guessing(recorded):
    recorded("readers_features_topfd.json")
    result = readers.read_features("sample_ms1.feature")

    with pytest.raises(_bridge.UsageError, match="no basis to say"):
        _ = result.retention_time_start_in_minutes


def test_the_expansion_caveat_is_present_for_ms1_feature(recorded):
    recorded("readers_features_topfd.json")

    result = readers.read_features("sample_ms1.feature")

    assert any("CHARGE STATE" in caveat for caveat in result.caveats), (
        "One row of this view is one charge state of one feature, not one line of the file. "
        "A caller comparing record_count to the file's line count must be told."
    )


# ---- read_matches --------------------------------------------------------------------------


def test_read_matches_parses_the_spectral_match_view(recorded):
    recorded("readers_matches_mspathfinder.json")

    result = readers.read_matches("results_IcTda.tsv")

    assert isinstance(result, readers.MatchRecords)
    assert "modifications" in result.column_names
    assert "accession" in result.column_names


def test_casanovo_is_decoy_is_none_not_false(recorded):
    recorded("readers_matches_casanovo.json")

    result = readers.read_matches("run.mztab")

    # De novo sequencing has no target/decoy label at all. False would let a caller filter on a
    # fabricated column - the same trap read_results already refuses for MSFragger.
    assert set(result.columns["is_decoy"]) == {None}
    assert any("de novo" in caveat for caveat in result.caveats)


def test_mspathfindert_reports_real_decoy_flags(recorded):
    recorded("readers_matches_mspathfinder.json")

    result = readers.read_matches("results_IcTda.tsv")

    assert all(value in (True, False) for value in result.columns["is_decoy"])


def test_read_matches_warns_that_nothing_is_fdr_filtered(recorded):
    recorded("readers_matches_casanovo.json")

    result = readers.read_matches("run.mztab")

    assert any("FDR" in caveat for caveat in result.caveats)


# ---- read_spectra --------------------------------------------------------------------------


def test_read_spectra_parses_scan_headers(recorded):
    recorded("readers_spectra_mzml.json")

    result = readers.read_spectra("run.mzML")

    assert isinstance(result, readers.ScanRecords)
    assert result.retention_time_unit == "minutes"
    assert "one_based_scan_number" in result.column_names
    assert result.peaks_included is False


def test_read_spectra_reports_the_files_total_alongside_the_filtered_count(recorded):
    # Deliberately the FILTERED payload. On an unfiltered one the two counts are equal by
    # construction, so `scan_count >= record_count` would pass for an implementation that simply
    # set scan_count = record_count - which is precisely the failure this field exists to prevent.
    recorded("readers_spectra_ms2.json")

    result = readers.read_spectra("run.mzML", ms_order=2)

    assert result.ms_order == 2
    assert result.scan_count > result.record_count, (
        "the file's real total must survive the filter, so a filter that matched nothing can "
        "never look like an empty file"
    )


def test_read_spectra_omits_peak_columns_by_default(recorded):
    recorded("readers_spectra_mzml.json")

    result = readers.read_spectra("run.mzML")

    assert "mz" not in result.columns
    assert "peak_count" in result.columns


def test_read_spectra_sends_ms_order_and_peaks(captured):
    readers.read_spectra("run.mzML", ms_order=2, peaks=True, limit=10)

    assert captured["args"] == (
        "readers", "read-spectra", "--path", "run.mzML",
        "--limit", "10", "--ms-order", "2", "--peaks",
    )


def test_read_spectra_omits_the_peaks_flag_when_false(captured):
    readers.read_spectra("run.mzML", peaks=False)

    assert "--peaks" not in captured["args"]


@pytest.mark.parametrize("bad", [0, -1, "2", 1.5, True])
def test_read_spectra_rejects_a_bad_ms_order_before_spawning_anything(bad, captured):
    with pytest.raises(_bridge.UsageError, match="ms_order"):
        readers.read_spectra("run.mzML", ms_order=bad)

    assert "args" not in captured, "Validation must happen before the bridge is invoked."


# ---- shared argument validation --------------------------------------------------------------

READ_VERBS = [
    ("read_records", "read-records"),
    ("read_features", "read-features"),
    ("read_matches", "read-matches"),
    ("read_spectra", "read-spectra"),
    ("read_results", "read-results"),
]


@pytest.mark.parametrize(("name", "verb"), READ_VERBS)
def test_every_read_verb_sends_its_own_verb(name, verb, captured):
    getattr(readers, name)("a.tsv")

    assert captured["args"][:2] == ("readers", verb)


@pytest.mark.parametrize(("name", "_verb"), READ_VERBS)
@pytest.mark.parametrize("bad", ["", "   ", None, 7])
def test_every_read_verb_rejects_a_blank_path(name, _verb, bad, captured):
    with pytest.raises(_bridge.UsageError, match="file path is required"):
        getattr(readers, name)(bad)

    assert "args" not in captured


@pytest.mark.parametrize(("name", "_verb"), READ_VERBS)
@pytest.mark.parametrize("bad", [-1, "5", 2.5, True])
def test_every_read_verb_rejects_a_bad_limit(name, _verb, bad, captured):
    with pytest.raises(_bridge.UsageError, match="limit"):
        getattr(readers, name)("a.tsv", limit=bad)

    assert "args" not in captured


@pytest.mark.parametrize(("name", "_verb"), READ_VERBS)
@pytest.mark.parametrize("bad", [-1, "5", 2.5, True])
def test_every_read_verb_rejects_a_bad_offset(name, _verb, bad, captured):
    with pytest.raises(_bridge.UsageError, match="offset"):
        getattr(readers, name)("a.tsv", offset=bad)

    assert "args" not in captured


@pytest.mark.parametrize(("name", "_verb"), READ_VERBS)
def test_every_read_verb_rejects_a_blank_out(name, _verb, captured):
    with pytest.raises(_bridge.UsageError, match="out must be"):
        getattr(readers, name)("a.tsv", out="  ")

    assert "args" not in captured


@pytest.mark.parametrize(("name", "_verb"), READ_VERBS)
def test_every_read_verb_omits_a_zero_offset(name, _verb, captured):
    getattr(readers, name)("a.tsv", offset=0)

    # A default that is sent explicitly is a default the bridge can disagree with later.
    assert "--offset" not in captured["args"]


# ---- the public surface ----------------------------------------------------------------------


def test_the_module_exports_every_read_verb():
    for name, _ in READ_VERBS:
        assert name in readers.__all__

    # Deliberately NOT `set(__all__) == set(dir(readers))`: `__dir__` is defined as
    # `sorted(__all__)`, so that comparison is a tautology that would pass for an `__all__` naming
    # symbols the module does not have.
    assert all(hasattr(readers, name) for name in readers.__all__), (
        "every exported name must resolve"
    )
    # Only names DEFINED here - `vars()` also holds this module's own imports (`Any`, `dataclass`),
    # which __all__ deliberately hides and which are not part of the public surface.
    defined_here = {
        name
        for name, value in vars(readers).items()
        if not name.startswith("_") and getattr(value, "__module__", None) == readers.__name__
    }
    assert defined_here - set(readers.__all__) == set(), (
        "a name defined here but missing from __all__ is invisible to dir() and help()"
    )


def test_the_view_names_are_exported_as_constants():
    assert readers.MS1_FEATURES == "ms1_features"
    assert readers.SPECTRAL_MATCH == "spectral_match"
    assert readers.SPECTRA == "spectra"


@pytest.mark.parametrize(
    ("name", "fixture"),
    [
        ("read_records", "readers_records_toppic.json"),
        ("read_features", "readers_features_topfd.json"),
        ("read_matches", "readers_matches_casanovo.json"),
        ("read_spectra", "readers_spectra_mzml.json"),
    ],
)
def test_converting_a_written_table_explains_itself(name, fixture, monkeypatch):
    data = payload(fixture)
    data = dict(data, columns=None, output={"path": "out.tsv", "format": "tsv", "row_count": 3})
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: data)

    result = getattr(readers, name)("a.tsv", out="out.tsv")

    assert result.records == [], "There is nothing to iterate when the table went to disk."
    assert result.output.path == "out.tsv"
