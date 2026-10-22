"""Tests for :mod:`pymzlib.library`, the projection of mzLib's ``.msl`` spectral-library reader.

Offline. Every fixture is recorded from the real bridge (scripts/record_fixture.py), from the two
libraries committed under ``fixtures/library/`` (written by the bridge tests' ``LibraryFixtures``):

* ``myPrositLib.msl``: mzLib's own ``myPrositLib.msp`` test library, five Prosit-predicted
  precursors, converted by mzLib. ``library_read_prosit.json``.
* ``synthetic_library.msl``: six precursors covering every field - target, decoy, entrapment
  target with calibrated minutes, entrapment decoy, a predicted precursor with every "not available"
  sentinel, a proteoform. ``library_read_synthetic.json``, ``library_read_fragments.json`` (one
  precursor's fragment ions), ``library_read_index_only.json``.
* ``library_read_many.json``: both, through ``--paths-stdin``.

What is pinned here is the call pyMzLib sends and that it projects mzLib's answer without
repairing it: a null stays ``None``, never a sentinel, and nothing is converted.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, library

FIXTURES = Path(__file__).parent / "fixtures"


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, stdin=None, timeout=None):
            seen["args"] = args
            seen["stdin"] = stdin
            seen["timeout"] = timeout
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: seen.setdefault("guard", (verb, since)))
        return seen

    return use


# ---- the call ---------------------------------------------------------------------------------


def test_a_plain_read_sends_only_the_path(replay):
    seen = replay("library_read_prosit.json")
    library.read("myPrositLib.msl")

    assert seen["args"] == ("library", "read", "--path", "myPrositLib.msl")
    assert seen["guard"] == ("library read", "0.6.0")
    assert seen["stdin"] is None


def test_every_option_is_sent_under_its_wire_name(replay):
    seen = replay("library_read_fragments.json")
    library.read(
        Path("synthetic_library.msl"), fragments=True, index_only_threshold=0, limit=1, offset=2,
        out="ions.tsv", timeout=30,
    )

    assert seen["args"] == (
        "library", "read", "--path", "synthetic_library.msl", "--limit", "1", "--offset", "2",
        "--out", "ions.tsv", "--fragments", "--index-only-threshold", "0",
    )
    assert seen["timeout"] == 30


def test_the_threshold_is_left_to_mzlib_when_not_given(replay):
    seen = replay("library_read_synthetic.json")
    result = library.read("synthetic_library.msl")

    assert "--index-only-threshold" not in seen["args"]
    # The echo is mzLib's own default, MslFileTypeHandler.DefaultIndexOnlyThresholdBytes.
    assert result.index_only_threshold == 1073741824


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"fragments": "yes"}, "fragments must be True or False"),
        ({"index_only_threshold": -1}, "index_only_threshold"),
        ({"index_only_threshold": True}, "index_only_threshold"),
        ({"index_only_threshold": 1.5e9}, "index_only_threshold"),
        ({"limit": -1}, "limit"),
        ({"offset": True}, "offset"),
    ],
)
def test_bad_arguments_are_refused_before_the_bridge_is_spawned(monkeypatch, kwargs, message):
    # Not even the version probe: require_verb runs only once the arguments are good.
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: pytest.fail("the bridge was probed"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))

    with pytest.raises(pymzlib.UsageError, match=message):
        library.read("synthetic_library.msl", **kwargs)


def test_a_list_given_to_read_points_at_read_many(monkeypatch):
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: pytest.fail("the bridge was probed"))
    with pytest.raises(pymzlib.UsageError, match=r"read_many\(\)"):
        library.read(["a.msl", "b.msl"])


def test_read_many_sends_the_list_on_stdin(replay):
    seen = replay("library_read_many.json")
    library.read_many(["synthetic_library.msl", Path("myPrositLib.msl")], threads=2, on_error="skip", fragments=True)

    assert seen["args"] == (
        "library", "read", "--paths-stdin", "--threads", "2", "--on-error", "skip", "--fragments",
    )
    assert seen["stdin"] == "synthetic_library.msl\nmyPrositLib.msl\n"


@pytest.mark.parametrize(
    "paths, kwargs, message",
    [
        ("one.msl", {}, "LIST of paths"),
        ([], {}, "empty"),
        (["a.msl", ""], {}, "Path 1"),
        (["a.msl"], {"threads": 0}, "threads"),
        (["a.msl"], {"on_error": "ignore"}, "on_error"),
        (["a.msl"], {"index_only_threshold": "1GB"}, "index_only_threshold"),
    ],
)
def test_read_many_refuses_bad_arguments_before_the_bridge(monkeypatch, paths, kwargs, message):
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: pytest.fail("the bridge was probed"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))

    with pytest.raises(pymzlib.UsageError, match=message):
        library.read_many(paths, **kwargs)


# ---- the projection ---------------------------------------------------------------------------


def test_the_envelope_is_projected_field_for_field(replay):
    replay("library_read_synthetic.json")
    result = library.read("synthetic_library.msl")
    wire = _load("library_read_synthetic.json")

    assert result.file_type is None
    assert result.reader == "MslLibrary"
    assert (result.format_version, result.is_compressed, result.is_index_only) == (3, False, False)
    assert (result.target_count, result.decoy_count, result.proteoform_count) == (4, 2, 1)
    assert result.record_count == result.returned_count == result.row_count == 6
    assert result.min_precursor_mz == wire["min_precursor_mz"]
    assert result.max_precursor_mz == wire["max_precursor_mz"]
    assert result.rows_not_read is None and result.retention_time_unit is None
    assert result.fragments_included is False
    assert result.truncated is False and result.offset == 0
    assert result.output is None
    assert result.caveats == wire["caveats"]
    assert result.column_names == wire["column_names"]
    assert list(result.columns) == result.column_names


def test_rows_are_the_stored_order_not_mz_order(replay):
    replay("library_read_synthetic.json")
    columns = library.read("synthetic_library.msl").columns

    assert columns["precursor_index"] == [0, 1, 2, 3, 4, 5]
    assert columns["precursor_mz"] != sorted(columns["precursor_mz"])


def test_labels_and_retention_time_units_are_per_row(replay):
    replay("library_read_synthetic.json")
    columns = library.read("synthetic_library.msl").columns

    assert columns["decoy_contaminant_target"] == ["T", "D", "ET", "ED", "T", "T"]
    assert columns["is_decoy"] == [False, True, False, True, False, False]
    assert columns["is_entrapment"] == [False, False, True, True, False, False]
    # Precursor 2 is calibrated minutes; every other is iRT. Nothing is converted.
    assert columns["retention_time_unit"] == ["iRT", "iRT", "minutes", "iRT", "iRT", "iRT"]
    assert columns["retention_time"][2] == pytest.approx(21.4, abs=1e-4)


def test_mzlibs_not_available_sentinels_arrive_as_none(replay):
    replay("library_read_synthetic.json")
    columns = library.read("synthetic_library.msl").columns

    predicted = 4
    for name in ("ion_mobility", "q_value", "nce", "protein_accession", "protein_name", "gene_name"):
        assert columns[name][predicted] is None, name
    assert columns["dissociation_type"][predicted] == "Unknown"
    assert columns["ion_mobility"][0] == pytest.approx(0.95, abs=1e-6)
    assert columns["q_value"][0] == pytest.approx(0.001)
    assert columns["nce"][0] == 28


def test_a_library_with_no_proteins_names_those_fields_absent(replay):
    replay("library_read_prosit.json")
    result = library.read("myPrositLib.msl")

    assert result.absent_fields == ["ion_mobility", "protein_accession", "protein_name", "gene_name"]
    assert result.columns["protein_accession"] == [None] * 5
    # Converted from a text library: mzLib's not-available defaults for what .msp cannot carry.
    assert result.columns["nce"] == [None] * 5
    assert result.columns["source_type"] == ["Predicted"] * 5


def test_fragment_rows_repeat_the_precursor_and_keep_mzlibs_values(replay):
    replay("library_read_fragments.json")
    result = library.read("synthetic_library.msl", fragments=True, limit=1)
    rows = result.records

    assert result.fragments_included is True
    assert (result.returned_count, result.row_count, result.truncated) == (1, 6, True)
    assert {r["full_sequence"] for r in rows} == {"PEPTIDEK"}
    assert [r["annotation"] for r in rows] == ["y1", "b2", "y2", "bIy[2-5]", "y3-18.0106", "y4-43.0058"]
    custom = rows[-1]
    assert custom["neutral_loss"] == pytest.approx(-43.005814)
    assert custom["residue_position"] is None       # mzLib reads 0 back for a custom loss
    internal = rows[3]
    assert (internal["secondary_product_type"], internal["secondary_fragment_number"]) == ("y", 5)
    assert all(r["secondary_fragment_number"] is None for r in rows if not r["is_internal"])
    assert rows[0]["neutral_loss"] == 0                # no loss is 0, not None
    assert max(r["fragment_intensity"] for r in rows) == 1
    assert {r["fragment_count"] for r in rows} == {6}


def test_an_index_only_read_says_so_and_gives_the_same_rows(replay):
    replay("library_read_index_only.json")
    index_only = library.read("synthetic_library.msl", index_only_threshold=0, limit=2)
    full = _load("library_read_synthetic.json")["columns"]

    assert index_only.is_index_only is True
    assert index_only.index_only_threshold == 0
    assert index_only.columns == {name: values[:2] for name, values in full.items()}


def test_a_batch_projects_one_report_per_library(replay):
    replay("library_read_many.json")
    batch = library.read_many(["synthetic_library.msl", "myPrositLib.msl"])

    assert (batch.file_count, batch.read_count, batch.failed_count) == (2, 2, 0)
    assert batch.record_count == batch.returned_count == batch.row_count == 11
    assert batch.column_names[:2] == ["source_index", "source_path"]
    assert batch.columns["source_index"] == [0] * 6 + [1] * 5
    assert [f.record_count for f in batch.files] == [6, 5]
    assert [f.reader for f in batch.files] == ["MslLibrary", "MslLibrary"]
    assert batch.files[0].proteoform_count == 1
    assert batch.files[1].absent_fields == ["ion_mobility", "protein_accession", "protein_name", "gene_name"]
    assert all(f.ok for f in batch.files) and batch.failed_files == []
    assert batch.records[6]["full_sequence"] == "ALAVDGAGKPGAEE"


def test_a_skipped_library_keeps_its_unknown_facts_as_none():
    report = library.LibraryFileReport._from_wire({
        "path": "C:/missing.msl", "file_type": None, "reader": None, "format_version": None,
        "is_compressed": None, "is_index_only": None, "target_count": None, "decoy_count": None,
        "proteoform_count": None, "min_precursor_mz": None, "max_precursor_mz": None,
        "record_count": None, "rows_not_read": None, "retention_time_unit": None, "caveats": [],
        "column_names": [], "absent_fields": [], "failed_fields": [], "excluded_fields": [],
        "error": {"kind": "usage", "type": "usage", "message": "File not found: 'C:/missing.msl'."},
    })

    assert not report.ok
    assert report.error.kind == "usage"
    assert (report.record_count, report.target_count, report.is_index_only) == (None, None, None)


def test_a_table_written_to_out_returns_no_columns():
    written = library.LibraryRecords._from_wire({
        **_load("library_read_synthetic.json"),
        "returned_count": 0, "row_count": 0, "columns": None,
        "output": {"path": "C:/out/library.tsv", "format": "tsv", "row_count": 6},
    })

    assert written.columns is None and written.records == []
    assert (written.output.path, written.output.row_count) == ("C:/out/library.tsv", 6)


def test_the_public_api_is_what_dir_lists():
    assert dir(library) == sorted(library.__all__)
    assert "library" in pymzlib.__all__
