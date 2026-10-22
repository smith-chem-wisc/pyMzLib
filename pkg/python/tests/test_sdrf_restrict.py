"""Tests for :func:`pymzlib.sdrf.restrict`, the projection of mzLib's ``SdrfSearchScope.Restrict``.

Offline except for the live tests at the end. The payloads are recorded from the real bridge over
mzLib's ``PXD008841.subset.sdrf.tsv`` (40 rows, 4 files, TMT10): one for a search that exercises
the reports (a calibrated derivative, the same run's raw file as well, a renamed file whose
acquisition is named, a file the SDRF lacks), one for a clean calibrated search written with
``out=``.

What the Python layer must get exactly right: the searched files become stdin lines, with
``searched<TAB>acquired`` only where an acquired file is named; ``out`` and the window are passed
through; and the restricted document is an :class:`SdrfDocument`, so its accessors work on it.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, sdrf

FIXTURES = Path(__file__).parent / "fixtures"
PXD008841 = FIXTURES / "PXD008841.subset.sdrf.tsv"
POOL1 = "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0"
POOL5 = "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0"
MESSY = {
    rf"D:\search\{POOL1}1-calib.mzML": None,
    f"{POOL1}2-calib.mzML": None,
    f"{POOL1}2.raw": None,
    "sample_A.mzML": f"{POOL5}1.raw",
    "not_in_sdrf.mzML": None,
}
CLEAN = [rf"D:\search\{p}{i}-calib.mzML" for p in (POOL1, POOL5) for i in (1, 2)]


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, **kwargs):
            if args == ("version",):
                return _load("bridge_version.json")
            seen["args"] = args
            seen["kwargs"] = kwargs
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        monkeypatch.setattr(_bridge, "_VERBS_SEEN", {})
        return seen

    return use


def test_a_mapping_sends_an_acquired_name_only_where_one_is_given(replay):
    seen = replay("sdrf_restrict_PXD008841.json")
    scoped = sdrf.restrict(PXD008841, MESSY, limit=2)

    assert seen["args"] == ("sdrf", "restrict", "--path", str(PXD008841), "--limit", "2")
    assert seen["kwargs"]["stdin"] == (
        rf"D:\search\{POOL1}1-calib.mzML" "\n"
        f"{POOL1}2-calib.mzML\n"
        f"{POOL1}2.raw\n"
        f"sample_A.mzML\t{POOL5}1.raw\n"
        "not_in_sdrf.mzML\n"
    )
    assert scoped.searched_files == list(MESSY)
    assert scoped.acquired_name_count == 1


def test_the_reports_are_mzlibs(replay):
    replay("sdrf_restrict_PXD008841.json")
    scoped = sdrf.restrict(PXD008841, MESSY, limit=2)

    assert scoped.input_row_count == 40 and scoped.row_count == 30
    assert scoped.returned_count == 2 and scoped.truncated and scoped.offset == 0
    assert scoped.searched_without_row == ["not_in_sdrf.mzML"]
    assert scoped.dropped_data_files == [f"{POOL5}2.raw"]
    assert scoped.ambiguous == [f"{POOL1}2: {POOL1}2-calib.mzML, {POOL1}2.raw"]
    assert scoped.shared_stem == []
    assert scoped.written is None
    assert len(scoped.caveats) == 6


def test_the_restricted_document_is_an_sdrf_document(replay):
    replay("sdrf_restrict_PXD008841.json")
    scoped = sdrf.restrict(PXD008841, MESSY, limit=2)

    assert isinstance(scoped, sdrf.SdrfDocument)
    data_file = scoped.index_of("comment[data file]")
    assert scoped.columns[data_file + 1] == sdrf.RestrictedSdrf.SEARCHED_DATA_FILE_COLUMN
    assert scoped.value("comment[data file]") == [f"{POOL1}1.raw"] * 2
    assert scoped.searched_data_files() == [f"{POOL1}1-calib.mzML"] * 2, "the name, no folder"
    assert all(len(row) == len(scoped.columns) for row in scoped.rows)


def test_a_list_sends_one_name_per_line_and_out_is_passed_through(replay, tmp_path):
    seen = replay("sdrf_restrict_PXD008841_out.json")
    target = tmp_path / "PXD008841_searched.sdrf.tsv"
    scoped = sdrf.restrict(PXD008841, [Path(CLEAN[0]), *CLEAN[1:]], out=target, limit=2)

    assert seen["args"][-4:] == ("--limit", "2", "--out", str(target))
    assert seen["kwargs"]["stdin"] == "\n".join(CLEAN) + "\n", "a PathLike crosses as its text"
    assert scoped.written == sdrf.WrittenSdrf(path="PXD008841_searched.sdrf.tsv", row_count=40)
    assert scoped.row_count == 40 and scoped.returned_count == 2


@pytest.mark.parametrize(
    "searched, message",
    [
        ("X-calib.mzML", "list of paths"),
        (Path("X-calib.mzML"), "list of paths"),
        ([], "at least one searched file"),
        ({}, "at least one searched file"),
        ([""], "non-blank"),
        (["a\tb"], "tab or newline"),
        ({"a.mzML": ""}, "acquired file for 'a.mzML'"),
        ({"a.mzML": "b\n.raw"}, "tab or newline"),
        (42, "list of paths or"),
    ],
)
def test_bad_searched_files_are_refused_before_the_bridge(monkeypatch, searched, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge must not be called"))
    with pytest.raises(_bridge.UsageError, match=message):
        sdrf.restrict(PXD008841, searched)


@pytest.mark.parametrize(
    "kwargs, message",
    [({"out": ""}, "ending in .tsv"), ({"limit": -1}, "limit"), ({"offset": -1}, "offset")],
)
def test_bad_options_are_refused_before_the_bridge(monkeypatch, kwargs, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge must not be called"))
    with pytest.raises(_bridge.UsageError, match=message):
        sdrf.restrict(PXD008841, CLEAN, **kwargs)


def test_a_blank_path_is_refused(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("called"))
    with pytest.raises(_bridge.UsageError):
        sdrf.restrict("", CLEAN)


# ---- live: the real bridge ----------------------------------------------------------------------


@pytest.fixture()
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so sdrf restrict went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")
    _bridge._VERBS_SEEN.clear()


def test_live_restriction_matches_the_recording(built_bridge):
    scoped = sdrf.restrict(PXD008841, MESSY, limit=2)
    recorded = _load("sdrf_restrict_PXD008841.json")

    assert scoped.rows == recorded["rows"]
    assert scoped.columns == recorded["column_names"]
    for name in ("searched_without_row", "dropped_data_files", "ambiguous", "shared_stem"):
        assert getattr(scoped, name) == recorded[name]


def test_live_out_writes_the_whole_document_where_mzlib_says(built_bridge, tmp_path):
    scoped = sdrf.restrict(PXD008841, CLEAN, out=tmp_path / "new" / "searched.tsv", limit=1)

    written = Path(scoped.written.path)
    assert written.name == "searched.tsv.sdrf.tsv", "mzLib adds .sdrf.tsv to a name not ending in it"
    assert written.exists() and not (tmp_path / "new" / "searched.tsv").exists()
    back = sdrf.read(written)
    assert back.row_count == 40 and scoped.returned_count == 1
    assert back.rows == sdrf.read(FIXTURES / "PXD008841_searched.sdrf.tsv").rows

    with pytest.raises(_bridge.UsageError, match=r"\.tsv"):
        sdrf.restrict(PXD008841, CLEAN, out=tmp_path / "searched.txt")


def test_live_restrict_then_design(built_bridge):
    d = sdrf.isobaric_design(
        FIXTURES / "PXD008841_searched.sdrf.tsv", kit="TMT10",
        plex_file_name_pattern=r"TMT_?pool(?<plex>\d+)", searched_files=CLEAN,
    )
    assert d.is_valid, d.report
    assert d.file_key_column == "comment[searched data file]"
    assert sorted({c.full_path for c in d.channels}) == sorted(CLEAN)

    unrestricted = sdrf.isobaric_design(
        PXD008841, kit="TMT10", plex_file_name_pattern=r"TMT_?pool(?<plex>\d+)", searched_files=CLEAN,
    )
    assert not unrestricted.is_valid
    assert any("has no SDRF row" in r for r in unrestricted.refusals)


def test_live_a_document_without_a_data_file_column_is_a_bridge_error(built_bridge, tmp_path):
    path = tmp_path / "bare.sdrf.tsv"
    path.write_text("source name\tassay name\nS1\trun 1\n", encoding="utf-8")
    with pytest.raises(_bridge.BridgeError) as caught:
        sdrf.restrict(path, ["X.raw"])
    assert caught.value.error_type == "ArgumentException"
    assert "comment[data file]" in str(caught.value)
