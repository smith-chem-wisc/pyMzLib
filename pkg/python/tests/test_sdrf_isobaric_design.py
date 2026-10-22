"""Tests for :func:`pymzlib.sdrf.isobaric_design`, the projection of mzLib's ``SdrfIsobaricDesign``.

Offline except for the live tests at the end. The payloads are recorded from the real bridge over
mzLib's own isobaric design fixtures, copied beside them: ``PXD008841.subset.sdrf.tsv`` (TMT10, two
plexes named only in the file names), ``PXD061609.subset.sdrf.tsv`` (TMTpro 18, three SAX arms of one
plex in ``comment[sample preparation batch]``) and ``PXD008841_searched.sdrf.tsv`` (PXD008841
restricted to a calibrated search by ``sdrf restrict --out``).

What the Python layer must get exactly right: every option is sent only when given, the condition
columns cross as ONE tab-joined option and searched files on stdin, a refusal stays a result, and
the numbers stay 1-based (nothing converts them).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, sdrf

FIXTURES = Path(__file__).parent / "fixtures"
PXD008841 = FIXTURES / "PXD008841.subset.sdrf.tsv"
PXD061609 = FIXTURES / "PXD061609.subset.sdrf.tsv"
SEARCHED = FIXTURES / "PXD008841_searched.sdrf.tsv"
POOL = r"TMT_?pool(?<plex>\d+)"
BATCH = "comment[sample preparation batch]"


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


def test_a_valid_design_projects_every_annotated_channel(replay):
    seen = replay("sdrf_isobaric_design_PXD008841.json")
    d = sdrf.isobaric_design(PXD008841, kit="TMT10", plex_file_name_pattern=POOL)

    assert seen["args"] == (
        "sdrf", "isobaric-design", "--path", str(PXD008841),
        "--kit", "TMT10", "--plex-file-name-pattern", POOL,
    )
    assert seen["kwargs"]["stdin"] is None
    assert d.is_valid and d.refusals == [] and d.notes == []
    assert d.kit == "TMT10" and d.kit_type == "TMT10" and len(d.kit_channels) == 10
    assert d.plex_file_name_pattern == POOL and d.plex_column is None and d.single_plex is None
    assert d.plex_ids == {"1": 1, "5": 2} and d.plex_count == 2
    assert d.file_count == 4 and d.record_count == 40 and len(d.channels) == 40
    assert d.file_key_column == "comment[data file]"
    assert d.condition_columns == ["factor value[disease]"] and not d.condition_columns_declared
    assert not d.searched_files_given and d.written is None
    assert d.column_names == [
        "full_path", "plex", "plex_id", "sample_name", "channel", "condition",
        "biological_replicate", "fraction", "technical_replicate",
    ]
    assert len(d.caveats) == 8 and d.caveats[0].startswith("1-BASED")
    assert "4 file(s), 2 plex(es), kit TMT10" in d.report


def test_the_numbers_stay_one_based(replay):
    replay("sdrf_isobaric_design_PXD008841.json")
    d = sdrf.isobaric_design(PXD008841, kit="TMT10", plex_file_name_pattern=POOL)

    assert {c.fraction for c in d.channels} == {1, 2}
    assert {c.technical_replicate for c in d.channels} == {1}
    assert {c.biological_replicate for c in d.channels} == {1}
    assert {c.plex_id for c in d.channels} == {1, 2}
    first_file = [c for c in d.channels if c.full_path == d.channels[0].full_path]
    assert [c.channel for c in first_file] == d.kit_channels, "reporter m/z order"
    assert first_file[-1].sample_name == "pool"


def test_a_refusal_is_a_result(replay):
    seen = replay("sdrf_isobaric_design_PXD008841_noplex.json")
    d = sdrf.isobaric_design(PXD008841, kit="TMT10")

    assert seen["args"] == ("sdrf", "isobaric-design", "--path", str(PXD008841), "--kit", "TMT10")
    assert not d.is_valid
    assert d.refusals[0].startswith("No plex source was declared")
    assert d.plex_source == "(not declared)"
    assert d.file_count == 0 and d.record_count == 0 and d.channels == [] and d.plex_ids == {}


def test_a_declared_batch_column_gives_three_plexes(replay):
    seen = replay("sdrf_isobaric_design_PXD061609.json")
    d = sdrf.isobaric_design(PXD061609, kit="TMT18", plex_column=BATCH)

    assert seen["args"][-4:] == ("--kit", "TMT18", "--plex-column", BATCH)
    assert d.is_valid and d.plex_ids == {"HighSalt": 1, "LowSalt": 2, "NoSAX": 3}
    assert d.record_count == 108 and len(d.kit_channels) == 18


def test_one_plex_over_three_arms_is_refused_every_reason_at_once(replay):
    replay("sdrf_isobaric_design_PXD061609_oneplex.json")
    d = sdrf.isobaric_design(PXD061609, kit="TMT18", single_plex="tube1")

    assert not d.is_valid and d.single_plex == "tube1"
    assert len(d.refusals) == 36
    assert all("more than once" in r for r in d.refusals)


def test_condition_columns_and_searched_files_cross_as_design_sends_them(replay, tmp_path):
    seen = replay("sdrf_isobaric_design_searched.json")
    runs = [tmp_path / "a-calib.mzML", "b-calib.mzML"]
    sdrf.isobaric_design(
        SEARCHED, kit="TMT10", plex_file_name_pattern=POOL,
        condition_columns=["factor value[disease]"], searched_files=runs, out=tmp_path / "design",
    )

    assert seen["args"][-5:] == (
        "--condition-columns", "factor value[disease]", "--searched-files-stdin",
        "--out", str(tmp_path / "design"),
    )
    assert seen["kwargs"]["stdin"] == f"{os.fspath(runs[0])}\nb-calib.mzML\n"


def test_a_restricted_document_keys_on_the_searched_name(replay):
    replay("sdrf_isobaric_design_searched.json")
    d = sdrf.isobaric_design(SEARCHED, kit="TMT10", plex_file_name_pattern=POOL, searched_files=["x"], out="design")

    assert d.file_key_column == "comment[searched data file]"
    assert d.searched_files_given
    assert {c.full_path for c in d.channels} == {
        rf"D:\search\HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0{i}-calib.mzML" for i in (1, 2)
    } | {rf"D:\search\HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0{i}-calib.mzML" for i in (1, 2)}
    assert d.written["record_count"] == 40 and d.written["file_count"] == 4
    assert d.written["path"].endswith("TmtDesign.txt")


def test_no_option_is_sent_unless_given(replay):
    seen = replay("sdrf_isobaric_design_PXD008841_noplex.json")
    sdrf.isobaric_design(PXD008841)
    assert seen["args"] == ("sdrf", "isobaric-design", "--path", str(PXD008841))
    assert seen["kwargs"]["stdin"] is None


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"kit": ""}, "kit must be a non-blank str"),
        ({"kit": 10}, "kit must be a non-blank str"),
        ({"plex_column": "  "}, "plex_column must be a non-blank str"),
        ({"plex_file_name_pattern": "a\nb"}, "newline"),
        ({"single_plex": ""}, "single_plex must be a non-blank str"),
        ({"condition_columns": "factor value[disease]"}, "list of column names"),
        ({"condition_columns": []}, "is empty"),
        ({"searched_files": "run.raw"}, "list of paths"),
        ({"searched_files": [""]}, "non-blank"),
        ({"out": ""}, "out must name a folder"),
    ],
)
def test_bad_arguments_are_refused_before_the_bridge(monkeypatch, kwargs, message):
    def no_call(*args, **kw):
        raise AssertionError("the bridge must not be called")

    monkeypatch.setattr(_bridge, "invoke", no_call)
    with pytest.raises(_bridge.UsageError, match=message):
        sdrf.isobaric_design(PXD008841, **kwargs)


def test_a_blank_path_is_refused(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("called"))
    with pytest.raises(_bridge.UsageError):
        sdrf.isobaric_design("", kit="TMT10")


def test_an_older_bridge_is_named_not_spawned(monkeypatch):
    calls = []

    def old_bridge(*args, **kwargs):
        calls.append(args)
        return {"verbs": ["sdrf design"]}

    monkeypatch.setattr(_bridge, "invoke", old_bridge)
    monkeypatch.setattr(_bridge, "_VERBS_SEEN", {})
    with pytest.raises(_bridge.UsageError, match="pyMzLib 0.6.0"):
        sdrf.isobaric_design(PXD008841, kit="TMT10")
    assert calls == [("version",)]


# ---- live: the real bridge ----------------------------------------------------------------------


@pytest.fixture()
def built_bridge():
    # As test_sdrf_design: skip in a source checkout with no built bridge, fail under CI, where the
    # bridge is always built and a skip would hide drift between the recording and mzLib.
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so sdrf isobaric-design went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")
    _bridge._VERBS_SEEN.clear()


def test_live_design_matches_the_recording_and_writes_tmt_design(built_bridge, tmp_path):
    folder = tmp_path / "new" / "design"
    d = sdrf.isobaric_design(PXD008841, kit="TMT10", plex_file_name_pattern=POOL, out=folder)
    recorded = _load("sdrf_isobaric_design_PXD008841.json")

    assert d.columns == recorded["columns"]
    assert d.plex_ids == recorded["plex_ids"]
    assert d.written == {"path": str(folder / "TmtDesign.txt"), "file_count": 4, "record_count": 40}
    lines = (folder / "TmtDesign.txt").read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "File", "Plex", "Sample Name", "TMT Channel", "Condition",
        "Biological Replicate", "Fraction", "Technical Replicate",
    ]
    first = lines[1].split("\t")
    c = d.channels[0]
    assert first == [c.full_path, c.plex, c.sample_name, c.channel, c.condition,
                     str(c.biological_replicate), str(c.fraction), str(c.technical_replicate)], (
        "the file and the wire carry the same 1-based numbers"
    )


def test_live_unknown_kit_is_the_isobaric_kits_usage_error(built_bridge):
    from pymzlib import isobaric

    with pytest.raises(_bridge.UsageError) as from_design:
        sdrf.isobaric_design(PXD008841, kit="TMT10plex", plex_file_name_pattern=POOL)
    with pytest.raises(_bridge.UsageError) as from_kits:
        isobaric.kits("TMT10plex")
    assert str(from_design.value) == str(from_kits.value)
