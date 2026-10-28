"""pymzlib.library.calibrate_rt against recorded bridge payloads.

The numbers are mzLib's, held to IrtCalibration.Fit by the bridge's C# tests
(LibraryCalibrateRtTests) and re-run live by test_library_calibrate_live.py. These tests hold the
Python layer: that it sends what the wire expects, parses what the wire returns, and refuses bad
input before starting a process.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from pymzlib import _bridge, library

FIXTURES = Path(__file__).parent / "fixtures"


def recorded(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture
def bridge(monkeypatch):
    """Answer every bridge call with ``state["payload"]``, recording the arguments and stdin."""
    state = {"calls": [], "payload": recorded("library_calibrate_rt_pbmc_linear.json")}

    def invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": list(args), "stdin": stdin, "timeout": timeout})
        return state["payload"]

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    return state


# ---- what is sent --------------------------------------------------------------------------------


def test_anchors_then_a_separator_then_queries_with_the_computed_time_blank(bridge):
    library.calibrate_rt(
        [("run1", 10, 20.5), {"group": "run1", "rt_minutes": 11.0, "irt": 22}],
        [("run1", None, 30), {"group": "run2", "rt_minutes": 12.5}, ("run1", float("nan"), 1.0)],
    )
    call = bridge["calls"][0]
    assert call["args"] == ["library", "calibrate-rt"]
    assert call["stdin"] == (
        "run1\t10.0\t20.5\nrun1\t11.0\t22.0\n--\nrun1\t\t30.0\nrun2\t12.5\t\nrun1\t\t1.0\n"
    )


def test_no_queries_sends_no_separator(bridge):
    library.calibrate_rt([("r", 1, 2)])
    library.calibrate_rt([("r", 1, 2)], [])
    assert bridge["calls"][0]["stdin"] == bridge["calls"][1]["stdin"] == "r\t1.0\t2.0\n"


def test_a_missing_anchor_time_crosses_blank_for_mzlib_to_refuse(bridge):
    library.calibrate_rt([("r", None, 2), ("r", 1, math.nan)])
    assert bridge["calls"][0]["stdin"] == "r\t\t2.0\nr\t1.0\t\n"


def test_options_are_sent_only_when_given(bridge):
    library.calibrate_rt([("r", 1, 2)], kind="lowess", minimum_anchors=5, bandwidth=0.5, knots=7,
                         robustness_iterations=0, on_error="skip", timeout=9)
    call = bridge["calls"][0]
    assert call["args"] == [
        "library", "calibrate-rt", "--kind", "lowess", "--minimum-anchors", "5", "--bandwidth", "0.5",
        "--knots", "7", "--robustness-iterations", "0", "--on-error", "skip",
    ]
    assert call["timeout"] == 9


@pytest.mark.parametrize(
    "anchors, queries, options, message",
    [
        ([], None, {}, "at least one anchor"),
        ("r\t1\t2", None, {}, "must be a list"),
        ([("r", 1)], None, {}, "must be (group, rt_minutes, irt)"),
        ([{"rt_minutes": 1, "irt": 2}], None, {}, "needs a 'group' key"),
        ([("", 1, 2)], None, {}, "non-empty group name"),
        ([("a\tb", 1, 2)], None, {}, "tab or newline"),
        ([("--", 1, 2)], None, {}, "separates anchors from queries"),
        ([("r", "early", 2)], None, {}, "must be a number or None"),
        ([("r", True, 2)], None, {}, "must be a number or None"),
        ([("r", 1, 2)], [("r", 1, 2)], {}, "exactly one of rt_minutes and irt"),
        ([("r", 1, 2)], [("r", None, None)], {}, "exactly one of rt_minutes and irt"),
        ([("r", 1, 2)], [{"group": "r"}], {}, "exactly one of rt_minutes and irt"),
        ([("r", 1, 2)], None, {"kind": "cubic"}, "kind must be one of linear, lowess"),
        ([("r", 1, 2)], None, {"minimum_anchors": 2.5}, "minimum_anchors must be a whole number"),
        ([("r", 1, 2)], None, {"knots": True}, "knots must be a whole number"),
        ([("r", 1, 2)], None, {"bandwidth": "wide"}, "bandwidth must be a number"),
        ([("r", 1, 2)], None, {"on_error": "ignore"}, "on_error must be 'fail' or 'skip'"),
    ],
)
def test_bad_input_is_refused_before_the_bridge_starts(bridge, anchors, queries, options, message):
    with pytest.raises(_bridge.UsageError, match=None) as raised:
        library.calibrate_rt(anchors, queries, **options)
    assert message in str(raised.value)
    assert bridge["calls"] == []


# ---- what is parsed ------------------------------------------------------------------------------


def test_parses_a_real_run(bridge):
    bridge["payload"] = recorded("library_calibrate_rt_hela.json")
    calibration = library.calibrate_rt([("x", 1, 2)])
    assert isinstance(calibration, library.RtCalibration)
    assert (calibration.kind, calibration.minimum_anchors, calibration.bandwidth) == ("lowess", 20, 0.3)
    assert (calibration.knots, calibration.robustness_iterations, calibration.on_error) == (100, 3, "fail")
    assert (calibration.group_count, calibration.fitted_count, calibration.failed_count) == (1, 1, 0)
    assert (calibration.anchor_count, calibration.query_count, calibration.row_count) == (418, 228, 646)
    assert calibration.column_names == ["group", "role", "rt_minutes", "irt", "fitted_irt", "is_extrapolated"]
    run = calibration.group("Fig2HeLa-0-5h_MHRM_R01_T0")
    assert run is calibration.groups[0] and run.ok and run.error is None
    assert run.kind == "lowess" and run.anchor_count == 418 and run.query_count == 228
    assert run.residual_sd == pytest.approx(5.0035, abs=1e-4)
    assert run.first_anchor_minutes == pytest.approx(15.988, abs=1e-3)
    assert run.last_anchor_minutes == pytest.approx(40.903, abs=1e-3)
    rows = calibration.records
    assert len(rows) == 646
    assert {r["role"] for r in rows} == {"anchor", "to_minutes", "to_irt"}
    assert all(r["fitted_irt"] is None for r in rows if r["role"] != "anchor")
    assert all(isinstance(r["fitted_irt"], float) for r in rows if r["role"] == "anchor")
    assert calibration.failed_groups == []
    with pytest.raises(KeyError):
        calibration.group("no such run")


def test_parses_refused_runs(bridge):
    bridge["payload"] = recorded("library_calibrate_rt_pbmc_skip.json")
    calibration = library.calibrate_rt([("x", 1, 2)], on_error="skip")
    assert calibration.on_error == "skip"
    assert calibration.failed_count == 3 and calibration.fitted_count == 0
    assert [g.anchor_count for g in calibration.groups] == [13, 11, 10]
    for group in calibration.groups:
        assert not group.ok and group.kind is None and group.residual_sd is None
        assert group.first_anchor_minutes is None and group.last_anchor_minutes is None
        assert isinstance(group.error, library.FitFailure)
        assert (group.error.kind, group.error.type) == ("usage", "ArgumentException")
        assert group.error.message.startswith("Calibration needs at least 20 anchors")
    assert calibration.failed_groups == calibration.groups
    assert calibration.columns["fitted_irt"] == [None] * calibration.row_count
    assert calibration.columns["is_extrapolated"] == [None] * calibration.row_count
    assert any(c.startswith("3 groups were not fitted") for c in calibration.caveats)


def test_parses_lines_per_run(bridge):
    calibration = library.calibrate_rt([("x", 1, 2)], kind="linear", minimum_anchors=10)
    assert calibration.kind == "linear" and calibration.minimum_anchors == 10
    assert [g.kind for g in calibration.groups] == ["linear"] * 3
    queries = [r for r in calibration.records if r["role"] == "to_minutes"]
    assert [round(r["rt_minutes"], 2) for r in queries] == [9.22, 9.45, 9.17]
    assert all(r["irt"] == 11.886096 for r in queries)


def test_kinds_are_mzlibs():
    assert library.CALIBRATION_KINDS == ("linear", "lowess")
