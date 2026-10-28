"""library calibrate-rt through the real, packaged bridge.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drift
apart. These re-run each recording's call on the files in ``fixtures/calibration/`` and check that
what the bridge returns today is what was recorded (so the docs' printed numbers are still true),
and that mzLib's refusals cross as usage errors. Local files only.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, library

FIXTURES = Path(__file__).parent / "fixtures"
CALIBRATION = FIXTURES / "calibration"


@pytest.fixture(autouse=True)
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so calibrate-rt went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def tsv(name: str) -> list[dict[str, str]]:
    with open(CALIBRATION / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def hela_call() -> library.RtCalibration:
    """The guide's call: q-value 0.01 or below are anchors; the rest are queried both ways."""
    rows = tsv("diann2_hela_rt.tsv")
    others = [r for r in rows if float(r["q_value"]) > 0.01]
    anchors = [(r["run"], float(r["retention_time"]), float(r["indexed_retention_time"]))
               for r in rows if float(r["q_value"]) <= 0.01]
    queries = ([(r["run"], None, float(r["indexed_retention_time"])) for r in others]
               + [(r["run"], float(r["retention_time"]), None) for r in others])
    return library.calibrate_rt(anchors, queries)


def pbmc_inputs():
    """The guide's and the docstring's: one precursor held out of each run's anchors and queried."""
    rows = tsv("diann_pbmc_rt.tsv")
    held_out = "AAAAAAAAGAFAGR2"
    anchors = [(r["run"], float(r["retention_time"]), float(r["indexed_retention_time"]))
               for r in rows if r["precursor_id"] != held_out]
    queries = [(r["run"], None, 11.886096) for r in rows if r["precursor_id"] == held_out]
    return anchors, queries


def assert_same(live: dict, recorded: dict) -> None:
    """Equal, except that floats may differ in their last bits across platforms (1e-9 relative)."""
    if isinstance(recorded, dict):
        assert list(live) == list(recorded)
        for key in recorded:
            assert_same(live[key], recorded[key])
    elif isinstance(recorded, list):
        assert len(live) == len(recorded)
        for a, b in zip(live, recorded):
            assert_same(a, b)
    elif isinstance(recorded, float) and not isinstance(live, bool):
        assert abs(live - recorded) <= 1e-9 * max(abs(recorded), 1e-300), (live, recorded)
    else:
        assert live == recorded


@pytest.mark.parametrize(
    "fixture, call",
    [
        ("library_calibrate_rt_hela.json", hela_call),
        ("library_calibrate_rt_pbmc_skip.json",
         lambda: library.calibrate_rt(*pbmc_inputs(), on_error="skip")),
        ("library_calibrate_rt_pbmc_linear.json",
         lambda: library.calibrate_rt(*pbmc_inputs(), kind="linear", minimum_anchors=10)),
    ],
)
def test_the_recordings_are_what_the_bridge_returns_today(fixture, call, monkeypatch):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    payloads = []
    real = _bridge.invoke

    def keep(*args, **kwargs):
        data = real(*args, **kwargs)
        payloads.append(data)
        return data

    monkeypatch.setattr(_bridge, "invoke", keep)
    call()
    assert_same(payloads[-1], recorded)


def test_refusals_cross_as_usage_errors():
    anchors, queries = pbmc_inputs()
    with pytest.raises(_bridge.UsageError, match=r"Group '20240316_PBMCs.*needs at least 20 anchors; got 13"):
        library.calibrate_rt(anchors, queries)
    with pytest.raises(_bridge.UsageError, match="Knots must be at least 3"):
        library.calibrate_rt(anchors, queries, knots=2, on_error="skip")
    with pytest.raises(_bridge.UsageError, match="Bandwidth must be in"):
        library.calibrate_rt(anchors, bandwidth=1.5)
    with pytest.raises(_bridge.UsageError, match="must be finite"):
        library.calibrate_rt([("r", 1, 1), ("r", 2, None), ("r", 3, 3)], minimum_anchors=3)
