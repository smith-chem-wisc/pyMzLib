"""The stats verbs through the real, packaged bridge, over the committed reference tables.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drifting
apart. These run the real executable on the files in ``fixtures/stats/`` and check two things:
that what it returns today is what was recorded (so the docs' printed numbers are still true), and
that it still agrees with limma and metafor, from Python, to 1e-8. Local files only, no network.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, stats

FIXTURES = Path(__file__).parent / "fixtures"
STATS = FIXTURES / "stats"


@pytest.fixture(autouse=True)
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so the stats verbs went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def tsv(name: str) -> list[dict[str, str]]:
    with open(STATS / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def worst_relative(ours, theirs) -> float:
    return max(abs(a - float(b)) / max(abs(float(b)), 1e-300) for a, b in zip(ours, theirs))


def test_fit_matches_limma_from_python():
    fit = stats.fit(STATS / "limma_reference_responses.tsv", STATS / "limma_reference_design.tsv",
                    ["age_decades"])
    limma = tsv("limma_ebayes_notrend.tsv")
    assert fit.feature_count == len(limma)
    for ours, theirs in [("t", "t_age"), ("p_value", "p_age"), ("bh_adjusted", "bh_age"),
                         ("posterior_variance", "s2_post"), ("df_total", "df_total"),
                         ("prior_variance", "s2_prior")]:
        assert worst_relative(fit.columns[ours], [row[theirs] for row in limma]) < 1e-8, ours


def test_meta_matches_metafor_from_python():
    studies = [(r["case"], float(r["yi"]), float(r["sei"])) for r in tsv("metafor_dl_inputs.tsv")]
    pooled = stats.meta(studies)
    results = {r["case"]: r for r in tsv("metafor_dl_results.tsv")}
    for row in pooled.records:
        expected = results[row["feature"]]
        for ours, theirs in [("estimate", "estimate"), ("standard_error", "se"), ("p_value", "pval"),
                             ("confidence_low", "ci_lb"), ("confidence_high", "ci_ub")]:
            assert worst_relative([row[ours]], [expected[theirs]]) < 1e-8, (row["feature"], ours)


@pytest.mark.parametrize(
    "fixture, call",
    [
        ("stats_fit_malat.json", lambda: stats.fit(
            STATS / "malat_dilution_log2_vs_fluc.tsv", STATS / "malat_dilution_design.tsv",
            ["malat_250ng", "malat_125ng"])),
        ("stats_fit_limma.json", lambda: stats.fit(
            STATS / "limma_reference_responses.tsv", STATS / "limma_reference_design.tsv",
            ["age_decades"])),
        ("stats_adjust.json", lambda: stats.adjust(
            [0.0002, 0.004, 0.019, None, 0.031, 0.2, float("nan"), 0.74])),
    ],
)
def test_the_recordings_are_what_the_bridge_returns_today(fixture, call):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    live = call()
    assert live.columns == recorded["columns"]
    assert live.caveats == recorded["caveats"]


def test_the_meta_recording_is_what_the_bridge_returns_today():
    recorded = json.loads((FIXTURES / "stats_meta_metafor.json").read_text(encoding="utf-8"))
    studies = [(r["case"], float(r["yi"]), float(r["sei"])) for r in tsv("metafor_dl_inputs.tsv")]
    assert stats.meta(studies).columns == recorded["columns"]


def test_a_trend_fit_runs_and_reports_its_basis():
    fit = stats.fit(STATS / "limma_reference_responses.tsv", STATS / "limma_reference_design.tsv",
                    ["age_decades"], trend=True)
    assert fit.prior.trended is True and fit.prior.scale is None
    assert 1 <= fit.prior.spline_basis_count <= 4


def test_bad_input_is_a_usage_error_from_the_bridge(tmp_path):
    design = tmp_path / "design.tsv"
    design.write_text("sample\tintercept\na\t1\nb\t1\n", encoding="utf-8")
    responses = tmp_path / "responses.tsv"
    responses.write_text("id\ta\tc\nf1\t1\t2\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="Not in the design: c"):
        stats.fit(responses, design, ["intercept"])
    with pytest.raises(_bridge.UsageError, match=r"outside \[0, 1\]"):
        stats.adjust([0.5, 2.0])
