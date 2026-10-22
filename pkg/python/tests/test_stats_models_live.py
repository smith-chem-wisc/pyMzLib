"""stats fit-logistic, fit-mixed and entrapment-fdp through the real, packaged bridge.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drift
apart. These run the real executable on the files in ``fixtures/stats/`` and check that what it
returns today is what was recorded (so the docs' printed numbers are still true), that it still
agrees with R's glm, nlme's lme and the entrapment reference script from Python, and that a refusal
crosses as a usage error. Local files only.
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


def relative(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


def test_fit_logistic_matches_glm_from_python():
    fit = stats.fit_logistic(STATS / "glm_reference_responses.tsv", STATS / "glm_reference_design.tsv",
                             ["Intercept", "x", "grp"])
    glm = tsv("glm_reference_fit.tsv")
    for j, coefficient in enumerate(["Intercept", "x", "grp"]):
        for ours, theirs in zip(fit.rows(coefficient), glm):
            if theirs["status"] == "separated":
                assert ours["status"] == "separated"
                continue
            assert ours["status"] == "fitted"
            assert relative(ours["estimate"], float(theirs[f"b{j}"])) < 1e-8
            assert relative(ours["standard_error"], float(theirs[f"se{j}"])) < 1e-6
    for ours, theirs in zip(fit.rows("x"), glm):
        if theirs["status"] == "fitted":
            assert relative(ours["p_value"], float(theirs["p1"])) < 1e-6
            assert relative(ours["deviance"], float(theirs["deviance"])) < 1e-8


@pytest.mark.parametrize("method, estimator", [("REML", None), ("ML", "ml")])
def test_fit_mixed_matches_nlme_from_python(method, estimator):
    fit = stats.fit_mixed(STATS / "lmm_reference_responses.tsv", STATS / "lmm_reference_design.tsv",
                          STATS / "lmm_reference_groups.tsv", ["Intercept", "age", "wp"], estimator=estimator)
    assert fit.estimator == method.lower()
    nlme = [r for r in tsv("lmm_reference_fit.tsv") if r["method"] == method]
    for j, coefficient in enumerate(["Intercept", "age", "wp"]):
        for ours, theirs in zip(fit.rows(coefficient), nlme):
            assert ours["status"] == "fitted" == theirs["status"]
            assert relative(ours["estimate"], float(theirs[f"b{j}"])) < 1e-5
            assert relative(ours["standard_error"], float(theirs[f"se{j}"])) < 1e-5
            assert ours["df"] == float(theirs[f"df{j}"])
    for ours, theirs in zip(fit.rows("age"), nlme):
        assert relative(ours["log_likelihood"], float(theirs["loglik"])) < 1e-5
        assert abs(ours["group_variance"] - float(theirs["tau2"])) < 1e-5 * float(theirs["sigma2"])


def test_entrapment_fdp_matches_the_reference_script_from_python():
    reference = tsv("entrapment_reference_r1.tsv")
    fdp = stats.entrapment_fdp(STATS / "entrapment_reference_targets.tsv",
                               STATS / "entrapment_reference_entrapments.tsv", ratio=1.0,
                               thresholds=[float(r["q_value_threshold"]) for r in reference])
    assert fdp.row_count == len(reference)
    for ours, theirs in zip(fdp.records, reference):
        assert ours["target_count"] == int(theirs["target_count"])
        assert ours["entrapment_count"] == int(theirs["entrapment_count"])
        for column in ("lower_bound", "combined", "paired"):
            assert abs(ours[column] - float(theirs[column])) < 1e-5, column


def assert_same_columns(live: dict, recorded: dict) -> None:
    """Equal, except that floats may differ in their last bits across platforms (1e-9 relative)."""
    assert list(live) == list(recorded)
    for name, values in recorded.items():
        assert len(live[name]) == len(values), name
        for i, (a, b) in enumerate(zip(live[name], values)):
            if isinstance(b, float) and isinstance(a, (int, float)) and not isinstance(a, bool):
                assert abs(a - b) <= 1e-9 * max(abs(b), 1e-300), (name, i, a, b)
            else:
                assert a == b, (name, i, a, b)


@pytest.mark.parametrize(
    "fixture, call",
    [
        ("stats_fit_logistic_glm.json", lambda: stats.fit_logistic(
            STATS / "glm_reference_responses.tsv", STATS / "glm_reference_design.tsv", ["Intercept", "x", "grp"])),
        ("stats_fit_logistic_malat.json", lambda: stats.fit_logistic(
            STATS / "malat_dilution_detected.tsv", STATS / "malat_dilution_design.tsv",
            ["malat_250ng", "malat_125ng"])),
        ("stats_fit_mixed_lmm_reml.json", lambda: stats.fit_mixed(
            STATS / "lmm_reference_responses.tsv", STATS / "lmm_reference_design.tsv",
            STATS / "lmm_reference_groups.tsv", ["Intercept", "age", "wp"])),
        ("stats_fit_mixed_lmm_ml.json", lambda: stats.fit_mixed(
            STATS / "lmm_reference_responses.tsv", STATS / "lmm_reference_design.tsv",
            STATS / "lmm_reference_groups.tsv", ["Intercept", "age", "wp"], estimator="ml")),
        ("stats_fit_mixed_malat.json", lambda: stats.fit_mixed(
            STATS / "malat_dilution_log2_vs_fluc.tsv", STATS / "malat_dilution_design.tsv",
            STATS / "malat_dilution_day.tsv", ["malat_250ng", "malat_125ng"])),
        ("stats_entrapment_fdp_reference.json", lambda: stats.entrapment_fdp(
            STATS / "entrapment_reference_targets.tsv", STATS / "entrapment_reference_entrapments.tsv",
            ratio=1.0)),
        ("stats_entrapment_fdp_reference_r2.json", lambda: stats.entrapment_fdp(
            STATS / "entrapment_reference_targets.tsv", STATS / "entrapment_reference_entrapments.tsv",
            ratio=2.0, thresholds=[0.01, 0.05])),
    ],
)
def test_the_recordings_are_what_the_bridge_returns_today(fixture, call):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    live = call()
    assert_same_columns(live.columns, recorded["columns"])
    assert live.caveats == recorded["caveats"]
    if "status_counts" in recorded:
        # Its keys are the bridge's reading of mzLib's code (mzLib#1487); a change must not slip past.
        assert live.status_counts == recorded["status_counts"]


def test_refusals_cross_as_usage_errors(tmp_path):
    proportions = tmp_path / "proportions.tsv"
    proportions.write_text("feature\tV1\tV2\nA\t0.5\t1\n", encoding="utf-8")
    design = tmp_path / "design.tsv"
    design.write_text("sample\tintercept\nV1\t1\nV2\t1\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="expected 0, 1 or missing.*feature 'A', sample 'V1'"):
        stats.fit_logistic(proportions, design, ["intercept"])
    with pytest.raises(_bridge.UsageError, match="refuses --max-iterations 0"):
        stats.fit_logistic(STATS / "glm_reference_responses.tsv", STATS / "glm_reference_design.tsv",
                           ["x"], max_iterations=0)
    groups = tmp_path / "groups.tsv"
    groups.write_text("sample\tday\tsite\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="needs exactly two"):
        stats.fit_mixed(STATS / "lmm_reference_responses.tsv", STATS / "lmm_reference_design.tsv", groups, ["age"])
    with pytest.raises(_bridge.UsageError, match="ratio must be a positive finite number"):
        stats.entrapment_fdp(STATS / "entrapment_reference_targets.tsv",
                             STATS / "entrapment_reference_entrapments.tsv", ratio=0.0)
    extra = tmp_path / "targets.tsv"
    extra.write_text("peptide\tscore\tq_value\nPEPTIDE\t10\t0.01\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="Not a column this file takes: peptide"):
        stats.entrapment_fdp(extra, STATS / "entrapment_reference_entrapments.tsv", ratio=1.0)
