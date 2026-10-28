"""stats fit-latent-mixed and fit-link through the real, packaged bridge.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drift
apart. These run the real executable on the files in ``fixtures/stats/`` and check that what it
returns today is what was recorded (so the docs' printed numbers are still true), that it still
agrees with lme4's glmer and with glmmTMB from Python, and that a refusal crosses as a usage error.
Local files only.
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

#: LogitNormalMixedModel finds its variances by numerical search, polished by Newton steps, so the
#: last digits of a recording made on one platform may move on another (as the mixed-model ones do).
CROSS_PLATFORM_OPTIMISED = 1e-6

#: FitLink's estimate is Nelder-Mead's stopping point, with no polishing step: its last digits follow
#: the platform's math library further than a polished optimum's.
CROSS_PLATFORM_SIMPLEX = 1e-4


@pytest.fixture(autouse=True)
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so the latent-scale verbs went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def tsv(name: str) -> list[dict[str, str]]:
    with open(STATS / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def close(a: object, b: object, relative: float) -> bool:
    if isinstance(b, float) and isinstance(a, (int, float)) and not isinstance(a, bool):
        return abs(a - b) <= relative * max(abs(b), 1e-300)
    return a == b


def scatter(**kwargs):
    return stats.fit_latent_mixed(STATS / "lnmm_reference_design.tsv", STATS / "lnmm_reference_groups.tsv",
                                  ["Intercept", "trait", "x"], link="identity",
                                  counts=STATS / "lnmm_scatter_counts.tsv", threads=-1, **kwargs)


def link_ml():
    return stats.fit_latent_mixed(STATS / "lnmm_reference_design.tsv", STATS / "lnmm_reference_groups.tsv",
                                  ["Intercept", "trait", "x"], link=(-0.4, 0.6, 0.3),
                                  counts=STATS / "lnmm_link_counts.tsv", depth=STATS / "lnmm_link_depth.tsv",
                                  residual_variance=STATS / "lnmm_link_residual_variance.tsv", estimator="ml",
                                  threads=-1)


def simulated():
    return stats.fit_latent_mixed(STATS / "lnmm_simulated_design.tsv", STATS / "lnmm_simulated_groups.tsv",
                                  ["Intercept", "trait"], link=(0.4, 0.6),
                                  responses=STATS / "lnmm_simulated_responses.tsv",
                                  counts=STATS / "lnmm_simulated_counts.tsv", threads=-1)


@pytest.mark.parametrize(
    "fixture, call",
    [
        ("stats_fit_latent_mixed_scatter.json", scatter),
        ("stats_fit_latent_mixed_link_ml.json", link_ml),
        ("stats_fit_latent_mixed_simulated.json", simulated),
    ],
)
def test_the_latent_mixed_recordings_are_what_the_bridge_returns_today(fixture, call):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    live = call()
    assert list(live.columns) == list(recorded["columns"])
    for name, values in recorded["columns"].items():
        bad = [(i, a, b) for i, (a, b) in enumerate(zip(live.columns[name], values))
               if not close(a, b, CROSS_PLATFORM_OPTIMISED)]
        assert len(live.columns[name]) == len(values) and not bad, (name, bad[:3])
    assert live.caveats == recorded["caveats"]
    # Its keys are the bridge's reading of mzLib's code (mzLib#1487); a change must not slip past.
    assert live.status_counts == recorded["status_counts"]
    assert live.link == stats.DetectionLink._from_wire(recorded["link"])


def test_the_link_recording_is_what_the_bridge_returns_today():
    recorded = json.loads((FIXTURES / "stats_fit_link_simulated.json").read_text(encoding="utf-8"))
    live = stats.fit_link(STATS / "link_cells.tsv", with_depth=True)
    for name, value in recorded["link"].items():
        assert close(getattr(live.link, name), float(value), CROSS_PLATFORM_SIMPLEX), name
    assert live.caveats == recorded["caveats"]
    assert (live.with_depth, live.with_overdispersion, live.detection_only) == (
        recorded["with_depth"], recorded["with_overdispersion"], recorded["detection_only"])


def test_fit_latent_mixed_matches_glmer_and_glmmtmb_from_python():
    ml = link_ml()
    glmer = [r for r in tsv("lnmm_counted_fit.tsv") if r["link"] == "link" and r["method"] == "ML"]
    for j, coefficient in enumerate(["Intercept", "trait", "x"]):
        for ours, theirs in zip(ml.rows(coefficient), glmer):
            if theirs["status"] != "fitted":
                continue
            se = float(theirs[f"se{j}"])
            assert abs(ours["estimate"] - float(theirs[f"b{j}"])) < 1e-5 * se
            assert abs(ours["standard_error"] - se) / se < 5e-5
    reml = scatter()
    tmb = [r for r in tsv("lnmm_scatter_fit.tsv") if r["method"] == "REML"]
    for ours, theirs in zip(reml.rows("trait"), tmb):
        tau2 = float(theirs["tau2"])
        assert abs(ours["group_variance"] - tau2) < 0.03 * tau2 + 0.003


def test_refusals_cross_as_usage_errors(tmp_path):
    design, groups = STATS / "lnmm_simulated_design.tsv", STATS / "lnmm_simulated_groups.tsv"
    counts = tmp_path / "counts.tsv"
    counts.write_text("feature\ts1\nA\t6/5\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="must list exactly the responses' features"):
        stats.fit_latent_mixed(design, groups, ["trait"], link="identity",
                               responses=STATS / "lnmm_simulated_responses.tsv", counts=counts)
    with pytest.raises(_bridge.UsageError, match="non-zero slope"):
        stats.fit_latent_mixed(design, groups, ["trait"], link=(0.0, 0.0),
                               counts=STATS / "lnmm_simulated_counts.tsv")
    with pytest.raises(_bridge.UsageError, match="supply one depth per sample"):
        stats.fit_latent_mixed(design, groups, ["trait"], link=(0.0, 1.0, 0.2),
                               counts=STATS / "lnmm_simulated_counts.tsv")
    cells = tmp_path / "cells.tsv"
    cells.write_text("latent_mean\tlatent_sd\tsuccesses\ttrials\tfeature\n0\t1\t1\t2\tA\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="Not a column this file takes: feature"):
        stats.fit_link(cells)
