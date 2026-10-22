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
        ("stats_fit_malat_marginal.json", lambda: stats.fit(
            STATS / "malat_dilution_log2_vs_fluc.tsv", STATS / "malat_dilution_design.tsv",
            ["malat_250ng", "malat_125ng"], estimator="marginal_likelihood")),
        ("stats_adjust.json", lambda: stats.adjust(
            [0.0002, 0.004, 0.019, None, 0.031, 0.2, float("nan"), 0.74])),
    ],
)
def test_the_recordings_are_what_the_bridge_returns_today(fixture, call):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    live = call()
    # The marginal-likelihood prior is a minimum found by Brent's search over log d0. Near a
    # minimum the likelihood is flat, so the last bits of Math.Log (MathNet's GammaLn is built on
    # it), which differ between platforms' math libraries, move it by about sqrt(eps) relative:
    # more than 1e-9, far less than the 1e-6 the mixed-model recordings are held to likewise.
    relative = 1e-6 if "marginal" in fixture else 1e-9
    assert_same_columns(live.columns, recorded["columns"], relative)
    assert live.caveats == recorded["caveats"]
    if "status_counts" in recorded:
        # Its keys come from an mzLib enum; a repin that grows it must not slip past the recording.
        assert live.status_counts == recorded["status_counts"]


def test_the_meta_recording_is_what_the_bridge_returns_today():
    recorded = json.loads((FIXTURES / "stats_meta_metafor.json").read_text(encoding="utf-8"))
    studies = [(r["case"], float(r["yi"]), float(r["sei"])) for r in tsv("metafor_dl_inputs.tsv")]
    assert_same_columns(stats.meta(studies).columns, recorded["columns"])


def assert_same_columns(live: dict, recorded: dict, relative: float = 1e-9) -> None:
    """Equal, except that floats may differ in their last bits.

    The recordings were made on one platform and CI replays them on others; the last digits of a
    log or a quantile can differ between them, and that is not drift. 1e-9 relative is far below
    the 1e-8 the reference comparisons promise, and far above floating-point noise.
    """
    assert list(live) == list(recorded)
    for name, values in recorded.items():
        assert len(live[name]) == len(values), name
        for i, (a, b) in enumerate(zip(live[name], values)):
            if isinstance(b, float) and isinstance(a, (int, float)) and not isinstance(a, bool):
                assert abs(a - b) <= relative * max(abs(b), 1e-300), (name, i, a, b)
            else:
                assert a == b, (name, i, a, b)


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


# ---- contrasts and the prior estimator ------------------------------------------------------------


def pairs_design(tmp_path: Path):
    """diff design on the dilution series with no reference level: every pair of loads a contrast."""
    from pymzlib import diff

    design = diff.design(FIXTURES / "diff" / "malat_dilution_samples.tsv",
                         factors=[diff.FactorSpec("malat_load")], out=tmp_path / "malat_pairs")
    return design.strata[0]


def test_the_contrasts_recording_is_what_the_chain_returns_today(tmp_path):
    """diff design --out, then stats fit --contrasts on what it wrote: the guide's chain, for real."""
    recorded = json.loads((FIXTURES / "stats_fit_malat_contrasts.json").read_text(encoding="utf-8"))
    stratum = pairs_design(tmp_path)
    live = stats.fit(STATS / "malat_dilution_log2_vs_fluc.tsv", stratum.design_file,
                     contrasts=stratum.contrasts_file)
    assert live.contrasts.tested_contrasts == recorded["contrasts"]["tested_contrasts"]
    assert live.contrasts.weights == recorded["contrasts"]["weights"]
    assert live.contrasts.confidence_level == recorded["contrasts"]["confidence_level"]
    assert_same_columns(live.contrasts.columns, recorded["contrasts"]["columns"])
    assert live.caveats == recorded["caveats"]


def test_a_contrast_of_one_coefficient_is_that_coefficients_test(tmp_path):
    stratum = pairs_design(tmp_path)
    fit = stats.fit(STATS / "malat_dilution_log2_vs_fluc.tsv", stratum.design_file,
                    ["malat_load=250ng"], contrasts=stratum.contrasts_file)
    assert fit.contrasts.weights["c1"] == [0, 1, 0]
    coefficient = fit.rows("malat_load=250ng")
    contrast = fit.contrasts.rows("c1")
    for ours, theirs in zip(coefficient, contrast):
        for column in ("estimate", "standard_error", "t", "p_value", "bh_adjusted"):
            assert ours[column] == theirs[column], column


def test_a_lower_confidence_gives_narrower_intervals_and_mzlib_refuses_one_outside_0_1(tmp_path):
    stratum = pairs_design(tmp_path)
    args = (STATS / "malat_dilution_log2_vs_fluc.tsv", stratum.design_file)
    wide = stats.fit(*args, contrasts=stratum.contrasts_file)
    narrow = stats.fit(*args, contrasts=stratum.contrasts_file, confidence=0.8)
    assert narrow.contrasts.confidence_level == 0.8
    for w, n in zip(wide.contrasts.records, narrow.contrasts.records):
        if w["status"] == "fitted":
            assert w["confidence_low"] < n["confidence_low"] < n["confidence_high"] < w["confidence_high"]
    with pytest.raises(_bridge.UsageError, match="strictly between 0 and 1"):
        stats.fit(*args, contrasts=stratum.contrasts_file, confidence=1.5)
    with pytest.raises(_bridge.UsageError, match="give --contrasts as well"):
        stats.fit(*args, ["malat_load=250ng"], confidence=0.8)


def test_the_default_estimator_is_the_one_mzlib_reports(tmp_path):
    args = (STATS / "malat_dilution_log2_vs_fluc.tsv", STATS / "malat_dilution_design.tsv", ["malat_250ng"])
    default = stats.fit(*args)
    legacy = stats.fit(*args, estimator="moments_legacy")
    assert default.prior.estimator == "moments_legacy"
    assert default.columns["t"] == legacy.columns["t"]
