"""pymzlib.stats against recorded bridge payloads: parsing, the stdin it writes, and its refusals.

The numbers themselves are mzLib's, and are held to limma and metafor by the bridge's C# tests and
by test_stats_live.py. These tests hold the Python layer: that it sends what the wire expects,
parses what the wire returns, and refuses bad input before starting a process.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, stats

FIXTURES = Path(__file__).parent / "fixtures"


def recorded(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture
def bridge(monkeypatch):
    """Answer every bridge call with ``state["payload"]``, recording the arguments and stdin."""
    state = {"calls": [], "payload": None}

    def invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": list(args), "stdin": stdin, "timeout": timeout})
        return state["payload"]

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    return state


# ---- fit -----------------------------------------------------------------------------------------


def test_fit_sends_the_files_and_one_coefficient_per_stdin_line(bridge, tmp_path):
    bridge["payload"] = recorded("stats_fit_malat.json")
    responses, design = tmp_path / "r.tsv", tmp_path / "d.tsv"
    stats.fit(responses, design, ["malat_250ng", "malat_125ng"], trend=True, spline_basis=3, threads=-1)

    call = bridge["calls"][0]
    assert call["args"][:2] == ["stats", "fit"]
    assert call["args"][call["args"].index("--responses") + 1] == str(responses)
    assert call["args"][call["args"].index("--design") + 1] == str(design)
    assert call["args"][call["args"].index("--threads") + 1] == "-1"
    assert call["args"][call["args"].index("--spline-basis") + 1] == "3"
    assert "--trend" in call["args"]
    assert call["stdin"] == "malat_250ng\nmalat_125ng\n"


def test_fit_omits_trend_options_it_was_not_given(bridge):
    bridge["payload"] = recorded("stats_fit_malat.json")
    stats.fit("r.tsv", "d.tsv", ["malat_250ng"])
    args = bridge["calls"][0]["args"]
    assert "--trend" not in args and "--spline-basis" not in args
    assert args[args.index("--threads") + 1] == "1"


def test_fit_parses_the_recorded_dilution(bridge):
    bridge["payload"] = recorded("stats_fit_malat.json")
    fit = stats.fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
                    ["malat_250ng", "malat_125ng"])

    assert isinstance(fit, stats.ModeratedFit)
    assert fit.tested_coefficients == ["malat_250ng", "malat_125ng"]
    assert fit.coefficient_names == ["intercept", "malat_250ng", "malat_125ng"]
    assert fit.row_count == fit.feature_count * 2 == 984
    assert sum(fit.status_counts.values()) == fit.feature_count
    assert set(fit.status_counts) == set(stats.FIT_STATUSES)
    assert fit.residual_df_differ is True
    assert isinstance(fit.prior, stats.VariancePrior)
    assert fit.prior.trended is False and fit.prior.df_infinite is False
    assert fit.prior.scale == pytest.approx(fit.columns["prior_variance"][0])
    assert any(c.startswith("residual_df_differ:") for c in fit.caveats)
    assert fit.column_names == list(fit.columns)


def test_unfitted_rows_carry_none_and_stay_out_of_the_family(bridge):
    bridge["payload"] = recorded("stats_fit_malat.json")
    fit = stats.fit("r.tsv", "d.tsv", ["malat_250ng", "malat_125ng"])
    for row in fit.records:
        fitted = row["status"] == "fitted"
        assert (row["p_value"] is not None) is fitted
        assert (row["bh_adjusted"] is not None) is fitted
        assert row["observed"] >= 0


def test_rows_selects_one_coefficient_in_feature_order(bridge):
    bridge["payload"] = recorded("stats_fit_malat.json")
    fit = stats.fit("r.tsv", "d.tsv", ["malat_250ng", "malat_125ng"])
    rows = fit.rows("malat_125ng")
    assert len(rows) == fit.feature_count
    assert {r["coefficient"] for r in rows} == {"malat_125ng"}
    assert [r["feature"] for r in rows] == fit.columns["feature"][: fit.feature_count]
    with pytest.raises(pymzlib.UsageError, match="was not tested"):
        fit.rows("intercept")


@pytest.mark.parametrize(
    "coefficients, message",
    [
        ("malat_250ng", "a list"),
        ([], "at least one coefficient"),
        ([""], "non-empty str"),
        ([3], "non-empty str"),
        (["a\nb"], "newline"),
    ],
)
def test_fit_refuses_bad_coefficients_without_starting_a_process(bridge, coefficients, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.fit("r.tsv", "d.tsv", coefficients)
    assert bridge["calls"] == []


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"trend": "yes"}, "trend must be True or False"),
        ({"threads": 1.5}, "threads must be an int"),
        ({"threads": True}, "threads must be an int"),
        ({"spline_basis": "4", "trend": True}, "spline_basis must be an int"),
    ],
)
def test_fit_refuses_bad_options_without_starting_a_process(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.fit("r.tsv", "d.tsv", ["x"], **kwargs)
    assert bridge["calls"] == []


# ---- fit: contrasts, confidence and the prior estimator -----------------------------------------


def test_fit_sends_contrasts_confidence_and_estimator_only_when_given(bridge, tmp_path):
    bridge["payload"] = recorded("stats_fit_malat_contrasts.json")
    contrasts = tmp_path / "contrasts_1.tsv"
    stats.fit("r.tsv", "d.tsv", contrasts=contrasts, confidence=0.9, estimator="marginal_likelihood")
    args = bridge["calls"][0]["args"]
    assert args[args.index("--contrasts") + 1] == str(contrasts)
    assert args[args.index("--confidence") + 1] == "0.9"
    assert args[args.index("--estimator") + 1] == "marginal_likelihood"
    assert bridge["calls"][0]["stdin"] == "", "contrasts alone: no coefficient on stdin"

    bridge["payload"] = recorded("stats_fit_malat.json")
    stats.fit("r.tsv", "d.tsv", ["malat_250ng"])
    args = bridge["calls"][1]["args"]
    assert not {"--contrasts", "--confidence", "--estimator"} & set(args), (
        "left out, they are not sent, so mzLib's defaults apply"
    )


def test_fit_parses_the_recorded_contrasts(bridge):
    bridge["payload"] = recorded("stats_fit_malat_contrasts.json")
    fit = stats.fit("malat_dilution_log2_vs_fluc.tsv", "design_1.tsv", contrasts="contrasts_1.tsv")

    assert fit.tested_coefficients == [] and fit.row_count == 0 and fit.records == []
    assert fit.contrasts_file.endswith("contrasts_1.tsv")
    pairs = fit.contrasts
    assert isinstance(pairs, stats.ModeratedContrasts)
    assert pairs.tested_contrasts == ["c1", "c2", "c3"]
    assert pairs.confidence_level == 0.95
    assert fit.coefficient_names == ["intercept", "malat_load=250ng", "malat_load=500ng"]
    assert pairs.weights == {"c1": [0, 1, 0], "c2": [0, 0, 1], "c3": [0, -1, 1]}
    assert pairs.row_count == 3 * fit.feature_count == len(pairs.records)
    assert pairs.column_names[:7] == [
        "feature", "contrast", "status", "estimate", "standard_error", "confidence_low", "confidence_high",
    ]
    rows = pairs.rows("c3")
    assert len(rows) == fit.feature_count
    assert [r["feature"] for r in rows] == pairs.columns["feature"][: fit.feature_count]
    for row in rows:
        fitted = row["status"] == "fitted"
        assert (row["confidence_low"] is not None) is fitted
        assert (row["bh_adjusted"] is not None) is fitted
        if fitted:
            assert row["confidence_low"] < row["estimate"] < row["confidence_high"]
    assert any(c.startswith("Contrasts are tested with the same variance prior") for c in fit.caveats)


def test_contrast_rows_refuses_a_contrast_that_was_not_tested(bridge):
    bridge["payload"] = recorded("stats_fit_malat_contrasts.json")
    pairs = stats.fit("r.tsv", "d.tsv", contrasts="c.tsv").contrasts
    with pytest.raises(pymzlib.UsageError, match="tested the contrasts: c1, c2, c3"):
        pairs.rows("c4")


def test_a_fit_without_contrasts_has_none(bridge):
    bridge["payload"] = recorded("stats_fit_malat.json")
    fit = stats.fit("r.tsv", "d.tsv", ["malat_250ng", "malat_125ng"])
    assert fit.contrasts is None and fit.contrasts_file is None
    assert fit.prior.estimator == "moments_legacy"


def test_fit_parses_the_estimator_mzlib_used(bridge):
    bridge["payload"] = recorded("stats_fit_malat_marginal.json")
    fit = stats.fit("r.tsv", "d.tsv", ["malat_250ng", "malat_125ng"], estimator="marginal_likelihood")
    assert fit.prior.estimator == "marginal_likelihood"
    assert fit.prior.estimator in stats.PRIOR_ESTIMATORS
    assert fit.caveats[0].startswith("marginal_likelihood:")


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"estimator": "reml"}, "estimator must be one of moments_legacy, marginal_likelihood"),
        ({"estimator": 1}, "estimator must be one of"),
        ({"confidence": "0.9"}, "confidence must be a number"),
        ({"confidence": True}, "confidence must be a number"),
        ({"contrasts": 3}, "contrasts must be a path"),
    ],
)
def test_fit_refuses_bad_contrast_options_without_starting_a_process(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.fit("r.tsv", "d.tsv", ["x"], **kwargs)
    assert bridge["calls"] == []


def test_fit_needs_a_coefficient_or_contrasts(bridge):
    with pytest.raises(pymzlib.UsageError, match="or a contrasts file"):
        stats.fit("r.tsv", "d.tsv")
    assert bridge["calls"] == []


def _older(name: str) -> dict:
    """A recording as a bridge from before pyMzLib 0.6.0 answers: no contrasts, no estimator."""
    data = recorded(name)
    data.pop("contrasts_file", None)
    data.pop("contrasts", None)
    data["prior"].pop("estimator", None)
    return data


def test_fit_says_so_when_an_older_bridge_ignored_contrasts(bridge):
    bridge["payload"] = _older("stats_fit_malat.json")
    with pytest.raises(pymzlib.UsageError, match=r"fit\(contrasts=\.\.\.\) needs the bridge from pyMzLib 0\.6\.0"):
        stats.fit("r.tsv", "d.tsv", ["malat_250ng"], contrasts="c.tsv")


def test_fit_says_so_when_an_older_bridge_ignored_the_estimator(bridge):
    bridge["payload"] = _older("stats_fit_malat.json")
    with pytest.raises(pymzlib.UsageError, match=r"fit\(estimator=\.\.\.\) needs the bridge from pyMzLib 0\.6\.0"):
        stats.fit("r.tsv", "d.tsv", ["malat_250ng"], estimator="moments_legacy")
    # Without the new options, an older bridge's answer is read as it always was.
    assert stats.fit("r.tsv", "d.tsv", ["malat_250ng"]).prior.estimator is None


def test_fit_names_the_cause_when_an_older_bridge_wanted_a_coefficient(monkeypatch):
    def invoke(*args, stdin=None, timeout=None):
        raise _bridge.UsageError("No coefficient to test was given on stdin. Give one name per line.")

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    with pytest.raises(pymzlib.UsageError, match="needs the bridge from pyMzLib 0.6.0"):
        stats.fit("r.tsv", "d.tsv", contrasts="c.tsv")


# ---- adjust --------------------------------------------------------------------------------------


def test_adjust_writes_one_line_per_value_with_untested_blank(bridge):
    bridge["payload"] = recorded("stats_adjust.json")
    adjusted = stats.adjust([0.0002, 0.004, 0.019, None, 0.031, 0.2, float("nan"), 0.74])

    assert bridge["calls"][0]["args"] == ["stats", "adjust"]
    assert bridge["calls"][0]["stdin"] == "0.0002\n0.004\n0.019\n\n0.031\n0.2\n\n0.74\n"
    assert adjusted.line_count == adjusted.row_count == 8
    assert adjusted.tested_count == 6
    assert adjusted.bh_adjusted[3] is None and adjusted.bh_adjusted[6] is None
    assert adjusted.bh_adjusted[0] == pytest.approx(0.0012)
    assert adjusted.bh_adjusted == adjusted.columns["bh_adjusted"]


def test_adjust_keeps_a_trailing_untested_entry(bridge):
    bridge["payload"] = recorded("stats_adjust.json")
    stats.adjust([0.1, None])
    assert bridge["calls"][0]["stdin"] == "0.1\n\n"


@pytest.mark.parametrize(
    "values, message",
    [([], "at least one"), ("0.1", "as a list"), ([0.1, "0.2"], r"p_values\[1\]"), ([True], r"p_values\[0\]")],
)
def test_adjust_refuses_bad_input_without_starting_a_process(bridge, values, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.adjust(values)
    assert bridge["calls"] == []


# ---- meta ----------------------------------------------------------------------------------------


def test_meta_writes_tab_separated_studies_from_tuples_and_mappings(bridge):
    bridge["payload"] = recorded("stats_meta_metafor.json")
    pooled = stats.meta(
        [("A", 0.5, 0.1), {"feature": "A", "estimate": 0.7, "standard_error": 0.2}, ["B", -1, 0.3]],
        confidence=0.9,
    )

    call = bridge["calls"][0]
    assert call["args"] == ["stats", "meta", "--confidence", "0.9"]
    assert call["stdin"] == "A\t0.5\t0.1\nA\t0.7\t0.2\nB\t-1.0\t0.3\n"
    assert isinstance(pooled, stats.MetaAnalysis)
    assert pooled.feature_count == pooled.row_count == len(pooled.columns["feature"])
    one_study = [r for r in pooled.records if r["studies"] == 1]
    assert all(r["q"] is None and r["i_squared"] is None for r in one_study)


def test_meta_parses_the_metafor_recording(bridge):
    bridge["payload"] = recorded("stats_meta_metafor.json")
    pooled = stats.meta([("x", 0.1, 0.1)])
    assert pooled.columns["feature"] == ["homogeneous", "heterogeneous", "two_studies", "random_eight"]
    assert pooled.study_count == 21
    assert pooled.confidence == 0.95
    assert all(not math.isnan(v) for v in pooled.columns["estimate"])


@pytest.mark.parametrize(
    "studies, message",
    [
        ([], "at least one study"),
        ("A\t1\t0.1", "a list of studies"),
        ({"feature": "A"}, "a list of studies"),
        ([("A", 1)], r"studies\[0\] must be"),
        ([{"feature": "A", "estimate": 1}], "standard_error"),
        ([("", 1, 0.1)], "non-empty feature"),
        ([("A\tB", 1, 0.1)], "tab"),
        ([("A\nB", 1, 0.1)], "newline"),
        ([("A", "1", 0.1)], "estimate must be a number"),
        ([("A", 1, None)], "standard_error must be a number"),
    ],
)
def test_meta_refuses_bad_studies_without_starting_a_process(bridge, studies, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.meta(studies)
    assert bridge["calls"] == []


def test_meta_refuses_a_non_numeric_confidence(bridge):
    with pytest.raises(pymzlib.UsageError, match="confidence must be a number"):
        stats.meta([("A", 1, 0.1)], confidence="high")


def test_module_lists_only_its_public_api():
    assert dir(stats) == sorted(stats.__all__)
    assert pymzlib.stats is stats
