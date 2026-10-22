"""pymzlib.stats.fit_logistic, fit_mixed and entrapment_fdp against recorded bridge payloads.

The numbers are mzLib's, held to R's glm, nlme's lme and the entrapment reference script by the
bridge's C# tests and by test_stats_models_live.py. These tests hold the Python layer: that it sends
what the wire expects (and nothing for an option left out, so mzLib's default applies), parses what
the wire returns, and refuses bad input before starting a process.
"""

from __future__ import annotations

import json
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
    state = {"calls": [], "payload": None, "verbs": []}

    def invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": list(args), "stdin": stdin, "timeout": timeout})
        return state["payload"]

    def require_verb(verb, since):
        state["verbs"].append((verb, since))

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", require_verb)
    return state


# ---- fit_logistic --------------------------------------------------------------------------------


def test_fit_logistic_sends_only_what_was_given(bridge, tmp_path):
    bridge["payload"] = recorded("stats_fit_logistic_glm.json")
    responses, design = tmp_path / "r.tsv", tmp_path / "d.tsv"
    stats.fit_logistic(responses, design, ["x", "grp"])
    call = bridge["calls"][0]
    assert call["args"] == ["stats", "fit-logistic", "--responses", str(responses),
                            "--design", str(design), "--threads", "1"]
    assert call["stdin"] == "x\ngrp\n"
    assert bridge["verbs"] == [("stats fit-logistic", "0.6.0")]
    stats.fit_logistic(responses, design, ["x"], max_iterations=25, tolerance=1e-8, threads=-1)
    assert bridge["calls"][1]["args"][-6:] == ["--threads", "-1", "--max-iterations", "25",
                                               "--tolerance", "1e-08"]


def test_fit_logistic_parses_the_glm_recording(bridge):
    bridge["payload"] = recorded("stats_fit_logistic_glm.json")
    fit = stats.fit_logistic("r.tsv", "d.tsv", ["Intercept", "x", "grp"])
    assert isinstance(fit, stats.LogisticModelFit)
    assert fit.tested_coefficients == ["Intercept", "x", "grp"]
    assert fit.max_iterations is None and fit.tolerance is None
    assert list(fit.status_counts) == list(stats.LOGISTIC_STATUSES)
    assert sum(fit.status_counts.values()) == fit.feature_count == 80
    assert fit.row_count == 3 * 80
    assert fit.column_names == ["feature", "coefficient", "status", "estimate", "standard_error", "z",
                                "p_value", "deviance", "iterations", "observed", "events"]
    for row in fit.rows("x"):
        assert row["status"] in stats.LOGISTIC_STATUSES
        fitted = row["status"] == "fitted"
        assert all((row[c] is not None) == fitted
                   for c in ("estimate", "standard_error", "z", "p_value", "deviance"))
        assert 0 <= row["events"] <= row["observed"] <= fit.sample_count
    with pytest.raises(pymzlib.UsageError, match="'age' was not reported by this fit"):
        fit.rows("age")


def test_fit_logistic_separated_rows_say_so_in_the_caveats(bridge):
    bridge["payload"] = recorded("stats_fit_logistic_malat.json")
    fit = stats.fit_logistic("r.tsv", "d.tsv", ["malat_250ng", "malat_125ng"])
    assert fit.status_counts["separated"] > 0
    assert any(c.startswith("separated: the design splits") for c in fit.caveats)
    assert any("were not fitted" in c and "separated" in c for c in fit.caveats)


@pytest.mark.parametrize(
    "args, kwargs, message",
    [
        (("r.tsv", "d.tsv", "x"), {}, r"for one, use \['x'\]"),
        (("r.tsv", "d.tsv", []), {}, "at least one coefficient"),
        (("r.tsv", "d.tsv", [""]), {}, "non-empty str"),
        (("r.tsv", "d.tsv", ["a\nb"]), {}, "newline"),
        ((5, "d.tsv", ["x"]), {}, "responses must be a path"),
        (("r.tsv", "d.tsv", ["x"]), {"max_iterations": 2.5}, "max_iterations must be an int"),
        (("r.tsv", "d.tsv", ["x"]), {"max_iterations": True}, "max_iterations must be an int"),
        (("r.tsv", "d.tsv", ["x"]), {"tolerance": "tight"}, "tolerance must be a number"),
        (("r.tsv", "d.tsv", ["x"]), {"threads": "all"}, "threads must be an int"),
    ],
)
def test_fit_logistic_refuses_bad_input_without_starting_a_process(bridge, args, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.fit_logistic(*args, **kwargs)
    assert bridge["calls"] == []


# ---- fit_mixed -----------------------------------------------------------------------------------


def test_fit_mixed_sends_the_estimator_only_when_given(bridge, tmp_path):
    bridge["payload"] = recorded("stats_fit_mixed_lmm_reml.json")
    r, d, g = tmp_path / "r.tsv", tmp_path / "d.tsv", tmp_path / "g.tsv"
    stats.fit_mixed(r, d, g, ["age"])
    assert bridge["calls"][0]["args"] == ["stats", "fit-mixed", "--responses", str(r), "--design", str(d),
                                          "--groups", str(g), "--threads", "1"]
    assert bridge["calls"][0]["stdin"] == "age\n"
    assert bridge["verbs"] == [("stats fit-mixed", "0.6.0")]
    stats.fit_mixed(r, d, g, ["age"], estimator="ml")
    assert bridge["calls"][1]["args"][-2:] == ["--estimator", "ml"]


def test_fit_mixed_parses_the_reference_recordings(bridge):
    bridge["payload"] = recorded("stats_fit_mixed_lmm_reml.json")
    fit = stats.fit_mixed("r.tsv", "d.tsv", "g.tsv", ["Intercept", "age", "wp"])
    assert isinstance(fit, stats.MixedModelFit)
    assert fit.estimator == "reml", "absent, mzLib's default, as mzLib reports it"
    assert fit.grouping == "group" and fit.group_count == 7
    assert fit.group_names == ["d1", "d2", "d3", "d4", "d5", "d6", "d7"]
    assert list(fit.status_counts) == list(stats.MIXED_STATUSES)
    assert sum(fit.status_counts.values()) == fit.feature_count == 120
    assert fit.column_names == ["feature", "coefficient", "status", "estimate", "standard_error", "df",
                                "t", "p_value", "residual_variance", "group_variance",
                                "log_likelihood", "observed", "groups"]
    for row in fit.rows("wp"):
        assert row["df"] in (4, 5) and row["groups"] in (6, 7)
        assert row["group_variance"] >= 0
    bridge["payload"] = recorded("stats_fit_mixed_lmm_ml.json")
    ml = stats.fit_mixed("r.tsv", "d.tsv", "g.tsv", ["Intercept", "age", "wp"], estimator="ml")
    assert ml.estimator == "ml"
    assert any(c.startswith("ml: residual_variance is the maximum-likelihood estimate") for c in ml.caveats)


def test_fit_mixed_null_df_rows_carry_t_but_no_p_value(bridge):
    bridge["payload"] = recorded("stats_fit_mixed_malat.json")
    fit = stats.fit_mixed("r.tsv", "d.tsv", "g.tsv", ["malat_250ng", "malat_125ng"])
    no_df = [r for r in fit.records if r["status"] == "fitted" and r["df"] is None]
    assert no_df, "the dilution series has group-level coefficients for some oligos"
    assert all(r["t"] is not None and r["p_value"] is None for r in no_df)
    assert any(c.endswith("estimate, standard_error and t are reported.") for c in fit.caveats)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"estimator": "REML"}, "estimator must be one of reml, ml"),
        ({"estimator": "lme4"}, "estimator must be one of reml, ml"),
        ({"threads": 1.0}, "threads must be an int"),
    ],
)
def test_fit_mixed_refuses_bad_input_without_starting_a_process(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.fit_mixed("r.tsv", "d.tsv", "g.tsv", ["age"], **kwargs)
    with pytest.raises(pymzlib.UsageError, match="groups must be a path"):
        stats.fit_mixed("r.tsv", "d.tsv", None, ["age"])  # type: ignore[arg-type]
    assert bridge["calls"] == []


# ---- entrapment_fdp ------------------------------------------------------------------------------


def test_entrapment_fdp_sends_the_ratio_and_thresholds_only_when_given(bridge, tmp_path):
    bridge["payload"] = recorded("stats_entrapment_fdp_reference.json")
    t, e = tmp_path / "t.tsv", tmp_path / "e.tsv"
    stats.entrapment_fdp(t, e, ratio=1)
    assert bridge["calls"][0]["args"] == ["stats", "entrapment-fdp", "--targets", str(t),
                                          "--entrapments", str(e), "--ratio", "1.0"]
    assert bridge["calls"][0]["stdin"] is None
    assert bridge["verbs"] == [("stats entrapment-fdp", "0.6.0")]
    stats.entrapment_fdp(t, e, ratio=2.0, thresholds=[0.01, 0.05])
    assert bridge["calls"][1]["args"][-4:] == ["--ratio", "2.0", "--thresholds", "0.01,0.05"]


def test_entrapment_fdp_parses_the_reference_recording(bridge):
    bridge["payload"] = recorded("stats_entrapment_fdp_reference.json")
    fdp = stats.entrapment_fdp("t.tsv", "e.tsv", ratio=1.0)
    assert isinstance(fdp, stats.EntrapmentFdpEstimates)
    assert fdp.ratio == 1.0 and fdp.thresholds_given is False
    assert fdp.thresholds == [i / 1000 for i in range(1, 101)], "mzLib's DefaultQValueThresholds"
    assert fdp.target_candidate_count == 7 and fdp.entrapment_candidate_count == 4
    assert fdp.partnered_count == 3 and fdp.foreign_count == 0
    assert fdp.row_count == fdp.threshold_count - fdp.empty_threshold_count == 100
    row = fdp.at(0.003)
    assert (row["target_count"], row["entrapment_count"]) == (2, 1)
    assert row["lower_bound"] == pytest.approx(1 / 3) and row["paired"] == pytest.approx(1 / 3)
    with pytest.raises(pymzlib.UsageError, match="not one of the thresholds"):
        fdp.at(0.2)


def test_entrapment_fdp_paired_is_none_off_a_ratio_of_one(bridge):
    bridge["payload"] = recorded("stats_entrapment_fdp_reference_r2.json")
    fdp = stats.entrapment_fdp("t.tsv", "e.tsv", ratio=2.0, thresholds=[0.01, 0.05])
    assert fdp.thresholds_given is True and fdp.thresholds == [0.01, 0.05]
    assert fdp.columns["paired"] == [None, None]
    assert fdp.at(0.01)["combined"] == pytest.approx(3 * 1.5 / 7)
    assert fdp.caveats[-1] == "paired is null at every threshold. Here the ratio is 2 and 0 of 4 entrapments are foreign."


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"ratio": "1"}, "ratio must be a number"),
        ({"ratio": True}, "ratio must be a number"),
        ({"ratio": 1.0, "thresholds": []}, "non-empty list of q-values"),
        ({"ratio": 1.0, "thresholds": "0.01"}, "non-empty list of q-values"),
        ({"ratio": 1.0, "thresholds": [0.01, "5%"]}, r"thresholds\[1\] must be a number"),
    ],
)
def test_entrapment_fdp_refuses_bad_input_without_starting_a_process(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.entrapment_fdp("t.tsv", "e.tsv", **kwargs)
    assert bridge["calls"] == []


def test_entrapment_fdp_has_no_default_ratio():
    with pytest.raises(TypeError, match="ratio"):
        stats.entrapment_fdp("t.tsv", "e.tsv")  # type: ignore[call-arg]
