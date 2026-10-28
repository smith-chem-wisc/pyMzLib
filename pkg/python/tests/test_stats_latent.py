"""pymzlib.stats.fit_latent_mixed and fit_link against recorded bridge payloads.

The numbers are mzLib's LogitNormalMixedModel, held to lme4's lmer and glmer and to glmmTMB by the
bridge's C# tests and by test_stats_latent_live.py. These tests hold the Python layer: that it sends
what the wire expects (and nothing for an option left out, so mzLib's default applies), parses what
the wire returns, and refuses bad input before starting a process. The replay bridge's matching of
these recordings is held here too, since every doctest of the two calls relies on it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, stats

import replay_bridge

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


# ---- fit_latent_mixed ----------------------------------------------------------------------------


def test_fit_latent_mixed_sends_only_what_was_given(bridge, tmp_path):
    bridge["payload"] = recorded("stats_fit_latent_mixed_scatter.json")
    d, g, c = tmp_path / "d.tsv", tmp_path / "g.tsv", tmp_path / "c.tsv"
    stats.fit_latent_mixed(d, g, ["trait"], link="identity", counts=c)
    call = bridge["calls"][0]
    assert call["args"] == ["stats", "fit-latent-mixed", "--design", str(d), "--groups", str(g),
                            "--link", "identity", "--threads", "1", "--counts", str(c)]
    assert call["stdin"] == "trait\n"
    assert bridge["verbs"] == [("stats fit-latent-mixed", "0.7.0")]

    r, depth, sigma = tmp_path / "r.tsv", tmp_path / "depth.tsv", tmp_path / "s.tsv"
    stats.fit_latent_mixed(d, g, ["trait", "x"], link=(-0.4, 0.6, 0.3), responses=r, counts=c, depth=depth,
                           residual_variance=sigma, estimator="ml", threads=-1)
    assert bridge["calls"][1]["args"][6:] == [
        "--link", "-0.4,0.6,0.3", "--threads", "-1", "--responses", str(r), "--counts", str(c),
        "--depth", str(depth), "--residual-variance", str(sigma), "--estimator", "ml"]


def test_a_link_sends_only_the_numbers_it_holds(bridge):
    bridge["payload"] = recorded("stats_fit_latent_mixed_scatter.json")

    def sent(link):
        stats.fit_latent_mixed("d.tsv", "g.tsv", ["trait"], link=link, counts="c.tsv")
        args = bridge["calls"][-1]["args"]
        return args[args.index("--link") + 1]

    assert sent("IDENTITY") == "identity"
    assert sent((0.4, 0.6)) == "0.4,0.6", "a tuple sends what it holds: the rest is mzLib's default"
    assert sent([1, 2, 0.5, 0.1]) == "1.0,2.0,0.5,0.1"
    assert sent(stats.DetectionLink(0.1, 0.9, 0.0, 0.0)) == "0.1,0.9,0.0,0.0"
    fitted = stats.LinkFit._from_wire(recorded("stats_fit_link_simulated.json"))
    assert sent(fitted) == ",".join(repr(v) for v in (fitted.link.intercept, fitted.link.slope,
                                                     fitted.link.depth_slope, fitted.link.overdispersion))


def test_fit_latent_mixed_parses_the_recordings(bridge):
    bridge["payload"] = recorded("stats_fit_latent_mixed_simulated.json")
    fit = stats.fit_latent_mixed("d.tsv", "g.tsv", ["Intercept", "trait"], link=(0.4, 0.6),
                                 responses="r.tsv", counts="c.tsv")
    assert isinstance(fit, stats.LatentMixedModelFit)
    assert fit.estimator == "reml", "absent, mzLib's default, as mzLib reports it"
    assert fit.link == stats.DetectionLink(0.4, 0.6, 0.0, 0.0)
    assert fit.grouping == "animal" and fit.group_count == 16
    assert (fit.measured_cells, fit.counted_cells, fit.missing_count) == (96, 96, 0)
    assert fit.depth_file is None and fit.residual_variance_file is None
    assert fit.responses_file is not None and fit.counts_file is not None
    assert list(fit.status_counts) == list(stats.LATENT_MIXED_STATUSES)
    assert fit.column_names == ["feature", "coefficient", "status", "estimate", "standard_error", "df", "t",
                                "p_value", "residual_variance", "residual_variance_given", "group_variance",
                                "log_likelihood", "measured", "counted", "groups"]
    for row in fit.rows("trait"):
        assert (row["measured"], row["counted"], row["df"]) == (32, 32, 14)
        assert row["residual_variance_given"] is False

    bridge["payload"] = recorded("stats_fit_latent_mixed_link_ml.json")
    ml = stats.fit_latent_mixed("d.tsv", "g.tsv", ["Intercept", "trait", "x"], link=(-0.4, 0.6, 0.3),
                                counts="c.tsv", depth="dp.tsv", residual_variance="s.tsv", estimator="ml")
    assert ml.estimator == "ml" and ml.responses_file is None and ml.measured_cells == 0
    assert ml.link == stats.DetectionLink(-0.4, 0.6, 0.3, 0.0)
    assert all(r["residual_variance_given"] is True for r in ml.records)
    assert any(c.startswith("16 of 16 features were fitted at the residual_variance you gave") for c in ml.caveats)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"link": "logit"}, "link takes 'identity'"),
        ({"link": (1.0,)}, "2 to 4 numbers"),
        ({"link": (1, 2, 3, 4, 5)}, "2 to 4 numbers"),
        ({"link": (0.0, "one")}, r"link\[1\] must be a number"),
        ({"link": None}, "link takes"),
        ({"link": "identity", "estimator": "REML"}, "estimator must be one of reml, ml"),
        ({"link": "identity", "threads": 1.0}, "threads must be an int"),
        ({"link": "identity", "depth": 3}, "depth must be a path"),
    ],
)
def test_fit_latent_mixed_refuses_bad_input_without_starting_a_process(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.fit_latent_mixed("d.tsv", "g.tsv", ["trait"], counts="c.tsv", **kwargs)
    assert bridge["calls"] == []


def test_fit_latent_mixed_needs_responses_or_counts(bridge):
    with pytest.raises(pymzlib.UsageError, match="needs responses= .* counts= .* or both"):
        stats.fit_latent_mixed("d.tsv", "g.tsv", ["trait"], link="identity")
    unconverged = stats.LinkFit("c.tsv", 1, False, False, False, True, None)
    with pytest.raises(pymzlib.UsageError, match="did not converge"):
        stats.fit_latent_mixed("d.tsv", "g.tsv", ["trait"], link=unconverged, counts="c.tsv")
    assert bridge["calls"] == []


# ---- fit_link ------------------------------------------------------------------------------------


def test_fit_link_sends_each_option_only_when_given(bridge, tmp_path):
    bridge["payload"] = recorded("stats_fit_link_simulated.json")
    cells = tmp_path / "cells.tsv"
    stats.fit_link(cells)
    assert bridge["calls"][0]["args"] == ["stats", "fit-link", "--cells", str(cells)]
    assert bridge["verbs"] == [("stats fit-link", "0.7.0")]
    stats.fit_link(cells, with_depth=True, with_overdispersion=False, detection_only=False)
    assert bridge["calls"][1]["args"][4:] == ["--with-depth", "true", "--with-overdispersion", "false",
                                             "--detection-only", "false"]


def test_fit_link_parses_the_recording(bridge):
    bridge["payload"] = recorded("stats_fit_link_simulated.json")
    fit = stats.fit_link("link_cells.tsv", with_depth=True)
    assert isinstance(fit, stats.LinkFit)
    assert (fit.cell_count, fit.depth_given, fit.with_depth) == (4000, True, True)
    assert (fit.with_overdispersion, fit.detection_only) == (False, True), "mzLib's defaults, reported"
    assert isinstance(fit.link, stats.DetectionLink)
    assert abs(fit.link.intercept + 0.5) < 0.15 and abs(fit.link.slope - 0.5) < 0.1
    assert fit.link.overdispersion == 0.0
    assert any(c.startswith("overdispersion was not fitted") for c in fit.caveats)
    bridge["payload"] = dict(recorded("stats_fit_link_simulated.json"), link=None)
    assert stats.fit_link("link_cells.tsv").link is None


@pytest.mark.parametrize("kwargs", [{"with_depth": 1}, {"detection_only": "false"}, {"with_overdispersion": 0}])
def test_fit_link_refuses_a_non_bool_option_without_starting_a_process(bridge, kwargs):
    with pytest.raises(pymzlib.UsageError, match="must be True, False or None"):
        stats.fit_link("cells.tsv", **kwargs)
    assert bridge["calls"] == []


# ---- the replay bridge, which the doctests of both calls run on ----------------------------------


def test_replay_matches_a_link_by_the_numbers_given():
    data = recorded("stats_fit_latent_mixed_link_ml.json")
    files = {"counts": "lnmm_link_counts.tsv", "depth": "lnmm_link_depth.tsv",
             "residual-variance": "lnmm_link_residual_variance.tsv"}
    assert replay_bridge.mismatch(data, {**files, "link": "-0.4,0.6,0.3"}) == ""
    assert replay_bridge.mismatch(data, {**files, "link": "-0.4,0.6"}) == ""
    assert "not --link" in replay_bridge.mismatch(data, {**files, "link": "-0.4,0.6,0.2"})
    assert "not --link" in replay_bridge.mismatch(data, {**files, "link": "identity"})
    scatter = recorded("stats_fit_latent_mixed_scatter.json")
    assert replay_bridge.mismatch(scatter, {"counts": "lnmm_scatter_counts.tsv", "link": "identity"}) == ""
    assert replay_bridge.mismatch(scatter, {"counts": "lnmm_scatter_counts.tsv", "link": "0.0,1.0"}) == ""


def test_replay_holds_the_optional_files_and_a_true_false_option():
    scatter = recorded("stats_fit_latent_mixed_scatter.json")
    assert "recorded without --responses" in replay_bridge.mismatch(scatter, {"responses": "r.tsv"})
    assert "recorded with --counts" in replay_bridge.mismatch(scatter, {"link": "identity"})
    assert "recorded without --depth" in replay_bridge.mismatch(
        scatter, {"counts": "lnmm_scatter_counts.tsv", "depth": "d.tsv"})
    link = recorded("stats_fit_link_simulated.json")
    assert replay_bridge.mismatch(link, {"cells": "link_cells.tsv", "with-depth": "true"}) == ""
    assert "with_depth=True" in replay_bridge.mismatch(link, {"cells": "link_cells.tsv", "with-depth": "false"})
