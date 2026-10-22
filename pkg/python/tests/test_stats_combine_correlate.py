"""pymzlib.stats.combine_p and correlate against recorded bridge payloads.

The numbers are mzLib's, held to R and ppcor by the bridge's C# tests and by
test_stats_combine_correlate_live.py. These tests hold the Python layer: that it sends what the
wire expects, parses what the wire returns, and refuses bad input before starting a process.
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
    state = {"calls": [], "payload": None}

    def invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": list(args), "stdin": stdin, "timeout": timeout})
        return state["payload"]

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    return state


# ---- combine_p -----------------------------------------------------------------------------------


def test_combine_p_writes_one_line_per_test_with_blank_missing_fields(bridge):
    bridge["payload"] = recorded("stats_combine_p_signed.json")
    stats.combine_p(
        [
            ("A", 0.01, 4, 0.5),
            ("A", None, 5, -1),
            {"feature": "B", "p_value": 0.2, "weight": 3, "effect": None},
            ["B", float("nan"), 2.5, 0.1],
        ],
        method="stouffer",
    )
    call = bridge["calls"][0]
    assert call["args"] == ["stats", "combine-p", "--method", "stouffer"]
    assert call["stdin"] == "A\t0.01\t4.0\t0.5\nA\t\t5.0\t-1.0\nB\t0.2\t3.0\t\nB\t\t2.5\t0.1\n"


def test_combine_p_sends_only_the_fields_each_entry_has(bridge):
    bridge["payload"] = recorded("stats_combine_p_fisher.json")
    stats.combine_p([("A", 0.1), {"feature": "B", "p_value": 0.2}], method="fisher")
    assert bridge["calls"][0]["stdin"] == "A\t0.1\nB\t0.2\n"
    # A mapping with an effect but no weight sends a blank weight, so the effect stays in field 4.
    stats.combine_p([{"feature": "A", "p_value": 0.1, "effect": -2}], method="fisher")
    assert bridge["calls"][1]["stdin"] == "A\t0.1\t\t-2.0\n"
    stats.combine_p([("A", 0.1, None)], method="stouffer")
    assert bridge["calls"][2]["stdin"] == "A\t0.1\t\n"


def test_combine_p_parses_the_fisher_recording(bridge):
    bridge["payload"] = recorded("stats_combine_p_fisher.json")
    combined = stats.combine_p([("x", 0.5)], method="fisher")
    assert isinstance(combined, stats.CombinedPValues)
    assert combined.method == "fisher" and combined.weighted is False and combined.signed is False
    assert combined.columns["feature"] == ["few", "tiny", "one", "many", "mixed", "tail"]
    assert combined.columns["studies"] == [3, 4, 1, 12, 5, 3]
    assert combined.columns["omitted"] == [0] * 6
    assert combined.line_count == 28 and combined.feature_count == combined.row_count == 6
    assert list(combined.p_values) == combined.columns["feature"]
    assert combined.p_values["one"] == pytest.approx(0.037)
    assert combined.column_names == ["feature", "studies", "omitted", "statistic", "p_value"]
    assert len(combined.caveats) == 3


def test_combine_p_parses_the_signed_recording(bridge):
    bridge["payload"] = recorded("stats_combine_p_signed.json")
    signed = stats.combine_p([("x", 0.5, 1, 1)], method="stouffer")
    assert signed.signed is True and signed.weighted is True and signed.method == "stouffer"
    assert signed.p_values["flips_sign"] > 0.4 > signed.p_values["up_everywhere"]
    assert any(c.startswith("Each p-value was read as two-sided") for c in signed.caveats)


@pytest.mark.parametrize(
    "tests, method, message",
    [
        ([("A", 0.1)], "Fisher", "method must be one of fisher, stouffer"),
        ([("A", 0.1)], None, "method must be one of"),
        ([], "fisher", "at least one test"),
        ("A\t0.1", "fisher", "a list of tests"),
        ({"feature": "A", "p_value": 0.1}, "fisher", "a list of tests"),
        ([("A",)], "fisher", r"tests\[0\] must be"),
        ([("A", 0.1, 1, 1, 1)], "stouffer", r"tests\[0\] must be"),
        ([{"feature": "A"}], "fisher", "needs 'feature' and 'p_value'"),
        ([("", 0.1)], "fisher", "non-empty feature"),
        ([("A\tB", 0.1)], "fisher", "tab"),
        ([("A\nB", 0.1)], "fisher", "newline"),
        ([("A", "0.1")], "fisher", r"tests\[0\] p_value must be a number"),
        ([("A", True)], "fisher", r"tests\[0\] p_value must be a number"),
        ([("A", 0.1, "heavy")], "stouffer", r"tests\[0\] weight must be a number"),
        ([("A", 0.1, None, "up")], "stouffer", r"tests\[0\] effect must be a number"),
    ],
)
def test_combine_p_refuses_bad_input_without_starting_a_process(bridge, tests, method, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.combine_p(tests, method=method)
    assert bridge["calls"] == []


def test_combine_p_has_no_default_method():
    with pytest.raises(TypeError, match="method"):
        stats.combine_p([("A", 0.1)])  # type: ignore[call-arg]


# ---- correlate -----------------------------------------------------------------------------------


def test_correlate_sends_the_files_and_given_only_when_set(bridge, tmp_path):
    bridge["payload"] = recorded("stats_correlate_malat.json")
    responses, trait, given = tmp_path / "r.tsv", tmp_path / "t.tsv", tmp_path / "g.tsv"
    stats.correlate(responses, trait)
    assert bridge["calls"][0]["args"] == ["stats", "correlate", "--responses", str(responses),
                                          "--trait", str(trait)]
    assert bridge["calls"][0]["stdin"] is None
    stats.correlate(responses, trait, given=given)
    assert bridge["calls"][1]["args"][-2:] == ["--given", str(given)]


def test_correlate_parses_the_malat_recording(bridge):
    bridge["payload"] = recorded("stats_correlate_malat.json")
    result = stats.correlate("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_load_ng.tsv")
    assert isinstance(result, stats.Correlations)
    assert result.trait_name == "malat_ng" and result.partial is False
    assert result.given_file is None and result.given_incomplete_count is None
    assert result.covariate_names == []
    assert set(result.method_counts) == set(stats.SPEARMAN_METHODS)
    assert sum(result.method_counts.values()) == result.feature_count == result.row_count == 492
    assert result.column_names == ["feature", "rho", "p_value", "n", "has_ties", "method"]
    for row in result.records:
        assert row["method"] in stats.SPEARMAN_METHODS
        estimable = row["method"] != "not_estimable"
        assert (row["rho"] is not None) == estimable and (row["p_value"] is not None) == estimable
        if row["method"] == "exact":
            assert row["n"] <= 9 and row["has_ties"] is False


def test_correlate_parses_the_partial_recording(bridge):
    bridge["payload"] = recorded("stats_correlate_limma_age_given_sex.json")
    result = stats.correlate("r.tsv", "t.tsv", given="g.tsv")
    assert result.partial is True and result.covariate_names == ["sex"]
    assert result.given_file.endswith("limma_reference_sex.tsv")
    assert result.given_incomplete_count == 0
    assert result.method_counts == {"exact": 0, "asymptotic": 400, "not_estimable": 0}


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"responses": 5, "trait": "t.tsv"}, "responses must be a path"),
        ({"responses": "r.tsv", "trait": b"t.tsv"}, "trait must be a path"),
        ({"responses": "r.tsv", "trait": "t.tsv", "given": ""}, "given must be a path"),
    ],
)
def test_correlate_refuses_a_non_path_without_starting_a_process(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        stats.correlate(**kwargs)
    assert bridge["calls"] == []
