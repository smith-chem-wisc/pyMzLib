"""pymzlib.ptmqtl against recorded bridge payloads.

The numbers are mzLib's, held to mzLib's own hand-computed cases by the bridge's C# tests
(pkg/bridge.tests/PtmQtlTests.cs) and by test_ptmqtl_live.py. These tests hold the Python layer:
that it sends what the wire expects (and nothing when the caller gives nothing, so mzLib's defaults
apply), parses what the wire returns, and refuses bad input before starting a process.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, ptmqtl

FIXTURES = Path(__file__).parent / "fixtures"


def recorded(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture
def bridge(monkeypatch):
    """Answer every bridge call with ``state["payload"]``, recording the arguments."""
    state = {"calls": [], "payload": None}

    def invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": list(args), "stdin": stdin, "timeout": timeout})
        return state["payload"]

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    return state


def test_module_is_exported():
    assert pymzlib.ptmqtl is ptmqtl
    assert "ptmqtl" in pymzlib.__all__
    assert set(dir(ptmqtl)) == set(ptmqtl.__all__)


# ---- occupancy -----------------------------------------------------------------------------------


def test_occupancy_sends_only_what_was_given(bridge):
    bridge["payload"] = recorded("ptmqtl_occupancy_cases.json")
    ptmqtl.occupancy("cases.tsv")
    assert bridge["calls"][0]["args"] == ["ptmqtl", "occupancy", "--observations", "cases.tsv"]
    ptmqtl.occupancy(Path("cases.tsv"), samples="s.tsv", categories=["Common Biological", "UniProt"], out="o.tsv")
    assert bridge["calls"][1]["args"] == [
        "ptmqtl", "occupancy", "--observations", "cases.tsv", "--samples", "s.tsv",
        "--categories", "Common Biological\tUniProt", "--out", "o.tsv",
    ]
    ptmqtl.occupancy("cases.tsv", modifications=["Common Biological:Phosphorylation on S"])
    assert bridge["calls"][2]["args"][-2:] == ["--modifications", "Common Biological:Phosphorylation on S"]


def test_occupancy_parses_the_cases_recording(bridge):
    bridge["payload"] = recorded("ptmqtl_occupancy_cases.json")
    occ = ptmqtl.occupancy("cases.tsv")
    assert isinstance(occ, ptmqtl.SiteOccupancy)
    assert occ.column_names == [
        "site_key", "protein_accession", "position", "residue", "modification", "category", "run",
        "state", "modified_intensity", "covering_intensity", "fraction", "unmodified_quantified",
        "covering_peptidoforms", "modified_peptidoforms",
    ]
    assert set(occ.state_counts) == set(ptmqtl.OCCUPANCY_STATES)
    assert occ.row_count == len(occ.records) == sum(occ.state_counts.values())
    assert occ.samples_file is None and occ.sample_count is None and occ.written is None
    assert occ.modifications is None and occ.categories is None and occ.ignored_columns == []
    by_key = {(r["site_key"], r["run"]): r for r in occ.records}
    ratio = by_key[("RATIO:S13:Common Biological:Phosphorylation on S", "r1")]
    assert ratio["fraction"] == pytest.approx(30 / 200, abs=1e-15)
    assert ratio["modified_intensity"] == 30 and ratio["covering_intensity"] == 200
    # A censored cell carries no fraction, never a 0.
    for record in occ.records:
        assert (record["fraction"] is None) == (record["state"] != "quantified")
    # mzLib's exclusions: no site from these proteins.
    proteins = {r["protein_accession"] for r in occ.records}
    assert proteins.isdisjoint({"EXCLUDED", "PEPTIDE_NTERM", "NTERM_AFTER_K"})
    assert by_key[("NTERM_AFTER_MET:S2:UniProt:N-acetylserine on S", "r1")]["position"] == 2
    assert any(c.startswith("position is 1-based") for c in occ.caveats)


def test_occupancy_carries_mzlibs_peptidoform_counts(bridge):
    # mzLib #1468: identifications, not intensity. The k of n of each cell.
    bridge["payload"] = recorded("ptmqtl_occupancy_cases.json")
    occ = ptmqtl.occupancy("cases.tsv")
    states = {r["run"]: (r["state"], r["covering_peptidoforms"], r["modified_peptidoforms"])
              for r in occ.records if r["protein_accession"] == "STATES"}
    assert states == {
        "countOnly": ("count_only", 1, 1),
        "floor": ("floor", 2, 1),
        "notDetected": ("not_detected", 1, 0),
        "quant": ("quantified", 2, 1),
    }
    ratio = next(r for r in occ.records if r["protein_accession"] == "RATIO")
    assert (ratio["covering_peptidoforms"], ratio["modified_peptidoforms"]) == (3, 1)
    assert any(c.startswith("covering_peptidoforms and modified_peptidoforms") for c in occ.caveats)


def test_occupancy_with_samples_reports_the_samples(bridge):
    bridge["payload"] = recorded("ptmqtl_occupancy_fractions.json")
    occ = ptmqtl.occupancy("fractions.tsv", samples="fractions_samples.tsv")
    assert occ.run_names == ["A", "B"]
    assert occ.run_count == 3 and occ.sample_count == 2
    assert occ.samples_file.endswith("fractions_samples.tsv")
    assert any(c.startswith("With --samples") for c in occ.caveats)


def test_occupancy_written_to_out_has_no_columns(bridge):
    bridge["payload"] = recorded("ptmqtl_occupancy_trait_out.json")
    occ = ptmqtl.occupancy("trait.tsv", out="trait_occupancy.tsv")
    assert occ.columns is None and occ.records == []
    assert isinstance(occ.written, pymzlib.readers.WrittenTable)
    assert occ.written.format == "tsv" and occ.written.row_count == occ.row_count == 24
    assert occ.written.path.endswith("trait_occupancy.tsv")


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"modifications": ["a"], "categories": ["b"]}, "not both"),
        ({"categories": "Common Biological"}, "takes a list of names"),
        ({"categories": []}, "non-empty list"),
        ({"categories": ["ok", ""]}, "non-empty str"),
        ({"modifications": ["a\tb"]}, "a tab or a line break"),
        ({"samples": 3}, "samples must be a path"),
        ({"out": b"x.tsv"}, "out must be a path"),
    ],
)
def test_occupancy_refuses_bad_arguments_before_any_process(bridge, kwargs, message):
    with pytest.raises(_bridge.UsageError, match=message):
        ptmqtl.occupancy("cases.tsv", **kwargs)
    assert bridge["calls"] == []


# ---- pairs ---------------------------------------------------------------------------------------


def test_pairs_sends_only_what_was_given(bridge):
    bridge["payload"] = recorded("ptmqtl_pairs_batch.json")
    ptmqtl.pairs("batch.tsv", type="covarying")
    assert bridge["calls"][0]["args"] == ["ptmqtl", "pairs", "--observations", "batch.tsv", "--type", "covarying"]
    ptmqtl.pairs(
        "batch.tsv", type="both", samples="s.tsv", min_quantified_fraction=0.5, keep_ceiling=True,
        remove_run_level=True, run_level_group="category",
    )
    assert bridge["calls"][1]["args"] == [
        "ptmqtl", "pairs", "--observations", "batch.tsv", "--type", "both", "--samples", "s.tsv",
        "--min-quantified-fraction", "0.5", "--keep-ceiling", "--remove-run-level",
        "--run-level-group", "category",
    ]
    ptmqtl.pairs("batch.tsv", type="covarying", remove_run_level=True, run_level_groups="groups.tsv")
    assert bridge["calls"][2]["args"][-2:] == ["--run-level-groups", "groups.tsv"]
    ptmqtl.pairs("depth.tsv", type="covarying", include_bounds=True, run_covariate=Path("psms.tsv"))
    assert bridge["calls"][3]["args"] == [
        "ptmqtl", "pairs", "--observations", "depth.tsv", "--type", "covarying",
        "--include-bounds", "--run-covariate", "psms.tsv",
    ]
    # Left out, neither is sent: mzLib's defaults (no bounds, no covariate) apply.
    assert "--include-bounds" not in bridge["calls"][0]["args"]
    assert "--run-covariate" not in bridge["calls"][0]["args"]


def test_pairs_parses_the_physical_recording(bridge):
    bridge["payload"] = recorded("ptmqtl_pairs_physical.json")
    found = ptmqtl.pairs("cases.tsv", type="physical")
    assert isinstance(found, ptmqtl.PtmPairs)
    assert found.type == "physical" and found.physical_count == 1 and found.covarying_count is None
    # Co-varying settings are not echoed for physical pairs: none applied.
    assert found.min_quantified_fraction is None and found.exclude_ceiling is None
    assert found.remove_run_level is None and found.covarying_run_count is None
    (pair,) = found.records
    assert pair["result_type"] == "physical" and pair["overlapping"] is True
    assert pair["statistic"] == pytest.approx(0.25, abs=1e-15)
    assert (pair["n"], pair["statistic_n"]) == (2, 1)
    assert pair["p_value"] is None and pair["q"] is None and pair["spearman_method"] is None
    assert pair["fdr_family"] == "P:intra"
    # The co-varying test fields: no test for physical pairs.
    assert pair["test_status"] is None
    assert pair["bounds_included"] is False and pair["run_covariate_removed"] is False
    assert found.include_bounds is None and found.min_untied_share is None
    assert found.run_covariate_file is None and found.run_covariate_runs_missing is None


def test_pairs_parses_mzlibs_defaults_and_the_run_level(bridge):
    bridge["payload"] = recorded("ptmqtl_pairs_batch.json")
    plain = ptmqtl.pairs("batch.tsv", type="covarying")
    # mzLib's excludeCeiling is nullable since #1468 and resolved inside CoVarying: absent stays None.
    assert plain.min_quantified_fraction == 0.7 and plain.exclude_ceiling is None
    assert plain.remove_run_level is False and plain.run_level_group is None
    assert plain.covarying_run_count == 12 and plain.physical_count is None
    assert plain.min_runs_for_run_level is None

    bridge["payload"] = recorded("ptmqtl_pairs_batch_run_level.json")
    given = ptmqtl.pairs("batch.tsv", type="covarying", remove_run_level=True)
    assert given.remove_run_level is True and given.run_level_group == "none"
    assert (given.min_runs_for_run_level, given.min_sites_for_run_level) == (10, 3)
    assert all(r["run_level_removed"] for r in given.records)
    not_estimable = [r for r in given.records if r["spearman_method"] == "not_estimable"]
    assert not_estimable and all(r["statistic"] is None and r["q"] is None for r in not_estimable)
    assert {r["test_status"] for r in not_estimable} == {"not_estimable"}
    assert any(c.startswith("With --remove-run-level") for c in given.caveats)
    assert any(c.startswith(f"{len(not_estimable)} co-varying pairs were not tested") for c in given.caveats)


def test_pairs_with_bounds_tests_what_the_guard_kept_out(bridge):
    # mzLib's BoundsEnterTypeAAsTiedRanks: six values in ten runs is below the 70% guard.
    bridge["payload"] = recorded("ptmqtl_pairs_bounds_plain.json")
    plain = ptmqtl.pairs("bounds.tsv", type="covarying")
    assert plain.row_count == 0 and plain.include_bounds is False and plain.min_untied_share is None

    bridge["payload"] = recorded("ptmqtl_pairs_bounds.json")
    bounded = ptmqtl.pairs("bounds.tsv", type="covarying", include_bounds=True)
    assert bounded.include_bounds is True and bounded.min_untied_share == 0.5
    assert bounded.exclude_ceiling is None
    (pair,) = bounded.records
    assert pair["n"] == 9 and pair["test_status"] == "tested" and pair["bounds_included"] is True
    assert pair["statistic"] is not None and pair["run_covariate_removed"] is False
    assert any(c.startswith("With --include-bounds, a pair is tested only if") for c in bounded.caveats)


def test_pairs_given_a_run_covariate(bridge):
    # mzLib's ARunCovariateIsPartialledOut: two sites that rise with each run's depth and a batch.
    def pair(result):
        return next(r for r in result.records
                    if r["site_a"].startswith("DEPTH_S") and r["site_b"].startswith("DEPTH_T"))

    bridge["payload"] = recorded("ptmqtl_pairs_depth.json")
    plain = ptmqtl.pairs("depth.tsv", type="covarying")
    assert plain.run_covariate_file is None and plain.run_covariate_name is None
    assert pair(plain)["run_covariate_removed"] is False

    bridge["payload"] = recorded("ptmqtl_pairs_depth_covariate.json")
    given = ptmqtl.pairs("depth.tsv", type="covarying", run_covariate="depth_psms.tsv")
    assert given.run_covariate_file.endswith("depth_psms.tsv")
    assert given.run_covariate_name == "target_psms" and given.run_covariate_runs_missing == []
    assert given.min_runs_for_run_level == 10 and given.min_sites_for_run_level is None
    assert pair(given)["run_covariate_removed"] is True and pair(given)["run_level_removed"] is False
    assert pair(given)["statistic"] != pair(plain)["statistic"]

    bridge["payload"] = recorded("ptmqtl_pairs_depth_covariate_run_level.json")
    both = ptmqtl.pairs("depth.tsv", type="covarying", run_covariate="depth_psms.tsv", remove_run_level=True)
    assert pair(both)["run_covariate_removed"] is True and pair(both)["run_level_removed"] is True
    assert abs(pair(both)["statistic"]) < 0.1


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"type": "all"}, "type must be one of"),
        ({"type": None}, "type must be one of"),
        ({"type": "covarying", "min_quantified_fraction": "0.7"}, "must be a number"),
        ({"type": "covarying", "min_quantified_fraction": True}, "must be a number"),
        ({"type": "covarying", "keep_ceiling": "yes"}, "keep_ceiling must be True or False"),
        ({"type": "covarying", "remove_run_level": 1}, "remove_run_level must be True or False"),
        ({"type": "covarying", "include_bounds": "yes"}, "include_bounds must be True or False"),
        ({"type": "covarying", "run_covariate": 3}, "run_covariate must be a path"),
        ({"type": "covarying", "remove_run_level": True, "run_level_group": "chemistry"}, "run_level_group must be one of"),
        ({"type": "covarying", "categories": ["a"], "modifications": ["b"]}, "not both"),
    ],
)
def test_pairs_refuses_bad_arguments_before_any_process(bridge, kwargs, message):
    with pytest.raises(_bridge.UsageError, match=message):
        ptmqtl.pairs("batch.tsv", **kwargs)
    assert bridge["calls"] == []


def test_pairs_requires_a_type():
    with pytest.raises(TypeError):
        ptmqtl.pairs("batch.tsv")  # type: ignore[call-arg]


# ---- trait_effects -------------------------------------------------------------------------------


def test_trait_effects_sends_only_what_was_given(bridge):
    bridge["payload"] = recorded("ptmqtl_trait_effects.json")
    ptmqtl.trait_effects("occ.tsv", "runs.tsv")
    assert bridge["calls"][0]["args"] == ["ptmqtl", "trait-effects", "--occupancy", "occ.tsv", "--runs", "runs.tsv"]
    ptmqtl.trait_effects("occ.tsv", "runs.tsv", min_runs=8, min_replicates=5, logit_epsilon=1e-6, keep_ceiling=True,
                         scale="log2_odds")
    assert bridge["calls"][1]["args"] == [
        "ptmqtl", "trait-effects", "--occupancy", "occ.tsv", "--runs", "runs.tsv",
        "--min-runs", "8", "--min-replicates", "5", "--logit-epsilon", "1e-06", "--keep-ceiling",
        "--scale", "log2_odds",
    ]
    assert "--scale" not in bridge["calls"][0]["args"]


def test_trait_effects_parses_the_recording(bridge):
    bridge["payload"] = recorded("ptmqtl_trait_effects.json")
    fit = ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv")
    assert isinstance(fit, ptmqtl.TraitEffects)
    assert (fit.min_runs, fit.min_replicates, fit.logit_epsilon, fit.exclude_ceiling) == (6, 4, 5e-05, True)
    assert fit.scale == "natural_log_odds" and fit.scale in ptmqtl.OCCUPANCY_SCALES
    assert list(fit.status_counts) == list(ptmqtl.TRAIT_STATUSES)
    assert fit.trait_name == "dose" and fit.covariate_names == ["injection"]
    assert fit.trait_run_count == fit.run_count == 12 and fit.occupancy_runs_not_used == []
    trait = next(r for r in fit.records if r["protein_accession"] == "TRAIT_SITE")
    assert trait["status"] == "fitted" and trait["effect"] == pytest.approx(0.5, abs=0.1)
    assert trait["df"] == 4 and (trait["runs"], trait["replicates"]) == (12, 6)
    assert len(trait["covariate_effects"]) == 1
    assert trait["covariate_effects"][0] == pytest.approx(0.02, abs=0.005)


def test_trait_effects_below_support_has_null_statistics(bridge):
    bridge["payload"] = recorded("ptmqtl_trait_effects_strict.json")
    fit = ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv", min_replicates=6)
    assert fit.min_replicates == 6 and fit.status_counts["below_support"] == 1
    censored = next(r for r in fit.records if r["protein_accession"] == "CENSORED_SITE")
    assert censored["status"] == "below_support"
    assert censored["effect"] is None and censored["p_value"] is None and censored["q"] is None
    assert censored["covariate_effects"] == [None]
    assert censored["median_fraction"] is not None and censored["replicates"] == 5


def test_trait_effects_on_log2_odds_is_the_natural_fit_over_ln2(bridge):
    # mzLib's TraitFitOnLog2OddsIsTheNaturalLogFitDividedByLn2: a change of unit, not a refit.
    import math

    bridge["payload"] = recorded("ptmqtl_trait_effects.json")
    natural = {r["site_key"]: r for r in ptmqtl.trait_effects("o.tsv", "r.tsv").records}
    bridge["payload"] = recorded("ptmqtl_trait_effects_log2.json")
    log2 = ptmqtl.trait_effects("o.tsv", "r.tsv", scale="log2_odds")
    assert log2.scale == "log2_odds"
    assert any(c.startswith("scale is log2_odds") for c in log2.caveats)
    for r in log2.records:
        n = natural[r["site_key"]]
        assert r["effect"] == pytest.approx(n["effect"] / math.log(2), rel=1e-12)
        assert r["standard_error"] == pytest.approx(n["standard_error"] / math.log(2), rel=1e-12)
        assert r["replicate_variance"] == pytest.approx(n["replicate_variance"] / math.log(2) ** 2, rel=1e-12)
        assert (r["df"], r["p_value"], r["q"]) == (n["df"], n["p_value"], n["q"])


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"min_runs": 6.0}, "min_runs must be an int"),
        ({"min_replicates": True}, "min_replicates must be an int"),
        ({"logit_epsilon": "small"}, "logit_epsilon must be a number"),
        ({"keep_ceiling": None}, "keep_ceiling must be True or False"),
        ({"scale": "log10_odds"}, "scale must be one of natural_log_odds, log2_odds"),
        ({"scale": 2}, "scale must be one of"),
    ],
)
def test_trait_effects_refuses_bad_arguments_before_any_process(bridge, kwargs, message):
    with pytest.raises(_bridge.UsageError, match=message):
        ptmqtl.trait_effects("occ.tsv", "runs.tsv", **kwargs)
    assert bridge["calls"] == []


def test_trait_effects_needs_both_paths(bridge):
    with pytest.raises(_bridge.UsageError, match="runs must be a path"):
        ptmqtl.trait_effects("occ.tsv", None)  # type: ignore[arg-type]
    assert bridge["calls"] == []
