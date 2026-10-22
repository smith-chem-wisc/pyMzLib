"""stats combine-p and stats correlate through the real, packaged bridge.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drift
apart. These run the real executable on the files in ``fixtures/stats/`` and check that what it
returns today is what was recorded (so the docs' printed numbers are still true), that it still
agrees with R from Python, and that a refusal crosses as a usage error. Local files only.
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


COHORTS = [
    ("up_everywhere", 0.003, 4, 0.8), ("up_everywhere", 0.02, 5, 0.5), ("up_everywhere", 0.3, 3, 0.3),
    ("flips_sign", 0.001, 4, 1.2), ("flips_sign", 0.004, 5, -0.9), ("flips_sign", 0.4, 3, 0.2),
    ("no_effect", 0.6, 4, -0.1), ("no_effect", 0.45, 5, 0.1), ("no_effect", 0.85, 3, -0.05),
]


def test_combine_p_matches_r_from_python():
    rows = tsv("combination_inputs.tsv")
    results = tsv("combination_results.tsv")
    fisher = stats.combine_p([(r["case"], float(r["p"])) for r in rows], method="fisher")
    plain = stats.combine_p([(r["case"], float(r["p"])) for r in rows], method="stouffer")
    weighted = stats.combine_p([(r["case"], float(r["p"]), float(r["weight"])) for r in rows],
                               method="stouffer")
    for result, column in [(fisher, "fisher_p"), (plain, "stouffer_p"), (weighted, "weighted_p")]:
        for ours, theirs in zip(result.columns["p_value"], results):
            assert abs(ours - float(theirs[column])) <= 1e-8 * float(theirs[column]), (column, theirs["case"])


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
        ("stats_combine_p_fisher.json", lambda: stats.combine_p(
            [(r["case"], float(r["p"])) for r in tsv("combination_inputs.tsv")], method="fisher")),
        ("stats_combine_p_weighted.json", lambda: stats.combine_p(
            [(r["case"], float(r["p"]), float(r["weight"])) for r in tsv("combination_inputs.tsv")],
            method="stouffer")),
        ("stats_combine_p_signed.json", lambda: stats.combine_p(COHORTS, method="stouffer")),
        ("stats_combine_p_unsigned.json", lambda: stats.combine_p(
            [c[:3] for c in COHORTS], method="stouffer")),
        ("stats_correlate_malat.json", lambda: stats.correlate(
            STATS / "malat_dilution_log2_vs_fluc.tsv", STATS / "malat_dilution_load_ng.tsv")),
        ("stats_correlate_limma_age.json", lambda: stats.correlate(
            STATS / "limma_reference_responses.tsv", STATS / "limma_reference_age.tsv")),
        ("stats_correlate_limma_age_given_sex.json", lambda: stats.correlate(
            STATS / "limma_reference_responses.tsv", STATS / "limma_reference_age.tsv",
            given=STATS / "limma_reference_sex.tsv")),
    ],
)
def test_the_recordings_are_what_the_bridge_returns_today(fixture, call):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    live = call()
    assert_same_columns(live.columns, recorded["columns"])
    assert live.caveats == recorded["caveats"]
    if "method_counts" in recorded:
        # Its keys come from an mzLib enum; a repin that grows it must not slip past the recording.
        assert live.method_counts == recorded["method_counts"]


def test_refusals_cross_as_usage_errors(tmp_path):
    with pytest.raises(_bridge.UsageError, match="every line or on none"):
        stats.combine_p([("A", 0.1, 2), ("A", 0.2)], method="stouffer")
    with pytest.raises(_bridge.UsageError, match="takes no weights"):
        stats.combine_p([("A", 0.1, 2)], method="fisher")
    with pytest.raises(_bridge.UsageError, match=r"Feature 'A'.*outside \[0, 1\]"):
        stats.combine_p([("A", 1.5)], method="fisher")
    trait = tmp_path / "trait.tsv"
    trait.write_text("sample\tage\tsex\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="has 2 value columns"):
        stats.correlate(STATS / "limma_reference_responses.tsv", trait)
