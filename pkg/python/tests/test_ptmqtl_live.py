"""ptmqtl occupancy, pairs and trait-effects through the real, packaged bridge.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drift
apart. These run the real executable on the files in ``fixtures/ptmqtl/`` (mzLib's own
hand-computed PtmQtl cases) and check that what it returns today is what was recorded, so the
docs' printed numbers are still true; that ``occupancy(out=)`` writes the table ``trait_effects``
reads; and that a refusal crosses as a usage error. Local files only.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, ptmqtl

FIXTURES = Path(__file__).parent / "fixtures"
PTMQTL = FIXTURES / "ptmqtl"


@pytest.fixture(autouse=True)
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so the ptmqtl verbs went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def assert_same_columns(live: dict, recorded: dict) -> None:
    """Equal, except that floats may differ in their last bits across platforms (1e-9 relative)."""
    assert list(live) == list(recorded)
    for name, values in recorded.items():
        assert len(live[name]) == len(values), name
        for i, (a, b) in enumerate(zip(live[name], values)):
            if isinstance(b, list):
                assert len(a) == len(b), (name, i)
                for x, y in zip(a, b):
                    assert (x is None) == (y is None) and (y is None or abs(x - y) <= 1e-9 * max(abs(y), 1e-12)), (name, i)
            elif isinstance(b, float) and isinstance(a, (int, float)) and not isinstance(a, bool):
                assert abs(a - b) <= 1e-9 * max(abs(b), 1e-300), (name, i, a, b)
            else:
                assert a == b, (name, i, a, b)


@pytest.mark.parametrize(
    "fixture, call",
    [
        ("ptmqtl_occupancy_cases.json", lambda: ptmqtl.occupancy(PTMQTL / "cases.tsv")),
        ("ptmqtl_occupancy_fractions.json", lambda: ptmqtl.occupancy(
            PTMQTL / "fractions.tsv", samples=PTMQTL / "fractions_samples.tsv")),
        ("ptmqtl_pairs_physical.json", lambda: ptmqtl.pairs(PTMQTL / "cases.tsv", type="physical")),
        ("ptmqtl_pairs_batch.json", lambda: ptmqtl.pairs(PTMQTL / "batch.tsv", type="covarying")),
        ("ptmqtl_pairs_batch_run_level.json", lambda: ptmqtl.pairs(
            PTMQTL / "batch.tsv", type="covarying", remove_run_level=True)),
        ("ptmqtl_trait_effects.json", lambda: ptmqtl.trait_effects(
            PTMQTL / "trait_occupancy.tsv", PTMQTL / "trait_runs.tsv")),
        ("ptmqtl_trait_effects_strict.json", lambda: ptmqtl.trait_effects(
            PTMQTL / "trait_occupancy.tsv", PTMQTL / "trait_runs.tsv", min_replicates=6)),
    ],
)
def test_the_recordings_are_what_the_bridge_returns_today(fixture, call):
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))
    live = call()
    assert_same_columns(live.columns, recorded["columns"])
    assert live.caveats == recorded["caveats"]
    for counts in ("state_counts", "status_counts"):
        # Keys from mzLib enums; a repin that grows one must not slip past the recording.
        if counts in recorded:
            assert getattr(live, counts) == recorded[counts]


def test_occupancy_out_writes_the_table_trait_effects_reads(tmp_path):
    out = tmp_path / "nested" / "occupancy.tsv"
    occ = ptmqtl.occupancy(PTMQTL / "trait.tsv", out=out)
    assert occ.columns is None and occ.written.row_count == occ.row_count
    assert out.read_text(encoding="utf-8") == (PTMQTL / "trait_occupancy.tsv").read_text(encoding="utf-8")
    fit = ptmqtl.trait_effects(out, PTMQTL / "trait_runs.tsv")
    recorded = json.loads((FIXTURES / "ptmqtl_trait_effects.json").read_text(encoding="utf-8"))
    assert_same_columns(fit.columns, recorded["columns"])


def test_filters_reach_mzlib(tmp_path):
    every = ptmqtl.occupancy(PTMQTL / "cases.tsv")
    biological = ptmqtl.occupancy(PTMQTL / "cases.tsv", categories=["Common Biological"])
    assert {r["category"] for r in biological.records} == {"Common Biological"}
    assert biological.categories == ["Common Biological"]
    assert biological.row_count < every.row_count
    phospho = ptmqtl.occupancy(PTMQTL / "cases.tsv", modifications=["Common Biological:Phosphorylation on S"])
    assert {r["modification"] for r in phospho.records} == {"Common Biological:Phosphorylation on S"}


def test_refusals_cross_as_usage_errors(tmp_path):
    with pytest.raises(_bridge.UsageError, match="shapes co-varying pairs only"):
        ptmqtl.pairs(PTMQTL / "batch.tsv", type="physical", remove_run_level=True)
    with pytest.raises(_bridge.UsageError, match="give --remove-run-level as well"):
        ptmqtl.pairs(PTMQTL / "batch.tsv", type="covarying", run_level_group="category")
    with pytest.raises(_bridge.UsageError, match="must end in|must name a .tsv"):
        ptmqtl.occupancy(PTMQTL / "cases.tsv", out=tmp_path / "occupancy.csv")
    repeated = tmp_path / "psms.tsv"
    lines = (PTMQTL / "cases.tsv").read_text(encoding="utf-8").splitlines()
    repeated.write_text("\n".join([lines[0], lines[1], lines[1]]) + "\n", encoding="utf-8")
    with pytest.raises(_bridge.UsageError, match="not one per PSM"):
        ptmqtl.occupancy(repeated)
    with pytest.raises(_bridge.UsageError, match="above 0 and below 0.5"):
        ptmqtl.trait_effects(PTMQTL / "trait_occupancy.tsv", PTMQTL / "trait_runs.tsv", logit_epsilon=0.5)
