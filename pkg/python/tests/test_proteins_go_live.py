"""``annotate_go`` through the real, packaged bridge, over the committed GO fixtures.

The offline tests replay a recorded payload, so they cannot see the bridge and this module drift
apart. These run the real executable on ``fixtures/proteins/`` and check the facts the offline tests
and the guide rely on. Local files only, no network: ``update_go`` is checked live by the bridge's
own ``GoLiveCanaryTests``, which skip when GO's server is down.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from pymzlib import _bridge, proteins

DATA = Path(__file__).parent / "fixtures" / "proteins"
GROUPS = DATA / "PXD036557_AllQuantifiedProteinGroups.tsv"
XML = DATA / "pxd036557_proteins.xml"
OBO = DATA / "go-pxd036557.obo"
MAP = DATA / "organelle_map.tsv"


@pytest.fixture(autouse=True)
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so annotate_go went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def test_the_recorded_payload_still_matches_the_bridge():
    live = proteins.annotate_go(GROUPS, XML, go_obo=OBO, category_map=MAP)

    assert (live.group_count, live.row_count, live.go.term_count) == (5, 563, 412)
    assert live.header["status_annotated"] == "4" and live.header["status_contaminant"] == "1"
    histones = [r for r in live.records if r["protein_group"] == "P0C0S5|Q71UI9"]
    assert sum(r["n_with"] == r["n_members"] for r in histones) == 57
    assert live.categories is not None and live.categories.row_count == 30


def test_out_holds_the_whole_table_while_the_wire_holds_none(tmp_path):
    out = tmp_path / "go.tsv"
    live = proteins.annotate_go(GROUPS, XML, go_obo=OBO, out=out, limit=0)

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "#!go_annotation_format 1"
    assert sum(not line.startswith("#!") for line in lines) == 564  # column header + 563 rows
    assert live.written == proteins.WrittenTable(str(out), 563)
    assert live.returned_count == 0 and live.records == [] and live.truncated


def test_a_wrong_extension_is_refused_before_anything_is_read(tmp_path):
    with pytest.raises(_bridge.UsageError, match=r"\.tsv"):
        proteins.annotate_go(GROUPS, XML, go_obo=OBO, out=tmp_path / "go.csv")


def test_a_missing_go_obo_points_at_update_go(tmp_path):
    with pytest.raises(_bridge.UsageError, match="update-go"):
        proteins.annotate_go(GROUPS, XML, go_obo=tmp_path / "go.obo")
    assert not (tmp_path / "go.obo").exists()
