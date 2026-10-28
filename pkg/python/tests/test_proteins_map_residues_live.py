"""``genes map-residues`` through the real, packaged bridge, over the committed cut of release 116.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drifting
apart. These run the real executable on ``fixtures/map_residues/`` (a verbatim cut of Ensembl release
116; see PROVENANCE.txt), re-run each recording's call and compare, and check the refusals the
recordings cannot show. Local files only, no network.
"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, proteins

FIXTURES = Path(__file__).parent / "fixtures"
INPUTS = FIXTURES / "map_residues"
GTFS = {
    "homo_sapiens": INPUTS / "Homo_sapiens.GRCh38.116.gtf",
    "mus_musculus": INPUTS / "Mus_musculus.GRCm39.116.gtf",
    "rattus_norvegicus": INPUTS / "Rattus_norvegicus.GRCr8.116.gtf",
}
ALIGNMENT = INPUTS / "Compara.116.protein_default.aa.fasta.gz"
TREES = INPUTS / "vertebrates.GeneTree_content.default.e116.txt.gz"
TP53, TRP53 = "ENSP00000269305", "ENSMUSP00000104298"


@pytest.fixture(autouse=True)
def built_bridge():
    # Same rule as test_proteins_live: skip in a source checkout with no built bridge, but fail under
    # CI, where the bridge is always built and a skip would hide the drift these tests exist for.
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so genes map-residues went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def _recorded(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def _queries(recording: dict) -> list[tuple]:
    c = recording["columns"]
    return list(zip(c["source_protein_id"], c["position"], c["target_protein_id"]))


def _map(queries, **kwargs):
    arguments = {"alignment": ALIGNMENT, "gene_trees": TREES, "gtfs": GTFS}
    arguments.update(kwargs)
    return proteins.map_residues(queries, **arguments)


@pytest.mark.parametrize(
    "fixture, with_proteins",
    [("genes_map_residues_tp53.json", True), ("genes_map_residues_chain.json", False)],
)
def test_each_recording_is_what_the_bridge_answers_now(fixture, with_proteins):
    recorded = _recorded(fixture)
    live = _map(_queries(recorded), proteins=with_proteins)

    assert live.columns == recorded["columns"]
    assert live.outcome_counts == recorded["outcome_counts"]
    assert live.alignment == proteins.ComparaAlignment._from_wire(recorded["alignment"]), "byte-exact inputs (-text)"
    assert live.gene_trees == proteins.GeneTrees._from_wire(recorded["gene_trees"])
    assert live.caveats == recorded["caveats"]
    assert (live.proteins.columns if live.proteins else None) == (
        recorded["proteins"]["columns"] if recorded["proteins"] else None
    )


def test_a_position_past_the_source_fails_the_whole_call_with_mzlibs_message():
    # mzLib's MapResidue throws for it (the bridge checks no length), so the first query's answer
    # is lost too: no row comes back.
    with pytest.raises(pymzlib.UsageError) as caught:
        _map([(TP53, 15, TRP53), (TP53, 394, TRP53)])
    message = str(caught.value)
    assert "Query 2 (ENSP00000269305, position 394, ENSMUSP00000104298)" in message
    assert "ENSP00000269305 has 393 residues." in message, "mzLib's own words"
    assert "no query of this call was answered" in message


def test_position_zero_is_mzlibs_refusal_too():
    with pytest.raises(pymzlib.UsageError, match="has 393 residues"):
        _map([(TP53, 0, TRP53)])


def test_a_source_not_kept_is_an_outcome_before_its_position_is_looked_at():
    # MapResidue checks the source before the position, so an unkept source with any position is
    # SourceNotInAlignment, not a refusal.
    mapped = _map([("ENSP00000352610", 100000, TRP53)])
    assert mapped.columns["outcome"] == ["SourceNotInAlignment"]


def test_without_gene_sets_mzlibs_default_keeps_every_gene_and_a_caveat_says_so():
    mapped = _map([(TP53, 15, TRP53)], gtfs=None)

    assert mapped.gene_sets == {}
    assert mapped.gene_trees.restricted_to_gene_sets is None
    assert mapped.columns["target_position"] == [18], "the same answer for the same proteins"
    assert any("No gene sets were given, so mzLib kept every gene" in c for c in mapped.caveats)


def test_a_ragged_alignment_is_mzlibs_correctness_failure(tmp_path):
    bad = tmp_path / "Compara.116.protein_default.aa.fasta.gz"
    with gzip.open(bad, "wt", encoding="utf-8") as f:
        f.write(">ENSDARP00000000001\nPEPT\n>ENSDARP00000000002\nPEP\n//\n")

    with pytest.raises(pymzlib.BridgeError) as caught:
        _map([(TP53, 15, TRP53)], alignment=bad)
    assert caught.value.error_type == "InvalidDataException"
    assert "is 3 columns wide, the alignment 4" in str(caught.value)


def test_a_missing_alignment_is_a_usage_error_before_any_input_is_read(tmp_path):
    with pytest.raises(pymzlib.UsageError, match="Compara gene-tree alignment not found"):
        _map([(TP53, 15, TRP53)], alignment=tmp_path / "absent.aa.fasta.gz")
