"""Tests for :func:`pymzlib.proteins.map_residues`, the projection of ``genes map-residues``.

Offline. The fixtures are recorded from the real bridge (scripts/record_fixture.py) over a real cut
of Ensembl release 116 in ``fixtures/map_residues/`` (every line verbatim; see PROVENANCE.txt):

* ``genes_map_residues_tp53.json``: human, mouse and rat gene sets, nine queries through the TP53
  and GAPDH gene trees' alignments, one per outcome and more, with ``--proteins``.
* ``genes_map_residues_chain.json``: the same inputs, one query (human p53 Ser15 to mouse), without
  ``--proteins``.

What is pinned here is what pyMzLib sends (the species lines, then ``--``, then one query per line,
and only the options it was given), what it refuses before spawning a process, and that it projects
mzLib's answer without repairing it: an outcome is mzLib's own name, a gap stays ``None``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymzlib import _bridge, proteins

FIXTURES = Path(__file__).parent / "fixtures"
HUMAN, MOUSE, RAT = "homo_sapiens", "mus_musculus", "rattus_norvegicus"
GTFS = {
    HUMAN: "Homo_sapiens.GRCh38.116.gtf",
    MOUSE: "Mus_musculus.GRCm39.116.gtf",
    RAT: "Rattus_norvegicus.GRCr8.116.gtf",
}
ALIGNMENT = "Compara.116.protein_default.aa.fasta.gz"
TREES = "vertebrates.GeneTree_content.default.e116.txt.gz"
TP53, TRP53, RAT_TP53 = "ENSP00000269305", "ENSMUSP00000104298", "ENSRNOP00000074031"
QUERIES = [
    (TP53, 15, TRP53), (TP53, 20, TRP53), (TP53, 46, TRP53), (TP53, 392, TRP53), (TRP53, 18, TP53),
    (TP53, 15, RAT_TP53), (TP53, 15, "ENSMUSP00000113213"), ("ENSP00000352610", 15, TRP53),
    (TP53, 15, "ENSMUSP00000005371"),
]


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, stdin=None, timeout=None):
            seen["args"] = args
            seen["stdin"] = stdin
            seen["timeout"] = timeout
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: seen.setdefault("guard", (verb, since)))
        return seen

    return use


def _map(queries=None, **kwargs):
    arguments = {"alignment": ALIGNMENT, "gene_trees": TREES, "gtfs": GTFS}
    arguments.update(kwargs)
    return proteins.map_residues(QUERIES if queries is None else queries, **arguments)


# ---- what is sent --------------------------------------------------------------------------------


def test_sends_the_species_lines_then_the_queries_on_stdin_and_no_option_it_was_not_given(replay):
    seen = replay("genes_map_residues_chain.json")
    _map([(TP53, 15, TRP53)])

    assert seen["args"] == ("genes", "map-residues", "--alignment", ALIGNMENT, "--gene-trees", TREES)
    assert seen["stdin"] == (
        "homo_sapiens\tgtf\tHomo_sapiens.GRCh38.116.gtf\n"
        "mus_musculus\tgtf\tMus_musculus.GRCm39.116.gtf\n"
        "rattus_norvegicus\tgtf\tRattus_norvegicus.GRCr8.116.gtf\n"
        "--\nENSP00000269305\t15\tENSMUSP00000104298\n"
    )
    assert seen["guard"] == ("genes map-residues", "0.7.0")
    assert seen["timeout"] is None


def test_proteins_and_timeout_are_sent_when_given(replay):
    seen = replay("genes_map_residues_tp53.json")
    _map(proteins=True, timeout=600, alignment=Path(ALIGNMENT))

    assert seen["args"][-1] == "--proteins"
    assert seen["args"][3] == ALIGNMENT
    assert seen["timeout"] == 600


def test_without_gene_sets_no_species_line_is_sent_and_mzlib_keeps_every_gene(replay):
    seen = replay("genes_map_residues_chain.json")
    _map([(TP53, 15, TRP53)], gtfs=None)

    assert seen["stdin"] == "--\nENSP00000269305\t15\tENSMUSP00000104298\n", "stdin starts with the separator"


def test_a_position_is_sent_as_given_never_checked_against_a_length(replay):
    # mzLib owns the rule (MapResidue throws for a position outside the source), so pyMzLib sends
    # 0, a negative and a huge position unchanged.
    seen = replay("genes_map_residues_chain.json")
    _map([(TP53, 0, TRP53), (TP53, -3, TRP53), (TP53, 10**6, TRP53)])

    assert seen["stdin"].split("--\n")[1].splitlines() == [
        f"{TP53}\t0\t{TRP53}", f"{TP53}\t-3\t{TRP53}", f"{TP53}\t1000000\t{TRP53}",
    ]


def test_gene_sets_are_sent_as_gene_set_lines(replay):
    seen = replay("genes_map_residues_chain.json")
    _map([(TP53, 15, TRP53)], gtfs={HUMAN: "human.gtf.gz"}, gene_sets={MOUSE: "mouse.genes.tsv"})

    assert seen["stdin"].split("--\n")[0].splitlines() == [
        "homo_sapiens\tgtf\thuman.gtf.gz", "mus_musculus\tgene-set\tmouse.genes.tsv",
    ]


# ---- what is refused before the bridge runs ----------------------------------------------------


@pytest.mark.parametrize(
    "queries, kwargs, match",
    [
        ([], {}, "queries is empty"),
        ("ENSP00000269305\t15\tENSMUSP00000104298", {}, "not one string"),
        ([(TP53, 15)], {}, "must be a \\(source_protein_id, position, target_protein_id\\) tuple"),
        ([(TP53, "15", TRP53)], {}, "position must be an int"),
        ([(TP53, 15.0, TRP53)], {}, "position must be an int"),
        ([(TP53, True, TRP53)], {}, "position must be an int"),
        ([("", 15, TRP53)], {}, "source must be a protein id"),
        ([(TP53, 15, None)], {}, "target must be a protein id"),
        ([(TP53, 15, "ENSMUSP\t1")], {}, "contains a tab or line break"),
        ([(TP53, 15, TRP53)], {"alignment": None}, "alignment is required"),
        ([(TP53, 15, TRP53)], {"gene_trees": ""}, "gene_trees must be a non-empty"),
        ([(TP53, 15, TRP53)], {"proteins": "yes"}, "proteins must be True or False"),
        ([(TP53, 15, TRP53)], {"gene_sets": {HUMAN: "h.tsv"}}, "'homo_sapiens' is in both gtfs and gene_sets"),
        ([(TP53, 15, TRP53)], {"gtfs": [(HUMAN, "x.gtf")]}, "gtfs must be a dict"),
    ],
)
def test_malformed_arguments_are_usage_errors_before_the_bridge_runs(monkeypatch, queries, kwargs, match):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge must not be called"))
    with pytest.raises(_bridge.UsageError, match=match):
        _map(queries, **kwargs)


# ---- what is projected ---------------------------------------------------------------------------


def test_every_outcome_is_mzlibs_name_and_a_gap_is_never_the_nearest_residue(replay):
    replay("genes_map_residues_tp53.json")
    mapped = _map(proteins=True)

    assert mapped.row_count == 9
    assert mapped.column_names == [
        "source_protein_id", "position", "target_protein_id", "outcome", "column", "target_position",
        "target_residue",
    ]
    assert mapped.columns["outcome"] == [
        "Aligned", "Aligned", "GapInTarget", "Aligned", "Aligned", "Aligned",
        "DifferentAlignments", "SourceNotInAlignment", "TargetNotInAlignment",
    ]
    assert list(mapped.outcome_counts) == list(proteins.COLUMN_OUTCOMES), "every outcome, mzLib's order"
    assert mapped.outcome_counts == {
        "Aligned": 5, "GapInTarget": 1, "DifferentAlignments": 1, "SourceNotInAlignment": 1,
        "TargetNotInAlignment": 1,
    }
    # Human S15, S20 and S392 are mouse S18, S23 and S389; S46 faces a gap in mouse.
    assert mapped.columns["target_position"] == [18, 23, None, 389, 15, 15, None, None, None]
    assert mapped.columns["target_residue"] == ["S", "S", None, "S", "S", "S", None, None, None]
    assert mapped.columns["column"] == [476, 547, 863, 3036, 476, 476, None, None, None]
    assert mapped.records[2] == {
        "source_protein_id": TP53, "position": 46, "target_protein_id": TRP53,
        "outcome": "GapInTarget", "column": 863, "target_position": None, "target_residue": None,
    }


def test_the_envelope_says_what_went_in(replay):
    replay("genes_map_residues_tp53.json")
    mapped = _map(proteins=True)

    a = mapped.alignment
    assert (a.source_file_name, a.release, a.collection) == (ALIGNMENT, "116", "protein_default")
    assert (a.alignment_count, a.protein_count) == (2, 5)
    assert a.restricted_to_gene_tree_content == mapped.gene_trees.sha256
    assert mapped.gene_trees.member_count == 5
    assert mapped.gene_trees.restricted_to_gene_sets == sorted(g.sha256 for g in mapped.gene_sets.values())
    assert isinstance(mapped.gene_sets[HUMAN], proteins.GeneSet), "the type resolve_genes reports"
    assert mapped.gene_sets[HUMAN].genome_build == "GRCh38.p14"
    assert list(mapped.gene_sets) == [HUMAN, MOUSE, RAT]
    assert any("2 queries name a protein the alignment did not keep" in c for c in mapped.caveats)


def test_proteins_hold_the_kept_canonical_proteins_with_their_ungapped_sequences(replay):
    replay("genes_map_residues_tp53.json")
    kept = _map(proteins=True).proteins

    assert kept.row_count == 5
    assert kept.columns["protein_id"] == sorted(kept.columns["protein_id"]), "mzLib's ordinal order"
    sequences = kept.sequences()
    assert len(sequences[TP53]) == 393 and len(sequences[TRP53]) == 390
    assert "-" not in "".join(sequences.values())
    assert sequences[TP53][15 - 1] == "S" and sequences[TRP53][18 - 1] == "S"
    trees = dict(zip(kept.columns["protein_id"], kept.columns["tree_id"]))
    assert trees[TP53] == trees[TRP53] == "ENSGT00950000183153"
    assert "ENSP00000352610" not in sequences, "a TP53 isoform Compara did not align"


def test_without_proteins_the_block_is_none(replay):
    replay("genes_map_residues_chain.json")
    site = _map([(TP53, 15, TRP53)])

    assert site.proteins_included is False and site.proteins is None
    assert site.records == [{
        "source_protein_id": TP53, "position": 15, "target_protein_id": TRP53, "outcome": "Aligned",
        "column": 476, "target_position": 18, "target_residue": "S",
    }]
    assert site.caveats == []
