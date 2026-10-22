"""Tests for :func:`pymzlib.proteins.build_orthology`, the projection of ``genes orthology``.

Offline. Both fixtures are recorded from the real bridge (scripts/record_fixture.py) over mzLib's
own three-species test world in ``fixtures/orthology/`` (its rows verbatim; see PROVENANCE.txt):

* ``genes_orthology_pair_status.json``: human, mouse and rat, unrestricted gene trees, with
  ``--pair-status homo_sapiens,mus_musculus``.
* ``genes_orthology_written.json``: the same inputs with ``--restrict-gene-trees`` and ``--out``.

What is pinned here is what pyMzLib sends (only the options it was given, the species lines and the
dumps on stdin, split by ``--``), what it refuses before spawning a process, and that it projects
mzLib's answer without repairing it: nulls stay ``None``, every status is mzLib's word.
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
DUMPS = [f"homologies/{s}/Compara.116.protein_default.homologies.tsv.gz" for s in GTFS]
TREES = "vertebrates.GeneTree_content.default.e116.txt.gz"


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


def _build(**kwargs):
    arguments = {"gtfs": GTFS, "gene_trees": TREES, "dumps": DUMPS}
    arguments.update(kwargs)
    return proteins.build_orthology("116", **arguments)


# ---- what is sent --------------------------------------------------------------------------------


def test_sends_the_species_lines_then_the_dumps_on_stdin_and_no_option_it_was_not_given(replay):
    seen = replay("genes_orthology_pair_status.json")
    _build()

    assert seen["args"] == ("genes", "orthology", "--release", "116", "--gene-trees", TREES)
    assert seen["stdin"] == (
        "homo_sapiens\tgtf\tHomo_sapiens.GRCh38.116.gtf\n"
        "mus_musculus\tgtf\tMus_musculus.GRCm39.116.gtf\n"
        "rattus_norvegicus\tgtf\tRattus_norvegicus.GRCr8.116.gtf\n"
        "--\n" + "\n".join(DUMPS) + "\n"
    )
    assert seen["guard"] == ("genes orthology", "0.6.0")
    assert seen["timeout"] is None


def test_every_option_given_is_sent_in_its_wire_spelling(replay):
    seen = replay("genes_orthology_written.json")
    _build(restrict_gene_trees=True, out=Path("compara-116"), pair_status=(MOUSE, HUMAN), timeout=600)

    args = seen["args"]
    assert "--restrict-gene-trees" in args
    assert args[args.index("--out") + 1] == "compara-116"
    assert args[args.index("--pair-status") + 1] == "mus_musculus,homo_sapiens", "from, then to"
    assert seen["timeout"] == 600


def test_gene_sets_are_sent_as_gene_set_lines_beside_gtf_lines(replay):
    seen = replay("genes_orthology_pair_status.json")
    _build(gtfs={HUMAN: "human.gtf.gz"}, gene_sets={MOUSE: "mouse.genes.tsv", RAT: Path("rat.genes.tsv")})

    species_lines = seen["stdin"].split("--\n")[0].splitlines()
    assert species_lines == [
        "homo_sapiens\tgtf\thuman.gtf.gz",
        "mus_musculus\tgene-set\tmouse.genes.tsv",
        "rattus_norvegicus\tgene-set\trat.genes.tsv",
    ]


def test_an_int_release_is_sent_as_its_digits(replay):
    seen = replay("genes_orthology_pair_status.json")
    proteins.build_orthology(116, gtfs=GTFS, gene_trees=TREES, dumps=DUMPS)

    assert seen["args"][2:4] == ("--release", "116")


# ---- what is refused before the bridge runs ----------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"gtfs": None}, "At least one species is required"),
        ({"gtfs": {}}, "At least one species is required"),
        ({"gtfs": [("homo_sapiens", "x.gtf")]}, "gtfs must be a dict"),
        ({"gene_sets": {HUMAN: "h.tsv"}}, "'homo_sapiens' is in both gtfs and gene_sets"),
        ({"gtfs": {"homo\tsapiens": "x.gtf"}}, "cannot hold a tab"),
        ({"gtfs": {"": "x.gtf"}}, "non-empty str"),
        ({"gtfs": {HUMAN: ""}}, "non-empty str or path"),
        ({"dumps": []}, "dumps is empty"),
        ({"dumps": ["--"]}, "is not a usable homology dump path"),
        ({"gene_trees": None}, "gene_trees is required"),
        ({"restrict_gene_trees": "yes"}, "restrict_gene_trees must be True or False"),
        ({"pair_status": "homo_sapiens,mus_musculus"}, "pair_status must be a"),
        ({"pair_status": (HUMAN,)}, "pair_status must be a"),
        ({"pair_status": (HUMAN, "mus,musculus")}, "must be a production name"),
        ({"out": ""}, "out must be a non-empty"),
    ],
)
def test_malformed_arguments_are_usage_errors_before_the_bridge_runs(monkeypatch, kwargs, match):
    def no_bridge(*args, **kwargs):
        raise AssertionError("the bridge must not be called")

    monkeypatch.setattr(_bridge, "invoke", no_bridge)
    with pytest.raises(_bridge.UsageError, match=match):
        _build(**kwargs)


@pytest.mark.parametrize("release", [None, True, "", "11 6", 116.0])
def test_a_release_that_is_not_one_is_a_usage_error(monkeypatch, release):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge must not be called"))
    with pytest.raises(_bridge.UsageError, match="release must be an Ensembl release"):
        proteins.build_orthology(release, gtfs=GTFS, gene_trees=TREES, dumps=DUMPS)


# ---- what is projected ---------------------------------------------------------------------------


def test_the_snapshot_says_what_went_in(replay):
    replay("genes_orthology_pair_status.json")
    snap = _build(pair_status=(HUMAN, MOUSE))

    assert snap.release == "116"
    assert snap.species == [HUMAN, MOUSE, RAT]
    assert (snap.format, snap.format_version) == ("ensembl-orthology-snapshot", 1)
    assert len(snap.snapshot_id) == 64
    assert snap.identical_duplicates_dropped == 0
    assert isinstance(snap.gene_sets[RAT], proteins.GeneSet), "the type resolve_genes reports"
    assert snap.gene_sets[HUMAN].gene_count == 4
    assert snap.gene_sets[HUMAN].genome_build is None, "the test GTFs carry no #!genome-build"
    assert snap.gene_trees.collection == "vertebrates"
    assert snap.gene_trees.member_count == 9
    assert snap.gene_trees.restricted_to_gene_sets is None
    assert {s: d.rows_kept for s, d in snap.dumps.items()} == {HUMAN: 1, MOUSE: 5, RAT: 1}
    assert {d.source_file_name for d in snap.dumps.values()} == {"Compara.116.protein_default.homologies.tsv.gz"}
    assert snap.written is None


def test_every_human_mouse_row_came_from_the_mouse_dump_and_none_was_lost(replay):
    replay("genes_orthology_pair_status.json")
    snap = _build(pair_status=(HUMAN, MOUSE))

    assert snap.dumps[HUMAN].rows_kept == 1, "the human dump holds only the human-rat row"
    assert snap.pair(HUMAN, MOUSE).ortholog_count == 3
    assert snap.pair(MOUSE, HUMAN) == snap.pair(HUMAN, MOUSE), "either order"
    assert snap.pair(MOUSE, MOUSE).paralog_count == 1
    assert [(p.species_a, p.species_b, p.row_count) for p in snap.pairs] == [
        (HUMAN, HUMAN, 0), (HUMAN, MOUSE, 3), (HUMAN, RAT, 1),
        (MOUSE, MOUSE, 1), (MOUSE, RAT, 2), (RAT, RAT, 0),
    ]


def test_pair_of_a_species_not_in_the_snapshot_is_a_usage_error(replay):
    replay("genes_orthology_pair_status.json")
    snap = _build(pair_status=(HUMAN, MOUSE))

    with pytest.raises(_bridge.UsageError, match="not a pair of this snapshot's species"):
        snap.pair(HUMAN, "danio_rerio")


def test_pair_status_is_mzlibs_status_for_every_gene(replay):
    replay("genes_orthology_pair_status.json")
    status = _build(pair_status=(HUMAN, MOUSE)).pair_status

    assert (status.from_species, status.to_species, status.row_count) == (HUMAN, MOUSE, 4)
    assert status.column_names == ["gene_id", "gene_name", "gene_biotype", "status"]
    assert status.columns["status"] == [
        "has_ortholog", "has_ortholog", "tree_lacks_target_species", "not_in_any_tree",
    ]
    assert status.columns["gene_name"] == [None] * 4, "no gene_name in the GTF: None, not ''"
    assert list(status.status_counts) == list(proteins.ORTHOLOGY_STATUSES), "every status, in mzLib's order"
    assert status.status_counts["not_in_gene_set"] == 0
    assert status.genes("has_ortholog") == ["ENSG00000000001", "ENSG00000000002"]
    assert status.records[3] == {
        "gene_id": "ENSG00000000004", "gene_name": None, "gene_biotype": "protein_coding",
        "status": "not_in_any_tree",
    }


def test_genes_with_an_unknown_status_is_a_usage_error(replay):
    replay("genes_orthology_pair_status.json")
    status = _build(pair_status=(HUMAN, MOUSE)).pair_status

    with pytest.raises(_bridge.UsageError, match="status must be one of"):
        status.genes("no_ortholog")


def test_the_caveats_name_the_shared_dump_name_and_the_unrestricted_trees(replay):
    replay("genes_orthology_pair_status.json")
    caveats = _build(pair_status=(HUMAN, MOUSE)).caveats

    assert any("share the file name Compara.116.protein_default.homologies.tsv.gz" in c for c in caveats)
    assert any("read unrestricted: 9 genes" in c for c in caveats)


def test_a_written_snapshot_lists_every_file_but_the_manifest(replay):
    replay("genes_orthology_written.json")
    snap = _build(restrict_gene_trees=True, out="compara-116")

    assert snap.pair_status is None
    assert snap.gene_trees.member_count == 8, "the zebrafish gene is not kept"
    assert snap.gene_trees.restricted_to_gene_sets == sorted(g.sha256 for g in snap.gene_sets.values())
    files = {f.path: f for f in snap.written.files}
    assert len(files) == 13 and "manifest.json" not in files
    assert files["views.sql"].rows is None, "not a table"
    assert files["pairs/homo_sapiens__homo_sapiens.parquet"].rows == 0, "an empty pair is still written"
    assert files["members/homo_sapiens.parquet"].rows == 3, "ENSG00000000004 is in no tree"
    assert all(len(f.sha256) == 64 and f.bytes > 0 for f in files.values())
    assert snap.written.manifest.endswith("manifest.json")
    assert snap.written.manifest.startswith(snap.written.directory)


def test_restricting_the_gene_trees_does_not_change_the_snapshot(replay):
    replay("genes_orthology_pair_status.json")
    unrestricted = _build(pair_status=(HUMAN, MOUSE))
    replay("genes_orthology_written.json")
    restricted = _build(restrict_gene_trees=True, out="compara-116")

    assert restricted.snapshot_id == unrestricted.snapshot_id
    assert restricted.pairs == unrestricted.pairs
    assert not any("unrestricted" in c for c in restricted.caveats)


def test_a_payload_without_optional_blocks_projects_to_none(replay, monkeypatch):
    data = _load("genes_orthology_written.json")
    data["written"] = None
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: data)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)

    snap = _build()
    assert snap.written is None and snap.pair_status is None
