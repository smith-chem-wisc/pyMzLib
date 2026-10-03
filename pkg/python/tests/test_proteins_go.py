"""Tests for Gene Ontology on protein groups: ``annotate_go`` and ``update_go``.

Offline. The Python layer's job is to send the right options and to turn the payload into objects
that keep a term-less row apart from a missing value. The biology - propagation, the union over
members, the status of a group with no term - is mzLib's, and the C# suite (GoTests) covers it,
including that ``out=`` is byte for byte the file mzLib's own writer produces.

Both payloads were recorded from the real bridge (not hand-written):

- ``proteins_annotate_go_pxd036557.json``: mzLib's PXD036557 MetaMorpheus 1.1.11 table, the five
  UniProt entries it names and GO release 2026-07-26 trimmed to the terms mzLib reached, with the
  organelle category map in ``fixtures/proteins/``.
- ``proteins_update_go.json``: a live download into an empty folder.

Only absolute paths were rewritten so the recording machine's layout is not committed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymzlib import _bridge, proteins

FIXTURES = Path(__file__).parent / "fixtures"
GROUPS = "PXD036557_AllQuantifiedProteinGroups.tsv"
XML = "pxd036557_proteins.xml"
OBO = "go-pxd036557.obo"
MAP = "organelle_map.tsv"


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def captured(monkeypatch):
    """Capture what a call would send, answering with a recorded payload."""
    seen: dict = {"payload": "proteins_annotate_go_pxd036557.json"}

    def fake_invoke(*args, **kwargs):
        seen["args"] = list(args)
        seen["stdin"] = kwargs.get("stdin")
        seen["timeout"] = kwargs.get("timeout")
        return _load(seen["payload"])

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    return seen


@pytest.fixture()
def go(captured) -> proteins.GoAnnotations:
    return proteins.annotate_go(GROUPS, XML, go_obo=OBO, category_map=MAP)


def _rows(go: proteins.GoAnnotations, group: str) -> list[dict]:
    return [r for r in go.records if r["protein_group"] == group]


# ---- the surface -------------------------------------------------------------------------------


def test_the_new_names_are_public():
    for name in ("annotate_go", "update_go", "GoAnnotations", "GoCategories", "GoRelease", "GoUpdate",
                 "AnnotationDatabase", "WrittenTable", "ANNOTATION_STATUSES"):
        assert name in proteins.__all__
        assert name in dir(proteins)


def test_the_statuses_are_mzlibs_four():
    assert proteins.ANNOTATION_STATUSES == ("annotated", "no_go_terms", "no_entry", "contaminant")


# ---- the payload -------------------------------------------------------------------------------


def test_every_non_decoy_group_is_annotated_and_the_decoy_skipped(go):
    assert (go.table_row_count, go.decoy_group_count, go.group_count) == (6, 1, 5)
    assert go.row_count == go.returned_count == len(go.records) == 563
    assert not go.truncated
    groups = list(dict.fromkeys(go.columns["protein_group"]))
    assert groups == ["P68363", "P05141", "P0C0S5|Q71UI9", "P02769", "P63104"]


def test_a_contaminant_is_one_term_less_row_that_says_why(go):
    (albumin,) = _rows(go, "P02769")
    assert albumin["annotation_status"] == "contaminant"
    assert albumin["go_id"] is None and albumin["aspect"] is None
    assert albumin["inherited"] is None and albumin["propagated"] is None
    assert albumin["n_with"] == 0 and albumin["accession_used"] == []


def test_the_union_keeps_both_histones_and_the_consensus_is_a_filter(go):
    histones = _rows(go, "P0C0S5|Q71UI9")
    consensus = [r for r in histones if r["n_with"] == r["n_members"]]
    assert (len(histones), len(consensus)) == (106, 57)
    only_one = next(r for r in histones if r["go_id"] == "GO:0000791")
    assert only_one["go_name"] == "euchromatin" and only_one["accession_used"] == ["P0C0S5"]


def test_evidence_is_kept_per_member(go):
    nucleosome = next(r for r in _rows(go, "P0C0S5|Q71UI9") if r["go_id"] == "GO:0000786")
    assert nucleosome["evidence_by_member"] == {"P0C0S5": ["ECO:0000353"], "Q71UI9": ["ECO:0000353"]}


def test_a_direct_annotation_is_not_propagated(go):
    inner = next(r for r in _rows(go, "P05141") if r["go_id"] == "GO:0005743")
    assert (inner["go_name"], inner["aspect"]) == ("mitochondrial inner membrane", "cellular_component")
    assert inner["propagated"] is False and inner["inherited"] is False
    assert any(r["propagated"] for r in _rows(go, "P05141")), "ancestors arrive as propagated rows"


def test_the_column_names_are_mzlibs_schema(go):
    assert go.column_names[0] == "protein_group"
    assert go.column_names[-3:] == ["go_release", "go_obo_sha256", "annotation_db_sha256"]
    assert len(go.column_names) == 19
    assert list(go.columns) == go.column_names


def test_provenance_pins_every_input(go):
    assert go.go == proteins.GoRelease("go-pxd036557.obo", go.go.sha256, "releases/2026-07-26", 412)
    assert go.annotation_database.file_type == "UniProtXml"
    assert go.annotation_database.protein_count == 5
    assert go.header["source_file_sha256"] == go.groups_file_sha256
    assert go.header["annotation_db_sha256"] == go.annotation_database.sha256
    assert set(go.columns["go_obo_sha256"]) == {go.go.sha256}


def test_the_header_counts_groups_not_rows(go):
    assert go.header["counter_q_value_max"] == "0.01"
    assert (go.header["status_annotated"], go.header["status_contaminant"]) == ("4", "1")
    assert go.header["n_multi_member_groups"] == "1"


def test_categories_join_on_go_id_and_say_the_subcategory(go):
    cats = go.categories
    assert isinstance(cats, proteins.GoCategories)
    assert (cats.map_name, cats.map_version, cats.anchor_count, cats.row_count) == ("organelle", "1", 9, 30)
    by_term = {r["go_id"]: r for r in cats.records}
    assert by_term["GO:0005743"]["subcategory"] == "mitochondrion:inner_membrane"
    assert by_term["GO:0005739"]["subcategory"] is None
    assert set(by_term) <= set(go.columns["go_id"])


def test_nothing_written_unless_asked(go):
    assert go.written is None and go.categories_written is None


def test_no_ids_were_dropped_in_a_strict_run(go):
    assert go.skip_unknown_go_ids is False and go.unresolved_go_ids == [] and go.caveats == []


# ---- arguments ---------------------------------------------------------------------------------


def test_the_required_inputs_go_on_argv(captured):
    proteins.annotate_go(GROUPS, Path(XML), go_obo=OBO, timeout=30)
    assert captured["args"] == ["proteins", "annotate-go", "--groups", GROUPS, "--database", XML,
                                "--go-obo", OBO, "--offset", "0"]
    assert captured["stdin"] is None and captured["timeout"] == 30


def test_the_optional_inputs_follow(captured):
    proteins.annotate_go(GROUPS, XML, go_obo=OBO, category_map=MAP, skip_unknown_go_ids=True,
                         out="go.tsv", categories_out="cats.tsv", limit=0, offset=5)
    args = captured["args"]
    for option, value in (("--category-map", MAP), ("--out", "go.tsv"), ("--categories-out", "cats.tsv"),
                          ("--limit", "0"), ("--offset", "5")):
        assert args[args.index(option) + 1] == value
    assert "--skip-unknown-go-ids" in args


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"limit": -1}, "limit must be an int"),
        ({"offset": -2}, "offset must be an int"),
        ({"limit": True}, "limit must be an int"),
        ({"limit": 1.5}, "limit must be an int"),
        ({"out": ""}, "out must be a non-empty"),
        ({"category_map": b"map.tsv"}, "category_map must be a non-empty"),
    ],
)
def test_bad_arguments_fail_before_a_process_starts(monkeypatch, kwargs, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match=message):
        proteins.annotate_go(GROUPS, XML, go_obo=OBO, **kwargs)


def test_go_obo_is_required(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match="go_obo is required"):
        proteins.annotate_go(GROUPS, XML, go_obo=None)  # type: ignore[arg-type]


def test_an_old_bridge_is_named_rather_than_spawned(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"verbs": ["version", "proteins read"]})
    monkeypatch.setattr(_bridge, "_VERBS_SEEN", {})
    with pytest.raises(_bridge.UsageError, match="needs the bridge from pyMzLib 0.3.0"):
        proteins.annotate_go(GROUPS, XML, go_obo=OBO)


# ---- update_go ---------------------------------------------------------------------------------


def test_update_go_reports_the_release_now_on_disk(captured):
    captured["payload"] = "proteins_update_go.json"
    update = proteins.update_go("go.obo")
    assert captured["args"] == ["proteins", "update-go", "--go-obo", "go.obo"]
    assert (update.existed_before, update.previous_sha256, update.changed) == (False, None, True)
    assert update.go.release == "releases/2026-07-26" and update.go.term_count == 48340
    assert update.url == "https://purl.obolibrary.org/obo/go.obo"


def test_update_go_needs_a_path(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match="go_obo must be a non-empty"):
        proteins.update_go("")


# ---- payload shapes the fixture does not cover -------------------------------------------------


def test_written_tables_and_an_absent_map_project():
    data = _load("proteins_annotate_go_pxd036557.json")
    data.update(written={"path": "go.tsv", "row_count": 563}, categories=None,
                categories_written=None, columns={}, returned_count=0, truncated=True)
    go = proteins.GoAnnotations._from_wire(data)
    assert go.written == proteins.WrittenTable("go.tsv", 563)
    assert go.categories is None and go.records == []
    assert proteins.WrittenTable._from_wire({"path": "c.tsv"}) == proteins.WrittenTable("c.tsv", None)
