"""Tests for the protein-database projection: ``read``, ``resolve_genes``, ``classify_peptides``.

Offline. The Python layer's job is to shape arguments - one path or many, contaminants marked, a
second list after a ``--`` line on stdin - and to turn the payload into objects that keep "this
format cannot say" apart from "there is nothing". The biology is mzLib's and the C# suite covers it.

The four payloads were recorded from the real bridge through this module (not hand-written), over
the databases in ``fixtures/proteins/``. Only absolute paths were rewritten to
``fixtures/proteins/...`` so the recording machine's layout is not committed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, proteins

FIXTURES = Path(__file__).parent / "fixtures"


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def captured(monkeypatch):
    """Capture what a call would send, answering with a chosen recorded payload."""
    seen: dict = {"payload": "proteins_read_human.json"}

    def fake_invoke(*args, **kwargs):
        seen["args"] = list(args)
        seen["stdin"] = kwargs.get("stdin")
        seen["timeout"] = kwargs.get("timeout")
        return _load(seen["payload"])

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


@pytest.fixture()
def human(monkeypatch) -> proteins.ProteinDatabase:
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load("proteins_read_human.json"))
    return proteins.read(["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
                         contaminants=["contaminants.fasta"], tables=proteins.TABLES)


# ---- the package surface -----------------------------------------------------------------------


def test_the_module_is_exported():
    assert "proteins" in pymzlib.__all__
    assert pymzlib.proteins is proteins


def test_dir_lists_only_the_public_api():
    assert "read" in dir(proteins)
    assert "Iterable" not in dir(proteins)


# ---- read: the payload -------------------------------------------------------------------------


def test_read_reports_organism_and_taxon_per_accession(human):
    taxa = human.taxonomy()
    assert taxa["P04406"] == "9606"           # UniProt XML dbReference
    assert taxa["Q9Z0X1"] == "10090"          # FASTA OX=
    assert human.organisms()["P02769"] == "Bos taurus"


def test_read_keeps_every_input_in_order_with_its_role(human):
    assert [f.file_type for f in human.files] == ["UniProtXml", "Fasta", "Fasta", "Fasta"]
    assert [f.contaminant for f in human.files] == [False, False, False, True]
    assert human.columns["source_index"] == sorted(human.columns["source_index"])
    contaminated = [a for a, c in zip(human.columns["accession"], human.columns["is_contaminant"]) if c]
    assert contaminated == ["P02769"]


def test_a_fasta_says_go_and_ensembl_are_absent_rather_than_empty(human):
    fasta = human.files[1]
    assert {"go_terms", "ensembl_genes"} <= set(fasta.absent_fields)
    assert any("FASTA carries no GO" in c for c in fasta.caveats)
    assert human.files[0].absent_fields == []
    go_sources = set(human.go_terms.columns["source_index"])
    assert go_sources == {0}, "only the XML can contribute GO rows"


def test_go_terms_are_a_long_table_with_aspect_and_evidence(human):
    go = human.go_terms
    assert go.row_count == len(go.columns["go_id"]) > 0
    first = go.records[0]
    assert first["go_id"].startswith("GO:")
    assert first["aspect"] in {"BiologicalProcess", "CellularComponent", "MolecularFunction", "Unknown"}
    assert isinstance(first["evidence_codes"], list)


def test_ensembl_links_carry_stable_and_versioned_ids(human):
    ens = human.ensembl_genes.records
    gapdh = [r for r in ens if r["accession"] == "P04406"]
    assert {r["gene_id"] for r in gapdh} == {"ENSG00000111640"}
    assert {r["versioned_gene_id"] for r in gapdh} == {"ENSG00000111640.15"}
    psca = [r for r in ens if r["accession"] == "O43653"]
    assert all(r["gene_version"] is None for r in psca), "unversioned is None, never 0"


def test_mass_and_length_are_numbers(human):
    row = human.records[human.columns["accession"].index("P04406")]
    assert row["length"] == 335
    assert row["monoisotopic_mass"] == pytest.approx(36030.398, abs=0.01)


def test_filtered_read_names_the_accessions_it_did_not_find(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load("proteins_read_filtered.json"))
    db = proteins.read("human_subset.xml", accessions=["O43653", "P04406-1"], tables=("proteins", "go_terms"))
    assert db.columns["accession"] == ["O43653"]
    assert db.accessions_not_found == ["P04406-1"], "exact matching: an isoform suffix does not find the entry"
    assert db.accession_filter_count == 2
    assert db.ensembl_genes is None


def test_taxonomy_needs_the_proteins_table():
    db = proteins.ProteinDatabase(file_count=1, read_count=1, failed_count=0, record_count=0)
    with pytest.raises(pymzlib.UsageError, match="proteins table was not requested"):
        db.taxonomy()


def test_a_skipped_file_carries_its_error():
    db = proteins.ProteinDatabase._from_wire({
        "file_count": 1, "read_count": 0, "failed_count": 1, "record_count": 0,
        "files": [{"source_index": 0, "path": "x.fasta", "error": {
            "kind": "usage", "type": "usage", "message": "Protein database not found: 'x.fasta'."}}],
    })
    assert db.files[0].error.kind == "usage"
    assert "x.fasta" in db.files[0].error.message


# ---- read: what is sent ------------------------------------------------------------------------


def test_one_database_goes_on_argv(captured):
    proteins.read("human.xml")
    args = captured["args"]
    assert args[:4] == ["proteins", "read", "--path", "human.xml"]
    assert "--contaminant" not in args
    assert args[args.index("--tables") + 1] == "proteins"
    assert args[args.index("--threads") + 1] == "1"
    assert captured["stdin"] is None


def test_a_lone_contaminant_database_is_marked(captured):
    proteins.read([], contaminants="crap.fasta")
    assert captured["args"][2:5] == ["--path", "crap.fasta", "--contaminant"]


def test_many_databases_go_on_stdin_with_contaminants_tagged(captured):
    proteins.read(["a.xml", Path("b.fasta")], contaminants=["c.fasta"], threads=-1, on_error="skip")
    assert "--paths-stdin" in captured["args"]
    assert captured["stdin"] == "a.xml\nb.fasta\nc.fasta\tcontaminant\n"
    assert captured["args"][captured["args"].index("--threads") + 1] == "-1"
    assert captured["args"][captured["args"].index("--on-error") + 1] == "skip"


def test_accessions_follow_the_paths_after_a_dash_dash_line(captured):
    proteins.read(["a.xml", "b.fasta"], accessions=[" P04406 ", "Q9Z0X1"], tables=proteins.TABLES,
                  sequences=True)
    assert captured["stdin"] == "a.xml\nb.fasta\n--\nP04406\nQ9Z0X1\n"
    assert "--accessions-stdin" in captured["args"]
    assert "--sequences" in captured["args"]
    assert captured["args"][captured["args"].index("--tables") + 1] == "proteins,go_terms,ensembl_genes"


def test_with_one_database_the_accessions_are_all_of_stdin(captured):
    proteins.read("a.xml", accessions=["P04406"])
    assert captured["stdin"] == "P04406\n"


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"databases": []}, "At least one protein database"),
        ({"databases": ""}, "non-empty"),
        ({"databases": ["a.xml", "a.xml"]}, "listed twice"),
        ({"databases": "a\tb.xml"}, "tab or newline"),
        ({"databases": "--"}, "not a usable"),
        ({"databases": "a.xml", "tables": ("proteins", "peptides")}, "tables must be"),
        ({"databases": "a.xml", "tables": ()}, "tables must be"),
        ({"databases": "a.xml", "threads": 0}, "threads must be"),
        ({"databases": "a.xml", "threads": True}, "threads must be"),
        ({"databases": "a.xml", "on_error": "ignore"}, "on_error"),
        ({"databases": "a.xml", "accessions": []}, "accessions is empty"),
        ({"databases": "a.xml", "accessions": "P04406"}, "not one string"),
        ({"databases": "a.xml", "accessions": ["P04406", " "]}, "blank entry"),
        ({"databases": "a.xml", "accessions": [42]}, "must be a str"),
        ({"databases": "a.xml", "accessions": ["P0\n4406"]}, "tab or newline"),
        ({"databases": [42]}, "non-empty str or path"),
    ],
)
def test_bad_read_arguments_fail_before_a_process_starts(monkeypatch, kwargs, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("no process should start"))
    with pytest.raises(pymzlib.UsageError, match=message):
        proteins.read(**kwargs)


def test_a_single_table_name_is_accepted_as_a_string(captured):
    proteins.read("a.xml", tables="go_terms")
    assert captured["args"][captured["args"].index("--tables") + 1] == "go_terms"


# ---- resolve_genes -----------------------------------------------------------------------------


@pytest.fixture()
def genes(monkeypatch) -> proteins.GeneResolutions:
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load("genes_resolve_human.json"))
    return proteins.resolve_genes(["human_subset.xml", "human_extra.fasta"],
                                  contaminants=["contaminants.fasta"],
                                  gtf="Homo_sapiens.GRCh38.116.gtf",
                                  xref="Homo_sapiens.GRCh38.116.uniprot.tsv")


def test_resolution_pins_the_release_and_every_input_by_hash(genes):
    assert genes.gene_set.release == "116"
    assert genes.gene_set.genome_build == "GRCh38.p14"
    assert len(genes.gene_set.sha256) == 64
    assert genes.xref.release == "116"
    assert all(len(f.search_database_sha256) == 64 for f in genes.files)
    assert set(genes.columns["gene_set_sha256"]) == {genes.gene_set.sha256}


def test_every_outcome_is_counted_and_each_protein_has_one(genes):
    assert set(genes.outcome_counts) == set(proteins.OUTCOMES)
    assert sum(genes.outcome_counts.values()) == genes.protein_count
    by_accession = {r["accession"]: r for r in genes.records}
    assert by_accession["P04406"]["outcome"] == "resolved"
    assert by_accession["P04406"]["gene_symbol"] == "GAPDH"
    assert by_accession["P04406"]["ensembl_xref_agrees"] is True
    assert by_accession["O43653"]["outcome"] == "off_primary_only"
    assert by_accession["P02768"]["outcome"] == "not_in_source"
    assert by_accession["P02769"]["outcome"] == "contaminant_not_mapped"


def test_resolution_warns_that_fasta_proteins_cannot_resolve(genes):
    assert any("came from a FASTA" in c for c in genes.caveats)


def test_resolve_sends_the_gene_set_and_xref(captured):
    captured["payload"] = "genes_resolve_human.json"
    proteins.resolve_genes("a.xml", gtf=Path("g.gtf.gz"), xref="x.tsv.gz")
    args = captured["args"]
    assert args[:2] == ["genes", "resolve"]
    assert args[args.index("--gtf") + 1] == "g.gtf.gz"
    assert args[args.index("--xref") + 1] == "x.tsv.gz"
    assert "--gene-set" not in args


def test_resolve_accepts_a_compact_gene_set(captured):
    captured["payload"] = "genes_resolve_human.json"
    proteins.resolve_genes(["a.xml", "b.xml"], gene_set="genes.tsv")
    assert captured["args"][captured["args"].index("--gene-set") + 1] == "genes.tsv"
    assert captured["stdin"] == "a.xml\nb.xml\n"


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({}, "exactly one of gtf"),
        ({"gtf": "g.gtf", "gene_set": "s.tsv"}, "exactly one of gtf"),
        ({"gtf": " "}, "exactly one of gtf"),
        ({"gtf": "g.gtf", "xref": ""}, "xref must be"),
    ],
)
def test_bad_resolve_arguments_fail_before_a_process_starts(monkeypatch, kwargs, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("no process should start"))
    with pytest.raises(pymzlib.UsageError, match=message):
        proteins.resolve_genes("a.xml", **kwargs)


# ---- classify_peptides -------------------------------------------------------------------------


@pytest.fixture()
def calls(monkeypatch) -> proteins.PeptideClassification:
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load("proteins_classify_peptides.json"))
    return proteins.classify_peptides(
        ["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
        ["human_subset.xml", "human_extra.fasta"], contaminants=["contaminants.fasta"])


def test_each_sharing_class_has_a_real_example(calls):
    sharing = calls.sharing_of()
    assert sharing["VGVNGFGR"] == proteins.UNIQUE
    assert sharing["ALSEQINIFFDYSGR"] == proteins.SHARED_WITHIN_GENE     # DYNC1I2 isoforms
    assert sharing["YLYEIAR"] == proteins.SHARED_ACROSS_GENES            # human and bovine albumin
    assert sharing["AEFVEVTK"] == proteins.UNIQUE                        # bovine only: a contaminant
    assert sharing["PEPTIDEK"] == proteins.NOT_IN_DATABASE
    assert sum(calls.sharing_counts.values()) == calls.peptide_count == 6


def test_i_and_l_are_one_residue_and_the_peptide_is_echoed_as_given(calls):
    row = calls.records[1]
    assert row["peptide"] == "LVLNGNPLTLFQER"          # GAPDH has LVINGNPITIFQER
    assert row["accessions"] == ["P04406"]
    assert calls.i_and_l_equivalent is True


def test_shared_across_genes_lists_both_proteins_and_no_common_gene(calls):
    row = calls.records[3]
    assert row["accessions"] == ["P02768", "P02769"]
    assert row["shared_gene_keys"] == []


def test_peptides_follow_the_paths_on_stdin(captured):
    captured["payload"] = "proteins_classify_peptides.json"
    proteins.classify_peptides(("PEPTIDEK", " VGVNGFGR"), ["a.xml"], contaminants="c.fasta", threads=2)
    assert captured["args"][:3] == ["proteins", "classify-peptides", "--paths-stdin"]
    assert captured["stdin"] == "a.xml\nc.fasta\tcontaminant\n--\nPEPTIDEK\nVGVNGFGR\n"
    assert "--on-error" not in captured["args"], "there is no skip mode for one search space"


def test_with_one_database_the_peptides_are_all_of_stdin(captured):
    captured["payload"] = "proteins_classify_peptides.json"
    proteins.classify_peptides(["PEPTIDEK"], "a.xml")
    assert captured["args"][2:4] == ["--path", "a.xml"]
    assert captured["stdin"] == "PEPTIDEK\n"


@pytest.mark.parametrize(
    "peptides, message",
    [([], "At least one peptide"), ("PEPTIDEK", "not one string"), (["PEPTIDEK", ""], "blank entry")],
)
def test_bad_peptides_fail_before_a_process_starts(monkeypatch, peptides, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("no process should start"))
    with pytest.raises(pymzlib.UsageError, match=message):
        proteins.classify_peptides(peptides, "a.xml")
