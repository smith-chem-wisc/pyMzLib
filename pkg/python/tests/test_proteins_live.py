"""The protein verbs through the real, packaged bridge, over the committed databases.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drifting
apart - a renamed wire key, a stdin section the bridge no longer splits, a column mzLib's
``GeneResolutionTsv`` schema gained. These run the real executable on the small databases in
``fixtures/proteins/`` and check the facts the offline tests rely on. Local files only, no network.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, proteins

DATABASES = Path(__file__).parent / "fixtures" / "proteins"
HUMAN_XML = DATABASES / "human_subset.xml"
HUMAN_FASTA = DATABASES / "human_extra.fasta"
MOUSE_FASTA = DATABASES / "mouse_aifm1.fasta"
BOVINE_ALBUMIN = DATABASES / "contaminants.fasta"
GTF = DATABASES / "Homo_sapiens.GRCh38.116.gtf"
XREF = DATABASES / "Homo_sapiens.GRCh38.116.uniprot.tsv"


@pytest.fixture(autouse=True)
def built_bridge():
    # Same rule as test_thermo_live: skip in a source checkout with no built bridge, but fail under
    # CI, where the bridge is always built and a skip would hide the drift these tests exist for.
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so the protein verbs went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def test_read_answers_accession_to_organism_and_taxon_across_formats():
    db = proteins.read([HUMAN_XML, MOUSE_FASTA], accessions=["P04406", "Q9Z0X1", "P99999"])

    assert db.taxonomy() == {"P04406": "9606", "Q9Z0X1": "10090"}
    assert db.accessions_not_found == ["P99999"]


def test_read_is_identical_at_any_thread_count():
    one = proteins.read([HUMAN_XML, HUMAN_FASTA, MOUSE_FASTA], tables=proteins.TABLES, threads=1)
    many = proteins.read([HUMAN_XML, HUMAN_FASTA, MOUSE_FASTA], tables=proteins.TABLES, threads=-1)

    assert one == many


def test_go_and_ensembl_tables_come_from_the_xml():
    db = proteins.read(HUMAN_XML, tables=("go_terms", "ensembl_genes"))

    assert db.columns is None
    assert "GO:0005737" in db.go_terms.columns["go_id"]
    assert "ENSG00000111640.15" in db.ensembl_genes.columns["versioned_gene_id"]


def test_resolve_genes_resolves_gapdh_and_pins_the_release():
    genes = proteins.resolve_genes(HUMAN_XML, gtf=GTF, xref=XREF)

    gapdh = next(r for r in genes.records if r["accession"] == "P04406")
    assert (gapdh["outcome"], gapdh["gene_id"], gapdh["gene_set_release"]) == (
        "resolved", "ENSG00000111640", "116")
    assert gapdh["ensembl_xref_agrees"] is True
    assert genes.column_names[2:] == [
        "accession", "entry_accession", "isoform", "namespace", "outcome", "n_genes", "gene_id",
        "versioned_gene_id", "gene_symbol", "gene_biotype", "off_primary_genes", "uniprot_gene_name",
        "source", "search_database_sha256", "gene_set_release", "gene_set_sha256",
        "ensembl_xref_agrees", "ensembl_xref_info_type", "ensembl_xref_sha256",
    ]


def test_classify_peptides_separates_bovine_from_human_albumin():
    calls = proteins.classify_peptides(["YLYEIAR", "AEFVEVTK", "LVLNGNPLTLFQER"],
                                       [HUMAN_XML, HUMAN_FASTA], contaminants=BOVINE_ALBUMIN)

    assert calls.columns["sharing"] == ["SharedAcrossGenes", "Unique", "Unique"]
    assert calls.columns["accessions"][1] == ["P02769"]


def test_a_modified_peptide_is_refused_by_mzlib_as_a_usage_error():
    with pytest.raises(pymzlib.UsageError, match="base sequence"):
        proteins.classify_peptides(["PEPT[Phospho]IDEK"], HUMAN_XML)
