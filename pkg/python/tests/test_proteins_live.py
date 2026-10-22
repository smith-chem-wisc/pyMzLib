"""The protein verbs through the real, packaged bridge, over the committed databases.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drifting
apart - a renamed wire key, a stdin section the bridge no longer splits, a column mzLib's
``GeneResolutionTsv`` schema gained. These run the real executable on the small databases in
``fixtures/proteins/`` and check the facts the offline tests rely on. Local files only, no network.
"""

from __future__ import annotations

import json
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


def _sequences() -> dict[str, str]:
    db = proteins.read(HUMAN_FASTA, contaminants=BOVINE_ALBUMIN, sequences=True)
    return dict(zip(db.columns["accession"], db.columns["sequence"]))


def test_align_leaves_every_option_to_mzlib_and_says_what_it_used():
    aln = proteins.align_many([("MLLLPEPTIDE", "PEPTIDE")])

    assert (aln.matrix, aln.gap_open, aln.gap_extend, aln.end_gaps, aln.max_cells) == (
        "BLOSUM62", 11, 1, "free", 200_000_000)
    assert aln.aligner_id == "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"
    assert aln.columns["score"] == [39]


def test_align_maps_canonical_dync1i2_onto_its_isoforms_as_recorded():
    seq = _sequences()
    live = proteins.align_many([(seq["Q13409"], seq["Q13409-2"]), (seq["Q13409"], seq["Q13409-3"])])
    recorded = json.loads((Path(__file__).parent / "fixtures" / "proteins_align_dync1i2.json").read_text("utf-8"))

    assert live.columns == recorded["columns"]
    assert live.residue_map.columns == recorded["residue_map"]["columns"]
    assert live.alignments[0].map_source_to_target(82) == 76


def test_align_options_reach_mzlib():
    aln = proteins.align("MLLLPEPTIDE", "PEPTIDE", gap_open=5, gap_extend=2, end_gaps="charged")

    assert aln.aligner_id == "global-affine;BLOSUM62;open=5;extend=2;end-gaps=charged"
    assert aln.score == 39 - (5 + 4 * 2)


@pytest.mark.parametrize(
    ("pairs", "kwargs", "match"),
    [
        ([("PEP-TIDE", "PEPTIDE")], {}, "gap character"),
        ([("PEPTIDE", "PEPTIDE"), ("AAAA", "AAAC")], {"max_cells": 1}, "Pair 1 was refused"),
        ([("PEPTIDE", "PEPTIDE")], {"gap_open": -1}, "--gap-open was refused by mzLib"),
        ([("PEPTIDE", "PEPTIDE")], {"gap_extend": 0}, "--gap-extend was refused by mzLib"),
    ],
)
def test_mzlibs_refusals_are_usage_errors(pairs, kwargs, match):
    with pytest.raises(pymzlib.UsageError, match=match):
        proteins.align_many(pairs, **kwargs)


def test_parse_accessions_never_repairs_what_it_is_given():
    parsed = proteins.parse_accessions(["P04637-2", " P04637", "P04637 ", "NP_000537.3_R72P", "P04637_2"])

    assert parsed.columns["input"] == ["P04637-2", " P04637", "P04637 ", "NP_000537.3_R72P", "P04637_2"]
    assert parsed.columns["namespace"] == ["uniprot", "unrecognized", "unrecognized", "refseq", "unrecognized"]
    assert parsed.columns["applied_variants"] == [None, None, None, "R72P", None]


def test_parse_accessions_and_resolve_genes_share_one_namespace_vocabulary():
    genes = proteins.resolve_genes(HUMAN_XML, gtf=GTF)
    parsed = proteins.parse_accessions(genes.columns["accession"])

    assert parsed.columns["namespace"] == genes.columns["namespace"]
    assert parsed.columns["entry_accession"] == genes.columns["entry_accession"]
    assert parsed.columns["isoform"] == genes.columns["isoform"]


# ---- genes orthology --------------------------------------------------------------------------

ORTHOLOGY = Path(__file__).parent / "fixtures" / "orthology"
ORTHOLOGY_GTFS = {
    "homo_sapiens": ORTHOLOGY / "Homo_sapiens.GRCh38.116.gtf",
    "mus_musculus": ORTHOLOGY / "Mus_musculus.GRCm39.116.gtf",
    "rattus_norvegicus": ORTHOLOGY / "Rattus_norvegicus.GRCr8.116.gtf",
}
ORTHOLOGY_DUMPS = [
    ORTHOLOGY / "homologies" / s / "Compara.116.protein_default.homologies.tsv.gz" for s in ORTHOLOGY_GTFS
]
ORTHOLOGY_TREES = ORTHOLOGY / "vertebrates.GeneTree_content.default.e116.txt.gz"


def _orthology(**kwargs):
    arguments = {"gtfs": ORTHOLOGY_GTFS, "gene_trees": ORTHOLOGY_TREES, "dumps": ORTHOLOGY_DUMPS}
    arguments.update(kwargs)
    return proteins.build_orthology("116", **arguments)


def test_orthology_matches_its_recording_over_the_committed_inputs():
    live = _orthology(pair_status=("homo_sapiens", "mus_musculus"))
    recorded = json.loads((Path(__file__).parent / "fixtures" / "genes_orthology_pair_status.json").read_text("utf-8"))

    assert live.snapshot_id == recorded["snapshot_id"], "the inputs are byte-exact (-text), so the id is too"
    assert live.pair_status.columns == recorded["pair_status"]["columns"]
    assert [p.row_count for p in live.pairs] == [p["row_count"] for p in recorded["pairs"]]
    assert live.caveats == recorded["caveats"]


def test_orthology_writes_parquet_once_and_refuses_to_overwrite(tmp_path):
    out = tmp_path / "compara-116"
    snap = _orthology(restrict_gene_trees=True, out=out)

    assert Path(snap.written.directory) == out.resolve()
    assert Path(snap.written.manifest).is_file()
    manifest = json.loads(Path(snap.written.manifest).read_text("utf-8"))
    assert manifest["snapshot_id"] == snap.snapshot_id
    assert [f["path"] for f in manifest["files"]] == [f.path for f in snap.written.files]
    assert all((out / f.path).stat().st_size == f.bytes for f in snap.written.files)
    with pytest.raises(pymzlib.UsageError, match="is not empty"):
        _orthology(out=out)


def test_orthology_without_every_dump_is_refused_by_mzlib():
    with pytest.raises(pymzlib.UsageError, match="No homology dump for mus_musculus"):
        _orthology(dumps=[ORTHOLOGY_DUMPS[0], ORTHOLOGY_DUMPS[2]])


def test_orthology_pair_status_outside_the_species_is_refused_before_reading():
    with pytest.raises(pymzlib.UsageError, match="not one of the species on stdin"):
        _orthology(pair_status=("homo_sapiens", "danio_rerio"))


def test_orthology_a_malformed_dump_is_mzlibs_correctness_failure(tmp_path):
    bad = tmp_path / "Compara.116.protein_default.homologies.tsv"
    bad.write_text("gene_stable_id\tspecies\n", encoding="utf-8")

    with pytest.raises(pymzlib.BridgeError) as caught:
        _orthology(dumps=[*ORTHOLOGY_DUMPS[:2], bad])
    assert caught.value.error_type == "InvalidDataException"
    assert "unexpected columns" in str(caught.value)
