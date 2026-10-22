"""Protein databases: what each protein *is*, which gene it belongs to, and whether a peptide
identifies it.

Questions a search result cannot answer on its own, each answered from the protein database you
searched, by mzLib:

**What is this accession?** :func:`read` loads a UniProt XML or a FASTA and returns one row per
protein - organism, NCBI taxonomy id, gene names, length, monoisotopic mass - and, on request, its
Gene Ontology terms and Ensembl gene links as long tables::

    >>> db = pymzlib.proteins.read(
    ...     ["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
    ...     contaminants=["contaminants.fasta"], tables=pymzlib.proteins.TABLES)
    >>> db.taxonomy()["Q9Z0X1"], db.organisms()["P02769"]
    ('10090', 'Bos taurus')
    >>> db.go_terms.columns["term_name"][:3]
    ['cytoplasm', 'cytosol', 'extracellular exosome']

**Which gene is it, reproducibly?** :func:`resolve_genes` resolves every protein to a stable
Ensembl gene id, counted against **a gene set you supply and pin** - an Ensembl GTF for one release -
and says, per protein, *how* it resolved (``resolved``, ``multi_gene``, ``off_primary_only``,
``not_in_source``, ...). Every row carries the sha256 of the database, the GTF and the optional
cross-reference table it was computed from::

    >>> genes = pymzlib.proteins.resolve_genes(
    ...     ["human_subset.xml", "human_extra.fasta"], contaminants=["contaminants.fasta"],
    ...     gtf="Homo_sapiens.GRCh38.116.gtf", xref="Homo_sapiens.GRCh38.116.uniprot.tsv")
    >>> genes.gene_set.release, genes.gene_set.genome_build, genes.outcome_counts["resolved"]
    ('116', 'GRCh38.p14', 1)

**Does this peptide identify one protein?** :func:`classify_peptides` classifies each peptide as
``Unique``, ``SharedWithinGene``, ``SharedAcrossGenes`` or ``NotInDatabase``, treating **I and L as
the same residue**, because a mass spectrometer cannot tell them apart::

    >>> calls = pymzlib.proteins.classify_peptides(
    ...     ["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
    ...     ["human_subset.xml", "human_extra.fasta"], contaminants=["contaminants.fasta"])
    >>> calls.columns["sharing"]
    ['Unique', 'Unique', 'SharedWithinGene', 'SharedAcrossGenes', 'Unique', 'NotInDatabase']

**What does each protein group do, every member kept?** :func:`annotate_go` reads the protein-group
table MetaMorpheus wrote, the UniProt XML it searched and a go.obo you pinned, and returns one row
per (group, GO term) that **any** member holds - directly or through an ancestor in the ontology -
naming the members that carry it. Consensus and direct-only views are filters on the rows.
:func:`update_go` is the one function that fetches a go.obo; :func:`annotate_go` never does::

    >>> go = pymzlib.proteins.annotate_go(
    ...     "PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
    ...     go_obo="go-pxd036557.obo", category_map="organelle_map.tsv")
    >>> go.group_count, go.row_count, go.header["status_contaminant"]
    (5, 563, '1')

**Where is this residue in the other sequence?** :func:`align` aligns two sequences with mzLib's
global aligner (BLOSUM62, affine gaps) and maps every residue to its counterpart in the other - an
isoform, an ortholog, an older release - or to ``None`` when it faces a gap. :func:`align_many`
aligns a list of pairs in one call::

    >>> aln = pymzlib.proteins.align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")
    >>> aln.map_source_to_target(12), aln.map_source_to_target(10)
    (9, None)

**What entry, isoform and variants does this accession name?** :func:`parse_accessions` reads
``P04637-2``, ``NP_000537.3_R72P`` or ``P04637_S70N_A80T`` back into entry, isoform or version, and
applied variants, with mzLib's grammars, and leaves anything else ``"unrecognized"`` rather than
guess::

    >>> parsed = pymzlib.proteins.parse_accessions(["P04637", "P04637-2", "Q13409-2",
    ...     "NP_000537.3_R72P", "P04637_S70N_A80T", "P04637_2", "DECOY_P02768",
    ...     "sp|P04637|P53_HUMAN", "P04637|P02768"])
    >>> parsed.columns["namespace"][4:7]
    ['uniprot', 'unrecognized', 'unrecognized']

**Which genes are orthologous across species?** :func:`build_orthology` builds mzLib's snapshot of
Ensembl Compara's homologies among your species, from their gene sets, the gene trees and **every**
species' homology dump (Compara splits each pair's relationships between the pair's two dumps), and
writes it as Parquet with ``out=``. With ``pair_status=`` it says, gene by gene, how far Compara's
inference got in one other species::

    >>> gtfs = {"homo_sapiens": "Homo_sapiens.GRCh38.116.gtf",
    ...         "mus_musculus": "Mus_musculus.GRCm39.116.gtf",
    ...         "rattus_norvegicus": "Rattus_norvegicus.GRCr8.116.gtf"}
    >>> snap = pymzlib.proteins.build_orthology(
    ...     "116", gtfs=gtfs, gene_trees="vertebrates.GeneTree_content.default.e116.txt.gz",
    ...     dumps=[f"homologies/{s}/Compara.116.protein_default.homologies.tsv.gz" for s in gtfs],
    ...     pair_status=("homo_sapiens", "mus_musculus"))
    >>> snap.pair_status.status_counts["has_ortholog"], snap.pair_status.row_count
    (2, 4)

These examples are the calls that recorded the payloads in pyMzLib's test fixtures, over the small
databases in ``pkg/python/tests/fixtures/proteins/`` (and, for orthology, mzLib's own three-species
test world in ``pkg/python/tests/fixtures/orthology/``). The GO example's table is a real MetaMorpheus
1.1.11 search (PXD036557).

:func:`read`, :func:`resolve_genes` and :func:`classify_peptides` take **one database or many**.
Pass a list and they are read in one bridge call, in order; ``threads`` says how many are read at
once (default 1), and the answer is identical at any value.
Mark contaminant databases with ``contaminants=`` - it changes answers: a contaminant is never mapped
to a gene, and a peptide it shares with a target is shared.

**What these functions will not tell you, and say so instead:**

- *A FASTA knows no GO terms and no Ensembl genes.* Its headers carry organism (``OS=``), taxonomy
  (``OX=``) and gene name (``GN=``) and nothing else, so it contributes no ``go_terms`` or
  ``ensembl_genes`` rows and :func:`resolve_genes` calls every FASTA protein ``not_in_source``. That
  is the format's silence, not a biological absence, and each result says so in
  :attr:`DatabaseFile.absent_fields` and ``caveats``. Read the UniProt XML of the same proteome.
- *No decoys are generated and no variants are expanded.* A database is read as written. A decoy
  already in the file (an accession starting ``DECOY``) is kept and flagged ``is_decoy``. A UniProt
  XML that records a **genotype** still has those variants applied by mzLib, which renames the
  accession (``P38936_C117Y``); the file's caveats say how many.
- *Masses are of the unmodified sequence as written* - the precursor, initiator methionine and
  signal peptide included - so they are not the mass of the mature protein.

Wire verbs: ``proteins read``, ``genes resolve``, ``genes orthology``,
``proteins classify-peptides``, ``proteins annotate-go``, ``proteins update-go``,
``proteins align`` and ``proteins parse-accession``; their language-neutral specs live in the bridge repository under
``design/verbs/``.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Union

from . import _bridge
from .readers import _Table

__all__ = [
    "ANNOTATION_STATUSES",
    "END_GAPS",
    "NAMESPACES",
    "NOT_IN_DATABASE",
    "ORTHOLOGY_STATUSES",
    "OUTCOMES",
    "SHARED_ACROSS_GENES",
    "SHARED_WITHIN_GENE",
    "TABLES",
    "UNIQUE",
    "Alignment",
    "Alignments",
    "AnnotationDatabase",
    "DatabaseFile",
    "FileError",
    "GeneResolutions",
    "GeneSet",
    "GeneTrees",
    "GoAnnotations",
    "GoCategories",
    "GoRelease",
    "GoUpdate",
    "HomologyDump",
    "OrthologySnapshot",
    "PairStatus",
    "ParsedAccessions",
    "PeptideClassification",
    "ProteinDatabase",
    "SnapshotFile",
    "SpeciesPair",
    "Table",
    "WrittenSnapshot",
    "WrittenTable",
    "XrefTable",
    "align",
    "align_many",
    "annotate_go",
    "build_orthology",
    "classify_peptides",
    "parse_accessions",
    "read",
    "resolve_genes",
    "update_go",
]

#: The tables :func:`read` can return. ``"proteins"`` alone is the default.
TABLES = ("proteins", "go_terms", "ensembl_genes")

#: Every protein containing the peptide has the same sequence. Identical sequences under two
#: accessions count as one, because no peptide can tell them apart.
UNIQUE = "Unique"
#: The peptide is in more than one distinct sequence, and all of them share a gene: it supports the
#: gene, not any one isoform.
SHARED_WITHIN_GENE = "SharedWithinGene"
#: The peptide is in sequences with no gene in common.
SHARED_ACROSS_GENES = "SharedAcrossGenes"
#: No target protein in the databases contains the peptide.
NOT_IN_DATABASE = "NotInDatabase"

#: Every outcome :func:`resolve_genes` can report, as written on the wire (mzLib's
#: ``GeneResolutionTsv.OutcomeName``). Every protein gets exactly one.
OUTCOMES = (
    "resolved",
    "multi_gene",
    "off_primary_only",
    "not_in_source",
    "unrecognized_accession",
    "contaminant_not_mapped",
)

#: Every ``annotation_status`` :func:`annotate_go` can report, as mzLib writes it
#: (``GoAnnotationTsv.StatusName``). A row with a GO term is always ``annotated``; the other three
#: mark the single term-less row of a group that has no term, and say why.
ANNOTATION_STATUSES = ("annotated", "no_go_terms", "no_entry", "contaminant")

#: The pyMzLib release whose bridge first dispatches ``proteins annotate-go`` and ``update-go``.
_GO_SINCE = "0.3.0"

#: The end-gap rules :func:`align` takes, as mzLib's aligner id spells them: ``"free"`` (mzLib's
#: default: an overhang at either end costs nothing) or ``"charged"`` (every gap costs).
END_GAPS = ("free", "charged")

#: Every ``namespace`` :func:`parse_accessions` can report: mzLib's ``AccessionNamespace``, lower
#: case, exactly as :func:`resolve_genes` writes it.
NAMESPACES = ("uniprot", "refseq", "unrecognized")

#: The pyMzLib release whose bridge first dispatches ``proteins align`` and ``parse-accession``.
_ALIGN_SINCE = "0.6.0"

#: Every status :func:`build_orthology` can report in :attr:`OrthologySnapshot.pair_status`, as
#: mzLib writes it (``OrthologySnapshot.StatusName``), in mzLib's order. Every gene gets exactly one,
#: so "Compara found nothing", "there was nothing to compare against" and "no call was attempted"
#: never share an empty result. ``not_in_gene_set`` is for a gene id from elsewhere, and never occurs
#: in a pair-status table, whose rows are the gene set itself.
ORTHOLOGY_STATUSES = (
    "has_ortholog",
    "no_edge_in_shared_tree",
    "tree_lacks_target_species",
    "not_in_any_tree",
    "not_in_gene_set",
)

#: The pyMzLib release whose bridge first dispatches ``genes orthology``.
_ORTHOLOGY_SINCE = "0.6.0"

_PathLike = Union[str, "os.PathLike[str]"]

#: The line that separates the database paths from a second list when both travel on stdin.
_SECTION = "--"


def __dir__() -> list[str]:
    """The public API only, so ``dir()`` is a map of what to call rather than of our imports."""
    return sorted(__all__)


# ---- result types ------------------------------------------------------------------------------


@dataclass(frozen=True)
class FileError:
    """Why one database could not be read, under ``on_error="skip"``.

    Attributes:
        kind: ``"usage"`` (the file is missing, or not a database this reads) or ``"correctness"``
            (mzLib could not parse it).
        type: The error type as the bridge classified it - ``"usage"``, or the .NET exception name.
        message: What went wrong, naming the file.
    """

    kind: str
    type: str
    message: str

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> FileError:
        return cls(kind=data.get("kind", ""), type=data.get("type", ""), message=data.get("message", ""))


@dataclass(frozen=True)
class DatabaseFile:
    """What happened to one input database.

    One per input, in input order, whether it was read or not - so a result for many databases
    always says which of them contributed.

    Attributes:
        source_index: Position in the input: databases first, then ``contaminants``. The
            ``source_index`` column of every table points back here.
        path: The absolute path.
        file_type: ``"UniProtXml"`` or ``"Fasta"``, chosen from the extension (``.xml``, or
            ``.fasta``/``.fa``/``.faa``/``.fas``, each optionally ``.gz``). ``None`` if the file
            failed before its type was known.
        reader: The mzLib loader used, e.g. ``"ProteinDbLoader.LoadProteinXML"``.
        contaminant: Whether this database was loaded as contaminants.
        protein_count: Proteins mzLib loaded from it, before any accession filter.
        decoy_count: Of those, how many the file itself marks as decoys (accession starts
            ``DECOY``). None are generated.
        record_count: Rows this database contributed to the main table. For :func:`read`, proteins
            that passed the accession filter; for :func:`resolve_genes`, resolution rows; for
            :func:`classify_peptides`, target proteins searched.
        search_database_sha256: For :func:`resolve_genes` only: the lower-case hex sha256 of the
            database's **decompressed** bytes, the key mzLib's resolver pins every row to. ``None``
            from the other two functions.
        caveats: What this file cannot tell you, or what mzLib did to it: a FASTA has no GO or
            Ensembl source; genotype variants were applied; FASTA lines the loader skipped.
        absent_fields: Tables and columns this **format** has no source for, so that an empty value
            reads as "this file cannot say" and never as "there is none". For a FASTA:
            ``["go_terms", "ensembl_genes", "ensembl_gene_ids"]``.
        error: ``None``, or why the file could not be read (``on_error="skip"`` only).
    """

    source_index: int
    path: str
    file_type: str | None
    reader: str | None
    contaminant: bool
    protein_count: int
    decoy_count: int
    record_count: int
    search_database_sha256: str | None = None
    caveats: list[str] = field(default_factory=list)
    absent_fields: list[str] = field(default_factory=list)
    error: FileError | None = None

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> DatabaseFile:
        error = data.get("error")
        return cls(
            source_index=int(data.get("source_index", 0)),
            path=data.get("path", ""),
            file_type=data.get("file_type"),
            reader=data.get("reader"),
            contaminant=bool(data.get("contaminant", False)),
            protein_count=int(data.get("protein_count", 0)),
            decoy_count=int(data.get("decoy_count", 0)),
            record_count=int(data.get("record_count", 0)),
            search_database_sha256=data.get("search_database_sha256"),
            caveats=list(data.get("caveats") or []),
            absent_fields=list(data.get("absent_fields") or []),
            error=FileError._from_wire(error) if error else None,
        )


@dataclass(frozen=True)
class Table(_Table):
    """A long table: one row per record, column-major.

    Attributes:
        row_count: Rows in the table.
        column_names: Column order.
        columns: Column name -> list of values, all the same length. Hand it straight to
            ``pandas.DataFrame(table.columns)``, or use :attr:`records` for a list of dicts.
    """

    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any] | None) -> Table | None:
        if data is None:
            return None
        return cls(
            row_count=int(data.get("row_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
        )


def _files(data: Mapping[str, Any]) -> list[DatabaseFile]:
    return [DatabaseFile._from_wire(f) for f in data.get("files") or []]


@dataclass(frozen=True)
class ProteinDatabase(_Table):
    """What :func:`read` returns: the proteins, and the long tables asked for.

    The ``proteins`` table (:attr:`columns`) has one row per protein, in database order:

    ======================  =============  ======================================================
    column                  type / unit    meaning; ``None`` means
    ======================  =============  ======================================================
    ``source_index``        int            which input (:attr:`files`) the row came from
    ``source_path``         str            that input's absolute path
    ``accession``           str            as the database wrote it (``P04406``, ``Q13409-2``)
    ``name``                str            entry name (``G3P_HUMAN``); ``None``: not recorded
    ``full_name``           str            protein name; ``None``: not recorded
    ``organism``            str            ``OS=`` / ``<organism>`` scientific name;
                                           ``None``: not recorded
    ``ncbi_taxonomy_id``    str            NCBI taxon (``"9606"``) from XML or FASTA ``OX=``;
                                           ``None``: the entry records none. Kept a string: it
                                           is an identifier, not a quantity.
    ``primary_gene_name``   str            the primary gene name (``GN=``); ``None``: none given
    ``gene_names``          list[str]      every gene name the entry gives, primary first,
                                           synonyms, ORF and locus names after; may be empty
    ``length``              int, residues  sequence length
    ``monoisotopic_mass``   float, Da      unmodified sequence plus one water; ``None``: the
                                           sequence holds a letter with no defined mass (X, B,
                                           Z, J)
    ``is_contaminant``      bool           loaded from a ``contaminants=`` database
    ``is_decoy``            bool           the file marks it a decoy (``DECOY`` prefix)
    ``is_entrapment``       bool           mzLib's entrapment flag (accession contains
                                           ``Random``)
    ``ensembl_gene_ids``    list[str]      distinct stable Ensembl gene ids, ordinal order;
                                           empty for FASTA (see ``absent_fields``)
    ``sequence``            str            only with ``sequences=True``
    ======================  =============  ======================================================

    Attributes:
        file_count: Database files given.
        read_count: Database files read.
        failed_count: Database files that failed (``on_error="skip"`` only; otherwise the call
            raises).
        record_count: Proteins that passed the accession filter - the rows of the ``proteins``
            table when it was requested.
        tables: The tables returned, in :data:`TABLES` order.
        accession_filter_count: Distinct accessions in the filter, or ``None`` with no filter.
        accessions_not_found: Filter accessions no database contained, in the order given; ``None``
            with no filter. **Check this**: matching is exact, so ``P04406-1`` does not find
            ``P04406``.
        caveats: Traps in this result as a whole, e.g. proteins with no computable mass.
        files: One :class:`DatabaseFile` per input, in input order.
        column_names: The ``proteins`` table's column order; ``None`` if it was not requested.
        columns: The ``proteins`` table; ``None`` if it was not requested.
        go_terms: One row per (protein, GO term), or ``None`` if not requested. Columns:
            ``source_index``, ``source_path``, ``accession``, ``go_id`` (``"GO:0005737"``),
            ``aspect`` (``"BiologicalProcess"``, ``"CellularComponent"``,
            ``"MolecularFunction"``, or ``"Unknown"`` when UniProt gave no ``C:``/``F:``/``P:``
            prefix - never guessed), ``term_name`` (prefix removed; ``None`` when the reference
            carried no name), ``evidence_codes`` (ECO ids, sorted; may be empty) and ``projects``
            (annotating sources such as ``UniProtKB``, sorted). A GO id repeated with several lines
            of evidence is **one** row with the evidence unioned.
        ensembl_genes: One row per (protein, Ensembl transcript), or ``None`` if not requested.
            Columns: ``source_index``, ``source_path``, ``accession``, ``transcript_id`` (as UniProt
            wrote it, versioned), ``protein_id`` (``None`` when UniProt gave none), ``gene_id``
            (stable - **join on this**), ``versioned_gene_id`` (as written, e.g.
            ``ENSG00000111640.15``) and ``gene_version`` (int; ``None`` when unversioned, which is
            not version 0).
    """

    file_count: int
    read_count: int
    failed_count: int
    record_count: int
    tables: list[str] = field(default_factory=list)
    accession_filter_count: int | None = None
    accessions_not_found: list[str] | None = None
    caveats: list[str] = field(default_factory=list)
    files: list[DatabaseFile] = field(default_factory=list)
    column_names: list[str] | None = None
    columns: dict[str, list[Any]] | None = None
    go_terms: Table | None = None
    ensembl_genes: Table | None = None

    def taxonomy(self) -> dict[str, str | None]:
        """Accession -> NCBI taxonomy id, for every protein in the table.

        Raises:
            UsageError: the ``proteins`` table was not requested.

        Examples:
            >>> db = read(["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
            ...           contaminants=["contaminants.fasta"], tables=TABLES)
            >>> db.taxonomy()["Q9Z0X1"]
            '10090'
        """
        columns = self._require_proteins()
        return dict(zip(columns["accession"], columns["ncbi_taxonomy_id"]))

    def organisms(self) -> dict[str, str | None]:
        """Accession -> organism name, for every protein in the table.

        Raises:
            UsageError: the ``proteins`` table was not requested.
        """
        columns = self._require_proteins()
        return dict(zip(columns["accession"], columns["organism"]))

    def _require_proteins(self) -> dict[str, list[Any]]:
        if self.columns is None:
            raise _bridge.UsageError(
                "The proteins table was not requested. Call read(..., tables=(\"proteins\", ...))."
            )
        return self.columns

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> ProteinDatabase:
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            record_count=int(data.get("record_count", 0)),
            tables=list(data.get("tables") or []),
            accession_filter_count=data.get("accession_filter_count"),
            accessions_not_found=(
                None if data.get("accessions_not_found") is None
                else list(data["accessions_not_found"])
            ),
            caveats=list(data.get("caveats") or []),
            files=_files(data),
            column_names=None if data.get("column_names") is None else list(data["column_names"]),
            columns=data.get("columns"),
            go_terms=Table._from_wire(data.get("go_terms")),
            ensembl_genes=Table._from_wire(data.get("ensembl_genes")),
        )


@dataclass(frozen=True)
class GeneSet:
    """The Ensembl gene set every resolution was counted against, and where it came from.

    Attributes:
        source_file_name: The GTF's file name. With ``gene_set=``, still the **GTF's** name: a
            compact table carries the provenance of the GTF it was made from.
        sha256: Lower-case hex sha256 of the GTF's bytes as read - the compressed bytes for a
            ``.gz``, which is what Ensembl publishes and checksums.
        release: The Ensembl release parsed from the file name (``"116"`` from
            ``Homo_sapiens.GRCh38.116.gtf.gz``); ``None`` when the name carries none. Keep
            Ensembl's file name.
        genome_build: The ``#!genome-build`` header, e.g. ``"GRCh38.p14"``; ``None`` if absent.
        genebuild_last_updated: The ``#!genebuild-last-updated`` header; ``None`` if absent.
        gene_count: Genes in the set.
    """

    source_file_name: str
    sha256: str
    release: str | None
    genome_build: str | None
    genebuild_last_updated: str | None
    gene_count: int

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> GeneSet:
        return cls(
            source_file_name=data.get("source_file_name", ""),
            sha256=data.get("sha256", ""),
            release=data.get("release"),
            genome_build=data.get("genome_build"),
            genebuild_last_updated=data.get("genebuild_last_updated"),
            gene_count=int(data.get("gene_count", 0)),
        )


@dataclass(frozen=True)
class XrefTable:
    """Ensembl's own accession-to-gene cross-references, when ``xref=`` was given.

    Attributes:
        source_file_name: The file's name (e.g. ``Homo_sapiens.GRCh38.116.uniprot.tsv.gz``).
        sha256: Lower-case hex sha256 of the file's bytes as read.
        release: The Ensembl release from the file name; ``None`` if absent.
        accession_count: Distinct accessions the table links.
    """

    source_file_name: str
    sha256: str
    release: str | None
    accession_count: int

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any] | None) -> XrefTable | None:
        if data is None:
            return None
        return cls(
            source_file_name=data.get("source_file_name", ""),
            sha256=data.get("sha256", ""),
            release=data.get("release"),
            accession_count=int(data.get("accession_count", 0)),
        )


@dataclass(frozen=True)
class GeneResolutions(_Table):
    """What :func:`resolve_genes` returns: one row per (protein, gene), or one outcome row per
    protein with no gene. Never a ``|``-joined cell, and never a pick among several genes.

    The column names after ``source_index`` and ``source_path`` are **mzLib's own**
    ``GeneResolutionTsv`` schema, so this table and one written by mzLib's TSV writer mean the same:

    ==========================  =========  ===========================================================
    column                      type       meaning; ``None`` means
    ==========================  =========  ===========================================================
    ``accession``               str        as the search reported it, variant suffix and all
    ``entry_accession``         str        the UniProt entry / unversioned RefSeq; verbatim otherwise
    ``isoform``                 int        UniProt isoform number; ``None``: no ``-N`` suffix
    ``namespace``               str        ``uniprot``, ``refseq`` or ``unrecognized``
    ``outcome``                 str        one of :data:`OUTCOMES`
    ``n_genes``                 int        genes on the gene set's assembly (0, 1, or more)
    ``gene_id``                 str        stable id (``ENSG...``); ``None``: no gene on the assembly
    ``versioned_gene_id``       str        the id as the database wrote it; ``None``: no gene
    ``gene_symbol``             str        the **release's** symbol (display only); ``None``: no gene,
                                           or the GTF gives it no name
    ``gene_biotype``            str        the release's biotype; ``None``: no gene
    ``off_primary_genes``       int        linked genes **outside** the set (ALT haplotypes, patches):
                                           dropped from ``n_genes``, but counted here
    ``uniprot_gene_name``       str        the database's own primary gene name, a label;
                                           ``None``: none
    ``source``                  str        ``search_database_dbreference``, or ``ensembl_xref`` for a
                                           gene only Ensembl links (adds an answer, never changes the
                                           outcome)
    ``search_database_sha256``  str        sha256 of the decompressed database
    ``gene_set_release``        str        Ensembl release; ``None``: the GTF name carries none
    ``gene_set_sha256``         str        sha256 of the GTF
    ``ensembl_xref_agrees``     bool       whether Ensembl's xref links this pair; ``None``: no gene on
                                           the row, or no ``xref=`` given (**unknown, not false**)
    ``ensembl_xref_info_type``  str        Ensembl's evidence (``DIRECT``, ``SEQUENCE_MATCH``,
                                           ``INFERRED_PAIR``); ``None``: no xref link
    ``ensembl_xref_sha256``     str        sha256 of the xref file; ``None``: no ``xref=``
    ==========================  =========  ===========================================================

    Attributes:
        file_count: Database files given.
        read_count: Database files read.
        failed_count: Database files that failed (``on_error="skip"`` only).
        protein_count: Proteins resolved.
        record_count: Rows in the table: at least one per protein, one per gene for
            ``multi_gene``, plus any ``ensembl_xref`` rows.
        gene_set: The gene set and its provenance.
        xref: The xref table and its provenance, or ``None``.
        outcome_counts: Outcome -> number of **proteins** with it. Every key in :data:`OUTCOMES` is
            present, zeros included.
        caveats: Traps in this result: FASTA proteins, a GTF with no release in its name, no xref.
        files: One :class:`DatabaseFile` per input, each with its ``search_database_sha256``.
        column_names: Column order.
        columns: Column name -> list of values.
    """

    file_count: int
    read_count: int
    failed_count: int
    protein_count: int
    record_count: int
    gene_set: GeneSet
    xref: XrefTable | None = None
    outcome_counts: dict[str, int] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)
    files: list[DatabaseFile] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> GeneResolutions:
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            protein_count=int(data.get("protein_count", 0)),
            record_count=int(data.get("record_count", 0)),
            gene_set=GeneSet._from_wire(data.get("gene_set") or {}),
            xref=XrefTable._from_wire(data.get("xref")),
            outcome_counts=dict(data.get("outcome_counts") or {}),
            caveats=list(data.get("caveats") or []),
            files=_files(data),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
        )


@dataclass(frozen=True)
class PeptideClassification(_Table):
    """What :func:`classify_peptides` returns: one row per peptide, in the order given.

    ====================  =========  ===============================================================
    column                type       meaning
    ====================  =========  ===============================================================
    ``peptide``           str        exactly as given - **not** I/L-folded
    ``sharing``           str        :data:`UNIQUE`, :data:`SHARED_WITHIN_GENE`,
                                     :data:`SHARED_ACROSS_GENES` or :data:`NOT_IN_DATABASE`
    ``accession_count``   int        target proteins containing it
    ``accessions``        list[str]  those proteins, distinct, ordinal order; empty when not found
    ``shared_gene_keys``  list[str]  the gene keys common to all of them, e.g.
                                     ``"ensembl:ENSG00000111640"``, ``"gene:Homo sapiens:GAPDH"``,
                                     ``"entry:P04406"``. Empty for ``NotInDatabase`` and
                                     ``SharedAcrossGenes``
    ====================  =========  ===============================================================

    Attributes:
        file_count: Database files searched.
        peptide_count: Rows: peptides, one per input peptide, duplicates included.
        target_protein_count: Proteins searched. Contaminants count: they are real sequences in
            the search space.
        decoy_proteins_ignored: Decoy proteins in the databases, which are never searched.
        i_and_l_equivalent: Always ``True``: I and L are one residue for matching. On the wire so
            that the rule travels with every result rather than living only in documentation.
        sharing_counts: Class -> peptides in it; every class present, zeros included.
        files: One :class:`DatabaseFile` per input.
        column_names: Column order.
        columns: Column name -> list of values.
    """

    file_count: int
    peptide_count: int
    target_protein_count: int
    decoy_proteins_ignored: int
    i_and_l_equivalent: bool
    sharing_counts: dict[str, int] = field(default_factory=dict)
    files: list[DatabaseFile] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)

    def sharing_of(self) -> dict[str, str]:
        """Peptide -> sharing class. A peptide given twice appears once (it classifies the same).

        Examples:
            >>> calls = classify_peptides(
            ...     ["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
            ...     ["human_subset.xml", "human_extra.fasta"], contaminants=["contaminants.fasta"])
            >>> calls.sharing_of()["YLYEIAR"]
            'SharedAcrossGenes'
        """
        return dict(zip(self.columns.get("peptide", []), self.columns.get("sharing", [])))

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> PeptideClassification:
        return cls(
            file_count=int(data.get("file_count", 0)),
            peptide_count=int(data.get("peptide_count", 0)),
            target_protein_count=int(data.get("target_protein_count", 0)),
            decoy_proteins_ignored=int(data.get("decoy_proteins_ignored", 0)),
            i_and_l_equivalent=bool(data.get("i_and_l_equivalent", True)),
            sharing_counts=dict(data.get("sharing_counts") or {}),
            files=_files(data),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
        )


@dataclass(frozen=True)
class GoRelease:
    """The go.obo a result was computed against. Record it with your results.

    Attributes:
        source_file_name: The file name read.
        sha256: Lower-case hex sha256 of the file's bytes. This, not the release, is what proves two
            runs used the same ontology.
        release: The file's ``data-version``, e.g. ``"releases/2026-07-26"``. ``None``: the file
            has no ``data-version`` header.
        term_count: Terms in the file, obsolete terms included (alternative ids are not counted).
    """

    source_file_name: str
    sha256: str
    release: str | None
    term_count: int

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> GoRelease:
        return cls(
            source_file_name=data.get("source_file_name", ""),
            sha256=data.get("sha256", ""),
            release=data.get("release"),
            term_count=int(data.get("term_count", 0)),
        )


@dataclass(frozen=True)
class AnnotationDatabase:
    """The protein database whose GO terms annotated the groups.

    Attributes:
        path: The absolute path.
        file_type: ``"UniProtXml"``, the only type :func:`annotate_go` accepts.
        reader: The mzLib loader used, ``"ProteinDbLoader.LoadProteinXML"``.
        protein_count: Entries loaded, as written (no decoys generated, no variants applied).
        sha256: Lower-case hex sha256 of the **decompressed** database, so ``.xml`` and
            ``.xml.gz`` agree. Every row's ``annotation_db_sha256`` is this.
        caveats: What mzLib did to the file while loading it, e.g. genotype variants applied.
    """

    path: str
    file_type: str | None
    reader: str | None
    protein_count: int
    sha256: str
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> AnnotationDatabase:
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type"),
            reader=data.get("reader"),
            protein_count=int(data.get("protein_count", 0)),
            sha256=data.get("sha256", ""),
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class WrittenTable:
    """A file :func:`annotate_go` wrote, when asked to.

    Attributes:
        path: The absolute path written.
        row_count: Rows written: the **whole** table, never the ``limit``/``offset`` window.
            ``None`` for the category table, whose rows are counted in :attr:`GoCategories.row_count`.
    """

    path: str
    row_count: int | None = None

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any] | None) -> WrittenTable | None:
        if not data:
            return None
        row_count = data.get("row_count")
        return cls(path=data.get("path", ""), row_count=None if row_count is None else int(row_count))


@dataclass(frozen=True)
class GoCategories(_Table):
    """Your category map applied to the annotated terms: one row per (term, category, subcategory).

    A term belongs to a category when one of that category's anchors is the term itself or one of
    its ancestors. A term under no anchor has **no row**, so absence means "outside your map". Join
    to :class:`GoAnnotations` on ``go_id``: both tables come from the same go.obo release.

    ===============  =====  ===================================================================
    column           type   meaning; ``None`` means
    ===============  =====  ===================================================================
    ``go_id``        str    an annotated term at or below one of the map's anchors
    ``category``     str    the category label, as your map writes it
    ``subcategory``  str    ``"category:subcategory"``, the most specific anchor within the
                            category; ``None``: the term reaches only the category's own anchor
    ===============  =====  ===================================================================

    Attributes:
        map_name: The map's declared ``map_name``.
        map_version: The map's declared ``map_version``. Yours to bump; mzLib only records it.
        source_file_name: The map file read.
        sha256: Lower-case hex sha256 of the map file. A version is a claim; the hash is proof.
        anchor_count: Rows (anchors) in the map.
        row_count: Rows in the table.
        column_names: Column order.
        columns: Column name -> list of values.
    """

    map_name: str
    map_version: str
    source_file_name: str
    sha256: str
    anchor_count: int
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any] | None) -> GoCategories | None:
        if data is None:
            return None
        return cls(
            map_name=data.get("map_name", ""),
            map_version=data.get("map_version", ""),
            source_file_name=data.get("source_file_name", ""),
            sha256=data.get("sha256", ""),
            anchor_count=int(data.get("anchor_count", 0)),
            row_count=int(data.get("row_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
        )


@dataclass(frozen=True)
class GoAnnotations(_Table):
    """What :func:`annotate_go` returns: one row per (protein group, GO term) that **any** member of
    the group holds, directly or through an ancestor.

    No member of a group is privileged - MetaMorpheus never picks a leading protein - so the table
    is the union, and every row says which members carry its term. The views people usually want
    are filters you apply: **consensus** is ``n_with == n_members``; **direct annotations only** is
    ``propagated is False``; **leave out isoform inheritance** is ``inherited is False``.

    The column names are **mzLib's own** ``GoAnnotationTsv`` schema, in its order, so this table and
    the file written by ``out=`` (or by mzLib inside a search) mean the same thing:

    =========================  ==========  ===========================================================
    column                     type        meaning; ``None`` means
    =========================  ==========  ===========================================================
    ``protein_group``          str         the group, as MetaMorpheus named it (``P0C0S5|Q71UI9``)
    ``accession_used``         list[str]   members carrying the term, directly or by propagation
    ``accession_direct``       list[str]   of those, members annotated to this exact term
    ``accession_inherited``    list[str]   of those, members whose terms were borrowed: an isoform
                                           (``P04406-2``) or sequence variant (``P04406_A20T``) absent
                                           from the database takes its entry's terms
    ``go_id``                  str         the term's primary id; ``None``: a term-less row
    ``go_name``                str         the term's name; ``None``: a term-less row
    ``aspect``                 str         ``biological_process``, ``cellular_component``,
                                           ``molecular_function`` or ``unknown``; ``None``: term-less
    ``evidence``               list[str]   evidence codes (ECO) pooled over the carrying members
    ``evidence_by_member``     dict        member -> its own evidence codes for this term
    ``inherited``              bool        every carrying member's terms were borrowed;
                                           ``None``: a term-less row
    ``propagated``             bool        no member is annotated to this exact term - it is implied
                                           by a more specific one; ``None``: a term-less row
    ``n_members``              int         members in the group
    ``n_with``                 int         members carrying the term (``0`` on a term-less row)
    ``entrapment_members``     list[str]   the group's entrapment members, on every row of the group
    ``annotation_status``      str         one of :data:`ANNOTATION_STATUSES`
    ``q_value``                float       the group's q-value from the table (never filtered on)
    ``go_release``             str         the go.obo ``data-version``; ``None``: the file has none
    ``go_obo_sha256``          str         sha256 of the go.obo
    ``annotation_db_sha256``   str         sha256 of the decompressed annotation database
    =========================  ==========  ===========================================================

    Attributes:
        groups_file: The absolute path of the protein-group table read.
        groups_file_sha256: Lower-case hex sha256 of that table, also in :attr:`header` as
            ``source_file_sha256``.
        table_row_count: Rows in the table, decoys included.
        decoy_group_count: Of those, decoy groups, which are skipped: decoys carry no GO.
        group_count: Groups annotated: every non-decoy group, contaminants included, each once even
            if MetaMorpheus wrote it twice.
        annotation_database: The database the terms came from.
        go: The ontology release the terms were propagated against.
        skip_unknown_go_ids: Whether ``skip_unknown_go_ids=True`` was given.
        unresolved_go_ids: GO ids the database cites that the ontology release lacks, dropped
            because ``skip_unknown_go_ids=True``. Always empty otherwise (the call fails instead).
        header: mzLib's ``#!key value`` header, as its writer produced it: the format version, the
            mzLib build, the release and every sha256, and five counters - ``n_multi_member_groups``
            and one ``status_<status>`` per status. **The counters count groups, not rows, and only
            groups at q <= ``counter_q_value_max`` (0.01)**; the rows themselves are not filtered.
            Every value is a string, as in the file.
        row_count: Rows in the whole table.
        returned_count: Rows returned here (the ``limit``/``offset`` window).
        offset: Rows skipped before the window: the ``offset`` applied.
        truncated: True when rows were left out of :attr:`columns` by ``limit`` or ``offset``.
        caveats: What this result cannot tell you: ids dropped from the database, a go.obo with no
            release, what mzLib did loading the database.
        written: Where the whole table was written, when ``out=`` was given.
        categories: Your category map applied to these terms, when ``category_map=`` was given.
        categories_written: Where the category table was written, when ``categories_out=`` was given.
        column_names: Column order: mzLib's schema.
        columns: Column name -> list of values, for the returned window.
    """

    groups_file: str
    groups_file_sha256: str
    table_row_count: int
    decoy_group_count: int
    group_count: int
    annotation_database: AnnotationDatabase
    go: GoRelease
    skip_unknown_go_ids: bool
    row_count: int
    returned_count: int
    offset: int = 0
    truncated: bool = False
    unresolved_go_ids: list[str] = field(default_factory=list)
    header: dict[str, str] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)
    written: WrittenTable | None = None
    categories: GoCategories | None = None
    categories_written: WrittenTable | None = None
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> GoAnnotations:
        return cls(
            groups_file=data.get("groups_file", ""),
            groups_file_sha256=data.get("groups_file_sha256", ""),
            table_row_count=int(data.get("table_row_count", 0)),
            decoy_group_count=int(data.get("decoy_group_count", 0)),
            group_count=int(data.get("group_count", 0)),
            annotation_database=AnnotationDatabase._from_wire(data.get("annotation_database") or {}),
            go=GoRelease._from_wire(data.get("go") or {}),
            skip_unknown_go_ids=bool(data.get("skip_unknown_go_ids", False)),
            row_count=int(data.get("row_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            unresolved_go_ids=list(data.get("unresolved_go_ids") or []),
            header=dict(data.get("header") or {}),
            caveats=list(data.get("caveats") or []),
            written=WrittenTable._from_wire(data.get("written")),
            categories=GoCategories._from_wire(data.get("categories")),
            categories_written=WrittenTable._from_wire(data.get("categories_written")),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
        )


@dataclass(frozen=True)
class GoUpdate:
    """What :func:`update_go` returns: which go.obo is now on disk, and whether it changed.

    Attributes:
        go_obo_file: The absolute path written.
        url: Where it was fetched from: GO's PURL, which always serves the **current** release.
        existed_before: Whether a file was already at the path.
        previous_sha256: That file's sha256; ``None``: there was none.
        changed: Whether the file on disk is now different. When it is and a file existed, the old
            one was kept beside it as ``go.obo.<yyyyMMdd-HHmmss-fff>`` (see :attr:`caveats`).
        go: The release now at the path.
        caveats: The backup kept, or a file with no ``data-version``.
    """

    go_obo_file: str
    url: str
    existed_before: bool
    previous_sha256: str | None
    changed: bool
    go: GoRelease
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> GoUpdate:
        return cls(
            go_obo_file=data.get("go_obo_file", ""),
            url=data.get("url", ""),
            existed_before=bool(data.get("existed_before", False)),
            previous_sha256=data.get("previous_sha256"),
            changed=bool(data.get("changed", False)),
            go=GoRelease._from_wire(data.get("go") or {}),
            caveats=list(data.get("caveats") or []),
        )


def _check_position(position: object, length: int, what: str) -> int:
    """A 1-based position in a sequence of ``length`` residues, or a :class:`UsageError`."""
    if isinstance(position, bool) or not isinstance(position, int):
        raise _bridge.UsageError(f"position must be an int; got {position!r}.")
    if not 1 <= position <= length:
        raise _bridge.UsageError(
            f"position {position} is outside the {what}, which has positions 1 to {length} "
            "(positions are 1-based, as in mzLib and UniProt)."
        )
    return position


@dataclass(frozen=True)
class Alignment:
    """One pair's alignment, from :func:`align` or :attr:`Alignments.alignments`.

    Positions are **1-based**, as in mzLib and UniProt. A residue facing a gap has no counterpart,
    and :meth:`map_source_to_target` returns ``None`` for it - never the nearest residue.

    Attributes:
        pair_index: Position of the pair in the input, from 0.
        aligner_id: mzLib's name for the method and every parameter, e.g.
            ``"global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"``. Store it with any position
            you map, so the mapping says how it was made.
        score: Substitution scores summed over the paired residues, less every charged gap, in
            BLOSUM62's half-bits. A raw score: it grows with length, and is not a bit score, an
            E-value or a percent identity.
        identical: Paired residues that are the same letter (case-insensitive).
        substituted: Paired residues that differ.
        source_residues_facing_gap: Source residues with no counterpart in the target.
        target_residues_facing_gap: Target residues with no counterpart in the source.
        source_length: The source's length, in residues.
        target_length: The target's length, in residues.
        aligned_source: The source with ``-`` wherever it faces a target residue across a gap.
        aligned_target: The target with ``-`` wherever it faces a source residue across a gap. As
            long as :attr:`aligned_source`; print the two one above the other to see the alignment.
        source_to_target: For each source position ``p``, item ``p - 1`` is the paired target
            position, or ``None`` when that residue faces a gap. ``None`` as a whole when the call
            asked for no residue map.
        target_to_source: The same, from each target position to the source.
        caveats: What the score and the mapping do and do not mean (the call's ``caveats``).
    """

    pair_index: int
    aligner_id: str
    score: int
    identical: int
    substituted: int
    source_residues_facing_gap: int
    target_residues_facing_gap: int
    source_length: int
    target_length: int
    aligned_source: str
    aligned_target: str
    source_to_target: tuple[int | None, ...] | None = None
    target_to_source: tuple[int | None, ...] | None = None
    caveats: list[str] = field(default_factory=list)

    def map_source_to_target(self, position: int) -> int | None:
        """The target position paired with source ``position``, or ``None`` if it faces a gap.

        Wraps nothing new: it reads the residue map mzLib's ``TryMapSourceToTarget`` filled.

        Args:
            position: A 1-based source position, 1 to :attr:`source_length`.

        Returns:
            The 1-based target position, or ``None`` when the source residue has no counterpart.

        Raises:
            UsageError: ``position`` is not an int in range, or the call asked for no residue map.

        Examples:
            >>> aln = align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")
            >>> aln.map_source_to_target(12), aln.map_source_to_target(10)
            (9, None)
        """
        return self._lookup(self.source_to_target, position, self.source_length, "source")

    def map_target_to_source(self, position: int) -> int | None:
        """The source position paired with target ``position``, or ``None`` if it faces a gap.

        Args:
            position: A 1-based target position, 1 to :attr:`target_length`.

        Returns:
            The 1-based source position, or ``None`` when the target residue has no counterpart.

        Raises:
            UsageError: ``position`` is not an int in range, or the call asked for no residue map.

        Examples:
            >>> aln = align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")
            >>> aln.map_target_to_source(9)
            12
        """
        return self._lookup(self.target_to_source, position, self.target_length, "target")

    @staticmethod
    def _lookup(table: tuple[int | None, ...] | None, position: int, length: int, what: str) -> int | None:
        if table is None:
            raise _bridge.UsageError(
                "This alignment has no residue map: it came from a call with residue_map=False."
            )
        return table[_check_position(position, length, what) - 1]


@dataclass(frozen=True)
class Alignments(_Table):
    """What :func:`align_many` returns: one row per pair, and every residue's counterpart.

    The pairs table (:attr:`columns`), one row per input pair, in input order:

    ==============================  ===========  ===================================================
    column                          type / unit  meaning
    ==============================  ===========  ===================================================
    ``pair_index``                  int          position of the pair in the input, from 0
    ``source_length``               residues     the source's length
    ``target_length``               residues     the target's length
    ``score``                       half-bits    BLOSUM62 score of the alignment (see the caveats)
    ``identical``                   residues     paired residues that are the same letter
    ``substituted``                 residues     paired residues that differ
    ``source_residues_facing_gap``  residues     source residues with no counterpart
    ``target_residues_facing_gap``  residues     target residues with no counterpart
    ``aligned_source``              str          the source, ``-`` where it faces a gap
    ``aligned_target``              str          the target, ``-`` where it faces a gap
    ==============================  ===========  ===================================================

    The residue map (:attr:`residue_map`), one row per residue of each sequence, pair by pair: a
    pair's source residues first, then its target residues:

    ===================  =====  =============================================================
    column               type   meaning; ``None`` means
    ===================  =====  =============================================================
    ``pair_index``       int    the pair
    ``sequence``         str    ``"source"`` or ``"target"``: which sequence ``position`` is in
    ``position``         int    1-based position in that sequence
    ``residue``          str    the residue there, as given
    ``mapped_position``  int    1-based position of its counterpart in the other sequence;
                                ``None``: it faces a gap (**never** the nearest residue)
    ``mapped_residue``   str    the counterpart, as given; ``None``: it faces a gap. A
                                different letter from ``residue`` is a substitution
    ===================  =====  =============================================================

    Attributes:
        aligner_id: mzLib's name for the method and every parameter, e.g.
            ``"global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"``. One aligner per call.
        matrix: The substitution matrix, always ``"BLOSUM62"``: mzLib's default and the only one
            offered.
        gap_open: The gap-open cost used, in BLOSUM62's half-bits: what you passed, or mzLib's
            default (11) when you passed nothing.
        gap_extend: The cost per gap residue used, in half-bits: yours, or mzLib's default (1).
        end_gaps: ``"free"`` or ``"charged"``: the end-gap rule used, yours or mzLib's default
            (``"free"``).
        max_cells: The table limit used, in cells: yours, or mzLib's default (200,000,000).
        record_count: Pairs aligned: the rows of :attr:`columns`.
        column_names: The pairs table's column order.
        columns: Column name -> one value per pair. Hand it to ``pandas.DataFrame(...)``.
        residue_map: The residue map, a :class:`Table`; ``None`` when ``residue_map=False``.
        caveats: What the score and the mapping do and do not mean.
    """

    aligner_id: str
    matrix: str
    gap_open: int
    gap_extend: int
    end_gaps: str
    max_cells: int
    record_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    residue_map: Table | None = None
    caveats: list[str] = field(default_factory=list)

    @property
    def alignments(self) -> list[Alignment]:
        """One :class:`Alignment` per pair, in input order, each with its residue map.

        Examples:
            >>> alns = align_many([("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"),
            ...                    ("MLLLPEPTIDE", "PEPTIDE")])
            >>> [a.map_target_to_source(1) for a in alns.alignments]
            [1, 5]
        """
        maps: dict[tuple[int, str], dict[int, int | None]] = {}
        if self.residue_map is not None:
            cols = self.residue_map.columns
            for pair, sequence, position, mapped in zip(
                cols.get("pair_index", []), cols.get("sequence", []),
                cols.get("position", []), cols.get("mapped_position", []),
            ):
                maps.setdefault((pair, sequence), {})[position] = mapped

        def table(pair: int, sequence: str, length: int) -> tuple[int | None, ...] | None:
            if self.residue_map is None:
                return None
            found = maps.get((pair, sequence), {})
            return tuple(found.get(p) for p in range(1, length + 1))

        out = []
        for row in self.records:
            pair = int(row["pair_index"])
            out.append(Alignment(
                pair_index=pair,
                aligner_id=self.aligner_id,
                score=int(row["score"]),
                identical=int(row["identical"]),
                substituted=int(row["substituted"]),
                source_residues_facing_gap=int(row["source_residues_facing_gap"]),
                target_residues_facing_gap=int(row["target_residues_facing_gap"]),
                source_length=int(row["source_length"]),
                target_length=int(row["target_length"]),
                aligned_source=row["aligned_source"],
                aligned_target=row["aligned_target"],
                source_to_target=table(pair, "source", int(row["source_length"])),
                target_to_source=table(pair, "target", int(row["target_length"])),
                caveats=list(self.caveats),
            ))
        return out

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> Alignments:
        return cls(
            aligner_id=data.get("aligner_id", ""),
            matrix=data.get("matrix", ""),
            gap_open=int(data.get("gap_open", 0)),
            gap_extend=int(data.get("gap_extend", 0)),
            end_gaps=data.get("end_gaps", ""),
            max_cells=int(data.get("max_cells", 0)),
            record_count=int(data.get("record_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
            residue_map=Table._from_wire(data.get("residue_map")),
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class ParsedAccessions(_Table):
    """What :func:`parse_accessions` returns: one row per accession, in the order given.

    ============================  ====  ==============================================================
    column                        type  meaning; ``None`` means
    ============================  ====  ==============================================================
    ``input``                     str   exactly as given
    ``entry_accession``           str   the UniProt entry or unversioned RefSeq accession
                                        (``P04637`` for ``P04637-2_S70N``); ``input`` when
                                        unrecognized. An **entry**, not a sequence
    ``isoform``                   int   UniProt ``-n`` isoform; ``None``: no suffix (always for
                                        RefSeq), or a number past 2,147,483,647
    ``version``                   int   RefSeq ``.n`` version; ``None``: no version (always for
                                        UniProt), or a number past 2,147,483,647
    ``namespace``                 str   one of :data:`NAMESPACES`
    ``parent_accession``          str   the entry the variants were applied to, isoform or version
                                        kept (``P04637-2`` for ``P04637-2_Q5*``); ``input`` when
                                        there are no applied variants
    ``applied_variants``          str   the variant suffix as written (``S70N_A80T``); ``None``: none
    ``has_applied_variants``      bool  ``applied_variants`` is not ``None``
    ============================  ====  ==============================================================

    Attributes:
        record_count: Accessions parsed: the rows, duplicates included.
        namespace_counts: Namespace -> accessions in it. Every key in :data:`NAMESPACES` is
            present, zeros included. A large ``"unrecognized"`` count is worth reading: see
            :attr:`caveats`.
        with_applied_variants_count: Accessions that name applied variants.
        column_names: Column order.
        columns: Column name -> one value per accession.
        caveats: What mzLib leaves unrecognized rather than repairs, and other traps.
    """

    record_count: int
    namespace_counts: dict[str, int] = field(default_factory=dict)
    with_applied_variants_count: int = 0
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    def entry_of(self) -> dict[str, str]:
        """Accession -> entry accession, to group isoforms and variants under their entry.

        An accession given twice appears once (it parses the same). An unrecognized one maps to
        itself.

        Examples:
            >>> parsed = parse_accessions(["P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P",
            ...     "P04637_S70N_A80T", "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"])
            >>> parsed.entry_of()["NP_000537.3_R72P"]
            'NP_000537'
        """
        return dict(zip(self.columns.get("input", []), self.columns.get("entry_accession", [])))

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> ParsedAccessions:
        return cls(
            record_count=int(data.get("record_count", 0)),
            namespace_counts=dict(data.get("namespace_counts") or {}),
            with_applied_variants_count=int(data.get("with_applied_variants_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class GeneTrees:
    """Compara's gene-tree file a snapshot was built with, and what was kept from it.

    Attributes:
        source_file_name: The file's name, e.g. ``vertebrates.GeneTree_content.default.e116.txt.gz``.
        sha256: Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a
            ``.gz``).
        release: The release in the file name (``.e116.txt``); ``None`` when the name carries none,
            and then mzLib could not hold the file to the release you gave - a caveat says so.
        collection: The Compara collection in the file name, e.g. ``"vertebrates"``; ``None`` when
            the name does not start ``<collection>.GeneTree_content.``.
        member_count: Genes held in memory from the file, in genes: every gene in it, or with
            ``restrict_gene_trees=True`` only those in your gene sets.
        restricted_to_gene_sets: The sha256s of the gene sets the members were restricted to;
            ``None`` when every gene in the file was kept (``restrict_gene_trees=False``).
    """

    source_file_name: str
    sha256: str
    release: str | None
    collection: str | None
    member_count: int
    restricted_to_gene_sets: list[str] | None = None

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> GeneTrees:
        restricted = data.get("restricted_to_gene_sets")
        return cls(
            source_file_name=data.get("source_file_name", ""),
            sha256=data.get("sha256", ""),
            release=data.get("release"),
            collection=data.get("collection"),
            member_count=int(data.get("member_count", 0)),
            restricted_to_gene_sets=list(restricted) if restricted is not None else None,
        )


@dataclass(frozen=True)
class HomologyDump:
    """One species' Compara homology dump, as read.

    Attributes:
        source_file_name: The file's name. Ensembl names every species' dump
            ``Compara.<release>.protein_default.homologies.tsv.gz``, so this alone does not say whose
            it is: the key of :attr:`OrthologySnapshot.dumps` does, read from the dump's rows.
        sha256: Lower-case hex sha256 of the file's bytes as read.
        release: The release in the file name; ``None`` when it carries none (then not checked).
        rows_read: Every data row read and checked, kept or not, in rows.
        rows_kept: Rows whose two species are both snapshot species, in rows.
    """

    source_file_name: str
    sha256: str
    release: str | None
    rows_read: int
    rows_kept: int

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> HomologyDump:
        return cls(
            source_file_name=data.get("source_file_name", ""),
            sha256=data.get("sha256", ""),
            release=data.get("release"),
            rows_read=int(data.get("rows_read", 0)),
            rows_kept=int(data.get("rows_kept", 0)),
        )


@dataclass(frozen=True)
class SpeciesPair:
    """How many relationships one species pair holds: one row each in ``pairs/<a>__<b>.parquet``.

    Attributes:
        species_a: The species that sorts first: side A of every row of the pair.
        species_b: The other; the same species for its paralogs, the pair ``(s, s)``.
        row_count: Relationships between the two, in rows. One row each, never merged or picked.
        ortholog_count: Of those, orthologs (``ortholog_one2one``, ``ortholog_one2many``,
            ``ortholog_many2many``), in rows. Always 0 for ``(s, s)``.
        paralog_count: Paralogs (``within_species_paralog``, ``other_paralog``, ``gene_split``),
            in rows. Always 0 between two species.
        homoeolog_count: Homoeologs, between subgenomes of one polyploid species, in rows.
    """

    species_a: str
    species_b: str
    row_count: int
    ortholog_count: int
    paralog_count: int
    homoeolog_count: int

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> SpeciesPair:
        return cls(
            species_a=data.get("species_a", ""),
            species_b=data.get("species_b", ""),
            row_count=int(data.get("row_count", 0)),
            ortholog_count=int(data.get("ortholog_count", 0)),
            paralog_count=int(data.get("paralog_count", 0)),
            homoeolog_count=int(data.get("homoeolog_count", 0)),
        )


@dataclass(frozen=True)
class SnapshotFile:
    """One file of a written snapshot.

    Attributes:
        path: Relative to :attr:`WrittenSnapshot.directory`, with forward slashes, e.g.
            ``"pairs/homo_sapiens__mus_musculus.parquet"``.
        rows: Rows in the Parquet table; ``None`` for a file that is not a table (``views.sql``).
            A pair with no relationship is still written, with 0 rows.
        bytes: The file's size in bytes.
        sha256: Lower-case hex sha256 of the file, as ``manifest.json`` lists it.
    """

    path: str
    rows: int | None
    bytes: int
    sha256: str

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> SnapshotFile:
        rows = data.get("rows")
        return cls(
            path=data.get("path", ""),
            rows=int(rows) if rows is not None else None,
            bytes=int(data.get("bytes", 0)),
            sha256=data.get("sha256", ""),
        )


@dataclass(frozen=True)
class WrittenSnapshot:
    """Where :func:`build_orthology` wrote the snapshot, with ``out=``.

    Attributes:
        directory: The ``out`` directory, absolute.
        files: Every file written but ``manifest.json``, in order of path: the Parquet tables
            ``genes/<species>.parquet``, ``members/<species>.parquet`` and
            ``pairs/<a>__<b>.parquet``, and ``views.sql``.
    """

    directory: str
    files: list[SnapshotFile] = field(default_factory=list)

    @property
    def manifest(self) -> str:
        """The path of ``manifest.json``, which lists every input and file with its sha256."""
        return os.path.join(self.directory, "manifest.json")

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any] | None) -> WrittenSnapshot | None:
        if data is None:
            return None
        return cls(
            directory=data.get("directory", ""),
            files=[SnapshotFile._from_wire(f) for f in data.get("files") or []],
        )


@dataclass(frozen=True)
class PairStatus(_Table):
    """Every gene of one species with how far Compara's inference got in another.

    One row per gene of :attr:`from_species`'s gene set, every biotype, in order of gene id:

    ================  ====  ==================================================================
    column            type  meaning; ``None`` means
    ================  ====  ==================================================================
    ``gene_id``       str   stable Ensembl gene id, unversioned
    ``gene_name``     str   the release's gene name, a label; ``None``: the GTF gives none
    ``gene_biotype``  str   the release's biotype (``unknown`` when the GTF omits it)
    ``status``        str   one of :data:`ORTHOLOGY_STATUSES` but ``not_in_gene_set``
    ================  ====  ==================================================================

    Attributes:
        from_species: Whose genes the rows are.
        to_species: The species each status is about.
        row_count: Genes: every gene of ``from_species``.
        status_counts: Status -> genes with it. Every key in :data:`ORTHOLOGY_STATUSES` is
            present, zeros included; ``not_in_gene_set`` is always 0 here.
        column_names: Column order.
        columns: Column name -> one value per gene.
    """

    from_species: str
    to_species: str
    row_count: int
    status_counts: dict[str, int] = field(default_factory=dict)
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)

    def genes(self, status: str) -> list[str]:
        """The gene ids with one status, in order of gene id.

        Raises:
            UsageError: ``status`` is not one of :data:`ORTHOLOGY_STATUSES`.

        See :func:`build_orthology` for an example.
        """
        if status not in ORTHOLOGY_STATUSES:
            raise _bridge.UsageError(
                f"status must be one of {', '.join(ORTHOLOGY_STATUSES)}; got {status!r}."
            )
        ids = self.columns.get("gene_id", [])
        return [g for g, s in zip(ids, self.columns.get("status", [])) if s == status]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any] | None) -> PairStatus | None:
        if data is None:
            return None
        return cls(
            from_species=data.get("from_species", ""),
            to_species=data.get("to_species", ""),
            row_count=int(data.get("row_count", 0)),
            status_counts=dict(data.get("status_counts") or {}),
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
        )


@dataclass(frozen=True)
class OrthologySnapshot:
    """What :func:`build_orthology` returns: an Ensembl Compara orthology snapshot, as built.

    The relationships themselves are in the Parquet files ``out=`` writes, not here: this says what
    went in, what came out, and, with ``pair_status=``, every gene's status in one other species.

    Attributes:
        release: The Ensembl release you gave.
        species: The snapshot's species, Ensembl production names, in sorted order.
        format: ``"ensembl-orthology-snapshot"``, as ``manifest.json`` names the format.
        format_version: The format's version (1 at mzLib 1.0.595).
        snapshot_id: mzLib's id for the snapshot: a sha256 over the format, the release, the
            species and every input's sha256. The same inputs give the same id, written or not,
            and it is ``manifest.json``'s ``snapshot_id``. Record it with any result drawn from the
            snapshot.
        identical_duplicates_dropped: Homology rows found in both of a pair's dumps, identical, and
            counted once, in rows. 0 in release 116, whose dumps partition.
        gene_sets: Species -> the :class:`GeneSet` it was counted against, the same type
            :func:`resolve_genes` reports. Its ``sha256`` is what ``manifest.json``'s
            ``gene_set_sha256`` and :func:`resolve_genes`'s ``gene_set_sha256`` column carry, so a
            resolution and a snapshot can be checked against one gene set.
        gene_trees: The gene-tree file, as :class:`GeneTrees`.
        dumps: Species -> its :class:`HomologyDump`. The species is read from the dump's rows.
        pairs: One :class:`SpeciesPair` per species pair ``(a, b)``, ``a <= b``, ``(s, s)``
            included. :meth:`pair` finds one in either order.
        written: Where the snapshot was written, as :class:`WrittenSnapshot`; ``None`` without
            ``out=``, when nothing was written.
        pair_status: Every gene's status, as :class:`PairStatus`; ``None`` without
            ``pair_status=``.
        caveats: This call's traps: a file name with no release (not checked against it), dumps
            sharing one file name, a GTF named for another release, duplicate rows, unrestricted
            gene trees.
    """

    release: str
    species: list[str]
    format: str
    format_version: int
    snapshot_id: str
    identical_duplicates_dropped: int
    gene_sets: dict[str, GeneSet] = field(default_factory=dict)
    gene_trees: GeneTrees | None = None
    dumps: dict[str, HomologyDump] = field(default_factory=dict)
    pairs: list[SpeciesPair] = field(default_factory=list)
    written: WrittenSnapshot | None = None
    pair_status: PairStatus | None = None
    caveats: list[str] = field(default_factory=list)

    def pair(self, a: str, b: str) -> SpeciesPair:
        """The counts for species ``a`` and ``b``, given in either order.

        Raises:
            UsageError: ``a`` or ``b`` is not one of :attr:`species`.

        See :func:`build_orthology` for an example.
        """
        first, second = sorted((a, b))
        for pair in self.pairs:
            if (pair.species_a, pair.species_b) == (first, second):
                return pair
        raise _bridge.UsageError(
            f"({a!r}, {b!r}) is not a pair of this snapshot's species: {', '.join(self.species)}."
        )

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> OrthologySnapshot:
        return cls(
            release=data.get("release", ""),
            species=list(data.get("species") or []),
            format=data.get("format", ""),
            format_version=int(data.get("format_version", 0)),
            snapshot_id=data.get("snapshot_id", ""),
            identical_duplicates_dropped=int(data.get("identical_duplicates_dropped", 0)),
            gene_sets={s: GeneSet._from_wire(g) for s, g in (data.get("gene_sets") or {}).items()},
            gene_trees=GeneTrees._from_wire(data.get("gene_trees") or {}),
            dumps={s: HomologyDump._from_wire(d) for s, d in (data.get("dumps") or {}).items()},
            pairs=[SpeciesPair._from_wire(p) for p in data.get("pairs") or []],
            written=WrittenSnapshot._from_wire(data.get("written")),
            pair_status=PairStatus._from_wire(data.get("pair_status")),
            caveats=list(data.get("caveats") or []),
        )


# ---- argument shaping --------------------------------------------------------------------------


def _paths(value: _PathLike | Iterable[_PathLike] | None, what: str) -> list[str]:
    """One path or many, as clean strings. A bare string is one path, never its characters."""
    if value is None:
        return []
    items = [value] if isinstance(value, (str, os.PathLike)) else list(value)
    paths = []
    for item in items:
        text = _bridge.path_text(item)
        if not text:
            raise _bridge.UsageError(f"Every {what} path must be a non-empty str or path; got {item!r}.")
        if "\t" in text or "\n" in text or "\r" in text:
            raise _bridge.UsageError(f"A {what} path contains a tab or newline: {text!r}.")
        if text == _SECTION:
            raise _bridge.UsageError(f"'{_SECTION}' is not a usable {what} path.")
        paths.append(text)
    return paths


def _lines(values: Iterable[str], what: str) -> list[str]:
    """A list of accessions or peptides for stdin, one per line."""
    if isinstance(values, str):
        raise _bridge.UsageError(
            f"{what} must be a list of strings, not one string - that would send one {what[:-1]} "
            f"per character. Wrap it: [{values!r}]."
        )
    out = []
    for value in values:
        if not isinstance(value, str):
            raise _bridge.UsageError(f"Every entry in {what} must be a str; got {value!r}.")
        text = value.strip()
        if not text:
            raise _bridge.UsageError(f"{what} contains a blank entry.")
        if "\t" in text or "\n" in text or "\r" in text:
            raise _bridge.UsageError(f"An entry in {what} contains a tab or newline: {value!r}.")
        out.append(text)
    return out


def _threads(threads: int) -> list[str]:
    if isinstance(threads, bool) or not isinstance(threads, int) or threads == 0 or threads < -1:
        raise _bridge.UsageError(
            f"threads must be 1 or more, or -1 for every core; got {threads!r}."
        )
    return ["--threads", str(threads)]


def _on_error(on_error: str) -> list[str]:
    if on_error not in ("fail", "skip"):
        raise _bridge.UsageError(f"on_error must be 'fail' or 'skip'; got {on_error!r}.")
    return ["--on-error", on_error]


def _database_args(
    databases: _PathLike | Iterable[_PathLike] | None,
    contaminants: _PathLike | Iterable[_PathLike] | None,
) -> tuple[list[str], list[str] | None]:
    """The database options, and the stdin path lines when there are several databases.

    One database travels as ``--path`` (plus ``--contaminant``); several go on stdin, one per line,
    a contaminant's line ending in a tab and ``contaminant``.
    """
    targets = _paths(databases, "database")
    extra = _paths(contaminants, "contaminant database")
    if not targets and not extra:
        raise _bridge.UsageError("At least one protein database is required, e.g. 'human.xml'.")

    everything = targets + extra
    if len({os.path.normcase(os.path.abspath(p)) for p in everything}) != len(everything):
        raise _bridge.UsageError("A database is listed twice.")

    if len(everything) == 1:
        args = ["--path", everything[0]]
        if extra:
            args.append("--contaminant")
        return args, None

    lines = targets + [f"{p}\tcontaminant" for p in extra]
    return ["--paths-stdin"], lines


def _stdin(path_lines: list[str] | None, second: list[str] | None) -> str | None:
    """Joins the path lines and a second list, separated by a ``--`` line when both are present."""
    if path_lines is None and second is None:
        return None
    if path_lines is None:
        return "\n".join(second or []) + "\n"
    if second is None:
        return "\n".join(path_lines) + "\n"
    return "\n".join([*path_lines, _SECTION, *second]) + "\n"


# ---- the verbs ---------------------------------------------------------------------------------


def read(
    databases: _PathLike | Iterable[_PathLike],
    *,
    tables: Sequence[str] = ("proteins",),
    accessions: Iterable[str] | None = None,
    contaminants: _PathLike | Iterable[_PathLike] | None = None,
    sequences: bool = False,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ProteinDatabase:
    """Read protein databases: one row per protein, and GO and Ensembl tables on request.

    Loads each database with mzLib's ``ProteinDbLoader``, exactly as a search would but with no
    decoys generated. A UniProt XML carries everything; a FASTA carries organism, taxon (``OX=``)
    and gene name (``GN=``) only - see :attr:`DatabaseFile.absent_fields`.

    Args:
        databases: A UniProt XML (``.xml``) or FASTA (``.fasta``, ``.fa``, ``.faa``, ``.fas``), each
            optionally gzipped, or a list of them. Read in the order given.
        tables: Which tables to return, from :data:`TABLES`. Default ``("proteins",)``. Ask for
            ``"go_terms"`` and ``"ensembl_genes"`` when you want them: a whole human proteome has
            about twenty GO rows per protein.
        accessions: Keep only proteins whose accession is one of these, **exact match**. Applies to
            every table. Misses are listed in :attr:`ProteinDatabase.accessions_not_found`. Sent on
            stdin, so tens of thousands are fine.
        contaminants: Databases to load as contaminants; ``is_contaminant`` is then ``True`` on
            their rows. They come after ``databases`` in ``source_index`` order.
        sequences: Add a ``sequence`` column to the proteins table.
        threads: Databases read at once. Default 1; ``-1`` means every core. A resource choice
            only: the result is identical at any value.
        on_error: ``"fail"`` (default) raises on the first unreadable database, in input order;
            ``"skip"`` records the failure in that file's :attr:`DatabaseFile.error` and reads the
            rest.
        timeout: Seconds to allow. ``None`` (default) waits: a proteome XML takes a while.

    Returns:
        A :class:`ProteinDatabase`.

    Raises:
        UsageError: no database, a blank or repeated path, an unknown table, a bad ``threads`` or
            ``on_error``, an empty ``accessions`` list; or (under ``on_error="fail"``) a missing
            file or one that is neither XML nor FASTA by name.
        BridgeError: mzLib could not parse a database (under ``on_error="fail"``).

    Examples:
        >>> db = read(["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
        ...           contaminants=["contaminants.fasta"], tables=TABLES)
        >>> db.record_count, db.go_terms.row_count, db.ensembl_genes.row_count
        (8, 47, 6)
        >>> db.files[1].file_type, db.files[1].absent_fields
        ('Fasta', ['go_terms', 'ensembl_genes', 'ensembl_gene_ids'])
    """
    if isinstance(tables, str):
        tables = (tables,)
    tables = list(tables)
    unknown = [t for t in tables if t not in TABLES]
    if not tables or unknown:
        raise _bridge.UsageError(f"tables must be a non-empty selection from {TABLES}; got {tables!r}.")

    db_args, path_lines = _database_args(databases, contaminants)
    args = ["proteins", "read", *db_args, "--tables", ",".join(tables), *_threads(threads), *_on_error(on_error)]
    if sequences:
        args.append("--sequences")

    filter_lines = None
    if accessions is not None:
        filter_lines = _lines(accessions, "accessions")
        if not filter_lines:
            raise _bridge.UsageError("accessions is empty; pass None to read every protein.")
        args.append("--accessions-stdin")

    data = _bridge.invoke(*args, stdin=_stdin(path_lines, filter_lines), timeout=timeout)
    return ProteinDatabase._from_wire(data)


def resolve_genes(
    databases: _PathLike | Iterable[_PathLike],
    *,
    gtf: _PathLike | None = None,
    gene_set: _PathLike | None = None,
    xref: _PathLike | None = None,
    contaminants: _PathLike | Iterable[_PathLike] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> GeneResolutions:
    """Resolve every protein to stable Ensembl gene ids, against a gene set you pin.

    Wraps mzLib's ``EnsemblGeneResolver``. The gene links come from the database itself (UniProt's
    ``<dbReference type="Ensembl">``), and each is **counted against the gene set you pass**: a gene
    in the set counts; a gene outside it (an ALT haplotype, a patch) is dropped from ``n_genes`` but
    counted in ``off_primary_genes``, and a protein with only such genes is ``off_primary_only``.

    **You supply the gene set, and you should pin it.** Gene ids, versions, symbols and the very
    set of genes on the primary assembly change between Ensembl releases, so a resolution means
    something only relative to one release. Use Ensembl's **primary-assembly** GTF (
    ``Species.Assembly.Release.gtf.gz``, not the ``chr_patch_hapl_scaff`` one): counting ALT
    haplotypes as separate genes made human release 116 look 6.99% multi-gene when it is 0.36%.
    Keep Ensembl's file name so the release is recorded, and record :attr:`GeneSet.sha256` with your
    results. Nothing is downloaded or defaulted here.

    Args:
        databases: The search database(s), UniProt XML for real answers. A FASTA carries no gene
            links, so every FASTA protein is ``not_in_source`` - a caveat says so.
        gtf: An Ensembl GTF, plain or ``.gz``. Only its ``gene`` rows are read.
        gene_set: Instead of ``gtf``, a compact gene table written by mzLib's
            ``EnsemblGeneSetWriter`` (about 0.5 MB for human, against a 141 MB GTF). It carries the
            GTF's provenance, so rows are keyed exactly as against the GTF. Give exactly one.
        xref: Optional. Ensembl's ``Species.Assembly.Release.uniprot.tsv.gz``, for a second
            opinion: each gene row then says whether Ensembl agrees, and a gene only Ensembl links
            gets its own row (``source == "ensembl_xref"``). Without it agreement is ``None`` -
            unknown, not false.
        contaminants: Databases to load as contaminants. Their proteins are
            ``contaminant_not_mapped`` - never mapped, and never silently dropped.
        threads: Databases read at once. Default 1; ``-1`` means every core. Same answer at any
            value.
        on_error: ``"fail"`` (default) or ``"skip"``, as for :func:`read`.
        timeout: Seconds to allow. ``None`` (default) waits: a gzipped GTF takes tens of seconds.

    Returns:
        A :class:`GeneResolutions`.

    Raises:
        UsageError: neither or both of ``gtf`` and ``gene_set``; a missing GTF, gene set, xref or
            database; any argument :func:`read` would refuse.
        BridgeError: a GTF gene row without ``gene_id``, an xref table with unexpected columns, or a
            database mzLib cannot parse.

    Examples:
        >>> genes = resolve_genes(
        ...     ["human_subset.xml", "human_extra.fasta"], contaminants=["contaminants.fasta"],
        ...     gtf="Homo_sapiens.GRCh38.116.gtf", xref="Homo_sapiens.GRCh38.116.uniprot.tsv")
        >>> genes.gene_set.release, genes.gene_set.genome_build
        ('116', 'GRCh38.p14')
        >>> genes.outcome_counts["contaminant_not_mapped"]
        1
    """
    gtf_text = _bridge.path_text(gtf) if gtf is not None else None
    set_text = _bridge.path_text(gene_set) if gene_set is not None else None
    if (gtf is None) == (gene_set is None) or gtf_text == "" or set_text == "":
        raise _bridge.UsageError(
            "Give exactly one of gtf= (an Ensembl GTF, e.g. 'Homo_sapiens.GRCh38.116.gtf.gz') or "
            "gene_set= (a compact table made from one). There is no default gene set: the release "
            "it pins is part of the answer."
        )

    db_args, path_lines = _database_args(databases, contaminants)
    args = ["genes", "resolve", *db_args, *_threads(threads), *_on_error(on_error)]
    args += ["--gtf", gtf_text] if gtf_text else ["--gene-set", set_text or ""]
    if xref is not None:
        xref_text = _bridge.path_text(xref)
        if not xref_text:
            raise _bridge.UsageError("xref must be a non-empty path or None.")
        args += ["--xref", xref_text]

    data = _bridge.invoke(*args, stdin=_stdin(path_lines, None), timeout=timeout)
    return GeneResolutions._from_wire(data)


def classify_peptides(
    peptides: Iterable[str],
    databases: _PathLike | Iterable[_PathLike],
    *,
    contaminants: _PathLike | Iterable[_PathLike] | None = None,
    threads: int = 1,
    timeout: float | None = None,
) -> PeptideClassification:
    """Classify peptides by how widely they are shared across the databases, with I = L.

    Wraps mzLib's ``PeptideUniquenessClassifier.Classify``. A protein contains a peptide when its
    sequence contains it **anywhere, whatever the protease** - deliberately conservative, so a
    peptide called ``Unique`` cannot be explained by another entry at a site the search's cleavage
    rules happened to skip. **I and L are the same residue** throughout: ``LVLNGNPLTLFQER`` finds
    GAPDH's ``LVINGNPITIFQER``.

    - ``Unique``: every protein containing it has the **same sequence** (identical entries under two
      accessions count as one; both are listed).
    - ``SharedWithinGene``: several distinct sequences, all sharing a gene - isoforms. The peptide
      supports the gene, not an isoform.
    - ``SharedAcrossGenes``: sequences with no gene in common.
    - ``NotInDatabase``: no target protein contains it.

    Two proteins share a gene when they share any key: an Ensembl gene id, organism plus primary
    gene name (so human and bovine ALB stay apart), or the UniProt entry (so ``P12345-2`` meets
    ``P12345``). Decoys are ignored; contaminants are included, because they are real sequences in
    the search space.

    Every database is part of one search space, so there is no ``on_error="skip"``: dropping a
    database that failed to load would report the peptides it contains as unique.

    Args:
        peptides: Unmodified base sequences, upper case (``"PEPTIDEK"``), one per result row in the
            same order, duplicates included. Sent on stdin. Strip modifications first: mzLib refuses
            ``"PEPT[Phospho]IDEK"`` and ``"peptidek"`` rather than guess.
        databases: The search database(s): UniProt XML or FASTA, optionally gzipped.
        contaminants: Contaminant databases searched alongside.
        threads: Databases read at once (the classification itself is one pass). Default 1;
            ``-1`` means every core. Same answer at any value.
        timeout: Seconds to allow. ``None`` (default) waits.

    Returns:
        A :class:`PeptideClassification`, one row per peptide.

    Raises:
        UsageError: no peptides, a blank one, one that is not upper-case A-Z (from mzLib, naming
            it), or any database problem :func:`read` would refuse.
        BridgeError: mzLib could not parse a database.

    Examples:
        >>> calls = classify_peptides(
        ...     ["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
        ...     ["human_subset.xml", "human_extra.fasta"], contaminants=["contaminants.fasta"])
        >>> calls.records[3]["peptide"], calls.records[3]["sharing"], calls.records[3]["accessions"]
        ('YLYEIAR', 'SharedAcrossGenes', ['P02768', 'P02769'])
    """
    peptide_lines = _lines(peptides, "peptides")
    if not peptide_lines:
        raise _bridge.UsageError("At least one peptide is required, e.g. ['PEPTIDEK'].")

    db_args, path_lines = _database_args(databases, contaminants)
    args = ["proteins", "classify-peptides", *db_args, *_threads(threads)]
    data = _bridge.invoke(*args, stdin=_stdin(path_lines, peptide_lines), timeout=timeout)
    return PeptideClassification._from_wire(data)


def _one_path(value: _PathLike | None, what: str, *, required: bool) -> str | None:
    """One path argument as clean text, or ``None`` when it is optional and absent."""
    if value is None:
        if required:
            raise _bridge.UsageError(f"{what} is required.")
        return None
    text = _bridge.path_text(value)
    if not text:
        raise _bridge.UsageError(f"{what} must be a non-empty str or path; got {value!r}.")
    return text


def annotate_go(
    groups: _PathLike,
    database: _PathLike,
    *,
    go_obo: _PathLike,
    category_map: _PathLike | None = None,
    skip_unknown_go_ids: bool = False,
    out: _PathLike | None = None,
    categories_out: _PathLike | None = None,
    limit: int | None = None,
    offset: int = 0,
    timeout: float | None = None,
) -> GoAnnotations:
    """Annotate MetaMorpheus protein groups with Gene Ontology terms, keeping every member.

    Wraps mzLib's ``GoGroupAnnotator`` over a stored protein-group table
    (``ProteinGroupFromTsv.ToGoAnnotationGroups``). For each group it returns one row per GO term
    that **any** member holds - directly, or by propagation up ``is_a`` and ``part_of`` - and each
    row names the members that carry it. Nothing is collapsed: the union, the consensus and the
    direct-only views are filters on the rows (see :class:`GoAnnotations`).

    **Every non-decoy group gets at least one row.** A group with no term gets a single row whose
    ``annotation_status`` says why: ``no_go_terms`` (its members have none), ``no_entry`` (a member
    is not in the database - annotate against the database the search used), or ``contaminant``.

    **Pin the ontology.** Terms and their ancestors change between GO releases, so a result means
    something only relative to one go.obo. This reads the file you name and never downloads one;
    fetch a release on purpose with :func:`update_go`, keep the file, and record
    :attr:`GoAnnotations.go` (its ``sha256``) with your results.

    Args:
        groups: A MetaMorpheus protein-group table: ``AllQuantifiedProteinGroups.tsv``,
            ``AllProteinGroups.tsv`` (a search without quantification), or one file's
            ``<file>_ProteinGroups.tsv``.
        database: The UniProt XML (``.xml`` or ``.xml.gz``) whose GO terms annotate the members:
            the proteome the search used. A FASTA is refused, because it carries no GO and every
            group would read ``no_go_terms`` whatever the proteins are.
        go_obo: A go.obo file. It must exist; nothing is fetched.
        category_map: Optional. Your own term-to-category map in mzLib's format
            (``#!category_map_format 1``, ``#!map_name``, ``#!map_version``, then ``category``,
            ``subcategory`` and ``anchor_go_id`` columns). mzLib ships no vocabulary. Adds
            :attr:`GoAnnotations.categories`.
        skip_unknown_go_ids: ``False`` (default) fails when the database cites a GO id the
            ontology release lacks - usually a UniProt release newer than the go.obo - and names
            every missing id. ``True`` drops each such id and lists it in
            :attr:`GoAnnotations.unresolved_go_ids`, so one new term does not cost the run.
        out: Write the **whole** table here with mzLib's own ``GoAnnotationTsv`` writer, provenance
            header included. Must end in ``.tsv``. For a large run, pair it with ``limit=0`` so only
            the summary comes back.
        categories_out: Write the category table here with mzLib's ``GoCategoryTsv`` writer. Must
            end in ``.tsv``, and needs ``category_map``.
        limit: At most this many rows (rows, not groups) in :attr:`GoAnnotations.columns`.
            ``None`` (default) returns them all. Never shortens ``out``.
        offset: Rows to skip before the window. Default 0.
        timeout: Seconds to allow. ``None`` (default) waits: a whole proteome XML takes a while.

    Returns:
        A :class:`GoAnnotations`.

    Raises:
        UsageError: a missing groups table, database, go.obo or category map; a FASTA database;
            an ``out`` or ``categories_out`` that is not ``.tsv`` or names an input;
            ``categories_out`` without ``category_map``; a negative ``limit`` or ``offset``. All of
            these are raised before any file is read.
        BridgeError: the database cites GO ids the release lacks (``InvalidDataException``, every id
            named) unless ``skip_unknown_go_ids=True``; a go.obo, category map or table mzLib
            cannot parse; a group MetaMorpheus wrote twice with different members.

    Examples:
        >>> go = annotate_go("PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
        ...                  go_obo="go-pxd036557.obo", category_map="organelle_map.tsv")
        >>> go.group_count, go.row_count, go.go.release
        (5, 563, 'releases/2026-07-26')
        >>> go.header["status_annotated"], go.header["status_contaminant"]
        ('4', '1')
        >>> histones = [r for r in go.records if r["protein_group"] == "P0C0S5|Q71UI9"]
        >>> len(histones), sum(r["n_with"] == r["n_members"] for r in histones)
        (106, 57)
    """
    args = [
        "proteins", "annotate-go",
        "--groups", _one_path(groups, "groups", required=True) or "",
        "--database", _one_path(database, "database", required=True) or "",
        "--go-obo", _one_path(go_obo, "go_obo", required=True) or "",
    ]
    for option, value in (("--category-map", category_map), ("--out", out), ("--categories-out", categories_out)):
        text = _one_path(value, option[2:].replace("-", "_"), required=False)
        if text is not None:
            args += [option, text]
    if skip_unknown_go_ids:
        args.append("--skip-unknown-go-ids")
    for name, value in (("limit", limit), ("offset", offset)):
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise _bridge.UsageError(f"{name} must be an int, zero or greater; got {value!r}.")
        args += [f"--{name}", str(value)]

    _bridge.require_verb("proteins annotate-go", since=_GO_SINCE)
    return GoAnnotations._from_wire(_bridge.invoke(*args, timeout=timeout))


def update_go(go_obo: _PathLike, *, timeout: float | None = None) -> GoUpdate:
    """Download the current Gene Ontology release to ``go_obo``, on purpose.

    Wraps mzLib's ``Loaders.UpdateGeneOntology``. The whole go.obo (about 37 MB) is streamed from
    GO's PURL, which always serves the **current** release. When a file is already at the path, it
    is kept beside the new one as ``go.obo.<yyyyMMdd-HHmmss-fff>`` if the download differs, and left
    alone if it is the same, so earlier runs stay reproducible. A failed download leaves any
    existing file untouched.

    This is the only function in pyMzLib that fetches go.obo. :func:`annotate_go` never does, so
    the release a result was computed against is always a file you chose to keep.

    Args:
        go_obo: Where to write, e.g. ``"go.obo"``. Its folder must exist.
        timeout: Seconds to allow. ``None`` (default) waits; mzLib itself gives up after two
            minutes without data.

    Returns:
        A :class:`GoUpdate`.

    Raises:
        UsageError: the folder does not exist.
        BridgeError: GO's server is unreachable or answered with an error. ``type`` is
            ``"ServiceUnavailable"`` for an outage, where retrying later makes sense.

    Examples:
        >>> update = update_go("go.obo")
        >>> update.go.release, update.go.term_count, update.changed
        ('releases/2026-07-26', 48340, True)
    """
    path = _one_path(go_obo, "go_obo", required=True) or ""
    _bridge.require_verb("proteins update-go", since=_GO_SINCE)
    return GoUpdate._from_wire(_bridge.invoke("proteins", "update-go", "--go-obo", path, timeout=timeout))


def _pair_lines(pairs: Iterable[Sequence[str]]) -> list[str]:
    """The pairs as stdin lines, ``source<TAB>target``, unchanged, refusing what cannot travel."""
    if isinstance(pairs, (str, bytes)):
        raise _bridge.UsageError(
            "pairs must be a list of (source, target) tuples, not one string. For one pair, call "
            "align(source, target)."
        )
    lines = []
    for index, pair in enumerate(pairs):
        if isinstance(pair, (str, bytes)) or not isinstance(pair, Sequence) or len(pair) != 2:
            raise _bridge.UsageError(f"Pair {index} must be a (source, target) tuple; got {pair!r}.")
        source, target = pair
        for name, value in (("source", source), ("target", target)):
            if not isinstance(value, str):
                raise _bridge.UsageError(f"Pair {index}'s {name} must be a str; got {value!r}.")
            # One pair is one stdin line, its two sequences split on the tab, so neither may hold
            # a tab or a line break: either would shift or split every later pair.
            if "\t" in value or "\n" in value or "\r" in value:
                raise _bridge.UsageError(f"Pair {index}'s {name} contains a tab or line break: {value!r}.")
        # The bridge skips a blank line, which would put every later result on the wrong pair.
        if not (source + target).strip():
            raise _bridge.UsageError(f"Pair {index} has two empty sequences; there is nothing to align.")
        lines.append(f"{source}\t{target}")
    return lines


def _aligner_args(
    gap_open: int | None, gap_extend: int | None, end_gaps: str | None, max_cells: int | None
) -> list[str]:
    """The aligner options that were given, and only those: mzLib's defaults apply to the rest."""
    args: list[str] = []
    for name, value in (("gap-open", gap_open), ("gap-extend", gap_extend), ("max-cells", max_cells)):
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, int):
            raise _bridge.UsageError(f"{name.replace('-', '_')} must be an int or None; got {value!r}.")
        args += [f"--{name}", str(value)]
    if end_gaps is not None:
        if end_gaps not in END_GAPS:
            raise _bridge.UsageError(f"end_gaps must be one of {END_GAPS} or None; got {end_gaps!r}.")
        args += ["--end-gaps", end_gaps]
    return args


def align_many(
    pairs: Iterable[Sequence[str]],
    *,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    end_gaps: str | None = None,
    max_cells: int | None = None,
    residue_map: bool = True,
    timeout: float | None = None,
) -> Alignments:
    """Align many pairs of protein sequences in one call, and map every residue across each pair.

    Wraps mzLib's ``PairwiseAligner.Align`` (mzLib #1428): a **global** alignment with affine gaps
    (Gotoh's form of Needleman-Wunsch) scored with **BLOSUM62**, and the residue map its
    ``TryMapSourceToTarget`` and ``TryMapTargetToSource`` give. Use it to carry a position - a
    modification site, a variant, a cleavage site - from one sequence to another: an isoform, an
    ortholog, an older release of the same entry.

    **Every option you leave out is mzLib's to set.** pyMzLib sends nothing for it, and the result
    says what the aligner used (:attr:`Alignments.gap_open` and the rest). At mzLib 1.0.595 the
    defaults are gap open 11, gap extend 1, end gaps free and a 200,000,000-cell limit: BLASTP's
    gap costs for BLOSUM62, and EMBOSS needle's free end gaps.

    Args:
        pairs: ``(source, target)`` tuples of one-letter sequences, one result row each, in order.
            Sent to mzLib exactly as given (case is kept, and compared case-insensitively). Either
            sequence may be empty, not both. mzLib refuses ``-`` and ``.``; it scores any other
            letter BLOSUM62 does not name (``U``, ``O``, a digit, a space) as ``X``, so strip
            modifications first.
        gap_open: Cost of opening a gap, in BLOSUM62's half-bits, charged once per gap: a gap of
            ``k`` residues costs ``gap_open + k * gap_extend``. ``None`` (default): mzLib's, 11.
            mzLib refuses a negative value.
        gap_extend: Cost of each residue in a gap, the first included, in half-bits. ``None``
            (default): mzLib's, 1. mzLib refuses less than 1.
        end_gaps: ``"free"`` - a gap at either end of either sequence costs nothing, so a signal
            peptide or a longer tail on one side is not penalised - or ``"charged"``. ``None``
            (default): mzLib's, ``"free"``. See :data:`END_GAPS`.
        max_cells: The largest table mzLib will build, ``(len(source) + 1) * (len(target) + 1)``
            cells of one byte each. A larger pair is refused rather than allowed to exhaust memory.
            ``None`` (default): mzLib's, 200,000,000 cells (about 200 MB). Two titin isoforms need
            about 1.2 billion; mzLib accepts up to 2,147,483,591.
        residue_map: ``True`` (default) returns :attr:`Alignments.residue_map`, one row per residue
            of both sequences. ``False`` sends ``--no-residue-map`` and returns only the pairs
            table, for a large batch where you want the scores and counts.
        timeout: Seconds to allow. ``None`` (default) waits.

    Returns:
        An :class:`Alignments`. Its :attr:`~Alignments.alignments` gives one :class:`Alignment`
        per pair, with :meth:`~Alignment.map_source_to_target`.

    Raises:
        UsageError: no pairs; a pair that is not two strings; a sequence with a tab or line
            break; two empty sequences; an option of the wrong type or ``end_gaps`` not in
            :data:`END_GAPS`. From mzLib, through the bridge: a gap cost or ``max_cells`` out of
            range, a sequence holding ``-`` or ``.``, or a pair whose table exceeds ``max_cells``
            (the message names the pair; the whole call fails).

    Examples:
        A deletion and an extra N-terminus, mzLib's own test cases:

        >>> alns = align_many([("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ"),
        ...                    ("MLLLPEPTIDE", "PEPTIDE")])
        >>> alns.aligner_id
        'global-affine;BLOSUM62;open=11;extend=1;end-gaps=free'
        >>> alns.columns["score"], alns.columns["source_residues_facing_gap"]
        ([80, 39], [3, 4])
        >>> alns.columns["aligned_target"]
        ['MKTAYIAK---ISFVKSHFSRQ', '----PEPTIDE']
    """
    lines = _pair_lines(pairs)
    if not lines:
        raise _bridge.UsageError("At least one pair is required, e.g. [('MKTAYIAK', 'MKTSYIAK')].")
    args = ["proteins", "align", *_aligner_args(gap_open, gap_extend, end_gaps, max_cells)]
    if not residue_map:
        args.append("--no-residue-map")

    _bridge.require_verb("proteins align", since=_ALIGN_SINCE)
    data = _bridge.invoke(*args, stdin="\n".join(lines) + "\n", timeout=timeout)
    return Alignments._from_wire(data)


def align(
    source: str,
    target: str,
    *,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    end_gaps: str | None = None,
    max_cells: int | None = None,
    timeout: float | None = None,
) -> Alignment:
    """Align two protein sequences with mzLib, and map residues from one onto the other.

    The one-pair form of :func:`align_many`: the same bridge call and the same mzLib defaults,
    returning one :class:`Alignment` with its residue map.

    Args:
        source: The sequence whose positions you have, one-letter codes, e.g. the canonical isoform.
        target: The sequence you want them in.
        gap_open: Gap-open cost in BLOSUM62's half-bits; ``None`` (default) is mzLib's, 11.
        gap_extend: Cost per gap residue in half-bits; ``None`` (default) is mzLib's, 1.
        end_gaps: ``"free"`` or ``"charged"``; ``None`` (default) is mzLib's, ``"free"``.
        max_cells: Table limit in cells; ``None`` (default) is mzLib's, 200,000,000.
        timeout: Seconds to allow. ``None`` (default) waits.

    Returns:
        An :class:`Alignment`.

    Raises:
        UsageError: as :func:`align_many`.

    Examples:
        mzLib's own deletion case: ``QRQ`` (source 9 to 11) is missing from the target.

        >>> aln = align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")
        >>> print(aln.aligned_source, aln.aligned_target, sep="\\n")
        MKTAYIAKQRQISFVKSHFSRQ
        MKTAYIAK---ISFVKSHFSRQ
        >>> aln.score, aln.identical, aln.source_residues_facing_gap
        (80, 19, 3)
        >>> [aln.map_source_to_target(p) for p in (8, 9, 10, 11, 12)]
        [8, None, None, None, 9]

        Charging the end gaps changes the score of an overhang, and the aligner id says so:

        >>> charged = align("MLLLPEPTIDE", "PEPTIDE", end_gaps="charged")
        >>> charged.score, charged.aligner_id
        (24, 'global-affine;BLOSUM62;open=11;extend=1;end-gaps=charged')
    """
    for name, value in (("source", source), ("target", target)):
        if not isinstance(value, str):
            raise _bridge.UsageError(f"{name} must be a str; got {value!r}.")
    result = align_many(
        [(source, target)], gap_open=gap_open, gap_extend=gap_extend, end_gaps=end_gaps,
        max_cells=max_cells, timeout=timeout,
    )
    return result.alignments[0]


def _accession_lines(accessions: Iterable[str]) -> list[str]:
    """The accessions as stdin lines, unchanged: mzLib parses and never repairs, so neither do we."""
    if isinstance(accessions, str):
        raise _bridge.UsageError(
            "accessions must be a list of strings, not one string: that would send one accession "
            f"per character. Wrap it: [{accessions!r}]."
        )
    lines = []
    for value in accessions:
        if not isinstance(value, str):
            raise _bridge.UsageError(f"Every accession must be a str; got {value!r}.")
        # The bridge skips blank lines and splits on line breaks, so either would put every later
        # row on the wrong accession. Refused here rather than silently re-aligned. A space around
        # an accession is NOT stripped: mzLib reads it as unrecognized, and that is the answer.
        if not value.strip():
            raise _bridge.UsageError("accessions contains a blank entry; every entry is one result row.")
        if "\n" in value or "\r" in value:
            raise _bridge.UsageError(f"An accession contains a line break: {value!r}.")
        lines.append(value)
    return lines


def parse_accessions(accessions: Iterable[str], *, timeout: float | None = None) -> ParsedAccessions:
    """Read accessions into entry, isoform or version, and applied variants, with mzLib.

    Wraps mzLib's ``VariantApplication.ParseAccession`` (mzLib #1382), and through it
    ``ProteinAccession.Parse``. It is the inverse of the names mzLib gives variant proteoforms:
    ``P04637_S70N`` is entry ``P04637`` with ``S70N`` applied, and ``NP_000537.3_R72P`` is RefSeq
    ``NP_000537``, version 3, with ``R72P``. Use it to link a search result's accession to anything
    you keep per entry - genes, orthologs, GO terms - without losing the isoform or the variants.

    **It parses and never repairs.** An accession outside UniProt's and RefSeq's protein grammars
    comes back ``"unrecognized"`` and verbatim - lower case, a stray space, a FASTA header, a
    ``DECOY_`` prefix, a ``|``-joined protein group, and ``P12345_2``, mzLib's loader counter for a
    *different* entry whose accession collided. pyMzLib sends each accession exactly as you give it.

    Args:
        accessions: One accession per result row, in order, duplicates included. Sent on stdin,
            unchanged: not stripped, not split. A blank entry or one with a line break is refused,
            because it would shift every later row.
        timeout: Seconds to allow. ``None`` (default) waits.

    Returns:
        A :class:`ParsedAccessions`, one row per accession.

    Raises:
        UsageError: no accessions; one that is not a str, is blank, or holds a line break; a single
            string instead of a list.

    Examples:
        >>> parsed = parse_accessions(["P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P",
        ...     "P04637_S70N_A80T", "P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"])
        >>> parsed.columns["entry_accession"][:5]
        ['P04637', 'P04637', 'Q13409', 'NP_000537', 'P04637']
        >>> parsed.columns["isoform"][:3], parsed.columns["version"][3]
        ([None, 2, 2], 3)
        >>> parsed.columns["applied_variants"][3:5]
        ['R72P', 'S70N_A80T']
        >>> parsed.namespace_counts
        {'unrecognized': 4, 'uniprot': 4, 'refseq': 1}
    """
    lines = _accession_lines(accessions)
    if not lines:
        raise _bridge.UsageError("At least one accession is required, e.g. ['P04637-2'].")
    _bridge.require_verb("proteins parse-accession", since=_ALIGN_SINCE)
    data = _bridge.invoke("proteins", "parse-accession", stdin="\n".join(lines) + "\n", timeout=timeout)
    return ParsedAccessions._from_wire(data)


def _species_lines(
    gtfs: Mapping[str, _PathLike] | None, gene_sets: Mapping[str, _PathLike] | None
) -> list[str]:
    """The stdin species lines, ``species<TAB>gtf|gene-set<TAB>path``, in the order given."""
    lines: list[str] = []
    seen: set[str] = set()
    for kind, mapping, what in (("gtf", gtfs, "gtfs"), ("gene-set", gene_sets, "gene_sets")):
        if mapping is None:
            continue
        if not isinstance(mapping, Mapping):
            raise _bridge.UsageError(
                f"{what} must be a dict from Ensembl production name to path, e.g. "
                f"{{'homo_sapiens': 'Homo_sapiens.GRCh38.116.gtf.gz'}}; got {mapping!r}."
            )
        for species, path in mapping.items():
            if not isinstance(species, str) or not species.strip():
                raise _bridge.UsageError(
                    f"Every species in {what} must be a non-empty str; got {species!r}."
                )
            name = species.strip()
            if any(c in name for c in "\t\r\n,"):
                raise _bridge.UsageError(
                    f"A species name cannot hold a tab, line break or comma: {species!r}."
                )
            if name in seen:
                raise _bridge.UsageError(
                    f"'{name}' is in both gtfs and gene_sets; give each species one gene set."
                )
            seen.add(name)
            lines.append(f"{name}\t{kind}\t{_paths(path, f'{what} {name}')[0]}")
    if not lines:
        raise _bridge.UsageError(
            "At least one species is required: gtfs={'homo_sapiens': 'Homo_sapiens.GRCh38.116.gtf.gz', "
            "...} (or gene_sets= for compact gene tables)."
        )
    return lines


def build_orthology(
    release: str | int,
    *,
    gene_trees: _PathLike,
    dumps: Iterable[_PathLike],
    gtfs: Mapping[str, _PathLike] | None = None,
    gene_sets: Mapping[str, _PathLike] | None = None,
    restrict_gene_trees: bool = False,
    out: _PathLike | None = None,
    pair_status: Sequence[str] | None = None,
    timeout: float | None = None,
) -> OrthologySnapshot:
    """Build an Ensembl Compara orthology snapshot for your species, and write it as Parquet.

    Wraps mzLib's ``OrthologySnapshot.Build`` and ``OrthologySnapshotWriter.WriteAsync`` (mzLib
    #1381). Which genes of one species are orthologous to which genes of another, kept exactly as
    Compara states it: one-to-many and many-to-many orthologies are several rows, paralogs are kept,
    and nothing is picked, merged or collapsed. Local files only; nothing is downloaded.

    **Every species' homology dump is required, and mzLib refuses a call without one.** Compara puts
    each of a pair's relationships in one of the pair's two dumps, arbitrarily, and in practice they
    partition: in release 116 every human-mouse orthology is in the *mouse* dump and none is in the
    human one. Reading the dump named after your species gives a complete-looking, badly incomplete
    answer.

    **Orthology is pairwise and not transitive.** A~B and B~C do not make A~C, so never chain pairs
    into groups. The written ``views.sql`` has a ``species_set3`` macro that forms a tuple only when
    every two members are orthologs.

    The relationships are in the Parquet files ``out=`` writes, which mzLib cannot read back and
    pyMzLib does not read (it has no third-party dependency): query them with DuckDB, pyarrow or
    R's arrow, as the proteins guide shows. Without ``out=``, the call
    checks that the inputs fit together and, with ``pair_status=``, says gene by gene how far
    Compara's inference got in one other species.

    Args:
        release: The Ensembl release of the Compara inputs, e.g. ``"116"`` (an int is sent as its
            digits). Required: mzLib holds the gene trees and every dump to it when their file names
            carry a release. A GTF's own number is not compared with it (Ensembl publishes some
            species with their own annotation release).
        gtfs: Species -> its Ensembl primary-assembly GTF, plain or ``.gz``, e.g.
            ``{"homo_sapiens": "Homo_sapiens.GRCh38.116.gtf.gz"}``. The keys are Ensembl production
            names, as the dumps' ``species`` column spells them, and they are the snapshot's species.
            Sent on stdin, one line per species.
        gene_sets: Species -> a compact gene table mzLib's ``EnsemblGeneSetWriter`` made from a GTF,
            instead of the GTF itself (it keeps the GTF's provenance, so the snapshot is the same).
            A species may be in ``gtfs`` or ``gene_sets``, not both.
        dumps: One Compara homology dump per species
            (``homologies/<species>/Compara.<release>.protein_default.homologies.tsv.gz``), in any
            order: each dump's species is read from its rows. Sent on stdin after the species lines.
        gene_trees: Compara's gene-tree file,
            ``<collection>.GeneTree_content.default.e<release>.txt.gz``. It tells a gene whose tree
            holds no gene of the other species from one whose tree does but got no ortholog call.
        restrict_gene_trees: Keep only the gene-tree members in your gene sets. ``False``
            (default) is mzLib's default and keeps every gene in the file - 4.2 million for the
            vertebrates of release 116. The snapshot is the same either way; only memory differs.
        out: A directory to write the snapshot to: Parquet tables (``genes/``, ``members/``,
            ``pairs/``), ``views.sql`` and ``manifest.json``. It must not exist yet or be empty:
            mzLib writes a snapshot once and never updates one in place. ``None`` (default): nothing
            is written.
        pair_status: ``(from_species, to_species)``: return every gene of ``from_species`` with its
            status in ``to_species`` (:class:`PairStatus`). Two different species of the snapshot,
            in that order: a status is not symmetric. ``None`` (default): not computed.
        timeout: Seconds to allow. ``None`` (default) waits: for human, mouse and rat at release
            116, mzLib reports about 20 s to build and 40 s more to write.

    Returns:
        An :class:`OrthologySnapshot`.

    Raises:
        UsageError: a missing input file; no species, no dumps, or a species in both ``gtfs`` and
            ``gene_sets``; ``out`` naming a file or a directory that is not empty; ``pair_status``
            not two different species of the snapshot (all checked before any input is read); and
            mzLib's refusal of inputs that do not fit together - a species without a dump, two
            dumps of one species, a dump of a species not given, or a file of another release.
        BridgeError: a file mzLib's reader refuses (``InvalidDataException``, naming the file and
            line), or data that breaks the snapshot's rules, such as a homology naming a gene its
            species' gene set lacks: the GTF and the dumps are from different releases.

    Examples:
        mzLib's own test world: three genes of human, mouse and rat, with every human-mouse row
        in the mouse dump, as in release 116.

        >>> gtfs = {"homo_sapiens": "Homo_sapiens.GRCh38.116.gtf",
        ...         "mus_musculus": "Mus_musculus.GRCm39.116.gtf",
        ...         "rattus_norvegicus": "Rattus_norvegicus.GRCr8.116.gtf"}
        >>> dumps = [f"homologies/{s}/Compara.116.protein_default.homologies.tsv.gz" for s in gtfs]
        >>> snap = build_orthology("116", gtfs=gtfs, dumps=dumps,
        ...     gene_trees="vertebrates.GeneTree_content.default.e116.txt.gz",
        ...     pair_status=("homo_sapiens", "mus_musculus"))
        >>> snap.dumps["mus_musculus"].rows_kept, snap.pair("mus_musculus", "homo_sapiens").row_count
        (5, 3)
        >>> snap.pair_status.columns["status"]
        ['has_ortholog', 'has_ortholog', 'tree_lacks_target_species', 'not_in_any_tree']
        >>> snap.pair_status.genes("not_in_any_tree")
        ['ENSG00000000004']

        With ``out=``, the snapshot is written, and its id is the same:

        >>> written = build_orthology("116", gtfs=gtfs, dumps=dumps,
        ...     gene_trees="vertebrates.GeneTree_content.default.e116.txt.gz",
        ...     restrict_gene_trees=True, out="compara-116")
        >>> written.snapshot_id == snap.snapshot_id
        True
        >>> [f.rows for f in written.written.files if f.path.startswith("pairs/")]
        [0, 3, 1, 1, 2, 0]
    """
    if isinstance(release, bool) or not isinstance(release, (str, int)):
        raise _bridge.UsageError(f"release must be an Ensembl release, e.g. '116'; got {release!r}.")
    release_text = str(release).strip()
    if not release_text or any(c.isspace() for c in release_text):
        raise _bridge.UsageError(f"release must be an Ensembl release, e.g. '116'; got {release!r}.")

    species_lines = _species_lines(gtfs, gene_sets)
    dump_lines = _paths(dumps, "homology dump")
    if not dump_lines:
        raise _bridge.UsageError(
            "dumps is empty. Give every species' Compara homology dump: each pair's relationships are "
            "split between the pair's two dumps, so mzLib requires them all."
        )

    args = [
        "genes", "orthology", "--release", release_text,
        "--gene-trees", _one_path(gene_trees, "gene_trees", required=True) or "",
    ]
    if not isinstance(restrict_gene_trees, bool):
        raise _bridge.UsageError(
            f"restrict_gene_trees must be True or False; got {restrict_gene_trees!r}."
        )
    if restrict_gene_trees:
        args.append("--restrict-gene-trees")
    out_text = _one_path(out, "out", required=False)
    if out_text is not None:
        args += ["--out", out_text]
    if pair_status is not None:
        if (
            isinstance(pair_status, (str, bytes))
            or not isinstance(pair_status, Sequence)
            or len(pair_status) != 2
        ):
            raise _bridge.UsageError(
                f"pair_status must be a (from_species, to_species) tuple, e.g. ('homo_sapiens', "
                f"'mus_musculus'); got {pair_status!r}."
            )
        for name in pair_status:
            if not isinstance(name, str) or not name.strip() or any(c in name for c in "\t\r\n,"):
                raise _bridge.UsageError(
                    f"Each species in pair_status must be a production name; got {name!r}."
                )
        args += ["--pair-status", ",".join(name.strip() for name in pair_status)]

    _bridge.require_verb("genes orthology", since=_ORTHOLOGY_SINCE)
    stdin = "\n".join([*species_lines, _SECTION, *dump_lines]) + "\n"
    data = _bridge.invoke(*args, stdin=stdin, timeout=timeout)
    return OrthologySnapshot._from_wire(data)
