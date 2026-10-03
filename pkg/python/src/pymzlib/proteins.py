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

These examples are the calls that recorded the payloads in pyMzLib's test fixtures, over the small
databases in ``pkg/python/tests/fixtures/proteins/``. The GO example's table is a real MetaMorpheus
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

Wire verbs: ``proteins read``, ``genes resolve``, ``proteins classify-peptides``,
``proteins annotate-go`` and ``proteins update-go``; their language-neutral specs live in the bridge
repository under ``design/verbs/``.
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
    "NOT_IN_DATABASE",
    "OUTCOMES",
    "SHARED_ACROSS_GENES",
    "SHARED_WITHIN_GENE",
    "TABLES",
    "UNIQUE",
    "AnnotationDatabase",
    "DatabaseFile",
    "FileError",
    "GeneResolutions",
    "GeneSet",
    "GoAnnotations",
    "GoCategories",
    "GoRelease",
    "GoUpdate",
    "PeptideClassification",
    "ProteinDatabase",
    "Table",
    "WrittenTable",
    "XrefTable",
    "annotate_go",
    "classify_peptides",
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
