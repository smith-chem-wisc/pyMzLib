"""Peptidoform-level questions: digest an annotated protein and fragment its peptides.

The question this answers is the one a mass spectrometrist actually asks — *what fragments would
I see for this protein's peptides?* — in one call:

    >>> import pymzlib
    >>> digest = pymzlib.peptidoform.fragments("P02768")
    >>> digest.accession, digest.sequence_length
    ('P02768', 609)

The defaults are opinions, not placeholders. Tryptic with the proline rule, two missed cleavages,
ETD, both termini, UniProt's annotated modifications applied. They are the choices this lab makes
when it does not have a reason to choose otherwise, so the common question needs no parameters —
and every one of them is reachable, because the point is to open the doors, not to hide them.

:func:`convert` rewrites full sequences from one notation to another with mzLib's
``SequenceConversionService``: a MetaMorpheus full sequence ``[UniProt:N-acetylserine on S]SEQK``
becomes ``[UNIMOD:1]SEQK`` in Unimod notation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping

from . import _bridge
from .readers import _Table

#: UniProtKB's own accession grammar (https://www.uniprot.org/help/accession_numbers). Checking
#: it here means a typo costs nothing instead of a network round trip and a puzzling HTTP 400.
_ACCESSION = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

__all__ = [
    "Fragment",
    "Peptide",
    "Digest",
    "ModificationCensus",
    "fragments",
    "CONVERSION_MODES",
    "ConvertedSequence",
    "SequenceConversions",
    "convert",
]


@dataclass(frozen=True)
class Fragment:
    """One backbone fragment ion.

    Attributes:
        product_type: The ion series, e.g. ``"c"`` or ``"zDot"`` for ETD, ``"b"``/``"y"`` for CID.
        fragment_number: Position in the series — ``c3`` is the third from the N-terminus.
        neutral_mass: Monoisotopic neutral mass in daltons. **Not** an m/z: no proton has been
            added and no charge assumed.

            Fragments deliberately expose no ``mz()``. Converting one correctly requires the fixed
            charge *within this fragment's span* — a c or z ion carries only the permanently
            charged modifications on the residues it contains, not the whole peptide's
            :attr:`Peptide.fixed_charges`. For an unmodified or neutrally-modified peptide,
            ``(neutral_mass + z * PROTON_MASS) / z`` is correct; for a fragment bearing a
            quaternary-ammonium modification (e.g. trimethyllysine) it is not, and per-fragment
            charge accounting is not yet provided.
        neutral_loss: Neutral loss in daltons, ``0.0`` when there is none.
        residue_position: One-based residue position in the peptide.
    """

    product_type: str
    fragment_number: int
    neutral_mass: float
    neutral_loss: float
    residue_position: int


@dataclass(frozen=True)
class Peptide:
    """One digested peptide, with its modifications and fragment ions.

    Attributes:
        base_sequence: The bare amino-acid sequence.
        full_sequence: The sequence with modifications written inline, as mzLib renders them.
        monoisotopic_mass: The neutral monoisotopic mass, modifications included.
        one_based_start / one_based_end: Position within the parent protein.
        missed_cleavages: How many cleavage sites the peptide spans.
        fixed_charges: Charges the peptide carries before any protonation, from modifications
            that leave a permanently charged residue. :meth:`mz` accounts for these.
        modifications: Each applied modification. ``one_based_residue`` indexes the peptide's own
            residues and is ``None`` for a terminal modification, which carries ``terminus``
            instead. (mzLib's internal dictionary reserves slot 1 for the N-terminus, so its keys
            are one past the residue they modify; that is corrected here rather than passed on.)
        fragments: The fragment ions for the requested dissociation type.
    """

    base_sequence: str
    full_sequence: str
    monoisotopic_mass: float
    length: int
    one_based_start: int
    one_based_end: int
    missed_cleavages: int
    fixed_charges: int = 0
    modifications: list[dict[str, Any]] = field(default_factory=list)
    fragments: list[Fragment] = field(default_factory=list)

    @property
    def is_modified(self) -> bool:
        """Whether this peptide carries at least one modification."""
        return bool(self.modifications)

    def mz(self, charge: int) -> float:
        """Return the m/z of the intact peptide at a given total charge.

        Two conventions are handled explicitly here, because getting either wrong is invisible in
        the answer.

        **The proton mass (1.007276), not the hydrogen atom (1.007825).** The difference is
        0.55 mDa — 1.1 ppm at m/z 500, which on an Orbitrap is a match versus a miss. Libraries
        differ on this and rarely say which they used.

        **Fixed charges are not double-counted.** Some modifications leave the residue permanently
        charged: trimethylation of a lysine ε-amine gives a quaternary ammonium, and UniProt
        records the delta as 43.054227 — C₃H₇ *minus an electron* — rather than the neutral
        43.054775. So :attr:`monoisotopic_mass` already carries that charge, and only
        ``charge - fixed_charges`` protons are added. Adding a full complement would put a 2+
        trimethylated peptide half a Thomson high, on the most important histone modification
        there is.

        A peptide with a fixed charge is therefore observable at that charge with no protonation
        at all, which is why ``charge`` may not be below :attr:`fixed_charges`.

        Args:
            charge: The total charge state, at least :attr:`fixed_charges` and at least 1.
        """
        if not isinstance(charge, int) or isinstance(charge, bool) or charge < 1:
            raise _bridge.UsageError(f"charge must be a positive whole number; got {charge!r}.")
        if charge < self.fixed_charges:
            raise _bridge.UsageError(
                f"This peptide already carries {self.fixed_charges} fixed charge(s) from its "
                f"modifications, so it cannot be observed at charge {charge}."
            )
        return (self.monoisotopic_mass + (charge - self.fixed_charges) * PROTON_MASS) / charge

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "Peptide":
        return cls(
            base_sequence=payload.get("base_sequence", ""),
            full_sequence=payload.get("full_sequence", ""),
            monoisotopic_mass=float(payload.get("monoisotopic_mass", 0.0)),
            length=int(payload.get("length", 0)),
            one_based_start=int(payload.get("one_based_start", 0)),
            one_based_end=int(payload.get("one_based_end", 0)),
            missed_cleavages=int(payload.get("missed_cleavages", 0)),
            fixed_charges=int(payload.get("fixed_charges", 0)),
            modifications=list(payload.get("modifications") or []),
            fragments=[
                Fragment(
                    product_type=f.get("product_type", ""),
                    fragment_number=int(f.get("fragment_number", 0)),
                    neutral_mass=float(f.get("neutral_mass", 0.0)),
                    neutral_loss=float(f.get("neutral_loss", 0.0)),
                    residue_position=int(f.get("residue_position", 0)),
                )
                for f in (payload.get("fragments") or [])
            ],
        )


#: The proton mass, in daltons. Stated here rather than buried, because which of the two nearby
#: constants a library used is invisible in its output and changes an answer by ~1 ppm.
PROTON_MASS = 1.00727646677


@dataclass(frozen=True)
class ModificationCensus:
    """What UniProt annotates, and what could actually be used.

    mzLib loads only ``modified residue`` and ``lipid moiety-binding region`` annotations; every
    other feature type is dropped **on feature type alone**, before any mass lookup. So the census
    sees the world at feature-*type* granularity — one entry per type in :attr:`by_type`, never per
    modification *name*. On serum albumin the 24 excluded features all sit under the single type
    ``glycosylation site``; at UniProt's finer name level 22 of those 24 are specifically
    ``N-linked (Glc) (glycation) lysine``, which *does* have a defined mass — but the census never
    surfaces that name, so "22" is a fact you confirm by reading the UniProt entry, not a number
    this class reports. Read the exclusion as "wrong feature type", **not** "no defined mass".

    The exclusion is still correct: glycation and glycosylation are labile, heterogeneous adducts,
    so assigning one an exact mass and a clean fragment ladder would describe a species you cannot
    observe. What this class exists for is that you should not have to *guess* it happened: for
    serum albumin, 14 sites are applied out of 38 annotated, and without this the 14 arrives with
    no indication that a rule was ever applied. See smith-chem-wisc/mzLib#1112.

    Attributes:
        sites: Distinct residue positions carrying at least one modification. A histone lists
            several alternatives at one residue — K9me1, K9me2, K9me3, K9ac are four
            modifications at one site — so this is always the smaller number and is not a
            modification count.
        applied: Modifications actually placed on the protein.
        annotated: Modification-like features UniProt lists.
        by_type: One entry per feature type, with ``count`` and whether it was ``loaded``.
        unresolved: Modification names UniProt annotated that could not be resolved to a mass —
            usually because the name is absent from UniProt's own ptmlist. These vanish silently
            otherwise: on histone H3.1, seven N6-lactoyllysine sites were dropped while the type
            summary still reported "modified residue … loaded".
    """

    sites: int
    applied: int
    annotated: int
    unresolved: list[str] = field(default_factory=list)
    by_type: list[dict[str, Any]] = field(default_factory=list)

    @property
    def excluded(self) -> int:
        """Annotated features mzLib did not apply, dropped on feature type (not for want of mass)."""
        return max(0, self.annotated - self.applied)

    def explain(self) -> str:
        """A one-paragraph, human-readable account of what was used and what was not.

        It names the excluded feature *types* and their counts — the only granularity the census
        has. It never reports a modification-*name*-level breakdown (e.g. "22 of 24 are glycation"),
        because :attr:`by_type` does not carry names; such a figure comes from reading the UniProt
        entry, not from this census.
        """
        if not self.excluded:
            return (
                f"All {self.annotated} annotated modifications were applied, across "
                f"{self.sites} residue positions."
            )
        sentences = [
            f"{self.applied} of {self.annotated} annotated modifications were applied, across "
            f"{self.sites} residue positions."
        ]
        excluded_types = [t for t in self.by_type if not t.get("loaded")]
        if excluded_types:
            named = ", ".join(f"{t['count']} × {t['type']}" for t in excluded_types)
            sentences.append(
                f"Excluded by type: {named} — mzLib loads only 'modified residue' "
                "and 'lipid moiety-binding region' annotations, so these were dropped on "
                "feature type alone. The exclusion is usually right: a glycation or "
                "glycosylation annotation describes a labile, heterogeneous adduct, so "
                "assigning it one exact mass and a clean fragment ladder would invent a "
                "species you cannot observe. But the reason is not reported, and the "
                "qualifier is not read: some annotations are marked 'in vitro' and some "
                "exist only in disease variants, which are different "
                "grounds for exclusion needing different judgements from you. Read the "
                "annotations on the UniProt entry before concluding anything about a "
                "specific site; this census can only tell you the count "
                "(smith-chem-wisc/mzLib#1112)."
            )
        if self.unresolved:
            sentences.append(
                f"Could not be resolved to a mass: {', '.join(self.unresolved)} — annotated by "
                "UniProt but absent from its own modification list, so they were dropped."
            )
        return " ".join(sentences)


@dataclass(frozen=True)
class Digest:
    """The result of digesting a protein and fragmenting its peptides.

    Attributes:
        accession: The accession UniProt returned.
        name: UniProt entry name, e.g. ``"ALBU_HUMAN"``.
        full_name: Recommended protein name, e.g. ``"Albumin"``.
        organism: Scientific name of the source organism.
        sequence_length: Length of the full precursor sequence, in residues.
        protease: The protease used, as mzLib names it.
        dissociation: The dissociation type used, which decides the fragment series.
        terminus: Which fragment termini were generated: ``"Both"``, ``"N"`` or ``"C"``.
        modifications_applied: ``False`` when ``modifications=False`` was passed.
        max_modifications: The ``max_modifications`` cap, echoed, in modifications per peptidoform.
        max_isoforms: The ``max_isoforms`` cap, echoed, in peptidoforms per peptide position.
        peptides_at_cap: Peptide positions whose peptidoforms reached ``max_isoforms``. Non-zero
            means the list is **truncated**; :attr:`truncated` is the same fact as a yes/no.
        modification_census: What UniProt annotates against what mzLib applied: distinct residue
            positions carrying a modification (``sites``), modifications loaded (``applied``),
            modification-like features annotated (``annotated``), names that resolved to no mass
            (``unresolved``), and a per-feature-type summary (``by_type``). See
            :class:`ModificationCensus`.
        peptides: One :class:`Peptide` per distinct peptidoform, in digestion order.
    """

    accession: str
    name: str
    full_name: str
    organism: str
    sequence_length: int
    protease: str
    dissociation: str
    terminus: str
    modifications_applied: bool
    max_modifications: int
    max_isoforms: int
    peptides_at_cap: int
    modification_census: ModificationCensus
    peptides: list[Peptide] = field(default_factory=list)

    @property
    def truncated(self) -> bool:
        """Whether any peptide hit the isoform cap, meaning the result is incomplete.

        A short answer and a truncated answer look identical from the outside. Check this before
        treating a Peptidoform list as exhaustive.
        """
        return self.peptides_at_cap > 0

    @property
    def modified_peptides(self) -> list[Peptide]:
        """Only the peptides carrying at least one modification."""
        return [p for p in self.peptides if p.is_modified]

    @property
    def fragment_count(self) -> int:
        """Total fragment ions across every peptide."""
        return sum(len(p.fragments) for p in self.peptides)


def fragments(
    accession: str,
    protease: str = "trypsin|P",
    dissociation: str = "ETD",
    modifications: bool = True,
    missed_cleavages: int = 2,
    min_length: int = 7,
    max_length: int | None = None,
    max_modifications: int = 2,
    max_isoforms: int = 1024,
    terminus: str = "Both",
    timeout: float | None = 300,
) -> Digest:
    """Fetch a UniProt entry, digest it, and fragment every peptide.

    Args:
        accession: A UniProtKB accession, e.g. ``"P02768"``.
        protease: **Read this if you are coming from MaxQuant or Mascot.** mzLib's ``"trypsin|P"``
            applies the classic Keil rule — cleave after K/R *except* before proline — and is the
            default here because it is what a mass spectrometrist usually means. mzLib's plain
            ``"trypsin"`` cleaves before proline too. That is the **reverse** of the MaxQuant and
            Mascot convention, where ``Trypsin/P`` denotes ignoring the proline rule. On serum
            albumin the two differ by 7 peptides out of about 200 (195 vs 202; tryptic, 2
            missed cleavages, min length 7) - a small count hiding a large semantic
            difference, since which peptides you get changes wherever a K/R precedes a
            proline.
        dissociation: ``"HCD"``/``"CID"`` (b and y ions), ``"ETD"``, and the rest of mzLib's
            dissociation types.

            ``"ETD"`` returns the c and zDot series (radical N-Ca cleavage yields c/z*, not the
            b/y of vibrational activation). Note that z* ions are suppressed N-terminal to
            proline while the complementary c ions are not, leaving about 4% of the c series
            unobservable (smith-chem-wisc/mzLib#1110).
        modifications: Apply UniProt's annotated modifications.

            Pass ``False`` for the bare sequences — a clean control, since the digest's
            distinct backbones are unchanged and only the modified variants of them go
            away. The peptidoform count drops a long way; on albumin from 303 to 195.

            This once carried a caveat saying ``False`` was *not* a clean control, because
            it discarded UniProt's whole feature table and with it the signal-peptide and
            propeptide boundaries mzLib digests at — costing albumin two peptides and
            changing the peptide list, not just the modifications on it (issue #8). That
            was true and is no longer: verified against the published bridge, both peptides
            are present either way and albumin gives 195 distinct base sequences with
            modifications on or off.
        missed_cleavages: Maximum missed cleavage sites per peptide.
        min_length: Shortest peptide to keep, in residues. The default of 7 silently discards shorter
            peptides — roughly a third of a histone digest — so pass ``min_length=1`` when you
            mean *every* peptide.
        max_length: Longest peptide to keep, in residues. ``None`` means unbounded.
        max_modifications: Maximum modifications considered per peptide. Modification isoforms
            are enumerated combinatorially: histone H3.1 yields 49 bare tryptic peptides, 2,563
            at two modifications and 7,040 at three.
        max_isoforms: Maximum modification isoforms per peptide position. mzLib's default of 1024
            **truncates silently** when it binds — on H3.1 at four modifications it discards
            about 30% of the Peptidoforms (13,700 down to 9,536). :attr:`Digest.peptides_at_cap`
            reports how many peptides hit it, so a truncated answer is visible rather than
            merely short.
        terminus: ``"Both"``, ``"N"`` or ``"C"``.
        timeout: Seconds to allow. Large proteins with many modification isoforms take longer.

    Returns:
        A :class:`Digest`. Check :attr:`Digest.modification_census` before trusting a
        modification count — it reports what was annotated as well as what was applied.

    Raises:
        UsageError: the accession, protease, dissociation type or terminus is not recognised.
        ServiceUnavailableError: UniProt was unreachable.

    Examples:
        >>> d = fragments("P02768")
        >>> print(d.modification_census.explain())                     # doctest: +ELLIPSIS
        14 of 38 annotated modifications were applied, across 14 residue positions. Excluded ...
    """
    if not isinstance(accession, str) or not accession.strip():
        raise _bridge.UsageError("A UniProt accession is required, e.g. 'P02768'.")
    canonical = accession.strip().upper()
    if not _ACCESSION.match(canonical):
        raise _bridge.UsageError(
            f"'{accession}' is not a valid UniProtKB accession. They look like 'P02768' or "
            "'A0A0B4J2D5' — see https://www.uniprot.org/help/accession_numbers."
        )
    # UniProt's accessions are upper-case and its API is case-sensitive, so the validated
    # canonical form is what gets sent, not the caller's original casing.
    if max_length is not None and (isinstance(max_length, bool) or not isinstance(max_length, int)):
        raise _bridge.UsageError(f"max_length must be a whole number or None; got {max_length!r}.")
    for name, value in (("missed_cleavages", missed_cleavages), ("min_length", min_length),
                        ("max_modifications", max_modifications), ("max_isoforms", max_isoforms)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise _bridge.UsageError(f"{name} must be a non-negative whole number; got {value!r}.")
    if max_isoforms < 1:
        raise _bridge.UsageError(f"max_isoforms must be at least 1; got {max_isoforms}.")

    args = [
        "peptidoform", "fragments",
        "--accession", canonical,
        "--protease", protease,
        "--dissociation", dissociation,
        "--terminus", terminus,
        "--missed-cleavages", str(missed_cleavages),
        "--min-length", str(min_length),
        "--max-length", str(max_length if max_length else 0),
        "--max-mods", str(max_modifications),
        "--max-isoforms", str(max_isoforms),
    ]
    if not modifications:
        args.append("--no-modifications")

    data = _bridge.invoke(*args, timeout=timeout)

    census = ModificationCensus(
        sites=int(data.get("annotated_modification_sites", 0)),
        applied=int(data.get("annotated_modifications_loaded", 0)),
        annotated=int(data.get("uniprot_annotated_features", 0)),
        by_type=list(data.get("uniprot_features_by_type") or []),
        unresolved=list(data.get("unresolved_modifications") or []),
    )

    return Digest(
        accession=data.get("accession", ""),
        name=data.get("name", "") or "",
        full_name=data.get("full_name", "") or "",
        organism=data.get("organism", "") or "",
        sequence_length=int(data.get("sequence_length", 0)),
        protease=data.get("protease", ""),
        dissociation=data.get("dissociation", ""),
        terminus=data.get("terminus", ""),
        modifications_applied=bool(data.get("modifications_applied", False)),
        max_modifications=int(data.get("max_modifications", 0)),
        max_isoforms=int(data.get("max_modification_isoforms", 0)),
        peptides_at_cap=int(data.get("peptides_at_isoform_cap", 0)),
        modification_census=census,
        peptides=[Peptide._from_wire(p) for p in (data.get("peptides") or [])],
    )


# ---- convert -----------------------------------------------------------------------------------

#: The first pyMzLib whose bridge has ``peptidoform convert`` (its spec's ``since.pymzlib``).
_CONVERT_SINCE = "0.4.0"

#: mzLib's ``SequenceConversionHandlingMode`` names, the values ``mode=`` accepts.
CONVERSION_MODES = ("ThrowException", "ReturnNull", "RemoveIncompatibleElements", "UsePrimarySequence")


@dataclass(frozen=True)
class ConvertedSequence:
    """One input sequence and what mzLib made of it.

    Attributes:
        input: The sequence exactly as given.
        output: The sequence in the target notation, or ``None`` when mzLib could not convert it
            (``status == "failed"``).
        status: mzLib's verdict. ``"converted"``: an output and nothing recorded against it.
            ``"converted_with_warnings"``: an output, but mzLib recorded a warning, an error or an
            incompatible item (a modification dropped, a character skipped); read
            :attr:`incompatible_items` and :attr:`warnings`. ``"failed"``: no output.
        failure_reason: mzLib's ``ConversionFailureReason`` (``"InvalidSequence"``,
            ``"IncompatibleModifications"``, ``"UnsupportedDirection"``, ``"UnknownFormat"``), or
            ``None`` when mzLib recorded none. A failed row can have none: the Unimod target under
            ``ReturnNull`` names only the :attr:`incompatible_items`.
        incompatible_items: The modifications (or other elements) mzLib could not write in the
            target, as mzLib describes them (``"Made Up:Not a modification on K @3(K)"``).
        warnings: mzLib's non-fatal messages for this input.
        errors: mzLib's error messages for this input.
    """

    input: str
    output: str | None
    status: str
    failure_reason: str | None
    incompatible_items: list[str]
    warnings: list[str]
    errors: list[str]

    @property
    def ok(self) -> bool:
        """Whether mzLib converted it with nothing recorded against it (``status == "converted"``)."""
        return self.status == "converted"


@dataclass(frozen=True)
class SequenceConversions(_Table):
    """What :func:`convert` returned: one row per input sequence, in input order.

    ``columns`` maps ``input``, ``output``, ``status``, ``failure_reason``, ``incompatible_items``,
    ``warnings`` and ``errors`` to one value per input, ready for
    ``pandas.DataFrame(result.columns)``. :attr:`sequences` is the same rows as
    :class:`ConvertedSequence` objects.

    Attributes:
        source_format: The source notation, spelled as mzLib registered it (``"mzLib"``).
        target_format: The target notation, spelled as mzLib registered it (``"Unimod"``).
        mode: The ``SequenceConversionHandlingMode`` used.
        source_formats: Every source notation mzLib has registered, sorted.
        target_formats: Every target notation mzLib has registered, sorted.
        record_count: Rows, one per input sequence, duplicates included.
        converted_count: Sequences with status ``"converted"``.
        warned_count: Sequences with status ``"converted_with_warnings"``.
        failed_count: Sequences with status ``"failed"``.
        column_names: The table's columns, in order.
        columns: Column name to one value per input.
        caveats: What a status does and does not promise, per target. Read these once.
    """

    source_format: str
    target_format: str
    mode: str
    source_formats: list[str]
    target_formats: list[str]
    record_count: int
    converted_count: int
    warned_count: int
    failed_count: int
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    @property
    def outputs(self) -> list[str | None]:
        """The converted sequences in input order, ``None`` where a row failed."""
        return list(self.columns.get("output") or [])

    @property
    def sequences(self) -> list[ConvertedSequence]:
        """Every row as a :class:`ConvertedSequence`, in input order."""
        return [
            ConvertedSequence(
                input=row["input"],
                output=row["output"],
                status=row["status"],
                failure_reason=row["failure_reason"],
                incompatible_items=list(row["incompatible_items"] or []),
                warnings=list(row["warnings"] or []),
                errors=list(row["errors"] or []),
            )
            for row in self.records
        ]

    @property
    def not_converted(self) -> list[ConvertedSequence]:
        """The rows whose status is not ``"converted"``: failed, or converted with warnings."""
        return [s for s in self.sequences if not s.ok]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SequenceConversions":
        return cls(
            source_format=data.get("source_format", ""),
            target_format=data.get("target_format", ""),
            mode=data.get("mode", ""),
            source_formats=list(data.get("source_formats") or []),
            target_formats=list(data.get("target_formats") or []),
            record_count=int(data.get("record_count", 0)),
            converted_count=int(data.get("converted_count", 0)),
            warned_count=int(data.get("warned_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns={k: list(v) for k, v in (data.get("columns") or {}).items()},
            caveats=list(data.get("caveats") or []),
        )


def _sequence_lines(sequences: Iterable[str]) -> list[str]:
    """The sequences as stdin lines, unchanged, refusing what would shift the rows."""
    if isinstance(sequences, str):
        raise _bridge.UsageError(
            "sequences must be a list of strings, not one string: that would send one sequence "
            f"per character. Wrap it: [{sequences!r}]."
        )
    lines = []
    for value in sequences:
        if not isinstance(value, str):
            raise _bridge.UsageError(f"Every sequence must be a str; got {value!r}.")
        # The bridge skips blank lines and splits on line breaks, so either would put every later
        # result on the wrong row. Refused here rather than silently re-aligned.
        if not value.strip():
            raise _bridge.UsageError("sequences contains a blank entry; every entry is one result row.")
        if "\n" in value or "\r" in value:
            raise _bridge.UsageError(f"A sequence contains a line break: {value!r}.")
        lines.append(value)
    return lines


def convert(
    sequences: Iterable[str],
    *,
    source: str = "mzLib",
    target: str = "Unimod",
    mode: str = "ReturnNull",
    threads: int = 1,
    timeout: float | None = 120,
) -> SequenceConversions:
    """Convert full sequences from one notation to another with mzLib, one result per input.

    Wraps mzLib's ``SequenceConversionService.Default`` (with ProForma registered). The main use is
    MetaMorpheus or mzLib full sequences to Unimod accessions: ``[UniProt:N-acetylserine on S]SEQK``
    becomes ``[UNIMOD:1]SEQK``. Nothing is converted in Python or in the bridge; every output and
    every status is mzLib's.

    **Choose Unimod, not ProForma, for UniProt-sourced modifications.** mzLib's ProForma target
    does not resolve them and writes them back under their mzLib name, status ``"converted"``. See
    the caveats.

    Args:
        sequences: Full sequences, one result row each, in order, duplicates included. Sent to
            mzLib exactly as given. A blank entry or one with a line break is refused, because it
            would shift every later row. Split an ambiguous MetaMorpheus full sequence
            (``|``-joined candidates) first: mzLib joins the candidates into one sequence.
        source: The notation the sequences are in. One of mzLib's registered source formats
            (:attr:`SequenceConversions.source_formats`), matched case-insensitively.
        target: The notation to write. One of mzLib's registered target formats
            (:attr:`SequenceConversions.target_formats`), matched case-insensitively.
        mode: mzLib's ``SequenceConversionHandlingMode``, one of :data:`CONVERSION_MODES`
            (case-insensitive).
            ``"ReturnNull"`` (default): a sequence mzLib cannot convert is a failed row.
            ``"RemoveIncompatibleElements"`` and ``"UsePrimarySequence"``: what the target cannot
            write is dropped and the row is ``"converted_with_warnings"``. ``"ThrowException"``:
            the first such sequence raises :class:`UsageError` naming it, and nothing is returned.
        threads: Sequences converted at once. Default 1; ``-1`` means every core. Same rows, in
            the same order, at any value.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`SequenceConversions` table, also as :attr:`SequenceConversions.sequences`.

    Raises:
        UsageError: no sequences; a blank one; ``source`` or ``target`` not a format mzLib
            registered (the message lists them); ``mode`` not one of :data:`CONVERSION_MODES`; or,
            under ``mode="ThrowException"``, a sequence mzLib could not convert.

    Examples:
        UniProt and MetaMorpheus modifications to Unimod accessions, and one mzLib cannot map:

        >>> result = convert([
        ...     "[UniProt:N-acetylserine on S]SEQK",
        ...     "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
        ...     "PEPM[Common Variable:Oxidation on M]K",
        ...     "PEPK[Made Up:Not a modification on K]R",
        ... ])
        >>> result.outputs
        ['[UNIMOD:1]SEQK', 'PEPK[UNIMOD:36]R', 'PEPM[UNIMOD:35]K', None]
        >>> failed = result.not_converted[0]
        >>> failed.status, failed.incompatible_items
        ('failed', ['Made Up:Not a modification on K @3(K)'])
    """
    lines = _sequence_lines(sequences)
    if not lines:
        raise _bridge.UsageError(
            "At least one sequence is required, e.g. ['[UniProt:N-acetylserine on S]SEQK']."
        )
    for name, value in (("source", source), ("target", target)):
        if not isinstance(value, str) or not value.strip():
            raise _bridge.UsageError(f"{name} must be a non-empty format name; got {value!r}.")
    canonical_mode = {m.lower(): m for m in CONVERSION_MODES}.get(mode.lower()) if isinstance(mode, str) else None
    if canonical_mode is None:
        raise _bridge.UsageError(
            f"mode must be one of {', '.join(CONVERSION_MODES)} (mzLib's "
            f"SequenceConversionHandlingMode); got {mode!r}."
        )
    if isinstance(threads, bool) or not isinstance(threads, int) or threads == 0 or threads < -1:
        raise _bridge.UsageError(f"threads must be 1 or more, or -1 for every core; got {threads!r}.")

    _bridge.require_verb("peptidoform convert", since=_CONVERT_SINCE)
    data = _bridge.invoke(
        "peptidoform", "convert",
        "--from", source, "--to", target, "--mode", canonical_mode, "--threads", str(threads),
        stdin="\n".join(lines) + "\n",
        timeout=timeout,
    )
    return SequenceConversions._from_wire(data)
