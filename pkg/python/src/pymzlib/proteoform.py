"""Proteoform-level questions: digest an annotated protein and fragment its peptides.

The question this answers is the one a mass spectrometrist actually asks — *what fragments would
I see for this protein's peptides?* — in one call:

    >>> import pymzlib
    >>> digest = pymzlib.proteoform.fragments("P02768")          # doctest: +SKIP
    >>> len(digest.peptides)                                      # doctest: +SKIP
    303

The defaults are opinions, not placeholders. Tryptic with the proline rule, two missed cleavages,
ETD, both termini, UniProt's annotated modifications applied. They are the choices this lab makes
when it does not have a reason to choose otherwise, so the common question needs no parameters —
and every one of them is reachable, because the point is to open the doors, not to hide them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Sequence

from . import _bridge

#: UniProtKB''s own accession grammar (https://www.uniprot.org/help/accession_numbers). Checking
#: it here means a typo costs nothing instead of a network round trip and a puzzling HTTP 400.
_ACCESSION = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

__all__ = ["Fragment", "Peptide", "Digest", "ModificationCensus", "fragments"]


@dataclass(frozen=True)
class Fragment:
    """One backbone fragment ion.

    Attributes:
        product_type: The ion series, e.g. ``"c"`` or ``"zDot"`` for ETD, ``"b"``/``"y"`` for CID.
        fragment_number: Position in the series — ``c3`` is the third from the N-terminus.
        neutral_mass: Monoisotopic neutral mass in daltons. **Not** an m/z: no proton has been
            added and no charge assumed. See :meth:`Peptide.mz` for that conversion, and read
            its note before putting a number in a figure.
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
        modifications: Each applied modification, with its one-based position and mass.
        fragments: The fragment ions for the requested dissociation type.
    """

    base_sequence: str
    full_sequence: str
    monoisotopic_mass: float
    length: int
    one_based_start: int
    one_based_end: int
    missed_cleavages: int
    modifications: list[dict[str, Any]] = field(default_factory=list)
    fragments: list[Fragment] = field(default_factory=list)

    @property
    def is_modified(self) -> bool:
        """Whether this peptide carries at least one modification."""
        return bool(self.modifications)

    def mz(self, charge: int) -> float:
        """Return the m/z of the intact peptide at a given charge.

        Uses the **proton** mass (1.007276), not the hydrogen atom mass (1.007825). The difference
        is 0.55 mDa, which at 500 m/z is 1.1 ppm — on an Orbitrap, the difference between a match
        and a miss. Libraries differ on this and rarely say which they used.

        Args:
            charge: A positive integer charge state.
        """
        if not isinstance(charge, int) or isinstance(charge, bool) or charge < 1:
            raise _bridge.UsageError(f"charge must be a positive whole number; got {charge!r}.")
        return (self.monoisotopic_mass + charge * PROTON_MASS) / charge

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

    A modification is only usable for mass spectrometry if it has a defined target **and** a
    defined mass. A glycosylation site annotated as "N-linked (GlcNAc...) asparagine" has neither
    a fixed composition nor a fixed mass, so there is nothing to add to a peptide, and a peptide
    carrying an ambiguous-mass PTM is not identifiable by mass anyway. Such annotations are
    therefore excluded.

    That exclusion is correct. What this class exists for is that you should not have to *guess*
    it happened: for serum albumin, 14 sites are applied out of 38 annotated, and without this
    the 14 arrives with no indication that a rule was ever applied.

    Attributes:
        sites: Distinct residue positions carrying at least one modification. A histone lists
            several alternatives at one residue — K9me1, K9me2, K9me3, K9ac are four
            modifications at one site — so this is always the smaller number and is not a
            modification count.
        applied: Modifications actually placed on the protein.
        annotated: Modification-like features UniProt lists.
        by_type: One entry per feature type, with ``count`` and whether it was ``loaded``.
    """

    sites: int
    applied: int
    annotated: int
    by_type: list[dict[str, Any]] = field(default_factory=list)

    @property
    def excluded(self) -> int:
        """Annotated features that could not be used because they have no defined mass."""
        return max(0, self.annotated - self.applied)

    def explain(self) -> str:
        """A one-paragraph, human-readable account of what was used and what was not."""
        if not self.excluded:
            return (
                f"All {self.annotated} annotated modifications were applied, across "
                f"{self.sites} residue positions."
            )
        skipped = ", ".join(
            f"{t['count']} × {t['type']}" for t in self.by_type if not t.get("loaded")
        )
        return (
            f"{self.applied} of {self.annotated} annotated modifications were applied, across "
            f"{self.sites} residue positions. "
            f"Excluded: {skipped} — these have no defined chemical composition, so no mass can "
            "be assigned and a peptide carrying one is not identifiable by mass spectrometry."
        )


@dataclass(frozen=True)
class Digest:
    """The result of digesting a protein and fragmenting its peptides."""

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
        treating a proteoform list as exhaustive.
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
            albumin the two differ by 37 peptides out of about 200.
        dissociation: ``"ETD"`` (c and z• ions), ``"HCD"``/``"CID"`` (b and y), and the rest of
            mzLib's dissociation types.
        modifications: Apply UniProt's annotated modifications. Pass ``False`` for the bare
            sequence — useful as a control, and the difference is usually large.
        missed_cleavages: Maximum missed cleavage sites per peptide.
        min_length / max_length: Peptide length bounds. ``max_length=None`` means unbounded.
        max_modifications: Maximum modifications considered per peptide. Modification isoforms
            are enumerated combinatorially: histone H3.1 yields 49 bare tryptic peptides, 2,563
            at two modifications and 7,040 at three.
        max_isoforms: Maximum modification isoforms per peptide position. mzLib's default of 1024
            **truncates silently** when it binds — on H3.1 at four modifications it discards
            about 30% of the proteoforms (13,700 down to 9,536). :attr:`Digest.peptides_at_cap`
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

    Example:
        >>> d = fragments("P02768")                                    # doctest: +SKIP
        >>> print(d.modification_census.explain())                     # doctest: +SKIP
        14 of 38 annotated modification sites were applied. Excluded: 24 × glycosylation site …
    """
    if not isinstance(accession, str) or not accession.strip():
        raise _bridge.UsageError("A UniProt accession is required, e.g. ''P02768''.")
    if not _ACCESSION.match(accession.strip().upper()):
        raise _bridge.UsageError(
            f"''{accession}'' is not a valid UniProtKB accession. They look like ''P02768'' or "
            "''A0A0B4J2D5'' — see https://www.uniprot.org/help/accession_numbers."
        )
    for name, value in (("missed_cleavages", missed_cleavages), ("min_length", min_length),
                        ("max_modifications", max_modifications)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise _bridge.UsageError(f"{name} must be a non-negative whole number; got {value!r}.")

    args = [
        "proteoform", "fragments",
        "--accession", accession.strip(),
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
