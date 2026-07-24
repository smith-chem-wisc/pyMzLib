"""Label-free quantification with FlashLFQ: quantify a search's peptides across mzML runs.

The question this answers is the one a quant workflow actually asks — *given these identifications
and these runs, how much of each peptide and protein is in each run?* — in one call:

    >>> import pymzlib
    >>> result = pymzlib.flashlfq.quantify(                          # doctest: +SKIP
    ...     psms="AllPSMs.psmtsv",
    ...     spectra=["run_3.mzML", "run_4.mzML"],
    ...     match_between_runs=True,
    ... )
    >>> result.peptide_count, result.protein_count                   # doctest: +SKIP
    (354, 943)

The whole pipeline is mzLib's: the result file is read by mzLib's ``Readers``, turned into FlashLFQ
identifications by mzLib's own converter, and quantified by ``FlashLfqEngine``. MetaMorpheus is not
involved — mzLib does it alone.

Names follow mzLib and FlashLFQ deliberately, so a value here means the same thing it does in the
FlashLFQ source, the MetaMorpheus output columns, and the FlashLFQ paper: ``match_between_runs``,
``ppm_tolerance``, ``mbr_ppm_tolerance``, ``sequence``, ``base_sequence``, ``protein_groups``,
``detection_types``, ``ProteinGroup``, ``FlashLfqResults``.

.. warning::

   **Do not quantify an MSFragger ``psm.tsv`` with this.** mzLib accepts one — it is one of only
   three formats implementing the record view FlashLFQ consumes — but the numbers that come back
   are wrong, and wrong silently.

   MSFragger writes retention time in **seconds**; MetaMorpheus writes **minutes**. mzLib's
   result-file readers pass the column through without converting it, and FlashLFQ then treats the
   value as minutes and searches a ±2-minute window around it. A peptide truly eluting at 60
   minutes is written as ``3600`` and then looked for at 3600 minutes - far past the end of any
   real gradient - so it is simply not found: quantification collapses toward zero rather than
   failing. This is an mzLib defect, not a pyMzLib one, and it affects any mzLib
   caller — it is reported upstream. Until it is fixed, quantify MetaMorpheus output only.

   ``pymzlib.readers.identify()`` will still report an MSFragger file as ``quantifiable``: that
   reports mzLib's interface, which the file genuinely implements. It is not an endorsement of the
   numbers.

**Three more limits worth knowing before you trust a number**, in the "surface it, don't hide it"
spirit of the rest of pyMzLib:

- **mzML only, for now.** Convert ``.raw``/``.d`` to mzML first; a non-mzML path is rejected up front.
- **A protein intensity can be ``None``.** FlashLFQ's median-polish protein quant marks a protein
  NaN when its peptide matrix is degenerate (too few peptides per run, or identical intensities
  across runs — a real artifact documented in mzLib's own tests). NaN is not valid JSON, so it
  arrives here as ``None`` — "could not be quantified" — rather than a silently wrong number. A
  *peptide* intensity, by contrast, is ``0.0`` when missing, never ``None``.
- **For match-between-runs, read the peaks, not the peptides.** The peptide roll-up
  (:attr:`FlashLfqResults.peptides`, mirroring ``QuantifiedPeptides.tsv``) reports far fewer MBR
  transfers than actually happened — a whole run's transfers can be absent. :attr:`FlashLfqResults.peaks`
  (and :meth:`FlashLfqResults.mbr_peaks`) is the complete surface; :attr:`FlashLfqResults.mbr_peak_count`
  is the number to trust.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence, Union

from . import _bridge

__all__ = [
    "SpectraFileInfo",
    "Peptide",
    "ProteinGroup",
    "Peak",
    "FlashLfqResults",
    "quantify",
]

#: One entry in the ``spectra`` argument: either a bare mzML path, or a mapping carrying the
#: experimental-design fields (``path`` plus any of ``condition``, ``biological_replicate``,
#: ``technical_replicate``, ``fraction``). The mapping keys are FlashLFQ's ``SpectraFileInfo`` names.
SpectraInput = Union[str, Mapping[str, Any]]


@dataclass(frozen=True)
class SpectraFileInfo:
    """One quantified run, mirroring mzLib's ``MassSpectrometry.SpectraFileInfo``.

    Attributes:
        file_name: The run's base name (no directory, no extension) — the key used everywhere else
            here to look up this run's intensity.
        full_path: The mzML path as provided.
        condition: The sample-group label, or ``""`` if none was given.
        biological_replicate / technical_replicate / fraction: The experimental-design coordinates.
        peak_count: Chromatographic peaks quantified in this run.
        mbr_peak_count: Of those, how many were transferred by match-between-runs — peaks quantified
            in this run for a peptide that was never identified *in* it. Zero unless
            ``match_between_runs`` was set.
    """

    file_name: str
    full_path: str
    condition: str
    biological_replicate: int
    technical_replicate: int
    fraction: int
    peak_count: int
    mbr_peak_count: int

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "SpectraFileInfo":
        return cls(
            file_name=payload.get("file_name", ""),
            full_path=payload.get("full_path", ""),
            condition=payload.get("condition", ""),
            biological_replicate=int(payload.get("biological_replicate", 0)),
            technical_replicate=int(payload.get("technical_replicate", 0)),
            fraction=int(payload.get("fraction", 0)),
            peak_count=int(payload.get("peak_count", 0)),
            mbr_peak_count=int(payload.get("mbr_peak_count", 0)),
        )


@dataclass(frozen=True)
class Peptide:
    """A quantified peptide, mirroring FlashLFQ's ``Peptide``.

    Attributes:
        sequence: The full (modified) sequence, as FlashLFQ renders it — the identity FlashLFQ
            quantifies. Two different modification states of one base sequence are two peptides.
        base_sequence: The bare amino-acid sequence.
        protein_groups: The protein group(s) this peptide belongs to, ``;``-joined.
        intensities: Run base name → intensity in that run. Missing is ``0.0``, **never ``None``** — a ``null``
            on the wire is resolved to ``0.0`` when parsed, so the invariant holds by
            construction rather than by hope (issue #7). Note too that where a peptide has
            several peaks in one run this roll-up reports **one** of them rather than their
            sum, so pivoting :attr:`FlashLfqResults.peaks` yourself will not reproduce these
            intensities exactly. Counting presence is unaffected; summing intensity is not.
            (unlike proteins). This roll-up mirrors FlashLFQ's ``QuantifiedPeptides.tsv`` and does
            **not** fully reflect match-between-runs — many transferred peptides read ``0.0`` /
            ``"NotDetected"`` here. For MBR-inclusive quantities use :attr:`FlashLfqResults.peaks`.
        detection_types: Run base name → how it was quantified there. Values FlashLFQ emits:
            ``"MSMS"`` (identified and quantified in this run), ``"MBR"`` (transferred by
            match-between-runs), ``"MSMSIdentifiedButNotQuantified"`` (an ID here but no usable
            peak), ``"MSMSAmbiguousPeakfinding"`` (more than one peptide fits the peak), and
            ``"NotDetected"``.
    """

    sequence: str
    base_sequence: str
    protein_groups: str
    intensities: dict[str, Any]
    detection_types: dict[str, str]

    def intensity(self, file_name: str) -> Any:
        """This peptide's intensity in the named run.

        **``0.0`` — not ``None`` — means "not quantified here."** (Only *protein* intensities are
        ever ``None``.) Treat ``0.0`` as missing, not as a measured absence: log-transforming it
        will mislead. And note a peptide that FlashLFQ *transferred* into this run by
        match-between-runs may still read ``0.0`` here — see :attr:`FlashLfqResults.peaks`. Returns
        ``0.0`` for a run that was never provided.
        """
        return self.intensities.get(file_name, 0.0)

    def detection_type(self, file_name: str) -> str:
        """How this peptide was quantified in the named run (``"NotDetected"`` if it was not)."""
        return self.detection_types.get(file_name, "NotDetected")

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "Peptide":
        return cls(
            sequence=payload.get("sequence", ""),
            base_sequence=payload.get("base_sequence", ""),
            protein_groups=payload.get("protein_groups", ""),
            # A null crossing the wire becomes 0.0 here, so the documented invariant --
            # "missing is 0.0, never None" -- is true by construction rather than by hope.
            # The bridge routes every double through a finite check, so a NaN peptide
            # intensity arrives as null; leaving it None would make intensity() return None
            # for a peptide, which is supposed to be a protein-only condition and is exactly
            # what callers branch on to find unresolvable proteins. See issue #7.
            intensities={
                run: (0.0 if value is None else value)
                for run, value in (payload.get("intensities") or {}).items()
            },
            detection_types=dict(payload.get("detection_types") or {}),
        )


@dataclass(frozen=True)
class ProteinGroup:
    """A quantified protein group, mirroring FlashLFQ's ``ProteinGroup``.

    Attributes:
        protein_group: The protein group name (accession, or ``;``-joined accessions).
        gene_name: The gene name, when the result file carried one.
        organism: The organism, when the result file carried one.
        intensities: Run base name → protein intensity in that run. **May be ``None``**: FlashLFQ's
            median-polish protein quant emits NaN (returned here as ``None``) when the peptide matrix
            for the protein is degenerate — too few peptides per run to resolve, or several runs
            reporting the same intensity — protecting you from a fabricated-looking number. Missing,
            as opposed to un-resolvable, is ``0.0``.
    """

    protein_group: str
    gene_name: str
    organism: str
    intensities: dict[str, Any]

    def intensity(self, file_name: str) -> Any:
        """This protein's intensity in the named run.

        ``None`` means FlashLFQ could not resolve a number (a degenerate peptide matrix); ``0.0``
        means simply not measured in this run.
        """
        return self.intensities.get(file_name, 0.0)

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "ProteinGroup":
        return cls(
            protein_group=payload.get("protein_group", ""),
            gene_name=payload.get("gene_name", ""),
            organism=payload.get("organism", ""),
            intensities=dict(payload.get("intensities") or {}),
        )


@dataclass(frozen=True)
class Peak:
    """One quantified chromatographic peak, mirroring FlashLFQ's ``ChromatographicPeak``.

    This is the surface to use for **match-between-runs**. Unlike the peptide roll-up
    (:attr:`Peptide.intensities`, which mirrors ``QuantifiedPeptides.tsv`` and drops most MBR
    transfers), the peaks fully represent every quantified peak, transferred or not. To build an
    MBR-inclusive peptide × run matrix, pivot these on ``(sequence, file_name)``.

    Attributes:
        file_name: The run this peak was measured in (base name).
        sequence: The full (modified) sequence of the peptide the peak was assigned to.
        base_sequence: The bare amino-acid sequence.
        intensity: The peak's intensity (``None`` only if FlashLFQ could not integrate it).
        detection_type: ``"MSMS"``, ``"MBR"``, ``"MSMSIdentifiedButNotQuantified"``,
            ``"MSMSAmbiguousPeakfinding"`` — the same vocabulary as :attr:`Peptide.detection_types`.
            Filter ``detection_type == "MBR"`` to see exactly the transferred peaks.
        retention_time: Apex retention time in minutes, or ``None`` if the peak has no apex.
        num_identifications: How many peptides could explain this peak. ``> 1`` means it is
            ambiguous and its intensity should be treated with care.
        protein_groups: The protein group(s) the assigned identification(s) belong to, ``;``-joined.
    """

    file_name: str
    sequence: str
    base_sequence: str
    intensity: Any
    detection_type: str
    retention_time: Any
    num_identifications: int
    protein_groups: str

    @property
    def is_mbr(self) -> bool:
        """Whether this peak was transferred by match-between-runs."""
        return self.detection_type == "MBR"

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "Peak":
        return cls(
            file_name=payload.get("file_name", ""),
            sequence=payload.get("sequence", ""),
            base_sequence=payload.get("base_sequence", ""),
            intensity=payload.get("intensity"),
            detection_type=payload.get("detection_type", ""),
            retention_time=payload.get("retention_time"),
            num_identifications=int(payload.get("num_identifications", 0)),
            protein_groups=payload.get("protein_groups", ""),
        )


@dataclass(frozen=True)
class FlashLfqResults:
    """The result of a quantification run, mirroring mzLib's ``FlashLfqResults``.

    Attributes:
        psm_file: The absolute path of the PSM result file that was quantified.
        identification_count: How many identifications were read from it.
        parameters: The FlashLFQ parameters actually used, echoed back with their mzLib names.
        spectra_files: One :class:`SpectraFileInfo` per run.
        peptides: One :class:`Peptide` per quantified modified sequence (mirrors
            ``QuantifiedPeptides.tsv``; does not fully carry MBR — use :attr:`peaks` for that).
        proteins: One :class:`ProteinGroup` per quantified protein group.
        peaks: Every quantified :class:`Peak` across all runs (mirrors ``QuantifiedPeaks.tsv``) —
            the complete surface, and the one to use for match-between-runs.
        output_directory: Where the FlashLFQ TSVs were written, or ``None`` if none were.
    """

    psm_file: str
    identification_count: int
    parameters: dict[str, Any]
    spectra_files: list[SpectraFileInfo]
    peptides: list[Peptide]
    proteins: list[ProteinGroup]
    peaks: list[Peak]
    output_directory: Any

    @property
    def peptide_count(self) -> int:
        """Number of quantified peptides."""
        return len(self.peptides)

    @property
    def protein_count(self) -> int:
        """Number of quantified protein groups."""
        return len(self.proteins)

    @property
    def mbr_peak_count(self) -> int:
        """Total match-between-runs peaks (transfers) across every run.

        This counts transferred **peaks**, not distinct peptides: one peptide rescued in two runs is
        two peaks here. For "how many peptides did MBR rescue," use :attr:`mbr_rescued_peptide_count`.
        Either way, do not count MBR from the peptide roll-up (:attr:`Peptide.detection_types`) — it
        under-counts. Zero unless ``match_between_runs`` was on.
        """
        return sum(f.mbr_peak_count for f in self.spectra_files)

    @property
    def mbr_peaks(self) -> list[Peak]:
        """Exactly the peaks transferred by match-between-runs (``detection_type == "MBR"``)."""
        return [p for p in self.peaks if p.is_mbr]

    @property
    def mbr_rescued_peptide_count(self) -> int:
        """Distinct peptides quantified in at least one run only by match-between-runs.

        Exactly: **the number of distinct ``sequence`` values among peaks whose
        ``detection_type`` is ``"MBR"``.** Stated in code terms because the prose version,
        "peptides quantified in at least one run *only* by match-between-runs", is subtly
        different and on real data the two diverge: peptides having **both** an MBR peak and
        a zero-intensity ``MSMS`` peak in the same run were identified there, so they are not
        rescues under the strict reading. On the K562 pair this returns 140 where the strict
        count is 135. Do not read ``mbr_rescued_peptide_count == mbr_peak_count`` as
        reassurance that nothing was double-counted; on that data both are 140, and they
        coincide only because every MBR peak happened to carry a distinct sequence.

        Distinct modified sequences among
        :attr:`mbr_peaks`. This equals :attr:`mbr_peak_count` only when no peptide was rescued in
        more than one run.
        """
        return len({p.sequence for p in self.mbr_peaks})

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "FlashLfqResults":
        return cls(
            psm_file=data.get("psm_file", ""),
            identification_count=int(data.get("identification_count", 0)),
            parameters=dict(data.get("parameters") or {}),
            spectra_files=[SpectraFileInfo._from_wire(f) for f in (data.get("spectra_files") or [])],
            peptides=[Peptide._from_wire(p) for p in (data.get("peptides") or [])],
            proteins=[ProteinGroup._from_wire(g) for g in (data.get("proteins") or [])],
            peaks=[Peak._from_wire(pk) for pk in (data.get("peaks") or [])],
            output_directory=data.get("output_directory"),
        )


def _spectra_stdin(spectra: Sequence[SpectraInput]) -> str:
    """Render the spectra argument as the bridge's tab-separated stdin lines.

    The wire format is one run per line: ``path[\\tcondition[\\tbiorep[\\ttechrep[\\tfraction]]]]``.
    Omitted trailing design fields are simply not written; the bridge applies the same defaults
    MetaMorpheus does with no experimental-design file (blank condition, each run its own
    biological replicate, fraction 0, technical replicate 0).
    """
    if not isinstance(spectra, Sequence) or isinstance(spectra, (str, bytes)):
        raise _bridge.UsageError("spectra must be a list of mzML paths or design mappings.")
    if len(spectra) == 0:
        raise _bridge.UsageError("At least one spectra file is required.")

    lines: list[str] = []
    for index, item in enumerate(spectra):
        if isinstance(item, str):
            path, design = item, {}
        elif isinstance(item, Mapping):
            path = item.get("path")
            if not path:
                raise _bridge.UsageError(f"spectra[{index}] is a mapping with no 'path'.")
            design = item
        else:
            raise _bridge.UsageError(
                f"spectra[{index}] must be an mzML path string or a mapping, got {type(item).__name__}."
            )

        path = str(path)
        if "\t" in path or "\n" in path:
            raise _bridge.UsageError(f"spectra[{index}] path may not contain a tab or newline.")

        fields = [
            path,
            str(design.get("condition", "")),
            _design_int(design, "biological_replicate", index),
            _design_int(design, "technical_replicate", index),
            _design_int(design, "fraction", index),
        ]
        # Drop trailing empties so a bare path stays a bare path on the wire.
        while len(fields) > 1 and fields[-1] == "":
            fields.pop()
        lines.append("\t".join(fields))

    return "\n".join(lines) + "\n"


def _design_int(design: Mapping[str, Any], key: str, index: int) -> str:
    value = design.get(key)
    if value is None:
        return ""
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise _bridge.UsageError(f"spectra[{index}] {key} must be a non-negative whole number; got {value!r}.")
    return str(value)


def quantify(
    psms: str,
    spectra: Sequence[SpectraInput],
    *,
    normalize: bool = False,
    ppm_tolerance: float = 10.0,
    isotope_ppm_tolerance: float = 5.0,
    integrate: bool = False,
    match_between_runs: bool = False,
    mbr_ppm_tolerance: float = 10.0,
    mbr_q_value_threshold: float = 0.05,
    use_shared_peptides_for_protein_quant: bool = False,
    bayesian_protein_quant: bool = False,
    use_pep_q_value: bool = False,
    max_threads: int = -1,
    output_directory: str | None = None,
    timeout: float | None = None,
) -> FlashLfqResults:
    """Quantify a search's peptides across mzML runs with FlashLFQ.

    Args:
        psms: Path to a PSM result file. **Use a MetaMorpheus ``.psmtsv`` or ``.osmtsv``.** mzLib
            also accepts an MSFragger ``psm.tsv`` here, but the result is wrong — see the warning
            in the module documentation above. Every run named in the file must have a matching
            mzML in ``spectra``; FlashLFQ matches identifications to runs by base file name.
        spectra: The mzML runs. Each entry is either a path (``"run_1.mzML"``) or a mapping
            carrying the experimental design (``{"path": "run_1.mzML", "condition": "treated",
            "biological_replicate": 1}``). Base file names must be unique.
        normalize: Normalize intensities across runs (FlashLFQ ``Normalize``).
        ppm_tolerance: Mass tolerance for peak-finding, in ppm (``PpmTolerance``).
        isotope_ppm_tolerance: Mass tolerance for isotope-envelope matching, in ppm.
        integrate: Integrate peak intensities rather than taking the apex. FlashLFQ recommends
            leaving this off.
        match_between_runs: Quantify a peptide in a run where it was not identified, by transferring
            the identification from a run where it was (``MatchBetweenRuns``). The reason to reach
            for FlashLFQ; off by default because it makes assumptions worth opting into. **Requires a
            complete, balanced design**: every condition and biological replicate must carry the same
            set of fractions — a missing replicate or fraction breaks the complementarity MBR relies
            on. Read the transferred peaks from :attr:`FlashLfqResults.peaks`, not the peptide table.
        mbr_ppm_tolerance: Mass tolerance for MBR transfers, in ppm.
        mbr_q_value_threshold: The q-value cutoff below which an MBR transfer is accepted.
        use_shared_peptides_for_protein_quant: Let peptides shared between protein groups contribute
            to protein quant (``UseSharedPeptidesForProteinQuant``).
        bayesian_protein_quant: Run FlashLFQ's Bayesian protein-fold-change engine.
        use_pep_q_value: Filter identifications on PEP q-value rather than q-value.
        max_threads: Worker threads; ``-1`` lets FlashLFQ choose.

            **This is not only a performance knob - it changes results.** With ``-1``,
            FlashLFQ's peptide roll-up nondeterministically drops some MBR intensities, so
            peptide and protein numbers vary between runs on byte-identical inputs. On the
            K562 pair, 6 peptides flip between ``0.0`` and a real intensity, which flips a
            borderline protein group between ``None`` and a number: unquantifiable in 5 of 6
            runs, quantified in the 6th. The peaks are stable throughout - only the roll-up
            wobbles. **Set ``max_threads=1`` for anything you intend to publish or
            reproduce.** See smith-chem-wisc/mzLib#1111.
        output_directory: If given, FlashLFQ also writes ``QuantifiedPeaks.tsv``,
            ``QuantifiedPeptides.tsv`` and ``QuantifiedProteins.tsv`` there.
        timeout: Seconds to allow. Large experiments legitimately take a while; ``None`` waits
            indefinitely.

    Returns:
        A :class:`FlashLfqResults`.

    Raises:
        UsageError: an argument is malformed, a run is not mzML, an mzML is missing, or the PSM
            file names a run with no mzML provided.
        BridgeError: FlashLFQ itself failed.
    """
    if not isinstance(psms, str) or not psms.strip():
        raise _bridge.UsageError("A PSM result file path is required, e.g. 'AllPSMs.psmtsv'.")

    stdin = _spectra_stdin(spectra)

    args: list[str] = ["quant", "flashlfq", "--psms", psms]
    if normalize:
        args.append("--normalize")
    args += ["--ppm", _number(ppm_tolerance, "ppm_tolerance")]
    args += ["--isotope-ppm", _number(isotope_ppm_tolerance, "isotope_ppm_tolerance")]
    if integrate:
        args.append("--integrate")
    if match_between_runs:
        args.append("--mbr")
    args += ["--mbr-ppm", _number(mbr_ppm_tolerance, "mbr_ppm_tolerance")]
    args += ["--mbr-q", _number(mbr_q_value_threshold, "mbr_q_value_threshold")]
    if use_shared_peptides_for_protein_quant:
        args.append("--shared-peptides")
    if bayesian_protein_quant:
        args.append("--bayesian")
    if use_pep_q_value:
        args.append("--use-pep-q")
    if isinstance(max_threads, bool) or not isinstance(max_threads, int):
        raise _bridge.UsageError(f"max_threads must be a whole number; got {max_threads!r}.")
    args += ["--threads", str(max_threads)]
    if output_directory is not None:
        if not isinstance(output_directory, str) or not output_directory.strip():
            raise _bridge.UsageError("output_directory must be a non-empty path or None.")
        args += ["--out", output_directory]

    data = _bridge.invoke(*args, stdin=stdin, timeout=timeout)
    return FlashLfqResults._from_wire(data)


def _number(value: Any, name: str) -> str:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _bridge.UsageError(f"{name} must be a number; got {value!r}.")
    if value != value or value in (float("inf"), float("-inf")):  # NaN / infinity
        raise _bridge.UsageError(f"{name} must be a finite number; got {value!r}.")
    return repr(float(value))
