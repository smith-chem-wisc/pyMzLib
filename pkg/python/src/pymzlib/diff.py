"""Differential abundance, set up by mzLib: the experiment's design, and every peptide in every sample.

Two questions come before any test, and mzLib answers both (mzLib #1450, #1455):

==============================================================  ======================  ===========================================
You want to know                                                Call                    mzLib computes it with
==============================================================  ======================  ===========================================
What is my design matrix, and which comparisons can it test?    :func:`design`          ``AnalysisDesign.ByStratum``,
                                                                                        ``DefaultContrasts``, ``Estimability``
What did each peptide measure in each sample, and if nothing,   :func:`observations`    ``FlashLfqObservations.FromPeptideTable``,
why? How far apart are the samples' loadings?                                           ``SharedPeptideNormalization.Measure``
==============================================================  ======================  ===========================================

**A sample sheet, not a matrix.** :func:`design` takes the table you already keep - one row per
sample, a column per condition, covariate, batch or tissue - and a model that says which column is
which. mzLib builds the matrix exactly as R's ``model.matrix`` would, says which comparisons the
samples cannot separate (a condition run in only one batch) and why, and writes the design table
:func:`pymzlib.stats.fit` reads::

    >>> design = pymzlib.diff.design("malat_dilution_samples.tsv",
    ...                              factors=[pymzlib.diff.FactorSpec("malat_load", "500ng")],
    ...                              out="malat_design")
    >>> design.strata[0].fitted_column_names
    ['intercept', 'malat_load=125ng', 'malat_load=250ng']

**Missing is never 0.** :func:`observations` reads a FlashLFQ or MetaMorpheus peptide table and
gives every peptide in every biological sample as a log2 intensity or ``None``, with the reason
in ``state`` (not detected, an ambiguous peak, a match-between-runs transfer left out). Fractions
are summed and technical replicates averaged by mzLib's rules, from the runs you list.

**What is not here yet.** mzLib does not yet roll peptides up to proteins for a test, choose which
rows a test uses, or assemble fits into its differential result table. pyMzLib does not do it
either: the peptide-by-sample table ``observations(responses_out=...)`` writes is the observation
table's own shape, nothing more.

Wire verbs: ``diff design`` and ``diff observations``; their language-neutral specs live in the
bridge repository under ``design/verbs/``.
"""

from __future__ import annotations

import math
import numbers
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Union

from . import _bridge
from .readers import WrittenTable, _Table

__all__ = [
    "ESTIMABILITIES",
    "OBSERVATION_STATES",
    "QUANT_BASES",
    "CovariateSpec",
    "Design",
    "DesignStratum",
    "FactorSpec",
    "NormalizationResult",
    "NotRunContrast",
    "ObservationRun",
    "Observations",
    "ProteinGroupInfo",
    "WrittenDesign",
    "WrittenResponses",
    "design",
    "observations",
]

#: The first pyMzLib whose bridge has the diff verbs (``since.pymzlib`` of their specs).
_SINCE = "0.6.0"

#: Every value of :func:`design`'s ``estimability`` column, as mzLib's ``AnalysisDesign`` writes
#: them. Only ``"estimable"`` rows have ``fitted_weights``.
ESTIMABILITIES = ("estimable", "not_estimable:confounded", "not_estimable:rank_deficient")

#: Every value of :func:`observations`' ``state`` column, mzLib's ``ObservationState`` from most to
#: least informative. Only ``"quantified"`` and ``"mbr_transferred"`` rows can carry a value.
OBSERVATION_STATES = (
    "quantified",
    "mbr_transferred",
    "ambiguous_peak",
    "identified_not_quantified",
    "not_detected",
)

#: The values :func:`observations` takes as ``basis=``: mzLib's ``QuantBasis`` machine names. There
#: is no default; you choose.
QUANT_BASES = ("msms_only", "mbr_kept")

_PathLike = Union[str, "os.PathLike[str]"]


def __dir__() -> list[str]:
    """The public API only, so ``dir()`` is a map of what to call rather than of our imports."""
    return sorted(__all__)


def _columns(data: Mapping[str, Any]) -> tuple[list[str], dict[str, list[Any]]]:
    names = list(data.get("column_names") or [])
    return names, dict(data.get("columns") or {})


# ---- inputs ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class FactorSpec:
    """A condition factor of the model: mzLib's ``FactorSpec``.

    Attributes:
        name: The sample sheet column holding each sample's level.
        reference_level: The level every other level is compared against, or ``None`` to compare
            every pair of levels (mzLib then warns, because the comparisons are many).
    """

    name: str
    reference_level: str | None = None


@dataclass(frozen=True)
class CovariateSpec:
    """A numeric covariate of the model, entered as ``(value - centre) / scale``: mzLib's ``CovariateSpec``.

    Leave a field ``None`` and nothing is sent for it, so mzLib's own default applies.

    Attributes:
        name: The sample sheet column holding each sample's value, in ``unit``.
        unit: The stored unit, e.g. ``"years"``.
        scale: Stored units per scaled unit, e.g. ``10`` years per decade. ``None``: mzLib's 1.
        scaled_unit: The unit a slope is per, e.g. ``"decade"``. ``None``: ``scale`` stored units
            (``"years"`` at scale 1, ``"10 years"`` at scale 10).
        centre: The value subtracted first, in ``unit``. ``None``: the stratum's mean.
    """

    name: str
    unit: str
    scale: float | None = None
    scaled_unit: str | None = None
    centre: float | None = None


@dataclass(frozen=True)
class ObservationRun:
    """One spectra file and its place in the design: mzLib's ``ObservationRun``.

    Attributes:
        file_name: The file as the peptide table's ``Intensity_`` columns name it, without its
            extension. A column with ``-calib`` after it (MetaMorpheus's calibrated files) matches.
        sample_id: The biological sample the run belongs to.
        fraction: The fraction; fractions of one technical replicate are summed. ``None`` (as an
            argument): not sent, so mzLib's 0.
        technical_replicate: The technical replicate; a sample's replicates are averaged in log2.
            ``None`` (as an argument): not sent, so mzLib's 0.
    """

    file_name: str
    sample_id: str
    fraction: int | None = None
    technical_replicate: int | None = None


def _name(value: object, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _bridge.UsageError(f"{what} must be a non-empty str; got {value!r}.")
    if any(c in value for c in "\t\r\n"):
        raise _bridge.UsageError(f"{what} contains a tab or a line break, which separate the fields: {value!r}.")
    return value.strip()


def _optional_name(value: object, what: str) -> str:
    return "" if value is None else _name(value, what)


def _number(value: object, what: str) -> str:
    if value is None:
        return ""
    if isinstance(value, bool) or not isinstance(value, numbers.Real) or not math.isfinite(float(value)):
        raise _bridge.UsageError(f"{what} must be a finite number; got {value!r}.")
    return repr(float(value))


def _integer(value: object, what: str) -> str:
    if value is None:
        return ""
    if isinstance(value, bool) or not isinstance(value, int):
        raise _bridge.UsageError(f"{what} must be an int; got {value!r}.")
    return str(value)


def _line(*fields: str) -> str:
    """A stdin line, with trailing fields left out when they were not given."""
    kept = list(fields)
    while kept and kept[-1] == "":
        kept.pop()
    return "\t".join(kept)


def _list(value: object, what: str) -> list:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise _bridge.UsageError(f"{what} takes a list; got {value!r}.")
    return list(value)


def _model_stdin(factors, covariates, interactions, strata, batch) -> str:
    lines = []
    for factor in _list(factors, "factors"):
        if isinstance(factor, str):
            factor = FactorSpec(factor)
        if not isinstance(factor, FactorSpec):
            raise _bridge.UsageError(f"Each factor is a FactorSpec or a column name; got {factor!r}.")
        lines.append(_line("factor", _name(factor.name, "A factor name"),
                           _optional_name(factor.reference_level, "A reference level")))
    for covariate in _list(covariates, "covariates"):
        if not isinstance(covariate, CovariateSpec):
            raise _bridge.UsageError(f"Each covariate is a CovariateSpec; got {covariate!r}.")
        lines.append(_line(
            "covariate",
            _name(covariate.name, "A covariate name"),
            _name(covariate.unit, f"The unit of covariate {covariate.name!r}"),
            _number(covariate.scale, f"The scale of covariate {covariate.name!r}"),
            _optional_name(covariate.scaled_unit, f"The scaled unit of covariate {covariate.name!r}"),
            _number(covariate.centre, f"The centre of covariate {covariate.name!r}"),
        ))
    for pair in _list(interactions, "interactions"):
        if isinstance(pair, str) or not isinstance(pair, Sequence) or len(pair) != 2:
            raise _bridge.UsageError(f"Each interaction is a pair of factor names; got {pair!r}.")
        lines.append(_line("interaction", _name(pair[0], "An interaction's factor"),
                           _name(pair[1], "An interaction's factor")))
    if isinstance(strata, str):
        strata = [strata]
    for column in _list(strata, "strata"):
        lines.append(_line("stratum", _name(column, "A stratum column")))
    if batch is not None:
        lines.append(_line("batch", _name(batch, "The batch column")))
    return "".join(line + "\n" for line in lines)


def _runs_stdin(runs) -> str:
    items = _list(runs, "runs")
    if not items:
        raise _bridge.UsageError("observations() needs the runs: one ObservationRun or (file, sample) tuple per file.")
    lines = []
    for run in items:
        if isinstance(run, ObservationRun):
            fields = (run.file_name, run.sample_id, run.fraction, run.technical_replicate)
        elif isinstance(run, Sequence) and not isinstance(run, str) and 2 <= len(run) <= 4:
            fields = tuple(run) + (None,) * (4 - len(run))
        else:
            raise _bridge.UsageError(
                f"Each run is an ObservationRun or a (file, sample[, fraction[, technical_replicate]]) tuple; got {run!r}."
            )
        lines.append(_line(
            _name(fields[0], "A run's file name"),
            _name(fields[1], "A run's sample"),
            _integer(fields[2], "A run's fraction"),
            _integer(fields[3], "A run's technical replicate"),
        ))
    return "".join(line + "\n" for line in lines)


def _path(value: object, what: str) -> str:
    text = _bridge.path_text(value)
    if not text:
        raise _bridge.UsageError(f"{what} must be a path (str or os.PathLike); got {value!r}.")
    return text


# ---- design results ----------------------------------------------------------------------------


@dataclass(frozen=True)
class NotRunContrast:
    """A contrast the model asks for that one stratum cannot run: mzLib's ``NotRunContrast``.

    Attributes:
        label: Which contrast, e.g. ``"condition (one level: A)"``.
        reason: Why: ``"not_in_stratum"`` (one side has no sample there).
    """

    label: str
    reason: str


@dataclass(frozen=True)
class WrittenDesign:
    """Where :func:`design` wrote its pivot files (``out=``).

    Attributes:
        directory: The absolute directory written into.
        file_count: Design and contrast files written (files).
    """

    directory: str
    file_count: int


@dataclass(frozen=True)
class DesignStratum:
    """One stratum's design, as mzLib's ``AnalysisDesign`` built it.

    Attributes:
        stratum: The stratum key: ``"all"``, or the stratum columns sorted by name as
            ``name=value`` joined by ``;`` (``%``, ``;`` and ``=`` in a value escaped).
        sample_count: Samples in this stratum (samples).
        sample_ids: The matrix rows, ordinal by sample id.
        column_names: mzLib's columns: ``intercept``; ``factor=level`` for each non-reference level;
            ``batch=level``; ``name (per unit)`` for a covariate; ``a=level:b=level`` for an
            interaction.
        fitted_column_names: The columns a fit takes: all of them at full rank, else each kept unless
            it is redundant with those before it. The design file's coefficient columns.
        rank: The numerical rank (columns), by mzLib's ``LinearModel`` full-rank test.
        aliased_terms: Terms that share a direction with another and cannot be told apart from it.
        warnings: What mzLib left out or changed here, in plain words.
        not_run: The contrasts this stratum cannot run, as :class:`NotRunContrast`.
        covariate_centres: The centre subtracted from each covariate here, in the covariate's stored
            unit: the one you gave, else the stratum's mean.
        matrix: The design matrix, one row per ``sample_ids`` entry and one value per
            ``column_names`` entry (coding: 0/1 indicators; covariate columns in scaled units).
            Equal to R's ``model.matrix`` for the same formula and level order.
        fitted_matrix: ``matrix`` with only the fitted columns (same coding): what ``stats.fit``
            fits.
        design_file: The design table written for this stratum (``design_<i>.tsv``), or ``None``
            without ``out=``. Give it to :func:`pymzlib.stats.fit` as ``design``.
        contrasts_file: The contrast weights written for this stratum (``contrasts_<i>.tsv``), one
            row per estimable contrast on the design file's columns; ``None`` without ``out=`` or
            when nothing here is estimable.
    """

    stratum: str
    sample_count: int
    sample_ids: list[str]
    column_names: list[str]
    fitted_column_names: list[str]
    rank: int
    aliased_terms: list[str]
    warnings: list[str]
    not_run: list[NotRunContrast]
    covariate_centres: dict[str, float]
    matrix: list[list[float]]
    fitted_matrix: list[list[float]]
    design_file: str | None
    contrasts_file: str | None

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> DesignStratum:
        return cls(
            stratum=str(data.get("stratum", "")),
            sample_count=int(data.get("sample_count", 0)),
            sample_ids=list(data.get("sample_ids") or []),
            column_names=list(data.get("column_names") or []),
            fitted_column_names=list(data.get("fitted_column_names") or []),
            rank=int(data.get("rank", 0)),
            aliased_terms=list(data.get("aliased_terms") or []),
            warnings=list(data.get("warnings") or []),
            not_run=[NotRunContrast(str(n.get("label", "")), str(n.get("reason", "")))
                     for n in data.get("not_run") or []],
            covariate_centres=dict(data.get("covariate_centres") or {}),
            matrix=[list(row) for row in data.get("matrix") or []],
            fitted_matrix=[list(row) for row in data.get("fitted_matrix") or []],
            design_file=data.get("design_file"),
            contrasts_file=data.get("contrasts_file"),
        )


@dataclass(frozen=True)
class Design(_Table):
    """What :func:`design` returns: each stratum's design, and every default contrast in one table.

    ``columns`` holds one row per (stratum, default contrast), strata in order:

    ==================  ============================  ============================================================
    column              unit                          meaning; ``None`` means
    ==================  ============================  ============================================================
    ``stratum``         -                             the stratum key (:attr:`DesignStratum.stratum`)
    ``contrast_id``     -                             ``c1``, ``c2``, ... numbered **within the stratum**
    ``label``           -                             ``factor=level vs factor=reference``, or a slope's column
    ``numerator``       -                             ``factor=level``; ``None`` = a slope
    ``denominator``     -                             ``factor=reference``; ``None`` = a slope
    ``covariate``       -                             the covariate of a slope; ``None`` = not a slope
    ``covariate_unit``  -                             the unit a slope is per, e.g. ``decade``
    ``covariate_scale`` -                             e.g. ``years/10, centred at 50 years``
    ``estimability``    -                             one of :data:`ESTIMABILITIES`
    ``weights``         weight per design column      one per :attr:`DesignStratum.column_names` entry
    ``fitted_weights``  weight per fitted column      one per :attr:`DesignStratum.fitted_column_names`
                                                      entry, the contrasts file's row; ``None`` = not
                                                      estimable
    ==================  ============================  ============================================================

    Attributes:
        samples_file: Absolute path of the sample sheet.
        sample_count: Samples read, across every stratum (samples).
        factors: The factors given, as :class:`FactorSpec`.
        covariates: The covariates given, as :class:`CovariateSpec`; a ``None`` field was not
            given, so mzLib's default applied.
        interactions: The interactions given, as ``(first, second)`` pairs.
        stratum_columns: The stratum columns given; empty for one stratum.
        batch_column: The batch column, or ``None`` when batch is not in the model.
        unused_columns: Sheet columns the model does not name, ignored.
        stratum_count: Designs built (strata).
        strata: One :class:`DesignStratum` per stratum, in mzLib's order (ordinal by key).
        estimability_counts: Default contrasts per verdict across strata (contrasts), keyed by
            :data:`ESTIMABILITIES`.
        written: Where ``out=`` wrote, as :class:`WrittenDesign`; ``None`` without ``out=``.
        row_count: Default contrasts across strata (rows).
        column_names: The table's columns, in order.
        columns: Column name -> list of values.
        caveats: How the matrix is coded and what the contrast ids mean, plus counts of contrasts
            that cannot be estimated or run.
    """

    samples_file: str
    sample_count: int
    factors: list[FactorSpec]
    covariates: list[CovariateSpec]
    interactions: list[tuple[str, str]]
    stratum_columns: list[str]
    batch_column: str | None
    unused_columns: list[str]
    stratum_count: int
    strata: list[DesignStratum]
    estimability_counts: dict[str, int]
    written: WrittenDesign | None
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    def stratum(self, key: str) -> DesignStratum:
        """One stratum's design, by its key.

        Args:
            key: A :attr:`DesignStratum.stratum` value, e.g. ``"all"`` or ``"organism part=liver"``.

        Raises:
            UsageError: there is no stratum with that key.
        """
        for stratum in self.strata:
            if stratum.stratum == key:
                return stratum
        raise _bridge.UsageError(
            f"No stratum '{key}'; the strata are: " + ", ".join(s.stratum for s in self.strata) + "."
        )

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> Design:
        names, columns = _columns(data)
        written = data.get("written")
        return cls(
            samples_file=str(data.get("samples_file", "")),
            sample_count=int(data.get("sample_count", 0)),
            factors=[FactorSpec(str(f.get("name", "")), f.get("reference_level")) for f in data.get("factors") or []],
            covariates=[
                CovariateSpec(str(c.get("name", "")), str(c.get("unit", "")), c.get("scale"),
                              c.get("scaled_unit"), c.get("centre"))
                for c in data.get("covariates") or []
            ],
            interactions=[(str(p.get("first", "")), str(p.get("second", ""))) for p in data.get("interactions") or []],
            stratum_columns=list(data.get("stratum_columns") or []),
            batch_column=data.get("batch_column"),
            unused_columns=list(data.get("unused_columns") or []),
            stratum_count=int(data.get("stratum_count", 0)),
            strata=[DesignStratum._from_wire(s) for s in data.get("strata") or []],
            estimability_counts=dict(data.get("estimability_counts") or {}),
            written=None if not written else WrittenDesign(str(written.get("directory", "")),
                                                           int(written.get("file_count", 0))),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


# ---- observation results -----------------------------------------------------------------------


@dataclass(frozen=True)
class ProteinGroupInfo:
    """A MetaMorpheus protein group, as its protein-group table describes it: mzLib's ``ProteinGroupInfo``.

    Attributes:
        name: The group's name as MetaMorpheus writes it: member accessions joined by ``|``.
        q_value: The group's ``Protein QValue`` (fraction 0 to 1), or ``None`` if not finite.
        is_decoy: The group is a decoy.
        is_contaminant: The group is a contaminant.
        is_entrapment: The group comes from an entrapment database.
        genes: The ``Gene`` cell as written, or ``None`` when the table has none.
    """

    name: str
    q_value: float | None
    is_decoy: bool
    is_contaminant: bool
    is_entrapment: bool
    genes: str | None


@dataclass(frozen=True)
class NormalizationResult:
    """mzLib's shared-peptide median normalization, as measured: its ``NormalizationResult``.

    Attributes:
        setting: ``"shared_peptide_median"`` (measured on the peptides with a value in every
            sample), ``"shared_peptide_median_half"`` (in at least half; used when too few are in
            all) or ``"none"`` (too few either way: no sample shifted).
        reference_set_size: Peptides in the reference set used (peptides), or in the last one tried
            when ``setting`` is ``"none"``. Decoy, contaminant and entrapment peptides never enter.
        minimum_reference_set_size: The floor used (peptides): ``min_peptides=``, else mzLib's 100.
        shifts: Each sample's shift in log2 units, subtracted with ``normalize=True``; ``None`` for
            a sample with no shift (not shifted).
        warnings: Which samples have no shift, and why.
    """

    setting: str
    reference_set_size: int
    minimum_reference_set_size: int
    shifts: dict[str, float | None]
    warnings: list[str]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> NormalizationResult:
        return cls(
            setting=str(data.get("setting", "")),
            reference_set_size=int(data.get("reference_set_size", 0)),
            minimum_reference_set_size=int(data.get("minimum_reference_set_size", 0)),
            shifts=dict(data.get("shifts") or {}),
            warnings=list(data.get("warnings") or []),
        )


@dataclass(frozen=True)
class WrittenResponses:
    """Where :func:`observations` wrote the peptide-by-sample responses table (``responses_out=``).

    Attributes:
        path: The absolute path written; give it to :func:`pymzlib.stats.fit` as ``responses``.
        format: Always ``"tsv"``.
        row_count: Peptide rows written (peptides).
        sample_count: Sample columns written (samples).
    """

    path: str
    format: str
    row_count: int
    sample_count: int


@dataclass(frozen=True)
class Observations(_Table):
    """What :func:`observations` returns: one row per (peptide, sample), peptides by full sequence.

    ========================  ================================  ==============================================
    column                    unit                              meaning; ``None`` means
    ========================  ================================  ==============================================
    ``full_sequence``         -                                 the engine's modified sequence: the key
    ``base_sequence``         -                                 the unmodified sequence
    ``protein_groups``        -                                 list of every group the peptide table names
                                                                (``UNDEFINED`` = in no surviving group)
    ``unique_protein_group``  -                                 the one group it is unique to; ``None`` =
                                                                shared, UNDEFINED, or none
    ``sample``                -                                 the biological sample
    ``log2_intensity``        log2 of FlashLFQ intensity units  fractions summed, replicates averaged in log2,
                                                                less the shift with ``normalize=True``;
                                                                ``None`` = **no value** (see ``state``), never 0
    ``state``                 -                                 one of :data:`OBSERVATION_STATES`
    ========================  ================================  ==============================================

    Attributes:
        peptides_file: Absolute path of the peptide table.
        protein_groups_file: Absolute path of the protein-group table, or ``None`` without one.
        samples_file: Absolute path of the ``samples=`` list, or ``None`` when every sample was kept.
        basis: ``"msms_only"`` or ``"mbr_kept"``, as given.
        normalize: Whether the shifts were subtracted.
        min_peptides: ``min_peptides=`` as given (peptides), or ``None``: mzLib's floor applied.
        run_count: Runs in the table (runs).
        runs: The runs as mzLib holds them, ordinal by file name, as :class:`ObservationRun`
            (fraction and technical replicate filled in).
        sample_count: Biological samples (samples).
        samples: The samples, ordinal by id, with or without values.
        samples_without_values: Samples in which no peptide has a value under this basis.
        peptide_count: Peptides (peptides): distinct full sequences.
        protein_group_count: Groups the protein-group table describes (protein groups).
        protein_groups: Those groups, as :class:`ProteinGroupInfo`.
        state_counts: Rows per state (cells), keyed by :data:`OBSERVATION_STATES`.
        normalization: The :class:`NormalizationResult`, measured whether or not it was applied.
        warnings: What mzLib noticed that did not stop the build, e.g. groups the peptides name that
            the protein-group table does not describe.
        written: Where ``out=`` wrote the table, as :class:`~pymzlib.readers.WrittenTable`; ``None``
            when the rows are in ``columns``.
        responses_written: Where ``responses_out=`` wrote, as :class:`WrittenResponses`, or ``None``.
        row_count: ``peptide_count`` x ``sample_count`` (rows), returned or written.
        column_names: The table's columns, in order.
        columns: Column name -> list of values; empty with ``out=``.
        caveats: How values are combined, what the basis does, what mzLib cannot check about the
            table, and what the normalization assumes.
    """

    peptides_file: str
    protein_groups_file: str | None
    samples_file: str | None
    basis: str
    normalize: bool
    min_peptides: int | None
    run_count: int
    runs: list[ObservationRun]
    sample_count: int
    samples: list[str]
    samples_without_values: list[str]
    peptide_count: int
    protein_group_count: int
    protein_groups: list[ProteinGroupInfo]
    state_counts: dict[str, int]
    normalization: NormalizationResult
    warnings: list[str]
    written: WrittenTable | None
    responses_written: WrittenResponses | None
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> Observations:
        names, columns = _columns(data)
        written = data.get("written")
        responses = data.get("responses_written")
        minimum = data.get("min_peptides")
        return cls(
            peptides_file=str(data.get("peptides_file", "")),
            protein_groups_file=data.get("protein_groups_file"),
            samples_file=data.get("samples_file"),
            basis=str(data.get("basis", "")),
            normalize=bool(data.get("normalize", False)),
            min_peptides=None if minimum is None else int(minimum),
            run_count=int(data.get("run_count", 0)),
            runs=[
                ObservationRun(str(r.get("file_name", "")), str(r.get("sample_id", "")),
                               int(r.get("fraction", 0)), int(r.get("technical_replicate", 0)))
                for r in data.get("runs") or []
            ],
            sample_count=int(data.get("sample_count", 0)),
            samples=list(data.get("samples") or []),
            samples_without_values=list(data.get("samples_without_values") or []),
            peptide_count=int(data.get("peptide_count", 0)),
            protein_group_count=int(data.get("protein_group_count", 0)),
            protein_groups=[
                ProteinGroupInfo(str(g.get("name", "")), g.get("q_value"), bool(g.get("is_decoy", False)),
                                 bool(g.get("is_contaminant", False)), bool(g.get("is_entrapment", False)),
                                 g.get("genes"))
                for g in data.get("protein_groups") or []
            ],
            state_counts=dict(data.get("state_counts") or {}),
            normalization=NormalizationResult._from_wire(data.get("normalization") or {}),
            warnings=list(data.get("warnings") or []),
            written=WrittenTable._from_wire(written) if written else None,
            responses_written=None if not responses else WrittenResponses(
                str(responses.get("path", "")), str(responses.get("format", "")),
                int(responses.get("row_count", 0)), int(responses.get("sample_count", 0))),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


# ---- calls -------------------------------------------------------------------------------------


def design(
    samples: _PathLike,
    *,
    factors: Sequence[FactorSpec | str] | None = None,
    covariates: Sequence[CovariateSpec] | None = None,
    interactions: Sequence[tuple[str, str]] | None = None,
    strata: Sequence[str] | str | None = None,
    batch: str | None = None,
    out: _PathLike | None = None,
    timeout: float | None = 60,
) -> Design:
    """Build each stratum's design matrix and default contrasts from a sample sheet, with mzLib.

    mzLib's ``AnalysisDesign.ByStratum`` splits the samples by stratum, and in each stratum builds
    one additive model: an intercept, each factor treatment-coded against its reference level,
    batch, each covariate as ``(value - centre) / scale``, and any interactions. The matrix equals
    R's ``model.matrix`` for the same formula. It then lists ``DefaultContrasts`` - each level
    against its reference, every pair when there is no reference, each covariate's slope - and
    judges whether each can be estimated at all from these samples.

    Args:
        samples: The sample sheet: a TSV with a header row and one row per biological sample. First
            column the sample id; the other columns are available by header. A column the model
            does not name is ignored (:attr:`Design.unused_columns`).
        factors: Condition factors, in column order, each a :class:`FactorSpec` or just a column
            name (no reference: every pair of levels is compared).
        covariates: Numeric covariates, in column order, each a :class:`CovariateSpec` with its
            stored unit.
        interactions: Two-factor interactions, as ``(first, second)`` pairs of factor names.
        strata: Columns whose values split the samples into strata, each designed on its own (one
            name, or a list).
        batch: The column holding each sample's batch, fitted as a fixed block. ``None``: no batch.
        out: A directory to write ``design_<i>.tsv`` (for :func:`pymzlib.stats.fit`) and
            ``contrasts_<i>.tsv`` into, per stratum. ``None``: nothing is written.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`Design`.

    Raises:
        UsageError: the sheet is missing or malformed; a named column is not in it; a column is
            given two roles; a sample lacks a factor or covariate value; or mzLib refuses the
            design (batch on some samples only, an interaction of an unknown factor).

    Examples:
        The MALAT1 dilution series :func:`pymzlib.stats.fit` analyses, described by each run's load
        rather than by a hand-made 0/1 design:

        >>> design = pymzlib.diff.design("malat_dilution_samples.tsv",
        ...                              factors=[pymzlib.diff.FactorSpec("malat_load", "500ng")],
        ...                              out="malat_design")
        >>> stratum = design.stratum("all")
        >>> stratum.sample_count, stratum.rank, stratum.fitted_column_names
        (27, 3, ['intercept', 'malat_load=125ng', 'malat_load=250ng'])
        >>> for row in design.records:
        ...     print(row["contrast_id"], row["label"], row["estimability"], row["fitted_weights"])
        c1 malat_load=125ng vs malat_load=500ng estimable [0, 1, 0]
        c2 malat_load=250ng vs malat_load=500ng estimable [0, 0, 1]

        Each contrast against the reference is one design column, so ``stats.fit`` tests it by
        that column's name in the design file this call wrote.
    """
    args = ["diff", "design", "--samples", _path(samples, "samples")]
    stdin = _model_stdin(factors, covariates, interactions, strata, batch)
    if out is not None:
        args += ["--out", _path(out, "out")]
    _bridge.require_verb("diff design", since=_SINCE)
    data = _bridge.invoke(*args, stdin=stdin, timeout=timeout)
    return Design._from_wire(data)


def observations(
    peptides: _PathLike,
    runs: Sequence[ObservationRun | tuple],
    *,
    basis: str,
    protein_groups: _PathLike | None = None,
    samples: _PathLike | None = None,
    normalize: bool = False,
    min_peptides: int | None = None,
    out: _PathLike | None = None,
    responses_out: _PathLike | None = None,
    timeout: float | None = 600,
) -> Observations:
    """Every peptide in every biological sample, as log2 intensity or no value with the reason, by mzLib.

    mzLib's ``FlashLfqObservations.FromPeptideTable`` reads a stored FlashLFQ or MetaMorpheus
    peptide table, maps each run's detection type to a state, sums each technical replicate's
    fractions on the linear scale and averages a sample's replicates in log2. Then
    ``SharedPeptideNormalization.Measure`` measures how far each sample sits from the others on
    the peptides they share; ``normalize=True`` subtracts it.

    Args:
        peptides: A FlashLFQ ``QuantifiedPeptides.tsv`` or a MetaMorpheus
            ``AllQuantifiedPeptides.tsv``.
        runs: Every run of the table, each an :class:`ObservationRun` or a
            ``(file, sample[, fraction[, technical_replicate]])`` tuple. ``file`` is the run's name
            as the ``Intensity_`` columns write it, without extension.
        basis: ``"msms_only"`` (no match-between-runs value enters) or ``"mbr_kept"``; see
            :data:`QUANT_BASES`. Required: mzLib has no default.
        protein_groups: A MetaMorpheus ``AllQuantifiedProteinGroups.tsv``, for each group's q-value
            and decoy, contaminant and entrapment flags. ``None``: none described.
        samples: A TSV whose first column (after a header row) lists the samples to keep, e.g. one
            stratum's design file from :func:`design`. ``None``: every sample.
        normalize: Subtract each sample's shift. The shifts are measured and reported either way.
        min_peptides: The floor on the normalization's reference set (peptides). ``None``: not
            sent, so mzLib's own floor, 100.
        out: A ``.tsv`` path to write the long table to instead of returning its rows.
        responses_out: A ``.tsv`` path to also write the peptide-by-sample log2 table that
            :func:`pymzlib.stats.fit` reads as ``responses``.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`Observations` table.

    Raises:
        UsageError: a file is missing or of the wrong type; ``basis`` is not one of
            :data:`QUANT_BASES`; the runs and the table's columns do not match one to one; or
            mzLib refuses the table (a duplicate peptide, a damaged intensity).

    Examples:
        Two runs of FlashLFQ's E. coli spike-in vignette, which their names say were spiked at 1x and
        at 2x:

        >>> table = pymzlib.diff.observations(
        ...     "EcoliSpikeIn_QuantifiedPeptides.tsv",
        ...     [("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"),
        ...      ("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x")],
        ...     basis="msms_only")
        >>> table.peptide_count, table.samples
        (243, ['spike_1x', 'spike_2x'])
        >>> table.state_counts
        {'quantified': 319, 'mbr_transferred': 43, 'ambiguous_peak': 9, 'identified_not_quantified': 6, 'not_detected': 109}
        >>> table.normalization.setting, table.normalization.reference_set_size
        ('none', 86)

        Only 86 peptides were measured in both runs, under mzLib's floor of 100, so mzLib measured
        no shift rather than trust too few.
    """
    if basis not in QUANT_BASES:
        raise _bridge.UsageError(
            f"basis must be one of {', '.join(QUANT_BASES)}; got {basis!r}. mzLib has no default: "
            "whether match-between-runs values enter is your choice."
        )
    if not isinstance(normalize, bool):
        raise _bridge.UsageError(f"normalize must be True or False; got {normalize!r}.")
    args = ["diff", "observations", "--peptides", _path(peptides, "peptides"), "--basis", basis]
    stdin = _runs_stdin(runs)
    if protein_groups is not None:
        args += ["--protein-groups", _path(protein_groups, "protein_groups")]
    if samples is not None:
        args += ["--samples", _path(samples, "samples")]
    if normalize:
        args.append("--normalize")
    if min_peptides is not None:
        args += ["--min-peptides", _integer(min_peptides, "min_peptides")]
    if out is not None:
        args += ["--out", _path(out, "out")]
    if responses_out is not None:
        args += ["--responses-out", _path(responses_out, "responses_out")]
    _bridge.require_verb("diff observations", since=_SINCE)
    data = _bridge.invoke(*args, stdin=stdin, timeout=timeout)
    return Observations._from_wire(data)
