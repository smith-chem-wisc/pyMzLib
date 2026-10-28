"""PTM site occupancy and PTM-QTL statistics, computed by mzLib's ``Quantification.PtmQtl``.

Three questions, each answered by mzLib (mzLib #1369, #1430, #1446, #1464, #1468):

==============================================================  ======================  ===============================================
You want to know                                                Call                    mzLib computes it with
==============================================================  ======================  ===============================================
What fraction of each site carries its modification, per run?  :func:`occupancy`       ``SiteOccupancyCalculator.Calculate``
Which modifications sit together, or rise and fall together?   :func:`pairs`           ``PtmPairEngine.Physical`` / ``.CoVarying``
Does a site's occupancy change with a per-run trait?           :func:`trait_effects`   ``SiteTraitEffects.Fit`` (a mixed model)
==============================================================  ======================  ===============================================

**Occupancy is a ratio of intensities.** For a site in a run, it is the summed intensity of the
quantified peptidoforms that carry the modification there, over the summed intensity of every
quantified peptidoform that covers the position - the estimator MetaMorpheus writes as
``IntensityOccupancy_``. mzLib holds every cell's **state** apart: a ``"floor"`` (the modified form
was identified but not quantified) is low but not zero, and ``"not_detected"`` (covered, and no
modified form identified) is a bound below what the run could detect. Neither has a ``fraction``,
so neither can be averaged in as 0 by accident::

    >>> occ = pymzlib.ptmqtl.occupancy("cases.tsv")
    >>> [(r["run"], r["state"], r["fraction"]) for r in occ.records
    ...  if r["protein_accession"] == "STATES"]
    [('countOnly', 'count_only', None), ('floor', 'floor', None), ('notDetected', 'not_detected', None), ('quant', 'quantified', 0.1)]

**Your data, as one table.** Every function reads the same plain table: one row per (run,
peptidoform, protein) identified by MS/MS, with its intensity in that run. That is mzLib's own
``PeptidoformObservation``, as columns. mzLib has no reader that builds it from a MetaMorpheus or
FlashLFQ result, so you assemble it - the guide shows how, from files pyMzLib already reads.

**One file between the steps.** :func:`occupancy` with ``out=`` writes the occupancy table that
:func:`trait_effects` reads. :func:`pairs` reads the observations themselves, because mzLib needs
the peptidoforms to find pairs on one molecule and to tell which pairs share a peptide.

Wire verbs: ``ptmqtl occupancy``, ``ptmqtl pairs`` and ``ptmqtl trait-effects``; their
language-neutral specs live in the bridge repository under ``design/verbs/``.
"""

from __future__ import annotations

import numbers
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Union

from . import _bridge
from .readers import WrittenTable, _Table

__all__ = [
    "OCCUPANCY_SCALES",
    "OCCUPANCY_STATES",
    "PAIR_TEST_STATUSES",
    "PAIR_TYPES",
    "RUN_LEVEL_GROUPS",
    "TRAIT_STATUSES",
    "PtmPairs",
    "SiteOccupancy",
    "TraitEffects",
    "occupancy",
    "pairs",
    "trait_effects",
]

#: The first pyMzLib whose bridge has the ptmqtl verbs (``since.pymzlib`` of their specs).
_SINCE = "0.6.0"

#: Every value of the occupancy ``state`` column, as mzLib's ``OccupancyState`` is written on the
#: wire. Only ``"quantified"`` cells have a ``fraction``.
OCCUPANCY_STATES = ("quantified", "floor", "count_only", "not_detected")

#: What :func:`pairs` takes as ``type=``. There is no default; you choose.
PAIR_TYPES = ("physical", "covarying", "both")

#: Every value of the :func:`pairs` ``test_status`` column (mzLib's ``PairTestStatus``), for
#: co-varying pairs; physical pairs have ``None``. Only ``"tested"`` pairs carry a statistic.
PAIR_TEST_STATUSES = ("tested", "not_estimable", "too_few_shared_runs", "both_sites_mostly_bounds")

#: What :func:`trait_effects` takes as ``scale=`` (mzLib's ``OccupancyScale``): the log odds of
#: occupancy in natural log (mzLib's default) or in log2, where an effect of +1 doubles the odds.
OCCUPANCY_SCALES = ("natural_log_odds", "log2_odds")

#: The named rules :func:`pairs` takes as ``run_level_group=``. For any other grouping, map each
#: modification to a group in a file and pass it as ``run_level_groups=``.
RUN_LEVEL_GROUPS = ("category", "modification")

#: Every value of the :func:`trait_effects` ``status`` column. Only ``"fitted"`` rows carry statistics.
TRAIT_STATUSES = (
    "fitted",
    "below_support",
    "too_few_observations",
    "rank_deficient",
    "too_few_groups",
    "variance_ratio_at_limit",
    "not_converged",
)

_PathLike = Union[str, "os.PathLike[str]"]


def __dir__() -> list[str]:
    """The public API only, so ``dir()`` is a map of what to call rather than of our imports."""
    return sorted(__all__)


def _columns(data: Mapping[str, Any]) -> tuple[list[str], dict[str, list[Any]] | None]:
    names = list(data.get("column_names") or [])
    columns = data.get("columns")
    return names, None if columns is None else dict(columns)


def _optional_list(value: Any) -> list[str] | None:
    return None if value is None else list(value)


def _optional_int(value: Any) -> int | None:
    return None if value is None else int(value)


# ---- result types ------------------------------------------------------------------------------


@dataclass(frozen=True)
class SiteOccupancy(_Table):
    """What :func:`occupancy` returns: one row per (site, run) where the site was identified or covered.

    Rows are in mzLib's order: by protein, then run, then position and modification.

    =========================  ============================  ===================================================
    column                     unit                          meaning; ``None`` means
    =========================  ============================  ===================================================
    ``site_key``               -                             ``accession:residue position:modification``, the
                                                             key :func:`pairs` and :func:`trait_effects` use
    ``protein_accession``      -                             the protein, as given
    ``position``               -                             1-based residue; a protein N-terminal modification
                                                             at the residue it sits on
    ``residue``                -                             one letter
    ``modification``           -                             ``Category:IdWithMotif``, as in the full sequence
    ``category``               -                             the text before the first ``:`` (``""`` if none)
    ``run``                    -                             the run, or with ``samples=`` the sample
    ``state``                  -                             one of :data:`OCCUPANCY_STATES`
    ``modified_intensity``     intensity (instrument units)  summed intensity of the quantified forms carrying
                                                             the modification here; 0 unless quantified
    ``covering_intensity``     intensity (instrument units)  summed intensity of every quantified form covering
                                                             the position; 0 for ``count_only``
    ``fraction``               fraction (0 to 1)             modified over covering intensity; ``None`` unless
                                                             ``state`` is ``"quantified"`` - censored, never 0
    ``unmodified_quantified``  -                             ``False`` with fraction 1 is a **ceiling**: the site
                                                             was seen only modified
    ``covering_peptidoforms``  peptidoforms                  distinct peptidoforms covering the position
                                                             identified in the run (quantified or not): the
                                                             n of k of n; ``None`` = not supplied
    ``modified_peptidoforms``  peptidoforms                  of those, the ones carrying this modification here:
                                                             the k; 0 for ``not_detected``
    =========================  ============================  ===================================================

    Attributes:
        observations_file: Absolute path of the observations table.
        samples_file: Absolute path of the run-to-sample table, or ``None`` without ``samples=``.
        observation_count: Observation rows read (rows).
        run_count: Distinct runs in the observations (runs).
        sample_count: Distinct samples (samples), or ``None`` without ``samples=``.
        run_names: The distinct values of the ``run`` column, in ordinal order: the runs, or the
            samples with ``samples=``.
        site_count: Distinct sites in the table (sites).
        state_counts: Cells per state (cells), keyed by :data:`OCCUPANCY_STATES`.
        modifications: The modification names the sites were limited to, or ``None`` for no filter.
        categories: The categories the sites were limited to, or ``None`` for no filter.
        ignored_columns: Columns of the observations file that are not fields of mzLib's
            observation, so were not read. A misspelled ``previous_residue`` shows up here.
        row_count: Cells in the table (rows), whether returned or written.
        column_names: The table's columns, in order; also the header ``out=`` writes.
        columns: Column name -> list of values, or ``None`` when the table was written to ``out=``.
        written: Where ``out=`` wrote the table, or ``None``.
        caveats: What the numbers are and are not, plus any that apply to this input.
    """

    observations_file: str
    samples_file: str | None
    observation_count: int
    run_count: int
    sample_count: int | None
    run_names: list[str]
    site_count: int
    state_counts: dict[str, int]
    modifications: list[str] | None
    categories: list[str] | None
    ignored_columns: list[str]
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] | None = None
    written: WrittenTable | None = None
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> SiteOccupancy:
        names, columns = _columns(data)
        written = data.get("written")
        return cls(
            observations_file=str(data.get("observations_file", "")),
            samples_file=data.get("samples_file"),
            observation_count=int(data.get("observation_count", 0)),
            run_count=int(data.get("run_count", 0)),
            sample_count=_optional_int(data.get("sample_count")),
            run_names=list(data.get("run_names") or []),
            site_count=int(data.get("site_count", 0)),
            state_counts=dict(data.get("state_counts") or {}),
            modifications=_optional_list(data.get("modifications")),
            categories=_optional_list(data.get("categories")),
            ignored_columns=list(data.get("ignored_columns") or []),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            written=WrittenTable._from_wire(written) if written else None,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class PtmPairs(_Table):
    """What :func:`pairs` returns: one row per pair of sites, physical pairs first.

    ==========================  ===================================  ==============================================
    column                      unit                                 meaning; ``None`` means
    ==========================  ===================================  ==============================================
    ``result_type``             -                                    ``"physical"`` (both on one peptidoform, so one
                                                                     molecule) or ``"covarying"`` (across runs)
    ``site_a``, ``site_b``      -                                    the two sites' keys, ``site_a`` first by
                                                                     ordinal order; a pair is written once
    ``same_protein``            -                                    both on one protein accession
    ``overlapping``             -                                    one peptide covers both; such a co-varying pair
                                                                     is correlated by construction
    ``statistic``               physical: fraction (0 to 1);         physical: median co-occupancy; covarying:
                                covarying: unitless (-1 to 1)        Spearman's rho. ``None`` = physical never
                                                                     quantified together, or not estimable
    ``p_value``                 fraction (0 to 1)                    two-sided Spearman; ``None`` = physical (no
                                                                     test) or not estimable
    ``q``                       fraction (0 to 1)                    Benjamini-Hochberg within ``fdr_family``;
                                                                     ``None`` = not tested (physical, overlapping,
                                                                     not estimable). Not a target-decoy q-value
    ``fdr_family``              -                                    mzLib's family: ``P`` or ``A``, then ``intra``
                                                                     or ``inter`` protein, e.g. ``"A:inter"``
    ``n``                       runs                                 physical: runs identified together; covarying:
                                                                     runs where both were quantified
    ``statistic_n``             runs                                 runs ``statistic`` was computed on
    ``spearman_method``         -                                    ``"exact"``, ``"asymptotic"`` or
                                                                     ``"not_estimable"``; ``None`` for physical
    ``run_level_removed``       -                                    whether rho is partial given the run level
    ``test_status``             -                                    one of :data:`PAIR_TEST_STATUSES`: tested, or
                                                                     why not; ``None`` for physical
    ``bounds_included``         -                                    whether one-form cells entered as tied bounds
    ``run_covariate_removed``   -                                    whether rho is partial given the run covariate
    ==========================  ===================================  ==============================================

    Attributes:
        observations_file: Absolute path of the observations table.
        samples_file: Absolute path of the run-to-sample table, or ``None`` without ``samples=``.
        type: ``"physical"``, ``"covarying"`` or ``"both"``, as you asked.
        observation_count: Observation rows read (rows).
        run_count: Distinct runs in the observations (runs).
        sample_count: Distinct samples (samples), or ``None`` without ``samples=``.
        run_names: The runs, or samples, the pairs were computed over, in ordinal order.
        covarying_run_count: The runs (or samples) where any site was identified or covered (runs):
            the denominator of ``min_quantified_fraction``. ``None`` for ``type="physical"``.
        min_quantified_fraction: The guard used (fraction of the runs, 0 to 1): yours, or mzLib's
            default. ``None`` for ``type="physical"``.
        exclude_ceiling: ``False`` with ``keep_ceiling=True``. ``None`` when you did not choose,
            so mzLib's own default applied (ceilings left out, or with ``include_bounds`` a high
            bound; ``PtmPairEngine.CoVarying``), and for ``type="physical"``.
        include_bounds: Whether one-form cells entered as tied bounds; ``None`` for
            ``type="physical"``.
        min_untied_share: mzLib's ``MinUntiedShare`` (fraction of the shared runs, 0 to 1): with
            bounds, a pair is tested only if one of its sites has at least this share of the shared
            runs outside its largest tie. ``None`` without bounds.
        run_covariate_file: Absolute path of the run covariate table, or ``None``.
        run_covariate_name: The covariate's name, the table's second header, or ``None``.
        run_covariate_runs_missing: Runs (or samples) with no covariate value, left out of every
            pair; empty when every run has one, ``None`` without ``run_covariate=``.
        remove_run_level: Whether the run level was removed; ``None`` for ``type="physical"``.
        run_level_group: ``"none"`` (one group), ``"category"``, ``"modification"`` or ``"file"``;
            ``None`` when the run level was not removed.
        run_level_groups_file: Absolute path of the modification-to-group table, or ``None``.
        min_runs_for_run_level: mzLib's ``MinRunsForRunLevel`` (runs): with the run level or a run
            covariate removed, a pair sharing fewer runs is ``"too_few_shared_runs"``. ``None`` when
            neither was removed.
        min_sites_for_run_level: mzLib's ``MinSitesForRunLevel`` (sites): a run whose group has fewer
            other sites there has no level. ``None`` when not removed.
        modifications: The modification names the sites were limited to, or ``None``.
        categories: The categories the sites were limited to, or ``None``.
        ignored_columns: Columns of the observations file that were not read.
        physical_count: Physical pairs found (pairs), or ``None`` for ``type="covarying"``.
        covarying_count: Co-varying pairs reported (pairs), or ``None`` for ``type="physical"``.
        row_count: Pairs in the table (rows).
        column_names: The table's columns, in order.
        columns: Column name -> one value per pair.
        caveats: What each pair type means, and the run level when it was removed, plus any that
            apply to this input.
    """

    observations_file: str
    samples_file: str | None
    type: str
    observation_count: int
    run_count: int
    sample_count: int | None
    run_names: list[str]
    covarying_run_count: int | None
    min_quantified_fraction: float | None
    exclude_ceiling: bool | None
    include_bounds: bool | None
    min_untied_share: float | None
    run_covariate_file: str | None
    run_covariate_name: str | None
    run_covariate_runs_missing: list[str] | None
    remove_run_level: bool | None
    run_level_group: str | None
    run_level_groups_file: str | None
    min_runs_for_run_level: int | None
    min_sites_for_run_level: int | None
    modifications: list[str] | None
    categories: list[str] | None
    ignored_columns: list[str]
    physical_count: int | None
    covarying_count: int | None
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> PtmPairs:
        names, columns = _columns(data)
        return cls(
            observations_file=str(data.get("observations_file", "")),
            samples_file=data.get("samples_file"),
            type=str(data.get("type", "")),
            observation_count=int(data.get("observation_count", 0)),
            run_count=int(data.get("run_count", 0)),
            sample_count=_optional_int(data.get("sample_count")),
            run_names=list(data.get("run_names") or []),
            covarying_run_count=_optional_int(data.get("covarying_run_count")),
            min_quantified_fraction=data.get("min_quantified_fraction"),
            exclude_ceiling=data.get("exclude_ceiling"),
            include_bounds=data.get("include_bounds"),
            min_untied_share=data.get("min_untied_share"),
            run_covariate_file=data.get("run_covariate_file"),
            run_covariate_name=data.get("run_covariate_name"),
            run_covariate_runs_missing=_optional_list(data.get("run_covariate_runs_missing")),
            remove_run_level=data.get("remove_run_level"),
            run_level_group=data.get("run_level_group"),
            run_level_groups_file=data.get("run_level_groups_file"),
            min_runs_for_run_level=_optional_int(data.get("min_runs_for_run_level")),
            min_sites_for_run_level=_optional_int(data.get("min_sites_for_run_level")),
            modifications=_optional_list(data.get("modifications")),
            categories=_optional_list(data.get("categories")),
            ignored_columns=list(data.get("ignored_columns") or []),
            physical_count=_optional_int(data.get("physical_count")),
            covarying_count=_optional_int(data.get("covarying_count")),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns or {},
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class TraitEffects(_Table):
    """What :func:`trait_effects` returns: one row per site with a usable response, by site key.

    ======================  ===============================  ==================================================
    column                  unit                             meaning; ``None`` means
    ======================  ===============================  ==================================================
    ``site_key``            -                                as :func:`occupancy` writes it
    ``protein_accession``,  -                                the site, as :func:`occupancy`
    ``position``,
    ``residue``,
    ``modification``
    ``status``              -                                one of :data:`TRAIT_STATUSES`
    ``effect``              log odds of occupancy per trait  change in the log odds of occupancy per unit of
                            unit, natural log or log2 as     the trait; ``exp(effect)`` (natural) or
                            ``scale`` says                   ``2 ** effect`` (log2) is the odds ratio.
                                                             ``None`` = not fitted
    ``standard_error``      as ``effect``                    its standard error
    ``df``                  degrees of freedom               containment rule: between-replicate for a trait
                                                             constant within each replicate
    ``p_value``             fraction (0 to 1)                two-sided, for the trait
    ``q``                   fraction (0 to 1)                Benjamini-Hochberg over this call's fitted sites
    ``runs``                runs                             runs contributing a response
    ``replicates``          replicates                       distinct replicates among them
    ``median_fraction``     fraction (0 to 1)                median occupancy over those runs
    ``covariate_effects``   log odds per covariate unit      one per covariate, in ``covariate_names`` order;
                            (as ``scale``)                   each ``None`` when not fitted
    ``replicate_variance``  squared log odds (as ``scale``)  between-replicate variance (random intercept)
    ``residual_variance``   squared log odds (as ``scale``)  within-replicate, run-to-run variance
    ======================  ===============================  ==================================================

    Attributes:
        occupancy_file: Absolute path of the occupancy table.
        runs_file: Absolute path of the runs table.
        trait_name: The runs table's third header: the trait.
        covariate_names: The covariate headers, in order; empty without covariates.
        cell_count: Occupancy rows read (cells).
        run_count: Runs in the runs table (runs).
        trait_run_count: Runs with a trait value (runs): the runs the design is built over.
        occupancy_runs_not_used: Occupancy runs with no trait value, so not used.
        min_runs: The support threshold used (runs): yours, or mzLib's default.
        min_replicates: The replicate threshold used (replicates): yours, or mzLib's default.
        logit_epsilon: The clamp used before the logit (fraction, 0 to 1).
        exclude_ceiling: Whether ceiling cells were left out.
        scale: The scale of every effect, standard error and variance: ``"natural_log_odds"``
            (mzLib's default) or ``"log2_odds"``, one of :data:`OCCUPANCY_SCALES`.
        status_counts: Sites per status (sites), keyed by :data:`TRAIT_STATUSES`.
        ignored_columns: Columns of the occupancy file that were not read.
        row_count: Sites in the table (rows).
        column_names: The table's columns, in order.
        columns: Column name -> one value per site.
        caveats: The model and what it leaves out, plus any that apply to this input.
    """

    occupancy_file: str
    runs_file: str
    trait_name: str
    covariate_names: list[str]
    cell_count: int
    run_count: int
    trait_run_count: int
    occupancy_runs_not_used: list[str]
    min_runs: int
    min_replicates: int
    logit_epsilon: float
    exclude_ceiling: bool
    scale: str
    status_counts: dict[str, int]
    ignored_columns: list[str]
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> TraitEffects:
        names, columns = _columns(data)
        return cls(
            occupancy_file=str(data.get("occupancy_file", "")),
            runs_file=str(data.get("runs_file", "")),
            trait_name=str(data.get("trait_name", "")),
            covariate_names=list(data.get("covariate_names") or []),
            cell_count=int(data.get("cell_count", 0)),
            run_count=int(data.get("run_count", 0)),
            trait_run_count=int(data.get("trait_run_count", 0)),
            occupancy_runs_not_used=list(data.get("occupancy_runs_not_used") or []),
            min_runs=int(data.get("min_runs", 0)),
            min_replicates=int(data.get("min_replicates", 0)),
            logit_epsilon=float(data.get("logit_epsilon", 0.0)),
            exclude_ceiling=bool(data.get("exclude_ceiling", True)),
            scale=str(data.get("scale", "")),
            status_counts=dict(data.get("status_counts") or {}),
            ignored_columns=list(data.get("ignored_columns") or []),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns or {},
            caveats=list(data.get("caveats") or []),
        )


# ---- calls -------------------------------------------------------------------------------------


def _path(value: object, what: str) -> str:
    text = _bridge.path_text(value)
    if not text:
        raise _bridge.UsageError(f"{what} must be a path (str or os.PathLike); got {value!r}.")
    return text


def _names(value: object, what: str) -> str:
    """A list of names as the wire's one tab-separated option."""
    if isinstance(value, str):
        raise _bridge.UsageError(f"{what} takes a list of names; for one, use [{value!r}].")
    if not isinstance(value, Sequence) or not value:
        raise _bridge.UsageError(f"{what} takes a non-empty list of names; got {value!r}.")
    for name in value:
        if not isinstance(name, str) or not name.strip():
            raise _bridge.UsageError(f"Every name in {what} must be a non-empty str; got {name!r}.")
        if any(c in name for c in "\t\r\n"):
            raise _bridge.UsageError(f"A name in {what} contains a tab or a line break: {name!r}.")
    return "\t".join(value)


def _flag(value: object, what: str) -> bool:
    if not isinstance(value, bool):
        raise _bridge.UsageError(f"{what} must be True or False; got {value!r}.")
    return value


def _filter_args(modifications: Sequence[str] | None, categories: Sequence[str] | None) -> list[str]:
    if modifications is not None and categories is not None:
        raise _bridge.UsageError(
            "Give modifications or categories, not both: mzLib takes one filter on the modification name."
        )
    if modifications is not None:
        return ["--modifications", _names(modifications, "modifications")]
    if categories is not None:
        return ["--categories", _names(categories, "categories")]
    return []


def occupancy(
    observations: _PathLike,
    *,
    samples: _PathLike | None = None,
    modifications: Sequence[str] | None = None,
    categories: Sequence[str] | None = None,
    out: _PathLike | None = None,
    timeout: float | None = 600,
) -> SiteOccupancy:
    """Compute every site's occupancy in every run, with mzLib's ``SiteOccupancyCalculator``.

    For a site (protein, position, modification) in a run, occupancy is the summed intensity of the
    quantified peptidoforms that carry the modification there, over the summed intensity of every
    quantified peptidoform that covers the position. Each cell keeps its state, so a censored cell
    is never mistaken for a 0: see :data:`OCCUPANCY_STATES`.

    mzLib always leaves out the categories ``Common Variable`` and ``Common Fixed``,
    peptide-terminal modifications and every C-terminal modification. An N-terminal modification is
    a site only on a peptide that starts at residue 1, or at 2 after the initiator Met is removed.

    Args:
        observations: A TSV with a header row naming its columns, one row per (run, peptidoform,
            protein) identified by MS/MS: ``run``, ``full_sequence`` (MetaMorpheus notation),
            ``protein_accession``, ``start_residue`` and ``end_residue`` (1-based, in the protein),
            ``intensity`` (instrument units; blank or ``NA`` for identified but not quantified), and
            optionally ``previous_residue`` (one letter). A shared peptide is one row per protein.
            Never one row per PSM, and never a match-between-runs transfer.
        samples: A TSV with a header row, then ``run``, ``sample``, to sum each sample's fractions
            before the ratio; the ``run`` column of the result then holds the sample. Map only
            fractions together, never technical replicates. Leave it out (``None``) for per-run
            occupancy.
        modifications: Count only these modifications as sites: exact names,
            ``Category:IdWithMotif``, as a list. Not with ``categories``. ``None``: no filter.
        categories: Count only modifications of these categories (the text before the first
            ``:``), as a list, e.g. ``["Common Biological"]``. Not with ``modifications``.
        out: Write the table to this ``.tsv`` instead of returning it - the file
            :func:`trait_effects` reads. Then :attr:`SiteOccupancy.columns` is ``None``.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`SiteOccupancy` table.

    Raises:
        UsageError: a file is missing or malformed; a required column is absent; ``out`` is not a
            ``.tsv`` or is the observations file; both filters were given; or mzLib refuses the
            observations - a repeated (run, peptidoform, protein), a span that does not fit the
            sequence, a negative intensity, a run with no sample.

    Examples:
        mzLib's own hand-computed cases, one protein each (``tests/fixtures/ptmqtl/README.md``):

        >>> occ = pymzlib.ptmqtl.occupancy("cases.tsv")
        >>> occ.site_count, occ.state_counts
        (9, {'quantified': 9, 'floor': 1, 'count_only': 3, 'not_detected': 1})
        >>> for r in occ.records:
        ...     if r["protein_accession"] in ("RATIO", "CEILING", "SHARED_POSITION"):
        ...         print(r["protein_accession"], r["modification"], r["fraction"],
        ...               r["unmodified_quantified"])
        CEILING Common Biological:Phosphorylation on S 1 False
        RATIO Common Biological:Phosphorylation on S 0.15 True
        SHARED_POSITION Common Biological:HexNAc on S 0.3 True
        SHARED_POSITION Common Biological:Phosphorylation on S 0.2 True

        ``RATIO`` is 30 / (30 + 70 + 100): a missed-cleavage form covering the site counts in the
        denominator. ``CEILING`` was seen only modified, so its 1 is a ceiling. The two
        modifications at one ``SHARED_POSITION`` share one denominator.

        A sample measured as two fractions is summed before the ratio, not averaged after:

        >>> occ = pymzlib.ptmqtl.occupancy("fractions.tsv", samples="fractions_samples.tsv")
        >>> [(r["protein_accession"], r["run"], r["state"], r["fraction"]) for r in occ.records]
        [('FLOOR_ACROSS', 'A', 'floor', None), ('SUMMED', 'A', 'quantified', 0.2), ('SUMMED', 'B', 'not_detected', None)]
    """
    args = ["ptmqtl", "occupancy", "--observations", _path(observations, "observations")]
    if samples is not None:
        args += ["--samples", _path(samples, "samples")]
    args += _filter_args(modifications, categories)
    if out is not None:
        args += ["--out", _path(out, "out")]
    _bridge.require_verb("ptmqtl occupancy", since=_SINCE)
    data = _bridge.invoke(*args, timeout=timeout)
    return SiteOccupancy._from_wire(data)


def pairs(
    observations: _PathLike,
    *,
    type: str,
    samples: _PathLike | None = None,
    modifications: Sequence[str] | None = None,
    categories: Sequence[str] | None = None,
    min_quantified_fraction: float | None = None,
    keep_ceiling: bool = False,
    include_bounds: bool = False,
    run_covariate: _PathLike | None = None,
    remove_run_level: bool = False,
    run_level_group: str | None = None,
    run_level_groups: _PathLike | None = None,
    timeout: float | None = 600,
) -> PtmPairs:
    """Find PTM-PTM pairs within one dataset, with mzLib's ``PtmPairEngine``.

    Two kinds of evidence, which you choose between with ``type``:

    * **physical** (``PtmPairEngine.Physical``): both modifications were identified on one
      peptidoform, so on one molecule. ``statistic`` is the median co-occupancy over the runs where
      both were quantified. There is no test.
    * **covarying** (``PtmPairEngine.CoVarying``): the two sites' occupancies rise and fall
      together (or oppositely) across runs. ``statistic`` is Spearman's rho on the runs where both
      are quantified, with a Benjamini-Hochberg ``q`` per family. Floors, undetected cells and, by
      default, ceilings are left out: two unrelated sites both read 1 in low-load runs, which would
      correlate them through detection alone. ``include_bounds=True`` lets them in as tied bounds
      instead, and ``run_covariate=`` correlates given each run's depth.

    Occupancy is computed inside by the same mzLib call :func:`occupancy` makes.

    Args:
        observations: :func:`occupancy`'s observations table, with the same rules.
        type: ``"physical"``, ``"covarying"`` or ``"both"`` (physical rows first), one of
            :data:`PAIR_TYPES`. Required: mzLib has no default between its two methods.
        samples: :func:`occupancy`'s run-to-sample table; pairs are then across samples.
        modifications: As :func:`occupancy`: only these modifications are sites.
        categories: As :func:`occupancy`: only these categories are sites.
        min_quantified_fraction: Co-varying only. A site enters pairs only if quantified in at least
            this fraction of the runs (fraction of the runs, above 0 and at most 1), rounded up.
            ``None``: mzLib's default, which :attr:`PtmPairs.min_quantified_fraction` reports.
        keep_ceiling: Co-varying only. Keep runs where a site was seen only modified. ``False``
            sends nothing, so mzLib's default (leave them out) applies.
        include_bounds: Co-varying only. Let a run where a site was seen in only one form enter as a
            tied bound, never an invented value: ``"not_detected"`` and ``"floor"`` cells rank below
            every fraction of the site, ceilings and count-only cells seen only modified above. A
            pair whose two sites are both mostly one bound is not tested
            (``"both_sites_mostly_bounds"``). ``False`` sends nothing: mzLib's default, no bounds.
        run_covariate: Co-varying only. A TSV with a header row whose second field names the
            covariate (e.g. ``target_psms``), then ``run``, value: each run's depth, in any unit.
            Pairs become the partial Spearman correlation given it; it is ranked, so log2 gives
            the same answer. Leave a run out for "unknown"; it is then left out of every pair.
        remove_run_level: Co-varying only. Correlate given each run's level - the median occupancy
            of the other tested sites of the pair's group - to remove a batch or loading shift.
            **Batch only**: where the level follows the design, it removes the effect you study.
        run_level_group: With ``remove_run_level``: ``"category"`` or ``"modification"``, one of
            :data:`RUN_LEVEL_GROUPS`, for one level per group of sites. ``None``: one level for every
            site, mzLib's default. Not with ``run_level_groups``.
        run_level_groups: With ``remove_run_level``: a TSV with a header row, then
            ``modification``, ``group``, for any grouping, such as chemistry. Every modification
            among the sites needs a group.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`PtmPairs` table.

    Raises:
        UsageError: ``type`` is not one of :data:`PAIR_TYPES`; a co-varying option was given with
            ``type="physical"``; a run-level grouping without ``remove_run_level``, or both
            groupings; ``min_quantified_fraction`` out of range; a file missing or malformed; or
            mzLib refuses the observations, as :func:`occupancy` lists, or a run covariate value
            that is not finite.

    Examples:
        Two modifications on one molecule, mzLib's ``PhysicalPairsCarryCoOccupancy`` case:

        >>> found = pymzlib.ptmqtl.pairs("cases.tsv", type="physical")
        >>> r = found.records[0]
        >>> r["site_a"]
        'TWO_ON_ONE:K14:Common Biological:Acetylation on K'
        >>> r["statistic"], r["n"], r["statistic_n"], r["p_value"]
        (0.25, 2, 1, None)

        Identified together in two runs, quantified in one, where 25 of the 100 intensity units
        covering both positions carried both marks.

        Twelve runs in two batches, every site 0.2 higher in the second (mzLib's
        ``RunLevelRemovalTakesOutABatchShift``). The batch alone correlates two sites:

        >>> plain = pymzlib.ptmqtl.pairs("batch.tsv", type="covarying")
        >>> given = pymzlib.ptmqtl.pairs("batch.tsv", type="covarying", remove_run_level=True)
        >>> def pair(result):
        ...     return next(r for r in result.records
        ...                 if r["site_a"].startswith("PAIR_S") and r["site_b"].startswith("PAIR_T"))
        >>> print(f"{pair(plain)['statistic']:.2f} {pair(given)['statistic']:.2f}")
        0.76 0.03

        Given each run's level, the pair is unrelated.

        Two sites with a value in six of ten runs, and a bound in the rest (mzLib's
        ``BoundsEnterTypeAAsTiedRanks``). Six is below mzLib's 70% guard, so no pair is tested
        until the bounds enter:

        >>> pymzlib.ptmqtl.pairs("bounds.tsv", type="covarying").row_count
        0
        >>> bounded = pymzlib.ptmqtl.pairs("bounds.tsv", type="covarying", include_bounds=True)
        >>> r = bounded.records[0]
        >>> r["n"], r["test_status"], r["bounds_included"]
        (9, 'tested', True)
    """
    if not isinstance(type, str) or type not in PAIR_TYPES:
        raise _bridge.UsageError(f"type must be one of {', '.join(PAIR_TYPES)}; got {type!r}.")
    args = ["ptmqtl", "pairs", "--observations", _path(observations, "observations"), "--type", type]
    if samples is not None:
        args += ["--samples", _path(samples, "samples")]
    args += _filter_args(modifications, categories)
    if min_quantified_fraction is not None:
        if isinstance(min_quantified_fraction, bool) or not isinstance(min_quantified_fraction, numbers.Real):
            raise _bridge.UsageError(
                f"min_quantified_fraction must be a number; got {min_quantified_fraction!r}."
            )
        args += ["--min-quantified-fraction", repr(float(min_quantified_fraction))]
    if _flag(keep_ceiling, "keep_ceiling"):
        args.append("--keep-ceiling")
    if _flag(include_bounds, "include_bounds"):
        args.append("--include-bounds")
    if run_covariate is not None:
        args += ["--run-covariate", _path(run_covariate, "run_covariate")]
    if _flag(remove_run_level, "remove_run_level"):
        args.append("--remove-run-level")
    if run_level_group is not None:
        if not isinstance(run_level_group, str) or run_level_group not in RUN_LEVEL_GROUPS:
            raise _bridge.UsageError(
                f"run_level_group must be one of {', '.join(RUN_LEVEL_GROUPS)}; got {run_level_group!r}. "
                "For any other grouping, pass a modification-to-group file as run_level_groups."
            )
        args += ["--run-level-group", run_level_group]
    if run_level_groups is not None:
        args += ["--run-level-groups", _path(run_level_groups, "run_level_groups")]
    _bridge.require_verb("ptmqtl pairs", since=_SINCE)
    data = _bridge.invoke(*args, timeout=timeout)
    return PtmPairs._from_wire(data)


def _count(value: object, what: str) -> str:
    if isinstance(value, bool) or not isinstance(value, int):
        raise _bridge.UsageError(f"{what} must be an int; got {value!r}.")
    return str(value)


def trait_effects(
    occupancy: _PathLike,
    runs: _PathLike,
    *,
    min_runs: int | None = None,
    min_replicates: int | None = None,
    logit_epsilon: float | None = None,
    keep_ceiling: bool = False,
    scale: str | None = None,
    timeout: float | None = 600,
) -> TraitEffects:
    """Fit, per site, the effect of a per-run trait on the log odds of occupancy, with mzLib's ``SiteTraitEffects``.

    The model, per site, is logit(occupancy) = b0 + b_trait x trait + covariates + a random
    intercept per biological replicate + error, fitted by REML (mzLib's ``MixedModel``). The
    random intercept is the point: repeated runs of one animal or donor are not counted as
    independent, so a trait that is constant within each replicate is tested at the
    between-replicate degrees of freedom.

    Only quantified cells are responses. Floors, count-only and undetected cells are left out,
    never entered as 0, so this answers "does the modified share change?", not "is the
    modification detected more often?". Ceilings are left out by default.

    Args:
        occupancy: The table :func:`occupancy` writes with ``out=`` (columns found by header name).
            A ``fraction`` cell is used as written, as mzLib does for a stored fraction.
        runs: A TSV with a header row, then one row per run: the run, its biological replicate,
            the trait (its header names it; any unit), then any covariates (each header names one).
            A run with a blank trait is not used.
        min_runs: Fewest runs with a response for a site to be fitted (runs). ``None``: mzLib's
            default, which :attr:`TraitEffects.min_runs` reports.
        min_replicates: Fewest distinct replicates among them (replicates). ``None``: mzLib's
            default.
        logit_epsilon: Fractions are clamped to [epsilon, 1 - epsilon] before the logit (fraction,
            above 0 and below 0.5). ``None``: mzLib's default.
        keep_ceiling: Use runs where the site was seen only modified. ``False`` sends nothing, so
            mzLib's default (leave them out) applies.
        scale: The unit of the effects and variances, one of :data:`OCCUPANCY_SCALES`.
            ``"natural_log_odds"``: ``exp(effect)`` is the odds ratio per trait unit.
            ``"log2_odds"``: an effect of +1 doubles the odds; mzLib divides the natural-log fit by
            ln 2 (variances by (ln 2)^2), and df, p and q do not change. ``None``: mzLib's default,
            which :attr:`TraitEffects.scale` reports.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`TraitEffects` table.

    Raises:
        UsageError: a file is missing or malformed, or repeats a (site, run); an option is out of
            range; or mzLib refuses the input - a run with a trait and no replicate, a covariate
            that is not finite, or a design that is not of full rank (a trait that never varies).

    Examples:
        mzLib's ``TraitDesign``: six replicates at three doses, two injections each, with
        logit(occupancy) = -2 + 0.5 x dose + a replicate offset. Occupancy is computed and written
        first; ``trait_runs.tsv`` gives each run's replicate, dose and an injection covariate:

        >>> occ = pymzlib.ptmqtl.occupancy("trait.tsv", out="trait_occupancy.tsv")
        >>> occ.written.row_count, occ.columns is None
        (24, True)
        >>> fit = pymzlib.ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv")
        >>> fit.trait_name, fit.covariate_names, fit.min_runs, fit.min_replicates
        ('dose', ['injection'], 6, 4)
        >>> for r in fit.records:
        ...     print(r["protein_accession"], r["status"], f"{r['effect']:.2f}", r["df"],
        ...           r["runs"], r["replicates"])
        CENSORED_SITE fitted 0.54 3 10 5
        TRAIT_SITE fitted 0.50 4 12 6

        The effect, 0.5 logit units per dose unit, comes back. ``CENSORED_SITE`` lost one
        replicate's two runs - a floor and a ceiling - so it rests on five replicates.

        The same fit in log2 odds, where +1 doubles the odds of the modified form:

        >>> log2 = pymzlib.ptmqtl.trait_effects("trait_occupancy.tsv", "trait_runs.tsv",
        ...                                     scale="log2_odds")
        >>> trait = next(r for r in log2.records if r["protein_accession"] == "TRAIT_SITE")
        >>> log2.scale, f"{trait['effect']:.2f}", f"{trait['p_value']:.5f}"
        ('log2_odds', '0.72', '0.00082')
    """
    args = ["ptmqtl", "trait-effects", "--occupancy", _path(occupancy, "occupancy"),
            "--runs", _path(runs, "runs")]
    if min_runs is not None:
        args += ["--min-runs", _count(min_runs, "min_runs")]
    if min_replicates is not None:
        args += ["--min-replicates", _count(min_replicates, "min_replicates")]
    if logit_epsilon is not None:
        if isinstance(logit_epsilon, bool) or not isinstance(logit_epsilon, numbers.Real):
            raise _bridge.UsageError(f"logit_epsilon must be a number; got {logit_epsilon!r}.")
        args += ["--logit-epsilon", repr(float(logit_epsilon))]
    if _flag(keep_ceiling, "keep_ceiling"):
        args.append("--keep-ceiling")
    if scale is not None:
        if not isinstance(scale, str) or scale not in OCCUPANCY_SCALES:
            raise _bridge.UsageError(f"scale must be one of {', '.join(OCCUPANCY_SCALES)}; got {scale!r}.")
        args += ["--scale", scale]
    _bridge.require_verb("ptmqtl trait-effects", since=_SINCE)
    data = _bridge.invoke(*args, timeout=timeout)
    return TraitEffects._from_wire(data)
