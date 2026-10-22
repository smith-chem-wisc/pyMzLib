"""Differential abundance with no R: limma's moderated t-test, Benjamini-Hochberg, random-effects
meta-analysis, p-value combination and Spearman correlation, computed by mzLib.

Five questions, each answered by mzLib's ``StatisticalModels`` (mzLib #1341, #1357, #1367, #1446):

==========================================================  =====================  =============================================
You want to know                                            Call                   mzLib computes it with
==========================================================  =====================  =============================================
Which features change with a coefficient, and how surely?   :func:`fit`            ``LinearModel.Fit`` + ``EmpiricalBayes.Moderate``
Which of these p-values survive a false discovery rate?     :func:`adjust`         ``MultipleTesting.BenjaminiHochberg``
What is one feature's effect, pooled across studies?        :func:`meta`           ``RandomEffectsMeta.Pool`` (DerSimonian-Laird)
What do several studies' p-values say together?             :func:`combine_p`      ``PValueCombination.Fisher`` / ``.Stouffer``
Which features rise or fall with a per-sample trait?        :func:`correlate`      ``SpearmanCorrelation.Correlate`` / ``.PartialCorrelate``
==========================================================  =====================  =============================================

**It is limma, checked.** :func:`fit` is limma's ``lmFit`` followed by ``eBayes(legacy = TRUE)``
(Smyth 2004), and mzLib holds it to limma 3.68.5 to 1e-8 relative on limma's own reference output.
pyMzLib ships that reference - limma's input, and limma's own ``eBayes`` output for it - so you
can run the comparison yourself::

    >>> fit = pymzlib.stats.fit("limma_reference_responses.tsv", "limma_reference_design.tsv",
    ...                         ["age_decades"])
    >>> import csv
    >>> with open("limma_ebayes_notrend.tsv", newline="") as fh:
    ...     limma = list(csv.DictReader(fh, delimiter="\\t"))
    >>> def worst(ours, theirs):
    ...     return max(abs(a - float(b)) / abs(float(b)) for a, b in zip(ours, theirs))
    >>> worst(fit.columns["t"], [row["t_age"] for row in limma]) < 1e-8
    True
    >>> worst(fit.columns["bh_adjusted"], [row["bh_age"] for row in limma]) < 1e-8
    True
    >>> fit.feature_count, fit.residual_df_differ
    (400, True)

Every one of the 400 moderated t-statistics and adjusted p-values agrees with limma to better than
one part in 10^8.

**What it is not.** It is limma's *legacy* estimator. When features have different residual
degrees of freedom - which omitting missing values causes, so in most label-free data - limma 3.61
and later default to a different prior estimator, and :attr:`ModeratedFit.residual_df_differ`
says when your input is such a case. ``robust = TRUE``, contrasts, the B-statistic and observation
weights are not implemented; write a contrast as its own design column instead.

**Your data, as files.** :func:`fit` reads two tab-separated files - a feature-by-sample table of
(already log-transformed) responses, and a sample-by-coefficient design - because that is the
shape a table library writes (``df.to_csv(path, sep="\\t")``) and it keeps a 10,000-protein table
off the command line. A blank, ``NA`` or ``NaN`` cell is missing and is left out of that feature's
fit, never imputed. **A 0 is an observation**, as in mzLib, so blank out any 0 that means "not
measured".

Wire verbs: ``stats fit``, ``stats adjust``, ``stats meta``, ``stats combine-p`` and
``stats correlate``; their language-neutral specs live in the bridge repository under
``design/verbs/``.
"""

from __future__ import annotations

import math
import numbers
import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Union

from . import _bridge
from .readers import _Table

__all__ = [
    "COMBINATION_METHODS",
    "FIT_STATUSES",
    "SPEARMAN_METHODS",
    "Adjusted",
    "CombinedPValues",
    "Correlations",
    "MetaAnalysis",
    "ModeratedFit",
    "VariancePrior",
    "adjust",
    "combine_p",
    "correlate",
    "fit",
    "meta",
]

#: The first pyMzLib whose bridge has the stats verbs (``since.pymzlib`` of their specs).
_SINCE = "0.3.0"

#: The first pyMzLib whose bridge has ``stats combine-p`` and ``stats correlate`` (mzLib 1.0.595).
_COMBINE_CORRELATE_SINCE = "0.6.0"

#: Every value of the ``status`` column, as mzLib's ``FeatureFitStatus`` is written on the wire.
#: Only ``"fitted"`` rows carry statistics or enter the Benjamini-Hochberg family.
FIT_STATUSES = ("fitted", "too_few_observations", "rank_deficient")

#: The rules :func:`combine_p` takes as ``method=``: mzLib's ``PValueCombinationMethod``. There is
#: no default; you choose.
COMBINATION_METHODS = ("fisher", "stouffer")

#: Every value of the ``method`` column of :func:`correlate`, as mzLib's ``SpearmanPValueMethod`` is
#: written on the wire. Only ``"not_estimable"`` rows have ``None`` for ``rho`` and ``p_value``.
SPEARMAN_METHODS = ("exact", "asymptotic", "not_estimable")

_PathLike = Union[str, "os.PathLike[str]"]


def __dir__() -> list[str]:
    """The public API only, so ``dir()`` is a map of what to call rather than of our imports."""
    return sorted(__all__)


def _columns(data: Mapping[str, Any]) -> tuple[list[str], dict[str, list[Any]]]:
    names = list(data.get("column_names") or [])
    return names, dict(data.get("columns") or {})


# ---- result types ------------------------------------------------------------------------------


@dataclass(frozen=True)
class VariancePrior:
    """The fitted prior of the per-feature residual variances, s²_g ~ s0² · χ²(d0)/d0.

    It depends only on the fit, never on which coefficient is tested, so one prior serves every
    coefficient in a :class:`ModeratedFit`.

    Attributes:
        df: Prior degrees of freedom d0, or ``None`` when it is infinite (see ``df_infinite``).
            Larger means the features' variances agree, so each is shrunk harder toward ``scale``.
            0 with fewer than 3 features: mzLib then fits no prior, as limma's ``squeezeVar``
            does, so each feature keeps its own variance.
        df_infinite: ``True`` when the residual variances are no more dispersed than sampling alone
            explains, so every feature takes the prior variance.
        trended: Whether s0² varies with each feature's average response (``trend=True``).
        spline_basis_count: Basis functions the trend spline used, intercept included; 1 without a
            trend.
        scale: s0², the prior variance, in squared response units, when every feature shares
            one. ``None`` when it differs by feature, which is then in the ``prior_variance``
            column: with a trend, or with fewer than 3 features, where no prior is fitted.
    """

    df: float | None
    df_infinite: bool
    trended: bool
    spline_basis_count: int
    scale: float | None

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> VariancePrior:
        return cls(
            df=data.get("df"),
            df_infinite=bool(data.get("df_infinite", False)),
            trended=bool(data.get("trended", False)),
            spline_basis_count=int(data.get("spline_basis_count", 1)),
            scale=data.get("scale"),
        )


@dataclass(frozen=True)
class ModeratedFit(_Table):
    """What :func:`fit` returns: one row per (tested coefficient, feature), with limma's moderated t.

    ``columns`` holds one block of ``feature_count`` rows per tested coefficient, in the order you
    named them, each block in the responses file's feature order:

    ======================  ==============================  ====================================================
    column                  unit                            meaning; ``None`` means
    ======================  ==============================  ====================================================
    ``feature``             -                               the id in the responses file's first column
    ``coefficient``         -                               the coefficient this row tests
    ``status``              -                               one of :data:`FIT_STATUSES`
    ``estimate``            response units per coefficient  least-squares coefficient: a **log2 fold change**
                            unit                            for a 0/1 group column on log2 intensities;
                                                            ``None`` = not fitted
    ``standard_error``      response units per coefficient  moderated standard error; ``None`` = not fitted
                            unit
    ``t``                   standard errors                 moderated t; ``None`` = not fitted
    ``df_total``            degrees of freedom              df_residual + prior df, capped at the pooled df
    ``p_value``             fraction (0 to 1)               two-sided; ``None`` = not fitted
    ``bh_adjusted``         fraction (0 to 1)               Benjamini-Hochberg over this coefficient's fitted
                                                            features; **not a target-decoy q-value**
    ``posterior_variance``  squared response units          the moderated residual variance
    ``prior_variance``      squared response units          s0² for this feature
    ``sigma``               response units                  residual standard deviation before moderation
    ``df_residual``         -                               observed samples minus coefficients
    ``observed``            samples                         samples with a finite response
    ``average_response``    response units                  the covariate ``trend=True`` fits on
    ======================  ==============================  ====================================================

    Attributes:
        responses_file: Absolute path of the responses table.
        design_file: Absolute path of the design table.
        feature_count: Feature rows read (features).
        sample_count: Sample columns read (samples).
        sample_names: In the responses header's order.
        coefficient_names: Every design column, in order: the model always fits all of them.
        tested_coefficients: The coefficients you asked to test, in order.
        trend: Whether the prior variance was allowed to trend with average response.
        threads: The thread count used (threads); the answer is identical at any count.
        missing_count: Response cells that were blank, NA or NaN (cells), left out of their fits.
        zero_count: Response cells that are exactly 0 (cells), fitted as observations.
        status_counts: Features per status (features), keyed by :data:`FIT_STATUSES`.
        residual_df_differ: ``True`` when the fitted features do not all have the same residual
            degrees of freedom. Then default limma (3.61 and later) would use a different prior
            estimator and report different moderated statistics; these are ``legacy = TRUE``.
        prior: The fitted :class:`VariancePrior`.
        row_count: Rows in the table (rows): ``feature_count`` x the number of tested coefficients.
        column_names: The table's columns, in order.
        columns: Column name -> list of values - the shape ``pandas.DataFrame`` takes as is.
        caveats: What the numbers are and are not, including any that apply to this input
            (unfitted features, zeros, ``residual_df_differ``).
    """

    responses_file: str
    design_file: str
    feature_count: int
    sample_count: int
    sample_names: list[str]
    coefficient_names: list[str]
    tested_coefficients: list[str]
    trend: bool
    threads: int
    missing_count: int
    zero_count: int
    status_counts: dict[str, int]
    residual_df_differ: bool
    prior: VariancePrior
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    def rows(self, coefficient: str) -> list[dict[str, Any]]:
        """The rows testing one coefficient, as dicts, in feature order.

        Args:
            coefficient: One of :attr:`tested_coefficients`.

        Raises:
            UsageError: ``coefficient`` was not tested in this fit.
        """
        if coefficient not in self.tested_coefficients:
            raise _bridge.UsageError(
                f"'{coefficient}' was not tested in this fit; it tested: "
                + ", ".join(self.tested_coefficients)
                + "."
            )
        return [r for r in self.records if r["coefficient"] == coefficient]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> ModeratedFit:
        names, columns = _columns(data)
        return cls(
            responses_file=str(data.get("responses_file", "")),
            design_file=str(data.get("design_file", "")),
            feature_count=int(data.get("feature_count", 0)),
            sample_count=int(data.get("sample_count", 0)),
            sample_names=list(data.get("sample_names") or []),
            coefficient_names=list(data.get("coefficient_names") or []),
            tested_coefficients=list(data.get("tested_coefficients") or []),
            trend=bool(data.get("trend", False)),
            threads=int(data.get("threads", 1)),
            missing_count=int(data.get("missing_count", 0)),
            zero_count=int(data.get("zero_count", 0)),
            status_counts=dict(data.get("status_counts") or {}),
            residual_df_differ=bool(data.get("residual_df_differ", False)),
            prior=VariancePrior._from_wire(data.get("prior") or {}),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class Adjusted(_Table):
    """What :func:`adjust` returns: one row per p-value given, in input order.

    ``columns["p_value"]`` is each p-value as read (a fraction, 0 to 1), ``None`` where you passed
    ``None`` or NaN; ``columns["bh_adjusted"]`` is its Benjamini-Hochberg adjusted value (a
    fraction, 0 to 1), ``None`` for an untested entry.

    Attributes:
        line_count: Values given (lines); the length of every column.
        tested_count: m, the finite p-values the adjustment is over (p-values).
        row_count: Rows in the table (rows), equal to ``line_count``.
        column_names: ``["p_value", "bh_adjusted"]``.
        columns: Column name -> one value per input.
        caveats: What the adjustment is and is not.
    """

    line_count: int
    tested_count: int
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @property
    def bh_adjusted(self) -> list[float | None]:
        """The adjusted values alone, aligned with the input (``None`` = not tested)."""
        return list(self.columns.get("bh_adjusted") or [])

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> Adjusted:
        names, columns = _columns(data)
        return cls(
            line_count=int(data.get("line_count", 0)),
            tested_count=int(data.get("tested_count", 0)),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class MetaAnalysis(_Table):
    """What :func:`meta` returns: one pooled row per feature, in order of first appearance.

    ===========================  ====================================  ===========================================
    column                       unit                                  meaning; ``None`` means
    ===========================  ====================================  ===========================================
    ``feature``                  -                                     as given
    ``studies``                  studies                               studies pooled
    ``estimate``                 units of the study estimates          random-effects pooled estimate
    ``standard_error``           units of the study estimates          its standard error
    ``confidence_low`` / _high   units of the study estimates          the interval at :attr:`confidence`
    ``p_value``                  fraction (0 to 1)                     two-sided, normal reference; not adjusted
                                                                       across features
    ``tau2``                     squared units of the study estimates  between-study variance; 0 for one study
    ``q``                        chi-square statistic                  Cochran's Q; ``None`` = one study
    ``i_squared``                fraction (0 to 1)                     share of variation from heterogeneity;
                                                                       ``None`` = one study
    ``direction_agree``          studies                               studies with the pooled estimate's sign
    ``direction_disagree``       studies                               studies with the opposite sign
    ``leave_one_out_max_delta``  units of the study estimates          largest change from dropping one study;
                                                                       ``None`` = one study
    ===========================  ====================================  ===========================================

    Attributes:
        study_count: Studies read, over every feature (studies).
        feature_count: Distinct features (features); the table's row count.
        confidence: The interval's coverage, as a fraction (0 to 1).
        row_count: Rows in the table (rows), equal to ``feature_count``.
        column_names: The table's columns, in order.
        columns: Column name -> one value per feature.
        caveats: What the pooling assumes.
    """

    study_count: int
    feature_count: int
    confidence: float
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> MetaAnalysis:
        names, columns = _columns(data)
        return cls(
            study_count=int(data.get("study_count", 0)),
            feature_count=int(data.get("feature_count", 0)),
            confidence=float(data.get("confidence", 0.95)),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class CombinedPValues(_Table):
    """What :func:`combine_p` returns: one combined row per feature, in order of first appearance.

    =============  ===============================  ====================================================
    column         unit                             meaning; ``None`` means
    =============  ===============================  ====================================================
    ``feature``    -                                as given
    ``studies``    tests                            tests combined: those with a finite p-value (after
                                                    the signed conversion, when effects were given)
    ``omitted``    tests                            this feature's tests that were not combined: no
                                                    p-value, or (signed) an effect of 0 or ``None``
    ``statistic``  chi-square on 2 x studies df     Fisher's X = -2 sum ln p, or Stouffer's
                   (fisher); standard normal        Z = sum w z / sqrt(sum w^2). ``None`` = not finite:
                   deviates (stouffer)              see below
    ``p_value``    fraction (0 to 1)                the upper tail of the statistic, so ONE direction
                                                    (effect > 0 when signed); ``None`` = nothing to
                                                    combine, or a p of 0 and a p of 1 together (stouffer)
    =============  ===============================  ====================================================

    A ``None`` statistic with a p-value is an infinite one: ``p_value`` 0 means plus infinity (a
    p-value of 0 among the tests), 1 means minus infinity (stouffer, a p-value of 1).

    Attributes:
        method: ``"fisher"`` or ``"stouffer"``, as you asked.
        weighted: Whether the tests carried weights (stouffer only).
        signed: Whether the tests carried effects, so every p-value was read as two-sided and
            converted to the one-sided p-value for effect > 0 before combining.
        line_count: Tests sent, over every feature (lines).
        feature_count: Distinct features (features); the table's row count.
        row_count: Rows in the table (rows), equal to ``feature_count``.
        column_names: The table's columns, in order.
        columns: Column name -> one value per feature.
        caveats: What the combination assumes, plus any that apply to this input (omitted tests,
            infinite statistics).
    """

    method: str
    weighted: bool
    signed: bool
    line_count: int
    feature_count: int
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @property
    def p_values(self) -> dict[str, float | None]:
        """Feature -> combined p-value (``None`` = nothing could be combined), in table order."""
        return dict(zip(self.columns.get("feature") or [], self.columns.get("p_value") or []))

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> CombinedPValues:
        names, columns = _columns(data)
        return cls(
            method=str(data.get("method", "")),
            weighted=bool(data.get("weighted", False)),
            signed=bool(data.get("signed", False)),
            line_count=int(data.get("line_count", 0)),
            feature_count=int(data.get("feature_count", 0)),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class Correlations(_Table):
    """What :func:`correlate` returns: one Spearman row per feature, in the responses file's order.

    ============  =========================  ===========================================================
    column        unit                       meaning; ``None`` means
    ============  =========================  ===========================================================
    ``feature``   -                          the id in the responses file's first column
    ``rho``       unitless (-1 to 1)         Spearman's rho, ties at their average rank; with
                                             ``given=``, the partial rho. ``None`` = not_estimable
    ``p_value``   fraction (0 to 1)          two-sided, against rho = 0; **not adjusted across
                                             features**. ``None`` = not_estimable
    ``n``         samples                    complete samples used: response, trait and every
                                             covariate finite
    ``has_ties``  -                          any variable tied among those samples (``False`` too
                                             when there were too few to rank: n below 3, or 4)
    ``method``    -                          one of :data:`SPEARMAN_METHODS`: ``"exact"`` (n <= 9, no
                                             ties, plain), ``"asymptotic"`` (t on n - 2 df, or
                                             n - 2 - k with k covariates kept) or ``"not_estimable"``
    ============  =========================  ===========================================================

    Attributes:
        responses_file: Absolute path of the responses table.
        trait_file: Absolute path of the trait table.
        given_file: Absolute path of the covariates table, or ``None`` without ``given=``.
        trait_name: The trait column's header.
        covariate_names: The covariate columns' headers, in order; empty without ``given=``.
        partial: ``True`` when ``given=`` was used, so every row is a partial correlation.
        feature_count: Feature rows read (features); the table's row count.
        sample_count: Sample columns read (samples).
        sample_names: In the responses header's order.
        missing_count: Response cells that were blank, NA or NaN (cells).
        trait_missing_count: Samples with no trait value (samples), so in no feature's correlation.
        given_incomplete_count: Samples missing a covariate (samples), so in no feature's
            correlation; ``None`` without ``given=``.
        method_counts: Features per method (features), keyed by :data:`SPEARMAN_METHODS`.
        row_count: Rows in the table (rows), equal to ``feature_count``.
        column_names: The table's columns, in order.
        columns: Column name -> one value per feature.
        caveats: How the p-values were computed and how they differ from R's defaults, plus a
            count of not_estimable features when there are any.
    """

    responses_file: str
    trait_file: str
    given_file: str | None
    trait_name: str
    covariate_names: list[str]
    partial: bool
    feature_count: int
    sample_count: int
    sample_names: list[str]
    missing_count: int
    trait_missing_count: int
    given_incomplete_count: int | None
    method_counts: dict[str, int]
    row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> Correlations:
        names, columns = _columns(data)
        incomplete = data.get("given_incomplete_count")
        return cls(
            responses_file=str(data.get("responses_file", "")),
            trait_file=str(data.get("trait_file", "")),
            given_file=data.get("given_file"),
            trait_name=str(data.get("trait_name", "")),
            covariate_names=list(data.get("covariate_names") or []),
            partial=bool(data.get("partial", False)),
            feature_count=int(data.get("feature_count", 0)),
            sample_count=int(data.get("sample_count", 0)),
            sample_names=list(data.get("sample_names") or []),
            missing_count=int(data.get("missing_count", 0)),
            trait_missing_count=int(data.get("trait_missing_count", 0)),
            given_incomplete_count=None if incomplete is None else int(incomplete),
            method_counts=dict(data.get("method_counts") or {}),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


# ---- calls -------------------------------------------------------------------------------------


def _one_line(text: str, what: str) -> str:
    if "\n" in text or "\r" in text:
        raise _bridge.UsageError(f"{what} contains a newline, which separates entries: {text!r}.")
    return text


def _number(value: object, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise _bridge.UsageError(f"{what} must be a number; got {value!r}.")
    return float(value)


def fit(
    responses: _PathLike,
    design: _PathLike,
    coefficients: Sequence[str],
    *,
    trend: bool = False,
    spline_basis: int | None = None,
    threads: int = 1,
    timeout: float | None = 600,
) -> ModeratedFit:
    """Fit one linear model per feature and test coefficients with limma's moderated t.

    Runs mzLib's ``LinearModel.Fit`` once - every design column is fitted - and then
    ``EmpiricalBayes.Moderate`` for each coefficient you name, which also applies
    Benjamini-Hochberg across that coefficient's fitted features. This is limma's
    ``lmFit(...)`` then ``eBayes(..., legacy = TRUE)``, reproduced by mzLib to 1e-8 relative.

    A feature is fitted on the samples where it was observed, and is reported - never dropped -
    when it cannot be: ``too_few_observations`` when it has no more observed samples than the
    design has coefficients, ``rank_deficient`` when its observed samples cannot separate the
    coefficients (a feature seen in only one group cannot have a group effect).

    Args:
        responses: A TSV with a header row. First column the feature id; every other column one
            sample, named by its header. Use the values as you want them modelled - log2
            intensities, usually. Blank, ``NA`` or ``NaN`` is missing; **0 is an observation**.
        design: A TSV with a header row and one row per sample: first column the sample name
            (matching the responses header exactly, in any order), then one numeric column per
            coefficient. Include an intercept column yourself if you want one.
        coefficients: Design column names to test, each with its own moderated test and its own
            Benjamini-Hochberg family.
        trend: Let the prior variance follow average response (limma's ``trend = TRUE``), for data
            whose low-abundance features are noisier.
        spline_basis: Basis functions of the trend spline, intercept included. Only with
            ``trend=True``. Leave it out to let mzLib choose, as limma does, for the number of
            features; mzLib may also use fewer than you ask for. Either way
            :attr:`VariancePrior.spline_basis_count` reports what it used.
        threads: Features fitted at once (threads), or ``-1`` for every core. The answer is
            identical at any value.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`ModeratedFit`.

    Raises:
        UsageError: a file is missing or malformed; the design's samples are not exactly the
            responses' samples; a coefficient is not in the design or is named twice; the design
            is not of full rank; or fewer than two features can be fitted, so no prior exists.

    Examples:
        A dilution series from mzLib's own RNA test data: MALAT1 loaded at 500, 250 and 125 ng
        against a constant 500 ng FLuc spike, log2 intensities normalised to FLuc in each run (see
        ``tests/fixtures/stats/README.md``). Halving MALAT1 should read as about -1:

        >>> fit = pymzlib.stats.fit("malat_dilution_log2_vs_fluc.tsv", "malat_dilution_design.tsv",
        ...                         ["malat_250ng", "malat_125ng"])
        >>> fit.feature_count, fit.sample_count, fit.status_counts
        (492, 27, {'fitted': 212, 'too_few_observations': 214, 'rank_deficient': 66})
        >>> fit.residual_df_differ, round(fit.prior.df, 1)
        (True, 19.6)
        >>> import statistics
        >>> malat = [r["estimate"] for r in fit.rows("malat_250ng")
        ...          if r["status"] == "fitted" and r["feature"].startswith("MALAT1:")]
        >>> len(malat), round(statistics.median(malat), 2)
        (116, -0.89)
    """
    if isinstance(coefficients, str):
        raise _bridge.UsageError(
            f"coefficients takes a list of design column names; for one, use [{coefficients!r}]."
        )
    if not isinstance(coefficients, Sequence) or not coefficients:
        raise _bridge.UsageError("fit() needs at least one coefficient to test, as a list.")
    for name in coefficients:
        if not isinstance(name, str) or not name.strip():
            raise _bridge.UsageError(f"Every coefficient must be a non-empty str; got {name!r}.")
        _one_line(name, "A coefficient name")
    if not isinstance(trend, bool):
        raise _bridge.UsageError(f"trend must be True or False; got {trend!r}.")
    if isinstance(threads, bool) or not isinstance(threads, int):
        raise _bridge.UsageError(f"threads must be an int; got {threads!r}.")

    args = ["stats", "fit", "--responses", _bridge.path_text(responses),
            "--design", _bridge.path_text(design), "--threads", str(threads)]
    if trend:
        args.append("--trend")
    if spline_basis is not None:
        if isinstance(spline_basis, bool) or not isinstance(spline_basis, int):
            raise _bridge.UsageError(f"spline_basis must be an int; got {spline_basis!r}.")
        args += ["--spline-basis", str(spline_basis)]

    _bridge.require_verb("stats fit", since=_SINCE)
    data = _bridge.invoke(*args, stdin="\n".join(coefficients) + "\n", timeout=timeout)
    return ModeratedFit._from_wire(data)


def adjust(p_values: Sequence[float | None], *, timeout: float | None = 60) -> Adjusted:
    """Benjamini-Hochberg adjust a list of p-values, keeping every position.

    For the i-th smallest of m p-values the adjusted value is min over j >= i of p_(j) m / j, capped
    at 1 (Benjamini and Hochberg 1995), computed by mzLib's ``MultipleTesting.BenjaminiHochberg``.
    ``None`` or NaN marks a feature that was not tested: it stays ``None`` and is **not counted in
    m**, which is the honest family. Do not write 1 for an untested feature; that enlarges m.

    :func:`fit` already reports ``bh_adjusted``; this is for p-values from anywhere else.

    Args:
        p_values: The p-values, each a fraction from 0 to 1, or ``None``/NaN for untested.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`Adjusted` table, row ``i`` for ``p_values[i]``.

    Raises:
        UsageError: ``p_values`` is empty or not a list, an entry is not a number, or a p-value is
            outside [0, 1].

    Examples:
        >>> adjusted = pymzlib.stats.adjust([0.0002, 0.004, 0.019, None, 0.031, 0.2, float("nan"), 0.74])
        >>> adjusted.tested_count
        6
        >>> [None if v is None else round(v, 4) for v in adjusted.bh_adjusted]
        [0.0012, 0.012, 0.038, None, 0.0465, 0.24, None, 0.74]
    """
    if isinstance(p_values, (str, bytes)) or not isinstance(p_values, Sequence) or not p_values:
        raise _bridge.UsageError("adjust() needs at least one p-value, as a list.")
    lines = []
    for i, value in enumerate(p_values):
        if value is None:
            lines.append("")
            continue
        number = _number(value, f"p_values[{i}]")
        lines.append("" if math.isnan(number) else repr(number))
    _bridge.require_verb("stats adjust", since=_SINCE)
    # A trailing newline ends the last line, so a final untested entry survives the trip.
    data = _bridge.invoke("stats", "adjust", stdin="\n".join(lines) + "\n", timeout=timeout)
    return Adjusted._from_wire(data)


def meta(
    studies: Iterable[Sequence[Any] | Mapping[str, Any]],
    *,
    confidence: float = 0.95,
    timeout: float | None = 60,
) -> MetaAnalysis:
    """Pool one effect size per study into a random-effects estimate per feature.

    Each feature is pooled on its own by mzLib's ``RandomEffectsMeta.Pool``, with the
    DerSimonian-Laird (1986) moment estimator of the between-study variance, as metafor's
    ``rma(method = "DL")`` does (mzLib matches metafor 5.2.1 to 1e-8 relative). Every row also says
    how many studies agree in direction and how far the estimate moves when any one is dropped -
    the two checks a reader of a pooled result asks for first.

    Args:
        studies: One entry per study: ``(feature, estimate, standard_error)``, or a mapping with
            those three keys. A feature's studies need not be adjacent. The standard error must be
            positive.
        confidence: Coverage of the reported interval, as a fraction (0 to 1).
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`MetaAnalysis` table, one row per feature in order of first appearance.

    Raises:
        UsageError: no studies; an entry without the three values; a feature name that is empty
            or holds a tab or newline; an estimate or standard error that is not a finite number,
            or a standard error that is not positive; ``confidence`` not strictly between 0 and 1.

    Examples:
        metafor's own DerSimonian-Laird reference cases, from mzLib's test data:

        >>> import csv
        >>> with open("metafor_dl_inputs.tsv", newline="") as fh:
        ...     studies = [(r["case"], float(r["yi"]), float(r["sei"]))
        ...                for r in csv.DictReader(fh, delimiter="\\t")]
        >>> pooled = pymzlib.stats.meta(studies)
        >>> for r in pooled.records:
        ...     print(r["feature"], r["studies"], round(r["estimate"], 4), round(r["tau2"], 4),
        ...           r["direction_agree"], r["direction_disagree"])
        homogeneous 5 0.3097 0 5 0
        heterogeneous 6 0.2513 0.0752 5 1
        two_studies 2 -0.2 0.06 1 1
        random_eight 8 0.107 0.0911 6 2
    """
    if isinstance(confidence, bool) or not isinstance(confidence, numbers.Real):
        raise _bridge.UsageError(f"confidence must be a number; got {confidence!r}.")
    if isinstance(studies, (str, bytes, Mapping)):
        raise _bridge.UsageError(
            "meta() takes a list of studies, each (feature, estimate, standard_error)."
        )
    lines = []
    for i, study in enumerate(studies):
        if isinstance(study, Mapping):
            try:
                feature, estimate, error = study["feature"], study["estimate"], study["standard_error"]
            except KeyError as missing:
                raise _bridge.UsageError(
                    f"studies[{i}] has no {missing.args[0]!r} key; a mapping needs feature, "
                    "estimate and standard_error."
                ) from None
        elif isinstance(study, Sequence) and not isinstance(study, (str, bytes)) and len(study) == 3:
            feature, estimate, error = study
        else:
            raise _bridge.UsageError(
                f"studies[{i}] must be (feature, estimate, standard_error); got {study!r}."
            )
        if not isinstance(feature, str) or not feature.strip():
            raise _bridge.UsageError(f"studies[{i}] needs a non-empty feature name; got {feature!r}.")
        if "\t" in feature:
            raise _bridge.UsageError(f"studies[{i}]'s feature name contains a tab: {feature!r}.")
        _one_line(feature, f"studies[{i}]'s feature name")
        lines.append(
            f"{feature}\t{_number(estimate, f'studies[{i}] estimate')!r}"
            f"\t{_number(error, f'studies[{i}] standard_error')!r}"
        )
    if not lines:
        raise _bridge.UsageError("meta() needs at least one study.")
    _bridge.require_verb("stats meta", since=_SINCE)
    data = _bridge.invoke(
        "stats", "meta", "--confidence", repr(float(confidence)),
        stdin="\n".join(lines) + "\n", timeout=timeout,
    )
    return MetaAnalysis._from_wire(data)


def _path(value: object, what: str) -> str:
    text = _bridge.path_text(value)
    if not text:
        raise _bridge.UsageError(f"{what} must be a path (str or os.PathLike); got {value!r}.")
    return text


def _cell(value: object, what: str) -> str:
    """A numeric field of a stdin line: ``None`` and NaN are blank (missing), a number its repr."""
    if value is None:
        return ""
    number = _number(value, what)
    return "" if math.isnan(number) else repr(number)


def combine_p(
    tests: Iterable[Sequence[Any] | Mapping[str, Any]],
    *,
    method: str,
    timeout: float | None = 60,
) -> CombinedPValues:
    """Combine several independent p-values per feature into one, by Fisher's or Stouffer's method.

    Each feature is combined on its own by mzLib's ``PValueCombination``: Fisher (1925),
    X = -2 sum ln p against chi-square on 2k degrees of freedom, or Stouffer (1949),
    Z = sum w z / sqrt(sum w^2) against the standard normal. mzLib holds both to R's ``pchisq`` and
    ``qnorm`` to 1e-8 relative. Use it when you have one p-value per study - per dataset, cohort or
    laboratory - for the same feature, and want one answer from all of them.

    **Direction.** Both rules combine *one-sided* p-values. A two-sided p-value from :func:`fit`
    does not say which way the effect went, so two studies that disagree look like two that agree.
    Give every test its ``effect`` (its estimate; only the sign is used) and mzLib first turns each
    p into the one-sided p-value for effect > 0 (``PValueCombination.OneSided``: p/2 for a positive
    effect, 1 - p/2 for a negative one), so opposite effects cancel. The combined ``p_value`` then
    tests effect > 0; for Stouffer, a ``p_value`` near 1 means a consistently negative effect.

    Neither rule estimates an effect size: :func:`meta` does that, from estimates and standard
    errors. Neither adjusts across features: send ``p_values`` to :func:`adjust`.

    Args:
        tests: One entry per test: ``(feature, p_value)``, ``(feature, p_value, weight)`` or
            ``(feature, p_value, weight, effect)``, or a mapping with ``feature`` and ``p_value``
            and optionally ``weight`` and ``effect``. A feature's tests need not be adjacent.
            ``p_value`` is a fraction from 0 to 1, or ``None``/NaN for a test that could not be run
            (it is omitted, and counted in the ``omitted`` column). ``weight`` (Stouffer only) is a
            positive number, for example the square root of the study's sample size; give one for
            every test or for none, and ``None`` for "no weight". ``effect`` is the test's estimate
            in any unit; give one for every test or for none, and ``None`` for a test with no
            direction (omitted). To give an effect without a weight, pass ``None`` as the weight.
        method: ``"fisher"`` or ``"stouffer"``, one of :data:`COMBINATION_METHODS`. Required:
            mzLib has no preferred rule. Fisher is driven by the strongest single result and takes
            no weights; Stouffer weighs every study, and takes weights.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`CombinedPValues` table, one row per feature in order of first appearance;
        :attr:`CombinedPValues.p_values` maps each feature to its combined p-value.

    Raises:
        UsageError: no tests; an entry that is not 2 to 4 values or a mapping with ``feature`` and
            ``p_value``; a feature name that is empty or holds a tab or newline; a value that is not
            a number; ``method`` not one of :data:`COMBINATION_METHODS`; a weight or an effect on
            some tests and not others; a weight with ``method="fisher"``; or, from mzLib, naming the
            feature, a p-value outside [0, 1] or a weight that is not finite and positive.

    Examples:
        mzLib's own reference cases (constructed p-values, not proteomics data), by both rules:

        >>> import csv
        >>> with open("combination_inputs.tsv", newline="") as fh:
        ...     rows = list(csv.DictReader(fh, delimiter="\\t"))
        >>> fisher = pymzlib.stats.combine_p([(r["case"], float(r["p"])) for r in rows],
        ...                                  method="fisher")
        >>> stouffer = pymzlib.stats.combine_p(
        ...     [(r["case"], float(r["p"]), float(r["weight"])) for r in rows], method="stouffer")
        >>> for f, s in zip(fisher.records, stouffer.records):
        ...     print(f["feature"], f["studies"], f"{f['p_value']:.3g}", f"{s['p_value']:.3g}")
        few 3 0.0871 0.15
        tiny 4 3.14e-11 0.00228
        one 1 0.037 0.037
        many 12 0.113 0.0323
        mixed 5 0.0191 0.713
        tail 3 2.14e-18 0.000836

        ``mixed`` shows the two rules apart. One p-value of 0.001 among five carries Fisher's to
        0.019; weighted Stouffer counts all five, and its heaviest test (p = 0.95) points the
        other way. Both agree with R:

        >>> with open("combination_results.tsv", newline="") as fh:
        ...     r_answers = list(csv.DictReader(fh, delimiter="\\t"))
        >>> def worst(ours, theirs):
        ...     return max(abs(a - float(b)) / abs(float(b)) for a, b in zip(ours, theirs))
        >>> worst(fisher.columns["p_value"], [x["fisher_p"] for x in r_answers]) < 1e-8
        True
        >>> worst(stouffer.columns["p_value"], [x["weighted_p"] for x in r_answers]) < 1e-8
        True
    """
    if not isinstance(method, str) or method not in COMBINATION_METHODS:
        raise _bridge.UsageError(
            f"method must be one of {', '.join(COMBINATION_METHODS)} (mzLib's "
            f"PValueCombinationMethod); got {method!r}."
        )
    if isinstance(tests, (str, bytes, Mapping)):
        raise _bridge.UsageError(
            "combine_p() takes a list of tests, each (feature, p_value[, weight[, effect]])."
        )
    names = ("p_value", "weight", "effect")
    lines = []
    for i, test in enumerate(tests):
        if isinstance(test, Mapping):
            if "feature" not in test or "p_value" not in test:
                raise _bridge.UsageError(
                    f"tests[{i}] needs 'feature' and 'p_value' keys, and optionally 'weight' and "
                    f"'effect'; got {sorted(map(str, test))}."
                )
            size = 4 if "effect" in test else 3 if "weight" in test else 2
            values = [test["feature"], test["p_value"], test.get("weight"), test.get("effect")][:size]
        elif isinstance(test, Sequence) and not isinstance(test, (str, bytes)) and 2 <= len(test) <= 4:
            values = list(test)
        else:
            raise _bridge.UsageError(
                f"tests[{i}] must be (feature, p_value[, weight[, effect]]); got {test!r}."
            )
        feature = values[0]
        if not isinstance(feature, str) or not feature.strip():
            raise _bridge.UsageError(f"tests[{i}] needs a non-empty feature name; got {feature!r}.")
        if "\t" in feature:
            raise _bridge.UsageError(f"tests[{i}]'s feature name contains a tab: {feature!r}.")
        _one_line(feature, f"tests[{i}]'s feature name")
        cells = [_cell(value, f"tests[{i}] {names[j]}") for j, value in enumerate(values[1:])]
        lines.append("\t".join([feature, *cells]))
    if not lines:
        raise _bridge.UsageError("combine_p() needs at least one test.")
    _bridge.require_verb("stats combine-p", since=_COMBINE_CORRELATE_SINCE)
    data = _bridge.invoke(
        "stats", "combine-p", "--method", method, stdin="\n".join(lines) + "\n", timeout=timeout
    )
    return CombinedPValues._from_wire(data)


def correlate(
    responses: _PathLike,
    trait: _PathLike,
    *,
    given: _PathLike | None = None,
    timeout: float | None = 600,
) -> Correlations:
    """Spearman-correlate every feature with a per-sample trait, optionally given covariates.

    For each feature (row) of the responses table, mzLib's ``SpearmanCorrelation.Correlate`` ranks
    the feature's values and the trait over the samples where both are known, and tests rho = 0,
    two-sided. With ``given``, ``SpearmanCorrelation.PartialCorrelate`` first removes from both sets
    of ranks what the covariates' ranks explain, as ``ppcor::pcor.test(method = "spearman")`` does.
    mzLib holds both to R's ``cor.test`` and to ppcor 1.1 to 1e-8 relative.

    Spearman uses ranks only. It needs no log transform, one extreme value cannot drive it, and it
    finds any relationship that only ever rises (or only ever falls), not just a straight line.

    Args:
        responses: :func:`fit`'s responses table: a TSV with a header row, the feature id first,
            then one column per sample, named by its header. Blank, ``NA`` or ``NaN`` is missing.
            **0 is an observation**, ranked like any other value.
        trait: A TSV with a header row and one row per sample: the sample name (exactly the
            responses' sample columns, in any order), then ONE value column, whose header names the
            trait - age, dose, loading. A blank or ``NA`` cell leaves that sample out.
        given: A TSV of the same shape with one or more covariate columns (batch, run order,
            loading) to partial out. Leave it out (``None``) for a plain correlation.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`Correlations` table, one row per feature in the responses file's order.

    Raises:
        UsageError: a file is missing or malformed; the trait file has other than one value column;
            the trait's or the covariates' samples are not exactly the responses' samples; a cell
            is not a number.

    Examples:
        The MALAT1 dilution series :func:`fit` analyses (500, 250 and 125 ng MALAT1 against a
        constant FLuc spike, log2 intensities normalised to FLuc), correlated with each run's
        MALAT1 load:

        >>> result = pymzlib.stats.correlate("malat_dilution_log2_vs_fluc.tsv",
        ...                                  "malat_dilution_load_ng.tsv")
        >>> result.trait_name, result.feature_count, result.sample_count
        ('malat_ng', 492, 27)
        >>> result.method_counts
        {'exact': 4, 'asymptotic': 308, 'not_estimable': 180}
        >>> import statistics
        >>> def median_rho(prefix):
        ...     return statistics.median(r["rho"] for r in result.records
        ...                              if r["rho"] is not None and r["feature"].startswith(prefix))
        >>> print(f"MALAT1 {median_rho('MALAT1:'):+.2f}, FLuc {median_rho('FLuc:'):+.2f}")
        MALAT1 +0.43, FLuc +0.00

        MALAT1's oligos rise with its load. FLuc's, the constant spike every run was normalised
        to, do not.
    """
    args = ["stats", "correlate", "--responses", _path(responses, "responses"),
            "--trait", _path(trait, "trait")]
    if given is not None:
        args += ["--given", _path(given, "given")]
    _bridge.require_verb("stats correlate", since=_COMBINE_CORRELATE_SINCE)
    data = _bridge.invoke(*args, timeout=timeout)
    return Correlations._from_wire(data)
