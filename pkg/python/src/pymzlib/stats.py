"""Differential abundance with no R: limma's moderated t-test, Benjamini-Hochberg, and random-effects
meta-analysis, computed by mzLib.

Three questions, each answered by mzLib's ``StatisticalModels`` (mzLib #1341 and #1357):

==========================================================  =====================  =============================================
You want to know                                            Call                   mzLib computes it with
==========================================================  =====================  =============================================
Which features change with a coefficient, and how surely?   :func:`fit`            ``LinearModel.Fit`` + ``EmpiricalBayes.Moderate``
Which of these p-values survive a false discovery rate?     :func:`adjust`         ``MultipleTesting.BenjaminiHochberg``
What is one feature's effect, pooled across studies?        :func:`meta`           ``RandomEffectsMeta.Pool`` (DerSimonian-Laird)
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

Wire verbs: ``stats fit``, ``stats adjust`` and ``stats meta``; their language-neutral specs live in
the bridge repository under ``design/verbs/``.
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
    "FIT_STATUSES",
    "Adjusted",
    "MetaAnalysis",
    "ModeratedFit",
    "VariancePrior",
    "adjust",
    "fit",
    "meta",
]

#: The first pyMzLib whose bridge has the stats verbs (``since.pymzlib`` of their specs).
_SINCE = "0.3.0"

#: Every value of the ``status`` column, as mzLib's ``FeatureFitStatus`` is written on the wire.
#: Only ``"fitted"`` rows carry statistics or enter the Benjamini-Hochberg family.
FIT_STATUSES = ("fitted", "too_few_observations", "rank_deficient")

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
        df_infinite: ``True`` when the residual variances are no more dispersed than sampling alone
            explains, so every feature takes the prior variance.
        trended: Whether s0² varies with each feature's average response (``trend=True``).
        spline_basis_count: Basis functions the trend spline used, intercept included; 1 without a
            trend.
        scale: s0², the prior variance, in squared response units; ``None`` when trended, in which
            case it is per feature, in the ``prior_variance`` column.
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
        spline_basis: Basis functions of the trend spline, intercept included (default 4). Only
            with ``trend=True``; mzLib may use fewer, and :attr:`VariancePrior.spline_basis_count`
            reports what it used.
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
