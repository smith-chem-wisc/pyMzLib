"""Predict peptide properties: retention time, fragment intensities, CCS, detectability.

mzLib ships clients for **37 published models** on the `Koina <https://koina.wilhelmlab.org/>`_
inference server, across five families. This module calls them::

    >>> import pymzlib
    >>> r = pymzlib.prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK", "ELVISLIVESK"])
    ... # doctest: +SKIP
    >>> r.columns["retention_time"], r.retention_time_unit                       # doctest: +SKIP
    ([5.5165, 129.7723], 'indexed_retention_time')

Start with :func:`models`. It is enumerated from mzLib rather than transcribed, and each entry
carries the constraints that decide whether a peptide can be sent at all - length bounds, allowed
charges, which UNIMOD modifications the model was trained on, and whether a collision energy or
instrument type is required.

Five things worth knowing before the first call
-----------------------------------------------

**Koina is someone else's GPU.** It is a public, shared, community-run service, free and
unauthenticated. The throttling defaults are mzLib's and are deliberately not raised here;
``max_batches`` and ``throttle_ms`` exist for a genuinely large job, not for going faster by
default.

**A prediction is an opinion, not a measurement.** Nothing here has been matched against a
spectrum, and no output is FDR-anything. That is obvious stated plainly and surprisingly easy to
forget once the numbers are in a DataFrame next to real ones.

**Retention time is not always in minutes.** Most of these models are trained on *indexed*
retention time - iRT, a dimensionless scale anchored to standard peptides - and only
``Chronologer_RT`` returns absolute minutes. :attr:`Predictions.retention_time_unit` tells you
which, as a value. Plotting an iRT against a gradient without fitting the iRT-to-minutes line
first is the commonest way to misread these numbers.

**Fragment arrays are ragged.** Koina returns a fixed-width grid with ``-1`` marking ions that
cannot exist for a given peptide, and mzLib drops those. So each row's three fragment arrays are
as long as *that peptide's* possible ions - 28 for a short tryptic peptide from a model whose
nominal count is 174 - and two rows in one call will differ. Index them per row, never as a
rectangle.

**A peptide that cannot be predicted still gets a row.** Too long, an unsupported modification, a
missing required collision energy: the value comes back ``None`` and ``warning`` says why, rather
than the row vanishing. Predictions therefore always line up with the peptides you sent.

.. note::

   **The local Chronologer predictor is deliberately not exposed.** mzLib also ships Chronologer as
   a local TorchSharp network: x64-only, extracting hundreds of megabytes of weights to a shared
   temporary path, and racing any concurrent process doing the same. pyMzLib publishes an arm64
   macOS wheel, so exposing it would either break that wheel or ship a function that fails on it.
   The same model is reachable over Koina as ``Chronologer_RT`` - but note the two report
   **different units**, absolute retention time over the network and % acetonitrile locally.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

from . import _bridge

__all__ = [
    "Constraint",
    "Model",
    "Predictions",
    "models",
    "retention_time",
    "fragments",
    "ccs",
    "detectability",
    "crosslink_fragments",
]


def __dir__() -> list[str]:
    """The public API only."""
    return sorted(__all__)


@dataclass(frozen=True)
class Constraint:
    """What a model requires of one optional input parameter.

    Attributes:
        requirement: One of ``"not_applicable"`` (do not send this parameter - the model has no
            such input), ``"any_value_required"`` (you must send one, any value), or ``"one_of"``
            (you must send one of :attr:`values`).
        values: The accepted values, or ``None`` when the requirement does not name any.

    Note:
        This is a **tri-state**, and mzLib expresses it as a nullable set whose emptiness means the
        opposite of what it looks like: ``null`` means *not applicable* and an *empty* set means
        *required, any value*. Reading the raw collection is how you conclude that
        ``Prosit_2020_intensity_CID`` accepts any collision energy (it accepts none - it is fixed
        at NCE 35) and that ``Prosit_2020_intensity_HCD`` accepts none (it requires one).
    """

    requirement: str
    values: Any = None

    @property
    def applicable(self) -> bool:
        """Whether this parameter should be sent at all."""
        return self.requirement != "not_applicable"

    def accepts(self, value: Any) -> bool:
        """Whether a particular value satisfies this constraint."""
        if self.requirement == "not_applicable":
            return value is None
        if self.requirement == "any_value_required":
            return value is not None
        return value in (self.values or [])

    @classmethod
    def _from_wire(cls, payload: Any) -> "Constraint":
        if not isinstance(payload, dict):
            return cls(requirement="not_applicable")
        return cls(
            requirement=payload.get("requirement") or "not_applicable",
            values=payload.get("values"),
        )


@dataclass(frozen=True)
class Model:
    """One Koina model mzLib can call, with the constraints that decide what you may send it.

    Attributes:
        model: The model's published Koina name - what you pass as ``model``.
        family: ``"retention_time"``, ``"fragment_intensity"``, ``"collisional_cross_section"``,
            ``"detectability"`` or ``"crosslink_intensity"``.
        verb: The bridge verb this family is called through, for cross-referencing.
        type: The mzLib class name, for cross-referencing the mzLib source.
        min_peptide_length: Shortest base sequence the model accepts.
        max_peptide_length: Longest. Frequently 30, which excludes a lot of real tryptic peptides.
        max_batch_size: Sequences per request, as the server accepts them. mzLib batches for you.
        allowed_unimod_ids: UNIMOD accessions the model was trained on. **Empty means the model
            accepts no modifications at all**, which is a real answer, not a missing one.
        precursor_charge: A :class:`Constraint`.
        collision_energy: A :class:`Constraint`.
        instrument_type: A :class:`Constraint`.
        fragmentation_type: A :class:`Constraint`.
        retention_time_unit: ``"indexed_retention_time"`` or ``"minutes"`` for the retention-time
            family; ``None`` for every other family.
        number_of_predicted_fragment_ions: The model's nominal ion count, or ``None`` when it is
            dynamic. **Not** the length of any row's fragment arrays - see the module docstring.
    """

    model: str
    family: str
    verb: Any = None
    type: Any = None
    min_peptide_length: Any = None
    max_peptide_length: Any = None
    max_batch_size: Any = None
    allowed_unimod_ids: list[int] = field(default_factory=list)
    precursor_charge: Constraint = field(default_factory=lambda: Constraint("not_applicable"))
    collision_energy: Constraint = field(default_factory=lambda: Constraint("not_applicable"))
    instrument_type: Constraint = field(default_factory=lambda: Constraint("not_applicable"))
    fragmentation_type: Constraint = field(default_factory=lambda: Constraint("not_applicable"))
    retention_time_unit: Any = None
    number_of_predicted_fragment_ions: Any = None
    error: Any = None

    @property
    def accepts_modifications(self) -> bool:
        """Whether the model was trained on any modifications at all."""
        return bool(self.allowed_unimod_ids)

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "Model":
        return cls(
            model=payload.get("model") or "",
            family=payload.get("family") or "unknown",
            verb=payload.get("verb"),
            type=payload.get("type"),
            min_peptide_length=payload.get("min_peptide_length"),
            max_peptide_length=payload.get("max_peptide_length"),
            max_batch_size=payload.get("max_batch_size"),
            allowed_unimod_ids=list(payload.get("allowed_unimod_ids") or []),
            precursor_charge=Constraint._from_wire(payload.get("precursor_charge")),
            collision_energy=Constraint._from_wire(payload.get("collision_energy")),
            instrument_type=Constraint._from_wire(payload.get("instrument_type")),
            fragmentation_type=Constraint._from_wire(payload.get("fragmentation_type")),
            retention_time_unit=payload.get("retention_time_unit"),
            number_of_predicted_fragment_ions=payload.get("number_of_predicted_fragment_ions"),
            error=payload.get("error"),
        )


@dataclass(frozen=True)
class Predictions:
    """A table of predictions, one row per peptide sent.

    Attributes:
        model: The model that produced them.
        row_count: Rows returned - always equal to the number of peptides sent, so predictions line
            up with inputs even where some could not be predicted.
        failed_row_count: Rows whose prediction is ``None`` with a ``warning`` explaining why. **Not
            an error**: too long, an unsupported modification, or a missing required parameter are
            normal outcomes.
        column_names: The field names, in order.
        columns: Field name -> list of values, the shape a DataFrame accepts directly. ``None``
            when ``out`` was given.
        output: Where the table was written, or ``None``.
        retention_time_unit: For :func:`retention_time` only: ``"indexed_retention_time"`` or
            ``"minutes"``. See the module docstring - this is not a formality.
        collisional_cross_section_unit: For :func:`ccs` only: ``"square_angstroms"``, never 1/K0.
        intensity_scale: For the fragment verbs: ``"relative"``.
        caveats: What these numbers cannot be trusted to mean.
    """

    model: str
    row_count: int
    failed_row_count: int
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None
    retention_time_unit: Any = None
    collisional_cross_section_unit: Any = None
    intensity_scale: Any = None
    caveats: list[str] = field(default_factory=list)

    @property
    def records(self) -> list[dict[str, Any]]:
        """The same data row-wise: one dict per peptide."""
        if not self.columns:
            return []
        names = [n for n in (self.column_names or list(self.columns)) if n in self.columns]
        length = min((len(self.columns[n]) for n in names), default=0)
        return [{name: self.columns[name][i] for name in names} for i in range(length)]

    @property
    def warnings(self) -> list[tuple[int, str]]:
        """``(row index, message)`` for every row that could not be predicted."""
        if not self.columns or "warning" not in self.columns:
            return []
        return [(i, w) for i, w in enumerate(self.columns["warning"]) if w]

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "Predictions":
        written = data.get("output")
        return cls(
            model=data.get("model") or "",
            row_count=int(data.get("row_count", 0)),
            failed_row_count=int(data.get("failed_row_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=written,
            retention_time_unit=data.get("retention_time_unit"),
            collisional_cross_section_unit=data.get("collisional_cross_section_unit"),
            intensity_scale=data.get("intensity_scale"),
            caveats=list(data.get("caveats") or []),
        )


def models(family: str | None = None, timeout: float | None = 60) -> list[Model]:
    """Every Koina model mzLib can call.

    Enumerated from mzLib itself, so it reflects the installed version and cannot go stale.

    Args:
        family: Restrict to one family, e.g. ``"retention_time"``. ``None`` returns all 37.
        timeout: Seconds to allow.

    Returns:
        One :class:`Model` per model. **This call does not touch the network** - it describes what
        mzLib can call, not what the server currently answers.

    Example:
        >>> [m.model for m in models("collisional_cross_section")]     # doctest: +SKIP
        ['AlphaPeptDeep_ccs_generic', 'IM2Deep']
    """
    args = ["predict", "models"]
    if family is not None:
        if not isinstance(family, str) or not family.strip():
            raise _bridge.UsageError("family must be a non-empty string or None.")
        args += ["--family", family.strip()]

    data = _bridge.invoke(*args, timeout=timeout)
    return [Model._from_wire(item) for item in (data.get("models") or [])]


def _rows_to_stdin(rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> str:
    """The input table as the tab-separated text the bridge reads from stdin.

    On stdin rather than argv because a real prediction run is thousands of peptides and argv has a
    hard ceiling of roughly 32 KB - the same reason PRIDE's explicit file selection goes this way.
    """
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if row.get(c) is None else str(row[c]) for c in columns))
    return "\n".join(lines) + "\n"


def _predict(
    verb: str,
    model: str,
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
    *,
    out: str | None,
    max_batches: int | None,
    throttle_ms: int | None,
    timeout: float | None,
) -> Predictions:
    """The shared body of every predict verb: validate, send, parse."""
    if not isinstance(model, str) or not model.strip():
        raise _bridge.UsageError(
            "A model name is required. prediction.models() lists them with their constraints."
        )
    if not rows:
        raise _bridge.UsageError("At least one peptide is required.")

    args = ["predict", verb, "--model", model.strip()]

    if out is not None:
        if not isinstance(out, str) or not out.strip():
            raise _bridge.UsageError("out must be a non-empty path or None.")
        args += ["--out", out.strip()]

    # Deliberately not defaulted to anything faster than mzLib's own values: Koina is a shared
    # community server, and a binding that maximised throughput by default would be spending
    # someone else's GPU time.
    for name, value, minimum in (("max-batches", max_batches, 1), ("throttle-ms", throttle_ms, 0)):
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise _bridge.UsageError(
                f"{name.replace('-', '_')} must be a whole number of {minimum} or more, or None; "
                f"got {value!r}."
            )
        args += [f"--{name}", str(value)]

    data = _bridge.invoke(*args, stdin=_rows_to_stdin(rows, columns), timeout=timeout)
    return Predictions._from_wire(data)


def _as_rows(
    peptides: Iterable[Any], sequence_column: str, **defaults: Any
) -> list[dict[str, Any]]:
    """Accept either bare sequences or dicts, so the simple call stays simple.

    ``["PEPTIDEK", "ELVISLIVESK"]`` and
    ``[{"sequence": "PEPTIDEK", "precursor_charge": 2}]`` are both natural, and requiring the
    second shape for a one-column model would be gratuitous ceremony.
    """
    rows: list[dict[str, Any]] = []
    for peptide in peptides:
        if isinstance(peptide, str):
            row: dict[str, Any] = {sequence_column: peptide}
        elif isinstance(peptide, Mapping):
            row = dict(peptide)
        else:
            raise _bridge.UsageError(
                f"Each peptide must be a sequence string or a mapping of columns; got {peptide!r}."
            )
        for key, value in defaults.items():
            row.setdefault(key, value)
        rows.append(row)

    return rows


def retention_time(
    model: str,
    peptides: Iterable[Any],
    *,
    out: str | None = None,
    max_batches: int | None = None,
    throttle_ms: int | None = None,
    timeout: float | None = None,
) -> Predictions:
    """Predict elution, one row per peptide.

    Args:
        model: A retention-time model name from :func:`models`, e.g. ``"Prosit_2019_irt"``.
        peptides: Sequences, in mzLib ``FullSequence`` notation.
        out: Write a tab-separated table here and return only a summary.
        max_batches: In-flight requests. Leave ``None`` for mzLib's polite default.
        throttle_ms: Delay between request chunks. Leave ``None`` for mzLib's default.
        timeout: Seconds to allow. A large batch legitimately takes a while.

    Returns:
        A :class:`Predictions`. **Check** :attr:`Predictions.retention_time_unit` - most of these
        models return iRT, not minutes.

    Example:
        >>> r = retention_time("Prosit_2019_irt", ["PEPTIDEK"])        # doctest: +SKIP
        >>> r.retention_time_unit                                      # doctest: +SKIP
        'indexed_retention_time'
    """
    return _predict(
        "retention-time", model, _as_rows(peptides, "sequence"), ["sequence"],
        out=out, max_batches=max_batches, throttle_ms=throttle_ms, timeout=timeout,
    )


def fragments(
    model: str,
    peptides: Iterable[Any],
    *,
    precursor_charge: int | None = None,
    collision_energy: int | None = None,
    instrument_type: str | None = None,
    fragmentation_type: str | None = None,
    out: str | None = None,
    max_batches: int | None = None,
    throttle_ms: int | None = None,
    timeout: float | None = None,
) -> Predictions:
    """Predict MS2 fragment m/z and relative intensity.

    Args:
        model: A fragment-intensity model name, e.g. ``"Prosit_2020_intensity_HCD"``.
        peptides: Sequences, or mappings carrying per-peptide ``precursor_charge`` and the rest.
        precursor_charge: Applied to every peptide that does not carry its own.
        collision_energy: Required by many models - check
            :attr:`Model.collision_energy`. A peptide missing a required one comes back with a
            warning rather than failing the call.
        instrument_type: Required by a few, e.g. ``"QE"``, ``"LUMOS"``.
        fragmentation_type: Required by a few, e.g. ``"HCD"``, ``"CID"``.
        out: Write a tab-separated table here; the arrays become ``;``-joined lists.
        max_batches: In-flight requests.
        throttle_ms: Delay between request chunks.
        timeout: Seconds to allow.

    Returns:
        A :class:`Predictions` whose ``fragment_annotations``, ``fragment_mz`` and
        ``fragment_intensity`` columns each hold **one array per row**, of differing lengths.

    Example:
        >>> f = fragments("Prosit_2020_intensity_HCD",                 # doctest: +SKIP
        ...               ["PEPTIDEK", "ELVISLIVESK"],
        ...               precursor_charge=2, collision_energy=28)
        >>> [len(a) for a in f.columns["fragment_mz"]]                 # doctest: +SKIP
        [28, 40]
    """
    rows = _as_rows(
        peptides, "sequence",
        precursor_charge=precursor_charge,
        collision_energy=collision_energy,
        instrument_type=instrument_type,
        fragmentation_type=fragmentation_type,
    )
    return _predict(
        "fragments", model, rows,
        ["sequence", "precursor_charge", "collision_energy", "instrument_type", "fragmentation_type"],
        out=out, max_batches=max_batches, throttle_ms=throttle_ms, timeout=timeout,
    )


def ccs(
    model: str,
    peptides: Iterable[Any],
    *,
    precursor_charge: int | None = None,
    out: str | None = None,
    max_batches: int | None = None,
    throttle_ms: int | None = None,
    timeout: float | None = None,
) -> Predictions:
    """Predict collisional cross-section, in square angstroms.

    Args:
        model: A CCS model name, e.g. ``"IM2Deep"``.
        peptides: Sequences, or mappings carrying per-peptide ``precursor_charge``.
        precursor_charge: Applied to every peptide that does not carry its own. **CCS depends on
            charge**, so this is not optional in practice.
        out: Write a tab-separated table here.
        max_batches: In-flight requests.
        throttle_ms: Delay between request chunks.
        timeout: Seconds to allow.

    Returns:
        A :class:`Predictions`. The unit is **square angstroms, never 1/K0** - converting to the
        reduced mobility a timsTOF reports needs drift-gas parameters mzLib does not carry.
    """
    rows = _as_rows(peptides, "sequence", precursor_charge=precursor_charge)
    return _predict(
        "ccs", model, rows, ["sequence", "precursor_charge"],
        out=out, max_batches=max_batches, throttle_ms=throttle_ms, timeout=timeout,
    )


def detectability(
    model: str,
    peptides: Iterable[Any],
    *,
    out: str | None = None,
    max_batches: int | None = None,
    throttle_ms: int | None = None,
    timeout: float | None = None,
) -> Predictions:
    """Predict flyability, as four class probabilities that sum to 1.

    Args:
        model: A detectability model name, e.g. ``"pfly_2024_fine_tuned"``.
        peptides: Sequences.
        out: Write a tab-separated table here.
        max_batches: In-flight requests.
        throttle_ms: Delay between request chunks.
        timeout: Seconds to allow.

    Returns:
        A :class:`Predictions` with ``not_detectable``, ``low_detectability``,
        ``intermediate_detectability`` and ``high_detectability`` columns. They are a distribution
        over classes, **not** an expected intensity and not a probability of detection.
    """
    return _predict(
        "detectability", model, _as_rows(peptides, "sequence"), ["sequence"],
        out=out, max_batches=max_batches, throttle_ms=throttle_ms, timeout=timeout,
    )


def crosslink_fragments(
    model: str,
    pairs: Iterable[Mapping[str, Any]],
    *,
    precursor_charge: int | None = None,
    collision_energy: int | None = None,
    out: str | None = None,
    max_batches: int | None = None,
    throttle_ms: int | None = None,
    timeout: float | None = None,
) -> Predictions:
    """Predict MS2 intensities for a crosslinked peptide pair.

    Args:
        model: A crosslink model name, e.g. ``"Prosit_2023_intensity_XL_CMS2"``.
        pairs: Mappings with ``alpha_sequence`` and, for most models, ``beta_sequence``.
        precursor_charge: Applied to every pair that does not carry its own.
        collision_energy: Applied to every pair that does not carry its own.
        out: Write a tab-separated table here.
        max_batches: In-flight requests.
        throttle_ms: Delay between request chunks.
        timeout: Seconds to allow.

    Returns:
        A :class:`Predictions`, with ragged fragment arrays as for :func:`fragments`.

    Warning:
        **This family takes a different sequence language from every other function here.** The
        others accept mzLib's ``FullSequence`` notation and convert it; the crosslink models reject
        it and require raw UNIMOD brackets - ``K[UNIMOD:1896]``. That is mzLib's constraint, not a
        choice made here, and it is repeated in :attr:`Predictions.caveats`.
    """
    rows = _as_rows(
        pairs, "alpha_sequence",
        precursor_charge=precursor_charge,
        collision_energy=collision_energy,
    )
    return _predict(
        "crosslink-fragments", model, rows,
        ["alpha_sequence", "beta_sequence", "precursor_charge", "collision_energy"],
        out=out, max_batches=max_batches, throttle_ms=throttle_ms, timeout=timeout,
    )
