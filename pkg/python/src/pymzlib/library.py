"""Read an mzLib spectral library (``.msl``): every precursor, and its fragment ions on request.

A spectral library holds, for each precursor - a peptide or proteoform at one charge - the m/z,
retention time and ion mobility a search expects to see it at, and the fragment ions it expects
in the MS2 spectrum. MetaMorpheus and mzLib write them in mzLib's binary ``.msl`` format.
:func:`read` returns one with mzLib's own reader, ``MslLibrary``::

    >>> import pymzlib
    >>> library = pymzlib.library.read("myPrositLib.msl")
    >>> library.record_count, library.target_count, library.decoy_count
    (5, 5, 0)
    >>> library.columns["full_sequence"][3]
    'QSQHM[Common Variable:Oxidation on M]TEVVR'
    >>> library.columns["retention_time_unit"][:2]
    ['iRT', 'iRT']

One row per precursor, in the library's stored order. With ``fragments=True`` the table is one row
per **fragment ion** instead, the precursor's columns repeated on each row, ready to filter by
product type or plot as a spectrum.

**Three things to know before you use the numbers:**

- **Retention time has a unit per row.** A library stores either an indexed retention time (iRT)
  or run minutes for each precursor (``MslLibraryEntry.RtIsCalibrated``), and one library can hold
  both. Read ``retention_time_unit`` on each row; there is no conversion between them without a
  calibration against a run, which :func:`calibrate_rt` fits.
- **"Not available" is ``None``, not a number.** mzLib stores an ion mobility of 0.0, a q-value of
  NaN and an NCE of 0 when a library does not have one, and documents each as "not available".
  They arrive as ``None``. A text library converted to ``.msl`` (from ``.msp``) has none of them.
- **The library is unfiltered.** Decoys (``D``, ``ED``) and entrapment precursors (``ET``,
  ``ED``) are rows; ``decoy_contaminant_target`` is mzLib's label for each.

Large libraries open **index-only**: precursors in memory, fragment blocks left on disk. mzLib
decides by file size (``MslFileTypeHandler.Open``, 1 GiB at mzLib 1.0.595); pass
``index_only_threshold`` to move the line. The rows are the same either way. Many libraries read
into one table with :func:`read_many`.

:func:`calibrate_rt` is the step between a library and a run. From anchors - precursors you
identified confidently, each with its observed retention time in minutes and its library iRT - it
fits mzLib's ``IrtCalibration`` for each run and converts the times you send with them: iRT to
minutes (where should this precursor elute?) or minutes to iRT (put this run on the library's
scale). The spectral-library guide works through all of it.
"""

from __future__ import annotations

import math
import numbers
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Iterable, Union

from . import _bridge
from .readers import ReadError, WrittenTable, _Table, _batch, _output, _window

__all__ = [
    "LibraryRecords",
    "LibraryFileReport",
    "LibraryBatch",
    "RtCalibration",
    "RtCalibrationGroup",
    "FitFailure",
    "CALIBRATION_KINDS",
    "read",
    "read_many",
    "calibrate_rt",
]

#: A list of paths for :func:`read_many`: each a ``str`` or ``pathlib.Path``.
PathList = Iterable[Union[str, "os.PathLike[str]"]]

#: The first pyMzLib whose bridge has ``library read``: its spec's ``since.pymzlib``.
_LIBRARY_SINCE = "0.6.0"

#: The first pyMzLib whose bridge has ``library calibrate-rt``: its spec's ``since.pymzlib``.
_CALIBRATE_SINCE = "0.7.0"


def __dir__() -> list[str]:
    """The public API only, not this module's own imports."""
    return sorted(__all__)


@dataclass(frozen=True)
class LibraryRecords(_Table):
    """An mzLib spectral library read into one table: what :func:`read` returns.

    Without fragments, one row per precursor. Its columns, under mzLib's ``MslLibraryEntry`` names:
    ``precursor_index`` (the 0-based position in the library, which is also the row order - not
    m/z order), ``full_sequence`` (mzLib bracket notation), ``base_sequence``, ``precursor_mz``
    (m/z), ``charge_state``, ``retention_time`` and ``retention_time_unit`` (``"iRT"`` or
    ``"minutes"``, per row), ``ion_mobility`` (1/K0, V s/cm2; ``None`` where mzLib stores 0.0, "not
    available"), ``is_decoy``, ``is_entrapment``, ``decoy_contaminant_target`` (``"T"``, ``"D"``,
    ``"ET"`` or ``"ED"``), ``protein_accession``, ``protein_name`` and ``gene_name`` (``None``
    where mzLib stores an empty string), ``is_proteotypic``, ``q_value`` (``None`` where mzLib
    stores NaN: a predicted library), ``elution_group_id``, ``molecule_type``, ``source_type``
    (``"Predicted"``, ``"Empirical"`` or ``"EmpiricalRefined"``), ``dissociation_type`` and ``nce``
    (``None`` where mzLib stores 0).

    With ``fragments=True``, one row per fragment ion: those columns, then ``fragment_count``,
    ``fragment_mz``, ``fragment_intensity`` (relative within the precursor, 0 to 1),
    ``product_type``, ``secondary_product_type`` and ``secondary_fragment_number`` (``None`` for a
    terminal ion), ``fragment_number``, ``residue_position`` (``None`` for a custom-loss fragment,
    whose stored position mzLib reads back as 0), ``fragment_charge``, ``neutral_loss`` (Da; 0 is
    no loss), ``is_internal``, ``is_diagnostic``, ``exclude_from_quant`` and ``annotation``. A
    precursor with no fragment ions is one row with those cells ``None``.

    Attributes:
        path: The absolute path that was read.
        file_type: Always ``None``: ``.msl`` is not one of mzLib's ``SupportedFileType`` members,
            so there is no type name to give. ``reader`` names the class that read it.
        reader: ``"MslLibrary"``.
        format_version: The file's ``.msl`` format version, 1 to 3 at mzLib 1.0.595.
        is_compressed: Whether the fragment section is zstd-compressed. A compressed library
            always loads whole.
        is_index_only: Whether mzLib opened it index-only: at or above ``index_only_threshold``
            and not compressed.
        target_count: Precursors that are not decoys, entrapment targets included.
        decoy_count: Decoy precursors, entrapment decoys included.
        proteoform_count: Precursors whose ``molecule_type`` is ``"Proteoform"``.
        min_precursor_mz: The smallest precursor m/z. ``None``: the library has no precursors.
        max_precursor_mz: The largest precursor m/z. ``None``: the library has no precursors.
        record_count: Precursors in the whole library, before ``offset`` and ``limit``.
        rows_not_read: Always ``None``: rows are not counted for a binary library - mzLib checks the
            whole file before returning any precursor, so every one is read or the call fails.
        retention_time_unit: Always ``None``: the unit varies by precursor, so it is the
            ``retention_time_unit`` column.
        caveats: What this table cannot be trusted to mean, each citing the mzLib source.
        column_names: The columns, in order.
        absent_fields: Columns this library has no source for, so every value is ``None``:
            ``ion_mobility`` when no precursor has one, ``protein_accession`` and ``protein_name``
            when there is no protein table, ``gene_name`` when there are no gene names.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows.
        excluded_fields: Always empty: every field crosses, the fragments with ``fragments=True``.
        error: Always ``None``: a library that cannot be read raises instead.
        fragments_included: Whether the rows are fragment ions.
        index_only_threshold: The size, in bytes, at and above which mzLib opened the library
            index-only: the one you passed, or mzLib's default (1073741824, 1 GiB).
        returned_count: Precursors returned in :attr:`columns` - precursors, not rows; the unit
            ``limit`` and ``offset`` count in. Zero when ``out`` was given.
        row_count: Rows in :attr:`columns`: ``returned_count`` without fragments, one per fragment
            ion with them. Zero when ``out`` was given.
        offset: Precursors skipped.
        truncated: Whether precursors were left behind, by either ``limit`` or ``offset``.
        columns: Column name -> list of values, one per row - the shape ``pandas.DataFrame``
            accepts. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
    """

    path: str
    reader: Any
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    row_count: int = 0
    file_type: Any = None
    format_version: int | None = None
    is_compressed: bool = False
    is_index_only: bool = False
    target_count: int = 0
    decoy_count: int = 0
    proteoform_count: int = 0
    min_precursor_mz: float | None = None
    max_precursor_mz: float | None = None
    rows_not_read: Any = None
    retention_time_unit: Any = None
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: WrittenTable | None = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None
    fragments_included: bool = False
    index_only_threshold: int | None = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "LibraryRecords":
        return cls(
            path=data.get("path", ""),
            reader=data.get("reader"),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            row_count=int(data.get("row_count", 0)),
            output=_output(data),
            fragments_included=bool(data.get("fragments_included", False)),
            index_only_threshold=data.get("index_only_threshold"),
            columns=data.get("columns"),
            **_library_block(data),
        )


def _library_block(data: dict[str, Any]) -> dict[str, Any]:
    """The per-library facts, shared by a single read and each ``files[]`` entry of a batch."""
    return dict(
        file_type=data.get("file_type"),
        format_version=data.get("format_version"),
        # Raw, not coerced: a skipped input's entry carries null for every fact it could not
        # establish, and a 0 or False here would read as a library that was opened and was empty.
        is_compressed=data.get("is_compressed"),
        is_index_only=data.get("is_index_only"),
        target_count=data.get("target_count"),
        decoy_count=data.get("decoy_count"),
        proteoform_count=data.get("proteoform_count"),
        min_precursor_mz=data.get("min_precursor_mz"),
        max_precursor_mz=data.get("max_precursor_mz"),
        rows_not_read=data.get("rows_not_read"),
        retention_time_unit=data.get("retention_time_unit"),
        caveats=list(data.get("caveats") or []),
        column_names=list(data.get("column_names") or []),
        absent_fields=list(data.get("absent_fields") or []),
        failed_fields=list(data.get("failed_fields") or []),
        excluded_fields=list(data.get("excluded_fields") or []),
        error=ReadError._from_wire(data.get("error")),
    )


@dataclass(frozen=True)
class LibraryFileReport:
    """The per-library facts of one input to :func:`read_many`: one entry of :attr:`LibraryBatch.files`.

    The fields a single :func:`read` reports at its top level, under the same names. For an input
    that failed under ``on_error="skip"``, :attr:`error` says why, and the facts that could not be
    established are ``None``, zero or empty.

    Attributes:
        path: The absolute path of the input.
        file_type: Always ``None``: ``.msl`` is not an mzLib ``SupportedFileType``.
        reader: ``"MslLibrary"``, or ``None`` for a failed input.
        format_version: The ``.msl`` format version. ``None``: a failed input.
        is_compressed: Whether the fragment section is compressed. ``None``: a failed input.
        is_index_only: Whether mzLib opened it index-only. ``None``: a failed input.
        target_count: Precursors that are not decoys, entrapment targets included. ``None``: a
            failed input.
        decoy_count: Decoy precursors, entrapment decoys included. ``None``: a failed input.
        proteoform_count: Proteoform precursors. ``None``: a failed input.
        min_precursor_mz: The smallest precursor m/z. ``None``: no precursors, or a failed input.
        max_precursor_mz: The largest precursor m/z. ``None``: no precursors, or a failed input.
        record_count: Precursors in the library. ``None``: a failed input - not zero, because
            nothing was counted.
        rows_not_read: Always ``None`` (see :attr:`LibraryRecords.rows_not_read`).
        retention_time_unit: Always ``None``: the unit is the ``retention_time_unit`` column.
        caveats: What this library's rows cannot be trusted to mean.
        column_names: This library's columns, without ``source_index``/``source_path``.
        absent_fields: Columns this library has no source for.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows.
        excluded_fields: Always empty.
        error: Why this input could not be read, under ``on_error="skip"``; ``None`` when it was.
    """

    path: str
    reader: Any = None
    record_count: int | None = None
    file_type: Any = None
    format_version: int | None = None
    is_compressed: bool | None = None
    is_index_only: bool | None = None
    target_count: int | None = None
    decoy_count: int | None = None
    proteoform_count: int | None = None
    min_precursor_mz: float | None = None
    max_precursor_mz: float | None = None
    rows_not_read: Any = None
    retention_time_unit: Any = None
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @property
    def ok(self) -> bool:
        """Whether this input was read (``error`` is ``None``)."""
        return self.error is None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "LibraryFileReport":
        count = data.get("record_count")
        return cls(
            path=data.get("path", ""),
            reader=data.get("reader"),
            record_count=None if count is None else int(count),
            **_library_block(data),
        )


@dataclass(frozen=True)
class LibraryBatch(_Table):
    """Many libraries read into ONE long table: what :func:`read_many` returns.

    The first two columns say where each row came from: ``source_index``, the library's 0-based
    position in the list you passed, and ``source_path``, its absolute path. The rest are
    :class:`LibraryRecords`' columns. Rows are grouped by library in input order, each library's
    rows in ``precursor_index`` order - **whatever** ``threads`` **was**.

    Attributes:
        file_count: Libraries in the list you passed.
        read_count: Libraries read.
        failed_count: Libraries that could not be read - always 0 unless ``on_error="skip"``.
        on_error: ``"fail"`` or ``"skip"``, as requested.
        record_count: Precursors, summed over the libraries read.
        returned_count: Precursors carried back in :attr:`columns`. Zero when ``out`` was given.
        row_count: Rows in :attr:`columns` - one per fragment ion with ``fragments=True``. Zero when
            ``out`` was given.
        fragments_included: Whether the rows are fragment ions.
        index_only_threshold: The size in bytes at and above which a library opened index-only.
        column_names: ``source_index``, ``source_path``, then the library columns, in order.
        columns: Column name -> list of values, one per row. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
        files: One :class:`LibraryFileReport` per input, in input order.
    """

    file_count: int
    read_count: int
    failed_count: int
    on_error: str
    record_count: int
    returned_count: int
    row_count: int = 0
    fragments_included: bool = False
    index_only_threshold: int | None = None
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: WrittenTable | None = None
    files: list[LibraryFileReport] = field(default_factory=list)

    @property
    def failed_files(self) -> list[LibraryFileReport]:
        """The inputs that could not be read, under ``on_error="skip"``. Empty otherwise."""
        return [report for report in self.files if report.error is not None]

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "LibraryBatch":
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            on_error=data.get("on_error") or "fail",
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            row_count=int(data.get("row_count", 0)),
            fragments_included=bool(data.get("fragments_included", False)),
            index_only_threshold=data.get("index_only_threshold"),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=_output(data),
            files=[LibraryFileReport._from_wire(entry) for entry in (data.get("files") or [])],
        )


def _options(fragments: bool, index_only_threshold: int | None) -> list[str]:
    """``--fragments`` and ``--index-only-threshold``, validated before the bridge is spawned."""
    if not isinstance(fragments, bool):
        raise _bridge.UsageError(f"fragments must be True or False; got {fragments!r}.")
    args = ["--fragments"] if fragments else []
    if index_only_threshold is not None:
        # bool first: isinstance(True, int) is True, so True would otherwise pass as 1 byte.
        if (
            isinstance(index_only_threshold, bool)
            or not isinstance(index_only_threshold, int)
            or index_only_threshold < 0
        ):
            raise _bridge.UsageError(
                "index_only_threshold must be a whole number of bytes, 0 or more, or None for "
                f"mzLib's default; got {index_only_threshold!r}."
            )
        args += ["--index-only-threshold", str(index_only_threshold)]
    return args


def read(
    path: str | os.PathLike[str],
    *,
    fragments: bool = False,
    index_only_threshold: int | None = None,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> LibraryRecords:
    """Read an mzLib spectral library (``.msl``) as one row per precursor, or per fragment ion.

    Wraps mzLib's ``MslFileTypeHandler.Open`` and ``MslLibrary.GetEntry``. Rows come in the
    library's stored order (``precursor_index``), not by m/z; sort on ``precursor_mz`` for that.
    mzLib's "not available" values - an ion mobility of 0.0, a NaN q-value, an NCE of 0, an empty
    protein accession - arrive as ``None``.

    Args:
        path: Path to a ``.msl`` file. Any other extension is a ``UsageError``: text libraries
            (``.msp``) are not read here.
        fragments: One row per **fragment ion** instead of one per precursor, with the precursor's
            columns repeated on each row. Off by default: a library has tens of fragments per
            precursor, so pair it with ``limit`` or ``out`` for a large one.
        index_only_threshold: A library this many **bytes** or larger opens index-only - precursors
            in memory, fragment blocks read from disk as needed; a smaller one loads whole. ``None``
            (the default) sends nothing, so mzLib's own default applies:
            ``MslFileTypeHandler.DefaultIndexOnlyThresholdBytes``, 1073741824 bytes (1 GiB) at
            mzLib 1.0.595. ``0`` opens every library index-only. The rows are the same either way.
        limit: Maximum **precursors** to return; with ``fragments`` each gives one row per fragment
            ion. ``None`` returns all of them.
        offset: **Precursors** to skip, in ``precursor_index`` order.
        out: Write a tab-separated table here and return only a summary.
            The path must end in ``.tsv`` (any case); any other is a ``UsageError`` before a read.
        timeout: Seconds to allow. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`LibraryRecords`.

    Raises:
        UsageError: the path is blank, missing, a directory or not ``.msl``; ``fragments`` is not
            a bool; ``index_only_threshold``, ``limit`` or ``offset`` is not a whole number of 0 or
            more; ``out`` is not a ``.tsv`` or names the input; or the bridge in use predates this
            function.
        BridgeError: mzLib could not read the library - ``error_type`` ``"FormatException"`` for
            a file that is not a valid ``.msl``, ``"InvalidDataException"`` for a checksum
            mismatch (a damaged file).

    Examples:
        >>> library = read("synthetic_library.msl")
        >>> library.columns["decoy_contaminant_target"]
        ['T', 'D', 'ET', 'ED', 'T', 'T']
        >>> library.columns["retention_time_unit"][:3]
        ['iRT', 'iRT', 'minutes']
        >>> library.columns["q_value"][4], library.columns["ion_mobility"][4]   # predicted: not available
        (None, None)

        One precursor's fragment ions:

        >>> ions = read("synthetic_library.msl", fragments=True, limit=1)
        >>> ions.returned_count, ions.row_count
        (1, 6)
        >>> ions.columns["annotation"]
        ['y1', 'b2', 'y2', 'bIy[2-5]', 'y3-18.0106', 'y4-43.0058']
    """
    # Every argument is checked before any process is started, the version probe included.
    extra = _options(fragments, index_only_threshold)
    args = _window("read", path, limit=limit, offset=offset, out=out, module="library")
    _bridge.require_verb("library read", since=_LIBRARY_SINCE)
    data = _bridge.invoke(*args, *extra, timeout=timeout)
    return LibraryRecords._from_wire(data)


def read_many(
    paths: PathList,
    *,
    fragments: bool = False,
    index_only_threshold: int | None = None,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> LibraryBatch:
    """Read many spectral libraries into one long table.

    :func:`read` for a list - a predicted library and the empirical one refined from it, say. Each
    library's facts (counts, load mode, format version) are in its :attr:`LibraryFileReport`.

    Args:
        paths: The ``.msl`` files to read, as a list (or any iterable) of ``str`` or
            ``pathlib.Path``. Read by ONE bridge process into ONE table, in this order; a file's
            position is its ``source_index``. A repeated path is refused rather than read twice.
        fragments: One row per fragment ion, in every library. Off by default.
        index_only_threshold: The size in **bytes** at and above which a library opens
            index-only, for every library. ``None`` (the default) leaves it to mzLib (1 GiB).
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`LibraryBatch.output`), one library at a time, in order. A batch that stops on an
            error removes its partial table.
            The path must end in ``.tsv`` (any case); any other is a ``UsageError`` before a read.
        threads: Libraries to read at once: ``1`` (the default), more, or ``-1`` for one per core.
            The result is **byte-identical at any value**, so this trades memory for speed and never
            changes an answer. The default is 1 because a library loaded whole is held in memory
            while its rows are written.
        on_error: ``"fail"`` (the default) stops at the first library that cannot be read and
            raises, naming it; ``"skip"`` records the failure in that library's
            :attr:`LibraryFileReport.error` and reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`LibraryBatch` with the columns of :class:`LibraryRecords`, after
        ``source_index`` and ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``fragments``, ``index_only_threshold``, ``threads`` or ``on_error`` is not a
            value listed above; or, with ``on_error="fail"``, a file is missing or is not ``.msl`` -
            the message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib failed to read a library; the message names
            the input the same way.

    Examples:
        >>> batch = read_many(["synthetic_library.msl", "myPrositLib.msl"])
        >>> batch.record_count, [f.record_count for f in batch.files]
        (11, [6, 5])
        >>> batch.files[1].absent_fields
        ['ion_mobility', 'protein_accession', 'protein_name', 'gene_name']
    """
    extra = _options(fragments, index_only_threshold)
    args, stdin = _batch("read", paths, threads=threads, on_error=on_error, out=out, module="library")
    _bridge.require_verb("library read", since=_LIBRARY_SINCE)
    data = _bridge.invoke(*args, *extra, stdin=stdin, timeout=timeout)
    return LibraryBatch._from_wire(data)


@dataclass(frozen=True)
class FitFailure:
    """Why mzLib did not fit one group of :func:`calibrate_rt`, under ``on_error="skip"``.

    Attributes:
        kind: Always ``"usage"``: mzLib refused the anchors you sent for this group.
        type: The mzLib exception name, ``"ArgumentException"``.
        message: mzLib's own message, e.g. ``"Calibration needs at least 20 anchors; got 14."``.
    """

    kind: str
    type: str
    message: str

    @classmethod
    def _from_wire(cls, payload: Mapping[str, Any] | None) -> "FitFailure | None":
        if not payload:
            return None
        return cls(
            kind=str(payload.get("kind", "")),
            type=str(payload.get("type", "")),
            message=str(payload.get("message", "")),
        )


@dataclass(frozen=True)
class RtCalibrationGroup:
    """One group's (one run's) calibration: one entry of :attr:`RtCalibration.groups`.

    Attributes:
        group: The group name, as you sent it.
        kind: ``"linear"`` or ``"lowess"``, the shape mzLib fitted (``IrtCalibrationModel.Kind``).
            ``None``: not fitted.
        anchor_count: Anchors this group was fitted from (anchors); the anchors you sent, if mzLib
            refused them.
        query_count: Times converted for this group (queries).
        residual_sd: How far the anchors that fit sit from the curve, in **iRT units**: mzLib's
            robust residual SD (``IrtCalibrationModel.ResidualSd``), wrong anchors set aside.
            Divide by the curve's slope (iRT per minute) to read it in minutes. ``None``: not
            fitted.
        first_anchor_minutes: The earliest anchor's time, in **minutes**: the map is fitted from
            here. ``None``: not fitted.
        last_anchor_minutes: The latest anchor's time, in **minutes**: the map is fitted up to
            here. ``None``: not fitted.
        error: Why mzLib refused this group's anchors, under ``on_error="skip"``; ``None`` when it
            was fitted.
    """

    group: str
    kind: str | None = None
    anchor_count: int = 0
    query_count: int = 0
    residual_sd: float | None = None
    first_anchor_minutes: float | None = None
    last_anchor_minutes: float | None = None
    error: FitFailure | None = None

    @property
    def ok(self) -> bool:
        """Whether mzLib fitted this group (``error`` is ``None``)."""
        return self.error is None

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "RtCalibrationGroup":
        return cls(
            group=str(data.get("group", "")),
            kind=data.get("kind"),
            anchor_count=int(data.get("anchor_count", 0)),
            query_count=int(data.get("query_count", 0)),
            residual_sd=data.get("residual_sd"),
            first_anchor_minutes=data.get("first_anchor_minutes"),
            last_anchor_minutes=data.get("last_anchor_minutes"),
            error=FitFailure._from_wire(data.get("error")),
        )


@dataclass(frozen=True)
class RtCalibration(_Table):
    """What :func:`calibrate_rt` returns: each run's calibration, and every time you sent, converted.

    One row per anchor, then one per query, in the order you sent them:

    ===================  ======  ======================================================================
    column               unit    meaning; ``None`` means
    ===================  ======  ======================================================================
    ``group``            -       as you sent it
    ``role``             -       ``"anchor"``; ``"to_irt"`` (you gave minutes); ``"to_minutes"`` (you
                                 gave iRT)
    ``rt_minutes``       min     as given, or ``ToRtMinutes(irt)`` for a ``to_minutes`` query;
                                 ``None`` = the group was not fitted and you gave none
    ``irt``              iRT     as given, or ``ToIrt(rt_minutes)`` for a ``to_irt`` query; ``None``
                                 = the group was not fitted and you gave none
    ``fitted_irt``       iRT     anchors only: the curve's iRT at the anchor's time, so
                                 ``irt - fitted_irt`` is its residual; ``None`` on a query, or the
                                 group was not fitted
    ``is_extrapolated``  -       ``rt_minutes`` lies outside the group's anchors, where the map is a
                                 straight-line continuation; ``None`` = the group was not fitted
    ===================  ======  ======================================================================

    Attributes:
        kind: ``"linear"`` or ``"lowess"``: the kind applied, yours or mzLib's default.
        minimum_anchors: The fewest anchors a group may be fitted from (anchors), yours or mzLib's.
        bandwidth: The fraction of anchors in each local LOWESS fit, in (0, 1], yours or mzLib's.
        knots: The points along the run the LOWESS curve is fitted at (points), yours or mzLib's.
        robustness_iterations: Rounds of bisquare reweighting (iterations), yours or mzLib's.
        on_error: ``"fail"`` or ``"skip"``, as requested.
        group_count: Groups (runs) sent, over anchors and queries (groups).
        fitted_count: Groups mzLib fitted (groups).
        failed_count: Groups mzLib refused (groups); non-zero only under ``on_error="skip"``.
        anchor_count: Anchors sent, over every group (anchors).
        query_count: Queries sent, over every group (queries).
        row_count: Rows in :attr:`columns` (rows): ``anchor_count + query_count``.
        groups: One :class:`RtCalibrationGroup` per group, in order of first appearance.
        column_names: The table's columns, in order.
        columns: Column name -> one value per row, the shape ``pandas.DataFrame`` accepts.
        caveats: What a calibration can and cannot tell you, plus any that apply to this call
            (queries outside the anchors, groups not fitted).
    """

    kind: str
    minimum_anchors: int
    bandwidth: float
    knots: int
    robustness_iterations: int
    on_error: str
    group_count: int
    fitted_count: int
    failed_count: int
    anchor_count: int
    query_count: int
    row_count: int
    groups: list[RtCalibrationGroup] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: dict[str, list[Any]] = field(default_factory=dict)
    caveats: list[str] = field(default_factory=list)

    def group(self, name: str) -> RtCalibrationGroup:
        """The calibration of the group called ``name``.

        Raises:
            KeyError: no group has that name.
        """
        for entry in self.groups:
            if entry.group == name:
                return entry
        raise KeyError(name)

    @property
    def failed_groups(self) -> list[RtCalibrationGroup]:
        """The groups mzLib refused, under ``on_error="skip"``. Empty otherwise."""
        return [entry for entry in self.groups if entry.error is not None]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "RtCalibration":
        return cls(
            kind=str(data.get("kind", "")),
            minimum_anchors=int(data.get("minimum_anchors", 0)),
            bandwidth=float(data.get("bandwidth", 0.0)),
            knots=int(data.get("knots", 0)),
            robustness_iterations=int(data.get("robustness_iterations", 0)),
            on_error=str(data.get("on_error") or "fail"),
            group_count=int(data.get("group_count", 0)),
            fitted_count=int(data.get("fitted_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            anchor_count=int(data.get("anchor_count", 0)),
            query_count=int(data.get("query_count", 0)),
            row_count=int(data.get("row_count", 0)),
            groups=[RtCalibrationGroup._from_wire(g) for g in (data.get("groups") or [])],
            column_names=list(data.get("column_names") or []),
            columns=dict(data.get("columns") or {}),
            caveats=list(data.get("caveats") or []),
        )


#: The values of ``kind`` in :func:`calibrate_rt`: mzLib's ``IrtCalibrationKind``.
CALIBRATION_KINDS = ("linear", "lowess")


def _calibration_lines(entries: Any, what: str, query: bool) -> list[str]:
    """One stdin line per anchor or query: group, rt_minutes, irt, tab-separated."""
    if entries is None:
        return []
    if isinstance(entries, (str, bytes, Mapping)) or not isinstance(entries, Iterable):
        raise _bridge.UsageError(
            f"{what} must be a list of (group, rt_minutes, irt) entries or mappings; got {entries!r}."
        )
    lines = []
    for i, entry in enumerate(entries):
        if isinstance(entry, Mapping):
            if "group" not in entry:
                raise _bridge.UsageError(
                    f"{what}[{i}] needs a 'group' key; got {sorted(map(str, entry))}."
                )
            values = [entry["group"], entry.get("rt_minutes"), entry.get("irt")]
        elif isinstance(entry, Sequence) and not isinstance(entry, (str, bytes)) and len(entry) == 3:
            values = list(entry)
        else:
            raise _bridge.UsageError(f"{what}[{i}] must be (group, rt_minutes, irt); got {entry!r}.")
        group = values[0]
        if not isinstance(group, str) or not group.strip():
            raise _bridge.UsageError(f"{what}[{i}] needs a non-empty group name; got {group!r}.")
        if any(c in group for c in "\t\n\r"):
            raise _bridge.UsageError(f"{what}[{i}]'s group name holds a tab or newline: {group!r}.")
        if group.strip() == "--":
            raise _bridge.UsageError(
                f"{what}[{i}]'s group name is '--', which separates anchors from queries on the wire."
            )
        cells = [
            _time_cell(value, f"{what}[{i}] {name}")
            for value, name in zip(values[1:], ("rt_minutes", "irt"))
        ]
        if query and (cells[0] == "") == (cells[1] == ""):
            raise _bridge.UsageError(
                f"{what}[{i}] must give exactly one of rt_minutes and irt, the other None: the "
                f"missing one is what the calibration computes; got {entry!r}."
            )
        lines.append("\t".join([group, *cells]))
    return lines


def _time_cell(value: object, what: str) -> str:
    """A time on a stdin line: ``None`` and NaN are blank, a number its repr."""
    if value is None:
        return ""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise _bridge.UsageError(f"{what} must be a number or None; got {value!r}.")
    number = float(value)
    return "" if math.isnan(number) else repr(number)


def _option(name: str, value: Any, *, integer: bool) -> list[str]:
    """``--name value`` when given; nothing when ``None``, so mzLib's default applies."""
    if value is None:
        return []
    wanted = numbers.Integral if integer else numbers.Real
    if isinstance(value, bool) or not isinstance(value, wanted):
        raise _bridge.UsageError(
            f"{name.replace('-', '_')} must be {'a whole number' if integer else 'a number'}, or "
            f"None for mzLib's default; got {value!r}."
        )
    return [f"--{name}", str(int(value)) if integer else repr(float(value))]


def calibrate_rt(
    anchors: Iterable[Sequence[Any] | Mapping[str, Any]],
    queries: Iterable[Sequence[Any] | Mapping[str, Any]] | None = None,
    *,
    kind: str | None = None,
    minimum_anchors: int | None = None,
    bandwidth: float | None = None,
    knots: int | None = None,
    robustness_iterations: int | None = None,
    on_error: str = "fail",
    timeout: float | None = 60,
) -> RtCalibration:
    """Fit each run's map from minutes onto a library's iRT scale, and convert times either way.

    A spectral library stores retention time on a run-independent scale, iRT (Escher et al. 2012).
    A run reports minutes. mzLib's ``IrtCalibration.Fit`` maps one onto the other from **anchors**:
    precursors you identified confidently in the run, each with its observed apex time in minutes
    and its library iRT. The fit is robust: it starts from binned medians and downweights anchors
    far from the curve (Tukey bisquare), so some wrong identifications do not bend it. It is
    strictly increasing, so it inverts. By default it is a LOWESS curve (Cleveland 1979) that
    follows a gradient that bends; ``kind="linear"`` fits one straight line.

    Each **group** (a run, typically) is fitted on its own anchors, and the **queries** you send for
    it are converted with that fit: give a time in minutes to get its iRT, or an iRT to get the
    minutes at which that precursor should elute in this run. mzLib's fitted curve cannot be
    returned, so send the queries in the same call as their anchors.

    Choosing anchors is yours: mzLib applies no q-value cutoff, and pyMzLib adds none.

    Args:
        anchors: One entry per anchor: ``(group, rt_minutes, irt)``, or a mapping with ``group``,
            ``rt_minutes`` and ``irt``. ``group`` is any name (the run); ``rt_minutes`` is the
            observed retention time in **minutes**; ``irt`` is the library's indexed retention time
            (**iRT**). mzLib refuses a ``None`` or NaN time.
        queries: Times to convert, one entry each: ``(group, rt_minutes, None)`` converts minutes to
            iRT, ``(group, None, irt)`` converts iRT to minutes; or a mapping with ``group`` and one
            of ``rt_minutes`` or ``irt``. Each group's queries use that group's fit. ``None`` (the
            default) converts nothing: you get each group's fit and its anchors' residuals.
        kind: ``"linear"`` or ``"lowess"`` (:data:`CALIBRATION_KINDS`, mzLib's
            ``IrtCalibrationKind``). ``None`` sends nothing, so mzLib's default, LOWESS, applies.
        minimum_anchors: A group with fewer anchors than this is refused, not fitted (anchors).
            ``None`` leaves mzLib's default, 20. mzLib requires at least 3.
        bandwidth: LOWESS only: the **fraction of anchors** in each local fit, in (0, 1]. Smaller
            follows the gradient, and its noise, more closely. ``None`` leaves mzLib's default,
            0.3.
        knots: LOWESS only: how many evenly spaced times between the first and last anchor the
            curve is fitted at (points). ``None`` leaves mzLib's default, 100. At least 3.
        robustness_iterations: Rounds of bisquare reweighting after the first fit (iterations).
            ``None`` leaves mzLib's default, 3. 0 or more.
        on_error: ``"fail"`` (the default) raises at the first group mzLib refuses, naming it;
            ``"skip"`` records the refusal in that group's :attr:`RtCalibrationGroup.error` and fits
            the rest.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`RtCalibration`: each group's fit in :attr:`RtCalibration.groups`, and one row
        per anchor and query.

    Raises:
        UsageError: no anchors; an entry that is not three values or a mapping with ``group``; an
            empty group name, or one holding a tab or newline; a time that is not a number or
            ``None``; a query with both times or neither; an option that is not a value listed
            above; or, from mzLib, an option out of range, or (with ``on_error="fail"``) a group it
            refuses: fewer anchors than ``minimum_anchors``, a ``None`` or NaN time, every anchor
            at one time, or iRT that does not increase with time.

    Examples:
        Three 24-minute DIA runs of PBMCs, from a DIA-NN report in mzLib's test data. One
        precursor, ``AAAAAAAAGAFAGR`` at charge 2, is held out of the anchors, and each run predicts
        where it should elute from its library iRT:

        >>> import csv
        >>> with open("diann_pbmc_rt.tsv", newline="") as fh:
        ...     rows = list(csv.DictReader(fh, delimiter="\\t"))
        >>> held_out = "AAAAAAAAGAFAGR2"
        >>> anchors = [(r["run"], float(r["retention_time"]), float(r["indexed_retention_time"]))
        ...            for r in rows if r["precursor_id"] != held_out]
        >>> observed = {r["run"]: float(r["retention_time"])
        ...             for r in rows if r["precursor_id"] == held_out}
        >>> library_irt = 11.886096
        >>> queries = [(run, None, library_irt) for run in observed]
        >>> calibration = pymzlib.library.calibrate_rt(anchors, queries, kind="linear",
        ...                                            minimum_anchors=10)
        >>> for row in calibration.records:
        ...     if row["role"] == "to_minutes":
        ...         print(row["group"][:8], f"{row['rt_minutes']:.2f}", f"{observed[row['group']]:.2f}")
        20240316 9.22 9.21
        20240317 9.45 9.23
        20240315 9.17 8.94

        Each prediction is within a quarter of a minute of where the precursor was seen.
    """
    if kind is not None and kind not in CALIBRATION_KINDS:
        raise _bridge.UsageError(
            f"kind must be one of {', '.join(CALIBRATION_KINDS)} (mzLib's IrtCalibrationKind), or "
            f"None for mzLib's default; got {kind!r}."
        )
    if on_error not in ("fail", "skip"):
        raise _bridge.UsageError(f"on_error must be 'fail' or 'skip'; got {on_error!r}.")
    args = ["library", "calibrate-rt"]
    if kind is not None:
        args += ["--kind", kind]
    args += _option("minimum-anchors", minimum_anchors, integer=True)
    args += _option("bandwidth", bandwidth, integer=False)
    args += _option("knots", knots, integer=True)
    args += _option("robustness-iterations", robustness_iterations, integer=True)
    if on_error == "skip":
        args += ["--on-error", "skip"]
    anchor_lines = _calibration_lines(anchors, "anchors", query=False)
    if not anchor_lines:
        raise _bridge.UsageError("calibrate_rt() needs at least one anchor.")
    query_lines = _calibration_lines(queries, "queries", query=True)
    stdin = "\n".join(anchor_lines + (["--", *query_lines] if query_lines else [])) + "\n"
    _bridge.require_verb("library calibrate-rt", since=_CALIBRATE_SINCE)
    data = _bridge.invoke(*args, stdin=stdin, timeout=timeout)
    return RtCalibration._from_wire(data)
