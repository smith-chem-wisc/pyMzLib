"""Read proteomics result files: what a file *is*, what you can do with it, and its records.

mzLib recognises 29 file types written by a dozen different search and deconvolution tools  - 
MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur,
FlashLFQ - and dispatches each to a parser it maintains. This module asks it what a path is::

    >>> import pymzlib
    >>> info = pymzlib.readers.identify("psm.tsv")     # doctest: +SKIP
    >>> info.file_type, info.views                     # doctest: +SKIP
    ('MsFraggerPsm', ['quantifiable'])

**Read** :attr:`FileInfo.views` **before assuming anything.** It is tempting to describe mzLib as
reading 29 formats into one uniform shape, and it does not. The formats fall into disjoint families
and several belong to no family at all:

- ``"quantifiable"`` - a cross-format record view (sequence, retention time, charge, mass, protein
  groups), and the input :func:`pymzlib.flashlfq.quantify` accepts. **Exactly three file types have
  it**: MetaMorpheus ``.psmtsv`` and ``.osmtsv``, and MSFragger ``psm.tsv``. This view reports what
  mzLib's *interface* offers, not that the numbers are comparable - see the warning below, and do
  not quantify the MSFragger one.
- ``"ms1_features"`` - deconvolved MS1 features (TopFD ``_ms1.feature``, Dinosaur).
- ``"spectra"`` - the file is spectra, not results (``.raw``, ``.mzML``, ``.mgf``, ``.d``, msalign).
- ``"spectral_match"`` - records are identifications but share no file-level interface
  (MsPathFinderT, Casanovo).
- ``[]`` - **an empty list is a real and common answer.** TopPIC, Crux, MSFragger's peptide/protein
  tables and the FlashDeconv formats each parse into their own record type with nothing in common.
  mzLib reads them; there is simply no uniform view to project them onto.

Call :func:`formats` for the whole table. It is enumerated from mzLib rather than transcribed, so it
cannot drift from what mzLib actually dispatches.

**Two things this module deliberately does not tell you**, in the "surface it, don't hide it" spirit
of the rest of pyMzLib:

- **Which tool wrote the file.** mzLib has a ``Software`` property that looks like the answer and is
  not: readers carry their software constant on a constructor that mzLib's own file factory does not
  use, so the value is unset for everything the factory returns - and it is not reliably set on the
  other constructor either. Rather than reconstruct a plausible answer, there is no ``software``
  field. :attr:`FileInfo.file_type` already names the tool.
- **Whether the numbers inside mean the same thing across formats.** They do not, and this is the
  trap most likely to produce a wrong result. mzLib's result-file readers pass through whatever the
  tool wrote, with no unit conversion: MetaMorpheus retention times are in **minutes**, MSFragger's
  and TopPIC's are in **seconds**, and TopFD changed from seconds to minutes between v1.6.2 and
  v1.7.0 *within the same file type*. Likewise ``is_decoy`` is hardcoded ``False`` for MSFragger,
  which means "mzLib cannot tell" rather than "target" - MSFragger's ``psm.tsv`` carries no
  target/decoy column at all. ``monoisotopic_mass`` is the *theoretical* peptide mass in **both**
  formats, never the observed precursor mass. Identifying a file is safe; comparing raw fields
  across formats is not.

- **Anything about confidence.** There is no q-value, PEP or score in this view, because
  ``IQuantifiableRecord`` carries only what FlashLFQ needs. **Nothing you get back is
  FDR-filtered**, even though every one of these files records confidence somewhere. Filter before
  you report.

.. warning::

   That units mismatch is not hypothetical. Passing an MSFragger ``psm.tsv`` to
   :func:`pymzlib.flashlfq.quantify` returns near-zero intensities, because FlashLFQ reads the
   seconds as minutes and searches for each peptide roughly sixty times too early in the gradient.
   ``identify()`` will still call the file ``quantifiable`` - that is mzLib's interface, honestly
   reported - but quantify MetaMorpheus output only until the upstream fix lands.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import _bridge

__all__ = [
    "Format",
    "FileInfo",
    "WrittenTable",
    "ResultRecords",
    "formats",
    "identify",
    "read_results",
]

#: The view name for the cross-format record shape :func:`pymzlib.flashlfq.quantify` consumes.
QUANTIFIABLE = "quantifiable"


def __dir__() -> list[str]:
    """The public API only.

    Without this, ``dir(pymzlib.readers)`` lists this module's own imports - ``Any``,
    ``annotations``, ``dataclass``, ``field`` - alongside the functions, which is four false leads
    out of eleven names on a package whose discovery story is ``dir()`` and ``help()``.
    """
    return sorted(__all__)


@dataclass(frozen=True)
class Format:
    """One file type mzLib can recognise.

    Attributes:
        file_type: mzLib's ``SupportedFileType`` name, e.g. ``"MsFraggerPsm"``, ``"psmtsv"``.
        extension: The extension or filename suffix mzLib dispatches on, e.g. ``"psm.tsv"``,
            ``"_ms1.feature"``. **Not unique across file types** - ``BrukerD`` and
            ``BrukerTimsTof`` are both ``.d`` (told apart by the directory's contents), and
            several formats share ``.tsv``.
        reader: The name of the mzLib class that parses it, for cross-referencing the mzLib source.
        views: The uniform views this format supports - see the module docstring. Often empty.
    """

    file_type: str
    extension: Any
    reader: Any
    views: list[str] = field(default_factory=list)

    @property
    def is_quantifiable(self) -> bool:
        """Whether this format offers the cross-format record view (and so feeds FlashLFQ)."""
        return QUANTIFIABLE in self.views

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "Format":
        return cls(
            file_type=payload.get("file_type", ""),
            extension=payload.get("extension"),
            reader=payload.get("reader"),
            views=list(payload.get("views") or []),
        )


@dataclass(frozen=True)
class FileInfo:
    """What a particular file is, and what can be done with it.

    Attributes:
        path: The absolute path that was identified.
        file_type: mzLib's ``SupportedFileType`` name.
        extension: The extension mzLib dispatched on.
        reader: The mzLib class that would parse it.
        views: The uniform views this file supports - see the module docstring. Often empty, which
            means mzLib can read the file but offers no cross-format projection of it.
    """

    path: str
    file_type: str
    extension: Any
    reader: Any
    views: list[str] = field(default_factory=list)

    @property
    def is_quantifiable(self) -> bool:
        """Whether this file offers the cross-format record view.

        When ``True``, the path can be passed straight to :func:`pymzlib.flashlfq.quantify` as
        ``psms``. When ``False``, mzLib can still read the file - it simply has no uniform view,
        so quantification would fail on it.
        """
        return QUANTIFIABLE in self.views

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "FileInfo":
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            extension=data.get("extension"),
            reader=data.get("reader"),
            views=list(data.get("views") or []),
        )


@dataclass(frozen=True)
class WrittenTable:
    """Where :func:`read_results` wrote a table, when asked to write one instead of returning it.

    Attributes:
        path: The absolute path written.
        format: Always ``"tsv"``. **Tab-separated, not comma-separated**, because these fields
            contain commas - MSFragger's mapped proteins are a comma-separated list inside a single
            field, and joined accessions look the same. It is also what every mzLib reader and
            writer uses. Read it with ``csv.reader(f, delimiter="\\t")`` or
            ``pandas.read_csv(path, sep="\\t")``.
        row_count: Rows written, excluding the header.
    """

    path: str
    format: str
    row_count: int

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "WrittenTable":
        return cls(
            path=payload.get("path", ""),
            format=payload.get("format", ""),
            row_count=int(payload.get("row_count", 0)),
        )


@dataclass(frozen=True)
class ResultRecords:
    """The uniform record view of a result file.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        record_count: Records in the **whole file**, regardless of ``limit`` or ``offset``.
        returned_count: Records actually carried back in :attr:`columns`. Zero when ``out`` was
            given, since the table went to disk instead.
        offset: The offset that was applied.
        truncated: **Whether records were left behind**, by either ``limit`` or ``offset``. A short
            answer and a complete one must never look alike, so check this rather than assuming.
        retention_time_unit: The unit :attr:`columns`' ``retention_time`` carries for this format -
            ``"minutes"``, ``"seconds"``, or ``"unknown"``. mzLib does not normalise it, so this
            differs per format and you must convert before comparing two files. Provided as a value
            so you can convert programmatically instead of hard-coding a table.
        rows_not_read: Data rows in the file that did not become records - mzLib drops a malformed
            row silently, so a non-zero value here means the file is partly unreadable and the
            table is incomplete. ``None`` when the count could not be established meaningfully.
        caveats: **What the uniform view cannot be trusted to mean for this format.** Empty for
            some formats, not for others; each entry cites the mzLib source it came from. Worth
            printing before comparing anything across formats - this is where you learn that
            MSFragger retention times are seconds while MetaMorpheus's are minutes.
        column_names: The field names, in order.
        columns: Field name -> list of values, one entry per record - the shape ``pandas.DataFrame``
            and ``polars.DataFrame`` both accept directly. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None`` if it was returned inline.
    """

    path: str
    file_type: str
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    retention_time_unit: str
    rows_not_read: Any
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None

    @property
    def retention_time_in_minutes(self) -> list[Any]:
        """``retention_time`` converted to minutes, whatever unit the format wrote.

        The conversion you would otherwise write by hand, using
        :attr:`retention_time_unit`. Raises if the unit is ``"unknown"`` rather than guessing -
        a silently unconverted axis is the specific mistake this module exists to prevent.
        """
        values = (self.columns or {}).get("retention_time") or []
        if self.retention_time_unit == "minutes":
            return list(values)
        if self.retention_time_unit == "seconds":
            return [None if v is None else v / 60.0 for v in values]
        raise _bridge.UsageError(
            f"Cannot convert retention time for '{self.file_type}': mzLib gives no basis to say "
            "what unit it is in. Inspect the values against scan numbers before comparing them."
        )

    @property
    def records(self) -> list[dict[str, Any]]:
        """The same data row-wise: one dict per record.

        A convenience for looping. If you are building a table, prefer :attr:`columns` - it is
        already the shape a DataFrame wants, and this rebuilds it. Empty when ``out`` was given.
        """
        if not self.columns:
            return []
        names = self.column_names or list(self.columns)
        return [
            {name: self.columns[name][i] for name in names}
            for i in range(self.returned_count)
        ]

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "ResultRecords":
        written = data.get("output")
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            retention_time_unit=data.get("retention_time_unit") or "unknown",
            rows_not_read=data.get("rows_not_read"),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=WrittenTable._from_wire(written) if written else None,
        )


def formats(timeout: float | None = 60) -> list[Format]:
    """Every file type mzLib can recognise.

    Enumerated from mzLib itself rather than from a list maintained here, so it reflects the
    installed version and cannot go stale.

    Args:
        timeout: Seconds to allow.

    Returns:
        One :class:`Format` per supported file type.

    Example:
        >>> quantifiable = [f.file_type for f in formats() if f.is_quantifiable]  # doctest: +SKIP
        >>> quantifiable                                                          # doctest: +SKIP
        ['psmtsv', 'osmtsv', 'MsFraggerPsm']
    """
    data = _bridge.invoke("readers", "formats", timeout=timeout)
    return [Format._from_wire(item) for item in (data.get("formats") or [])]


def identify(path: str, timeout: float | None = 60) -> FileInfo:
    """Identify a result file without parsing its contents.

    Cheap by design: mzLib resolves the type and stops, so identifying a million-row file costs no
    more than identifying an empty one. It is not, however, *pure* - mzLib disambiguates a bare
    ``.tsv`` by reading its first line, a ``.mztab`` by its first five, and a Bruker ``.d`` by which
    analysis file the directory holds. An unreadable file will therefore raise.

    Args:
        path: Path to a result or spectra file. A Bruker ``.d`` directory is also accepted.
        timeout: Seconds to allow.

    Returns:
        A :class:`FileInfo` naming the format and the views it supports.

    Raises:
        UsageError: the path is blank, does not exist, or is not a file type mzLib recognises.
            mzLib has no "unknown" result - a file is dispatchable or it is an error - so use
            :func:`formats` to see what is supported, or catch this to test a file.

    Example:
        >>> info = identify("AllPSMs.psmtsv")                      # doctest: +SKIP
        >>> info.file_type, info.is_quantifiable                   # doctest: +SKIP
        ('psmtsv', True)
    """
    if not isinstance(path, str) or not path.strip():
        raise _bridge.UsageError("A file path is required, e.g. 'AllPSMs.psmtsv'.")

    data = _bridge.invoke("readers", "identify", "--path", path.strip(), timeout=timeout)
    return FileInfo._from_wire(data)


def read_results(
    path: str,
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | None = None,
    timeout: float | None = None,
) -> ResultRecords:
    """Read a result file into the uniform record view.

    Only the three file types offering the ``"quantifiable"`` view can be read this way - check
    :func:`identify` first, or catch the error. A file without the view is rejected with a message
    naming the views it does have.

    **There is no default row limit.** A result file can carry a million rows, and truncating by
    default would mean the ordinary call returns a table that looks complete and is not. For a large
    file use ``out`` rather than paging: see the note on ``offset`` below.

    Args:
        path: Path to a MetaMorpheus ``.psmtsv`` / ``.osmtsv`` or an MSFragger ``psm.tsv``.
        limit: Maximum records to return. ``None`` (the default) returns all of them.
            :attr:`ResultRecords.truncated` reports whether anything was left behind.
        offset: Records to skip. **This is a window, not a cursor.** mzLib materializes the whole
            file on every call - its readers look lazy and are not - so paging re-reads and
            re-parses the file once per page. For a large file, one call with ``out`` is right and
            a paging loop is quadratic.
        out: Write the records to this path as a **tab-separated** table and return only a summary,
            instead of carrying them back in the envelope. The intended path for large files, not
            an escape hatch. Tab-separated because these fields contain commas.
        timeout: Seconds to allow. A large file legitimately takes a while; ``None`` waits
            indefinitely.

    Returns:
        A :class:`ResultRecords`. Read :attr:`ResultRecords.caveats` before trusting a field across
        formats.

    Raises:
        UsageError: the path is blank, missing, not a recognised format, or has no quantifiable view.

    Example:
        >>> r = read_results("AllPSMs.psmtsv")                       # doctest: +SKIP
        >>> r.record_count, r.truncated                              # doctest: +SKIP
        (8, False)
        >>> import pandas as pd                                      # doctest: +SKIP
        >>> pd.DataFrame(r.columns)                                  # doctest: +SKIP
    """
    if not isinstance(path, str) or not path.strip():
        raise _bridge.UsageError("A file path is required, e.g. 'AllPSMs.psmtsv'.")

    args: list[str] = ["readers", "read-results", "--path", path.strip()]

    if limit is not None:
        if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
            raise _bridge.UsageError(f"limit must be a non-negative whole number or None; got {limit!r}.")
        args += ["--limit", str(limit)]

    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise _bridge.UsageError(f"offset must be a non-negative whole number; got {offset!r}.")
    if offset:
        args += ["--offset", str(offset)]

    if out is not None:
        if not isinstance(out, str) or not out.strip():
            raise _bridge.UsageError("out must be a non-empty path or None.")
        args += ["--out", out.strip()]

    data = _bridge.invoke(*args, timeout=timeout)
    return ResultRecords._from_wire(data)
