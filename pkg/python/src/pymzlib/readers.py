"""Read proteomics result files: what a file *is*, what you can do with it, and its records.

mzLib recognises 29 file types written by a dozen different search and deconvolution tools  - 
MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur,
FlashLFQ - and dispatches each to a parser it maintains. This module asks it what a path is::

    >>> import pymzlib
    >>> info = pymzlib.readers.identify("psm.tsv")     # doctest: +SKIP
    >>> info.file_type, info.views                     # doctest: +SKIP
    ('MsFraggerPsm', ['quantifiable'])

...and reads it, whatever it turns out to be::

    >>> table = pymzlib.readers.read_records("toppic_prsm.tsv")   # doctest: +SKIP
    >>> table.record_type, len(table.column_names)                # doctest: +SKIP
    ('ToppicPrsm', 36)

**Every one of the 29 formats is readable** - :func:`read_records` reads any of them. What differs
between formats is not *whether* you can read them but *what the columns mean*, and that is what
:attr:`FileInfo.views` tells you. It is tempting to describe mzLib as reading 29 formats into one
uniform shape; it does not. They fall into disjoint families, and several belong to no family at
all:

+---------------------+-------+------------------------+---------------------------------------+
| view                | types | function               | columns                               |
+=====================+=======+========================+=======================================+
| ``"quantifiable"``  | 3     | :func:`read_results`   | uniform: sequence, RT, charge, mass,  |
|                     |       |                        | protein groups. What                  |
|                     |       |                        | :func:`pymzlib.flashlfq.quantify`     |
|                     |       |                        | consumes.                             |
+---------------------+-------+------------------------+---------------------------------------+
| ``"ms1_features"``  | 2     | :func:`read_features`  | uniform: m/z, charge, RT range,       |
|                     |       |                        | intensity, isotope count.             |
+---------------------+-------+------------------------+---------------------------------------+
| ``"spectral_match"``| 4     | :func:`read_matches`   | uniform: scan, sequences, accession,  |
|                     |       |                        | decoy flag, modifications.            |
+---------------------+-------+------------------------+---------------------------------------+
| ``"spectra"``       | 7     | :func:`read_spectra`   | uniform: scan headers, and peaks on   |
|                     |       |                        | request.                              |
+---------------------+-------+------------------------+---------------------------------------+
| *(any)*             | 29    | :func:`read_records`   | **this format's own fields**, under   |
|                     |       |                        | mzLib's names. Not uniform.           |
+---------------------+-------+------------------------+---------------------------------------+

``views == []`` is a real and common answer - thirteen types have it. TopPIC, Crux, MSFragger's
peptide and protein tables and the FlashDeconv formats each parse into their own record type with
nothing in common. mzLib reads them and so does :func:`read_records`; there is simply no uniform
view to project them onto, and inventing one here would mean publishing a schema mzLib does not
have.

So: **use a typed view when you need numbers that mean the same thing across files, and**
:func:`read_records` **when you need everything a format has.** A ``.psmtsv`` read through
:func:`read_results` gives 10 comparable columns; the same file through :func:`read_records` gives
73, including the q-values and scores the uniform view does not carry.

Call :func:`formats` for the whole table. It is enumerated from mzLib rather than transcribed, so it
cannot drift from what mzLib actually dispatches.

**Three things this module deliberately does not tell you**, in the "surface it, don't hide it" spirit
of the rest of pyMzLib:

- **Which tool wrote the file.** mzLib has a ``Software`` property that looks like the answer and is
  not: readers carry their software constant on a constructor that mzLib's own file factory does not
  use, so the value is unset for everything the factory returns - and it is not reliably set on the
  other constructor either. Rather than reconstruct a plausible answer, there is no ``software``
  field. :attr:`FileInfo.file_type` already names the tool.
- **Whether the numbers inside mean the same thing across formats.** They do not, and this is the
  trap most likely to produce a wrong result. mzLib's result-file readers pass through whatever the
  tool wrote: MetaMorpheus retention times are in **minutes** and MSFragger's are too (mzLib PR
  #1116 converts them at the reader), but TopPIC's are still in **seconds**, and TopFD changed from
  seconds to minutes between v1.6.2 and v1.7.0 *within the same file type*. Likewise ``is_decoy`` is hardcoded ``False`` for MSFragger,
  which means "mzLib cannot tell" rather than "target" - MSFragger's ``psm.tsv`` carries no
  target/decoy column at all - so ``is_decoy`` arrives as **``None``** for that format rather than
  a fabricated ``False``. ``monoisotopic_mass`` is the *theoretical* peptide mass in **both**
  formats, never the observed precursor mass. Identifying a file is safe; comparing raw fields
  across formats is not.

- **Anything about confidence.** There is no q-value, PEP or score in this view, because
  ``IQuantifiableRecord`` carries only what FlashLFQ needs. **Nothing you get back is
  FDR-filtered**, even though every one of these files records confidence somewhere. Filter before
  you report.

.. note::

   That units mismatch was not hypothetical, and the fix shows where such things belong. Passing an
   MSFragger ``psm.tsv`` to :func:`pymzlib.flashlfq.quantify` used to return near-zero intensities,
   because FlashLFQ read the seconds as minutes and searched for each peptide about sixty times too
   late in the gradient. It was fixed **upstream in mzLib** (`#1116
   <https://github.com/smith-chem-wisc/mzLib/pull/1116>`_, converting at the reader) rather than
   papered over here, so every mzLib consumer benefits and this library's caveat and
   ``retention_time_unit`` changed with it. That is the standing rule: a value whose meaning or
   availability is wrong is repaired in the core contract, and a binding discloses rather than
   repairs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import _bridge

__all__ = [
    "QUANTIFIABLE",
    "MS1_FEATURES",
    "SPECTRAL_MATCH",
    "SPECTRA",
    "Format",
    "FileInfo",
    "WrittenTable",
    "ResultRecords",
    "NativeRecords",
    "FeatureRecords",
    "MatchRecords",
    "ScanRecords",
    "formats",
    "identify",
    "read_results",
    "read_records",
    "read_features",
    "read_matches",
    "read_spectra",
]

#: The view name for the cross-format record shape :func:`pymzlib.flashlfq.quantify` consumes.
QUANTIFIABLE = "quantifiable"

#: The view name for deconvolved MS1 features - :func:`read_features`.
MS1_FEATURES = "ms1_features"

#: The view name for records that are identifications - :func:`read_matches`.
SPECTRAL_MATCH = "spectral_match"

#: The view name for files that are spectra rather than results - :func:`read_spectra`.
SPECTRA = "spectra"


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


class _Table:
    """The row-wise view every columnar result shares.

    Not a dataclass and not public: it carries no fields of its own, only the two derivations that
    would otherwise be copy-pasted into five result types and drift. Each concrete class declares
    its own fields, because they genuinely differ - a spectra read has ``scan_count``, a feature
    read has ``retention_time_unit``, a native read has ``excluded_fields``.
    """

    @property
    def records(self) -> list[dict[str, Any]]:
        """The same data row-wise: one dict per record.

        A convenience for looping. If you are building a table, prefer ``columns`` - it is already
        the shape a DataFrame wants, and this rebuilds it. Empty when ``out`` was given.
        """
        columns = getattr(self, "columns", None)
        if not columns:
            return []
        names = [n for n in (getattr(self, "column_names", None) or list(columns)) if n in columns]
        # Row count from the columns themselves, not from returned_count: the two come from
        # different wire fields, and trusting the count would raise IndexError - or silently drop
        # rows - if they ever disagreed.
        length = min((len(columns[n]) for n in names), default=0)
        return [{name: columns[name][i] for name in names} for i in range(length)]

    def _require_columns(self, verb: str) -> dict[str, Any]:
        columns = getattr(self, "columns", None)
        if columns is None:
            where = getattr(getattr(self, "output", None), "path", None)
            raise _bridge.UsageError(
                "The records were written to "
                + (f"'{where}'" if where else "disk")
                + f" rather than returned, so there is nothing here to convert. Read that file, or "
                f"call {verb}() without out=."
            )
        return columns


@dataclass(frozen=True)
class ResultRecords(_Table):
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
            printing before comparing anything across formats - this is where you learn that, e.g.,
            TopPIC retention times are seconds while MetaMorpheus's and MSFragger's are minutes.
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
        self._require_columns("read_results")
        if "retention_time" not in self.columns:
            raise _bridge.UsageError(
                "This result has no retention_time column, so it cannot be converted."
            )

        values = self.columns["retention_time"]
        if self.retention_time_unit == "minutes":
            return list(values)
        if self.retention_time_unit == "seconds":
            return [None if v is None else v / 60.0 for v in values]
        raise _bridge.UsageError(
            f"Cannot convert retention time for '{self.file_type}': mzLib gives no basis to say "
            "what unit it is in. Inspect the values against scan numbers before comparing them."
        )

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


@dataclass(frozen=True)
class NativeRecords(_Table):
    """A result file read into **its own** fields, whatever format it is.

    What :func:`read_records` returns. Unlike every other result type in this module, the columns
    here are **not uniform**: they are the fields of this format's own mzLib record type, under
    mzLib's own names in snake_case. A TopPIC file gives you TopPIC's thirty-six columns; a Crux
    file gives you Crux's twenty-three. Always read :attr:`column_names` rather than assuming.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        reader: The mzLib class that parsed it.
        record_type: The mzLib record class the columns came from, e.g. ``"ToppicPrsm"``. Cross-
            reference it against the mzLib source to find out what a column means.
        views: The uniform views this file *also* supports, if any - see the module docstring.
        record_count: Records in the **whole file**, regardless of ``limit`` or ``offset``.
        returned_count: Records carried back in :attr:`columns`. Zero when ``out`` was given.
        offset: The offset that was applied.
        truncated: Whether records were left behind, by either ``limit`` or ``offset``.
        excluded_fields: **Fields of the record type that could not become columns**, each with the
            reason. Nested objects and dictionaries have no faithful column shape, and inventing
            one would mean publishing a schema mzLib does not have. Listed rather than dropped, so
            an absent column is never mistaken for an absent field.
        failed_fields: Fields that **raised** while being read, with the exception type. Several
            mzLib properties are computed and assume a UniProt-style FASTA header - Crux's and
            MsPathFinderT's ``accession`` are ``protein_id.split("|")[1]`` - so on other databases
            they throw. Those cells arrive ``None`` instead of taking the whole read down, but a
            failure must not look like missing data, so it is named here.
        column_names: The field names, in order - base-class fields first, then declared ones.
        columns: Field name -> list of values, one entry per record. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None`` if it was returned inline.
    """

    path: str
    file_type: str
    reader: Any
    record_type: str
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    views: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "NativeRecords":
        written = data.get("output")
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            reader=data.get("reader"),
            record_type=data.get("record_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            views=list(data.get("views") or []),
            excluded_fields=list(data.get("excluded_fields") or []),
            failed_fields=list(data.get("failed_fields") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=WrittenTable._from_wire(written) if written else None,
        )


@dataclass(frozen=True)
class FeatureRecords(_Table):
    """Deconvolved MS1 features, in the cross-format ``ms1_features`` view.

    What :func:`read_features` returns. Columns are ``mz``, ``charge``, ``retention_time_start``,
    ``retention_time_end``, ``intensity`` and ``number_of_isotopes`` - the same for every format
    that offers the view, so they *are* comparable across files, subject to
    :attr:`retention_time_unit`.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        record_count: Features in the whole file. **For ``_ms1.feature`` this exceeds the file's
            line count**, because mzLib expands each deconvolved feature into one row per charge
            state. See :attr:`caveats`.
        returned_count: Features carried back in :attr:`columns`. Zero when ``out`` was given.
        offset: The offset that was applied.
        truncated: Whether features were left behind, by either ``limit`` or ``offset``.
        retention_time_unit: ``"minutes"``, ``"seconds"``, or ``"unknown"``. **It is genuinely
            ``"unknown"`` for ``_ms1.feature``**: TopFD wrote seconds through v1.6.2 and minutes
            from v1.7.0 without changing the file type. That is not a gap in this library; it is
            the honest state of the format.
        caveats: What this view cannot be trusted to mean for this format, each citing the mzLib
            source it came from.
        column_names: The field names, in order.
        columns: Field name -> list of values. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
    """

    path: str
    file_type: str
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    retention_time_unit: str
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None

    @property
    def retention_time_start_in_minutes(self) -> list[Any]:
        """``retention_time_start`` in minutes, or a raised error if the unit is unknown.

        Raises rather than guessing. For ``_ms1.feature`` the unit genuinely is unknown, and a
        silently unconverted time axis is the specific mistake this module exists to prevent -
        mzLib's own deconvolution code guesses here, and this will not.
        """
        return self._converted("retention_time_start")

    @property
    def retention_time_end_in_minutes(self) -> list[Any]:
        """``retention_time_end`` in minutes, or a raised error if the unit is unknown."""
        return self._converted("retention_time_end")

    def _converted(self, name: str) -> list[Any]:
        columns = self._require_columns("read_features")
        if name not in columns:
            raise _bridge.UsageError(f"This result has no {name} column, so it cannot be converted.")
        values = columns[name]
        if self.retention_time_unit == "minutes":
            return list(values)
        if self.retention_time_unit == "seconds":
            return [None if v is None else v / 60.0 for v in values]
        raise _bridge.UsageError(
            f"Cannot convert retention time for '{self.file_type}': mzLib gives no basis to say what "
            "unit it is in. TopFD changed from seconds to minutes at v1.7.0 without changing the "
            "file type, so check the values against your gradient length before comparing them."
        )

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "FeatureRecords":
        written = data.get("output")
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            retention_time_unit=data.get("retention_time_unit") or "unknown",
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=WrittenTable._from_wire(written) if written else None,
        )


@dataclass(frozen=True)
class MatchRecords(_Table):
    """Identifications, in the cross-format ``spectral_match`` view.

    What :func:`read_matches` returns. Columns are ``file_name_without_extension``,
    ``one_based_scan_number``, ``base_sequence``, ``full_sequence``, ``accession``, ``is_decoy``,
    ``modifications`` and ``modification_count``.

    **Nothing here is FDR-filtered**, and unlike the quantifiable view there is not even a hint of
    confidence to filter on - mzLib's ``ISpectralMatch`` carries identity fields only. Every format
    that offers this view records an E-value or q-value in columns :func:`read_records` will give
    you. Filter before you report.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        record_count: Matches in the whole file.
        returned_count: Matches carried back in :attr:`columns`. Zero when ``out`` was given.
        offset: The offset that was applied.
        truncated: Whether matches were left behind, by either ``limit`` or ``offset``.
        caveats: What this view cannot be trusted to mean for this format - that MsPathFinderT
            infers decoys from an ``XXX`` name prefix, that Casanovo's scan numbers are mzTab
            indices, and that Casanovo's ``is_decoy`` is ``None`` because de novo sequencing has no
            target/decoy label at all.
        column_names: The field names, in order.
        columns: Field name -> list of values. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
    """

    path: str
    file_type: str
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "MatchRecords":
        written = data.get("output")
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=WrittenTable._from_wire(written) if written else None,
        )


@dataclass(frozen=True)
class ScanRecords(_Table):
    """Scan headers - and optionally peaks - from a spectra file.

    What :func:`read_spectra` returns. Retention times here **are** in minutes for every format:
    mzLib's spectra readers convert at the boundary, unlike its result-file readers, which pass the
    tool's own unit through untouched.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        reader: The mzLib class that parsed it, e.g. ``"Mzml"``, ``"ThermoRawFileReader"``.
        scan_count: Scans in the **whole file**, before any ``ms_order`` filter. Reported alongside
            :attr:`record_count` so a filter that matched nothing can never look like an empty file.
        ms_order: The MS level filtered to, or ``None`` if unfiltered.
        record_count: Scans that passed the ``ms_order`` filter.
        returned_count: Scans carried back in :attr:`columns`. Zero when ``out`` was given.
        offset: The offset that was applied.
        truncated: Whether scans were left behind, by either ``limit`` or ``offset``.
        peaks_included: Whether ``mz`` and ``intensity`` are present. When ``False``, ``peak_count``
            still tells you how many peaks each scan has.
        retention_time_unit: Always ``"minutes"`` for this view.
        caveats: What this view cannot be trusted to mean for this format - that msalign holds
            deconvolved neutral masses rather than m/z, that MGF scan numbers come from a title
            line, that Bruker needs Windows-x64 native libraries.
        column_names: The field names, in order.
        columns: Field name -> list of values. When :attr:`peaks_included`, ``mz`` and ``intensity``
            are each a **list of lists** - one array per scan. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
    """

    path: str
    file_type: str
    reader: Any
    scan_count: int
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    peaks_included: bool
    retention_time_unit: str
    ms_order: Any = None
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None

    @property
    def total_ion_current(self) -> list[Any]:
        """The ``total_ion_current`` column, for the commonest plot there is."""
        columns = self._require_columns("read_spectra")
        return list(columns.get("total_ion_current") or [])

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "ScanRecords":
        written = data.get("output")
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            reader=data.get("reader"),
            scan_count=int(data.get("scan_count", 0)),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            peaks_included=bool(data.get("peaks_included", False)),
            retention_time_unit=data.get("retention_time_unit") or "minutes",
            ms_order=data.get("ms_order"),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=WrittenTable._from_wire(written) if written else None,
        )


def _window(
    verb: str,
    path: str,
    *,
    limit: int | None,
    offset: int,
    out: str | None,
) -> list[str]:
    """The ``--path``/``--limit``/``--offset``/``--out`` argument list, validated.

    Shared by all five read verbs so they cannot drift on what counts as a usable argument. Every
    check raises before the bridge is spawned: a caller who passed ``limit=-1`` wants an error, not
    a process launch.
    """
    if not isinstance(path, str) or not path.strip():
        raise _bridge.UsageError("A file path is required, e.g. 'AllPSMs.psmtsv'.")

    args = ["readers", verb, "--path", path.strip()]

    if limit is not None:
        # bool first: `isinstance(True, int)` is True in Python, so limit=True would otherwise sail
        # through as limit=1 and silently return one row.
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

    return args


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
    args = _window("read-results", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return ResultRecords._from_wire(data)


def read_records(
    path: str,
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | None = None,
    timeout: float | None = None,
) -> NativeRecords:
    """Read **any** file mzLib recognises, into that format's own fields.

    This is the exhaustive verb: if :func:`identify` succeeds on a path, this reads it. All
    twenty-nine file types, including the thirteen that belong to no cross-format view at all -
    TopPIC, Crux, MSFragger's peptide and protein tables, the FlashDeconv formats - which no other
    function here can touch.

    **The columns are not uniform, by design.** They are this format's own mzLib record fields,
    under mzLib's own names in snake_case: a TopPIC file gives thirty-six columns, a Crux file
    twenty-three, an experiment annotation five. Read :attr:`NativeRecords.column_names`, and use
    :func:`read_results`, :func:`read_features` or :func:`read_matches` when you need columns that
    mean the same thing across formats.

    Nothing is silently dropped. A field that could not become a column is named in
    :attr:`NativeRecords.excluded_fields`, and one that raised while being read is named in
    :attr:`NativeRecords.failed_fields` - so a missing column never has to be guessed at.

    Args:
        path: Path to any file mzLib recognises. A Bruker ``.d`` directory is also accepted.
        limit: Maximum records to return. ``None`` (the default) returns all of them.
        offset: Records to skip. A window, not a cursor - see :func:`read_results`.
        out: Write a **tab-separated** table here and return only a summary.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`NativeRecords`.

    Raises:
        UsageError: the path is blank, missing, or not a file type mzLib recognises.

    Example:
        >>> r = read_records("toppic_prsm.tsv")                    # doctest: +SKIP
        >>> r.record_type, len(r.column_names)                     # doctest: +SKIP
        ('ToppicPrsm', 36)
        >>> import pandas as pd                                    # doctest: +SKIP
        >>> pd.DataFrame(r.columns)[["e_value", "q_value_spectrum_level"]]   # doctest: +SKIP
    """
    args = _window("read-records", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return NativeRecords._from_wire(data)


def read_features(
    path: str,
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | None = None,
    timeout: float | None = None,
) -> FeatureRecords:
    """Read deconvolved MS1 features, in the cross-format ``ms1_features`` view.

    Two file types offer it: TopFD/FLASHDeconv ``_ms1.feature`` and Dinosaur ``.feature.tsv``. A
    file without the view is rejected with a message naming the views it does have.

    **One row is not one line of the file for ``_ms1.feature``.** mzLib expands each deconvolved
    feature into one single-charge feature per charge in its recorded range, so a hundred-feature
    file can read as a thousand rows. Dinosaur is one-for-one. Both facts are in
    :attr:`FeatureRecords.caveats`, and :func:`read_records` gives the file's own rows either way.

    Args:
        path: Path to a ``_ms1.feature`` or Dinosaur ``.feature.tsv``.
        limit: Maximum features to return. ``None`` returns all of them.
        offset: Features to skip.
        out: Write a tab-separated table here and return only a summary.
        timeout: Seconds to allow.

    Returns:
        A :class:`FeatureRecords`. Check :attr:`FeatureRecords.retention_time_unit` before
        comparing times - it is ``"unknown"`` for ``_ms1.feature`` and that is the honest answer.

    Raises:
        UsageError: the path is blank, missing, unrecognised, or has no ``ms1_features`` view.

    Example:
        >>> f = read_features("sample_ms1.feature")                # doctest: +SKIP
        >>> f.record_count, f.retention_time_unit                  # doctest: +SKIP
        (25, 'unknown')
    """
    args = _window("read-features", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return FeatureRecords._from_wire(data)


def read_matches(
    path: str,
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | None = None,
    timeout: float | None = None,
) -> MatchRecords:
    """Read identifications, in the cross-format ``spectral_match`` view.

    Four file types offer it: MsPathFinderT's targets, decoys and combined results, and Casanovo's
    ``.mztab``. These are the identification formats that share no *file*-level interface, so
    :func:`read_results` cannot reach them.

    **Nothing here is FDR-filtered, and there is no confidence column to filter on** - mzLib's
    ``ISpectralMatch`` carries identity fields only. Every one of these formats records an E-value
    or q-value that :func:`read_records` will give you. Filter before you report.

    Args:
        path: Path to an MsPathFinderT ``_IcTarget.tsv`` / ``_IcDecoy.tsv`` / ``_IcTDA.tsv``, or a
            Casanovo ``.mztab``.
        limit: Maximum matches to return. ``None`` returns all of them.
        offset: Matches to skip.
        out: Write a tab-separated table here and return only a summary.
        timeout: Seconds to allow.

    Returns:
        A :class:`MatchRecords`. Read :attr:`MatchRecords.caveats` before trusting ``is_decoy`` -
        it is inferred from a name prefix for MsPathFinderT and is ``None`` for Casanovo.

    Example:
        >>> m = read_matches("results_IcTda.tsv")                  # doctest: +SKIP
        >>> m.record_count, m.columns["modifications"]             # doctest: +SKIP
        (6, ['', '12:Oxidation on M', '', '', '4:Acetylation on K', ''])
    """
    args = _window("read-matches", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return MatchRecords._from_wire(data)


def read_spectra(
    path: str,
    *,
    limit: int | None = None,
    offset: int = 0,
    ms_order: int | None = None,
    peaks: bool = False,
    out: str | None = None,
    timeout: float | None = None,
) -> ScanRecords:
    """Read the scans of a spectra file: headers always, peaks on request.

    Seven file types offer the ``spectra`` view: ``.mzML``, ``.mgf``, ``_ms1.msalign``,
    ``_ms2.msalign``, Thermo ``.raw``, Bruker ``.d`` and timsTOF ``.d``.

    **Peaks are opt-in and should stay that way unless you need them.** A scan header is tens of
    bytes; its peak list is thousands, and a mid-size mzML holds tens of thousands of scans. With
    ``peaks=True`` the ``mz`` and ``intensity`` columns each become a list of arrays, one per scan -
    so pair it with ``limit``, ``ms_order`` or ``out``.

    Args:
        path: Path to a spectra file. A Bruker ``.d`` directory is also accepted.
        limit: Maximum scans to return. ``None`` returns all of them.
        offset: Scans to skip, applied **after** ``ms_order``.
        ms_order: Keep only scans at this MS level - ``1`` for survey scans, ``2`` for fragment
            scans. Applied before ``offset`` and ``limit``, so ``ms_order=2, limit=10`` means the
            first ten MS2 scans rather than the MS2 scans among the first ten.
        peaks: Include the ``mz`` and ``intensity`` arrays. Off by default.
        out: Write a tab-separated table here and return only a summary. With ``peaks=True`` each
            cell holds a ``;``-joined list.
        timeout: Seconds to allow. Reading a large ``.raw`` legitimately takes a while.

    Returns:
        A :class:`ScanRecords`. Retention times are in minutes for every format.

    Raises:
        UsageError: the path is blank, missing, unrecognised, has no ``spectra`` view, or
            ``ms_order`` is less than 1.

    Example:
        >>> s = read_spectra("run.mzML", ms_order=2, limit=5)      # doctest: +SKIP
        >>> s.scan_count, s.record_count                           # doctest: +SKIP
        (14238, 11902)
        >>> s.columns["selected_ion_mz"]                           # doctest: +SKIP
        [447.7391, 551.2903, 638.8215, 712.3344, 805.9012]
    """
    args = _window("read-spectra", path, limit=limit, offset=offset, out=out)

    if ms_order is not None:
        if isinstance(ms_order, bool) or not isinstance(ms_order, int) or ms_order < 1:
            raise _bridge.UsageError(
                f"ms_order must be a whole number of 1 or more, or None; got {ms_order!r}."
            )
        args += ["--ms-order", str(ms_order)]

    if peaks:
        args += ["--peaks"]

    data = _bridge.invoke(*args, timeout=timeout)
    return ScanRecords._from_wire(data)
