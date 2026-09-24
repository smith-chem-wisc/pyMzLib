"""Read mass-spectrometry data files and proteomics search results: what a file *is*, what you
can do with it, and its records.

**Spectra files are read here too, not just search output.** :func:`read_spectra` reads
**mzML**, Thermo ``.raw``, Bruker ``.d``, timsTOF ``.d``, MGF and msalign - scan headers always,
peaks on request::

    >>> import pymzlib
    >>> scans = pymzlib.readers.read_spectra("sliced_ethcd.mzML", limit=3)
    >>> scans.scan_count, scans.columns["retention_time"][:2]
    (6, [38.92571663975, 38.92606115255])

mzLib recognises 36 file types in all - the instrument and deconvolution formats above, plus the
output of a dozen search tools: MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux,
Casanovo, FlashDeconv, Dinosaur, DIA-NN, FlashLFQ, Pytheas, and any mzIdentML writer - and
dispatches each to a parser it maintains. This module asks it what a path is::

    >>> # Not run in CI: no "readers identify" recording to replay yet.
    >>> info = pymzlib.readers.identify("psm.tsv")     # doctest: +SKIP
    >>> info.file_type, info.views                     # doctest: +SKIP
    ('MsFraggerPsm', ['quantifiable'])

...and reads it, whatever it turns out to be::

    >>> table = pymzlib.readers.read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv")
    >>> table.record_type, len(table.column_names)
    ('ToppicPrsm', 36)

**Every one of the 36 formats is readable** - :func:`read_records` reads any of them. What differs
between formats is not *whether* you can read them but *what the columns mean*, and that is what
:attr:`FileInfo.views` tells you. It is tempting to describe mzLib as reading 36 formats into one
uniform shape; it does not. They fall into disjoint families, and several belong to no family at
all:

+---------------------+-------+------------------------+---------------------------------------+
| view                | types | function               | columns                               |
+=====================+=======+========================+=======================================+
| ``"quantifiable"``  | 4     | :func:`read_results`   | uniform: sequence, RT, charge, mass,  |
|                     |       |                        | protein groups. What                  |
|                     |       |                        | :func:`pymzlib.flashlfq.quantify`     |
|                     |       |                        | consumes.                             |
+---------------------+-------+------------------------+---------------------------------------+
| ``"ms1_features"``  | 2     | :func:`read_features`  | uniform: m/z, charge, RT range,       |
|                     |       |                        | intensity, isotope count.             |
+---------------------+-------+------------------------+---------------------------------------+
| ``"spectral_match"``| 6     | :func:`read_matches`   | uniform: scan, sequences, accession,  |
|                     |       |                        | decoy flag, modifications.            |
+---------------------+-------+------------------------+---------------------------------------+
| ``"spectra"``       | 7     | :func:`read_spectra`   | uniform: scan headers, and peaks on   |
|                     |       |                        | request.                              |
+---------------------+-------+------------------------+---------------------------------------+
| *(any)*             | 36    | :func:`read_records`   | **this format's own fields**, under   |
|                     |       |                        | mzLib's names. Not uniform.           |
+---------------------+-------+------------------------+---------------------------------------+

``views == []`` is a real and common answer - seventeen types have it. TopPIC, Crux, MSFragger's
peptide and protein tables, the FlashDeconv formats and MetaMorpheus's protein-group and peptide
tables each parse into their own record type with nothing in common. mzLib reads them and so does :func:`read_records`; there is simply no uniform
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

**Two quantification tables have typed functions of their own** (mzLib 1.0.592, #1347). They
belong to no uniform view - :func:`read_records` reads them too - but their per-sample values are
dictionaries that :func:`read_records` cannot project, so these return them as **long** tables,
one row per record per sample:

+-------------------------------------+------------------------------------------------------------+
| function                            | reads                                                      |
+=====================================+============================================================+
| :func:`read_protein_groups`         | MetaMorpheus ``AllQuantifiedProteinGroups.tsv``: one row   |
|                                     | per group per sample group - spectral count, intensity.    |
+-------------------------------------+------------------------------------------------------------+
| :func:`read_quantified_peptides`    | FlashLFQ ``QuantifiedPeptides.tsv`` (and MetaMorpheus's    |
|                                     | ``AllQuantifiedPeptides.tsv``): one row per peptide per    |
|                                     | sample - intensity, detection type.                        |
+-------------------------------------+------------------------------------------------------------+
| :func:`read_occupancy`              | the same protein-group table's PTM site occupancy: one row |
|                                     | per group, sample group, basis and modified site.          |
+-------------------------------------+------------------------------------------------------------+

**Many files at once.** Every read function has a ``_many`` twin - :func:`read_spectra_many`,
:func:`read_records_many`, :func:`identify_many` and so on - that takes a list of paths and
returns ONE long table, whose first two columns say which file each row came from
(``source_index``, ``source_path``), plus one :class:`FileReport` per file in
:attr:`ReadBatch.files`. The list is read by one bridge process, ``threads`` files at a time, so
two hundred runs cost one .NET start-up rather than two hundred. The answer is identical at any
``threads``; the default is 1 because each file in flight is held whole in memory. There is
deliberately no Python-side thread or process pool anywhere in this module: parallelism is the
bridge's to express, once, where it can be counted.

**Four ways a field can have no value, and each is named.** A ``None`` in a table can mean four
different things, and conflating them is how a missing column turns into a measured zero:

- ``absent_fields`` - the function defines the field, but **this file's format has no column for
  it**: MBR Score in a current FlashLFQ peaks table (mzLib #1345), ``q_value`` in an MsPathFinderT
  targets file, apex intensity in a FLASHDeconv feature file. Every value in the column is
  ``None``, whatever default mzLib filled in.
- ``failed_fields`` - the field exists, but **reading it threw** for some rows (a UniProt-shaped
  accession parsed from a plain FASTA header). ``None`` in those rows.
- ``excluded_fields`` - the field has **no column shape** (a dictionary, a nested object), so it
  does not cross in this function at all; each entry names the function that does carry it.
- a ``None`` in one row, with the field in none of those lists - the value is **genuinely missing
  for that row**: the precursor of an MS1 scan, a blank intensity cell.

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

import os
from dataclasses import dataclass, field
from typing import Any, Iterable, Union

from . import _bridge

__all__ = [
    "QUANTIFIABLE",
    "MS1_FEATURES",
    "SPECTRAL_MATCH",
    "SPECTRA",
    "Format",
    "FileInfo",
    "WrittenTable",
    "ReadError",
    "SpectraSource",
    "SkippedMatch",
    "FileReport",
    "ResultRecords",
    "NativeRecords",
    "FeatureRecords",
    "MatchRecords",
    "ScanRecords",
    "ProteinGroupRecords",
    "QuantifiedPeptideRecords",
    "OccupancyRecords",
    "ReadBatch",
    "IdentifyBatch",
    "formats",
    "identify",
    "identify_many",
    "read_results",
    "read_results_many",
    "read_records",
    "read_records_many",
    "read_features",
    "read_features_many",
    "read_matches",
    "read_matches_many",
    "read_spectra",
    "read_spectra_many",
    "read_protein_groups",
    "read_protein_groups_many",
    "read_quantified_peptides",
    "read_quantified_peptides_many",
    "read_occupancy",
    "read_occupancy_many",
]

#: The view name for the cross-format record shape :func:`pymzlib.flashlfq.quantify` consumes.
QUANTIFIABLE = "quantifiable"

#: The view name for deconvolved MS1 features - :func:`read_features`.
MS1_FEATURES = "ms1_features"

#: The view name for records that are identifications - :func:`read_matches`.
SPECTRAL_MATCH = "spectral_match"

#: The view name for files that are spectra rather than results - :func:`read_spectra`.
SPECTRA = "spectra"

#: A list of paths for a ``_many`` function: each a ``str`` or ``pathlib.Path``.
PathList = Iterable[Union[str, "os.PathLike[str]"]]

#: The pyMzLib release whose bridge first has the mzLib 1.0.592 quantification views - the
#: ``since.pymzlib`` of their specs. See :func:`pymzlib._bridge.require_verb`.
_QUANT_VIEWS_SINCE = "0.2.0"


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
class ReadError:
    """Why one file of a ``_many`` read could not be read, under ``on_error="skip"``.

    Attributes:
        kind: ``"usage"`` - the input was wrong: not found, not a type mzLib recognises, or not one
            this function reads; ``"correctness"`` - mzLib recognised the file and failed to parse
            it; or ``"service_unavailable"``. The same three-way split every pyMzLib error makes
            (:class:`~pymzlib.UsageError`, :class:`~pymzlib.BridgeError`,
            :class:`~pymzlib.ServiceUnavailableError`).
        type: What the failure would have been raised as on its own: ``"usage"``, or the .NET
            exception type, e.g. ``"MzLibException"``, ``"HeaderValidationException"``.
        message: mzLib's own message.
    """

    kind: str
    type: str
    message: str

    @classmethod
    def _from_wire(cls, payload: dict[str, Any] | None) -> "ReadError | None":
        if not payload:
            return None
        return cls(
            kind=payload.get("kind", ""),
            type=payload.get("type", ""),
            message=payload.get("message", ""),
        )


@dataclass(frozen=True)
class SpectraSource:
    """What a spectra file records about the run it came from (mzLib ``SourceFile``, #1349).

    The facts a batch or instrument confound is built from: which instrument, which unit of that
    model, and when. Every field is ``None`` when the file does not record it - MGF and msalign
    record none of them.

    Attributes:
        instrument_model: The instrument model's name, e.g. ``"Orbitrap Fusion Lumos"``. ``None``:
            the file does not record it.
        instrument_model_accession: Its PSI-MS accession, e.g. ``"MS:1002416"``. mzML carries one;
            Thermo ``.raw`` records only the name, so it is ``None`` there rather than looked up.
            mzLib's rule is to **match on the accession, never the name**.
        instrument_serial_number: The instrument's serial number, verbatim, e.g. ``"FSN10189"``.
            Together with the model it tells two instruments of one model apart. Free text: a
            converter's placeholder such as ``"Serial Number N/A"`` is reported as written.
            ``None``: the file does not record it.
        acquisition_start_time: When acquisition started, as ISO-8601 text, e.g.
            ``"2021-03-16T17:09:07Z"``. It ends in ``Z`` only when the source fixed the instant;
            otherwise it is the acquisition computer's wall-clock time with no offset, because the
            instant is unknown. ``None``: the file does not record it.
        acquisition_start_time_is_utc: ``True`` when :attr:`acquisition_start_time` is a UTC
            instant; ``False`` when it is local instrument-PC time (always, for Thermo ``.raw``)
            or absent. A ``.raw`` and ProteoWizard's mzML of it can differ by the site's UTC
            offset, because ProteoWizard assumes the converting machine's time zone.
    """

    instrument_model: str | None = None
    instrument_model_accession: str | None = None
    instrument_serial_number: str | None = None
    acquisition_start_time: str | None = None
    acquisition_start_time_is_utc: bool = False

    @property
    def acquired_at(self) -> Any:
        """:attr:`acquisition_start_time` as a :class:`datetime.datetime`, or ``None``.

        Timezone-aware (UTC) when :attr:`acquisition_start_time_is_utc`, and **naive** otherwise:
        a local clock reading with no known offset must not be given one.
        """
        if not self.acquisition_start_time:
            return None
        import datetime

        text = self.acquisition_start_time
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        return datetime.datetime.fromisoformat(text)

    @classmethod
    def _from_wire(cls, payload: dict[str, Any] | None) -> "SpectraSource | None":
        if not payload:
            return None
        return cls(
            instrument_model=payload.get("instrument_model"),
            instrument_model_accession=payload.get("instrument_model_accession"),
            instrument_serial_number=payload.get("instrument_serial_number"),
            acquisition_start_time=payload.get("acquisition_start_time"),
            acquisition_start_time_is_utc=bool(payload.get("acquisition_start_time_is_utc", False)),
        )


@dataclass(frozen=True)
class SkippedMatch:
    """An mzIdentML identification item mzLib did not turn into a row, and why (mzLib #1313).

    Attributes:
        spectrum_identification_item_id: The item's ``id`` in the document.
        spectrum_id: The spectrum's nativeID, as written.
        reason: Why it was skipped, e.g. ``"crosslink identification"``, or a modification that
            does not resolve to a Unimod entry.
    """

    spectrum_identification_item_id: str
    spectrum_id: str
    reason: str

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "SkippedMatch":
        return cls(
            spectrum_identification_item_id=payload.get("spectrum_identification_item_id", ""),
            spectrum_id=payload.get("spectrum_id", ""),
            reason=payload.get("reason", ""),
        )


def _skipped(payload: Any) -> list[SkippedMatch] | None:
    return None if payload is None else [SkippedMatch._from_wire(item) for item in payload]


@dataclass(frozen=True)
class FileReport:
    """The per-file facts of one input to a ``_many`` read: one entry of :attr:`ReadBatch.files`.

    One per input, **in input order**; a file's rows in the batch table are those whose
    ``source_index`` is its position. The fields are the ones a single-file read reports at its top
    level, under the same names, so what you learn about one reads the same for many.

    The fields every reader reports come first; the rest belong to particular readers and are
    ``None`` elsewhere. For an input that failed under ``on_error="skip"``, :attr:`error` says why
    and the facts that could not be established are ``None`` or empty.

    Attributes:
        path: The absolute path of the input.
        file_type: mzLib's ``SupportedFileType`` name. ``None``: the type could not be determined
            (a failed input).
        reader: The mzLib class that parsed the file. ``None``: a failed input.
        record_count: Records mzLib read from the whole file, in the reader's own unit - matches,
            scans passing ``ms_order``, groups, peptides. ``None``: a failed input - not zero,
            because nothing was counted.
        rows_not_read: Data lines that did not become records, counted where one line is one
            record (the psmtsv family and MSFragger); mzLib drops a malformed line silently.
            ``None``: not established for this format, or a failed input.
        retention_time_unit: The unit of the reader's time columns for this file: ``"minutes"``,
            ``"seconds"`` or ``"unknown"``. ``None``: the reader has no time column, or it is
            :func:`read_records_many`, whose columns are the format's own.
        caveats: What the reader cannot be trusted to mean for this file, each citing the mzLib
            source.
        column_names: This file's columns, without the ``source_index``/``source_path`` pair the
            batch table adds.
        absent_fields: Columns the reader defines but **this file has no source for** - its
            format has no such column - so every value in them is ``None``. See the module
            docstring's "four ways a field can have no value".
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows;
            those cells are ``None``.
        excluded_fields: Fields that have no column shape and so do not cross, each as
            ``{field, type, reason, verb}``, where ``verb`` names the bridge command that does carry
            it (``None`` if none does).
        error: Why this input could not be read, under ``on_error="skip"``; ``None`` when it was.
        record_type: :func:`read_records_many` only - the mzLib record class the columns belong to.
        views: :func:`read_records_many` only - the uniform views this file also offers.
        scan_count: :func:`read_spectra_many` only - scans in the whole file, before ``ms_order``.
        source: :func:`read_spectra_many` only - the instrument and acquisition facts.
        skipped_count: mzIdentML items mzLib did not represent (:func:`read_matches_many`,
            :func:`read_records_many`). ``None``: the format has no such list.
        skipped: The same items, each a :class:`SkippedMatch` with its reason.
        sample_labels: The quantification readers only - the file's sample labels, in header
            order.
        truncated_cell_count: :func:`read_occupancy_many` only - occupancy cells the writer cut
            short.
    """

    path: str
    file_type: str | None = None
    reader: str | None = None
    record_count: int | None = None
    rows_not_read: int | None = None
    retention_time_unit: str | None = None
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None
    record_type: str | None = None
    views: list[str] | None = None
    scan_count: int | None = None
    source: SpectraSource | None = None
    skipped_count: int | None = None
    skipped: list[SkippedMatch] | None = None
    sample_labels: list[str] | None = None
    truncated_cell_count: int | None = None

    @property
    def ok(self) -> bool:
        """Whether this input was read (``error`` is ``None``)."""
        return self.error is None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "FileReport":
        count = data.get("record_count")
        views = data.get("views")
        labels = data.get("sample_labels")
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type"),
            reader=data.get("reader"),
            record_count=None if count is None else int(count),
            rows_not_read=data.get("rows_not_read"),
            retention_time_unit=data.get("retention_time_unit"),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            absent_fields=list(data.get("absent_fields") or []),
            failed_fields=list(data.get("failed_fields") or []),
            excluded_fields=list(data.get("excluded_fields") or []),
            error=ReadError._from_wire(data.get("error")),
            record_type=data.get("record_type"),
            views=None if views is None else list(views),
            scan_count=data.get("scan_count"),
            source=SpectraSource._from_wire(data.get("source")),
            skipped_count=data.get("skipped_count"),
            skipped=_skipped(data.get("skipped")),
            sample_labels=None if labels is None else list(labels),
            truncated_cell_count=data.get("truncated_cell_count"),
        )


@dataclass(frozen=True)
class FileInfo:
    """What a particular file is, and what can be done with it.

    Attributes:
        path: The absolute path that was identified.
        file_type: mzLib's ``SupportedFileType`` name. Empty for an input of
            :func:`identify_many` that failed under ``on_error="skip"``.
        extension: The extension mzLib dispatched on. ``None`` for a failed input.
        reader: The mzLib class that would parse it. ``None`` for a failed input.
        views: The uniform views this file supports - see the module docstring. Often empty, which
            means mzLib can read the file but offers no cross-format projection of it.
        error: Always ``None`` from :func:`identify`, which raises instead. From
            :func:`identify_many` with ``on_error="skip"``, why this file could not be identified.
    """

    path: str
    file_type: str
    extension: Any
    reader: Any
    views: list[str] = field(default_factory=list)
    error: ReadError | None = None

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
            file_type=data.get("file_type") or "",
            extension=data.get("extension"),
            reader=data.get("reader"),
            views=list(data.get("views") or []),
            error=ReadError._from_wire(data.get("error")),
        )


@dataclass(frozen=True)
class WrittenTable:
    """Where a read wrote its table, when asked to write one instead of returning it.

    Attributes:
        path: The absolute path written.
        format: Always ``"tsv"``. **Tab-separated, not comma-separated**, because these fields
            contain commas - MSFragger's mapped proteins are a comma-separated list inside a single
            field, and joined accessions look the same. It is also what every mzLib reader and
            writer uses. Read it with ``csv.reader(f, delimiter="\\t")`` or
            ``pandas.read_csv(path, sep="\\t")``.
        row_count: Rows written, excluding the header - for a ``_many`` read, every file's rows.
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


def _output(data: dict[str, Any]) -> WrittenTable | None:
    written = data.get("output")
    return WrittenTable._from_wire(written) if written else None


def _block(data: dict[str, Any]) -> dict[str, Any]:
    """The per-file fields every single-file result carries, from its payload."""
    return dict(
        reader=data.get("reader"),
        absent_fields=list(data.get("absent_fields") or []),
        failed_fields=list(data.get("failed_fields") or []),
        excluded_fields=list(data.get("excluded_fields") or []),
        error=ReadError._from_wire(data.get("error")),
    )


class _Table:
    """The row-wise view, and the unit conversion, every columnar result shares.

    Not a dataclass and not public: it carries no fields of its own, only the derivations that
    would otherwise be copy-pasted into nine result types and drift. Each concrete class declares
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

    def _to_minutes(self, verb: str, name: str) -> list[Any]:
        """A time column in minutes, converting each file by its own ``retention_time_unit``.

        Raises rather than guessing when a unit is ``"unknown"``: a silently unconverted axis is
        the specific mistake this module exists to prevent.
        """
        columns = self._require_columns(verb)
        if name not in columns:
            raise _bridge.UsageError(f"This result has no {name} column, so it cannot be converted.")

        values = columns[name]
        files = getattr(self, "files", None)
        if files is not None:
            sources = [files[i] for i in columns["source_index"]]
        else:
            sources = [self] * len(values)

        converted = []
        for value, source in zip(values, sources):
            unit = getattr(source, "retention_time_unit", None)
            if unit == "minutes":
                converted.append(value)
            elif unit == "seconds":
                converted.append(None if value is None else value / 60.0)
            else:
                raise _bridge.UsageError(
                    f"Cannot convert {name} for '{getattr(source, 'file_type', '')}' "
                    f"({getattr(source, 'path', '')}): mzLib gives no basis to say what unit it is "
                    "in. TopFD, for one, changed from seconds to minutes at v1.7.0 without changing "
                    "the file type, so check the values against your gradient length before "
                    "comparing them."
                )
        return converted


@dataclass(frozen=True)
class ResultRecords(_Table):
    """The uniform record view of a result file.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        record_count: Records in the **whole file**, regardless of ``limit`` or ``offset``.
        returned_count: Records actually carried back in :attr:`columns`. Zero when ``out`` was
            given, since the table went to disk instead.
        offset: Records skipped: the offset that was applied.
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
        reader: The mzLib class that parsed the file, e.g. ``"PsmFromTsvFile"``.
        absent_fields: Columns this view defines but **this file's format has no column for**, so
            every value is ``None``: ``["is_decoy"]`` for MSFragger and DIA-NN, which write no
            target/decoy label. See the module docstring's "four ways a field can have no value".
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows;
            those cells are ``None``.
        excluded_fields: Fields with no column shape. Always empty for this view.
        error: Always ``None``: a file that cannot be read raises instead. Present because it is
            part of the per-file facts every reader reports (see :class:`FileReport`).
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
    reader: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @property
    def retention_time_in_minutes(self) -> list[Any]:
        """``retention_time`` converted to minutes, whatever unit the format wrote.

        The conversion you would otherwise write by hand, using
        :attr:`retention_time_unit`. Raises if the unit is ``"unknown"`` rather than guessing -
        a silently unconverted axis is the specific mistake this module exists to prevent.
        """
        return self._to_minutes("read_results", "retention_time")

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "ResultRecords":
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
            output=_output(data),
            **_block(data),
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
        record_count: Records in the **whole file**, regardless of ``limit`` or ``offset``.
        returned_count: Records carried back in :attr:`columns`. Zero when ``out`` was given.
        offset: Records skipped: the offset that was applied.
        truncated: Whether records were left behind, by either ``limit`` or ``offset``.
        views: The uniform views this file *also* supports, if any - see the module docstring.
        excluded_fields: **Fields of the record type that could not become columns**, each as
            ``{field, type, reason, verb}``. Nested objects and dictionaries have no faithful
            column shape, and inventing one would mean publishing a schema mzLib does not have.
            Listed rather than dropped, so an absent column is never mistaken for an absent field;
            ``verb`` names the bridge command that does carry it - ``sample_groups`` points at
            ``readers read-protein-groups`` (:func:`read_protein_groups`), mzIdentML's ``scores``
            at ``readers read-matches`` (:func:`read_matches` with ``scores=True``).
        failed_fields: Fields that **raised** while being read, with the exception type. Several
            mzLib properties are computed and assume a UniProt-style FASTA header - Crux's and
            MsPathFinderT's ``accession`` are ``protein_id.split("|")[1]`` - so on other databases
            they throw. Those cells arrive ``None`` instead of taking the whole read down, but a
            failure must not look like missing data, so it is named here.
        absent_fields: Columns **this file has no source for**: the record type reads them from an
            optional column the file does not have, so every value is mzLib's default, and is
            ``None`` here. ``["mbr_score"]`` on a current FlashLFQ peaks table (mzLib #1345);
            ``["q_value", "pep_q_value"]`` on an MsPathFinderT targets file, where mzLib would
            otherwise report a q-value of 0 for every row.
        column_names: The field names, in order - base-class fields first, then declared ones.
        columns: Field name -> list of values, one entry per record. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None`` if it was returned inline.
        rows_not_read: Data rows that did not become records, counted where one line is one
            record (the psmtsv family and MSFragger), since mzLib drops a malformed line silently.
            ``None``: not established for this format.
        retention_time_unit: Always ``None``: the columns are this format's own, and several
            formats carry more than one time column in different units. Use a typed view when
            units matter.
        caveats: Always empty for this function.
        skipped_count: mzIdentML only - identification items mzLib did not represent as records
            (crosslinks, unresolvable modifications, substitutions), so ``record_count +
            skipped_count`` is the items in the file. ``None`` for every other format, which keeps
            no such list.
        skipped: The same items, each a :class:`SkippedMatch` with its reason. ``None`` with
            :attr:`skipped_count`.
        error: Always ``None``: a file that cannot be read raises instead.
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
    absent_fields: list[str] = field(default_factory=list)
    rows_not_read: Any = None
    retention_time_unit: Any = None
    caveats: list[str] = field(default_factory=list)
    skipped_count: int | None = None
    skipped: list[SkippedMatch] | None = None
    error: ReadError | None = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "NativeRecords":
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
            output=_output(data),
            absent_fields=list(data.get("absent_fields") or []),
            rows_not_read=data.get("rows_not_read"),
            retention_time_unit=data.get("retention_time_unit"),
            caveats=list(data.get("caveats") or []),
            skipped_count=data.get("skipped_count"),
            skipped=_skipped(data.get("skipped")),
            error=ReadError._from_wire(data.get("error")),
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
        offset: Features skipped: the offset that was applied.
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
        reader: The mzLib class that parsed the file, e.g. ``"Ms1FeatureFile"``.
        rows_not_read: Always ``None``: rows that did not become records are not counted for this
            view, because one line of an ``_ms1.feature`` is not one feature.
        absent_fields: Columns this file has no source for, so every value is ``None``:
            ``number_of_isotopes`` for every ``_ms1.feature`` (mzLib never sets it there), and
            ``intensity`` for a FLASHDeconv file, whose schema has no apex intensity - mzLib would
            otherwise report zero for every feature.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows.
        excluded_fields: Always empty for this view.
        error: Always ``None``: a file that cannot be read raises instead.
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
    reader: Any = None
    rows_not_read: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @property
    def retention_time_start_in_minutes(self) -> list[Any]:
        """``retention_time_start`` in minutes, or a raised error if the unit is unknown.

        Raises rather than guessing. For ``_ms1.feature`` the unit genuinely is unknown, and a
        silently unconverted time axis is the specific mistake this module exists to prevent -
        mzLib's own deconvolution code guesses here, and this will not.
        """
        return self._to_minutes("read_features", "retention_time_start")

    @property
    def retention_time_end_in_minutes(self) -> list[Any]:
        """``retention_time_end`` in minutes, or a raised error if the unit is unknown."""
        return self._to_minutes("read_features", "retention_time_end")

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "FeatureRecords":
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
            output=_output(data),
            rows_not_read=data.get("rows_not_read"),
            **_block(data),
        )


@dataclass(frozen=True)
class MatchRecords(_Table):
    """Identifications, in the cross-format ``spectral_match`` view.

    What :func:`read_matches` returns. Columns are ``file_name_without_extension``,
    ``one_based_scan_number``, ``base_sequence``, ``full_sequence``, ``accession``, ``is_decoy``,
    ``modifications``, ``modification_count``, ``q_value``, ``rank`` and ``pass_threshold`` - and,
    with ``scores=True``, ``match_index``, ``score_name`` and ``score_value``.

    **Nothing here is FDR-filtered.** ``q_value`` is filled only by mzIdentML and by an
    MsPathFinderT file that carries a ``QValue`` column; :attr:`absent_fields` names it for every
    other file. Filter before you report.

    Attributes:
        path: The absolute path that was read.
        file_type: mzLib's ``SupportedFileType`` name.
        record_count: Matches in the whole file.
        returned_count: Matches carried back in :attr:`columns` - the unit ``limit`` and
            ``offset`` count in. Zero when ``out`` was given.
        row_count: Rows in :attr:`columns`: one per match, or with ``scores=True`` one per match
            and score, so it can exceed :attr:`returned_count`. Zero when ``out`` was given.
        offset: Matches skipped: the offset that was applied.
        truncated: Whether matches were left behind, by either ``limit`` or ``offset``.
        caveats: What this view cannot be trusted to mean for this format - that MsPathFinderT
            infers decoys from an ``XXX`` name prefix, that Casanovo's scan numbers are mzTab
            indices, that mzIdentML lists every candidate rank.
        column_names: The field names, in order.
        columns: Field name -> list of values. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
        reader: The mzLib class that parsed the file, e.g. ``"MzIdentMLResultFile"``.
        scores_included: Whether ``scores=True`` was asked for, making the table long.
        skipped_count: mzIdentML only - identification items mzLib did not turn into rows
            (crosslinks, unresolvable modifications, substitutions); ``record_count +
            skipped_count`` is the items in the file. ``None`` for every other format, which keeps
            no such list.
        skipped: The same items, each a :class:`SkippedMatch` with its reason. ``None`` with
            :attr:`skipped_count`.
        rows_not_read: Always ``None``: rows that did not become records are not counted for this
            view.
        retention_time_unit: Always ``None``: this view has no time column.
        absent_fields: Columns this file's format has no source for, so every value is ``None``:
            ``is_decoy`` for every format but MsPathFinderT; ``q_value`` for Casanovo and an
            MsPathFinderT file without a ``QValue`` column (mzLib would otherwise report 0, a
            perfect q-value); ``rank`` and ``pass_threshold`` for everything but mzIdentML; with
            ``scores=True``, ``score_name`` and ``score_value`` for everything but mzIdentML.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows -
            MsPathFinderT's ``accession`` on a FASTA without UniProt-style headers.
        excluded_fields: Always empty for this view.
        error: Always ``None``: a file that cannot be read raises instead.
    """

    path: str
    file_type: str
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    row_count: int = 0
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None
    reader: Any = None
    scores_included: bool = False
    skipped_count: int | None = None
    skipped: list[SkippedMatch] | None = None
    rows_not_read: Any = None
    retention_time_unit: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "MatchRecords":
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            row_count=int(data.get("row_count", 0)),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=_output(data),
            scores_included=bool(data.get("scores_included", False)),
            skipped_count=data.get("skipped_count"),
            skipped=_skipped(data.get("skipped")),
            rows_not_read=data.get("rows_not_read"),
            retention_time_unit=data.get("retention_time_unit"),
            **_block(data),
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
        offset: Scans skipped: the offset that was applied, counted after ``ms_order``.
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
        source: What the file records about the run - instrument model and serial number, and
            when acquisition started (mzLib #1349) - as a :class:`SpectraSource`. ``None`` only
            when the reader built no description at all; a field the file does not record is
            ``None`` inside it.
        rows_not_read: Always ``None``: rows that did not become records are not counted for this
            view.
        absent_fields: Always empty for this view: a scan field a file does not record is ``None``
            per scan, as for the precursor of an MS1 scan.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some scans.
        excluded_fields: Always empty for this view.
        error: Always ``None``: a file that cannot be read raises instead.
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
    source: SpectraSource | None = None
    rows_not_read: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @property
    def total_ion_current(self) -> list[Any]:
        """The ``total_ion_current`` column, for the commonest plot there is."""
        columns = self._require_columns("read_spectra")
        return list(columns.get("total_ion_current") or [])

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "ScanRecords":
        block = _block(data)
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
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
            output=_output(data),
            source=SpectraSource._from_wire(data.get("source")),
            rows_not_read=data.get("rows_not_read"),
            **block,
        )


@dataclass(frozen=True)
class ProteinGroupRecords(_Table):
    """A MetaMorpheus protein-group table in **long** form: one row per group per sample group.

    What :func:`read_protein_groups` returns. Columns: ``protein_group_name``, ``gene``,
    ``organism``, ``decoy_contaminant_target``, ``is_decoy``, ``is_contaminant``,
    ``is_entrapment``, ``q_value``, then ``sample_label``, ``spectral_count`` and ``intensity``.

    **Unfiltered.** Decoys, contaminants and groups above 1% FDR are all rows, as MetaMorpheus
    writes them. Filter on ``q_value`` and ``decoy_contaminant_target`` before counting.

    Attributes:
        path: The absolute path that was read.
        file_type: ``"MetaMorpheusQuantifiedProteinGroups"``.
        reader: ``"ProteinGroupFromTsvFile"``.
        sample_labels: The file's sample-group labels, in header order, verbatim - e.g.
            ``"QE-002106_GM1_a-calib"``. Condition and replicate cannot be recovered from them;
            map them to your design yourself.
        record_count: Protein groups in the whole file - groups, not rows.
        returned_count: Groups carried back in :attr:`columns` - groups, not rows; the unit
            ``limit`` and ``offset`` count in. Zero when ``out`` was given.
        row_count: Rows in :attr:`columns`: returned groups times sample groups. Zero when ``out``
            was given.
        offset: Groups skipped. ``limit`` and ``offset`` count groups, not rows.
        truncated: Whether groups were left behind, by either ``limit`` or ``offset``.
        rows_not_read: Always ``None``: rows that did not become records are not counted for this
            table - every data line is a record, or the file fails to read.
        retention_time_unit: Always ``None``: this table has no time column.
        caveats: What this table cannot be trusted to mean, each citing the mzLib source.
        column_names: The field names, in order.
        columns: Field name -> list of values. ``intensity`` is ``None`` where the cell was blank
            (not quantified in that sample group), which is not zero. ``None`` when ``out`` was
            given.
        output: Where the table was written, or ``None``.
        absent_fields: Columns this file has no column for, so every value is ``None`` - ``gene``
            or ``organism`` when MetaMorpheus did not write them; ``sample_label``,
            ``spectral_count`` and ``intensity`` when the file has no per-sample columns, in which
            case each group is one row.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows.
        excluded_fields: ``count_occupancy`` and ``intensity_occupancy``, whose ``verb`` is
            ``readers read-occupancy`` - :func:`read_occupancy` has them, one row per site.
        error: Always ``None``: a file that cannot be read raises instead.
    """

    path: str
    file_type: str
    reader: Any
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    row_count: int = 0
    sample_labels: list[str] = field(default_factory=list)
    rows_not_read: Any = None
    retention_time_unit: Any = None
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "ProteinGroupRecords":
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            row_count=int(data.get("row_count", 0)),
            sample_labels=list(data.get("sample_labels") or []),
            rows_not_read=data.get("rows_not_read"),
            retention_time_unit=data.get("retention_time_unit"),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=_output(data),
            **_block(data),
        )


@dataclass(frozen=True)
class QuantifiedPeptideRecords(_Table):
    """A FlashLFQ peptide table in **long** form: one row per peptide per sample.

    What :func:`read_quantified_peptides` returns. Columns: ``sequence`` (the full, modified
    sequence), ``base_sequence``, ``peak_order``, ``protein_groups``, ``gene_names``,
    ``organism``, then ``sample_label``, ``intensity``, ``detection_type`` and ``retention_time``.

    **An intensity of 0 is not a measurement.** FlashLFQ writes a literal 0 for a peptide it did
    not quantify in a sample; ``detection_type`` (``"MSMS"``, ``"MBR"``, ``"NotDetected"``, ...)
    is what tells the two apart.

    Attributes:
        path: The absolute path that was read.
        file_type: ``"FlashLFQQuantifiedPeptide"``.
        reader: ``"QuantifiedPeptideFile"``.
        sample_labels: The file's sample labels, in header order, verbatim.
        record_count: Peptides in the whole file - peptides, not rows.
        returned_count: Peptides carried back in :attr:`columns` - peptides, not rows. Zero when
            ``out`` was given.
        row_count: Rows in :attr:`columns`: returned peptides times samples. Zero when ``out`` was
            given.
        offset: Peptides skipped. ``limit`` and ``offset`` count peptides, not rows.
        truncated: Whether peptides were left behind, by either ``limit`` or ``offset``.
        rows_not_read: Always ``None``: rows that did not become records are not counted for this
            table - every data line is a record, or the file fails to read.
        retention_time_unit: ``"minutes"``: the header says ``RetentionTime (min)``. The column is
            filled only by IsoTracker output; :attr:`absent_fields` names it otherwise.
        caveats: What this table cannot be trusted to mean, each citing the mzLib source.
        column_names: The field names, in order.
        columns: Field name -> list of values. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
        absent_fields: Columns this file has no column for, so every value is ``None`` -
            ``peak_order`` and ``retention_time`` for anything but IsoTracker output.
        failed_fields: ``"field: ExceptionType"`` for each column whose read threw on some rows.
        excluded_fields: Always empty for this table.
        error: Always ``None``: a file that cannot be read raises instead.
    """

    path: str
    file_type: str
    reader: Any
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    row_count: int = 0
    sample_labels: list[str] = field(default_factory=list)
    rows_not_read: Any = None
    retention_time_unit: Any = "minutes"
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "QuantifiedPeptideRecords":
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            row_count=int(data.get("row_count", 0)),
            sample_labels=list(data.get("sample_labels") or []),
            rows_not_read=data.get("rows_not_read"),
            retention_time_unit=data.get("retention_time_unit") or "minutes",
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=_output(data),
            **_block(data),
        )


@dataclass(frozen=True)
class OccupancyRecords(_Table):
    """PTM site occupancy from a MetaMorpheus protein-group table: one row per modified site.

    What :func:`read_occupancy` returns. Each row is one site in one occupancy cell - one group,
    one sample group, one ``basis``. Columns: ``protein_group_name``, ``sample_label``, ``basis``
    (``"count"`` or ``"intensity"``), ``entity_index``, ``position``, ``is_n_terminus``,
    ``modification``, ``fraction``, ``numerator``, ``denominator`` and ``cell_is_truncated``.

    **Which number to trust depends on the basis.** A ``"count"`` cell prints its fraction to two
    decimals, so use ``numerator / denominator`` (PSMs modified / PSMs covering the site); an
    ``"intensity"`` cell prints its numerator and denominator to four significant digits, so use
    ``fraction``.

    Attributes:
        path: The absolute path that was read.
        file_type: ``"MetaMorpheusQuantifiedProteinGroups"``.
        reader: ``"ProteinGroupFromTsvFile"``.
        sample_labels: The file's sample-group labels, in header order.
        truncated_cell_count: Occupancy cells the writer cut short or replaced with ``"Output too
            long for Excel"``. Their complete sites are rows with ``cell_is_truncated`` set; a cut
            cell with no complete site gives no rows and is counted only here.
        record_count: Protein groups in the whole file - groups, not rows.
        returned_count: Groups carried back - groups, not rows. Zero when ``out`` was given.
        row_count: Site rows in :attr:`columns`: one per group, sample group, basis and modified
            site. Zero when ``out`` was given.
        offset: Groups skipped. ``limit`` and ``offset`` count groups, not rows.
        truncated: Whether groups were left behind by ``limit`` or ``offset`` - not the same thing
            as a cell the writer cut short; see :attr:`truncated_cell_count`.
        rows_not_read: Always ``None``: rows that did not become records are not counted for this
            table - every data line is a record, or the file fails to read.
        retention_time_unit: Always ``None``: this table has no time column.
        caveats: What this table cannot be trusted to mean.
        column_names: The field names, in order.
        columns: Field name -> list of values. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None``.
        absent_fields: Every site column, when the file has no occupancy columns at all (then
            there are no rows).
        failed_fields: ``"count_occupancy: FormatException"`` (or ``intensity_occupancy``) when a
            cell was not an occupancy cell at all: mzLib refuses to guess, and that cell gives no
            rows.
        excluded_fields: Always empty for this table.
        error: Always ``None``: a file that cannot be read raises instead.
    """

    path: str
    file_type: str
    reader: Any
    record_count: int
    returned_count: int
    offset: int
    truncated: bool
    row_count: int = 0
    sample_labels: list[str] = field(default_factory=list)
    truncated_cell_count: int = 0
    rows_not_read: Any = None
    retention_time_unit: Any = None
    caveats: list[str] = field(default_factory=list)
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None
    absent_fields: list[str] = field(default_factory=list)
    failed_fields: list[str] = field(default_factory=list)
    excluded_fields: list[dict[str, Any]] = field(default_factory=list)
    error: ReadError | None = None

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "OccupancyRecords":
        return cls(
            path=data.get("path", ""),
            file_type=data.get("file_type", ""),
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            row_count=int(data.get("row_count", 0)),
            sample_labels=list(data.get("sample_labels") or []),
            truncated_cell_count=int(data.get("truncated_cell_count") or 0),
            rows_not_read=data.get("rows_not_read"),
            retention_time_unit=data.get("retention_time_unit"),
            caveats=list(data.get("caveats") or []),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=_output(data),
            **_block(data),
        )


@dataclass(frozen=True)
class ReadBatch(_Table):
    """Many files read into ONE long table: what every ``_many`` read function returns.

    The table's first two columns say where each row came from: ``source_index``, the file's
    0-based position in the list you passed, and ``source_path``, its absolute path. The rest are
    the single-file function's columns. Rows are grouped by file in input order, each file's rows
    in the file's own order - **whatever** ``threads`` **was**, so the same list always gives the
    same table. The per-file facts are in :attr:`files`, one :class:`FileReport` per input.

    Attributes:
        verb: The bridge command that produced it, e.g. ``"readers read-spectra"``.
        file_count: Files in the list you passed.
        read_count: Files read.
        failed_count: Files that could not be read - always 0 unless ``on_error="skip"``, since
            otherwise the first failure raises. Each one's :attr:`FileReport.error` says why.
        on_error: ``"fail"`` or ``"skip"``, as requested.
        record_count: Records read, summed over the files read, in the reader's own unit (scans,
            matches, groups, ...). Not always the row count: a long reader gives several rows per
            record.
        returned_count: Records carried back in :attr:`columns`, summed over the files, in the
            reader's own unit. Zero when ``out`` was given.
        row_count: Rows in :attr:`columns` - more than :attr:`returned_count` for a long reader
            (protein groups, peptides, occupancy, matches with scores). Zero when ``out`` was
            given.
        column_names: ``source_index``, ``source_path``, then the reader's columns, in order.
        columns: Column name -> list of values, one per row - the shape ``pandas.DataFrame`` and
            ``polars.DataFrame`` accept directly. ``None`` when ``out`` was given.
        output: Where the table was written, or ``None`` if it was returned inline. Written one
            file at a time, so memory held at most ``threads`` files.
        files: One :class:`FileReport` per input, in input order.
        ms_order: :func:`read_spectra_many` - the MS level every file was filtered to, or ``None``.
        peaks_included: :func:`read_spectra_many` - whether ``mz`` and ``intensity`` are present.
        scores_included: :func:`read_matches_many` - whether the table is long by score.
    """

    verb: str
    file_count: int
    read_count: int
    failed_count: int
    on_error: str
    record_count: int
    returned_count: int
    row_count: int = 0
    column_names: list[str] = field(default_factory=list)
    columns: Any = None
    output: Any = None
    files: list[FileReport] = field(default_factory=list)
    ms_order: Any = None
    peaks_included: bool = False
    scores_included: bool = False

    @property
    def failed_files(self) -> list[FileReport]:
        """The inputs that could not be read, under ``on_error="skip"``. Empty otherwise."""
        return [report for report in self.files if report.error is not None]

    def in_minutes(self, column: str) -> list[Any]:
        """A time column in minutes, converting each file's rows by that file's own unit.

        The multi-file form of :attr:`ResultRecords.retention_time_in_minutes`: every file's
        :attr:`FileReport.retention_time_unit` is applied to that file's rows, so a batch that mixes
        formats comes back on one axis.

        Args:
            column: The time column, e.g. ``"retention_time"`` or ``"retention_time_start"``.

        Returns:
            The column's values in minutes, one per row.

        Raises:
            UsageError: the table went to ``out``, has no such column, or a file's unit is
                ``"unknown"`` - raised rather than guessed.
        """
        return self._to_minutes(self.verb, column)

    @classmethod
    def _from_wire(cls, verb: str, data: dict[str, Any]) -> "ReadBatch":
        return cls(
            verb=verb,
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            on_error=data.get("on_error") or "fail",
            record_count=int(data.get("record_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            row_count=int(data.get("row_count", 0)),
            column_names=list(data.get("column_names") or []),
            columns=data.get("columns"),
            output=_output(data),
            files=[FileReport._from_wire(entry) for entry in (data.get("files") or [])],
            ms_order=data.get("ms_order"),
            peaks_included=bool(data.get("peaks_included", False)),
            scores_included=bool(data.get("scores_included", False)),
        )


@dataclass(frozen=True)
class IdentifyBatch:
    """What :func:`identify_many` returns: one :class:`FileInfo` per path, in the order given.

    Attributes:
        file_count: Paths given.
        read_count: Paths identified.
        failed_count: Paths that could not be identified - always 0 unless ``on_error="skip"``.
        on_error: ``"fail"`` or ``"skip"``, as requested.
        files: One :class:`FileInfo` per path, in input order; a failed one has
            :attr:`FileInfo.error` set.
    """

    file_count: int
    read_count: int
    failed_count: int
    on_error: str
    files: list[FileInfo] = field(default_factory=list)

    @property
    def failed_files(self) -> list[FileInfo]:
        """The paths that could not be identified, under ``on_error="skip"``."""
        return [info for info in self.files if info.error is not None]

    @classmethod
    def _from_wire(cls, data: dict[str, Any]) -> "IdentifyBatch":
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            on_error=data.get("on_error") or "fail",
            files=[FileInfo._from_wire(item) for item in (data.get("files") or [])],
        )


# ---- argument assembly --------------------------------------------------------------------------


def _window(
    verb: str,
    path: str | os.PathLike[str],
    *,
    limit: int | None,
    offset: int,
    out: str | os.PathLike[str] | None,
) -> list[str]:
    """The ``--path``/``--limit``/``--offset``/``--out`` argument list, validated.

    Shared by every single-file read so they cannot drift on what counts as a usable argument.
    Every check raises before the bridge is spawned: a caller who passed ``limit=-1`` wants an
    error, not a process launch.
    """
    if isinstance(path, (list, tuple, set, frozenset)):
        raise _bridge.UsageError(
            f"{verb.replace('-', '_')}() reads ONE file; for a list of {len(path)} paths call "
            f"{verb.replace('-', '_')}_many(), which reads them in one bridge process."
        )
    path = _bridge.path_text(path)
    if not path:
        raise _bridge.UsageError("A file path is required, e.g. 'AllPSMs.psmtsv'.")

    args = ["readers", verb, "--path", path]

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

    args += _out(out)
    return args


def _out(out: str | os.PathLike[str] | None) -> list[str]:
    if out is None:
        return []
    text = _bridge.path_text(out)
    if not text:
        raise _bridge.UsageError("out must be a non-empty path or None.")
    return ["--out", text]


def _batch(
    verb: str,
    paths: PathList,
    *,
    threads: int,
    on_error: str,
    out: str | os.PathLike[str] | None = None,
) -> tuple[list[str], str]:
    """The arguments and stdin of a ``_many`` read: the list goes on stdin, one path per line.

    The one implementation behind every ``_many`` function. **No loop over files happens here,
    and none may be added.** The whole list is handed to ONE bridge process, which reads
    ``threads`` files at once and returns one table. A Python-side thread or process pool would pay
    .NET start-up once per file and leave nobody owning the total thread count - see the bridge's
    ``design/PARALLELISM.md``.
    """
    if isinstance(paths, (str, bytes, bytearray, os.PathLike)) or not isinstance(paths, Iterable):
        raise _bridge.UsageError(
            f"{verb.split()[-1].replace('-', '_')}_many takes a LIST of paths; got "
            f"{type(paths).__name__}. For one file, call the function without _many."
        )
    texts = []
    for index, item in enumerate(paths):
        text = _bridge.path_text(item)
        if not text:
            raise _bridge.UsageError(
                f"Path {index} of the list is blank or not a path ({item!r}); every entry must be a "
                "str or pathlib.Path."
            )
        if "\n" in text or "\r" in text:
            # The list travels one path per line, so a line break inside a path would split it.
            raise _bridge.UsageError(f"Path {index} of the list contains a line break: {text!r}.")
        texts.append(text)
    if not texts:
        raise _bridge.UsageError("The list of paths is empty; give at least one file.")

    if isinstance(threads, bool) or not isinstance(threads, int) or (threads < 1 and threads != -1):
        raise _bridge.UsageError(
            f"threads must be a whole number of 1 or more, or -1 for one per core; got {threads!r}."
        )
    if on_error not in ("fail", "skip"):
        raise _bridge.UsageError(f"on_error must be 'fail' or 'skip'; got {on_error!r}.")

    args = ["readers", verb, "--paths-stdin", "--threads", str(threads), "--on-error", on_error]
    args += _out(out)
    return args, "\n".join(texts) + "\n"


def _read_many(
    verb: str,
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None,
    threads: int,
    on_error: str,
    timeout: float | None,
    extra: tuple[str, ...] = (),
) -> ReadBatch:
    args, stdin = _batch(verb, paths, threads=threads, on_error=on_error, out=out)
    data = _bridge.invoke(*args, *extra, stdin=stdin, timeout=timeout)
    return ReadBatch._from_wire(f"readers {verb}", data)


def _ms_order_args(ms_order: int | None, peaks: bool) -> tuple[str, ...]:
    extra: list[str] = []
    if ms_order is not None:
        if isinstance(ms_order, bool) or not isinstance(ms_order, int) or ms_order < 1:
            raise _bridge.UsageError(
                f"ms_order must be a whole number of 1 or more, or None; got {ms_order!r}."
            )
        extra += ["--ms-order", str(ms_order)]
    if peaks:
        extra += ["--peaks"]
    return tuple(extra)


# ---- identify -----------------------------------------------------------------------------------


def formats(timeout: float | None = 60) -> list[Format]:
    """Every file type mzLib can recognise.

    Enumerated from mzLib itself rather than from a list maintained here, so it reflects the
    installed version and cannot go stale.

    Args:
        timeout: Seconds to allow.

    Returns:
        One :class:`Format` per supported file type.

    Examples:
        >>> quantifiable = [f.file_type for f in formats() if f.is_quantifiable]
        >>> quantifiable
        ['psmtsv', 'osmtsv', 'MsFraggerPsm', 'DiaNnReport']
    """
    data = _bridge.invoke("readers", "formats", timeout=timeout)
    return [Format._from_wire(item) for item in (data.get("formats") or [])]


def identify(path: str | os.PathLike[str], timeout: float | None = 60) -> FileInfo:
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

    Examples:
        >>> info = identify("PXD078927_msgf_1_1_0.mzid")
        >>> info.file_type, info.views
        ('MzIdentML', ['spectral_match'])
    """
    path = _bridge.path_text(path)
    if not path:
        raise _bridge.UsageError("A file path is required, e.g. 'AllPSMs.psmtsv'.")

    data = _bridge.invoke("readers", "identify", "--path", path, timeout=timeout)
    return FileInfo._from_wire(data)


def identify_many(
    paths: PathList,
    *,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = 60,
) -> IdentifyBatch:
    """Identify many files in one call, without parsing their contents.

    :func:`identify` for a list: one bridge process, the answers in the order given. The cheap way
    to sort a directory of unknown files before reading any of them.

    Args:
        paths: The files to identify, as a list of ``str`` or ``pathlib.Path``.
        threads: Files identified at once: ``1`` (the default), more, or ``-1`` for one per core.
            The answer is the same at any value.
        on_error: ``"fail"`` (the default) raises on the first file that cannot be identified;
            ``"skip"`` answers it with :attr:`FileInfo.error` set and carries on.
        timeout: Seconds to allow for the whole list.

    Returns:
        An :class:`IdentifyBatch`, whose :attr:`IdentifyBatch.files` has one :class:`FileInfo` per
        path, in order.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = identify_many(["run.mzML", "AllPSMs.psmtsv", "notes.txt"],
        ...                       on_error="skip")                             # doctest: +SKIP
        >>> [f.file_type for f in batch.files]                                 # doctest: +SKIP
        ['MzML', 'psmtsv', '']
        >>> batch.failed_files[0].error.kind                                   # doctest: +SKIP
        'usage'
    """
    args, stdin = _batch("identify", paths, threads=threads, on_error=on_error)
    data = _bridge.invoke(*args, stdin=stdin, timeout=timeout)
    return IdentifyBatch._from_wire(data)


# ---- the table readers --------------------------------------------------------------------------


def read_results(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> ResultRecords:
    """Read a result file into the uniform record view.

    Only the four file types offering the ``"quantifiable"`` view can be read this way - check
    :func:`identify` first, or catch the error. A file without the view is rejected with a message
    naming the views it does have. :func:`read_results_many` reads a list of them into one table.

    **There is no default row limit.** A result file can carry a million rows, and truncating by
    default would mean the ordinary call returns a table that looks complete and is not. For a large
    file use ``out`` rather than paging: see the note on ``offset`` below.

    Args:
        path: Path to a MetaMorpheus ``.psmtsv`` / ``.osmtsv``, an MSFragger ``psm.tsv`` or a
            DIA-NN ``report.tsv``.
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

    Examples:
        >>> r = read_results("FraggerPsm_FragPipev21.1_psm.tsv", limit=2)
        >>> r.record_count, r.returned_count, r.absent_fields
        (5, 2, ['is_decoy'])

        pandas is not a dependency of pyMzLib, so these lines are not executed:

        >>> import pandas as pd                                      # doctest: +SKIP
        >>> pd.DataFrame(r.columns)                                  # doctest: +SKIP
    """
    args = _window("read-results", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return ResultRecords._from_wire(data)


def read_results_many(
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read many result files into one long table of the uniform record view.

    :func:`read_results` for a list, which may mix the four quantifiable formats: the columns are
    the same for all of them. Their units need not be - use :meth:`ReadBatch.in_minutes` for
    ``retention_time`` - and each file's caveats are in its :class:`FileReport`.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`ReadBatch` with the columns of :class:`ResultRecords`, after ``source_index`` and
        ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = read_results_many(["a/AllPSMs.psmtsv", "b/psm.tsv"], threads=2)   # doctest: +SKIP
        >>> [(f.file_type, f.record_count) for f in batch.files]                      # doctest: +SKIP
        [('psmtsv', 8), ('MsFraggerPsm', 6)]
    """
    return _read_many("read-results", paths, out=out, threads=threads, on_error=on_error, timeout=timeout)


def read_records(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> NativeRecords:
    """Read **any** file mzLib recognises, into that format's own fields.

    This is the exhaustive verb: if :func:`identify` succeeds on a path, this reads it. All
    thirty-six file types, including the seventeen that belong to no cross-format view at all -
    TopPIC, Crux, MSFragger's peptide and protein tables, the FlashDeconv formats, MetaMorpheus's
    protein-group and peptide tables - which no other
    function here can touch.

    **The columns are not uniform, by design.** They are this format's own mzLib record fields,
    under mzLib's own names in snake_case: a TopPIC file gives thirty-six columns, a Crux file
    twenty-three, an experiment annotation five. Read :attr:`NativeRecords.column_names`, and use
    :func:`read_results`, :func:`read_features` or :func:`read_matches` when you need columns that
    mean the same thing across formats.

    Nothing is silently dropped. A field that could not become a column is named in
    :attr:`NativeRecords.excluded_fields`, one that raised while being read in
    :attr:`NativeRecords.failed_fields`, and one this file's format has no column for in
    :attr:`NativeRecords.absent_fields` - so a missing value never has to be guessed at.

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

    Examples:
        >>> r = read_records("ToppicPrsm_TopPICv1.6.2_prsm.tsv")
        >>> r.record_type, len(r.column_names), r.record_count
        ('ToppicPrsm', 36, 4)

        pandas is not a dependency of pyMzLib, so this line is not executed:

        >>> import pandas as pd                                    # doctest: +SKIP
        >>> pd.DataFrame(r.columns)[["e_value", "q_value_spectrum_level"]]   # doctest: +SKIP
    """
    args = _window("read-records", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return NativeRecords._from_wire(data)


def read_records_many(
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read many files of ONE record type into one long table of their own fields.

    :func:`read_records` for a list. Because the columns are the record type's own, **every file
    must have the same record type** - all mzIdentML, all FlashLFQ peaks tables, all TopPIC PrSMs.
    A mixed list is refused with a ``UsageError`` before anything is parsed, whatever ``on_error``
    says, with the groups named so you can split it.
    Each file's ``record_type``, ``absent_fields`` and mzIdentML ``skipped`` list are in its
    :class:`FileReport`.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`ReadBatch` with the shared record type's columns after ``source_index`` and
        ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = read_records_many(sorted(Path("runs").glob("*.mzid.gz")), threads=4)   # doctest: +SKIP
        >>> batch.files[0].record_type, batch.record_count                                 # doctest: +SKIP
        ('MzIdentMLRecord', 48213)
    """
    return _read_many("read-records", paths, out=out, threads=threads, on_error=on_error, timeout=timeout)


def read_features(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
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

    Examples:
        >>> f = read_features("Ms1Feature_TopFDv1.6.2_ms1.feature", limit=5)
        >>> f.record_count, f.returned_count, f.retention_time_unit
        (25, 5, 'unknown')
        >>> f.absent_fields
        ['number_of_isotopes']
    """
    args = _window("read-features", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return FeatureRecords._from_wire(data)


def read_features_many(
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read many feature files into one long table of the ``ms1_features`` view.

    :func:`read_features` for a list. Units differ per file - ``_ms1.feature`` times are genuinely
    of unknown unit - so :meth:`ReadBatch.in_minutes` raises rather than put two units on one axis.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`ReadBatch` with the columns of :class:`FeatureRecords`, after ``source_index``
        and ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = read_features_many(["a.feature.tsv", "b.feature.tsv"])      # doctest: +SKIP
        >>> batch.file_count, batch.read_count                                  # doctest: +SKIP
        (2, 2)
    """
    return _read_many("read-features", paths, out=out, threads=threads, on_error=on_error, timeout=timeout)


def read_matches(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    scores: bool = False,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> MatchRecords:
    """Read identifications, in the cross-format ``spectral_match`` view.

    Six file types offer it: MsPathFinderT's targets, decoys and combined results, Casanovo's
    ``.mztab``, and mzIdentML (``.mzid``, and ``.mzid.gz`` read without unpacking it, mzLib #1313).
    These are the identification formats that share no *file*-level interface, so
    :func:`read_results` cannot reach them.

    **Nothing here is FDR-filtered.** ``q_value`` is the only confidence column, filled by
    mzIdentML (the PSM-level q-value, mzLib #1306) and by an MsPathFinderT ``_IcTda.tsv``; for
    every other file it is named in :attr:`MatchRecords.absent_fields` and ``None``. mzIdentML
    also lists every candidate, not only accepted ones: filter on ``rank == 1`` and
    ``pass_threshold`` before counting.

    Args:
        path: Path to an MsPathFinderT ``_IcTarget.tsv`` / ``_IcDecoy.tsv`` / ``_IcTda.tsv``, a
            Casanovo ``.mztab``, or an mzIdentML ``.mzid`` / ``.mzid.gz``.
        limit: Maximum matches to return. ``None`` returns all of them.
        offset: Matches to skip.
        scores: Add each search engine's own scores as **long** rows: one row per match and
            score, with ``match_index`` (the match's 0-based position among the file's matches),
            ``score_name`` (e.g. ``"MS-GF:SpecEValue"``) and ``score_value``. Only mzIdentML
            records them (mzLib #1306); other files keep one row per match and name the two
            columns in ``absent_fields``. ``limit`` and ``offset`` still count matches.
        out: Write a tab-separated table here and return only a summary.
        timeout: Seconds to allow.

    Returns:
        A :class:`MatchRecords`. Read :attr:`MatchRecords.caveats` before trusting ``is_decoy`` -
        it is inferred from a name prefix for MsPathFinderT and absent for Casanovo and mzIdentML.

    Raises:
        UsageError: the path is blank, missing, unrecognised, or has no ``spectral_match`` view.

    Examples:
        >>> m = read_matches("MsPathFinderT_WithMods_IcTda.tsv")
        >>> m.record_count, m.columns["modifications"][:3]
        (5, ['37:Oxidation on M', '', '1:Acetylation on X;27:Acetylation on K'])
        >>> s = read_matches("PXD078927_msgf_1_1_0.mzid", limit=1, scores=True)
        >>> s.columns["score_name"][:3], s.columns["q_value"][0]
        (['MS-GF:RawScore', 'MS-GF:DeNovoScore', 'MS-GF:SpecEValue'], 0)
    """
    args = _window("read-matches", path, limit=limit, offset=offset, out=out)
    if scores:
        args += ["--scores"]
    data = _bridge.invoke(*args, timeout=timeout)
    return MatchRecords._from_wire(data)


def read_matches_many(
    paths: PathList,
    *,
    scores: bool = False,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read many identification files into one long table of the ``spectral_match`` view.

    :func:`read_matches` for a list - a directory of mzIdentML submissions, say. Each file's
    ``skipped`` items and ``absent_fields`` are in its :class:`FileReport`.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.
        scores: As for :func:`read_matches`: one row per match and score, for the files that
            record scores (mzIdentML).

    Returns:
        A :class:`ReadBatch` with the columns of :class:`MatchRecords`, after ``source_index`` and
        ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = read_matches_many(["a.mzid", "b.mzid.gz"], scores=True)     # doctest: +SKIP
        >>> sorted(set(batch.columns["score_name"]))[:2]                        # doctest: +SKIP
        ['MS-GF:DeNovoScore', 'MS-GF:EValue']
    """
    return _read_many("read-matches", paths, out=out, threads=threads, on_error=on_error, timeout=timeout,
                      extra=("--scores",) if scores else ())


def read_spectra(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    ms_order: int | None = None,
    peaks: bool = False,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> ScanRecords:
    """Read the scans of a spectra file: headers always, peaks on request.

    Seven file types offer the ``spectra`` view: ``.mzML``, ``.mgf``, ``_ms1.msalign``,
    ``_ms2.msalign``, Thermo ``.raw``, Bruker ``.d`` and timsTOF ``.d``. Each read also reports
    what the file records about the run - instrument model and serial number, and when acquisition
    started (:attr:`ScanRecords.source`, mzLib #1349).

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

    Examples:
        >>> s = read_spectra("sliced_ethcd.mzML", ms_order=2, limit=2)
        >>> s.scan_count, s.record_count, s.returned_count
        (6, 5, 2)
        >>> s.columns["selected_ion_mz"]
        [548.453918457031, 796.765197753906]
        >>> s.source.instrument_model, s.source.acquisition_start_time
        ('Orbitrap Fusion', '2021-03-16T17:09:07Z')
    """
    args = _window("read-spectra", path, limit=limit, offset=offset, out=out)
    args += list(_ms_order_args(ms_order, peaks))
    data = _bridge.invoke(*args, timeout=timeout)
    return ScanRecords._from_wire(data)


def read_spectra_many(
    paths: PathList,
    *,
    ms_order: int | None = None,
    peaks: bool = False,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read the scans of many spectra files into one long table.

    :func:`read_spectra` for a list - every run of an experiment, with each file's instrument,
    serial number and acquisition start in its :attr:`FileReport.source`. Retention times are
    minutes for every format, so the table is on one time axis.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.
        ms_order: Keep only scans at this MS level, in every file.
        peaks: Include the ``mz`` and ``intensity`` arrays. Off by default, and worth keeping off
            for a list: with peaks, write the batch to ``out``.

    Returns:
        A :class:`ReadBatch` with the columns of :class:`ScanRecords`, after ``source_index`` and
        ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> runs = read_spectra_many(["run1.raw", "run2.raw"], ms_order=1, threads=2)   # doctest: +SKIP
        >>> [f.source.instrument_serial_number for f in runs.files]                   # doctest: +SKIP
        ['FSN20410', 'FSN20410']
    """
    return _read_many("read-spectra", paths, out=out, threads=threads, on_error=on_error, timeout=timeout,
                      extra=_ms_order_args(ms_order, peaks))


def read_protein_groups(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> ProteinGroupRecords:
    """Read a MetaMorpheus protein-group table as one row per group per sample group.

    ``AllQuantifiedProteinGroups.tsv`` keeps each group's per-sample spectral counts and
    intensities in columns named after the samples, which mzLib 1.0.592 reads into a dictionary
    (mzLib #1347) that :func:`read_records` cannot project. This is that dictionary as a **long**
    table - ``sample_label``, ``spectral_count``, ``intensity`` - beside the fields you filter on.
    A wide matrix is one ``pivot`` away; the readers guide shows it.

    **The table is unfiltered**, as MetaMorpheus writes it: decoys, contaminants and groups above
    1% FDR are all rows. The group's other fields - coverage, masses, member counts - are in
    :func:`read_records`, joined on ``protein_group_name``; PTM site occupancy is
    :func:`read_occupancy`.

    Args:
        path: Path to a MetaMorpheus ``AllQuantifiedProteinGroups.tsv`` (:func:`identify` reports
            ``MetaMorpheusQuantifiedProteinGroups``).
        limit: Maximum groups to return; each group gives one row per sample group. ``None``
            returns all of them.
        offset: Groups to skip.
        out: Write a tab-separated table here and return only a summary.
        timeout: Seconds to allow.

    Returns:
        A :class:`ProteinGroupRecords`.

    Raises:
        UsageError: the path is blank or missing, or is not a protein-group table; or the bridge
            in use predates this function.

    Examples:
        >>> g = read_protein_groups("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv", limit=1)
        >>> g.record_count, g.returned_count, len(g.sample_labels), g.row_count
        (6, 1, 18, 18)
        >>> g.records[0]["sample_label"], g.records[0]["intensity"]
        ('QE-002106_GM1_a-calib', 1838317.1674502683)
    """
    _bridge.require_verb("readers read-protein-groups", since=_QUANT_VIEWS_SINCE)
    args = _window("read-protein-groups", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return ProteinGroupRecords._from_wire(data)


def read_protein_groups_many(
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read many MetaMorpheus protein-group tables into one long table.

    :func:`read_protein_groups` for a list - one search per condition, say. Each file's sample
    labels are in its :attr:`FileReport.sample_labels`.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`ReadBatch` with the columns of :class:`ProteinGroupRecords`, after
        ``source_index`` and ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = read_protein_groups_many(["search1/AllQuantifiedProteinGroups.tsv",
        ...                                   "search2/AllQuantifiedProteinGroups.tsv"])   # doctest: +SKIP
        >>> batch.record_count, batch.row_count                                          # doctest: +SKIP
        (12, 216)
    """
    _bridge.require_verb("readers read-protein-groups", since=_QUANT_VIEWS_SINCE)
    return _read_many("read-protein-groups", paths, out=out, threads=threads, on_error=on_error, timeout=timeout)


def read_quantified_peptides(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> QuantifiedPeptideRecords:
    """Read a FlashLFQ peptide table as one row per peptide per sample.

    FlashLFQ's ``QuantifiedPeptides.tsv`` - and MetaMorpheus's ``AllQuantifiedPeptides.tsv``,
    which is the same format - keeps per-sample intensities and detection types in columns named
    after the samples, which mzLib 1.0.592 reads into a dictionary (mzLib #1347). This is that
    dictionary as a **long** table: ``sample_label``, ``intensity``, ``detection_type`` and, for
    IsoTracker output, ``retention_time``.

    **An intensity of 0 is not a measured zero.** FlashLFQ writes 0 for a peptide it did not
    quantify in a sample; filter on ``detection_type`` before taking a mean or a log.

    Args:
        path: Path to a FlashLFQ ``QuantifiedPeptides.tsv`` or MetaMorpheus
            ``AllQuantifiedPeptides.tsv`` (:func:`identify` reports ``FlashLFQQuantifiedPeptide``).
        limit: Maximum peptides to return; each gives one row per sample. ``None`` returns all.
        offset: Peptides to skip.
        out: Write a tab-separated table here and return only a summary.
        timeout: Seconds to allow.

    Returns:
        A :class:`QuantifiedPeptideRecords`.

    Raises:
        UsageError: the path is blank or missing, or is not a peptide table; or the bridge in use
            predates this function.

    Examples:
        >>> p = read_quantified_peptides("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", limit=1)
        >>> p.returned_count, p.row_count, p.absent_fields
        (1, 18, ['peak_order', 'retention_time'])
        >>> p.records[0]["intensity"], p.records[0]["detection_type"]
        (0, 'NotDetected')
    """
    _bridge.require_verb("readers read-quantified-peptides", since=_QUANT_VIEWS_SINCE)
    args = _window("read-quantified-peptides", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return QuantifiedPeptideRecords._from_wire(data)


def read_quantified_peptides_many(
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read many FlashLFQ peptide tables into one long table.

    :func:`read_quantified_peptides` for a list.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`ReadBatch` with the columns of :class:`QuantifiedPeptideRecords`, after
        ``source_index`` and ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> batch = read_quantified_peptides_many(["a/QuantifiedPeptides.tsv",
        ...                                        "b/QuantifiedPeptides.tsv"])   # doctest: +SKIP
        >>> batch.read_count                                                     # doctest: +SKIP
        2
    """
    _bridge.require_verb("readers read-quantified-peptides", since=_QUANT_VIEWS_SINCE)
    return _read_many("read-quantified-peptides", paths, out=out, threads=threads, on_error=on_error,
                      timeout=timeout)


def read_occupancy(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    out: str | os.PathLike[str] | None = None,
    timeout: float | None = None,
) -> OccupancyRecords:
    """Read the PTM site occupancy of a MetaMorpheus protein-group table, one row per site.

    MetaMorpheus writes two occupancy cells per group per sample group - one from PSM counts
    (``CountOccupancy_``) and one from intensities (``IntensityOccupancy_``) - each a list of
    modified sites encoded as text. mzLib 1.0.592 parses them (``ModificationOccupancyCell``,
    #1347); this returns every site of every cell as a row, with ``basis`` saying which cell it
    came from.

    Args:
        path: Path to a MetaMorpheus ``AllQuantifiedProteinGroups.tsv``.
        limit: Maximum groups to return. A group with no modified sites gives no rows.
        offset: Groups to skip.
        out: Write a tab-separated table here and return only a summary.
        timeout: Seconds to allow.

    Returns:
        An :class:`OccupancyRecords`. Check :attr:`OccupancyRecords.truncated_cell_count` and
        :attr:`OccupancyRecords.failed_fields` - a cut or malformed cell shortens the table.

    Raises:
        UsageError: the path is blank or missing, or is not a protein-group table; or the bridge
            in use predates this function.

    Examples:
        >>> o = read_occupancy("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
        >>> o.record_count, o.returned_count, o.row_count, o.truncated_cell_count
        (6, 6, 95, 0)
        >>> row = o.records[0]
        >>> row["basis"], row["position"], row["modification"], row["numerator"], row["denominator"]
        ('count', 329, 'Deamidation on N', 1, 2)
    """
    _bridge.require_verb("readers read-occupancy", since=_QUANT_VIEWS_SINCE)
    args = _window("read-occupancy", path, limit=limit, offset=offset, out=out)
    data = _bridge.invoke(*args, timeout=timeout)
    return OccupancyRecords._from_wire(data)


def read_occupancy_many(
    paths: PathList,
    *,
    out: str | os.PathLike[str] | None = None,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> ReadBatch:
    """Read the PTM site occupancy of many MetaMorpheus protein-group tables into one table.

    :func:`read_occupancy` for a list. Each file's ``truncated_cell_count`` is in its
    :class:`FileReport`.

    Args:
        paths: The files to read, as a list (or any iterable) of ``str`` or ``pathlib.Path``.
            Read by ONE bridge process into ONE table, in this order; a file's position is its
            ``source_index``. A repeated path is refused rather than read twice.
        out: Write the table to this path as **tab-separated** text and return only a summary
            (:attr:`ReadBatch.output`). Files are written one at a time, in order, so memory holds
            at most ``threads`` files however long the list is - the way to read hundreds. A
            batch that stops on an error removes its partial table.
        threads: Files to read at once: ``1`` (the default), more, or ``-1`` for one per core. The
            result is **byte-identical at any value** - files are always returned in input order -
            so this trades memory for speed and never changes an answer. The default is 1 because
            every reader holds a whole file in memory, so ``threads=8`` can mean eight whole files
            at once.
        on_error: ``"fail"`` (the default) stops at the first file that cannot be read and raises,
            naming it; ``"skip"`` records the failure in that file's :attr:`FileReport.error` and
            reads the rest.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        A :class:`ReadBatch` with the columns of :class:`OccupancyRecords`, after ``source_index``
        and ``source_path``.

    Raises:
        UsageError: ``paths`` is not a list, is empty, or holds a blank entry; a path is
            repeated; ``threads`` or ``on_error`` is not a value listed above; or, with
            ``on_error="fail"``, a file is missing or is not one this function reads - the
            message then starts ``Input <i> ('<path>')``.
        BridgeError: with ``on_error="fail"``, mzLib recognised a file and failed to parse it;
            the message names the input the same way.

    Examples:
        >>> # Not run in CI: the replay bridge answers single-path calls only.
        >>> sites = read_occupancy_many(["s1/AllQuantifiedProteinGroups.tsv",
        ...                              "s2/AllQuantifiedProteinGroups.tsv"])   # doctest: +SKIP
        >>> [f.truncated_cell_count for f in sites.files]                       # doctest: +SKIP
        [0, 0]
    """
    _bridge.require_verb("readers read-occupancy", since=_QUANT_VIEWS_SINCE)
    return _read_many("read-occupancy", paths, out=out, threads=threads, on_error=on_error, timeout=timeout)
