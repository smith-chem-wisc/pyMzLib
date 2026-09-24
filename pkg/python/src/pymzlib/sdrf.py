"""Read SDRF-Proteomics experimental-design files, and pool several into one table.

Every other reader in pyMzLib answers *what did the search find*. SDRF answers *what was
searched* - which sample, which organism part, which replicate, which instrument settings - and
that is the half you need to group results across experiments::

    >>> import pymzlib
    >>> doc = pymzlib.sdrf.read("PXD000070.sdrf.tsv")
    >>> doc.row_count, len(doc.columns)
    (6, 31)
    >>> doc.value("characteristics[organism]")               # doctest: +ELLIPSIS
    ['plasmodium falciparum', 'plasmodium falciparum', ...]

Pool several experiments into one analysis table, giving each a name you choose::

    >>> pooled = pymzlib.sdrf.pool({
    ...     "PXD000070.sdrf.tsv": "malaria",
    ...     "PXD026824.sdrf.tsv": "colon",
    ... }, limit=4)
    >>> pooled.document_count, pooled.row_count, pooled.labels
    (2, 24, ['malaria', 'colon'])

**This module is row-major, and every other reader here is columnar.** That is not a style
choice. :func:`pymzlib.readers.read_records` and friends hand back ``columns``, a
name-to-values dict, because their column names are a schema. SDRF's are *data*, and they
**repeat**: 649 files in the curated corpus carry ``comment[modification parameters]`` more than
once, up to eight times in one file, and one file repeats an *empty* name 23 times. A dict keyed
by name would silently keep one occurrence and drop the rest. So :attr:`SdrfDocument.columns` is
a **list that may contain duplicates**, :attr:`SdrfDocument.rows` is a list of cell lists, and
position is what links them. Use :meth:`SdrfDocument.value` for the first cell under a name and
:meth:`SdrfDocument.all` for every one.

**Three things worth knowing before you index anything:**

*Rows are ragged.* ``len(row)`` may be less than ``len(columns)``. Real files are like this -
PXD059974 in mzLib's own fixtures has a 46-column header with 17 of its 23 rows carrying 42
cells - and mzLib preserves it rather than padding, so the file round-trips byte for byte.
:meth:`SdrfDocument.value` returns ``None`` for a position a row does not reach.

*Cells are raw strings, never interpreted.* The SDRF key=value grammar
(``"NT=Oxidation;AC=UNIMOD:35"``) arrives exactly as written. It is not decoded, because it
cannot be told apart from a cell that merely contains ``=`` and ``;`` - ``comment[file uri]``
routinely carries pre-signed download URLs whose query strings contain ``Signature=`` and
``Expires=``.

*A reserved word is a real value.* ``"not available"`` and ``"not applicable"`` mean *the
experiment stated an absence*, which is not the same as a column the document does not have.
``None`` means the latter. Do not collapse the two.

**Three questions this module answers about a document, and one about ages.** Each is mzLib's
own answer, projected once in the bridge for all three bindings:

- :func:`validate` - *is this file well-formed?* Structural findings from mzLib's
  ``SdrfValidator``: required columns, column order and casing, cell shape, row-key uniqueness.
- :func:`lint` - *do these files write the same thing the same way?* Cross-document drift from
  ``SdrfDriftLint``: one accession under two names, a column that is a CV term in one file and
  free text in another, ``"Homo sapiens"`` against ``"homo sapiens"``.
- :func:`assess` - *does this file describe its samples at all?* ``SdrfSampleInformativeness``'s
  verdict - ``Informative``, ``Partial`` or ``Skeleton`` - with the counts behind it.
- :func:`samples` - *what does it say about each sample?* ``SdrfSampleBlock``: one sample per
  ``source name``, its characteristics and factor values, and ages read into years.
- :func:`parse_ages` - ``SdrfAge.TryParse`` over any age cells, SDRF or not.

The first three are blind in different places, which is why there are three. A file of
``"not available"`` validates cleanly and lints clean, and only :func:`assess` sees that it says
nothing. Twenty individually valid files can still disagree, and only :func:`lint` looks across
them. Each result's ``caveats`` names the blind spot it leaves.

:func:`validate`, :func:`assess` and :func:`samples` each have a ``*_many`` form that reads a whole
corpus in one bridge call, with an explicit ``threads`` count - the batch is parallelised inside
the bridge, never by a Python loop.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from . import _bridge
from .readers import _Table

__all__ = [
    "SdrfDocument",
    "PooledSdrf",
    "WrittenSdrf",
    "read",
    "pool",
    "FileError",
    "ValidationMessage",
    "SdrfValidation",
    "ValidatedFile",
    "SdrfValidationBatch",
    "validate",
    "validate_many",
    "SdrfDrift",
    "lint",
    "SdrfAssessment",
    "AssessedFile",
    "SdrfAssessmentBatch",
    "assess",
    "assess_many",
    "SdrfSamples",
    "SampledFile",
    "SdrfSamplesBatch",
    "samples",
    "samples_many",
    "ParsedAges",
    "parse_ages",
    "VERDICTS",
    "AGE_PRECISIONS",
    "AGE_REFUSALS",
]


@dataclass(frozen=True)
class WrittenSdrf:
    """Where :func:`pool` wrote the merged document, when asked to write one.

    Attributes:
        path: The path written.
        row_count: Rows written - the **whole** merged document, not the windowed slice. ``limit``
            and ``offset`` shape what comes back over the wire; they never shorten the file.
    """

    path: str
    row_count: int

    @classmethod
    def _from_wire(cls, payload: Mapping[str, Any]) -> "WrittenSdrf":
        return cls(path=payload.get("path", ""), row_count=int(payload.get("row_count", 0)))


@dataclass(frozen=True)
class SdrfDocument:
    """One SDRF-Proteomics document: an ordered header, and rows of raw cells.

    Attributes:
        path: The path that was read.
        columns: The column names, **verbatim and in document order**. Names may repeat, and the
            order is part of the document, so this is a list rather than a set. Names are never
            case-normalised: the corpus contains ``"comment[MS min charge]"`` and
            ``"Material Type"``, and rewriting them would silently alter someone's file.
        rows: One list of cells per row. **Ragged**: a row may be shorter than ``columns``.
        row_count: Rows in the **whole document**, regardless of ``limit`` or ``offset``.
        returned_count: Rows actually carried back in :attr:`rows`.
        offset: The offset that was applied.
        truncated: Whether rows were left behind, by either ``limit`` or ``offset``. A short answer
            and a complete one must never look alike, so check this rather than assuming.
        caveats: What this document's data cannot tell you about itself - raggedness, repeated
            names, reserved words. Worth printing the first time you read an unfamiliar file.
    """

    path: str
    columns: list[str]
    rows: list[list[str]]
    row_count: int
    returned_count: int
    offset: int
    truncated: bool
    caveats: list[str] = field(default_factory=list)

    def index_of(self, column: str) -> int:
        """The position of the first column with this name, or ``-1`` if absent.

        Comparison is exact and case-**sensitive**, matching the SDRF specification and mzLib.
        """
        try:
            return self.columns.index(column)
        except ValueError:
            return -1

    def indexes_of(self, column: str) -> list[int]:
        """Every position carrying this name, in document order. Empty when the column is absent."""
        return [i for i, name in enumerate(self.columns) if name == column]

    def value(self, column: str) -> list[str | None]:
        """The first cell under ``column``, one entry per returned row.

        ``None`` means **the document does not have this column, or the row is too short to reach
        it**. It does not mean "empty": the SDRF reserved words ``"not available"`` and
        ``"not applicable"`` are real values that an experiment chose to write, and they come back
        as themselves.

        Examples:
            >>> doc = read("PXD000070.sdrf.tsv")
            >>> doc.value("characteristics[disease]")     # doctest: +ELLIPSIS
            ['not applicable', 'not applicable', ...]
        """
        i = self.index_of(column)
        if i < 0:
            return [None] * len(self.rows)
        return [row[i] if i < len(row) else None for row in self.rows]

    def all(self, column: str) -> list[list[str]]:
        """Every cell under ``column``, one list per returned row.

        The accessor for a multi-cardinality column such as ``comment[modification parameters]``,
        which legitimately repeats - up to eight times in one corpus file. Positions a row is too
        short to reach are skipped rather than reported as ``None``, so each inner list holds only
        cells that exist.
        """
        indexes = self.indexes_of(column)
        return [[row[i] for i in indexes if i < len(row)] for row in self.rows]

    @property
    def records(self) -> list[dict[str, str]]:
        """The rows as dicts, for the common case of a document with no repeated column names.

        **Lossy when a name repeats** - later positions overwrite earlier ones - which is exactly
        why it is not the primary shape. :attr:`has_repeated_columns` says whether that applies to
        this document. For a repeating document use :meth:`all` instead.
        """
        return [
            {name: row[i] for i, name in enumerate(self.columns) if i < len(row)}
            for row in self.rows
        ]

    @property
    def has_repeated_columns(self) -> bool:
        """Whether any column name appears more than once - see :attr:`records`."""
        return len(self.columns) != len(set(self.columns))

    @property
    def ragged_row_count(self) -> int:
        """How many returned rows carry fewer cells than there are columns."""
        return sum(1 for row in self.rows if len(row) < len(self.columns))

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfDocument":
        return cls(
            path=data.get("path", ""),
            columns=list(data.get("column_names") or []),
            rows=[list(r) for r in (data.get("rows") or [])],
            row_count=int(data.get("row_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class PooledSdrf(SdrfDocument):
    """Several SDRF documents merged into one table.

    Everything :class:`SdrfDocument` offers, plus where the rows came from. ``path`` is empty:
    a pooled table is not a file that was read.

    Attributes:
        document_count: How many documents were pooled.
        paths: The paths pooled, in the order given.
        labels: The provenance label used for each, in the same order - either what you supplied
            or mzLib's ``containing-folder/file-stem`` default.
        written: Where the merged document was written, when ``out`` was given. ``None`` otherwise.
    """

    document_count: int = 0
    paths: list[str] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)
    written: WrittenSdrf | None = None

    #: The column :meth:`pool` adds to record which document each row came from.
    SOURCE_DOCUMENT_COLUMN = "comment[source document]"

    def source_documents(self) -> list[str | None]:
        """The provenance label of each returned row - which document it came from."""
        return self.value(self.SOURCE_DOCUMENT_COLUMN)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "PooledSdrf":
        written = data.get("written")
        return cls(
            path="",
            columns=list(data.get("column_names") or []),
            rows=[list(r) for r in (data.get("rows") or [])],
            row_count=int(data.get("row_count", 0)),
            returned_count=int(data.get("returned_count", 0)),
            offset=int(data.get("offset", 0)),
            truncated=bool(data.get("truncated", False)),
            caveats=list(data.get("caveats") or []),
            document_count=int(data.get("document_count", 0)),
            paths=list(data.get("paths") or []),
            labels=list(data.get("labels") or []),
            written=WrittenSdrf._from_wire(written) if written else None,
        )


def _window(limit: int | None, offset: int) -> list[str]:
    """Validates and renders the shared limit/offset options. Mirrors :mod:`pymzlib.readers`."""
    args: list[str] = []
    if limit is not None:
        if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
            raise _bridge.UsageError(
                f"limit must be a non-negative whole number or None; got {limit!r}."
            )
        args += ["--limit", str(limit)]
    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise _bridge.UsageError(f"offset must be a non-negative whole number; got {offset!r}.")
    if offset:
        args += ["--offset", str(offset)]
    return args


def read(
    path: str | os.PathLike[str],
    *,
    limit: int | None = None,
    offset: int = 0,
    timeout: float | None = 60,
) -> SdrfDocument:
    """Read one SDRF-Proteomics file.

    Args:
        path: Path to a ``.sdrf.tsv`` file.
        limit: Maximum rows to return. ``None`` (the default) returns all of them.
        offset: Rows to skip.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`SdrfDocument`.

    Raises:
        UsageError: the path is blank, or the file is missing or unreadable as SDRF.

    Examples:
        >>> doc = read("PXD000070.sdrf.tsv")
        >>> doc.value("characteristics[organism part]")[0]
        'human erythrocytes'
        >>> doc.all("comment[modification parameters]")[0]     # doctest: +ELLIPSIS
        ['NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed', 'NT=Oxidation;AC=UNIMOD:35;...]
    """
    path = _bridge.path_text(path)
    if not path:
        raise _bridge.UsageError("A file path is required, e.g. 'PXD000070.sdrf.tsv'.")

    args = ["sdrf", "read", "--path", path, *_window(limit, offset)]
    data = _bridge.invoke(*args, timeout=timeout)
    return SdrfDocument._from_wire(data)


def pool(
    documents: Sequence[str] | Mapping[str, str],
    *,
    out: str | os.PathLike[str] | None = None,
    limit: int | None = None,
    offset: int = 0,
    timeout: float | None = 60,
) -> PooledSdrf:
    """Merge several SDRF documents into one analysis table.

    Columns are the union of every document's, ordered by SDRF's own block structure, and a name
    that repeats is carried at the highest multiplicity any single document used, so nothing is
    dropped. A cell a document did not have is filled with the reserved word ``"not available"``,
    and a ``comment[source document]`` column records which document each row came from.

    **Give your documents names.** Pass a mapping to choose them. With a plain sequence, mzLib
    falls back to ``containing-folder/file-stem``, which depends on where the files happen to sit -
    so the same two documents pooled from a different directory produce a different table. The
    returned :attr:`PooledSdrf.caveats` says so when that fallback was used.

    **The result is an analysis table, not something to deposit.** ``source name`` + ``assay
    name`` + ``comment[label]`` is unique within one document, but two experiments may both have a
    ``"Sample 1"``, so a pooled table will usually violate SDRF's uniqueness rule. Use the
    source-document column as part of any key.

    Args:
        documents: The files to pool - a sequence of paths, or a ``{path: label}`` mapping.
        out: Write the merged document here as SDRF. The **whole** document is written regardless
            of ``limit`` and ``offset``.
        limit: Maximum rows to return over the wire. ``None`` (the default) returns all of them.
        offset: Rows to skip.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`PooledSdrf`.

    Raises:
        UsageError: no documents were given, a path is missing, or a label is blank.

    Examples:
        >>> pooled = pool({"PXD000070.sdrf.tsv": "malaria", "PXD026824.sdrf.tsv": "colon"}, limit=4)
        >>> pooled.returned_count, pooled.truncated, set(pooled.source_documents())
        (4, True, {'malaria'})
    """
    lines = _document_lines(documents, verb="pool")

    args = ["sdrf", "pool", *_window(limit, offset)]
    if out is not None:
        out = _bridge.path_text(out)
        if not out:
            raise _bridge.UsageError("out must be a non-empty path or None.")
        args += ["--out", out]

    data = _bridge.invoke(*args, stdin="\n".join(lines), timeout=timeout)
    return PooledSdrf._from_wire(data)



def _document_lines(documents: Sequence[str] | Mapping[str, str], *, verb: str) -> list[str]:
    """Renders ``documents`` as the ``path[<TAB>label]`` stdin lines ``pool`` and ``lint`` share."""
    if isinstance(documents, Mapping):
        pairs = [(str(p), str(label)) for p, label in documents.items()]
    elif isinstance(documents, str):
        # A bare string is a sequence of characters, so this would otherwise send one document per
        # letter and fail with a baffling "file not found: P".
        raise _bridge.UsageError(
            f"{verb}() takes several documents; pass a list or a {{path: label}} mapping. "
            f"For one document use read({documents!r})."
        )
    elif isinstance(documents, Sequence):
        pairs = [(str(p), "") for p in documents]
    else:
        raise _bridge.UsageError(
            f"documents must be a sequence of paths or a {{path: label}} mapping; got {documents!r}."
        )

    if not pairs:
        raise _bridge.UsageError("At least one SDRF document is required.")

    lines = []
    for path, label in pairs:
        if not path.strip():
            raise _bridge.UsageError("A document path may not be blank.")
        if isinstance(documents, Mapping) and not label.strip():
            raise _bridge.UsageError(
                f"The label for '{path}' is blank. Give every document a label, or pass a plain "
                "list to accept mzLib's path-derived default for all of them."
            )
        if "\t" in path or "\t" in label or "\n" in path or "\n" in label:
            raise _bridge.UsageError(
                f"A path or label contains a tab or newline, which the bridge uses to separate "
                f"them: {path!r} / {label!r}."
            )
        lines.append(f"{path.strip()}\t{label.strip()}" if label.strip() else path.strip())
    return lines


# ---- shared by the per-document verbs ------------------------------------------------------------

#: The three verdicts :func:`assess` can return, best first. Enum member names, as mzLib spells them.
VERDICTS: tuple[str, ...] = ("Informative", "Partial", "Skeleton")

#: How much of an age a cell pins down - mzLib's ``SdrfAgePrecision``. See :func:`parse_ages`.
AGE_PRECISIONS: tuple[str, ...] = ("Exact", "Range", "LowerBound", "UpperBound")

#: Why an age cell was refused. mzLib's ``TryParse`` answers only yes or no; the bridge names the
#: reason from the cell: ``empty``; ``reserved_word`` (``"not available"`` and friends, any case);
#: ``no_unit`` (a bare number such as ``63`` - years or days cannot be told apart); ``unreadable``.
AGE_REFUSALS: tuple[str, ...] = ("empty", "reserved_word", "no_unit", "unreadable")


@dataclass(frozen=True)
class FileError:
    """Why one document in a ``*_many`` call produced no result, under ``on_error="skip"``.

    Attributes:
        kind: ``"usage"`` (the file is missing) or ``"correctness"`` (mzLib could not read it -
            ``message`` then starts with the .NET exception type, e.g. ``"MzLibException: ..."``).
        message: What went wrong, naming the file.
    """

    kind: str
    message: str

    @classmethod
    def _from_wire(cls, payload: Mapping[str, Any] | None) -> "FileError | None":
        if not payload:
            return None
        return cls(kind=payload.get("kind", ""), message=payload.get("message", ""))


def _one_path(path: str | os.PathLike[str]) -> str:
    text = _bridge.path_text(path)
    if not text:
        raise _bridge.UsageError("A file path is required, e.g. 'PXD000070.sdrf.tsv'.")
    return text


def _bulk_args(
    paths: Sequence[str | os.PathLike[str]], threads: int, on_error: str, *, verb: str
) -> tuple[list[str], str]:
    """Validates a ``*_many`` call and renders its options and its stdin."""
    if isinstance(paths, (str, os.PathLike)):
        raise _bridge.UsageError(
            f"{verb}_many() takes a list of paths. For one document use {verb}({paths!r})."
        )
    if not isinstance(paths, Sequence):
        raise _bridge.UsageError(f"paths must be a list of file paths; got {paths!r}.")
    if not paths:
        raise _bridge.UsageError("At least one SDRF document is required.")
    lines = []
    for p in paths:
        text = _bridge.path_text(p)
        if not text:
            raise _bridge.UsageError(f"Every path must be a non-blank str or PathLike; got {p!r}.")
        if "\n" in text or "\r" in text:
            raise _bridge.UsageError(
                f"A path contains a newline, which separates paths on the wire: {text!r}."
            )
        lines.append(text)
    if (
        isinstance(threads, bool)
        or not isinstance(threads, int)
        or not (threads >= 1 or threads == -1)
    ):
        raise _bridge.UsageError(
            f"threads must be a positive whole number, or -1 for every core; got {threads!r}."
        )
    if on_error not in ("fail", "skip"):
        raise _bridge.UsageError(f"on_error must be 'fail' or 'skip'; got {on_error!r}.")
    args = ["--paths-stdin", "--threads", str(threads), "--on-error", on_error]
    return args, "\n".join(lines)


def _columns(data: Mapping[str, Any]) -> tuple[list[str], dict[str, list[Any]]]:
    names = list(data.get("column_names") or [])
    columns = {k: list(v) for k, v in (data.get("columns") or {}).items()}
    return names, columns


# ---- validate ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ValidationMessage:
    """One finding from mzLib's ``SdrfValidator``, located as precisely as its rule allows.

    Attributes:
        severity: ``"Error"`` (the document cannot be reliably consumed - a ragged row, a missing
            required column, two indistinguishable rows) or ``"Warning"`` (it deviates from the
            specification but can still be joined - casing, an empty cell, a missing recommended
            column).
        rule: The stable rule name, e.g. ``"RequiredColumn"``, ``"RowWidth"``,
            ``"RowKeyUniqueness"``, ``"ReservedWordCase"``. Stable so that you can count or
            suppress by it.
        message: A human-readable description including the offending value.
        row_index: 0-based index into the document's data rows, or ``None`` for a finding about
            the whole document (a missing column, an empty header).
        line_number: The 1-based line in the file - ``row_index + 2``, since the header is line 1.
            ``None`` exactly when ``row_index`` is.
        column_name: The column involved, or ``None`` when the finding is not about one column.
    """

    severity: str
    rule: str
    message: str
    row_index: int | None
    line_number: int | None
    column_name: str | None


def _messages(columns: Mapping[str, list[Any]]) -> list[ValidationMessage]:
    return [
        ValidationMessage(
            severity=columns["severity"][i],
            rule=columns["rule"][i],
            message=columns["message"][i],
            row_index=columns["row_index"][i],
            line_number=columns["line_number"][i],
            column_name=columns["column_name"][i],
        )
        for i in range(len(columns.get("rule", [])))
    ]


@dataclass(frozen=True)
class SdrfValidation(_Table):
    """The structural findings for one SDRF document.

    The findings are a table - ``columns`` maps ``severity``, ``rule``, ``message``, ``row_index``,
    ``line_number`` and ``column_name`` to one value per finding, ready for
    ``pandas.DataFrame(result.columns)``. :attr:`messages` gives the same rows as objects.

    Attributes:
        path: The absolute path validated.
        is_valid: ``True`` when there is no ``Error``. **Warnings never make a document invalid**:
            mzLib calibrated every severity against the 1,236-file curated corpus, and a rule that
            fired on most curated files was judged wrong rather than the files.
        error_count: Findings of severity ``Error``.
        warning_count: Findings of severity ``Warning``.
        message_count: All findings; the length of every column.
        row_count: Data rows in the document (the header is not a row).
        column_names: The table's columns, in order.
        columns: Column name to one value per finding.
        caveats: What a clean result does *not* mean - read these once.
    """

    path: str
    is_valid: bool
    error_count: int
    warning_count: int
    message_count: int
    row_count: int
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    @property
    def messages(self) -> list[ValidationMessage]:
        """Every finding as a :class:`ValidationMessage`, in mzLib's order."""
        return _messages(self.columns)

    @property
    def errors(self) -> list[ValidationMessage]:
        """The findings that make the document invalid."""
        return [m for m in self.messages if m.severity == "Error"]

    @property
    def warnings(self) -> list[ValidationMessage]:
        """The findings worth fixing that do not make the document invalid."""
        return [m for m in self.messages if m.severity == "Warning"]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfValidation":
        names, columns = _columns(data)
        return cls(
            path=data.get("path", ""),
            is_valid=bool(data.get("is_valid", False)),
            error_count=int(data.get("error_count", 0)),
            warning_count=int(data.get("warning_count", 0)),
            message_count=int(data.get("message_count", 0)),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class ValidatedFile:
    """One document's summary in a :func:`validate_many` result.

    Every count is ``None`` when the document could not be read, and :attr:`error` then says why.

    Attributes:
        path: The absolute path.
        is_valid: No ``Error`` findings; ``None`` if unread.
        error_count: ``Error`` findings; ``None`` if unread.
        warning_count: ``Warning`` findings; ``None`` if unread.
        message_count: All findings; ``None`` if unread.
        row_count: Data rows in the document; ``None`` if unread.
        error: Why it was not read, or ``None`` when it was.
    """

    path: str
    is_valid: bool | None
    error_count: int | None
    warning_count: int | None
    message_count: int | None
    row_count: int | None
    error: FileError | None

    @classmethod
    def _from_wire(cls, f: Mapping[str, Any]) -> "ValidatedFile":
        return cls(
            path=f.get("path", ""),
            is_valid=f.get("is_valid"),
            error_count=f.get("error_count"),
            warning_count=f.get("warning_count"),
            message_count=f.get("message_count"),
            row_count=f.get("row_count"),
            error=FileError._from_wire(f.get("error")),
        )


@dataclass(frozen=True)
class SdrfValidationBatch(_Table):
    """The findings for many SDRF documents, as one long table.

    ``columns`` holds the same six columns as :class:`SdrfValidation`, preceded by
    ``source_index`` (the document's 0-based position in the list you passed) and ``source_path``.
    Rows are in input order, then mzLib's order within a document, whatever ``threads`` was.

    Attributes:
        file_count: Documents given.
        read_count: Documents validated.
        failed_count: Documents that could not be read (only non-zero under ``on_error="skip"``).
        valid_count: Documents read with no ``Error``.
        message_count: Findings over every document read.
        files: One :class:`ValidatedFile` per input, in input order.
        column_names: The table's columns, ``source_index`` and ``source_path`` first.
        columns: Column name to one value per finding.
        caveats: As for :class:`SdrfValidation`.
    """

    file_count: int
    read_count: int
    failed_count: int
    valid_count: int
    message_count: int
    files: list[ValidatedFile]
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    @property
    def messages(self) -> list[ValidationMessage]:
        """Every finding as a :class:`ValidationMessage`; use ``columns`` to keep the source."""
        return _messages(self.columns)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfValidationBatch":
        names, columns = _columns(data)
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            valid_count=int(data.get("valid_count", 0)),
            message_count=int(data.get("message_count", 0)),
            files=[ValidatedFile._from_wire(f) for f in (data.get("files") or [])],
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


def validate(path: str | os.PathLike[str], *, timeout: float | None = 60) -> SdrfValidation:
    """Check one SDRF file against the SDRF-Proteomics specification's structural rules.

    Calls mzLib's ``SdrfValidator.Validate``. Structure only - required and recommended columns,
    column order, casing, malformed names, ragged rows, integer replicate and fraction columns,
    reserved-word casing, and the one hard row rule: ``source name`` + ``assay name`` +
    ``comment[label]`` must be unique. Controlled-vocabulary accessions are **not** resolved.

    Args:
        path: Path to a ``.sdrf.tsv`` file.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`SdrfValidation`: ``is_valid``, the counts, and one row per finding.

    Raises:
        UsageError: the path is blank or the file does not exist.
        BridgeError: mzLib could not read the file at all (for example, it is empty).

    Example:
        >>> result = validate("PXD000070.sdrf.tsv")            # doctest: +SKIP
        >>> result.is_valid, result.error_count                # doctest: +SKIP
        (True, 0)
        >>> for m in result.warnings:                          # doctest: +SKIP
        ...     print(m.line_number, m.rule, m.column_name)
    """
    data = _bridge.invoke("sdrf", "validate", "--path", _one_path(path), timeout=timeout)
    return SdrfValidation._from_wire(data)


def validate_many(
    paths: Sequence[str | os.PathLike[str]],
    *,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> SdrfValidationBatch:
    """Validate many SDRF files in one bridge call.

    One process for the whole corpus instead of one per file - roughly 120 ms saved per document -
    with the parallelism inside the bridge, where it belongs. The result does not depend on
    ``threads``: rows are always in input order.

    Args:
        paths: The files, in the order you want them reported. Each may be given once.
        threads: Documents validated at once. ``1`` (the default) keeps one document in memory at
            a time; ``-1`` uses every core.
        on_error: ``"fail"`` (the default) raises on the first unreadable document in input order.
            ``"skip"`` records the failure in that document's :attr:`ValidatedFile.error` and
            carries on.
        timeout: Seconds to allow for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        An :class:`SdrfValidationBatch`.

    Raises:
        UsageError: no paths, a blank or repeated path, a bad ``threads`` or ``on_error``, or -
            under ``"fail"`` - a file that does not exist.
        BridgeError: under ``"fail"``, mzLib could not read one of the files.

    Example:
        >>> paths = sorted(Path("corpus").glob("*.sdrf.tsv"))                 # doctest: +SKIP
        >>> batch = validate_many(paths, threads=-1, on_error="skip")          # doctest: +SKIP
        >>> batch.valid_count, batch.read_count                                # doctest: +SKIP
    """
    args, stdin = _bulk_args(paths, threads, on_error, verb="validate")
    data = _bridge.invoke("sdrf", "validate", *args, stdin=stdin, timeout=timeout)
    return SdrfValidationBatch._from_wire(data)


# ---- lint ----------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SdrfDrift(_Table):
    """Concepts a set of SDRF documents wrote inconsistently, from mzLib's ``SdrfDriftLint``.

    ``columns`` is a long table, **one row per finding x variant**: a finding with three spellings
    is three rows sharing a ``finding_index``. Its columns:

    - ``finding_index`` - which finding the row belongs to, 0-based, most impactful first.
    - ``kind`` - ``AccessionNameConflict`` (one accession, several names), ``NameAccessionConflict``
      (one name, several accessions - the serious one), ``MixedTermAndFreeText``,
      ``ColumnNameVariant`` (names differing only by case or spacing), ``ValueCaseVariant``.
    - ``concept`` - what was written inconsistently: an accession, a name, a normalised value, or a
      normalised column name, depending on ``kind``.
    - ``column_name`` - the column involved; ``None`` for ``ColumnNameVariant``, whose finding is
      about a name rather than a column's values.
    - ``variant_rank`` - 0 for the majority spelling, then by frequency.
    - ``value`` - this spelling exactly as written.
    - ``occurrences`` - how widely it was used (a count of documents, at mzLib 1.0.592).
    - ``documents`` - the labels of the documents that used it, as a list.

    Attributes:
        document_count: Documents linted.
        paths: The paths, in the order given.
        labels: The label each document is named by in ``documents`` - yours, or mzLib's
            ``containing-folder/file-stem`` default.
        finding_count: Distinct findings (not rows).
        column_names: The table's columns.
        columns: Column name to one value per finding x variant.
        caveats: Including why the majority spelling is not advice.
    """

    document_count: int
    paths: list[str]
    labels: list[str]
    finding_count: int
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    def findings(self) -> list[list[dict[str, Any]]]:
        """The rows grouped by finding: one list of variant dicts per finding, majority first."""
        grouped: dict[int, list[dict[str, Any]]] = {}
        for record in self.records:
            grouped.setdefault(record["finding_index"], []).append(record)
        return [grouped[i] for i in sorted(grouped)]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfDrift":
        names, columns = _columns(data)
        return cls(
            document_count=int(data.get("document_count", 0)),
            paths=list(data.get("paths") or []),
            labels=list(data.get("labels") or []),
            finding_count=int(data.get("finding_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


def lint(
    documents: Sequence[str] | Mapping[str, str], *, timeout: float | None = 60
) -> SdrfDrift:
    """Find the concepts a set of SDRF documents annotated inconsistently.

    Validity is a property of one file; comparability is a property of the relationship between
    files. Every document can pass :func:`validate` and the pooled table still be unusable because
    one file writes ``"Homo sapiens"`` and another ``"homo sapiens"``. This is mzLib's check for that.

    Columns that are unique per row by construction are never compared - ``source name``,
    ``assay name``, ``comment[data file]``, ``comment[searched data file]`` (since mzLib 1.0.592),
    ``comment[file uri]`` and ``comment[source document]`` - so differing file names are not
    reported as drift. Reserved words and empty cells are skipped too, so a set of documents that
    say nothing lints clean: use :func:`assess` for that.

    Args:
        documents: The files - a sequence of paths, or a ``{path: label}`` mapping, exactly as for
            :func:`pool`. **Label them**: labels are how ``documents`` names each file in the
            result, and mzLib's default depends on where the files sit.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`SdrfDrift`, one row per finding x variant.

    Raises:
        UsageError: no documents, a missing file, or a blank label.

    Example:
        >>> drift = lint({"a.sdrf.tsv": "cohort", "b.sdrf.tsv": "partner"})   # doctest: +SKIP
        >>> for variants in drift.findings():                                  # doctest: +SKIP
        ...     print(variants[0]["kind"], [v["value"] for v in variants])
        ValueCaseVariant ['Homo sapiens', 'homo sapiens']
    """
    lines = _document_lines(documents, verb="lint")
    data = _bridge.invoke("sdrf", "lint", stdin="\n".join(lines), timeout=timeout)
    return SdrfDrift._from_wire(data)


# ---- assess --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SdrfAssessment(_Table):
    """Whether one SDRF document's sample half describes an experimental design.

    mzLib's ``SdrfSampleInformativeness`` asks three questions, and the verdict is how many pass:
    all three is ``Informative``, none is ``Skeleton``, anything between is ``Partial``.

    ``columns`` is the evidence: one row per column a check read, with ``role``
    (``factor_value``, ``sample_characteristic`` or ``biological_replicate``), ``column_name``,
    ``rows``, ``filled`` (rows with a real answer), ``absent`` (empty or a reserved word),
    ``distinct_values`` (compared ignoring case and surrounding space) and ``fill_rate``
    (``filled / rows``, a fraction from 0 to 1).

    Attributes:
        path: The absolute path assessed.
        verdict: ``"Informative"``, ``"Partial"`` or ``"Skeleton"`` - see :data:`VERDICTS`.
        factor_value_varies: Some ``factor value[...]`` column holds two or more real answers.
        sample_is_described: Some ``characteristics[...]`` column other than organism and
            biological replicate holds a real answer. Organism is left out because a search can
            fill it without a human.
        biological_replicate_varies: ``characteristics[biological replicate]`` holds two or more
            real answers. ``False`` also when the column is missing (there is then no
            ``biological_replicate`` row in ``columns``).
        row_count: Data rows in the document.
        column_names: The evidence table's columns.
        columns: Column name to one value per evidence row.
        caveats: Including when ``Partial`` is perfectly legitimate.
    """

    path: str
    verdict: str
    factor_value_varies: bool
    sample_is_described: bool
    biological_replicate_varies: bool
    row_count: int
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfAssessment":
        names, columns = _columns(data)
        return cls(
            path=data.get("path", ""),
            verdict=data.get("verdict", ""),
            factor_value_varies=bool(data.get("factor_value_varies", False)),
            sample_is_described=bool(data.get("sample_is_described", False)),
            biological_replicate_varies=bool(data.get("biological_replicate_varies", False)),
            row_count=int(data.get("row_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class AssessedFile:
    """One document's verdict in an :func:`assess_many` result; ``None`` fields mean unread.

    Attributes:
        path: The absolute path.
        verdict: ``"Informative"``, ``"Partial"`` or ``"Skeleton"``; ``None`` if unread.
        factor_value_varies: As :attr:`SdrfAssessment.factor_value_varies`; ``None`` if unread.
        sample_is_described: As :attr:`SdrfAssessment.sample_is_described`; ``None`` if unread.
        biological_replicate_varies: As :attr:`SdrfAssessment.biological_replicate_varies`;
            ``None`` if unread.
        row_count: Data rows; ``None`` if unread.
        error: Why it was not read, or ``None``.
    """

    path: str
    verdict: str | None
    factor_value_varies: bool | None
    sample_is_described: bool | None
    biological_replicate_varies: bool | None
    row_count: int | None
    error: FileError | None

    @classmethod
    def _from_wire(cls, f: Mapping[str, Any]) -> "AssessedFile":
        return cls(
            path=f.get("path", ""),
            verdict=f.get("verdict"),
            factor_value_varies=f.get("factor_value_varies"),
            sample_is_described=f.get("sample_is_described"),
            biological_replicate_varies=f.get("biological_replicate_varies"),
            row_count=f.get("row_count"),
            error=FileError._from_wire(f.get("error")),
        )


@dataclass(frozen=True)
class SdrfAssessmentBatch(_Table):
    """Verdicts for many SDRF documents, with the evidence as one long table.

    Attributes:
        file_count: Documents given.
        read_count: Documents assessed.
        failed_count: Documents that could not be read (only under ``on_error="skip"``).
        verdict_counts: ``{"informative": n, "partial": n, "skeleton": n}`` over the documents read.
        files: One :class:`AssessedFile` per input, in input order.
        column_names: The evidence columns of :class:`SdrfAssessment`, with ``source_index`` and
            ``source_path`` first.
        columns: Column name to one value per evidence row.
        caveats: As for :class:`SdrfAssessment`.
    """

    file_count: int
    read_count: int
    failed_count: int
    verdict_counts: dict[str, int]
    files: list[AssessedFile]
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    def paths_with(self, *verdicts: str) -> list[str]:
        """The paths whose verdict is one of ``verdicts``, in input order.

        Raises:
            UsageError: a verdict that is not one of :data:`VERDICTS` - a typo would otherwise
                silently select nothing.

        Example:
            >>> keep = batch.paths_with("Informative", "Partial")   # doctest: +SKIP
        """
        unknown = [v for v in verdicts if v not in VERDICTS]
        if unknown:
            raise _bridge.UsageError(f"Unknown verdict(s) {unknown}; expected some of {VERDICTS}.")
        return [f.path for f in self.files if f.verdict in verdicts]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfAssessmentBatch":
        names, columns = _columns(data)
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            verdict_counts={k: int(v) for k, v in (data.get("verdict_counts") or {}).items()},
            files=[AssessedFile._from_wire(f) for f in (data.get("files") or [])],
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


def assess(path: str | os.PathLike[str], *, timeout: float | None = 60) -> SdrfAssessment:
    """Decide whether an SDRF file actually describes its samples, or is a valid skeleton.

    A file generated from a list of data files - reserved words in every sample column, one
    replicate number, no factor - passes :func:`validate`, because reserved words are the
    specification's correct way to say nothing. For grouping results by biology it is the same as
    having no SDRF at all. mzLib's ``SdrfSampleInformativeness.Assess`` is the gate for that.

    Args:
        path: Path to a ``.sdrf.tsv`` file.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`SdrfAssessment`: the verdict, the three checks, and the per-column counts.

    Raises:
        UsageError: the path is blank or the file does not exist.
        BridgeError: mzLib could not read the file.

    Example:
        >>> a = assess("PXD000070.sdrf.tsv")                                   # doctest: +SKIP
        >>> a.verdict, a.factor_value_varies, a.sample_is_described            # doctest: +SKIP
        ('Partial', False, True)
    """
    data = _bridge.invoke("sdrf", "assess", "--path", _one_path(path), timeout=timeout)
    return SdrfAssessment._from_wire(data)


def assess_many(
    paths: Sequence[str | os.PathLike[str]],
    *,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> SdrfAssessmentBatch:
    """Assess many SDRF files in one bridge call - the gate for a corpus.

    Args:
        paths: The files, in the order you want them reported. Each may be given once.
        threads: Documents assessed at once. ``1`` by default; ``-1`` uses every core. The result
            is the same at any value.
        on_error: ``"fail"`` (default) or ``"skip"``, as for :func:`validate_many`.
        timeout: Seconds for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        An :class:`SdrfAssessmentBatch`.

    Raises:
        UsageError: as for :func:`validate_many`.
        BridgeError: under ``"fail"``, mzLib could not read one of the files.

    Example:
        >>> batch = assess_many(paths, threads=-1, on_error="skip")   # doctest: +SKIP
        >>> batch.verdict_counts                                      # doctest: +SKIP
        {'informative': 1, 'partial': 1, 'skeleton': 1}
        >>> usable = batch.paths_with("Informative")                  # doctest: +SKIP
    """
    args, stdin = _bulk_args(paths, threads, on_error, verb="assess")
    data = _bridge.invoke("sdrf", "assess", *args, stdin=stdin, timeout=timeout)
    return SdrfAssessmentBatch._from_wire(data)


# ---- samples -------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SdrfSamples(_Table):
    """Every sample in one SDRF document and what its rows agree on, as a long table.

    A sample is a ``source name``, matched ignoring case and surrounding space, exactly as mzLib's
    ``SdrfSampleBlock.BySourceName`` keys it: one sample measured in three fractions is one sample.
    ``columns`` has one row per sample x column x position, samples in document order:

    - ``source_name`` - the sample, as its first row spelled it.
    - ``sample_row_count`` - how many rows of the document carry it.
    - ``column_name`` - the header's own spelling (``characteristics[Age]`` stays itself).
    - ``column_kind`` - ``source_name``, ``characteristic`` or ``factor_value``.
    - ``position`` - 0 for the first occurrence of a repeated column, 1 for the second, and so on;
      ``None`` on a conflicting row.
    - ``status`` - ``"agreed"``, or ``"conflicting"`` when the sample's rows disagree about the
      column: mzLib then **withholds** it rather than guess, and ``value`` is ``None``.
    - ``value`` - the cell, verbatim. Reserved words are real answers and come back as themselves.
    - ``age_years``, ``age_min_years``, ``age_max_years``, ``age_precision``,
      ``age_follows_specification``, ``age_refusal`` - filled only on ``characteristics[age]``
      rows, with exactly the meanings :class:`ParsedAges` documents. On any other row all six are
      ``None``.

    Attributes:
        path: The absolute path read.
        sample_count: Distinct samples (source names).
        row_count: Data rows in the document.
        conflict_count: Sample x column pairs withheld as conflicting.
        problems: Rows mzLib could not place - no ``source name`` column, or a blank source name.
            Empty for a well-formed document.
        column_names: The table's columns.
        columns: Column name to one value per row of the table.
        caveats: Including how to key samples across documents.
    """

    path: str
    sample_count: int
    row_count: int
    conflict_count: int
    problems: list[str]
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    def ages(self) -> list[dict[str, Any]]:
        """Only the parsed-age rows: one per sample whose ``characteristics[age]`` was agreed."""
        return [
            r for r in self.records if r["age_refusal"] is not None or r["age_years"] is not None
        ]

    def conflicts(self) -> list[tuple[str, str]]:
        """``(source_name, column_name)`` for every column mzLib withheld as conflicting."""
        return [
            (r["source_name"], r["column_name"])
            for r in self.records
            if r["status"] == "conflicting"
        ]

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfSamples":
        names, columns = _columns(data)
        return cls(
            path=data.get("path", ""),
            sample_count=int(data.get("sample_count", 0)),
            row_count=int(data.get("row_count", 0)),
            conflict_count=int(data.get("conflict_count", 0)),
            problems=list(data.get("problems") or []),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


@dataclass(frozen=True)
class SampledFile:
    """One document's summary in a :func:`samples_many` result; ``None`` fields mean unread.

    Attributes:
        path: The absolute path.
        sample_count: Distinct samples; ``None`` if unread.
        row_count: Data rows; ``None`` if unread.
        conflict_count: Columns withheld as conflicting; ``None`` if unread.
        problems: Rows mzLib could not place; ``None`` if unread.
        error: Why it was not read, or ``None``.
    """

    path: str
    sample_count: int | None
    row_count: int | None
    conflict_count: int | None
    problems: list[str] | None
    error: FileError | None

    @classmethod
    def _from_wire(cls, f: Mapping[str, Any]) -> "SampledFile":
        problems = f.get("problems")
        return cls(
            path=f.get("path", ""),
            sample_count=f.get("sample_count"),
            row_count=f.get("row_count"),
            conflict_count=f.get("conflict_count"),
            problems=list(problems) if problems is not None else None,
            error=FileError._from_wire(f.get("error")),
        )


@dataclass(frozen=True)
class SdrfSamplesBatch(_Table):
    """The samples of many SDRF documents, as one long table.

    ``columns`` has the columns of :class:`SdrfSamples`, preceded by ``source_index`` and
    ``source_path``. **A source name is only unique within one document** - two studies may both
    have a ``"Sample 1"`` - so key a sample on ``(source_index, source_name)``.

    Attributes:
        file_count: Documents given.
        read_count: Documents read.
        failed_count: Documents that could not be read (only under ``on_error="skip"``).
        sample_count: Samples over every document read.
        files: One :class:`SampledFile` per input, in input order.
        column_names: The table's columns.
        columns: Column name to one value per row.
        caveats: As for :class:`SdrfSamples`.
    """

    file_count: int
    read_count: int
    failed_count: int
    sample_count: int
    files: list[SampledFile]
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "SdrfSamplesBatch":
        names, columns = _columns(data)
        return cls(
            file_count=int(data.get("file_count", 0)),
            read_count=int(data.get("read_count", 0)),
            failed_count=int(data.get("failed_count", 0)),
            sample_count=int(data.get("sample_count", 0)),
            files=[SampledFile._from_wire(f) for f in (data.get("files") or [])],
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


def samples(path: str | os.PathLike[str], *, timeout: float | None = 60) -> SdrfSamples:
    """Lift each sample's characteristics and factor values out of an SDRF file, ages parsed.

    Calls mzLib's ``SdrfSampleBlock.BySourceName``: the sample half of the document (``source
    name``, every ``characteristics[...]``, every ``factor value[...]``) per sample, merged over
    the sample's rows. Where those rows **disagree** - one fraction says ``normal`` and another
    ``COVID-19`` - mzLib names the column and withholds it rather than let the first row win.

    Args:
        path: Path to a ``.sdrf.tsv`` file.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`SdrfSamples` long table.

    Raises:
        UsageError: the path is blank or the file does not exist.
        BridgeError: mzLib could not read the file.

    Example:
        >>> s = samples("cohort.sdrf.tsv")                                     # doctest: +SKIP
        >>> s.sample_count, s.conflicts()                                      # doctest: +SKIP
        (6, [('S6', 'characteristics[disease]')])
        >>> [(a["source_name"], a["age_years"]) for a in s.ages()][:3]         # doctest: +SKIP
        [('S1', 58), ('S2', 62.5), ('S3', 90)]
    """
    data = _bridge.invoke("sdrf", "samples", "--path", _one_path(path), timeout=timeout)
    return SdrfSamples._from_wire(data)


def samples_many(
    paths: Sequence[str | os.PathLike[str]],
    *,
    threads: int = 1,
    on_error: str = "fail",
    timeout: float | None = None,
) -> SdrfSamplesBatch:
    """Lift the samples of many SDRF files in one bridge call.

    Args:
        paths: The files, in the order you want them reported. Each may be given once.
        threads: Documents read at once. ``1`` by default; ``-1`` uses every core. The result is
            the same at any value.
        on_error: ``"fail"`` (default) or ``"skip"``, as for :func:`validate_many`.
        timeout: Seconds for the whole batch. ``None`` (the default) waits indefinitely.

    Returns:
        An :class:`SdrfSamplesBatch`.

    Raises:
        UsageError: as for :func:`validate_many`.
        BridgeError: under ``"fail"``, mzLib could not read one of the files.
    """
    args, stdin = _bulk_args(paths, threads, on_error, verb="samples")
    data = _bridge.invoke("sdrf", "samples", *args, stdin=stdin, timeout=timeout)
    return SdrfSamplesBatch._from_wire(data)


# ---- parse_ages ----------------------------------------------------------------------------------


@dataclass(frozen=True)
class ParsedAges(_Table):
    """Age cells read into years by mzLib's ``SdrfAge.TryParse``, one row per cell given.

    ``columns`` has one row per input cell, **in input order** (blank cells included), with:

    - ``cell`` - the cell exactly as you gave it.
    - ``years`` - the single figure to place the sample on an age axis, in **years**: the age,
      a range's **midpoint**, or a bound. ``None`` when refused.
    - ``min_years`` - the youngest the cell allows, in years; ``0`` for an upper bound.
      ``None`` when refused.
    - ``max_years`` - the oldest the cell allows, in years. ``None`` when refused **or** when there
      is no upper bound (``precision == "LowerBound"``, e.g. ``>=90Y``): mzLib's value is +infinity,
      which JSON cannot carry. ``precision`` tells the two apart.
    - ``precision`` - one of :data:`AGE_PRECISIONS`. ``None`` when refused.
    - ``follows_specification`` - ``True`` for the specification's own ``nYnMnD`` / ``nW`` grammar,
      ``False`` for unambiguous words (``"3 year"``, ``"6-8 weeks"``). ``None`` when refused.
    - ``refusal`` - ``None`` when the cell was read, else one of :data:`AGE_REFUSALS`.

    A month is 1/12 year, a week 7/365.25, a day 1/365.25, an hour 1/(365.25 x 24).

    Attributes:
        cell_count: Cells given; the length of every column.
        parsed_count: Cells read as an age.
        column_names: The table's columns.
        columns: Column name to one value per cell.
        caveats: The unit conventions and refusal rules.
    """

    cell_count: int
    parsed_count: int
    column_names: list[str]
    columns: dict[str, list[Any]]
    caveats: list[str] = field(default_factory=list)

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "ParsedAges":
        names, columns = _columns(data)
        return cls(
            cell_count=int(data.get("cell_count", 0)),
            parsed_count=int(data.get("parsed_count", 0)),
            column_names=names,
            columns=columns,
            caveats=list(data.get("caveats") or []),
        )


def parse_ages(cells: Sequence[str | None], *, timeout: float | None = 60) -> ParsedAges:
    """Read ``characteristics[age]`` cells into years, refusing anything that would need a guess.

    Calls mzLib's ``SdrfAge.TryParse`` on each cell. It reads the specification's grammar
    (``58Y``, ``30Y6M``, ``16W``), ranges (``40Y-85Y``, ``6-8 weeks``), bounds (``>=90Y``,
    ``<1Y``) and unambiguous words (``3 year``, ``4 hour``). It **refuses** a bare number - ``63``
    is 11% of real age cells, and 63 years and 63 days are both plausible in one study - as well as
    reserved words and free text. A refused cell has no age; ``refusal`` says why.

    One bridge call for any number of cells, so pass them all at once rather than looping.

    Args:
        cells: The cells, e.g. ``doc.value("characteristics[age]")``. ``None`` is sent as an empty
            cell and comes back refused as ``"empty"``, so the result stays aligned with the input.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        A :class:`ParsedAges` table, row ``i`` for cell ``i``.

    Raises:
        UsageError: ``cells`` is a bare string or empty, or a cell is not a string or contains a
            newline.

    Example:
        >>> ages = parse_ages(["58Y", "40Y-85Y", ">=90Y", "63"])            # doctest: +SKIP
        >>> ages.columns["years"], ages.columns["refusal"]                 # doctest: +SKIP
        ([58, 62.5, 90, None], [None, None, None, 'no_unit'])
    """
    if isinstance(cells, str):
        raise _bridge.UsageError(
            f"parse_ages() takes a list of cells; for one cell use parse_ages([{cells!r}])."
        )
    if not isinstance(cells, Sequence) or not cells:
        raise _bridge.UsageError("parse_ages() needs at least one cell, as a list.")
    lines = []
    for c in cells:
        text = "" if c is None else c
        if not isinstance(text, str):
            raise _bridge.UsageError(f"Every cell must be a str or None; got {c!r}.")
        if "\n" in text or "\r" in text:
            raise _bridge.UsageError(f"A cell contains a newline, which separates cells: {text!r}.")
        lines.append(text)
    # A trailing newline ends the last cell, so a final blank cell survives the trip.
    data = _bridge.invoke("sdrf", "parse-age", stdin="\n".join(lines) + "\n", timeout=timeout)
    return ParsedAges._from_wire(data)
