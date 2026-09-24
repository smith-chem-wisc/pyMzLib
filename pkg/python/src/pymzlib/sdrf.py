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

**What this module does not do yet is validate.** mzLib models SDRF's structural rules in
``SdrfValidator`` and its vocabulary-drift rules in ``SdrfDriftLint``. Both have been public
since mzLib #1207, which the pinned mzLib (1.0.589) includes, so the bridge can call them. They
are not exposed yet. When they are, the rules will be projected once, in the bridge, for all
three bindings - a second implementation here is exactly the per-binding repair pyMzLib exists
to avoid. Until then, this module reads, pools and reports honestly, and makes no claim about
whether a document is *correct*.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from . import _bridge

__all__ = ["SdrfDocument", "PooledSdrf", "WrittenSdrf", "read", "pool"]


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
    if isinstance(documents, Mapping):
        pairs = [(str(p), str(label)) for p, label in documents.items()]
    elif isinstance(documents, str):
        # A bare string is a sequence of characters, so this would otherwise pool one document per
        # letter and fail with a baffling "file not found: P".
        raise _bridge.UsageError(
            "pool() takes several documents; pass a list or a {path: label} mapping. "
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

    args = ["sdrf", "pool", *_window(limit, offset)]
    if out is not None:
        out = _bridge.path_text(out)
        if not out:
            raise _bridge.UsageError("out must be a non-empty path or None.")
        args += ["--out", out]

    data = _bridge.invoke(*args, stdin="\n".join(lines), timeout=timeout)
    return PooledSdrf._from_wire(data)
