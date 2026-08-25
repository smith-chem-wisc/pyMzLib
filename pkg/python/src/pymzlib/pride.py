"""PRIDE Archive access, backed by mzLib's ``PrideArchiveClient``.

The PRIDE Archive (https://www.ebi.ac.uk/pride/archive/) is EBI's public proteomics data
repository. This module lets a Python user list what is in a project and pull files down,
using the same paging, URL-resolution, and safe-download logic that mzLib uses in C#.

    >>> import pymzlib
    >>> files = pymzlib.pride.list_files("PXD000001")
    >>> len(files)
    8
    >>> files[0].file_name
    'PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz'
    >>> raw = [f for f in files if f.category == "RAW"]
    >>> pymzlib.pride.download("PXD000001", "downloads", category="RAW")
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import _bridge

__all__ = [
    "PrideFile",
    "PrideFtpFile",
    "PrideProjectSearchResult",
    "ProjectNotFoundError",
    "list_files",
    "search",
    "list_ftp_files",
    "download",
    "download_files",
    "total_size_bytes",
    "approximate_total_size_bytes",
]

#: A PRIDE-style repository accession: a short letter prefix and a run of digits, e.g. PXD000001.
_ACCESSION_PATTERN = re.compile(r"^[A-Z]{2,4}[0-9]{4,}$")


class ProjectNotFoundError(_bridge.PyMzLibError):
    """No project with that accession exists, or it has no files.

    PRIDE answers an unknown accession with an empty result rather than a 404, so earlier versions
    of pyMzLib returned an empty list. That was a mistake: an empty list is indistinguishable from
    "this project genuinely has nothing matching", so a typo'd accession produced a script that
    reported "0 files, done" and moved on. A wrong answer that looks like a right answer is worse
    than an error.
    """


def _normalise_accession(accession: object) -> str:
    """Validate and canonicalise an accession, failing loudly rather than returning nothing.

    Accessions are upper-cased, because PRIDE's API is case-sensitive on the accession while
    pyMzLib's own category matching is case-insensitive — two rules pointing opposite ways is a
    trap, and this is the one that can be fixed without surprising anybody.
    """
    if not isinstance(accession, str):
        raise _bridge.UsageError(
            f"accession must be a string like 'PXD000001'; got {type(accession).__name__} ({accession!r})."
        )

    candidate = accession.strip().upper()
    if not candidate:
        raise _bridge.UsageError("A PRIDE project accession is required, e.g. 'PXD000001'.")
    if not _ACCESSION_PATTERN.match(candidate):
        raise _bridge.UsageError(
            f"'{accession}' is not a valid repository accession. Expected a short letter prefix "
            "followed by digits, e.g. 'PXD000001'."
        )
    return candidate


def _normalise_destination(destination: object) -> Path:
    """Reject a blank destination instead of quietly writing into the current directory.

    ``Path("")`` is ``Path(".")``, so ``dest = config.get("outdir", "")`` used to spray a
    multi-gigabyte project across the working directory. The docstring had always promised this
    raised; now it does.
    """
    if not isinstance(destination, (str, Path)):
        raise _bridge.UsageError(
            f"destination must be a path or string; got {type(destination).__name__} ({destination!r})."
        )

    # Check the ORIGINAL text, not the constructed Path. `Path("")` is `Path(".")`, whose str() is
    # "." — non-blank and truthy — so a guard applied after construction lets `Path("")` through
    # and writes into the working directory, which is the very thing this function exists to stop.
    raw = str(destination)
    if not raw.strip():
        raise _bridge.UsageError("A destination directory is required; got an empty path.")

    return Path(destination)


def _normalise_extensions(extensions: object) -> list[str]:
    """Accept a sequence of extensions, and refuse a bare string.

    ``extensions=".raw"`` satisfies the ``Sequence[str]`` annotation and iterates as four
    characters, so it used to match nothing, download nothing, and exit successfully. In a batch
    script that is zero files and a green exit code.
    """
    if extensions is None:
        return []
    if isinstance(extensions, str):
        raise _bridge.UsageError(
            f"extensions must be a list of extensions, not a single string. "
            f"Did you mean [{extensions!r}]?"
        )
    try:
        values = list(extensions)
    except TypeError as exc:
        raise _bridge.UsageError(
            f"extensions must be a list of extensions; got {type(extensions).__name__}."
        ) from exc

    for value in values:
        if not isinstance(value, str):
            raise _bridge.UsageError(
                f"Each extension must be a string; got {type(value).__name__} ({value!r})."
            )
        if "," in value:
            raise _bridge.UsageError(
                f"An extension may not contain a comma; got {value!r}. Pass separate list items."
            )

    kept = [v.strip() for v in values if v.strip()]
    # Fail-open again, one layer up: a caller who asked for extensions and whose list normalises to
    # nothing would have had `--ext` omitted entirely, which the bridge reads as "no filter" and
    # downloads the whole project. Asking for a filter and getting everything is never right.
    if values and not kept:
        raise _bridge.UsageError(
            f"extensions was given but names no extensions; got {list(values)!r}. "
            "Omit it to download every file type."
        )
    return kept


def _reject_flag_like(name: str, value: str) -> str:
    """Refuse a value that would be parsed as another option by the bridge.

    The bridge's parser treats ``--a --b`` as two flags, so a value beginning with ``--`` silently
    discards the option it belonged to — and can smuggle in a flag the caller never intended.
    """
    if value.startswith("-"):
        raise _bridge.UsageError(
            f"{name} may not begin with '-'; got {value!r}. That would be read as another option."
        )
    return value


def _parse_timestamp(value: str | None) -> datetime | None:
    """Convert an ISO-8601 timestamp from the bridge into a timezone-aware ``datetime``."""
    if not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


@dataclass(frozen=True)
class PrideFile:
    """One file belonging to a PRIDE Archive project.

    Attributes:
        file_name: The file's name, e.g. ``"run1.raw"``.
        file_size_bytes: Size in bytes as reported by PRIDE.
        checksum: The repository's checksum, or ``""`` if it provides none.
        category: The file category, e.g. ``"RAW"``, ``"PEAK"``, ``"SEARCH"``, ``"OTHER"``.
        https_url: A direct HTTPS download URL, or ``None`` when the file is only reachable
            by a protocol that cannot be fetched over HTTPS (Aspera-only files).
        locations: Every published location as ``{"accession", "name", "value"}`` dicts, for
            callers that want the raw controlled-vocabulary terms.
        submission_date / publication_date / updated_date: Repository timestamps.
    """

    file_name: str
    file_size_bytes: int
    checksum: str
    category: str
    category_accession: str
    https_url: str | None
    locations: list[dict[str, str]] = field(default_factory=list)
    submission_date: datetime | None = None
    publication_date: datetime | None = None
    updated_date: datetime | None = None
    project_accession: str = ""

    @property
    def size_mb(self) -> float:
        """The file size in megabytes, for the common case of eyeballing a manifest."""
        return self.file_size_bytes / 1_000_000

    @property
    def extension(self) -> str:
        """The file's lowercase extension including the dot, e.g. ``".raw"``. Empty if none."""
        return Path(self.file_name).suffix.lower()

    @property
    def downloadable(self) -> bool:
        """Whether this file can be fetched by :func:`download` (i.e. has an HTTPS location)."""
        return self.https_url is not None

    def as_dict(self) -> dict[str, Any]:
        """Return every attribute, **including the computed ones**, as a plain dict.

        Use this rather than ``vars(f)`` when building a table. ``size_mb``, ``extension`` and
        ``downloadable`` are properties, so ``vars()`` and ``dataclasses.asdict()`` both skip
        them — which silently produced a DataFrame missing the three attributes the
        documentation pushes hardest, including the ``downloadable`` flag used to filter out
        files that cannot be fetched.

        Example:
            >>> import pandas as pd                                    # doctest: +SKIP
            >>> df = pd.DataFrame([f.as_dict() for f in files])        # doctest: +SKIP
        """
        record = asdict(self)
        record["size_mb"] = self.size_mb
        record["extension"] = self.extension
        record["downloadable"] = self.downloadable
        return record

    @classmethod
    def _from_wire(cls, payload: dict[str, Any], project_accession: str = "") -> "PrideFile":
        return cls(
            file_name=payload.get("file_name", ""),
            file_size_bytes=int(payload.get("file_size_bytes", 0)),
            checksum=payload.get("checksum", ""),
            category=payload.get("category", ""),
            category_accession=payload.get("category_accession", ""),
            https_url=payload.get("https_url"),
            locations=list(payload.get("locations") or []),
            submission_date=_parse_timestamp(payload.get("submission_date")),
            publication_date=_parse_timestamp(payload.get("publication_date")),
            updated_date=_parse_timestamp(payload.get("updated_date")),
            project_accession=project_accession,
        )


@dataclass(frozen=True)
class PrideFtpFile:
    """One file found by walking a PRIDE project's **FTP directory tree** — the complete listing.

    This is what :func:`list_ftp_files` returns, and the difference from :class:`PrideFile` is the
    whole point: the FTP walk sees everything the project holds, including the files PRIDE's REST
    manifest omits and files nested in subdirectories. The trade-off is the size: PRIDE's directory
    index rounds it to about three significant figures, so it is ``approximate_size_bytes`` — a good
    project-size estimate, not the exact number of bytes you will transfer.

    Attributes:
        relative_path: Path relative to the project's FTP root, e.g. ``"run1.raw"`` or, for a file
            in a subdirectory, ``"generated/summary.mztab"``.
        file_name: The bare file name — the last segment of ``relative_path``.
        url: The HTTPS URL the file can be downloaded from.
        approximate_size_bytes: PRIDE's rounded index size in bytes. For the exact transfer size of
            one file, issue an HTTP HEAD against ``url`` and read its ``Content-Length``.
        project_accession: The accession this file was listed under.
    """

    relative_path: str
    file_name: str
    url: str
    approximate_size_bytes: int
    project_accession: str = ""

    @property
    def approximate_size_mb(self) -> float:
        """The approximate size in megabytes, for eyeballing a project's footprint."""
        return self.approximate_size_bytes / 1_000_000

    @property
    def extension(self) -> str:
        """The file's lowercase extension including the dot, e.g. ``".raw"``. Empty if none."""
        return Path(self.file_name).suffix.lower()

    def as_dict(self) -> dict[str, Any]:
        """Return every attribute, **including the computed ones**, as a plain dict for a DataFrame."""
        record = asdict(self)
        record["approximate_size_mb"] = self.approximate_size_mb
        record["extension"] = self.extension
        return record

    @classmethod
    def _from_wire(cls, payload: dict[str, Any], project_accession: str = "") -> "PrideFtpFile":
        return cls(
            relative_path=payload.get("relative_path", ""),
            file_name=payload.get("file_name", ""),
            url=payload.get("url", ""),
            approximate_size_bytes=int(payload.get("approximate_size_bytes", 0)),
            project_accession=project_accession,
        )


def list_files(accession: str, page_size: int = 100, timeout: float | None = 300) -> list[PrideFile]:
    """Return the file manifest of a PRIDE Archive project.

    **This is what PRIDE's REST API publishes, which is not always everything in the
    project.** For PXD000001 the API returns **8** files while the FTP tree holds **13**, and
    the five it omits include the two largest: ``...60min_01-20141210.mzML`` (450 MB) and the
    matching ``.mzXML`` (472 MB), exactly the modern open-format conversions most people want.
    The omission is PRIDE's, not mzLib's. **If completeness matters, use :func:`list_ftp_files`,**
    which walks the project's FTP directory and returns everything it actually holds.

    Paging is handled for you: however many pages the project spans, you get one list.

    Args:
        accession: The project accession, e.g. ``"PXD000001"``.
        page_size: How many files to request per underlying API call. Only affects how the
            manifest is fetched, never what you get back.
        timeout: Seconds to allow for the whole fetch.

    Returns:
        Every file in the project, in repository order. An unknown accession yields an empty
        list rather than an error — that is PRIDE's own behavior, preserved here.

    Raises:
        UsageError: the accession is blank or the page size is not positive.
        BridgeError: PRIDE returned an error status or was unreachable.
    """
    canonical = _normalise_accession(accession)
    if isinstance(page_size, bool) or not isinstance(page_size, int):
        raise _bridge.UsageError(
            f"page_size must be a whole number; got {type(page_size).__name__} ({page_size!r})."
        )
    if page_size <= 0:
        raise _bridge.UsageError(f"page_size must be positive; got {page_size}.")
    if page_size > 2_147_483_647:
        raise _bridge.UsageError(f"page_size is larger than the API allows; got {page_size}.")

    data = _bridge.invoke(
        "pride", "files",
        "--accession", canonical,
        "--page-size", str(page_size),
        timeout=timeout,
    )
    files = [PrideFile._from_wire(item, canonical) for item in data.get("files", [])]

    if not files:
        raise ProjectNotFoundError(
            f"PRIDE returned no files for '{canonical}'. Either the accession does not exist "
            "(check for a typo) or the project is private. PRIDE does not distinguish the two, "
            "so neither can pyMzLib."
        )
    return files


def list_ftp_files(accession: str, timeout: float | None = 300) -> list[PrideFtpFile]:
    """Return the **complete** file list of a PRIDE project, read from its FTP directory tree.

    This is the authoritative counterpart to :func:`list_files`. Where ``list_files`` returns
    PRIDE's REST manifest — which is knowingly incomplete, omitting for PXD000001 five of the
    project's 13 files, including the two largest — this walks the FTP directory (subdirectories
    included) and returns everything the project actually holds. Reach for it whenever completeness
    or a true project size matters; use :func:`list_files` when you want the rich metadata (category,
    checksum, controlled-vocabulary locations) that the REST manifest carries and the directory
    index does not.

    This is a **listing** surface only. :func:`download` and :func:`download_files` operate on the
    REST manifest, so a file that appears *only* here — the whole point of this function — is not
    accepted by them; fetch it directly from its :attr:`PrideFtpFile.url` with an ordinary HTTPS
    client (e.g. ``urllib.request.urlretrieve(f.url, f.file_name)``).

    The sizes are approximate: PRIDE's directory index rounds them (see
    :attr:`PrideFtpFile.approximate_size_bytes`), so :func:`approximate_total_size_bytes` is an
    estimate — but an estimate over the *whole* project, unlike :func:`total_size_bytes`.

        >>> ftp = pymzlib.pride.list_ftp_files("PXD000001")           # doctest: +SKIP
        >>> len(ftp)                                                  # doctest: +SKIP
        13
        >>> nested = [f.relative_path for f in ftp if "/" in f.relative_path]   # doctest: +SKIP

    Args:
        accession: The project accession, e.g. ``"PXD000001"``.
        timeout: Seconds to allow for the whole walk, which spans one request per directory.

    Returns:
        Every file under the project's FTP root, subdirectories included, in the order the walk
        encounters them. Never empty — an empty result is raised as an error (see below).

    Raises:
        UsageError: the accession is blank or malformed.
        ProjectNotFoundError: no project has that accession (or it lacks the publication date that
            locates its FTP directory), or the directory listed no files. Same "no such project"
            signal :func:`list_files` raises, so one ``except`` catches both — the difference from
            ``list_files`` is only *how* it is detected (mzLib resolves the project before walking,
            so a typo fails loudly rather than looking like an empty project).
        ServiceUnavailableError / BridgeError: PRIDE was unreachable, or a directory fetch failed.
    """
    canonical = _normalise_accession(accession)

    try:
        data = _bridge.invoke(
            "pride", "ftp-files",
            "--accession", canonical,
            timeout=timeout,
        )
    except _bridge.BridgeError as exc:
        # mzLib resolves the project (and its publication date) before walking, so an unknown
        # accession comes back as an MzLibException rather than an empty list. Re-map it to the same
        # ProjectNotFoundError that list_files() raises so callers catch one type for "not there".
        # ServiceUnavailableError (a BridgeError subclass) has a different error_type and is left to
        # propagate — an outage is not a missing project.
        if exc.error_type == "MzLibException":
            raise ProjectNotFoundError(
                f"PRIDE has no project '{canonical}' (or it lacks the publication date needed to "
                "locate its FTP directory). Check for a typo — a private project looks the same."
            ) from exc
        raise

    files = [PrideFtpFile._from_wire(item, canonical) for item in data.get("files", [])]
    if not files:
        # The project resolved but its FTP directory listed nothing. Almost always PRIDE changing
        # its autoindex format, never a real published project. Raise rather than hand back [], the
        # same "0 files, done" trap list_files() defends against.
        raise ProjectNotFoundError(
            f"The FTP directory for '{canonical}' listed no files. Either the project is genuinely "
            "empty or PRIDE's directory-index format has changed; list_files() may still return "
            "its REST manifest."
        )
    return files


def download(
    accession: str,
    destination: str | Path,
    category: str | None = None,
    extensions: Sequence[str] | None = None,
    overwrite: bool = True,
    timeout: float | None = None,
) -> list[Path]:
    """Download a project's files, optionally filtered, and return where they landed.

    Files are streamed to a temporary name and moved into place only once complete, so an
    interrupted download never leaves a truncated file behind.

    Args:
        accession: The project accession, e.g. ``"PXD000001"``.
        destination: Directory to write into. Created if it does not exist.
        category: Keep only files of this category, e.g. ``"RAW"``. ``None`` keeps all.
        extensions: Keep only files with these extensions, e.g. ``[".raw", ".mzML"]``.
            ``None`` keeps all. Combined with ``category`` as AND.
        overwrite: When ``False``, a file already present at the destination is left alone
            and not re-fetched — a cheap resume for a large project.
        timeout: Seconds to allow. ``None`` (the default) waits as long as it takes, which
            is usually what you want for multi-gigabyte projects.

    Returns:
        The paths written, in manifest order.

    Raises:
        UsageError: the accession or destination is blank.
        BridgeError: a request failed, or a selected file has no HTTPS location.
    """
    canonical = _normalise_accession(accession)
    target = _normalise_destination(destination)
    wanted = _normalise_extensions(extensions)

    args = ["pride", "download", "--accession", canonical, "--dest", str(target)]
    if category is not None:
        if not isinstance(category, str):
            raise _bridge.UsageError(
                f"category must be a string like 'RAW'; got {type(category).__name__} ({category!r})."
            )
        if not category.strip():
            raise _bridge.UsageError(
                "category is empty. Omit it to download every category, rather than passing a "
                "blank value — a filter that selects nothing must not silently select everything."
            )
        args += ["--category", _reject_flag_like("category", category.strip())]
    if wanted:
        args += ["--ext", _reject_flag_like("extensions", ",".join(wanted))]
    if not overwrite:
        args.append("--no-overwrite")

    data = _bridge.invoke(*args, timeout=timeout)
    written = [Path(p) for p in data.get("paths", [])]

    # Same doctrine as list_files and download_files, which this function was inconsistent with:
    # a filter that matched nothing is nearly always a filter that does not mean what its author
    # thought, and reporting success with an empty list lets a batch script carry on as though
    # the work had been done.
    if not written and (category is not None or wanted):
        raise _bridge.UsageError(
            f"No file in {canonical} matched "
            f"{'category ' + repr(category) if category else ''}"
            f"{' and ' if category and wanted else ''}"
            f"{'extensions ' + repr(wanted) if wanted else ''}. "
            "Use list_files() to see what the project actually contains — note that "
            "compressed files such as 'x.mgf.gz' have the extension '.gz'."
        )
    return written


def download_files(
    files: Iterable[PrideFile],
    destination: str | Path,
    overwrite: bool = True,
    timeout: float | None = None,
) -> list[Path]:
    """Download exactly the files you selected, and nothing else.

    This is the counterpart to :func:`list_files`, and usually the one you want. Filter the
    manifest however you like — in Python, with the full expressiveness of Python — and hand the
    result straight back:

        >>> files = list_files("PXD000001")                             # doctest: +SKIP
        >>> small = [f for f in files if f.size_mb < 5 and f.downloadable]
        >>> download_files(small, "downloads")                          # doctest: +SKIP

    :func:`download`'s ``category`` and ``extensions`` filters can only express what they were
    built to express; "under 5 MB", "the three newest", or "everything except the MGF" cannot be
    said in that vocabulary at all. They can all be said in a list comprehension.

    Args:
        files: The :class:`PrideFile` objects to fetch, from one project.
        destination: Directory to write into. Created if it does not exist.
        overwrite: When ``False``, files already present are left alone and not re-fetched.
        timeout: Seconds to allow. ``None`` waits as long as it takes.

    Returns:
        The paths written, in the order the repository lists them.

    Raises:
        UsageError: the selection is empty, spans several projects, or includes a file with no
            HTTPS location.
    """
    target = _normalise_destination(destination)
    selected = list(files)

    if not selected:
        raise _bridge.UsageError(
            "No files selected. An empty selection is almost always a filter that did not match "
            "what you expected, so pyMzLib refuses it rather than reporting success."
        )
    for item in selected:
        if not isinstance(item, PrideFile):
            raise _bridge.UsageError(
                f"download_files expects PrideFile objects from list_files(); got "
                f"{type(item).__name__} ({item!r})."
            )

    unreachable = [f.file_name for f in selected if not f.downloadable]
    if unreachable:
        raise _bridge.UsageError(
            f"{len(unreachable)} of {len(selected)} selected files have no HTTPS location and "
            f"cannot be downloaded (e.g. {unreachable[0]!r}). Filter on `.downloadable` first."
        )

    accessions = {f.project_accession for f in selected if f.project_accession}
    if len(accessions) > 1:
        raise _bridge.UsageError(
            f"All files must come from one project; got {sorted(accessions)}."
        )
    if not accessions:
        raise _bridge.UsageError(
            "These PrideFile objects carry no project accession, so pyMzLib cannot tell which "
            "project to fetch from. Obtain them from list_files()."
        )

    args = [
        "pride", "download",
        "--accession", accessions.pop(),
        "--dest", str(target),
        "--names-from-stdin",
    ]
    if not overwrite:
        args.append("--no-overwrite")

    # The selection travels on stdin rather than argv: a few thousand names would blow the ~32 KB
    # command-line ceiling. The framing is newline-delimited, which is *almost* general — a POSIX
    # file name may legally contain a newline, so such a name would split into two and silently
    # select the wrong files. PRIDE has never published one, but "never seen it" is not a contract,
    # so it is refused explicitly rather than mis-parsed quietly.
    embedded_newline = [f.file_name for f in selected if "\n" in f.file_name or "\r" in f.file_name]
    if embedded_newline:
        raise _bridge.UsageError(
            f"Cannot select {embedded_newline[0]!r}: the file name contains a line break, which the "
            "selection format cannot represent. Please open an issue — this is a limitation worth "
            "fixing properly if a real repository ever publishes such a name."
        )

    payload = "\n".join(f.file_name for f in selected)
    data = _bridge.invoke(*args, stdin=payload, timeout=timeout)
    return [Path(p) for p in data.get("paths", [])]


def total_size_bytes(files: Iterable[PrideFile]) -> int:
    """Sum the sizes of some files.

    **This is the size PRIDE reports, which is not the number of bytes you will transfer.**
    For compressed files PRIDE frequently reports the *decompressed* size: in PXD000001 the
    reported size of ``PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz`` is 16,448,103 bytes, exactly
    what ``gzip -l`` gives as its uncompressed length, while the actual download is 5,984,662
    bytes, 2.75x smaller.

    **It is also a sum over an incomplete manifest.** For PXD000001 this returns 0.51 GB; the
    project on disk is 1.44 GB, because PRIDE's API omits five files including the two largest
    (see :func:`list_files`). The two errors run in opposite directions and do **not** cancel. For
    a size that covers the whole project, use :func:`approximate_total_size_bytes` over
    :func:`list_ftp_files`.

    >>> files = list_files("PXD000001")           # doctest: +SKIP
    >>> total_size_bytes(f for f in files if f.category == "RAW") / 1e9   # doctest: +SKIP
    0.51
    """
    return sum(f.file_size_bytes for f in files)


def approximate_total_size_bytes(files: Iterable[PrideFtpFile]) -> int:
    """Sum the approximate sizes of some FTP files.

    This is the honest project-size number, and the counterpart to :func:`total_size_bytes` with the
    trade-offs reversed. It sums over the **complete** FTP listing (:func:`list_ftp_files`), so no
    files are missing — but each size is PRIDE's directory-index value, rounded to about three
    significant figures, so the total is an **estimate**, not an exact byte count. For PXD000001 it
    lands near the true 1.44 GB, where :func:`total_size_bytes` reports 0.51 GB over the incomplete
    REST manifest. When you need the exact bytes for one file, HTTP HEAD its
    :attr:`PrideFtpFile.url` and read ``Content-Length``.

    >>> ftp = list_ftp_files("PXD000001")                       # doctest: +SKIP
    >>> approximate_total_size_bytes(ftp) / 1e9                 # doctest: +SKIP
    1.44
    """
    return sum(f.approximate_size_bytes for f in files)


#: PRIDE answers a longer keyword with HTTP 500, which is indistinguishable from an outage, so
#: mzLib refuses one first and so do we. Mirrors ``PrideArchiveClient.MaxKeywordLength``.
MAX_KEYWORD_LENGTH = 1000


def _parse_date(value: str | None) -> date | None:
    """Convert a bare calendar date from the bridge into a :class:`datetime.date`.

    Deliberately **not** :func:`_parse_timestamp`, and the distinction is the whole point. The
    search endpoint sends a date with no time and no offset — ``"2026-08-15"`` — so mzLib types
    these ``DateTime`` where :class:`PrideFile`'s are ``DateTimeOffset``, precisely to avoid
    attaching an offset the wire never carried. ``datetime.fromisoformat`` on a bare date returns
    a naive **midnight**, which is a time PRIDE never reported; a :class:`date` says exactly what
    was sent and nothing more.
    """
    if not value:
        return None
    try:
        return date.fromisoformat(value[:10])
    except ValueError:
        return None


@dataclass(frozen=True)
class PrideProjectSearchResult:
    """One hit from :func:`search`.

    **This is not a project's full metadata, and the two are not interchangeable.** PRIDE serves
    search from a separate Elasticsearch projection in which every controlled-vocabulary field has
    been **flattened to a display string**: the same project reports its instruments as
    ``["Q Exactive"]`` here and as structured terms with accessions from the metadata endpoint,
    contacts collapse from ten-field objects to a display name, and publications to a single
    pre-formatted citation string. That is a property of PRIDE's wire, not a simplification chosen
    here — the accessions are simply not sent. Follow :attr:`accession` when you need the
    vocabulary.

    **A zero or an empty list means "not reported", never a measured zero.** PRIDE omits nothing as
    null, so absence arrives as ``0``, ``""`` or ``[]``, and several fields are genuinely sparse —
    sampled across 1,600 hits, :attr:`project_tags` was populated on 2.6%, :attr:`sdrf` on 2.4%,
    :attr:`other_omics_links` on 18%, and the bot/hub/organic trio on under half. Do not read a
    ``download_count`` of 0 as "nobody downloaded it".

    Attributes:
        accession: The project accession — the key to everything else in this module. Usually a
            ``PXD`` accession, but search also returns legacy ``PRD`` and affinity ``PAD`` ones, so
            treat it as an opaque key rather than assuming a prefix.
        title: The project title.
        project_description: The submitter's free-text description.
        sample_processing_protocol: How the sample was prepared, as free text.
        data_processing_protocol: How the data were searched and processed, as free text.
        doi: The dataset DOI, or ``""`` if PRIDE has not minted one.
        submission_type: ``"COMPLETE"``, ``"PARTIAL"``, ``"AFFINITY"``, or ``"PRIDE"`` for legacy
            submissions. Not a closed set — PRIDE has added values before.
        sdrf: The project's SDRF metadata as a single space-joined bag of term *values*, flattened
            by the search index. **Not a file, filename or URL** — nothing can be fetched with it
            and the row/column structure is gone. For a real SDRF see :mod:`pymzlib.sdrf`.
        submission_date / publication_date / updated_date: Calendar :class:`~datetime.date` values,
            not timestamps — see :func:`_parse_date`. ``None`` when PRIDE reported none.
        project_tags: PRIDE's coarse classification tags.
        keywords: The submitter's free-text keywords. **May contain empty and whitespace-only
            strings** — PRIDE ships them on roughly 9% of hits. They are passed through rather than
            filtered, because dropping them here would make this module disagree with mzLib and with
            the Rust and R bindings about what a project's keywords are. Filter before joining.
        submitters / lab_pis / affiliations: Display names and affiliations, flattened from the
            structured contact objects the metadata endpoint returns.
        instruments / softwares / quantification_methods: Display names.
        sample_attributes: Sample characteristics by display *value* (e.g. ``"liver"``). Flattened:
            which characteristic each value describes is **not recoverable** from a search hit.
        organisms / organism_parts / diseases: Display names.
        references: Publications, each a single pre-formatted citation string. A PubMed ID or DOI
            cannot be read out of one without parsing the string PRIDE assembled.
        experiment_types: e.g. ``"Data-independent acquisition"``.
        project_file_names: File *names* only — a search convenience, **not the manifest**. It
            carries no sizes, categories or download locations. Use :func:`list_files`, or
            :func:`list_ftp_files` for the complete list, to act on the files.
        other_omics_links: Links to related datasets in other omics repositories.
        highlights: Why this project matched, keyed by the field each match was found in, with the
            matched terms wrapped in ``<em>`` markup. The keys vary per hit and per query. This is
            the one thing search returns that the metadata endpoint cannot.
        yearly_downloads: ``{"year", "count"}`` dicts.
        download_count / avg_downloads_per_file / percentile: Download popularity.
        bot_count / hub_count / organic_count: Downloads split by traffic kind.
    """

    accession: str
    title: str = ""
    project_description: str = ""
    sample_processing_protocol: str = ""
    data_processing_protocol: str = ""
    doi: str = ""
    submission_type: str = ""
    sdrf: str = ""
    submission_date: date | None = None
    publication_date: date | None = None
    updated_date: date | None = None
    project_tags: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    submitters: list[str] = field(default_factory=list)
    lab_pis: list[str] = field(default_factory=list)
    affiliations: list[str] = field(default_factory=list)
    instruments: list[str] = field(default_factory=list)
    softwares: list[str] = field(default_factory=list)
    quantification_methods: list[str] = field(default_factory=list)
    sample_attributes: list[str] = field(default_factory=list)
    organisms: list[str] = field(default_factory=list)
    organism_parts: list[str] = field(default_factory=list)
    diseases: list[str] = field(default_factory=list)
    references: list[str] = field(default_factory=list)
    experiment_types: list[str] = field(default_factory=list)
    project_file_names: list[str] = field(default_factory=list)
    other_omics_links: list[str] = field(default_factory=list)
    highlights: dict[str, list[str]] = field(default_factory=dict)
    yearly_downloads: list[dict[str, Any]] = field(default_factory=list)
    download_count: int = 0
    avg_downloads_per_file: float = 0.0
    percentile: int = 0
    bot_count: int = 0
    hub_count: int = 0
    organic_count: int = 0

    @property
    def matched_fields(self) -> list[str]:
        """Which PRIDE fields the query hit, from :attr:`highlights`. Empty if PRIDE reported none."""
        return sorted(self.highlights)

    def as_dict(self) -> dict[str, Any]:
        """Every attribute, **including the computed ones**, as a plain dict.

        Use this rather than ``vars()`` when building a table: ``matched_fields`` is a property, so
        ``dataclasses.asdict()`` silently omits it. Same reasoning as :meth:`PrideFile.as_dict`.
        """
        record = asdict(self)
        record["matched_fields"] = self.matched_fields
        return record

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "PrideProjectSearchResult":
        def strings(name: str) -> list[str]:
            return [str(x) for x in (payload.get(name) or [])]

        return cls(
            accession=payload.get("accession", ""),
            title=payload.get("title", ""),
            project_description=payload.get("project_description", ""),
            sample_processing_protocol=payload.get("sample_processing_protocol", ""),
            data_processing_protocol=payload.get("data_processing_protocol", ""),
            doi=payload.get("doi", ""),
            submission_type=payload.get("submission_type", ""),
            sdrf=payload.get("sdrf", ""),
            submission_date=_parse_date(payload.get("submission_date")),
            publication_date=_parse_date(payload.get("publication_date")),
            updated_date=_parse_date(payload.get("updated_date")),
            project_tags=strings("project_tags"),
            keywords=strings("keywords"),
            submitters=strings("submitters"),
            lab_pis=strings("lab_pis"),
            affiliations=strings("affiliations"),
            instruments=strings("instruments"),
            softwares=strings("softwares"),
            quantification_methods=strings("quantification_methods"),
            sample_attributes=strings("sample_attributes"),
            organisms=strings("organisms"),
            organism_parts=strings("organism_parts"),
            diseases=strings("diseases"),
            references=strings("references"),
            experiment_types=strings("experiment_types"),
            project_file_names=strings("project_file_names"),
            other_omics_links=strings("other_omics_links"),
            highlights={
                str(k): [str(v) for v in (vs or [])]
                for k, vs in (payload.get("highlights") or {}).items()
            },
            yearly_downloads=list(payload.get("yearly_downloads") or []),
            download_count=int(payload.get("download_count", 0)),
            avg_downloads_per_file=float(payload.get("avg_downloads_per_file", 0.0)),
            percentile=int(payload.get("percentile", 0)),
            bot_count=int(payload.get("bot_count", 0)),
            hub_count=int(payload.get("hub_count", 0)),
            organic_count=int(payload.get("organic_count", 0)),
        )


def search(
    keyword: str,
    page_size: int = 100,
    timeout: float | None = 300,
) -> list[PrideProjectSearchResult]:
    """Find PRIDE projects by keyword.

    **The discovery entry point.** Every other function here takes an accession you already have;
    this is the one that produces them, so you can go from a subject to a dataset without leaving
    Python.

    Paging is handled for you: however many pages the result set spans, you get one list, with no
    accession repeated.

    Args:
        keyword: What to search for, e.g. ``"phosphoproteome"``. Matched across titles,
            descriptions, keywords, organisms and more — :attr:`~PrideProjectSearchResult.highlights`
            on each hit says which fields actually matched.
        page_size: How many hits to request per underlying API call. Only affects how the result
            set is fetched, never what you get back.
        timeout: Seconds to allow for the whole fetch.

    Returns:
        Every matching project. **An empty list is a real answer** — PRIDE reports no hits as an
        empty result rather than an error, so unlike :func:`list_files` this does not raise
        :class:`ProjectNotFoundError`: there is no accession here that could have been a typo.

    Raises:
        UsageError: the keyword is blank, over ``MAX_KEYWORD_LENGTH``, or begins with ``-``; or the
            page size is not positive.
        BridgeError: PRIDE returned an error status or was unreachable.

    Note:
        PRIDE pages a **live index** with no stable cursor, so a result set that changes during a
        multi-page fetch shifts its own paging. A project published mid-fetch is served on two pages
        and deduplicated, so it comes back once; a project *removed* mid-fetch can fall between two
        pages and be missed. A search whose hits fit on one page cannot be affected.

    Example:
        >>> hits = search("plasmodium falciparum schizont")        # doctest: +SKIP
        >>> hits[0].accession, hits[0].organisms                   # doctest: +SKIP
        ('PXD070842', ['Homo sapiens (human)', 'Plasmodium falciparum (isolate 3d7)'])
        >>> hits[0].matched_fields                                 # doctest: +SKIP
        ['references', 'title']
        >>> files = list_files(hits[0].accession)                  # doctest: +SKIP
    """
    if not isinstance(keyword, str) or not keyword.strip():
        raise _bridge.UsageError("A search keyword is required, e.g. 'phosphoproteome'.")

    canonical = _reject_flag_like("keyword", keyword.strip())
    if len(canonical) > MAX_KEYWORD_LENGTH:
        raise _bridge.UsageError(
            f"keyword may be at most {MAX_KEYWORD_LENGTH} characters; got {len(canonical)}. "
            "PRIDE answers a longer keyword with HTTP 500, which cannot be told apart from the "
            "service being down."
        )

    if isinstance(page_size, bool) or not isinstance(page_size, int):
        raise _bridge.UsageError(
            f"page_size must be a whole number; got {type(page_size).__name__} ({page_size!r})."
        )
    if page_size <= 0:
        raise _bridge.UsageError(f"page_size must be positive; got {page_size}.")
    if page_size > 2_147_483_647:
        raise _bridge.UsageError(f"page_size is larger than the API allows; got {page_size}.")

    data = _bridge.invoke(
        "pride", "search",
        "--keyword", canonical,
        "--page-size", str(page_size),
        timeout=timeout,
    )
    return [PrideProjectSearchResult._from_wire(item) for item in data.get("results", [])]
