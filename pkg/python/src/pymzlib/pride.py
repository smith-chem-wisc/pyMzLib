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
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import _bridge

__all__ = [
    "PrideFile",
    "ProjectNotFoundError",
    "list_files",
    "download",
    "download_files",
    "total_size_bytes",
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


def list_files(accession: str, page_size: int = 100, timeout: float | None = 300) -> list[PrideFile]:
    """Return the complete file manifest of a PRIDE Archive project.

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
    """Sum the sizes of some files — e.g. to see what a download will cost before starting it.

    >>> files = list_files("PXD000001")           # doctest: +SKIP
    >>> total_size_bytes(f for f in files if f.category == "RAW") / 1e9   # doctest: +SKIP
    0.51
    """
    return sum(f.file_size_bytes for f in files)
