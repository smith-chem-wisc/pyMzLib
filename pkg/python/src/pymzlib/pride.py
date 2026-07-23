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

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import _bridge

__all__ = ["PrideFile", "list_files", "download", "total_size_bytes"]


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

    @classmethod
    def _from_wire(cls, payload: dict[str, Any]) -> "PrideFile":
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
    if not accession or not accession.strip():
        raise _bridge.UsageError("A PRIDE project accession is required, e.g. 'PXD000001'.")
    if page_size <= 0:
        raise _bridge.UsageError(f"page_size must be positive; got {page_size}.")

    data = _bridge.invoke(
        "pride", "files",
        "--accession", accession.strip(),
        "--page-size", str(page_size),
        timeout=timeout,
    )
    return [PrideFile._from_wire(item) for item in data.get("files", [])]


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
    if not accession or not accession.strip():
        raise _bridge.UsageError("A PRIDE project accession is required, e.g. 'PXD000001'.")
    destination = Path(destination)

    args = ["pride", "download", "--accession", accession.strip(), "--dest", str(destination)]
    if category:
        args += ["--category", category]
    if extensions:
        args += ["--ext", ",".join(extensions)]
    if not overwrite:
        args.append("--no-overwrite")

    data = _bridge.invoke(*args, timeout=timeout)
    return [Path(p) for p in data.get("paths", [])]


def total_size_bytes(files: Iterable[PrideFile]) -> int:
    """Sum the sizes of some files — e.g. to see what a download will cost before starting it.

    >>> files = list_files("PXD000001")           # doctest: +SKIP
    >>> total_size_bytes(f for f in files if f.category == "RAW") / 1e9   # doctest: +SKIP
    0.51
    """
    return sum(f.file_size_bytes for f in files)
