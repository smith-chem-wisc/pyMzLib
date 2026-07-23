"""Identify proteomics result files: what a file *is*, and what you can do with it.

mzLib recognises 29 file types written by a dozen different search and deconvolution tools —
MetaMorpheus, MSFragger, TopPIC, TopFD, MsPathFinderT, Crux, Casanovo, FlashDeconv, Dinosaur,
FlashLFQ — and dispatches each to a parser it maintains. This module asks it what a path is::

    >>> import pymzlib
    >>> info = pymzlib.readers.identify("psm.tsv")     # doctest: +SKIP
    >>> info.file_type, info.views                     # doctest: +SKIP
    ('MsFraggerPsm', ['quantifiable'])

**Read** :attr:`FileInfo.views` **before assuming anything.** It is tempting to describe mzLib as
reading 29 formats into one uniform shape, and it does not. The formats fall into disjoint families
and several belong to no family at all:

- ``"quantifiable"`` — a cross-format record view (sequence, retention time, charge, mass, protein
  groups), and the input :func:`pymzlib.flashlfq.quantify` accepts. **Exactly three file types have
  it**: MetaMorpheus ``.psmtsv`` and ``.osmtsv``, and MSFragger ``psm.tsv``.
- ``"ms1_features"`` — deconvolved MS1 features (TopFD ``_ms1.feature``, Dinosaur).
- ``"spectra"`` — the file is spectra, not results (``.raw``, ``.mzML``, ``.mgf``, ``.d``, msalign).
- ``"spectral_match"`` — records are identifications but share no file-level interface
  (MsPathFinderT, Casanovo).
- ``[]`` — **an empty list is a real and common answer.** TopPIC, Crux, MSFragger's peptide/protein
  tables and the FlashDeconv formats each parse into their own record type with nothing in common.
  mzLib reads them; there is simply no uniform view to project them onto.

Call :func:`formats` for the whole table. It is enumerated from mzLib rather than transcribed, so it
cannot drift from what mzLib actually dispatches.

**Two things this module deliberately does not tell you**, in the "surface it, don't hide it" spirit
of the rest of pyMzLib:

- **Which tool wrote the file.** mzLib has a ``Software`` property that looks like the answer and is
  not: readers carry their software constant on a constructor that mzLib's own file factory does not
  use, so the value is unset for everything the factory returns — and it is not reliably set on the
  other constructor either. Rather than reconstruct a plausible answer, there is no ``software``
  field. :attr:`FileInfo.file_type` already names the tool.
- **Whether the numbers inside mean the same thing across formats.** They do not, and this is the
  trap most likely to produce a wrong result. mzLib's result-file readers pass through whatever the
  tool wrote, with no unit conversion: MetaMorpheus retention times are in **minutes**, MSFragger's
  and TopPIC's are in **seconds**, and TopFD changed from seconds to minutes between v1.6.2 and
  v1.7.0 *within the same file type*. Likewise ``is_decoy`` is hardcoded ``False`` for MSFragger
  (mzLib does not read its decoys), and MSFragger's "monoisotopic mass" is the *theoretical* peptide
  mass while MetaMorpheus's is the observed one. Identifying a file is safe; comparing raw fields
  across formats is not.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import _bridge

__all__ = [
    "Format",
    "FileInfo",
    "formats",
    "identify",
]

#: The view name for the cross-format record shape :func:`pymzlib.flashlfq.quantify` consumes.
QUANTIFIABLE = "quantifiable"


@dataclass(frozen=True)
class Format:
    """One file type mzLib can recognise.

    Attributes:
        file_type: mzLib's ``SupportedFileType`` name, e.g. ``"MsFraggerPsm"``, ``"psmtsv"``.
        extension: The extension or filename suffix mzLib dispatches on, e.g. ``"psm.tsv"``,
            ``"_ms1.feature"``. **Not unique across file types** — ``BrukerD`` and
            ``BrukerTimsTof`` are both ``.d`` (told apart by the directory's contents), and
            several formats share ``.tsv``.
        reader: The name of the mzLib class that parses it, for cross-referencing the mzLib source.
        views: The uniform views this format supports — see the module docstring. Often empty.
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
        views: The uniform views this file supports — see the module docstring. Often empty, which
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
        ``psms``. When ``False``, mzLib can still read the file — it simply has no uniform view,
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
    more than identifying an empty one. It is not, however, *pure* — mzLib disambiguates a bare
    ``.tsv`` by reading its first line, a ``.mztab`` by its first five, and a Bruker ``.d`` by which
    analysis file the directory holds. An unreadable file will therefore raise.

    Args:
        path: Path to a result or spectra file. A Bruker ``.d`` directory is also accepted.
        timeout: Seconds to allow.

    Returns:
        A :class:`FileInfo` naming the format and the views it supports.

    Raises:
        UsageError: the path is blank, does not exist, or is not a file type mzLib recognises.
            mzLib has no "unknown" result — a file is dispatchable or it is an error — so use
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
