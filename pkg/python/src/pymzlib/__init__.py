"""pyMzLib — mzLib for Python.

mzLib (https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry and proteomics
library written in C#. pyMzLib makes its functionality callable from Python, with no .NET
installation and no third-party Python dependencies: everything needed ships inside the
package.

Read a mass-spectrometry data file - **mzML**, Thermo ``.raw``, Bruker ``.d``, timsTOF
``.d``, MGF or msalign - with :func:`pymzlib.readers.read_spectra`::

    import pymzlib

    scans = pymzlib.readers.read_spectra("run.mzML", ms_order=2, limit=5, peaks=True)
    print(scans.scan_count, scans.columns["retention_time"])

The same module identifies and reads all 38 file types mzLib knows, search results included.
The PRIDE Archive is covered too::

    files = pymzlib.pride.list_files("PXD000001")
    print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")

    pymzlib.pride.download("PXD000001", "downloads", category="RAW")

And a differential-abundance test is limma's moderated t, with no R, in
:func:`pymzlib.stats.fit`::

    fit = pymzlib.stats.fit("log2_intensities.tsv", "design.tsv", ["treated"])
    print(fit.status_counts, fit.rows("treated")[0]["bh_adjusted"])
"""

from . import flashlfq, isobaric, peptidoform, pride, proteins, readers, sdrf, stats
from .pride import ProjectNotFoundError
from ._bridge import (
    BridgeError,
    ServiceUnavailableError,
    BridgeTimeoutError,
    BridgeNotFoundError,
    PyMzLibError,
    UsageError,
    bridge_path,
    bridge_version,
    BridgeVersion,
)

__version__ = "0.2.0"

__all__ = [
    "flashlfq",
    "isobaric",
    "peptidoform",
    "pride",
    "proteins",
    "readers",
    "sdrf",
    "stats",
    "PyMzLibError",
    "BridgeError",
    "ServiceUnavailableError",
    "BridgeTimeoutError",
    "ProjectNotFoundError",
    "BridgeNotFoundError",
    "UsageError",
    "bridge_path",
    "bridge_version",
    "BridgeVersion",
    "__version__",
]
