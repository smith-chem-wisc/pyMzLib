"""pyMzLib — mzLib for Python.

mzLib (https://github.com/smith-chem-wisc/mzLib) is a mass-spectrometry and proteomics
library written in C#. pyMzLib makes its functionality callable from Python, with no .NET
installation and no third-party Python dependencies: everything needed ships inside the
package.

The first area covered is the PRIDE Archive::

    import pymzlib

    files = pymzlib.pride.list_files("PXD000001")
    print(f"{len(files)} files, {pymzlib.pride.total_size_bytes(files) / 1e9:.2f} GB")

    pymzlib.pride.download("PXD000001", "downloads", category="RAW")
"""

from . import flashlfq, peptidoform, pride
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
)

__version__ = "0.1.0.dev0"

__all__ = [
    "flashlfq",
    "peptidoform",
    "pride",
    "PyMzLibError",
    "BridgeError",
    "ServiceUnavailableError",
    "BridgeTimeoutError",
    "ProjectNotFoundError",
    "BridgeNotFoundError",
    "UsageError",
    "bridge_path",
    "bridge_version",
    "__version__",
]
