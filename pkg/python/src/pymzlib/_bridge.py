"""Locating and invoking the bundled mzLib bridge executable.

This module is the only place in pyMzLib that knows the bridge exists. Everything above it
sees ordinary Python functions and objects. That boundary is deliberate: the transport
(today, a self-contained .NET executable invoked per call) can be replaced by an in-process
binding or a long-lived local server without any public API changing.

Nothing here imports a third-party package. pyMzLib declares no runtime dependencies, so it
cannot participate in a dependency conflict inside anyone's environment.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

__all__ = [
    "PyMzLibError",
    "BridgeError",
    "ServiceUnavailableError",
    "UsageError",
    "BridgeNotFoundError",
    "bridge_path",
    "invoke",
]

#: Wire-format version this Python package understands. The bridge reports its own; a
#: mismatch means the two halves were built from different sources.
PROTOCOL_VERSION = 1

#: Environment variable pointing at a bridge executable, overriding the bundled one.
#: Used during development, before a self-contained binary has been staged into the package.
BRIDGE_ENV_VAR = "PYMZLIB_BRIDGE"

#: The error type the bridge uses for availability failures. Must match
#: Program.ServiceUnavailableType on the C# side — the two halves agree by this string.
SERVICE_UNAVAILABLE_TYPE = "ServiceUnavailable"


class PyMzLibError(Exception):
    """Base class for every error pyMzLib raises."""


class BridgeNotFoundError(PyMzLibError):
    """The bridge executable could not be located.

    In a released wheel this should be impossible: the executable ships inside the package.
    It normally means pyMzLib is being run from a source checkout where the bridge has not
    been built yet (see ``pkg/build/publish-bridge.ps1``).
    """


class UsageError(PyMzLibError, ValueError):
    """A call was malformed — a missing or invalid argument. Raised before any work happens."""


class BridgeError(PyMzLibError):
    """mzLib reported a failure.

    Attributes:
        error_type: The .NET exception type name, e.g. ``HttpRequestException``. Useful for
            distinguishing a network failure from a bad accession without parsing prose.
    """

    def __init__(self, error_type: str, message: str) -> None:
        super().__init__(message)
        self.error_type = error_type


class ServiceUnavailableError(BridgeError):
    """An external service is unavailable — down, rate-limited, timing out, or unreachable.

    This is deliberately a distinct type, because the difference between "the repository is
    having a bad morning" and "something is broken" is the difference between retrying later and
    filing a bug. Catch it to back off and retry::

        try:
            files = pymzlib.pride.list_files("PXD000001")
        except pymzlib.ServiceUnavailableError:
            ...   # EBI's problem; try again later
        except pymzlib.BridgeError:
            ...   # ours

    The classification happens in the bridge rather than here, so every consumer of the wire
    format gets it and not only Python. HTTP 408, 429, and 5xx count as unavailable; 404 and 400
    do not, because a wrong URL or a malformed request is our problem and excusing it as an
    outage would hide a real bug.
    """


def _platform_tag() -> str:
    """The subdirectory name a bridge build is staged under, matching .NET runtime identifiers."""
    system = platform.system()
    machine = platform.machine().lower()
    arch = {"amd64": "x64", "x86_64": "x64", "arm64": "arm64", "aarch64": "arm64"}.get(machine, machine)
    prefix = {"Windows": "win", "Linux": "linux", "Darwin": "osx"}.get(system)
    if prefix is None:  # pragma: no cover - no other platform is supported
        raise BridgeNotFoundError(f"Unsupported platform: {system} {machine}")
    return f"{prefix}-{arch}"


def bridge_path() -> Path:
    """Return the path of the bridge executable that will be used.

    Resolution order: the ``PYMZLIB_BRIDGE`` environment variable, then the copy staged
    inside this package for the current platform.

    Raises:
        BridgeNotFoundError: if neither exists.
    """
    override = os.environ.get(BRIDGE_ENV_VAR)
    if override:
        candidate = Path(override)
        if not candidate.is_file():
            raise BridgeNotFoundError(f"{BRIDGE_ENV_VAR} points at '{override}', which is not a file.")
        return candidate

    executable = "mzlib-bridge.exe" if sys.platform == "win32" else "mzlib-bridge"
    candidate = Path(__file__).parent / "_dotnet" / _platform_tag() / executable
    if not candidate.is_file():
        raise BridgeNotFoundError(
            f"No mzLib bridge for this platform at '{candidate}'. "
            f"In a source checkout, build one and set {BRIDGE_ENV_VAR} to its path."
        )
    return candidate


def invoke(*args: str, timeout: float | None = None) -> Any:
    """Run one bridge command and return the decoded ``data`` payload.

    Args:
        *args: The command and its options, e.g. ``"pride", "files", "--accession", "PXD000001"``.
        timeout: Seconds to wait before giving up. ``None`` waits indefinitely, which is the
            right default for a large download.

    Returns:
        Whatever the command's ``data`` field contains — usually a dict.

    Raises:
        UsageError: the command or its arguments were malformed.
        BridgeError: mzLib itself failed (network, bad accession, disk).
        PyMzLibError: the bridge produced output this version cannot interpret.
    """
    command = [str(bridge_path()), *args]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ServiceUnavailableError(
            "Timeout", f"mzLib bridge timed out after {timeout}s."
        ) from exc

    stdout = completed.stdout.strip()
    if not stdout:
        # A silent non-zero exit means the process died before it could report anything —
        # surface stderr, which is the only evidence left.
        raise PyMzLibError(
            f"mzLib bridge exited with code {completed.returncode} and no output. "
            f"stderr: {completed.stderr.strip() or '(empty)'}"
        )

    try:
        envelope = json.loads(stdout)
    except json.JSONDecodeError as exc:
        raise PyMzLibError(f"mzLib bridge returned output that is not JSON: {stdout[:400]}") from exc

    if envelope.get("ok"):
        return envelope.get("data")

    error = envelope.get("error") or {}
    error_type = error.get("type", "Unknown")
    message = error.get("message", "mzLib reported a failure with no message.")
    if error_type == "usage":
        raise UsageError(message)
    if error_type == SERVICE_UNAVAILABLE_TYPE:
        raise ServiceUnavailableError(error_type, message)
    raise BridgeError(error_type, message)


def bridge_version() -> dict[str, Any]:
    """Return the bridge's own version information, and check protocol compatibility.

    Raises:
        PyMzLibError: if the bridge speaks a different wire format than this package.
    """
    info = invoke("version", timeout=60)
    reported = info.get("protocol")
    if reported != PROTOCOL_VERSION:
        raise PyMzLibError(
            f"mzLib bridge speaks protocol {reported}, but this pyMzLib expects {PROTOCOL_VERSION}. "
            "The Python package and the bridge were built from different sources."
        )
    return info
