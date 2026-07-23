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
    "BridgeTimeoutError",
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


class BridgeTimeoutError(PyMzLibError):
    """The bridge process did not finish within the timeout.

    Deliberately **not** a :class:`ServiceUnavailableError`, and the distinction is the whole
    point. A subprocess timeout has several possible causes and only one of them is a slow
    service: the bridge may be wedged, the executable may be corrupt, antivirus may be holding
    it, or the caller may simply have passed a timeout that was too short. Reporting all of that
    as "the repository is down" is how a real bug gets skipped by every test suite and never
    seen again.

    If the caller wants a slow network to be treated as an outage, they can catch this
    explicitly — but the library will not guess on their behalf.
    """


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


def _validate_timeout(timeout: float | None) -> None:
    """Reject timeouts that cannot mean anything, before spawning a process.

    ``subprocess`` accepts ``0``, negatives, ``inf`` and ``nan`` and then fails in ways that
    point nowhere near the caller's mistake — ``inf`` raises ``OverflowError`` from deep inside
    the platform's clock, and ``0`` looks exactly like a service that never answered.
    """
    if timeout is None:
        return
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise UsageError(f"timeout must be a number of seconds or None; got {timeout!r}.")
    if timeout != timeout:  # NaN
        raise UsageError("timeout must be a number of seconds or None; got nan.")
    if timeout <= 0:
        raise UsageError(f"timeout must be greater than zero; got {timeout!r}.")
    if timeout == float("inf"):
        raise UsageError("timeout must be finite; pass None to wait indefinitely.")


def invoke(*args: str, stdin: str | None = None, timeout: float | None = None) -> Any:
    """Run one bridge command and return the decoded ``data`` payload.

    Args:
        *args: The command and its options, e.g. ``"pride", "files", "--accession", "PXD000001"``.
        stdin: Text to write to the process's standard input. Used when a payload would not fit
            on the command line — argv has a hard ceiling of roughly 32 KB, and a few thousand
            file names exceed it.
        timeout: Seconds to wait before giving up. ``None`` waits indefinitely, which is the
            right default for a large download.

    Returns:
        Whatever the command's ``data`` field contains — usually a dict.

    Raises:
        UsageError: the command or its arguments were malformed.
        BridgeError: mzLib itself failed (network, bad accession, disk).
        PyMzLibError: the bridge produced output this version cannot interpret.
    """
    _validate_timeout(timeout)
    for arg in args:
        if not isinstance(arg, str):
            raise UsageError(
                f"Internal error: bridge arguments must be strings, got {type(arg).__name__} "
                f"({arg!r}). This is a bug in pyMzLib, not in your call."
            )
        if "\x00" in arg:
            raise UsageError("Arguments may not contain a null character.")
    command = [str(bridge_path()), *args]
    try:
        completed = subprocess.run(
            command,
            input=stdin,
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise BridgeTimeoutError(
            f"mzLib bridge did not finish within {timeout}s. This may mean the service is slow, "
            "but it can equally mean the bridge is wedged or the timeout was too short — "
            "pyMzLib will not guess which."
        ) from exc
    except OSError as exc:
        # A missing execute bit, a quarantined binary, or a file that is not an executable at
        # all. Without this the caller sees a bare OSError, contradicting the promise that
        # PyMzLibError is the base class of everything this package raises.
        raise PyMzLibError(
            f"Could not run the mzLib bridge at '{command[0]}': {exc}"
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
