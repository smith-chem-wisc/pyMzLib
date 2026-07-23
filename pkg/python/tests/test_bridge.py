"""Tests for the transport layer itself.

`_bridge` is where pyMzLib meets the outside world, so its failure paths are the ones a user
actually hits on a bad day: the executable is missing, the process dies, the output isn't JSON,
the version doesn't match. Those paths are invisible to the PRIDE tests, which only ever see the
happy path, and they are exactly the code most likely to be wrong — a wrong error message wastes
somebody's afternoon.

None of these tests touch the network or the real executable.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from pymzlib import _bridge


class _Completed:
    """Stands in for `subprocess.CompletedProcess` without running anything."""

    def __init__(self, stdout: str = "", stderr: str = "", returncode: int = 0) -> None:
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode


@pytest.fixture()
def fake_bridge(monkeypatch, tmp_path):
    """Point `bridge_path()` at a file that exists but is never actually executed."""
    stub = tmp_path / ("mzlib-bridge.exe" if sys.platform == "win32" else "mzlib-bridge")
    stub.write_text("not a real executable")
    monkeypatch.setattr(_bridge, "bridge_path", lambda: stub)
    return stub


def _run_returns(monkeypatch, completed):
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: completed)


# --------------------------------------------------------------------------- locating the payload


def test_env_override_takes_precedence(monkeypatch, tmp_path):
    stub = tmp_path / "custom-bridge"
    stub.write_text("x")
    monkeypatch.setenv(_bridge.BRIDGE_ENV_VAR, str(stub))
    assert _bridge.bridge_path() == stub


def test_env_override_pointing_at_nothing_says_so(monkeypatch, tmp_path):
    monkeypatch.setenv(_bridge.BRIDGE_ENV_VAR, str(tmp_path / "does-not-exist"))
    with pytest.raises(_bridge.BridgeNotFoundError, match="not a file"):
        _bridge.bridge_path()


def test_missing_payload_names_the_path_and_the_way_out(monkeypatch):
    """The message has to be actionable: a bare 'not found' leaves the user with nowhere to go."""
    monkeypatch.delenv(_bridge.BRIDGE_ENV_VAR, raising=False)
    monkeypatch.setattr(Path, "is_file", lambda self: False)
    with pytest.raises(_bridge.BridgeNotFoundError) as caught:
        _bridge.bridge_path()
    message = str(caught.value)
    assert "_dotnet" in message
    assert _bridge.BRIDGE_ENV_VAR in message


@pytest.mark.parametrize(
    ("system", "machine", "expected"),
    [
        ("Windows", "AMD64", "win-x64"),
        ("Linux", "x86_64", "linux-x64"),
        ("Linux", "aarch64", "linux-arm64"),
        ("Darwin", "arm64", "osx-arm64"),
        ("Darwin", "x86_64", "osx-x64"),
    ],
)
def test_platform_tags_match_dotnet_runtime_identifiers(monkeypatch, system, machine, expected):
    """These strings must equal the RIDs publish-bridge.ps1 stages under, or nothing is found."""
    monkeypatch.setattr("platform.system", lambda: system)
    monkeypatch.setattr("platform.machine", lambda: machine)
    assert _bridge._platform_tag() == expected


def test_unsupported_platform_is_reported_not_guessed(monkeypatch):
    monkeypatch.setattr("platform.system", lambda: "Plan9")
    monkeypatch.setattr("platform.machine", lambda: "x86_64")
    with pytest.raises(_bridge.BridgeNotFoundError, match="Unsupported platform"):
        _bridge._platform_tag()


# --------------------------------------------------------------------------- the envelope


def test_success_returns_only_the_data(monkeypatch, fake_bridge):
    _run_returns(monkeypatch, _Completed(stdout=json.dumps({"ok": True, "data": {"a": 1}})))
    assert _bridge.invoke("anything") == {"a": 1}


def test_usage_failure_maps_to_usage_error(monkeypatch, fake_bridge):
    envelope = {"ok": False, "error": {"type": "usage", "message": "Missing required option --x."}}
    _run_returns(monkeypatch, _Completed(stdout=json.dumps(envelope), returncode=2))
    with pytest.raises(_bridge.UsageError, match="Missing required option"):
        _bridge.invoke("anything")


def test_runtime_failure_preserves_the_dotnet_error_type(monkeypatch, fake_bridge):
    """`error_type` is what lets a caller tell a network blip from a bad accession."""
    envelope = {"ok": False, "error": {"type": "HttpRequestException", "message": "503"}}
    _run_returns(monkeypatch, _Completed(stdout=json.dumps(envelope), returncode=1))
    with pytest.raises(_bridge.BridgeError) as caught:
        _bridge.invoke("anything")
    assert caught.value.error_type == "HttpRequestException"


def test_failure_with_no_message_still_raises_something_readable(monkeypatch, fake_bridge):
    _run_returns(monkeypatch, _Completed(stdout=json.dumps({"ok": False}), returncode=1))
    with pytest.raises(_bridge.BridgeError, match="no message"):
        _bridge.invoke("anything")


# --------------------------------------------------------------------------- when things go wrong


def test_silent_death_surfaces_stderr(monkeypatch, fake_bridge):
    """A process that dies before writing anything leaves stderr as the only evidence."""
    _run_returns(monkeypatch, _Completed(stdout="", stderr="Segmentation fault", returncode=139))
    with pytest.raises(_bridge.PyMzLibError) as caught:
        _bridge.invoke("anything")
    assert "139" in str(caught.value)
    assert "Segmentation fault" in str(caught.value)


def test_silent_death_with_no_stderr_does_not_produce_an_empty_message(monkeypatch, fake_bridge):
    _run_returns(monkeypatch, _Completed(stdout="", stderr="", returncode=1))
    with pytest.raises(_bridge.PyMzLibError, match=r"\(empty\)"):
        _bridge.invoke("anything")


def test_non_json_output_is_quoted_back_not_swallowed(monkeypatch, fake_bridge):
    """If the bridge prints something unexpected, the user needs to see what it was."""
    _run_returns(monkeypatch, _Completed(stdout="Unhandled exception. System.Whatever"))
    with pytest.raises(_bridge.PyMzLibError, match="not JSON"):
        _bridge.invoke("anything")


def test_timeout_becomes_a_typed_error_not_a_subprocess_artifact(monkeypatch, fake_bridge):
    """Nothing above _bridge should ever have to import subprocess to handle a failure."""

    def timing_out(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="mzlib-bridge", timeout=5)

    monkeypatch.setattr(subprocess, "run", timing_out)
    with pytest.raises(_bridge.BridgeTimeoutError):
        _bridge.invoke("anything", timeout=5)


def test_a_timeout_is_not_reported_as_a_service_outage(monkeypatch, fake_bridge):
    """The regression this guards is subtle and was live: every subprocess timeout used to raise
    ServiceUnavailableError, which the canary suites turn into a skip. A wedged bridge, a corrupt
    binary, or a caller passing too small a timeout all reported "EBI is down" and the live tests
    passed green — precisely the failure this project's testing convention exists to prevent."""

    def timing_out(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="mzlib-bridge", timeout=1)

    monkeypatch.setattr(subprocess, "run", timing_out)
    with pytest.raises(_bridge.PyMzLibError) as caught:
        _bridge.invoke("anything", timeout=1)
    assert not isinstance(caught.value, _bridge.ServiceUnavailableError)


@pytest.mark.parametrize("bad", [0, -1, -5.0, float("inf"), float("nan"), "5", [], True])
def test_unusable_timeouts_are_rejected_before_spawning_anything(bad, fake_bridge, monkeypatch):
    """subprocess accepts these and then fails somewhere unrecognisable — inf raises OverflowError
    from the platform clock, and 0 is indistinguishable from a service that never answered."""
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: pytest.fail("should not have run"))
    with pytest.raises(_bridge.UsageError):
        _bridge.invoke("anything", timeout=bad)


def test_an_unlaunchable_bridge_is_a_pymzlib_error(monkeypatch, fake_bridge):
    """A missing execute bit or a quarantined binary must not escape as a bare OSError."""

    def not_executable(*args, **kwargs):
        raise OSError(8, "Exec format error")

    monkeypatch.setattr(subprocess, "run", not_executable)
    with pytest.raises(_bridge.PyMzLibError, match="Could not run the mzLib bridge"):
        _bridge.invoke("anything")


def test_null_bytes_are_rejected_rather_than_raising_from_subprocess(fake_bridge):
    with pytest.raises(_bridge.UsageError, match="null character"):
        _bridge.invoke("pride", "files", "--accession", "PX\x00D")


def test_every_error_is_catchable_as_one_type():
    """A user should be able to write `except PyMzLibError` and be done."""
    for error in (_bridge.UsageError, _bridge.BridgeNotFoundError, _bridge.BridgeError):
        assert issubclass(error, _bridge.PyMzLibError)


# --------------------------------------------------------------------------- version handshake


def test_matching_protocol_returns_the_version_info(monkeypatch, fake_bridge):
    info = {"bridge": "1.0.0.0", "protocol": _bridge.PROTOCOL_VERSION, "runtime": "8.0.27"}
    _run_returns(monkeypatch, _Completed(stdout=json.dumps({"ok": True, "data": info})))
    assert _bridge.bridge_version() == info


def test_mismatched_protocol_fails_loudly(monkeypatch, fake_bridge):
    """Halves built from different sources must not silently produce wrong results."""
    info = {"bridge": "1.0.0.0", "protocol": _bridge.PROTOCOL_VERSION + 1}
    _run_returns(monkeypatch, _Completed(stdout=json.dumps({"ok": True, "data": info})))
    with pytest.raises(_bridge.PyMzLibError, match="different sources"):
        _bridge.bridge_version()
