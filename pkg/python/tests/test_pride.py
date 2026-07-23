"""Tests for the PRIDE surface.

Two tiers, deliberately separated:

* **Offline** — the bridge is replaced by a recorded payload, so these run anywhere, in
  milliseconds, with no network. They test the Python layer's own behavior.
* **Network** (``-m network``) — real calls to the live PRIDE Archive through the real
  bridge. These are the ones that would catch mzLib or PRIDE changing under us.

Run the fast set with ``pytest -m "not network"``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, pride

FIXTURE = Path(__file__).parent / "fixtures" / "pride_PXD000001_files.json"


@pytest.fixture()
def recorded_manifest(monkeypatch):
    """Serve a recorded PXD000001 manifest instead of calling the bridge."""
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))

    def fake_invoke(*args, timeout=None):
        return payload

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return payload


# --------------------------------------------------------------------------- offline


def test_list_files_parses_every_file(recorded_manifest):
    files = pride.list_files("PXD000001")
    assert len(files) == len(recorded_manifest["files"])
    assert all(isinstance(f, pride.PrideFile) for f in files)


def test_file_fields_are_typed(recorded_manifest):
    first = pride.list_files("PXD000001")[0]
    assert isinstance(first.file_size_bytes, int)
    assert first.submission_date is not None
    assert first.submission_date.tzinfo is not None, "timestamps must be timezone-aware"


def test_derived_properties(recorded_manifest):
    files = pride.list_files("PXD000001")
    fasta = next(f for f in files if f.file_name.endswith(".fasta"))
    assert fasta.extension == ".fasta"
    assert fasta.size_mb == pytest.approx(fasta.file_size_bytes / 1_000_000)
    assert fasta.downloadable is True


def test_total_size_matches_sum(recorded_manifest):
    files = pride.list_files("PXD000001")
    assert pride.total_size_bytes(files) == sum(f.file_size_bytes for f in files)


def test_aspera_only_file_is_not_downloadable():
    """A file with no HTTPS location must report itself as not downloadable, not crash."""
    aspera_only = pride.PrideFile._from_wire(
        {"file_name": "x.raw", "file_size_bytes": 1, "https_url": None}
    )
    assert aspera_only.downloadable is False


@pytest.mark.parametrize("accession", ["", "   ", None])
def test_blank_accession_rejected_before_any_work(accession):
    with pytest.raises(pymzlib.UsageError):
        pride.list_files(accession)


def test_nonpositive_page_size_rejected():
    with pytest.raises(pymzlib.UsageError):
        pride.list_files("PXD000001", page_size=0)


def test_bridge_failure_surfaces_as_bridge_error(monkeypatch):
    """A .NET-side failure must arrive as a typed Python exception, not a subprocess artifact."""

    def failing_invoke(*args, timeout=None):
        raise _bridge.BridgeError("HttpRequestException", "PRIDE Archive request failed with status 503")

    monkeypatch.setattr(_bridge, "invoke", failing_invoke)
    with pytest.raises(pymzlib.BridgeError) as caught:
        pride.list_files("PXD000001")
    assert caught.value.error_type == "HttpRequestException"


def test_malformed_timestamp_becomes_none_rather_than_crashing():
    """PRIDE metadata is not ours to control; one odd date must not sink a whole manifest."""
    odd = pride.PrideFile._from_wire(
        {"file_name": "x.raw", "submission_date": "not a date", "publication_date": None}
    )
    assert odd.submission_date is None
    assert odd.publication_date is None


# --------------------------------------------------------------------------- download, offline
#
# The live download test covers the real path, but it is slow and routinely deselected. These
# cover what is actually easy to get wrong and invisible when it is: how the call is assembled.


@pytest.fixture()
def captured_args(monkeypatch):
    """Capture the argument list `download` builds instead of running anything."""
    seen = {}

    def fake_invoke(*args, timeout=None):
        seen["args"] = list(args)
        seen["timeout"] = timeout
        return {"paths": ["out/a.raw", "out/b.raw"]}

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


def test_download_returns_paths(captured_args, tmp_path):
    written = pride.download("PXD000001", tmp_path)
    assert [p.name for p in written] == ["a.raw", "b.raw"]
    assert all(isinstance(p, Path) for p in written)


def test_download_passes_accession_and_destination(captured_args, tmp_path):
    pride.download("  PXD000001  ", tmp_path)
    args = captured_args["args"]
    assert args[:2] == ["pride", "download"]
    assert args[args.index("--accession") + 1] == "PXD000001", "accession must be trimmed"
    assert args[args.index("--dest") + 1] == str(tmp_path)


def test_download_omits_filters_when_not_asked_for(captured_args, tmp_path):
    pride.download("PXD000001", tmp_path)
    assert "--category" not in captured_args["args"]
    assert "--ext" not in captured_args["args"]


def test_download_passes_category(captured_args, tmp_path):
    pride.download("PXD000001", tmp_path, category="RAW")
    args = captured_args["args"]
    assert args[args.index("--category") + 1] == "RAW"


def test_download_joins_extensions_with_commas(captured_args, tmp_path):
    pride.download("PXD000001", tmp_path, extensions=[".raw", ".mzML"])
    args = captured_args["args"]
    assert args[args.index("--ext") + 1] == ".raw,.mzML"


def test_overwrite_flag_is_not_inverted(captured_args, tmp_path):
    """The flag is negative on the wire and positive in Python; getting this backwards would
    silently re-download entire projects, or silently skip files the user wanted refreshed."""
    pride.download("PXD000001", tmp_path, overwrite=True)
    assert "--no-overwrite" not in captured_args["args"]

    pride.download("PXD000001", tmp_path, overwrite=False)
    assert "--no-overwrite" in captured_args["args"]


def test_download_defaults_to_no_timeout(captured_args, tmp_path):
    """Multi-gigabyte transfers legitimately take hours; a default timeout would truncate them."""
    pride.download("PXD000001", tmp_path)
    assert captured_args["timeout"] is None


@pytest.mark.parametrize("accession", ["", "   ", None])
def test_download_rejects_blank_accession_before_touching_the_network(accession, tmp_path):
    with pytest.raises(pymzlib.UsageError):
        pride.download(accession, tmp_path)


# --------------------------------------------------------------------------- network


@pytest.mark.network
def test_bridge_protocol_matches():
    info = pymzlib.bridge_version()
    assert info["protocol"] == _bridge.PROTOCOL_VERSION


@pytest.mark.network
def test_live_manifest_matches_recording():
    """The live manifest should still contain what we recorded — a canary on PRIDE and mzLib."""
    live = pride.list_files("PXD000001")
    assert len(live) >= 8
    assert any(f.category == "RAW" for f in live)


@pytest.mark.network
def test_unknown_accession_returns_empty_not_error():
    """PRIDE's own behavior for an unknown accession is an empty result; preserve it."""
    assert pride.list_files("PXD999999999") == []


@pytest.mark.network
@pytest.mark.slow
def test_download_writes_the_file(tmp_path):
    written = pride.download("PXD000001", tmp_path, extensions=[".fasta"])
    assert len(written) == 1
    assert written[0].is_file()
    assert written[0].stat().st_size > 0
    assert not list(tmp_path.glob("*.partial")), "no partial file may survive a successful download"
