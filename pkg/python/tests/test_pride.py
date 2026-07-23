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


# --------------------------------------------------------------------------- input that used to
# be accepted silently
#
# Every test below corresponds to a way the library previously reported success while doing
# something other than what was asked. Silent wrong answers are the worst class of bug in a
# library like this, because the caller has no reason to look.


@pytest.mark.parametrize("accession", ["pxd000001", "  PXD000001  ", "Pxd000001"])
def test_accession_case_and_whitespace_are_normalised(recorded_manifest, accession):
    """PRIDE's API is case-sensitive on the accession while our category matching is not. A
    lowercase accession used to return [] — indistinguishable from an empty project."""
    assert pride.list_files(accession)


@pytest.mark.parametrize("accession", ["banana", "PXD", "12345", "PXD00", "PXD000001x", "-PXD1"])
def test_malformed_accessions_are_rejected_not_silently_empty(accession):
    with pytest.raises(pymzlib.UsageError):
        pride.list_files(accession)


@pytest.mark.parametrize("accession", [123, None, ["PXD000001"], b"PXD000001"])
def test_non_string_accession_gives_a_usage_error_not_an_attribute_error(accession):
    with pytest.raises(pymzlib.UsageError):
        pride.list_files(accession)


def test_a_valid_but_unknown_accession_raises_rather_than_returning_empty(monkeypatch):
    """PRIDE answers an unknown accession with an empty result. Passing that through as [] meant
    a typo produced '0 files, done' and a script that carried on."""
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"files": []})
    with pytest.raises(pride.ProjectNotFoundError, match="PXD999999999"):
        pride.list_files("PXD999999999")


@pytest.mark.parametrize("destination", ["", "   ", None, 123])
def test_blank_destination_is_refused_instead_of_writing_to_the_cwd(destination):
    """Path("") is Path("."), so dest = cfg.get("outdir", "") sprayed a project into os.getcwd()."""
    with pytest.raises(pymzlib.UsageError):
        pride.download("PXD000001", destination)


def test_a_bare_string_of_extensions_is_refused(tmp_path):
    """`.raw` iterates as four characters, matched nothing, and exited successfully."""
    with pytest.raises(pymzlib.UsageError, match=r"\[\'\.raw\'\]"):
        pride.download("PXD000001", tmp_path, extensions=".raw")


def test_a_blank_category_is_refused_rather_than_selecting_everything(tmp_path):
    """A degenerate filter used to leave the bridge filter null and download the whole project."""
    with pytest.raises(pymzlib.UsageError):
        pride.download("PXD000001", tmp_path, category="   ")


@pytest.mark.parametrize("value", ["--no-overwrite", "-x"])
def test_flag_like_filter_values_are_refused(tmp_path, value):
    """The bridge parser reads `--category --no-overwrite` as two flags, dropping the category and
    silently enabling a flag the caller never asked for."""
    with pytest.raises(pymzlib.UsageError, match="another option"):
        pride.download("PXD000001", tmp_path, category=value)


@pytest.mark.parametrize("page_size", ["100", None, 2.5, 2_147_483_648])
def test_bad_page_sizes_are_usage_errors(page_size):
    with pytest.raises(pymzlib.UsageError):
        pride.list_files("PXD000001", page_size=page_size)


# --------------------------------------------------------------------------- download_files
#
# The workflow now closes: what list_files() returns can be handed straight to download_files().


def test_as_dict_includes_the_computed_properties(recorded_manifest):
    """`vars(f)` silently omits size_mb / extension / downloadable — the three attributes the
    documentation pushes hardest, including the one the FAQ tells you to filter on."""
    record = pride.list_files("PXD000001")[0].as_dict()
    for key in ("size_mb", "extension", "downloadable", "file_name"):
        assert key in record


def test_files_know_which_project_they_came_from(recorded_manifest):
    assert all(f.project_accession == "PXD000001" for f in pride.list_files("PXD000001"))


def test_download_files_sends_the_selection_on_stdin(recorded_manifest, monkeypatch, tmp_path):
    seen = {}

    def fake_invoke(*args, stdin=None, timeout=None):
        seen["args"] = list(args)
        seen["stdin"] = stdin
        return {"paths": ["out/a"]}

    files = pride.list_files("PXD000001")
    chosen = [f for f in files if f.size_mb < 5 and f.downloadable]
    monkeypatch.setattr(_bridge, "invoke", fake_invoke)

    pride.download_files(chosen, tmp_path)

    assert "--names-from-stdin" in seen["args"]
    assert seen["stdin"].splitlines() == [f.file_name for f in chosen]
    # The names must not be on argv: a few thousand of them would exceed the ~32 KB ceiling.
    assert not any(f.file_name in " ".join(seen["args"]) for f in chosen)


def test_download_files_refuses_an_empty_selection(tmp_path):
    """An empty selection is nearly always a filter that did not match what the caller expected."""
    with pytest.raises(pymzlib.UsageError, match="No files selected"):
        pride.download_files([], tmp_path)


def test_download_files_refuses_files_that_cannot_be_fetched(recorded_manifest, tmp_path):
    aspera_only = pride.PrideFile._from_wire(
        {"file_name": "x.raw", "file_size_bytes": 1, "https_url": None}, "PXD000001"
    )
    with pytest.raises(pymzlib.UsageError, match="no HTTPS location"):
        pride.download_files([aspera_only], tmp_path)


def test_download_files_refuses_a_mixed_project_selection(tmp_path):
    a = pride.PrideFile._from_wire({"file_name": "a", "https_url": "https://x/a"}, "PXD000001")
    b = pride.PrideFile._from_wire({"file_name": "b", "https_url": "https://x/b"}, "PXD000002")
    with pytest.raises(pymzlib.UsageError, match="one project"):
        pride.download_files([a, b], tmp_path)


def test_download_files_rejects_things_that_are_not_pride_files(tmp_path):
    with pytest.raises(pymzlib.UsageError, match="PrideFile"):
        pride.download_files(["a.raw"], tmp_path)


# The live canaries now live in test_pride_live.py, where each one routes through
# external_service() so a PRIDE outage skips with an explanatory message instead of
# failing. Keeping them here, bare, made a red build ambiguous.