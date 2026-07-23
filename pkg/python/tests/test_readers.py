"""Tests for result-file identification.

Offline throughout. The Python layer's job here is small and worth pinning exactly: assemble the
argument list, parse the wire payload into typed objects, and answer "can this file be quantified"
without the caller having to know that the answer is a string in a list. The dispatch itself is
mzLib's and is covered by the bridge's C# suite.

The ``readers formats`` fixture is recorded from the real bridge rather than hand-written, so it
reflects what mzLib actually dispatches.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, readers

FORMATS_FIXTURE = Path(__file__).parent / "fixtures" / "readers_formats.json"


@pytest.fixture()
def recorded_formats(monkeypatch):
    """Serve the recorded ``readers formats`` payload instead of running the bridge."""
    payload = json.loads(FORMATS_FIXTURE.read_text(encoding="utf-8"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: payload)
    return payload


@pytest.fixture()
def captured(monkeypatch):
    """Capture the args a call would send, without running anything."""
    seen: dict = {}

    def fake_invoke(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return {"path": "C:/x/AllPSMs.psmtsv", "file_type": "psmtsv", "extension": ".psmtsv",
                "reader": "PsmFromTsvFile", "views": ["quantifiable"]}

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# ---- identify ------------------------------------------------------------------------------


def test_identify_parses_the_payload_into_a_typed_result(captured):
    info = readers.identify("AllPSMs.psmtsv")

    assert isinstance(info, readers.FileInfo)
    assert info.file_type == "psmtsv"
    assert info.reader == "PsmFromTsvFile"
    assert info.views == ["quantifiable"]


def test_identify_sends_the_expected_verb_and_option(captured):
    readers.identify("  AllPSMs.psmtsv  ")

    assert captured["args"] == ("readers", "identify", "--path", "AllPSMs.psmtsv")


def test_is_quantifiable_is_true_only_for_the_quantifiable_view(captured):
    assert readers.identify("AllPSMs.psmtsv").is_quantifiable is True


def test_a_file_with_no_uniform_view_is_not_quantifiable(monkeypatch):
    # TopPIC is read perfectly well by mzLib and has no cross-format view. A caller must be able to
    # find that out here rather than discover it as a failure inside flashlfq.quantify().
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "path": "C:/x/run_prsm.tsv", "file_type": "ToppicPrsm", "extension": "_prsm.tsv",
        "reader": "ToppicSearchResultFile", "views": [],
    })

    info = readers.identify("run_prsm.tsv")

    assert info.views == []
    assert info.is_quantifiable is False


def test_a_spectra_file_is_not_quantifiable(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {
        "path": "C:/x/run.mzML", "file_type": "MzML", "extension": ".mzML",
        "reader": "MsDataFileToResultFileAdapter", "views": ["spectra"],
    })

    assert readers.identify("run.mzML").is_quantifiable is False


def test_missing_wire_fields_do_not_crash_the_parser(monkeypatch):
    # The wire contract is versioned but a consumer should degrade rather than raise on a payload
    # that predates a field.
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"file_type": "psmtsv"})

    info = readers.identify("x.psmtsv")

    assert info.views == []
    assert info.is_quantifiable is False
    assert info.extension is None


@pytest.mark.parametrize("bad", ["", "   ", None, 42])
def test_a_blank_or_non_string_path_is_rejected_without_starting_a_process(bad):
    with pytest.raises(pymzlib.UsageError):
        readers.identify(bad)


# ---- formats -------------------------------------------------------------------------------


def test_formats_parses_every_entry(recorded_formats):
    formats = readers.formats()

    assert len(formats) == recorded_formats["format_count"]
    assert all(isinstance(f, readers.Format) for f in formats)
    assert all(f.file_type for f in formats)


def test_exactly_three_formats_offer_the_quantifiable_view(recorded_formats):
    # The headline fact about this module, pinned on the Python side too: mzLib reads 29 formats
    # and only these three can feed flashlfq.quantify(). Documented in the module docstring, so if
    # mzLib widens the set this test fails and the docs get corrected rather than quietly lying.
    quantifiable = [f.file_type for f in readers.formats() if f.is_quantifiable]

    assert sorted(quantifiable) == sorted(["psmtsv", "osmtsv", "MsFraggerPsm"])


def test_most_formats_have_no_uniform_view(recorded_formats):
    viewless = [f for f in readers.formats() if not f.views]

    assert len(viewless) > 10, "an empty view list is the common case, not an edge case"


def test_extensions_are_not_unique_across_formats(recorded_formats):
    # Documented in Format.extension and worth pinning: both Bruker types are ".d", so extension is
    # not a key. Anyone tempted to build a dict keyed on it should find this test instead of a bug.
    extensions = [f.extension for f in readers.formats()]

    assert len(extensions) != len(set(extensions))


def test_formats_sends_no_options(monkeypatch):
    seen = {}
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: seen.update(args=a) or {"formats": []})

    readers.formats()

    assert seen["args"] == ("readers", "formats")


# ---- module wiring -------------------------------------------------------------------------


def test_readers_is_exported_from_the_package():
    # A module missing from __all__ is invisible to `from pymzlib import *` and to the docs build.
    assert "readers" in pymzlib.__all__
    assert pymzlib.readers is readers
