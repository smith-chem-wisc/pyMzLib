"""``pymzlib.library`` through the real, built bridge.

``test_library.py`` runs against recorded payloads. These run the same calls against the bridge
itself, which is the only way to catch the recordings drifting from what the bridge now says, a
batch whose table changes with ``threads``, or a refusal that does not reach Python as the error
type the docstrings promise.

Gated like ``test_readers_bulk_live.py``: skipped in a source checkout with no built bridge, failed
under CI, where the bridge is always built. Offline: the libraries are the two committed under
``fixtures/library/``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, library

HERE = Path(__file__).parent
FIXTURES = HERE / "fixtures"
SYNTHETIC = FIXTURES / "library" / "synthetic_library.msl"
PROSIT = FIXTURES / "library" / "myPrositLib.msl"


@pytest.fixture()
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so library read went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")
    _bridge._VERBS_SEEN.clear()


def without_paths(data: dict) -> dict:
    """A payload with its machine-specific absolute paths removed, for comparing recordings."""
    trimmed = {key: value for key, value in data.items() if key not in ("path", "files")}
    if isinstance(trimmed.get("columns"), dict) and "source_path" in trimmed["columns"]:
        trimmed["columns"] = {k: v for k, v in trimmed["columns"].items() if k != "source_path"}
    if "files" in data:
        trimmed["files"] = [{k: v for k, v in f.items() if k != "path"} for f in data["files"]]
    return trimmed


def recorded(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def test_the_bridge_lists_library_read(built_bridge):
    assert "library read" in pymzlib.bridge_version()["verbs"]


@pytest.mark.parametrize(
    "fixture, args",
    [
        ("library_read_prosit.json", ["--path", str(PROSIT)]),
        ("library_read_synthetic.json", ["--path", str(SYNTHETIC)]),
        ("library_read_fragments.json", ["--path", str(SYNTHETIC), "--fragments", "--limit", "1"]),
        ("library_read_index_only.json", ["--path", str(SYNTHETIC), "--index-only-threshold", "0", "--limit", "2"]),
    ],
)
def test_the_recordings_are_what_the_bridge_says_now(built_bridge, fixture, args):
    live = _bridge.invoke("library", "read", *args)

    assert without_paths(live) == without_paths(recorded(fixture))


def test_the_batch_recording_is_what_the_bridge_says_now(built_bridge):
    live = _bridge.invoke("library", "read", "--paths-stdin", stdin=f"{SYNTHETIC}\n{PROSIT}\n")

    assert without_paths(live) == without_paths(recorded("library_read_many.json"))


def test_a_batch_is_identical_at_one_and_four_threads(built_bridge):
    one = library.read_many([SYNTHETIC, PROSIT], fragments=True, threads=1)
    four = library.read_many([SYNTHETIC, PROSIT], fragments=True, threads=4)

    assert one == four
    assert one.row_count > one.returned_count == 11


def test_index_only_and_full_load_give_the_same_fragment_rows(built_bridge):
    full = library.read(SYNTHETIC, fragments=True)
    index_only = library.read(SYNTHETIC, fragments=True, index_only_threshold=0)

    assert (full.is_index_only, index_only.is_index_only) == (False, True)
    assert index_only.columns == full.columns


def test_out_writes_the_table(built_bridge, tmp_path):
    target = tmp_path / "nested" / "ions.tsv"
    written = library.read(SYNTHETIC, fragments=True, out=target)

    lines = target.read_text(encoding="utf-8").splitlines()
    assert written.columns is None
    assert lines[0].split("\t") == written.column_names
    assert len(lines) - 1 == written.output.row_count == 14


def test_a_text_library_is_a_usage_error(built_bridge, tmp_path):
    text = tmp_path / "library.msp"
    text.write_text("Name: PEPTIDE/2\n", encoding="utf-8")

    with pytest.raises(pymzlib.UsageError, match=r"\.msl"):
        library.read(text)


def test_a_damaged_library_is_mzlibs_correctness_error(built_bridge, tmp_path):
    damaged = tmp_path / "damaged.msl"
    data = bytearray(SYNTHETIC.read_bytes())
    data[200] ^= 0xFF                                  # inside the data the CRC-32 covers
    damaged.write_bytes(bytes(data))

    with pytest.raises(pymzlib.BridgeError) as caught:
        library.read(damaged)
    assert not isinstance(caught.value, pymzlib.UsageError)
    assert caught.value.error_type in ("InvalidDataException", "FormatException")


def test_a_skipped_library_is_reported_and_the_rest_read(built_bridge, tmp_path):
    batch = library.read_many([SYNTHETIC, tmp_path / "missing.msl"], on_error="skip")

    assert (batch.read_count, batch.failed_count) == (1, 1)
    (failed,) = batch.failed_files
    assert failed.error.kind == "usage"
    assert failed.record_count is None and failed.target_count is None
