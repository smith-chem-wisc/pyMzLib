"""The readers batch through the real, built bridge.

Everything in ``test_readers_bulk.py`` runs against recorded payloads. These run the same calls
against the bridge itself, which is the only way to catch the recordings drifting from what the
bridge now says, a ``_many`` call whose table changes with ``threads``, or a stale bridge that no
longer lists the verbs pyMzLib calls.

Gated like ``test_thermo_live.py``: skipped in a source checkout with no built bridge, failed under
CI, where the bridge is always built. The quantification checks also need mzLib's own 1.0.592 test
files from the pinned ``code/mzLib`` worktree, and skip without them.
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, readers

HERE = Path(__file__).parent
FIXTURES = HERE / "fixtures"
THERMO = FIXTURES / "thermo_scan_descriptions.raw"
MZLIB_TEST = HERE.parents[2] / "code" / "mzLib" / "mzLib" / "Test"
PROTEIN_GROUPS = MZLIB_TEST / "FileReadingTests" / "ExternalFileTypes" / "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"
PEPTIDES = MZLIB_TEST / "FileReadingTests" / "ExternalFileTypes" / "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"


@pytest.fixture()
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so the readers batch went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")
    _bridge._VERBS_SEEN.clear()


def mzlib_file(path: Path) -> str:
    if not path.is_file():
        pytest.skip(f"mzLib test file not in this checkout: {path}")
    return str(path)


def without_paths(data: dict) -> dict:
    """A payload with its machine-specific absolute paths removed, for comparing recordings."""
    return {key: value for key, value in data.items() if key != "path"}


def test_the_bridge_lists_the_verbs_this_package_calls(built_bridge):
    verbs = pymzlib.bridge_version()["verbs"]

    for verb in ("readers read-protein-groups", "readers read-quantified-peptides", "readers read-occupancy"):
        assert verb in verbs


@pytest.mark.parametrize("fixture, verb, path, extra", [
    ("readers_protein_groups.json", "read-protein-groups", PROTEIN_GROUPS, ["--limit", "1"]),
    ("readers_quantified_peptides.json", "read-quantified-peptides", PEPTIDES, ["--limit", "1"]),
    ("readers_occupancy.json", "read-occupancy", PROTEIN_GROUPS, []),
])
def test_the_recorded_quant_fixtures_still_match_the_live_bridge(built_bridge, fixture, verb, path, extra):
    live = _bridge.invoke("readers", verb, "--path", mzlib_file(path), *extra, timeout=120)
    recorded = json.loads((FIXTURES / fixture).read_text(encoding="utf-8"))

    assert without_paths(live) == without_paths(recorded), (
        f"{fixture} no longer matches the bridge: re-record it from the live bridge")


def test_a_list_reads_the_same_at_one_and_four_threads(built_bridge, tmp_path):
    # Two copies under different names: a batch refuses a repeated path, not a repeated file.
    copy = tmp_path / "second.raw"
    shutil.copyfile(THERMO, copy)
    paths = [THERMO, copy]

    one = readers.read_spectra_many(paths, threads=1, timeout=300)
    four = readers.read_spectra_many(paths, threads=4, timeout=300)

    assert one.columns == four.columns
    assert [f.record_count for f in one.files] == [54, 54]
    assert one.columns["source_index"].count(1) == 54
    assert one.files[0].source.instrument_model == "Orbitrap Fusion Lumos"
    assert one.files[0].source.acquisition_start_time_is_utc is False


def test_a_list_skips_what_it_cannot_read_when_asked(built_bridge, tmp_path):
    missing = tmp_path / "gone.raw"

    batch = readers.read_spectra_many([missing, THERMO], on_error="skip", timeout=300)

    assert (batch.read_count, batch.failed_count) == (1, 1)
    assert batch.failed_files[0].error.kind == "usage"
    assert set(batch.columns["source_index"]) == {1}


def test_a_list_stops_on_what_it_cannot_read_by_default(built_bridge, tmp_path):
    with pytest.raises(pymzlib.UsageError, match=r"^Input 1 \("):
        readers.read_spectra_many([THERMO, tmp_path / "gone.raw"], timeout=300)


def test_identify_many_answers_in_order(built_bridge, tmp_path):
    batch = readers.identify_many([THERMO, tmp_path / "gone.mzML"], on_error="skip")

    assert [f.file_type for f in batch.files] == ["ThermoRaw", ""]
    assert batch.files[1].error.kind == "usage"
