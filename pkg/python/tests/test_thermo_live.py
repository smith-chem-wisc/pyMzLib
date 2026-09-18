"""Read a real Thermo ``.raw`` file through the real, packaged bridge.

Every other readers test is offline against recorded payloads, and the bridge's C# suite proves
each file type is reachable. Neither could see the failure this guards: the C# suite runs mzLib
in-process from loose assemblies, but the wheel ships mzLib inside a single-file executable, and
Thermo's RawFileReader cannot run from inside one. 0.1.0 shipped with every ``.raw`` read failing
("Method invocation failed on Method[ThreadedFileFactory]") while every test passed. Only a read
through the packaged executable, on each platform CI tests the wheel on, can catch that.

The fixture is ``ScanDescriptionTestData.raw`` from mzLib's own test data: 54 scans, 0.5 MB.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from pymzlib import _bridge, readers

THERMO_FIXTURE = Path(__file__).parent / "fixtures" / "thermo_scan_descriptions.raw"


@pytest.fixture()
def built_bridge():
    # Same rule as the live check in test_readers.py: skip in a source checkout with no built
    # bridge, but fail under CI, where the bridge is always built and a skip would hide the one
    # test that reads a vendor file through the packaged executable.
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so .raw reading went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def test_a_thermo_raw_file_reads_through_the_packaged_bridge(built_bridge):
    scans = readers.read_spectra(THERMO_FIXTURE, limit=3, peaks=True)

    assert scans.file_type == "ThermoRaw"
    assert scans.scan_count == 54
    assert len(scans.records) == 3
    # Peaks prove the reader got past opening the file and into the scan data itself. (The third
    # scan in this file is genuinely empty, so not every scan has peaks.)
    assert all(len(record["mz"]) == record["peak_count"] for record in scans.records)
    assert scans.records[0]["peak_count"] == 37
