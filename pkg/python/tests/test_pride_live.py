"""Live canary against the real PRIDE Archive.

Everything else in this suite runs against a recorded fixture, which means it would keep passing
indefinitely even if EBI changed its API tomorrow. These are the tests that would notice.

Every one routes through :func:`external_service`, so a PRIDE outage **skips** with an explanatory
message while a genuine contract break still **fails**. A red build here should always mean "we
broke something", never "EBI is having a bad morning" — the convention is borrowed from mzLib,
where these carry ``[Category("ExternalService")]`` and run in a dedicated CI job.

Run with ``pytest -m network``.
"""

from __future__ import annotations

import pytest

import pymzlib
from conftest import external_service

pytestmark = pytest.mark.network

#: A small, long-stable public project. Changing this changes what the canary proves.
CANARY_ACCESSION = "PXD000001"


def test_the_bridge_and_the_package_agree_on_the_protocol():
    """Not a network test in spirit, but it needs the real executable rather than a stub."""
    with external_service():
        info = pymzlib.bridge_version()
    assert info["protocol"] == pymzlib._bridge.PROTOCOL_VERSION


def test_the_api_still_answers_and_the_manifest_still_parses():
    with external_service():
        files = pymzlib.pride.list_files(CANARY_ACCESSION)

    assert files, (
        "PRIDE answered but reported no files for a project known to have them — "
        "the response shape has probably changed."
    )


def test_the_fields_the_python_layer_reads_are_still_populated():
    """Each assertion here corresponds to an attribute users actually touch."""
    with external_service():
        files = pymzlib.pride.list_files(CANARY_ACCESSION)

    first = files[0]
    assert first.file_name
    assert first.category
    assert first.file_size_bytes > 0
    assert first.submission_date is not None
    assert first.submission_date.tzinfo is not None, "timestamps must stay timezone-aware"


def test_at_least_one_file_is_still_reachable_over_https():
    """The FTP-to-HTTPS upgrade is an assumption about EBI's publishing, not a guarantee."""
    with external_service():
        files = pymzlib.pride.list_files(CANARY_ACCESSION)

    assert any(f.downloadable for f in files), (
        "No file exposed an HTTPS location — the FTP-to-HTTPS upgrade assumption may no longer hold, "
        "which would break every download."
    )


def test_an_unknown_accession_still_returns_empty_rather_than_erroring():
    """PRIDE's own behavior, which the Python layer deliberately preserves."""
    with external_service():
        assert pymzlib.pride.list_files("PXD999999999") == []


@pytest.mark.slow
def test_a_real_download_still_works_end_to_end(tmp_path):
    with external_service():
        written = pymzlib.pride.download(CANARY_ACCESSION, tmp_path, extensions=[".fasta"])

    assert len(written) == 1
    assert written[0].is_file()
    assert written[0].stat().st_size > 0
    assert not list(tmp_path.glob("*.partial")), "no partial file may survive a successful download"
