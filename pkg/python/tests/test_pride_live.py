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


def test_an_unknown_accession_raises_rather_than_reporting_an_empty_project():
    """PRIDE answers an unknown accession with an empty result rather than a 404. pyMzLib no
    longer passes that through: an empty list is indistinguishable from a project that genuinely
    has no matching files, so a typo used to produce '0 files, done' and a script that carried on."""
    with external_service():
        with pytest.raises(pymzlib.pride.ProjectNotFoundError):
            pymzlib.pride.list_files("PXD999999999")


def test_a_real_selection_can_be_downloaded_directly(tmp_path):
    """The end-to-end shape the API is for: list, filter in Python, download exactly that."""
    with external_service():
        files = pymzlib.pride.list_files(CANARY_ACCESSION)
        chosen = [f for f in files if f.downloadable and f.size_mb < 2]
        assert chosen, "expected at least one small file in the canary project"
        written = pymzlib.pride.download_files(chosen[:1], tmp_path)

    assert len(written) == 1
    assert written[0].is_file() and written[0].stat().st_size > 0


@pytest.mark.slow
def test_a_real_download_still_works_end_to_end(tmp_path):
    with external_service():
        written = pymzlib.pride.download(CANARY_ACCESSION, tmp_path, extensions=[".fasta"])

    assert len(written) == 1
    assert written[0].is_file()
    assert written[0].stat().st_size > 0
    assert not list(tmp_path.glob("*.partial")), "no partial file may survive a successful download"


# ------------------------------------------------------------------- ftp-files


def test_the_ftp_listing_is_more_complete_than_the_rest_manifest():
    """The whole reason list_ftp_files exists: the FTP tree holds files the REST manifest omits.

    Asserted as a comparison (FTP > REST), not a hard count, so it survives PRIDE re-curating the
    project — while still catching the two ways this could break: the FTP walk finding nothing, or
    the REST manifest quietly becoming complete (in which case the extra verb has lost its purpose).
    """
    with external_service():
        rest = pymzlib.pride.list_files(CANARY_ACCESSION)
        ftp = pymzlib.pride.list_ftp_files(CANARY_ACCESSION)

    assert ftp, "the FTP walk returned nothing for a project known to have files"
    assert len(ftp) > len(rest), (
        f"the FTP listing ({len(ftp)}) should exceed the REST manifest ({len(rest)}); if PRIDE's "
        "REST API has become complete, list_ftp_files no longer earns its place."
    )


def test_the_ftp_fields_the_python_layer_reads_are_still_populated():
    with external_service():
        ftp = pymzlib.pride.list_ftp_files(CANARY_ACCESSION)

    first = ftp[0]
    assert first.relative_path
    assert first.file_name
    assert first.url.startswith("https://"), "the download URL must be the HTTPS location"
    assert first.approximate_size_bytes > 0


def test_the_ftp_total_size_covers_more_than_the_rest_manifest_total():
    """The complete-but-approximate total should exceed the incomplete REST total — the two size
    numbers fail in opposite directions, and this pins that the FTP one is the larger, fuller one."""
    with external_service():
        rest = pymzlib.pride.list_files(CANARY_ACCESSION)
        ftp = pymzlib.pride.list_ftp_files(CANARY_ACCESSION)

    assert pymzlib.pride.approximate_total_size_bytes(ftp) > pymzlib.pride.total_size_bytes(rest)
