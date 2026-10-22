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

import hashlib

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


def test_the_ftp_total_size_is_populated_and_consistent():
    """The complete-project size is a positive estimate that sums the listing.

    Deliberately *not* compared against the REST total: PRIDE reports decompressed sizes for `.gz`
    files, so the incomplete REST total is inflated by an unbounded amount and the two numbers "fail
    in opposite directions" (as the docs say) — an FTP-vs-REST inequality is not durable. What is
    durable: the total is positive and at least the largest single file it sums.
    """
    with external_service():
        ftp = pymzlib.pride.list_ftp_files(CANARY_ACCESSION)

    total = pymzlib.pride.approximate_total_size_bytes(ftp)
    assert total >= max(f.approximate_size_bytes for f in ftp) > 0


# ------------------------------------------------------------------- checksums


def test_the_checksum_list_still_parses_and_covers_the_deposited_files():
    """mzLib #1427. Asserted as relations to the FTP tree rather than counts, so it survives PRIDE
    re-curating the project: every row names a file at the project's root, and none names one of
    PRIDE's generated/ files."""
    with external_service():
        checksums = pymzlib.pride.file_checksums(CANARY_ACCESSION)
        ftp = pymzlib.pride.list_ftp_files(CANARY_ACCESSION)

    assert checksums, "PRIDE has always listed checksums for the canary; an empty list means the route changed"
    root = {f.relative_path for f in ftp if "/" not in f.relative_path}
    generated = {f.relative_path for f in ftp if f.relative_path.startswith("generated/")}
    assert set(checksums) <= root, f"rows outside the FTP root: {sorted(set(checksums) - root)}"
    assert generated, "the canary should still hold generated/ files, which have no checksum row"
    assert all(len(c.md5) == 32 and c.size_bytes > 0 for c in checksums.values())


def test_a_checksum_row_matches_the_bytes_a_download_delivers(tmp_path):
    """The point of the list: its size and MD5 are those of the file a download writes. A small
    file in both the REST manifest and the checksum list keeps this fast."""
    with external_service():
        checksums = pymzlib.pride.file_checksums(CANARY_ACCESSION)
        files = pymzlib.pride.list_files(CANARY_ACCESSION)
        chosen = min(
            (f for f in files if f.downloadable and f.file_name in checksums),
            key=lambda f: checksums[f.file_name].size_bytes,
        )
        (written,) = pymzlib.pride.download_files([chosen], tmp_path)

    expected = checksums[chosen.file_name]
    assert written.stat().st_size == expected.size_bytes
    assert hashlib.md5(written.read_bytes()).hexdigest() == expected.md5.lower()


def test_an_unknown_accession_raises_rather_than_reporting_no_checksums():
    """PRIDE answers an unknown accession with an empty list, as it does a project with no
    checksums. mzLib asks for the project to tell them apart, so a typo is not '0 files to check'."""
    with external_service():
        with pytest.raises(pymzlib.pride.ProjectNotFoundError):
            pymzlib.pride.file_checksums("PXD999999999")
