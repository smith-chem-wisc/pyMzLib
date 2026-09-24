"""Tests for the PRIDE surface.

Two tiers, deliberately separated:

* **Offline** — the bridge is replaced by a recorded payload, so these run anywhere, in
  milliseconds, with no network. They test the Python layer's own behavior.
* **Network** (``-m network``) — real calls to the live PRIDE Archive through the real
  bridge. These are the ones that would catch mzLib or PRIDE changing under us.

Run the fast set with ``pytest -m "not network"``.
"""

from __future__ import annotations

import copy
import json
from datetime import date, datetime
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


# --------------------------------------------------------------- offline: ftp-files

# The complete FTP listing the bridge returns for `pride ftp-files` (mzLib #1121): three top-level
# files — one of which the REST manifest omits — plus a file nested in a subdirectory.
FTP_LISTING = {
    "accession": "PXD000001",
    "file_count": 4,
    "approximate_total_size_bytes": 1638 + 210 * 1024 * 1024 + 429 * 1024 * 1024 + 864256,
    "files": [
        {"relative_path": "README.txt", "file_name": "README.txt",
         "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/README.txt",
         "approximate_size_bytes": 1638},
        {"relative_path": "run1.raw", "file_name": "run1.raw",
         "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/run1.raw",
         "approximate_size_bytes": 210 * 1024 * 1024},
        {"relative_path": "hidden_from_rest.mzML", "file_name": "hidden_from_rest.mzML",
         "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/hidden_from_rest.mzML",
         "approximate_size_bytes": 429 * 1024 * 1024},
        {"relative_path": "generated/summary.mztab", "file_name": "summary.mztab",
         "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/generated/summary.mztab",
         "approximate_size_bytes": 864256},
    ],
}


@pytest.fixture()
def recorded_ftp_listing(monkeypatch):
    """Serve a recorded PXD000001 FTP listing instead of calling the bridge.

    Hands each test a deep copy, the way ``recorded_manifest`` re-parses its JSON per use, so no
    test can leak a mutation into another through the shared module-level dict.
    """
    monkeypatch.setattr(_bridge, "invoke", lambda *args, timeout=None: copy.deepcopy(FTP_LISTING))
    return FTP_LISTING


def test_list_ftp_files_parses_every_file(recorded_ftp_listing):
    files = pride.list_ftp_files("PXD000001")
    # Pin len(files) to the payload's declared file_count, so a bridge/fixture drift in either is caught.
    assert len(files) == recorded_ftp_listing["file_count"] == 4
    assert all(isinstance(f, pride.PrideFtpFile) for f in files)


def test_ftp_file_fields_are_typed(recorded_ftp_listing):
    first = pride.list_ftp_files("PXD000001")[0]
    assert first.relative_path == "README.txt"
    assert isinstance(first.approximate_size_bytes, int)
    assert first.url.startswith("https://")
    assert first.project_accession == "PXD000001"


def test_ftp_nested_file_keeps_subdirectory_path_but_bare_leaf_name(recorded_ftp_listing):
    nested = next(f for f in pride.list_ftp_files("PXD000001") if "/" in f.relative_path)
    assert nested.relative_path == "generated/summary.mztab"
    assert nested.file_name == "summary.mztab"
    assert nested.extension == ".mztab"


def test_ftp_derived_properties_and_as_dict(recorded_ftp_listing):
    run = next(f for f in pride.list_ftp_files("PXD000001") if f.file_name == "run1.raw")
    assert run.extension == ".raw"
    assert run.approximate_size_mb == pytest.approx(run.approximate_size_bytes / 1_000_000)
    # as_dict must carry the computed properties, which asdict()/vars() skip.
    record = run.as_dict()
    assert record["approximate_size_mb"] == run.approximate_size_mb
    assert record["extension"] == ".raw"


def test_approximate_total_size_sums_over_the_complete_listing(recorded_ftp_listing):
    files = pride.list_ftp_files("PXD000001")
    assert pride.approximate_total_size_bytes(files) == sum(f.approximate_size_bytes for f in files)


def test_list_ftp_files_invokes_the_ftp_files_verb(monkeypatch):
    # The verb name is the contract with the bridge; a typo would only surface against the live
    # service, so pin it offline.
    seen: dict = {}

    def capturing_invoke(*args, timeout=None):
        seen["args"] = args
        seen["timeout"] = timeout
        return FTP_LISTING

    monkeypatch.setattr(_bridge, "invoke", capturing_invoke)
    pride.list_ftp_files("pxd000001", timeout=42)

    args = seen["args"]
    assert args[:2] == ("pride", "ftp-files")
    # Positional (flag/value adjacency), not mere membership: an accession under the wrong flag,
    # or a valueless --accession, must fail. Matches test_download_passes_accession_and_destination.
    assert args[args.index("--accession") + 1] == "PXD000001"  # normalised, upper-cased
    assert seen["timeout"] == 42


@pytest.mark.parametrize("accession", ["", "   ", "not-an-accession", 12345])
def test_list_ftp_files_rejects_a_bad_accession_before_any_work(accession, monkeypatch):
    # Install a raising sentinel so a validation regression fails loudly here instead of falling
    # through to the real bridge subprocess — an earlier mock-less test downloaded 480 MB of
    # PXD000001 into the source tree. This test's whole point is "we never reach the bridge".
    def must_not_run(*a, **k):
        raise AssertionError("validation was bypassed — list_ftp_files reached the bridge")

    monkeypatch.setattr(_bridge, "invoke", must_not_run)
    with pytest.raises(_bridge.UsageError):
        pride.list_ftp_files(accession)


def test_list_ftp_files_unknown_project_raises_project_not_found(monkeypatch):
    # mzLib resolves the project before walking, so an unknown accession comes back as an
    # MzLibException. list_ftp_files re-maps that to ProjectNotFoundError — the same "no such
    # project" signal list_files raises — so a caller catches one type across both functions.
    def failing_invoke(*args, timeout=None):
        raise _bridge.BridgeError(
            "MzLibException", "PRIDE Archive has no project with accession 'PXD999999'."
        )

    monkeypatch.setattr(_bridge, "invoke", failing_invoke)
    with pytest.raises(pride.ProjectNotFoundError):
        pride.list_ftp_files("PXD999999")  # valid form, so it reaches the bridge, then fails there


def test_list_ftp_files_no_publication_date_raises_project_not_found(monkeypatch):
    def failing_invoke(*args, timeout=None):
        raise _bridge.BridgeError(
            "MzLibException",
            "PRIDE project 'PXD000001' has no publication date, so its FTP directory cannot be located.",
        )

    monkeypatch.setattr(_bridge, "invoke", failing_invoke)
    with pytest.raises(pride.ProjectNotFoundError):
        pride.list_ftp_files("PXD000001")


def test_list_ftp_files_cyclic_listing_is_not_project_not_found(monkeypatch):
    # Also an MzLibException, but the project exists: its listing is broken. Reporting it as
    # "no such project, check for a typo" would send the caller looking for the wrong mistake.
    def failing_invoke(*args, timeout=None):
        raise _bridge.BridgeError(
            "MzLibException",
            "PRIDE FTP directory nesting exceeded 64 levels at 'https://x/'; the listing may be cyclic.",
        )

    monkeypatch.setattr(_bridge, "invoke", failing_invoke)
    with pytest.raises(_bridge.BridgeError) as caught:
        pride.list_ftp_files("PXD000001")
    assert not isinstance(caught.value, pride.ProjectNotFoundError)


def test_list_ftp_files_empty_listing_raises_rather_than_returning_empty(monkeypatch):
    # The project resolved but the directory parsed to nothing (e.g. PRIDE changed its autoindex).
    # That must not come back as [] — the "0 files, done" trap list_files also defends against.
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"files": []})
    with pytest.raises(pride.ProjectNotFoundError):
        pride.list_ftp_files("PXD000001")


def test_list_ftp_files_transport_failure_stays_a_bridge_error(monkeypatch):
    # Only "not found" (MzLibException) is re-mapped; a genuine transport failure keeps its
    # BridgeError so a caller can still tell an outage/HTTP error from a missing project.
    def failing_invoke(*args, timeout=None):
        raise _bridge.BridgeError("HttpRequestException", "PRIDE listing failed with status 500.")

    monkeypatch.setattr(_bridge, "invoke", failing_invoke)
    with pytest.raises(_bridge.BridgeError) as caught:
        pride.list_ftp_files("PXD000001")
    assert not isinstance(caught.value, pride.ProjectNotFoundError)


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
def test_accession_case_and_whitespace_are_normalised(recorded_manifest, monkeypatch, accession):
    """PRIDE's API is case-sensitive on the accession while our category matching is not. A
    lowercase accession used to return [] — indistinguishable from an empty project.

    Asserting only that the call succeeds would pass even if normalisation were removed, because
    the fixture ignores its arguments. The point is what reaches the bridge, so that is what is
    checked."""
    seen = {}

    def capturing_invoke(*args, timeout=None):
        seen["args"] = list(args)
        return recorded_manifest

    monkeypatch.setattr(_bridge, "invoke", capturing_invoke)
    pride.list_files(accession)

    assert seen["args"][seen["args"].index("--accession") + 1] == "PXD000001"


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


def test_download_raises_when_a_filter_matches_nothing(monkeypatch, tmp_path):
    """The doctrine has to hold in all three functions or it is not a doctrine. list_files and
    download_files both refuse to report success on nothing; download used to return [] happily,
    which in a batch script is zero files and a green exit code."""
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"paths": []})
    with pytest.raises(pymzlib.UsageError, match="matched"):
        pride.download("PXD000001", tmp_path, category="NOSUCHCATEGORY")


def test_download_without_a_filter_may_legitimately_write_nothing(monkeypatch, tmp_path):
    """The counterpart: no filter means 'everything', and only then is an empty result the
    repository's answer rather than a mistake in the call."""
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"paths": []})
    assert pride.download("PXD000001", tmp_path) == []


def test_a_path_of_empty_string_is_the_current_directory_and_cannot_be_distinguished(captured_args):
    """`Path("")` and `Path(".")` are the *same object* — pathlib collapses the empty string at
    construction, so by the time a guard sees it there is nothing left to detect.

    This is a real limitation and it is documented rather than papered over: the common mistake,
    `dest = config.get("outdir", "")`, passes a **string** and is caught. Someone who writes
    `Path(cfg.get("outdir", ""))` gets the current directory, exactly as if they had written
    `Path(".")` — which is a legitimate destination we must not refuse.

    (An earlier version of this test asserted a raise and had no mock, so it ran the real bridge
    and downloaded 480 MB of PXD000001 into the source tree — demonstrating the hazard rather
    more vividly than intended.)
    """
    assert Path("") == Path(".")
    pride.download("PXD000001", Path(""))
    assert captured_args["args"][captured_args["args"].index("--dest") + 1] == "."


def test_an_extension_list_that_names_nothing_is_refused(tmp_path):
    """It would otherwise normalise to [], omit --ext entirely, and download the whole project —
    the same fail-open the bridge was just fixed for, one layer up."""
    with pytest.raises(pymzlib.UsageError, match="names no extensions"):
        pride.download("PXD000001", tmp_path, extensions=["", "   "])


def test_a_file_name_containing_a_newline_is_refused(tmp_path):
    """The stdin framing is newline-delimited, and a POSIX file name may legally contain one.
    Splitting it silently would select the wrong files."""
    odd = pride.PrideFile._from_wire(
        {"file_name": "two\nlines.raw", "https_url": "https://x/a"}, "PXD000001"
    )
    with pytest.raises(pymzlib.UsageError, match="line break"):
        pride.download_files([odd], tmp_path)


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


# ---- search ----------------------------------------------------------------------------------
#
# The fixture is a real recorded `pride search` payload, chosen because it carries the three things
# that are easy to get wrong and impossible to invent by hand: bare calendar dates, a `highlights`
# map with PRIDE's own field names as keys, and a genuinely blank string inside `keywords`.

SEARCH_FIXTURE = Path(__file__).parent / "fixtures" / "pride_search_plasmodium.json"


@pytest.fixture()
def recorded_search(monkeypatch):
    payload = json.loads(SEARCH_FIXTURE.read_text(encoding="utf-8"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: payload)
    return payload


@pytest.fixture()
def hits(recorded_search):
    return pride.search("plasmodium falciparum schizont")


def test_search_returns_one_result_per_wire_entry(hits, recorded_search):
    assert len(hits) == recorded_search["result_count"] == 6
    assert hits[0].accession == "PXD070842"


def test_dates_are_calendar_dates_not_timestamps(hits):
    """The single most important thing this projection gets right.

    PRIDE's search endpoint sends "2025-11-17" — no time, no offset — which is why mzLib types
    these `DateTime` where PrideFile's are `DateTimeOffset`. `datetime.fromisoformat` would hand
    back a naive midnight, a time PRIDE never reported, and `pride._parse_timestamp` (used by
    PrideFile) does exactly that. A `date` says what was sent and nothing more.
    """
    assert isinstance(hits[0].submission_date, date)
    assert not isinstance(hits[0].submission_date, datetime)
    assert hits[0].submission_date == date(2025, 11, 17)


def test_an_absent_date_is_none_not_year_one():
    """`default(DateTime)` is 0001-01-01, and crossing that as a real date would be a trap.

    The DTO's date members are non-nullable, so a field PRIDE omitted arrives as year 1 rather than
    null. A caller doing arithmetic on that gets a two-thousand-year interval and no warning.
    """
    assert pride._parse_date(None) is None
    assert pride._parse_date("") is None
    assert pride._parse_date("not a date") is None


def test_highlights_keys_survive_unchanged(hits):
    """They are PRIDE's field names, not ours, so the snake_case policy must not have touched them.

    `JsonNamingPolicy.SnakeCaseLower` renames properties; `DictionaryKeyPolicy` is unset, so these
    pass through. Pinned because a later options change could quietly start rewriting them, and the
    keys are the only thing that says WHY a project matched.
    """
    assert hits[0].matched_fields == ["references", "title"]
    assert set(hits[0].highlights) <= {"references", "title", "project_description", "keywords",
                                       "organisms", "instruments", "sample_attributes",
                                       "data_processing_protocol", "sample_processing_protocol",
                                       "accession", "diseases", "organismsPart", "submitters"}


def test_blank_keywords_are_passed_through_not_filtered(hits):
    """PRIDE ships empty and whitespace-only keywords on ~9% of hits, and this fixture has one.

    Filtering them here would be a repair on the island: pyMzLib would then disagree with mzLib,
    and with the Rust and R bindings, about what a project's keywords are. The docstring tells the
    caller to filter; the projection does not do it for them.
    """
    assert any(not k.strip() for h in hits for k in h.keywords)


def test_a_zero_count_is_reported_as_zero_not_none(hits):
    """Zero means "not reported", but mapping it to None would assert an absence mzLib does not."""
    assert all(isinstance(h.download_count, int) for h in hits)
    assert all(isinstance(h.percentile, int) for h in hits)


def test_cv_fields_arrive_flattened_to_strings(hits):
    """The reason this is a separate type from a project's metadata, pinned as a fact about the wire."""
    assert hits[0].organisms == ["Homo sapiens (human)", "Plasmodium falciparum (isolate 3d7)"]
    assert all(isinstance(x, str) for h in hits for x in h.instruments)


def test_as_dict_carries_the_computed_property(hits):
    """`dataclasses.asdict` skips properties — the same trap PrideFile.as_dict already documents."""
    record = hits[0].as_dict()
    assert record["matched_fields"] == hits[0].matched_fields
    assert record["accession"] == "PXD070842"


def test_search_sends_the_keyword_and_page_size(monkeypatch):
    seen: dict = {}
    monkeypatch.setattr(
        _bridge, "invoke", lambda *a, **k: seen.update(args=a) or {"results": []}
    )

    pride.search("phospho", page_size=25)

    assert seen["args"] == ("pride", "search", "--keyword", "phospho", "--page-size", "25")


def test_no_hits_is_an_empty_list_not_an_error(monkeypatch):
    """Unlike list_files, which raises: there is no accession here that could have been a typo.

    PRIDE reports no matches as an empty result rather than a 404, and for a keyword that is the
    honest answer - "nothing matched" is a real finding, not a mistake to protect the caller from.
    """
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: {"keyword": "x", "result_count": 0,
                                                            "results": []})

    assert pride.search("qwertyuiop") == []


@pytest.mark.parametrize("keyword", ["", "   ", None, 7])
def test_search_requires_a_keyword(keyword):
    with pytest.raises(pymzlib.UsageError, match="search keyword is required"):
        pride.search(keyword)


def test_an_over_long_keyword_is_refused_before_the_subprocess():
    """PRIDE answers a very long keyword with HTTP 500, which reads as an outage rather than a bug.

    ExternalServiceTestHelper and the bridge's own classifier both treat a 5xx as "the service is
    down" and SKIP, so without this a caller's mistake would be reported as EBI having a bad day.
    Refused here, before any process is spawned, and again in the bridge.
    """
    with pytest.raises(pymzlib.UsageError, match="at most 1000 characters"):
        pride.search("x" * (pride.MAX_KEYWORD_LENGTH + 1))


def test_a_keyword_that_looks_like_an_option_is_refused():
    """The bridge parses `--a --b` as two flags, so a leading dash would discard the option."""
    with pytest.raises(pymzlib.UsageError, match="may not begin with"):
        pride.search("--page-size")


@pytest.mark.parametrize("page_size", [0, -1, 2.5, "10", True])
def test_search_rejects_a_bad_page_size(page_size):
    with pytest.raises(pymzlib.UsageError):
        pride.search("phospho", page_size=page_size)
