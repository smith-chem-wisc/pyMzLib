"""Tests for the SDRF-Proteomics projection.

Offline except for the two canaries at the end. The Python layer's job here is small and worth
pinning exactly: validate arguments before spending a subprocess, render the pool stdin, and turn
a row-major payload into objects that do not quietly lie about ragged rows, repeated column names,
or the difference between a reserved word and an absent column. The reading itself is mzLib's and
is covered by the bridge's C# suite.

The fixtures are recorded from the real bridge rather than hand-written, so they carry the shape
mzLib actually produces - including PXD059974's raggedness and PXD000070's eight repetitions of
``comment[modification parameters]``, both of which are the point rather than an accident. Only
the ``path``/``paths`` fields are rewritten, so the recording machine's directory layout does not
end up committed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, sdrf

FIXTURES = Path(__file__).parent / "fixtures"
READ_FIXTURE = FIXTURES / "sdrf_read_PXD000070.json"
RAGGED_FIXTURE = FIXTURES / "sdrf_read_ragged.json"
POOL_FIXTURE = FIXTURES / "sdrf_pool_two.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture()
def document(monkeypatch) -> sdrf.SdrfDocument:
    """PXD000070 as the bridge really returned it."""
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load(READ_FIXTURE))
    return sdrf.read("PXD000070.sdrf.tsv")


@pytest.fixture()
def ragged(monkeypatch) -> sdrf.SdrfDocument:
    """PXD059974: a 46-column header whose rows carry 42 cells."""
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load(RAGGED_FIXTURE))
    return sdrf.read("PXD059974.sdrf.tsv")


@pytest.fixture()
def pooled(monkeypatch) -> sdrf.PooledSdrf:
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: _load(POOL_FIXTURE))
    return sdrf.pool({"a.sdrf.tsv": "malaria", "b.sdrf.tsv": "colon"})


@pytest.fixture()
def captured(monkeypatch):
    """Capture what a call would send, without running anything."""
    seen: dict = {}

    def fake_invoke(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return _load(POOL_FIXTURE)

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# ---- the shape ------------------------------------------------------------------------------


def test_columns_is_a_list_because_sdrf_names_repeat(document):
    # The whole reason this module is row-major. PXD000070 carries
    # comment[modification parameters] eight times; a dict keyed by name would keep one.
    assert isinstance(document.columns, list)
    assert document.has_repeated_columns
    assert len(document.indexes_of("comment[modification parameters]")) == 8


def test_all_returns_every_occurrence_where_value_returns_the_first(document):
    every = document.all("comment[modification parameters]")[0]
    first = document.value("comment[modification parameters]")[0]

    assert len(every) == 8
    assert every[0] == first
    assert first.startswith("NT=Carbamidomethyl")


def test_a_cv_cell_crosses_intact_rather_than_being_split_on_its_semicolons(document):
    # The defect this module exists to close: readers.read_records() joins an SdrfRow's cells with
    # ";" and the SDRF key=value grammar is itself semicolon-delimited, so the result cannot be
    # split back apart. Here one cell is one string, semicolons and all.
    cell = document.value("comment[modification parameters]")[0]

    assert cell == "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"
    assert ";" in cell


def test_an_absent_column_is_none_and_a_reserved_word_is_itself(document):
    # These must not collapse. None means "this document has no such column"; "not applicable"
    # means "the experiment said there is none", which is a real answer a curator wrote.
    assert document.value("characteristics[nonesuch]") == [None] * document.returned_count
    assert document.value("characteristics[disease]")[0] == "not applicable"


def test_a_short_row_reports_none_rather_than_shifting_later_columns(ragged):
    # Raggedness is data: mzLib preserves it so the file round-trips. Padding here would invent
    # cells, and shifting would silently misattribute every value past the gap.
    assert ragged.ragged_row_count == 17
    short = next(i for i, row in enumerate(ragged.rows) if len(row) < len(ragged.columns))
    full = next(i for i, row in enumerate(ragged.rows) if len(row) == len(ragged.columns))

    # The same column, reached from a row that has it and one that does not. The short row must
    # say None rather than borrowing the value from an earlier position.
    last = ragged.columns[-1]
    assert ragged.value(last)[full] is not None
    assert ragged.value(last)[short] is None

    # And a position both rows DO reach still lines up, which is what "no shifting" means.
    assert ragged.value("characteristics[organism]")[short] is not None


def test_caveats_name_the_raggedness_rather_than_leaving_it_to_be_discovered(ragged):
    assert any("SHORT" in c for c in ragged.caveats)


def test_truncated_says_rows_were_left_behind(pooled):
    # A short answer and a complete one must never look alike. The pooled fixture was recorded
    # with limit=4 against a 24-row merge, so it is the windowed case.
    assert pooled.truncated
    assert pooled.returned_count == 4 < pooled.row_count == 24


def test_a_complete_read_is_not_marked_truncated(ragged):
    assert not ragged.truncated
    assert ragged.returned_count == ragged.row_count


def test_records_is_offered_but_flagged_lossy_when_names_repeat(document):
    # It is the convenient shape and the wrong one for this file; has_repeated_columns is how a
    # caller finds that out before trusting it.
    assert document.has_repeated_columns
    assert len(document.records[0]) < len(document.columns)


# ---- pool -----------------------------------------------------------------------------------


def test_pool_carries_provenance_and_the_labels_it_was_given(pooled):
    assert pooled.document_count == 2
    assert pooled.labels == ["malaria", "colon"]
    assert sdrf.PooledSdrf.SOURCE_DOCUMENT_COLUMN in pooled.columns
    assert set(pooled.source_documents()) <= {"malaria", "colon"}


def test_pool_renders_one_tab_separated_stdin_line_per_document(captured):
    sdrf.pool({"a.sdrf.tsv": "malaria", "b.sdrf.tsv": "colon"})

    assert captured["kwargs"]["stdin"] == "a.sdrf.tsv\tmalaria\nb.sdrf.tsv\tcolon"
    assert captured["args"][:2] == ("sdrf", "pool")


def test_pool_without_labels_sends_bare_paths(captured):
    sdrf.pool(["a.sdrf.tsv", "b.sdrf.tsv"])

    assert captured["kwargs"]["stdin"] == "a.sdrf.tsv\nb.sdrf.tsv"


def test_pool_warns_when_provenance_fell_back_to_the_path(pooled, monkeypatch):
    # The recorded fixture was pooled WITH labels, so it must not carry the warning...
    assert not any("reproducible" in c for c in pooled.caveats)

    # ...and the bridge adds it when they are absent. Pinned here as the contract the Python side
    # passes through, since a caller who never sees this warning will publish an unreproducible
    # table without knowing.
    payload = _load(POOL_FIXTURE)
    payload["caveats"] = payload["caveats"] + ["... not reproducible on another machine ..."]
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: payload)

    assert any("reproducible" in c for c in sdrf.pool(["a.sdrf.tsv", "b.sdrf.tsv"]).caveats)


def test_pool_writes_the_whole_document_not_the_window(captured):
    sdrf.pool(["a.sdrf.tsv", "b.sdrf.tsv"], out="merged.sdrf.tsv", limit=4)

    args = captured["args"]
    assert "--out" in args and args[args.index("--out") + 1] == "merged.sdrf.tsv"
    # The fixture's own numbers say it: 4 rows came back, the file holds all 24.
    assert _load(POOL_FIXTURE)["row_count"] == 24


# ---- argument validation, before a subprocess is spent ---------------------------------------


def test_a_bare_string_is_refused_rather_than_pooled_letter_by_letter():
    # A str IS a Sequence, so without this pool("a.sdrf.tsv") would send one line per character
    # and fail with "file not found: a".
    with pytest.raises(pymzlib.UsageError, match="several documents"):
        sdrf.pool("a.sdrf.tsv")


def test_an_empty_selection_is_refused():
    with pytest.raises(pymzlib.UsageError, match="At least one"):
        sdrf.pool([])


def test_a_blank_label_is_refused_rather_than_silently_defaulted():
    # Half-labelled input is the trap: the bridge refuses it too, but failing here costs no
    # subprocess and names the offending path.
    with pytest.raises(pymzlib.UsageError, match="label for"):
        sdrf.pool({"a.sdrf.tsv": "malaria", "b.sdrf.tsv": "   "})


def test_a_tab_in_a_label_is_refused_because_it_is_the_field_separator():
    with pytest.raises(pymzlib.UsageError, match="tab or newline"):
        sdrf.pool({"a.sdrf.tsv": "mal\taria"})


@pytest.mark.parametrize("path", ["", "   ", None, 7])
def test_read_requires_a_path(path):
    with pytest.raises(pymzlib.UsageError, match="file path is required"):
        sdrf.read(path)


@pytest.mark.parametrize("limit", [-1, 1.5, "3", True])
def test_read_rejects_a_bad_limit(limit):
    with pytest.raises(pymzlib.UsageError, match="limit must be"):
        sdrf.read("a.sdrf.tsv", limit=limit)


@pytest.mark.parametrize("offset", [-1, 2.5, "0", True])
def test_read_rejects_a_bad_offset(offset):
    with pytest.raises(pymzlib.UsageError, match="offset must be"):
        sdrf.read("a.sdrf.tsv", offset=offset)


def test_the_window_options_reach_the_bridge(monkeypatch):
    seen: dict = {}
    monkeypatch.setattr(
        _bridge, "invoke", lambda *a, **k: seen.update(args=a) or _load(READ_FIXTURE)
    )

    sdrf.read("a.sdrf.tsv", limit=2, offset=1)

    assert seen["args"] == ("sdrf", "read", "--path", "a.sdrf.tsv", "--limit", "2", "--offset", "1")


# ---- live canaries ---------------------------------------------------------------------------


def _bridge_is_built() -> bool:
    try:
        return pymzlib.bridge_path().exists()
    except pymzlib.BridgeNotFoundError:
        return False


MZLIB_SDRF = (
    Path(__file__).parents[3]
    / "code" / "mzLib" / "mzLib" / "Test" / "FileReadingTests" / "ExternalFileTypes"
    / "PXD000070.sdrf.tsv"
)


@pytest.mark.skipif(not _bridge_is_built(), reason="bridge binary not built")
@pytest.mark.skipif(not MZLIB_SDRF.exists(), reason="mzLib worktree fixture not present")
def test_the_recorded_read_fixture_still_matches_the_live_bridge():
    """The recorded payload must not drift from what the bridge actually emits.

    A fixture that has quietly diverged makes every offline test above prove something about a
    payload shape nothing produces any more.
    """
    live = sdrf.read(str(MZLIB_SDRF))
    recorded = sdrf.SdrfDocument._from_wire(_load(READ_FIXTURE))

    assert live.columns == recorded.columns
    assert live.rows == recorded.rows
    assert live.row_count == recorded.row_count
    assert live.caveats == recorded.caveats


@pytest.mark.skipif(not _bridge_is_built(), reason="bridge binary not built")
@pytest.mark.skipif(not MZLIB_SDRF.exists(), reason="mzLib worktree fixture not present")
def test_pooling_one_document_with_itself_keeps_both_copies_apart():
    """Two labels over the same file: provenance is what distinguishes otherwise identical rows.

    Also the cheapest live proof that the stdin rendering survives a round trip, since a pooled
    table of one repeated document is the case where a dropped line is most visible.
    """
    pooled = sdrf.pool({str(MZLIB_SDRF): "first"})
    twice = sdrf.pool([str(MZLIB_SDRF), str(MZLIB_SDRF)])

    assert pooled.document_count == 1
    assert twice.document_count == 2
    assert twice.row_count == 2 * pooled.row_count
