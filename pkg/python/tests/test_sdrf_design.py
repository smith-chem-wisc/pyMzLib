"""Tests for :func:`pymzlib.sdrf.design`, the projection of mzLib's ``SdrfLabelFreeDesign``.

Offline except for the live canary at the end. The payloads are recorded from the real bridge over
mzLib's own design fixtures, copied beside them: ``PXD067622.sdrf.tsv`` (24 runs, two factors),
``PXD049018.sdrf.tsv`` (treatment ``not available``, refused) and ``PXD067622_studywide.sdrf.tsv``
(PXD067622 numbered by the study-wide index in each file name, which mzLib ranks back per
condition).

What the Python layer must get exactly right: the condition columns cross as ONE tab-joined
option, searched files cross on stdin with the flag, ``out`` is passed through, a refusal stays a
result (only ``spectra()``/``run_design()`` raise), and the 0-based coordinates reach FlashLFQ's
shapes unchanged.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, flashlfq, sdrf

FIXTURES = Path(__file__).parent / "fixtures"
VALID = FIXTURES / "PXD067622.sdrf.tsv"
REFUSED = FIXTURES / "PXD049018.sdrf.tsv"
STUDYWIDE = FIXTURES / "PXD067622_studywide.sdrf.tsv"
BOTH = ["factor value[genotype]", "factor value[treatment]"]


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, **kwargs):
            seen["args"] = args
            seen["kwargs"] = kwargs
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        return seen

    return use


def test_a_valid_design_projects_every_run(replay):
    seen = replay("sdrf_design_PXD067622.json")
    d = sdrf.design(VALID, condition_columns=BOTH)

    assert seen["args"] == (
        "sdrf", "design", "--path", str(VALID),
        "--condition-columns", "factor value[genotype]\tfactor value[treatment]",
    )
    assert seen["kwargs"]["stdin"] is None
    assert d.is_valid and d.file_count == 24 and len(d.files) == 24
    assert d.refusals == [] and d.notes == []
    assert d.condition_columns == BOTH and d.condition_columns_declared
    assert not d.searched_files_given
    assert d.file_key_column == "comment[data file]"
    assert d.written is None
    assert len(set(d.columns["condition"])) == 8
    assert {f.biological_replicate for f in d.files} == {0, 1, 2}
    assert all(f.full_path == f.file_name + ".raw" for f in d.files)
    assert "24 file(s), 8 condition(s)" in d.report
    assert len(d.caveats) == 6


def test_the_design_feeds_flashlfq_unchanged(replay):
    replay("sdrf_design_PXD067622.json")
    d = sdrf.design(VALID, condition_columns=BOTH)

    spectra = d.spectra()
    assert spectra[0] == {
        "path": "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4.raw",
        "condition": "SPRTN-TurboID CA_DMSO (vehicle)",
        "biological_replicate": 0,
        "technical_replicate": 0,
        "fraction": 0,
    }
    # The run design is keyed by file_name, as median_polish wants it.
    assert d.run_design()[0]["file_name"] == "20240830_HF_LC3_MAA_RK_12032_CA_DMSO4"
    # And FlashLFQ's own stdin renderer accepts every row: 0-based, non-negative whole numbers.
    lines = flashlfq._spectra_stdin(spectra).splitlines()
    assert len(lines) == 24
    assert lines[0].split("\t")[1:] == ["SPRTN-TurboID CA_DMSO (vehicle)", "0", "0", "0"]


def test_a_refusal_is_a_result_and_blocks_only_the_hand_off(replay):
    replay("sdrf_design_PXD049018.json")
    d = sdrf.design(REFUSED, condition_columns=BOTH)

    assert not d.is_valid
    assert len(d.refusals) == 20 and d.file_count == 0 and d.files == []
    assert "MSB67868ABand_01.raw" in d.refusals[0] and "'not available'" in d.refusals[0]
    with pytest.raises(_bridge.UsageError, match="refused") as caught:
        d.spectra()
    assert d.refusals[0] in str(caught.value), "the reasons travel with the error"
    with pytest.raises(_bridge.UsageError, match="refused"):
        d.run_design()


def test_the_renumbering_is_kept_in_notes(replay):
    replay("sdrf_design_studywide.json")
    ranked = sdrf.design(STUDYWIDE, condition_columns=BOTH)
    replay("sdrf_design_PXD067622.json")
    reference = sdrf.design(VALID, condition_columns=BOTH)

    assert ranked.is_valid and len(ranked.notes) == 7
    assert any("22 -> 1, 23 -> 2, 24 -> 3" in n for n in ranked.notes)
    assert ranked.columns == reference.columns


def test_no_condition_columns_sends_no_option(replay):
    seen = replay("sdrf_design_PXD067622.json")
    sdrf.design(VALID)
    assert seen["args"] == ("sdrf", "design", "--path", str(VALID))


def test_searched_files_go_on_stdin_with_the_flag(replay, tmp_path):
    seen = replay("sdrf_design_PXD067622.json")
    runs = [tmp_path / "a.raw", "b.raw"]
    sdrf.design(VALID, condition_columns=BOTH, searched_files=runs)

    assert seen["args"][-1] == "--searched-files-stdin"
    assert seen["kwargs"]["stdin"] == f"{os.fspath(runs[0])}\nb.raw\n"


def test_out_is_passed_through(replay, tmp_path):
    seen = replay("sdrf_design_PXD067622.json")
    target = tmp_path / "ExperimentalDesign.tsv"
    sdrf.design(VALID, condition_columns=BOTH, out=target)
    assert seen["args"][-2:] == ("--out", str(target))


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"condition_columns": "factor value[treatment]"}, "list of column names"),
        ({"condition_columns": []}, "is empty"),
        ({"condition_columns": ["  "]}, "non-blank"),
        ({"condition_columns": ["a\tb"]}, "tab or newline"),
        ({"searched_files": "run.raw"}, "list of paths"),
        ({"searched_files": [""]}, "non-blank"),
        ({"searched_files": ["a\nb"]}, "newline"),
        ({"out": ""}, "ending in .tsv"),
    ],
)
def test_bad_arguments_are_refused_before_the_bridge(monkeypatch, kwargs, message):
    def no_call(*args, **kw):
        raise AssertionError("the bridge must not be called")

    monkeypatch.setattr(_bridge, "invoke", no_call)
    with pytest.raises(_bridge.UsageError, match=message):
        sdrf.design(VALID, **kwargs)


def test_a_blank_path_is_refused(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("called"))
    with pytest.raises(_bridge.UsageError):
        sdrf.design("")


# ---- live: the real bridge ----------------------------------------------------------------------


@pytest.fixture()
def built_bridge():
    # As test_proteins_live: skip in a source checkout with no built bridge, fail under CI, where
    # the bridge is always built and a skip would hide drift between the recording and mzLib.
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so sdrf design went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def test_live_design_matches_the_recording_and_writes_a_one_based_file(built_bridge, tmp_path):
    out = tmp_path / "ExperimentalDesign.tsv"
    d = sdrf.design(VALID, condition_columns=BOTH, out=out)
    recorded = _load("sdrf_design_PXD067622.json")

    assert d.columns == recorded["columns"]
    assert d.written == {"path": str(out), "file_count": 24}
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "FileName\tCondition\tBiorep\tFraction\tTechrep"
    first = lines[1].split("\t")
    assert int(first[2]) == d.files[0].biological_replicate + 1

    with pytest.raises(_bridge.UsageError, match=r"\.tsv"):
        sdrf.design(VALID, out=tmp_path / "design.txt")
