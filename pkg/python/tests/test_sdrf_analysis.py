"""Tests for validate, lint, assess, samples and parse_ages in :mod:`pymzlib.sdrf`.

Offline except for the live canaries at the end. Every payload is recorded from the real bridge
(``sdrf_validate_*.json`` and friends), over three small documents that sit beside them:
``sdrf_cohort.sdrf.tsv`` (six samples, one age per precision, a bare number, a mis-cased reserved
word, and a sample whose two rows disagree about its disease), ``sdrf_cohort_partner.sdrf.tsv``
(one drift of every kind against the cohort) and ``sdrf_skeleton.sdrf.tsv``. Only path fields are
rewritten to base names, so the recording machine's layout is not committed.

The Python layer's job is small and worth pinning exactly: refuse bad arguments before spending a
process, render stdin and options faithfully, and turn the payload into objects that keep every
null meaning apart - a document-level finding from a row-level one, an unbounded age from a refused
one, a withheld column from an absent one.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, sdrf

FIXTURES = Path(__file__).parent / "fixtures"
COHORT = FIXTURES / "sdrf_cohort.sdrf.tsv"
PARTNER = FIXTURES / "sdrf_cohort_partner.sdrf.tsv"
SKELETON = FIXTURES / "sdrf_skeleton.sdrf.tsv"


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    """Answer the next bridge call with a recorded payload, and remember what was sent."""
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, **kwargs):
            seen["args"] = args
            seen["kwargs"] = kwargs
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        return seen

    return use


# ---- validate -------------------------------------------------------------------------------


def test_validate_projects_mzlibs_findings(replay):
    seen = replay("sdrf_validate_cohort.json")
    result = sdrf.validate(COHORT)

    assert seen["args"] == ("sdrf", "validate", "--path", str(COHORT))
    assert result.is_valid and result.error_count == 0
    assert result.warning_count == result.message_count == 2
    assert [m.rule for m in result.warnings] == ["ReservedWordCase", "ReservedWordCase"]
    assert result.errors == []
    # 0-based rows 8 and 9 are lines 10 and 11: the header is line 1.
    assert [(m.row_index, m.line_number) for m in result.messages] == [(8, 10), (9, 11)]


def test_a_document_level_finding_has_no_row_or_line(replay):
    # None, not 0: a missing column is about the document, and 0 would read as "the first row".
    replay("sdrf_validate_skeleton.json")
    result = sdrf.validate(SKELETON)

    assert not result.is_valid
    required = [m for m in result.errors if m.rule == "RequiredColumn"]
    assert required, "the skeleton lacks most required columns"
    assert all(m.row_index is None and m.line_number is None for m in required)
    assert all(m.column_name for m in required)


def test_validate_columns_are_dataframe_ready(replay):
    replay("sdrf_validate_skeleton.json")
    result = sdrf.validate(SKELETON)

    assert list(result.columns) == result.column_names
    assert {len(v) for v in result.columns.values()} == {result.message_count}
    assert result.records[0].keys() == set(result.column_names)


def test_validate_many_sends_paths_on_stdin_with_explicit_options(replay):
    seen = replay("sdrf_validate_bulk.json")
    sdrf.validate_many([SKELETON, COHORT], threads=4, on_error="skip")

    assert seen["args"] == (
        "sdrf", "validate", "--paths-stdin", "--threads", "4", "--on-error", "skip",
    )
    assert seen["kwargs"]["stdin"] == f"{SKELETON}\n{COHORT}"


def test_validate_many_keeps_a_failed_file_visible(replay):
    replay("sdrf_validate_bulk.json")
    batch = sdrf.validate_many(["sdrf_skeleton.sdrf.tsv", "sdrf_cohort.sdrf.tsv", "missing.sdrf.tsv"],
                               on_error="skip")

    assert (batch.file_count, batch.read_count, batch.failed_count) == (3, 2, 1)
    assert batch.valid_count == 1
    skeleton, cohort, missing = batch.files
    assert skeleton.is_valid is False and cohort.is_valid is True
    assert missing.is_valid is None and missing.message_count is None
    assert missing.error == sdrf.FileError(kind="usage", message=missing.error.message)
    assert "missing.sdrf.tsv" in missing.error.message
    # The long table starts with its provenance and is in input order.
    assert batch.column_names[:2] == ["source_index", "source_path"]
    assert batch.columns["source_index"] == sorted(batch.columns["source_index"])
    assert 2 not in batch.columns["source_index"]
    assert len(batch.messages) == batch.message_count


# ---- lint -----------------------------------------------------------------------------------


def test_lint_renders_documents_exactly_as_pool_does(replay):
    seen = replay("sdrf_lint_cohort.json")
    sdrf.lint({"a.sdrf.tsv": "cohort", "b.sdrf.tsv": "partner"})

    assert seen["args"] == ("sdrf", "lint")
    assert seen["kwargs"]["stdin"] == "a.sdrf.tsv\tcohort\nb.sdrf.tsv\tpartner"


def test_lint_is_one_row_per_variant_grouped_by_finding(replay):
    replay("sdrf_lint_cohort.json")
    drift = sdrf.lint({"a.sdrf.tsv": "cohort", "b.sdrf.tsv": "partner"})
    findings = drift.findings()

    assert drift.finding_count == len(findings) == 4
    assert {f[0]["kind"] for f in findings} == {
        "AccessionNameConflict", "MixedTermAndFreeText", "ColumnNameVariant", "ValueCaseVariant",
    }
    for variants in findings:
        assert [v["variant_rank"] for v in variants] == list(range(len(variants)))
        assert all(isinstance(v["documents"], list) for v in variants)


def test_a_column_name_finding_has_no_column(replay):
    replay("sdrf_lint_cohort.json")
    drift = sdrf.lint(["a.sdrf.tsv", "b.sdrf.tsv"])
    names = next(f for f in drift.findings() if f[0]["kind"] == "ColumnNameVariant")

    assert {v["column_name"] for v in names} == {None}
    assert {v["value"] for v in names} == {"characteristics[sex]", "Characteristics[sex]"}


def test_lint_never_reports_searched_data_file_names(replay):
    # mzLib #1335. The fixtures write S1_F1.mzML and s1_f1.mzML on purpose.
    replay("sdrf_lint_cohort.json")
    drift = sdrf.lint(["a.sdrf.tsv", "b.sdrf.tsv"])

    assert "comment[searched data file]" not in drift.columns["column_name"]


def test_lint_refuses_what_pool_refuses():
    with pytest.raises(pymzlib.UsageError, match=r"lint\(\) takes several documents"):
        sdrf.lint("a.sdrf.tsv")
    with pytest.raises(pymzlib.UsageError, match="label for"):
        sdrf.lint({"a.sdrf.tsv": "x", "b.sdrf.tsv": " "})


# ---- assess ---------------------------------------------------------------------------------


def test_assess_carries_the_verdict_and_its_evidence(replay):
    seen = replay("sdrf_assess_cohort.json")
    result = sdrf.assess(COHORT)

    assert seen["args"] == ("sdrf", "assess", "--path", str(COHORT))
    assert result.verdict == "Informative" and result.verdict in sdrf.VERDICTS
    assert result.factor_value_varies and result.sample_is_described
    assert result.biological_replicate_varies
    assert result.columns["role"][0] == "factor_value"
    assert result.columns["role"][-1] == "biological_replicate"
    assert all(0 <= r <= 1 for r in result.columns["fill_rate"])
    # characteristics[cell type] is all "not applicable": present, but absent in every row.
    cell_type = result.columns["column_name"].index("characteristics[cell type]")
    assert result.columns["filled"][cell_type] == 0
    assert result.columns["absent"][cell_type] == result.columns["rows"][cell_type]


def test_assess_many_counts_verdicts_and_selects_by_them(replay):
    replay("sdrf_assess_bulk.json")
    batch = sdrf.assess_many(["cohort", "skeleton", "PXD000070"])

    assert batch.verdict_counts == {"informative": 1, "partial": 1, "skeleton": 1}
    assert [f.verdict for f in batch.files] == ["Informative", "Skeleton", "Partial"]
    assert batch.paths_with("Informative", "Partial") == [batch.files[0].path, batch.files[2].path]
    with pytest.raises(pymzlib.UsageError, match="Unknown verdict"):
        batch.paths_with("informative")


# ---- samples --------------------------------------------------------------------------------


def test_samples_is_one_sample_per_source_name_in_document_order(replay):
    replay("sdrf_samples_cohort.json")
    result = sdrf.samples(COHORT)

    assert result.sample_count == 6 and result.row_count == 12
    assert list(dict.fromkeys(result.columns["source_name"])) == ["S1", "S2", "S3", "S4", "S5", "S6"]
    assert set(result.columns["sample_row_count"]) == {2}
    assert result.problems == []


def test_a_conflicting_column_is_withheld_but_named(replay):
    replay("sdrf_samples_cohort.json")
    result = sdrf.samples(COHORT)

    assert result.conflict_count == 1
    assert result.conflicts() == [("S6", "characteristics[disease]")]
    row = next(r for r in result.records if r["status"] == "conflicting")
    assert row["value"] is None and row["position"] is None


def test_ages_are_in_years_and_every_null_is_explained(replay):
    replay("sdrf_samples_cohort.json")
    ages = {r["source_name"]: r for r in sdrf.samples(COHORT).ages()}

    assert ages["S1"]["age_years"] == 58 and ages["S1"]["age_precision"] == "Exact"
    assert ages["S2"]["age_years"] == 62.5 and ages["S2"]["age_precision"] == "Range"
    assert (ages["S2"]["age_min_years"], ages["S2"]["age_max_years"]) == (40, 85)
    # >=90Y: max is null because it is UNBOUNDED, and precision is what says so.
    assert ages["S3"]["age_max_years"] is None and ages["S3"]["age_precision"] == "LowerBound"
    assert ages["S3"]["age_refusal"] is None
    # Refused cells: no age, and a reason.
    assert ages["S4"]["age_years"] is None and ages["S4"]["age_refusal"] == "no_unit"
    assert ages["S5"]["age_refusal"] == "reserved_word"
    # mzLib #1333: a degenerate range pins one age.
    assert ages["S6"]["age_precision"] == "Exact"
    assert all(r["age_refusal"] in (None, *sdrf.AGE_REFUSALS) for r in ages.values())


def test_non_age_rows_carry_no_age_at_all(replay):
    replay("sdrf_samples_cohort.json")
    rows = [r for r in sdrf.samples(COHORT).records if r["column_name"] != "characteristics[age]"]

    assert rows
    assert all(r[f"age_{f}"] is None for r in rows for f in
               ("years", "min_years", "max_years", "precision", "follows_specification", "refusal"))


def test_samples_many_keys_samples_by_source(replay):
    seen = replay("sdrf_samples_bulk.json")
    batch = sdrf.samples_many([COHORT, PARTNER], threads=-1)

    assert seen["args"][-4:] == ("--threads", "-1", "--on-error", "fail")
    assert batch.sample_count == 8
    assert [f.sample_count for f in batch.files] == [6, 2]
    assert batch.files[0].problems == [] and batch.files[0].error is None
    keys = set(zip(batch.columns["source_index"], batch.columns["source_name"]))
    assert (0, "S1") in keys and (1, "P1") in keys


# ---- parse_ages -----------------------------------------------------------------------------


def test_parse_ages_is_one_row_per_cell_in_order(replay):
    seen = replay("sdrf_parse_age.json")
    cells = ["58Y", "30Y6M", "40Y-85Y", "40Y-40Y", ">=90Y", "<1Y", "6-8 weeks", "63",
             "not available", None, "about forty"]
    ages = sdrf.parse_ages(cells)

    assert seen["args"] == ("sdrf", "parse-age")
    # None goes as an empty cell, and a trailing newline keeps a final blank cell alive.
    assert seen["kwargs"]["stdin"] == "\n".join(c or "" for c in cells) + "\n"
    assert ages.cell_count == len(cells) and ages.parsed_count == 7
    assert ages.columns["cell"] == [c or "" for c in cells]
    assert ages.columns["precision"][:7] == ["Exact", "Exact", "Range", "Exact", "LowerBound",
                                             "UpperBound", "Range"]
    assert ages.columns["refusal"][7:] == ["no_unit", "reserved_word", "empty", "unreadable"]
    assert ages.columns["years"][1] == 30.5
    assert ages.columns["min_years"][5] == 0 and ages.columns["max_years"][5] == 1
    assert ages.columns["max_years"][4] is None
    assert ages.columns["follows_specification"][6] is False
    assert ages.columns["years"][6] == pytest.approx(7 * 7 / 365.25)


def test_parse_ages_accepts_a_document_column_directly(replay):
    # The intended pipeline: read() then parse_ages(doc.value(...)). value() yields None for a row
    # too short to reach the column, which must not shift later ages.
    seen = replay("sdrf_parse_age.json")
    sdrf.parse_ages(["58Y", None, "63"])
    assert seen["kwargs"]["stdin"] == "58Y\n\n63\n"


@pytest.mark.parametrize("cells, match", [
    ("58Y", "takes a list of cells"),
    ([], "at least one cell"),
    ([58], "str or None"),
    (["58Y\n63"], "newline"),
])
def test_parse_ages_rejects_bad_input_without_a_process(cells, match):
    with pytest.raises(pymzlib.UsageError, match=match):
        sdrf.parse_ages(cells)


# ---- argument validation shared by the *_many forms -----------------------------------------


@pytest.mark.parametrize("verb", ["validate_many", "assess_many", "samples_many"])
def test_a_single_path_is_pointed_at_the_single_form(verb):
    with pytest.raises(pymzlib.UsageError, match=r"takes a list of paths"):
        getattr(sdrf, verb)("a.sdrf.tsv")
    with pytest.raises(pymzlib.UsageError, match=r"takes a list of paths"):
        getattr(sdrf, verb)(Path("a.sdrf.tsv"))


@pytest.mark.parametrize("kwargs, match", [
    ({"threads": 0}, "threads must be"),
    ({"threads": -2}, "threads must be"),
    ({"threads": True}, "threads must be"),
    ({"threads": 1.5}, "threads must be"),
    ({"on_error": "ignore"}, "on_error must be"),
])
def test_bad_bulk_options_are_refused(kwargs, match):
    with pytest.raises(pymzlib.UsageError, match=match):
        sdrf.validate_many(["a.sdrf.tsv"], **kwargs)


@pytest.mark.parametrize("paths, match", [
    ([], "At least one"),
    (["a.sdrf.tsv", "  "], "non-blank"),
    (["a\nb.sdrf.tsv"], "newline"),
    ({"a.sdrf.tsv"}, "list of file paths"),
])
def test_bad_path_lists_are_refused(paths, match):
    with pytest.raises(pymzlib.UsageError, match=match):
        sdrf.assess_many(paths)


@pytest.mark.parametrize("verb", ["validate", "assess", "samples"])
@pytest.mark.parametrize("path", ["", "   ", None, 7])
def test_the_single_forms_require_a_path(verb, path):
    with pytest.raises(pymzlib.UsageError, match="file path is required"):
        getattr(sdrf, verb)(path)


def test_file_error_is_none_for_a_file_that_was_read():
    assert sdrf.FileError._from_wire(None) is None
    assert sdrf.FileError._from_wire({}) is None


# ---- live canaries --------------------------------------------------------------------------


def _bridge_is_built() -> bool:
    try:
        return pymzlib.bridge_path().exists()
    except pymzlib.BridgeNotFoundError:
        return False


live = pytest.mark.skipif(not _bridge_is_built(), reason="bridge binary not built")


def _same(live_data, recorded_name: str, cls) -> None:
    """A recorded fixture must still be what the bridge emits, paths aside."""
    recorded = cls._from_wire(_load(recorded_name))
    assert live_data.column_names == recorded.column_names
    for name in recorded.columns:
        if name != "source_path":
            assert live_data.columns[name] == recorded.columns[name], name
    assert live_data.caveats == recorded.caveats


@live
def test_the_recorded_validate_fixture_still_matches_the_live_bridge():
    _same(sdrf.validate(COHORT), "sdrf_validate_cohort.json", sdrf.SdrfValidation)
    _same(sdrf.validate(SKELETON), "sdrf_validate_skeleton.json", sdrf.SdrfValidation)


@live
def test_the_recorded_assess_and_samples_fixtures_still_match_the_live_bridge():
    _same(sdrf.assess(COHORT), "sdrf_assess_cohort.json", sdrf.SdrfAssessment)
    _same(sdrf.samples(COHORT), "sdrf_samples_cohort.json", sdrf.SdrfSamples)
    _same(sdrf.samples_many([COHORT, PARTNER]), "sdrf_samples_bulk.json", sdrf.SdrfSamplesBatch)


@live
def test_the_recorded_lint_and_age_fixtures_still_match_the_live_bridge():
    drift = sdrf.lint({str(COHORT): "cohort", str(PARTNER): "partner"})
    _same(drift, "sdrf_lint_cohort.json", sdrf.SdrfDrift)
    cells = ["58Y", "30Y6M", "40Y-85Y", "40Y-40Y", ">=90Y", "<1Y", "6-8 weeks", "63",
             "not available", "", "about forty"]
    _same(sdrf.parse_ages(cells), "sdrf_parse_age.json", sdrf.ParsedAges)


@live
def test_a_live_bulk_call_is_the_same_at_any_thread_count():
    paths = [COHORT, SKELETON, PARTNER]
    one = sdrf.assess_many(paths, threads=1)
    every = sdrf.assess_many(paths, threads=-1)
    assert one == every


@live
def test_live_skip_reports_a_missing_file_and_fail_raises():
    missing = FIXTURES / "no-such.sdrf.tsv"
    batch = sdrf.validate_many([COHORT, missing], on_error="skip")
    assert batch.files[1].error is not None and batch.files[1].error.kind == "usage"
    with pytest.raises(pymzlib.UsageError, match="no-such"):
        sdrf.validate_many([COHORT, missing])
