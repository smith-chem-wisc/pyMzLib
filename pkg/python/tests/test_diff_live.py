"""diff design and diff observations through the real, packaged bridge.

The offline tests replay recorded payloads, so they cannot see the bridge and this module drift
apart. These run the real executable on the files in ``fixtures/diff/`` and check that what it
returns today is what was recorded (so the docs' printed numbers are still true), that the files
it writes compose with stats.fit from Python, and that a refusal crosses as a usage error. Local
files only.
"""

from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath

import pytest

import pymzlib
from pymzlib import _bridge, diff, stats

FIXTURES = Path(__file__).parent / "fixtures"
DIFF = FIXTURES / "diff"
STATS = FIXTURES / "stats"

#: Fields that hold a path, which differs between the recording machine and this one.
PATHS = {"samples_file", "peptides_file", "protein_groups_file", "design_file", "contrasts_file", "directory", "path"}


@pytest.fixture(autouse=True)
def built_bridge():
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so the diff verbs went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def file_name(path: str) -> str:
    """The last part of a path written on either platform: a recording made on Windows is replayed on Linux."""
    return PureWindowsPath(path).name if "\\" in path else PurePosixPath(path).name


def without_paths(value):
    if isinstance(value, dict):
        return {k: (file_name(v) if k in PATHS and isinstance(v, str) else without_paths(v)) for k, v in value.items()}
    if isinstance(value, list):
        return [without_paths(v) for v in value]
    return value


def recorded(name: str) -> dict:
    return without_paths(json.loads((FIXTURES / name).read_text(encoding="utf-8")))


def mm_runs() -> list[tuple[str, str]]:
    header = (DIFF / "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv").read_text(encoding="utf-8").splitlines()[0]
    names = [c[len("Intensity_"):].removesuffix("-calib") for c in header.split("\t") if c.startswith("Intensity_")]
    return [(n, n) for n in names]


SPIKE_IN = [("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"), ("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x")]


def test_the_design_recordings_are_still_what_the_bridge_says(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cases = {
        "diff_design_malat.json": dict(samples=DIFF / "malat_dilution_samples.tsv",
                                       factors=[diff.FactorSpec("malat_load", "500ng")], out="malat_design"),
        "diff_design_factor_batch_covariate.json": dict(
            samples=DIFF / "design_factor_batch_covariate_samples.tsv", factors=[diff.FactorSpec("condition", "ctrl")],
            covariates=[diff.CovariateSpec("age", "years", scale=10, scaled_unit="decade", centre=50)], batch="batch"),
        "diff_design_confounded.json": dict(samples=DIFF / "design_confounded_samples.tsv",
                                            factors=[diff.FactorSpec("condition", "A")], batch="batch"),
        "diff_design_strata.json": dict(samples=DIFF / "stratified_samples.tsv",
                                        factors=[diff.FactorSpec("condition", "A")], strata=["organism part"]),
    }
    for name, kwargs in cases.items():
        live = diff.design(**kwargs)
        assert without_paths(_as_wire(live)) == _as_wire(diff.Design._from_wire(recorded(name))), name


def _as_wire(design: diff.Design) -> dict:
    """Enough of a Design to compare two answers field by field."""
    return {
        "strata": [vars(s) for s in design.strata],
        "columns": design.columns,
        "caveats": design.caveats,
        "estimability_counts": design.estimability_counts,
        "unused_columns": design.unused_columns,
    }


def test_the_observation_recordings_are_still_what_the_bridge_says(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cases = {
        "diff_observations_mm_mbr_kept.json": ("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", mm_runs(), dict(
            basis="mbr_kept", protein_groups=DIFF / "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")),
        "diff_observations_mm_msms_only.json": ("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", mm_runs(), dict(
            basis="msms_only", protein_groups=DIFF / "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")),
        "diff_observations_spike_in.json": ("EcoliSpikeIn_QuantifiedPeptides.tsv", SPIKE_IN, dict(basis="msms_only")),
        "diff_observations_spike_in_normalized.json": ("EcoliSpikeIn_QuantifiedPeptides.tsv", SPIKE_IN, dict(
            basis="msms_only", min_peptides=50, normalize=True, responses_out="spike_in_responses.tsv")),
    }
    for name, (peptides, runs, kwargs) in cases.items():
        live = diff.observations(DIFF / peptides, runs, **kwargs)
        expected = diff.Observations._from_wire(recorded(name))
        assert live.columns == expected.columns, name
        assert live.normalization == expected.normalization, name
        assert live.state_counts == expected.state_counts, name
        assert live.caveats == expected.caveats, name


def test_design_then_fit_is_the_committed_fit(tmp_path):
    """The design file diff design writes is the same model as the committed 0/1 design."""
    design = diff.design(DIFF / "malat_dilution_samples.tsv", factors=[diff.FactorSpec("malat_load", "500ng")],
                         out=tmp_path)
    ours = stats.fit(STATS / "malat_dilution_log2_vs_fluc.tsv", design.strata[0].design_file,
                     ["malat_load=250ng", "malat_load=125ng"])
    committed = stats.fit(STATS / "malat_dilution_log2_vs_fluc.tsv", STATS / "malat_dilution_design.tsv",
                          ["malat_250ng", "malat_125ng"])
    assert ours.status_counts == committed.status_counts
    for a, b in zip(ours.columns["estimate"], committed.columns["estimate"]):
        assert (a is None) == (b is None)
        if a is not None:
            assert a == pytest.approx(b, rel=1e-9, abs=1e-12)


def test_the_responses_table_is_what_stats_fit_reads(tmp_path):
    responses = tmp_path / "responses.tsv"
    table = diff.observations(DIFF / "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", mm_runs(), basis="mbr_kept",
                              responses_out=responses)
    design = tmp_path / "design.tsv"
    design.write_text("sample\tintercept\n" + "".join(f"{s}\t1\n" for s in table.samples), encoding="utf-8")
    fit = stats.fit(responses, design, ["intercept"])
    assert fit.feature_count == table.peptide_count and fit.sample_names == table.samples


def test_a_refusal_crosses_as_a_usage_error(tmp_path):
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tcondition\ns1\tA\ns2\t\n", encoding="utf-8")
    with pytest.raises(pymzlib.UsageError, match="has no value for factor 'condition'"):
        diff.design(sheet, factors=[diff.FactorSpec("condition", "A")])
    with pytest.raises(pymzlib.UsageError, match="not in the experimental design"):
        diff.observations(DIFF / "EcoliSpikeIn_QuantifiedPeptides.tsv", SPIKE_IN[:1], basis="msms_only")
