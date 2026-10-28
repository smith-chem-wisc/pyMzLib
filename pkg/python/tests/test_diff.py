"""pymzlib.diff.design and observations against recorded bridge payloads.

The numbers are mzLib's, held to R's model.matrix and to mzLib's own tables by the bridge's C#
tests (pkg/bridge.tests/DifferentialTests.cs) and by test_diff_live.py. These tests hold the Python
layer: that it sends what the wire expects, sends nothing for what the caller left out, parses what
the wire returns, and refuses bad input before starting a process.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, diff

FIXTURES = Path(__file__).parent / "fixtures"
DIFF = FIXTURES / "diff"


def recorded(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture
def bridge(monkeypatch):
    """Answer every bridge call with ``state["payload"]``, recording the arguments and stdin."""
    state = {"calls": [], "payload": None}

    def invoke(*args, stdin=None, timeout=None):
        state["calls"].append({"args": list(args), "stdin": stdin, "timeout": timeout})
        return state["payload"]

    monkeypatch.setattr(_bridge, "invoke", invoke)
    monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: None)
    return state


def test_the_module_is_exported():
    assert pymzlib.diff is diff
    assert "diff" in pymzlib.__all__
    assert set(diff.__all__) <= set(dir(diff))


# ---- design: what is sent --------------------------------------------------------------------------


def test_design_writes_one_stdin_line_per_term_and_leaves_out_what_was_not_given(bridge):
    bridge["payload"] = recorded("diff_design_factor_batch_covariate.json")
    diff.design(
        "sheet.tsv",
        factors=[diff.FactorSpec("condition", "ctrl"), "sex"],
        covariates=[diff.CovariateSpec("age", "years", scale=10, scaled_unit="decade", centre=50),
                    diff.CovariateSpec("bmi", "kg/m2"),
                    diff.CovariateSpec("dose", "mg", centre=0)],
        interactions=[("condition", "sex")],
        strata="tissue",
        batch="plex",
    )
    call = bridge["calls"][0]
    assert call["args"] == ["diff", "design", "--samples", "sheet.tsv"]
    assert call["stdin"].split("\n") == [
        "factor\tcondition\tctrl",
        "factor\tsex",
        "covariate\tage\tyears\t10.0\tdecade\t50.0",
        "covariate\tbmi\tkg/m2",
        "covariate\tdose\tmg\t\t\t0.0",
        "interaction\tcondition\tsex",
        "stratum\ttissue",
        "batch\tplex",
        "",
    ]


def test_design_sends_out_only_when_given(bridge):
    bridge["payload"] = recorded("diff_design_malat.json")
    diff.design("s.tsv", factors=["malat_load"], out=Path("dir"))
    diff.design("s.tsv")
    assert bridge["calls"][0]["args"][-2:] == ["--out", "dir"]
    assert "--out" not in bridge["calls"][1]["args"]
    assert bridge["calls"][1]["stdin"] == "", "no model: nothing on stdin, mzLib's empty DesignSpec"


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"factors": "condition"}, "takes a list"),
        ({"factors": [3]}, "FactorSpec or a column name"),
        ({"factors": [diff.FactorSpec("a\tb")]}, "tab or a line break"),
        ({"factors": [diff.FactorSpec("")]}, "non-empty str"),
        ({"covariates": [("age", "years")]}, "CovariateSpec"),
        ({"covariates": [diff.CovariateSpec("age", "years", scale=float("inf"))]}, "finite number"),
        ({"covariates": [diff.CovariateSpec("age", "years", centre=True)]}, "finite number"),
        ({"interactions": [("a", "b", "c")]}, "pair of factor names"),
        ({"interactions": ["ab"]}, "pair of factor names"),
        ({"batch": ""}, "non-empty str"),
    ],
)
def test_design_refuses_a_malformed_model_before_calling(bridge, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        diff.design("s.tsv", **kwargs)
    assert bridge["calls"] == []


def test_design_refuses_a_missing_path(bridge):
    with pytest.raises(pymzlib.UsageError, match="samples must be a path"):
        diff.design(None)
    assert bridge["calls"] == []


# ---- design: what is parsed ------------------------------------------------------------------------


def test_design_parses_the_malat_recording(bridge):
    bridge["payload"] = recorded("diff_design_malat.json")
    design = diff.design("malat_dilution_samples.tsv", factors=[diff.FactorSpec("malat_load", "500ng")], out="d")
    assert isinstance(design, diff.Design)
    assert design.factors == [diff.FactorSpec("malat_load", "500ng")]
    assert design.covariates == [] and design.interactions == [] and design.batch_column is None
    assert design.stratum_count == 1 and design.sample_count == 27
    assert design.estimability_counts == {e: n for e, n in zip(diff.ESTIMABILITIES, (2, 0, 0))}
    assert isinstance(design.written, diff.WrittenDesign) and design.written.file_count == 2
    stratum = design.stratum("all")
    assert isinstance(stratum, diff.DesignStratum)
    assert stratum.column_names == stratum.fitted_column_names == ["intercept", "malat_load=125ng", "malat_load=250ng"]
    assert len(stratum.matrix) == 27 and stratum.matrix == stratum.fitted_matrix
    assert stratum.design_file.endswith("design_1.tsv") and stratum.contrasts_file.endswith("contrasts_1.tsv")
    assert [r["fitted_weights"] for r in design.records] == [[0, 1, 0], [0, 0, 1]]
    assert design.column_names == [
        "stratum", "contrast_id", "label", "numerator", "denominator", "covariate", "covariate_unit",
        "covariate_scale", "estimability", "weights", "fitted_weights",
    ]
    with pytest.raises(pymzlib.UsageError, match="No stratum 'liver'"):
        design.stratum("liver")


def test_design_matrix_is_rs_model_matrix():
    payload = recorded("diff_design_factor_batch_covariate.json")
    design = diff.Design._from_wire(payload)
    with open(DIFF / "design_factor_batch_covariate_matrix.tsv", newline="", encoding="utf-8") as fh:
        r = [[float(v) for v in row] for row in list(csv.reader(fh, delimiter="\t"))[1:]]
    stratum = design.strata[0]
    assert [[round(v, 12) for v in row] for row in stratum.matrix] == [[round(v, 12) for v in row] for row in r]
    assert stratum.covariate_centres == {"age": 50}
    assert design.covariates == [diff.CovariateSpec("age", "years", 10, "decade", 50)]
    slope = design.records[1]
    assert (slope["covariate"], slope["covariate_unit"], slope["covariate_scale"]) == (
        "age", "decade", "years/10, centred at 50 years")


def test_design_confounded_recording_says_why():
    design = diff.Design._from_wire(recorded("diff_design_confounded.json"))
    stratum = design.strata[0]
    assert stratum.rank == 2 and stratum.aliased_terms == ["condition", "batch"]
    assert stratum.fitted_column_names == ["intercept", "condition=B"]
    assert design.records[0]["estimability"] == "not_estimable:confounded"
    assert design.records[0]["fitted_weights"] is None
    assert any("not estimable" in c for c in design.caveats)


def test_design_strata_recording_says_what_cannot_run():
    design = diff.Design._from_wire(recorded("diff_design_strata.json"))
    assert [s.stratum for s in design.strata] == ["organism part=kidney", "organism part=liver", "organism part=muscle"]
    kidney = design.stratum("organism part=kidney")
    # mzLib #1476: kidney skips the analysis-wide c1, listed under its own label and id.
    assert kidney.not_run == [diff.NotRunContrast("condition=B vs condition=A", "not_in_stratum", "c1")]
    assert kidney.column_names == ["intercept"]
    assert design.columns["contrast_id"] == ["c1", "c1"]
    assert design.stratum_columns == ["organism part"]


# ---- observations: what is sent --------------------------------------------------------------------


def test_observations_write_one_line_per_run_and_send_only_what_was_given(bridge):
    bridge["payload"] = recorded("diff_observations_spike_in.json")
    diff.observations("p.tsv", [("a", "S1"), ("b", "S1", 1), diff.ObservationRun("c", "S1", technical_replicate=1),
                                ["d", "S2", None, 2]], basis="msms_only")
    call = bridge["calls"][0]
    assert call["args"] == ["diff", "observations", "--peptides", "p.tsv", "--basis", "msms_only"]
    assert call["stdin"] == "a\tS1\nb\tS1\t1\nc\tS1\t\t1\nd\tS2\t\t2\n"


def test_observations_send_every_option_given(bridge):
    bridge["payload"] = recorded("diff_observations_spike_in_normalized.json")
    diff.observations("p.tsv", [("a", "S1")], basis="mbr_kept", protein_groups="g.tsv", samples="s.tsv",
                      normalize=True, min_peptides=50, out="o.tsv", responses_out=Path("r.tsv"))
    assert bridge["calls"][0]["args"] == [
        "diff", "observations", "--peptides", "p.tsv", "--basis", "mbr_kept", "--protein-groups", "g.tsv",
        "--samples", "s.tsv", "--normalize", "--min-peptides", "50", "--out", "o.tsv", "--responses-out", "r.tsv",
    ]


@pytest.mark.parametrize(
    "runs, kwargs, message",
    [
        ([("a", "S1")], {"basis": "mbr"}, "basis must be one of msms_only, mbr_kept"),
        ([("a", "S1")], {"basis": None}, "basis must be one of"),
        ([], {"basis": "msms_only"}, "needs the runs"),
        ([("a",)], {"basis": "msms_only"}, "Each run is"),
        ("a\tS1", {"basis": "msms_only"}, "takes a list"),
        ([("a", "S1", 1.5)], {"basis": "msms_only"}, "fraction must be an int"),
        ([("a", "S1", 0, True)], {"basis": "msms_only"}, "technical replicate must be an int"),
        ([("a b\t", "S1")], {"basis": "msms_only"}, "tab or a line break"),
        ([("a", "S1")], {"basis": "msms_only", "normalize": "yes"}, "normalize must be True or False"),
        ([("a", "S1")], {"basis": "msms_only", "min_peptides": 50.0}, "min_peptides must be an int"),
    ],
)
def test_observations_refuse_bad_input_before_calling(bridge, runs, kwargs, message):
    with pytest.raises(pymzlib.UsageError, match=message):
        diff.observations("p.tsv", runs, **kwargs)
    assert bridge["calls"] == []


# ---- observations: what is parsed ------------------------------------------------------------------


def test_observations_parse_the_metamorpheus_recording():
    table = diff.Observations._from_wire(recorded("diff_observations_mm_mbr_kept.json"))
    assert table.basis == "mbr_kept" and table.normalize is False and table.min_peptides is None
    assert table.peptide_count == 3 and table.sample_count == table.run_count == 18
    assert table.row_count == 54 == len(table.records)
    assert list(table.state_counts) == list(diff.OBSERVATION_STATES)
    assert table.runs[0] == diff.ObservationRun("QE-002106_GM1_a", "QE-002106_GM1_a", 0, 0)
    alb = next(g for g in table.protein_groups if g.name == "P02769")
    assert alb == diff.ProteinGroupInfo("P02769", 0, False, True, False, "ALB")
    assert [g.name for g in table.protein_groups if g.is_decoy] == ["DECOY_P62750"]
    assert table.protein_groups_file.endswith("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv")
    assert table.samples_file is None and table.written is None and table.responses_written is None
    assert isinstance(table.normalization, diff.NormalizationResult)
    assert table.normalization.setting == "none"
    assert set(table.normalization.shifts.values()) == {None}
    # Missing is None, never 0, and every None has its reason.
    assert all(v is None or v > 0 for v in table.columns["log2_intensity"])
    assert all(r["state"] in diff.OBSERVATION_STATES for r in table.records)
    assert all(r["state"] in ("quantified", "mbr_transferred") for r in table.records if r["log2_intensity"] is not None)
    assert table.records[0]["protein_groups"] == ["Q16643"] == [table.records[0]["unique_protein_group"]]


def test_the_two_bases_differ_only_where_a_transfer_is():
    kept = diff.Observations._from_wire(recorded("diff_observations_mm_mbr_kept.json")).records
    only = diff.Observations._from_wire(recorded("diff_observations_mm_msms_only.json")).records
    differing = [(a, b) for a, b in zip(kept, only) if a["log2_intensity"] != b["log2_intensity"]]
    assert len(differing) == 1
    a, b = differing[0]
    assert a["state"] == b["state"] == "mbr_transferred" and b["log2_intensity"] is None
    assert a["sample"] == "QE-002118_GM6_a"


def test_the_spike_in_normalization_takes_away_the_spike():
    plain = diff.Observations._from_wire(recorded("diff_observations_spike_in.json"))
    shifted = diff.Observations._from_wire(recorded("diff_observations_spike_in_normalized.json"))
    assert shifted.normalize is True and shifted.min_peptides == 50
    assert shifted.normalization.setting == "shared_peptide_median"
    assert shifted.normalization.reference_set_size == 86
    shifts = shifted.normalization.shifts
    assert shifts["spike_1x"] == pytest.approx(-shifts["spike_2x"])
    assert isinstance(shifted.responses_written, diff.WrittenResponses)
    assert shifted.responses_written.row_count == 243 and shifted.responses_written.sample_count == 2
    for a, b in zip(plain.records, shifted.records):
        if a["log2_intensity"] is None:
            assert b["log2_intensity"] is None
        else:
            assert b["log2_intensity"] == pytest.approx(a["log2_intensity"] - shifts[a["sample"]])


def test_the_example_inputs_are_vendored():
    for name in ("malat_dilution_samples.tsv", "design_factor_batch_covariate_samples.tsv",
                 "design_confounded_samples.tsv", "stratified_samples.tsv",
                 "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
                 "EcoliSpikeIn_QuantifiedPeptides.tsv"):
        assert (DIFF / name).is_file(), name
