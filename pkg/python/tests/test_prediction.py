"""Tests for peptide-property prediction.

Offline apart from the canaries at the end. The payloads are **recorded from the real bridge
against the real Koina server**, so a wire-shape change shows up here as a parse failure rather
than as a fixture that agrees with a Python file and with nothing else. What is under test is the
Python layer's job: build the input table, validate arguments before spawning anything, and parse
the payload into typed objects.

Two of these tests exist because of traps rather than features. ``Constraint`` encodes a tri-state
mzLib expresses as a nullable set whose emptiness means the opposite of what it looks like, and the
fragment arrays are ragged, so a caller who indexes them as a rectangle gets a wrong answer rather
than an exception.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymzlib import _bridge, prediction

FIXTURES = Path(__file__).parent / "fixtures"


def payload(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def recorded(monkeypatch):
    """Serve a recorded payload, capturing the args and stdin the call would have sent."""
    seen: dict = {}

    def serve(name: str):
        data = payload(name)

        def fake_invoke(*args, **kwargs):
            seen["args"] = args
            seen["kwargs"] = kwargs
            return data

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        return seen

    return serve


@pytest.fixture()
def captured(monkeypatch):
    seen: dict = {}

    def fake_invoke(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return {}

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# ---- the catalogue ---------------------------------------------------------------------------


def test_every_model_mzlib_can_call_is_listed(recorded):
    recorded("prediction_models.json")

    catalogue = prediction.models()

    assert len(catalogue) == 37
    assert all(isinstance(m, prediction.Model) for m in catalogue)
    # Enumerated from mzLib, so a model added upstream appears without editing this package.
    assert all(m.error is None for m in catalogue), "every model must be constructible"


def test_the_five_families_are_all_represented(recorded):
    recorded("prediction_models.json")

    families = {m.family for m in prediction.models()}

    assert families == {
        "retention_time",
        "fragment_intensity",
        "collisional_cross_section",
        "detectability",
        "crosslink_intensity",
    }


def test_every_model_names_the_verb_that_calls_it(recorded):
    recorded("prediction_models.json")

    assert all(m.verb for m in prediction.models())


def test_the_retention_time_unit_is_per_model_not_per_family(recorded):
    recorded("prediction_models.json")
    catalogue = {m.model: m for m in prediction.models()}

    # The distinction that makes a bare IsIndexed boolean insufficient: one model in the family
    # returns absolute minutes and the rest return a dimensionless index.
    assert catalogue["Prosit_2019_irt"].retention_time_unit == "indexed_retention_time"
    assert catalogue["Chronologer_RT"].retention_time_unit == "minutes"
    # ...and it is meaningless outside that family.
    assert catalogue["IM2Deep"].retention_time_unit is None


def test_a_constraint_distinguishes_not_applicable_from_required_any(recorded):
    recorded("prediction_models.json")
    catalogue = {m.model: m for m in prediction.models()}

    # The trap this type exists for. mzLib expresses both as a nullable set: null means "this model
    # has no such input", empty means "required, any value". Reading the raw collection makes CID
    # look permissive and HCD look impossible, which is backwards for both.
    hcd = catalogue["Prosit_2020_intensity_HCD"].collision_energy
    cid = catalogue["Prosit_2020_intensity_CID"].collision_energy

    assert hcd.requirement == "any_value_required"
    assert hcd.applicable is True
    assert cid.requirement == "not_applicable"
    assert cid.applicable is False


def test_a_constraint_lists_its_values_when_it_has_them(recorded):
    recorded("prediction_models.json")
    catalogue = {m.model: m for m in prediction.models()}

    altimeter = catalogue["Altimeter_2024_intensities"].collision_energy
    assert altimeter.requirement == "one_of"
    assert min(altimeter.values) == 20 and max(altimeter.values) == 40
    assert altimeter.accepts(30) is True
    assert altimeter.accepts(45) is False

    unispec = catalogue["UniSpec"].instrument_type
    assert unispec.requirement == "one_of"
    assert "LUMOS" in unispec.values


def test_a_model_that_accepts_no_modifications_says_so(recorded):
    recorded("prediction_models.json")
    catalogue = {m.model: m for m in prediction.models()}

    # An empty allowed-UNIMOD list is a real answer, not a missing one.
    assert catalogue["pfly_2024_fine_tuned"].allowed_unimod_ids == []
    assert catalogue["pfly_2024_fine_tuned"].accepts_modifications is False
    assert catalogue["Prosit_2019_irt"].accepts_modifications is True


def test_models_sends_a_family_filter_when_asked(captured):
    prediction.models(" retention_time ")
    assert captured["args"] == ("predict", "models", "--family", "retention_time")


def test_models_does_not_send_an_empty_family(captured):
    prediction.models()
    assert "--family" not in captured["args"]


# ---- retention time --------------------------------------------------------------------------


def test_retention_time_parses_a_prediction_table(recorded):
    recorded("prediction_rt_irt.json")

    result = prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK", "ELVISLIVESK"])

    assert isinstance(result, prediction.Predictions)
    assert result.model == "Prosit_2019_irt"
    assert result.row_count == 2
    assert len(result.columns["retention_time"]) == 2


def test_an_irt_model_is_not_reported_as_minutes(recorded):
    recorded("prediction_rt_irt.json")
    result = prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK"])

    assert result.retention_time_unit == "indexed_retention_time"
    assert any("iRT" in caveat for caveat in result.caveats)


def test_chronologer_over_the_network_is_reported_as_minutes(recorded):
    recorded("prediction_rt_minutes.json")
    result = prediction.retention_time("Chronologer_RT", ["PEPTIDEK"])

    assert result.retention_time_unit == "minutes"
    # ...and says so, because the LOCAL Chronologer returns % acetonitrile from the same weights.
    assert any("acetonitrile" in caveat for caveat in result.caveats)


def test_retention_time_sends_the_peptides_on_stdin(captured):
    prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK", "ELVISLIVESK"])

    assert captured["args"] == ("predict", "retention-time", "--model", "Prosit_2019_irt")
    # On stdin rather than argv: argv has a hard ceiling of roughly 32 KB and a real run is
    # thousands of peptides.
    assert captured["kwargs"]["stdin"] == "sequence\nPEPTIDEK\nELVISLIVESK\n"


# ---- fragments -------------------------------------------------------------------------------


def test_fragment_arrays_are_ragged_and_that_is_correct(recorded):
    recorded("prediction_fragments.json")

    result = prediction.fragments(
        "Prosit_2020_intensity_HCD", ["PEPTIDEK", "ELVISLIVESK"],
        precursor_charge=2, collision_energy=28,
    )

    lengths = [len(a) for a in result.columns["fragment_mz"]]
    # Koina returns a fixed-width grid with -1 for impossible ions and mzLib drops those, so each
    # row is as long as ITS peptide's possible ions. Indexing these as a rectangle is wrong.
    assert len(set(lengths)) > 1, "two peptides of different lengths must give different arrays"
    for row in result.records:
        assert len(row["fragment_mz"]) == len(row["fragment_intensity"])
        assert len(row["fragment_mz"]) == len(row["fragment_annotations"])


def test_the_nominal_ion_count_is_not_the_array_length(recorded):
    recorded("prediction_fragments.json")
    result = prediction.fragments("Prosit_2020_intensity_HCD", ["PEPTIDEK"], precursor_charge=2)

    # The model's published count is 174; a short tryptic peptide gets a fraction of that.
    assert len(result.columns["fragment_mz"][0]) < 174


def test_intensities_are_declared_relative(recorded):
    recorded("prediction_fragments.json")
    result = prediction.fragments("Prosit_2020_intensity_HCD", ["PEPTIDEK"], precursor_charge=2)

    assert result.intensity_scale == "relative"
    assert any("RELATIVE" in caveat for caveat in result.caveats)


def test_a_peptide_that_cannot_be_predicted_still_gets_a_row(recorded):
    # Prosit_2020_intensity_HCD requires a collision energy; omitting it must not lose the row,
    # or predictions would no longer line up with the peptides that were sent.
    recorded("prediction_fragments_warned.json")

    result = prediction.fragments("Prosit_2020_intensity_HCD", ["PEPTIDEK"], precursor_charge=2)

    assert result.row_count == 1
    assert result.failed_row_count == 1
    assert result.columns["fragment_mz"][0] is None
    assert result.warnings and "CollisionEnergy" in result.warnings[0][1]


def test_fragments_sends_every_input_column(captured):
    prediction.fragments(
        "Prosit_2020_intensity_HCD", ["PEPTIDEK"],
        precursor_charge=2, collision_energy=28, instrument_type="QE", fragmentation_type="HCD",
    )

    stdin = captured["kwargs"]["stdin"]
    assert stdin.splitlines()[0] == (
        "sequence\tprecursor_charge\tcollision_energy\tinstrument_type\tfragmentation_type"
    )
    assert stdin.splitlines()[1] == "PEPTIDEK\t2\t28\tQE\tHCD"


def test_a_per_peptide_mapping_overrides_the_shared_default(captured):
    prediction.fragments(
        "Prosit_2020_intensity_HCD",
        ["PEPTIDEK", {"sequence": "ELVISLIVESK", "precursor_charge": 3}],
        precursor_charge=2,
    )

    rows = captured["kwargs"]["stdin"].splitlines()
    assert rows[1].startswith("PEPTIDEK\t2")
    assert rows[2].startswith("ELVISLIVESK\t3")


def test_an_absent_optional_column_is_an_empty_cell_not_the_string_none(captured):
    prediction.fragments("Prosit_2020_intensity_HCD", ["PEPTIDEK"], precursor_charge=2)

    # "None" would reach the bridge as a value and fail to parse as an integer, reporting a
    # confusing error about a column the caller never set.
    assert "None" not in captured["kwargs"]["stdin"]
    assert captured["kwargs"]["stdin"].splitlines()[1] == "PEPTIDEK\t2\t\t\t"


# ---- ccs and detectability -------------------------------------------------------------------


def test_ccs_is_in_square_angstroms_not_reduced_mobility(recorded):
    recorded("prediction_ccs.json")
    result = prediction.ccs("IM2Deep", ["PEPTIDEK"], precursor_charge=2)

    assert result.collisional_cross_section_unit == "square_angstroms"
    assert any("1/K0" in caveat for caveat in result.caveats)


def test_detectability_returns_four_classes_that_sum_to_one(recorded):
    recorded("prediction_detectability.json")
    result = prediction.detectability("pfly_2024_fine_tuned", ["PEPTIDEK"])

    classes = [
        "not_detectable", "low_detectability",
        "intermediate_detectability", "high_detectability",
    ]
    assert all(name in result.columns for name in classes)
    assert abs(sum(result.columns[name][0] for name in classes) - 1.0) < 1e-5


# ---- validation ------------------------------------------------------------------------------

VERBS = [
    ("retention_time", ["PEPTIDEK"]),
    ("fragments", ["PEPTIDEK"]),
    ("ccs", ["PEPTIDEK"]),
    ("detectability", ["PEPTIDEK"]),
    ("crosslink_fragments", [{"alpha_sequence": "PEPTIDEK"}]),
]


@pytest.mark.parametrize(("name", "peptides"), VERBS)
@pytest.mark.parametrize("bad", ["", "   ", None, 7])
def test_every_verb_rejects_a_blank_model_before_spawning_anything(name, peptides, bad, captured):
    with pytest.raises(_bridge.UsageError, match="model name is required"):
        getattr(prediction, name)(bad, peptides)

    assert "args" not in captured


@pytest.mark.parametrize(("name", "peptides"), VERBS)
def test_every_verb_rejects_an_empty_peptide_list(name, peptides, captured):
    with pytest.raises(_bridge.UsageError, match="At least one peptide"):
        getattr(prediction, name)("some_model", [])

    assert "args" not in captured


@pytest.mark.parametrize(("name", "peptides"), VERBS)
@pytest.mark.parametrize("bad", [0, -1, "5", 2.5, True])
def test_every_verb_rejects_a_bad_max_batches(name, peptides, bad, captured):
    with pytest.raises(_bridge.UsageError, match="max_batches"):
        getattr(prediction, name)("some_model", peptides, max_batches=bad)

    assert "args" not in captured


@pytest.mark.parametrize(("name", "peptides"), VERBS)
def test_every_verb_leaves_the_politeness_defaults_alone(name, peptides, captured):
    getattr(prediction, name)("some_model", peptides)

    # Koina is a shared community server. A binding that raised the throughput defaults would be
    # spending someone else's GPU time without being asked.
    assert "--max-batches" not in captured["args"]
    assert "--throttle-ms" not in captured["args"]


def test_a_peptide_that_is_neither_a_string_nor_a_mapping_is_refused(captured):
    with pytest.raises(_bridge.UsageError, match="sequence string or a mapping"):
        prediction.retention_time("Prosit_2019_irt", [42])

    assert "args" not in captured


def test_the_module_exports_what_it_documents():
    for name, _ in VERBS:
        assert name in prediction.__all__
    assert all(hasattr(prediction, name) for name in prediction.__all__)
    defined_here = {
        name
        for name, value in vars(prediction).items()
        if not name.startswith("_") and getattr(value, "__module__", None) == prediction.__name__
    }
    assert defined_here - set(prediction.__all__) == set()


# ---- live canaries ---------------------------------------------------------------------------


@pytest.mark.network
def test_the_koina_catalogue_still_matches_what_this_package_describes():
    """The recorded catalogue must not drift from the installed mzLib."""
    from .conftest import external_service

    with external_service("Koina"):
        live = {m.model for m in prediction.models()}

    recorded = {m["model"] for m in payload("prediction_models.json")["models"] if "model" in m}
    assert live == recorded, (
        "the recorded model catalogue has drifted from mzLib; re-record the fixture"
    )


@pytest.mark.network
def test_koina_still_answers_and_the_fields_this_package_reads_are_populated():
    from .conftest import external_service

    with external_service("Koina"):
        result = prediction.retention_time("Prosit_2019_irt", ["PEPTIDEK", "ELVISLIVESK"])

    assert result.row_count == 2
    assert result.failed_row_count == 0
    assert all(value is not None for value in result.columns["retention_time"])
    assert result.retention_time_unit == "indexed_retention_time"
