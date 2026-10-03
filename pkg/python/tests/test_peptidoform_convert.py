"""Tests for :func:`pymzlib.peptidoform.convert`, the projection of mzLib's sequence conversion.

Offline. Every fixture is recorded from the real bridge (scripts/record_fixture.py):

* ``peptidoform_convert_unimod.json``: four probe sequences to Unimod under ReturnNull - two UniProt
  modifications, one MetaMorpheus one, and one mzLib cannot map.
* ``peptidoform_convert_psmtsv.json`` / ``..._proforma.json``: the full sequences of mzLib's own
  ``BottomUpExample.psmtsv`` to Unimod and to ProForma.

The accessions are mzLib's. What is pinned here is that pyMzLib sends the call mzLib needs and
projects mzLib's answer row for row, without repairing it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymzlib import _bridge, peptidoform

FIXTURES = Path(__file__).parent / "fixtures"

PROBE = [
    "[UniProt:N-acetylserine on S]SEQK",
    "PEPK[UniProt:N6,N6-dimethyllysine on K]R",
    "PEPM[Common Variable:Oxidation on M]K",
    "PEPK[Made Up:Not a modification on K]R",
]


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, stdin=None, timeout=None):
            seen["args"] = args
            seen["stdin"] = stdin
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: seen.setdefault("guard", (verb, since)))
        return seen

    return use


def test_the_probe_cases_map_to_their_unimod_accessions(replay):
    seen = replay("peptidoform_convert_unimod.json")
    result = peptidoform.convert(PROBE)

    assert seen["args"] == (
        "peptidoform", "convert", "--from", "mzLib", "--to", "Unimod",
        "--mode", "ReturnNull", "--threads", "1",
    )
    assert seen["stdin"] == "\n".join(PROBE) + "\n"
    assert seen["guard"] == ("peptidoform convert", "0.4.0")

    rows = result.sequences
    assert [r.input for r in rows] == PROBE
    assert rows[0].output == "[UNIMOD:1]SEQK"        # N-acetylserine
    assert rows[1].output == "PEPK[UNIMOD:36]R"      # N6,N6-dimethyllysine, not Ethyl (UNIMOD:280)
    assert rows[2].output == "PEPM[UNIMOD:35]K"      # Oxidation on M
    assert [r.status for r in rows[:3]] == ["converted"] * 3
    assert all(r.ok for r in rows[:3])


def test_an_unmappable_modification_fails_its_row_and_names_itself(replay):
    replay("peptidoform_convert_unimod.json")
    result = peptidoform.convert(PROBE)

    (bad,) = result.not_converted
    assert bad.input == "PEPK[Made Up:Not a modification on K]R"
    assert bad.output is None
    assert bad.status == "failed"
    assert not bad.ok
    assert bad.incompatible_items == ["Made Up:Not a modification on K @3(K)"]
    # mzLib's Unimod serializer under ReturnNull records no reason code; projected as None, not filled in.
    assert bad.failure_reason is None
    assert result.outputs == ["[UNIMOD:1]SEQK", "PEPK[UNIMOD:36]R", "PEPM[UNIMOD:35]K", None]
    assert (result.record_count, result.converted_count, result.warned_count, result.failed_count) == (4, 3, 0, 1)


def test_the_envelope_names_the_formats_mzlib_registered(replay):
    replay("peptidoform_convert_unimod.json")
    result = peptidoform.convert(PROBE)

    assert (result.source_format, result.target_format, result.mode) == ("mzLib", "Unimod", "ReturnNull")
    assert result.source_formats == ["MassShift", "Modomics", "ProForma", "mzLib"]
    assert result.target_formats == ["Chronologer", "Essential", "MassShift", "ProForma", "Unimod", "mzLib"]
    assert result.column_names == [
        "input", "output", "status", "failure_reason", "incompatible_items", "warnings", "errors"
    ]
    assert len(result.caveats) == 6
    assert any("ProForma" in c and "UniProt" in c for c in result.caveats)


def test_real_psmtsv_sequences_to_unimod(replay):
    replay("peptidoform_convert_psmtsv.json")
    psms = _load("readers_results_psmtsv.json")
    result = peptidoform.convert(psms["columns"]["full_sequence"])

    assert [r.input for r in result.sequences] == psms["columns"]["full_sequence"]
    assert result.outputs[0] == "YPIEH[UNIMOD:34]GIVTNWDDMEK"
    assert result.outputs[4] == "YPIEH[UNIMOD:34]GIVTNWDDM[UNIMOD:35]EK"
    assert result.converted_count == result.record_count == 8
    assert not any("UniProt:" in o for o in result.outputs)


def test_proforma_leaves_uniprot_modifications_unresolved(replay):
    """The mzLib gap the caveats describe, pinned so a fix upstream is noticed when re-recorded."""
    seen = replay("peptidoform_convert_psmtsv_proforma.json")
    psms = _load("readers_results_psmtsv.json")
    result = peptidoform.convert(psms["columns"]["full_sequence"], target="ProForma")

    assert seen["args"][5] == "ProForma"
    assert result.outputs[0] == "YPIEH[UniProt:Tele-methylhistidine on H]GIVTNWDDMEK"
    assert result.outputs[2] == "AYHEQLSVAEITNAC[UNIMOD:4]FEPANQMVK"
    assert result.sequences[0].status == "converted"


def test_mode_and_threads_are_passed_through(replay):
    seen = replay("peptidoform_convert_unimod.json")
    peptidoform.convert(PROBE, mode="removeincompatibleelements", threads=-1)
    assert seen["args"][-4:] == ("--mode", "RemoveIncompatibleElements", "--threads", "-1")


@pytest.mark.parametrize(
    "kwargs, sequences, message",
    [
        ({}, "PEPTIDE", "not one string"),
        ({}, [], "At least one sequence"),
        ({}, ["PEPTIDE", "  "], "blank entry"),
        ({}, ["PEP\nTIDE"], "line break"),
        ({}, ["PEPTIDE", 7], "must be a str"),
        ({"mode": "Strict"}, ["PEPTIDE"], "mode must be one of"),
        ({"mode": None}, ["PEPTIDE"], "mode must be one of"),
        ({"source": " "}, ["PEPTIDE"], "source must be"),
        ({"target": None}, ["PEPTIDE"], "target must be"),
        ({"threads": 0}, ["PEPTIDE"], "threads"),
        ({"threads": True}, ["PEPTIDE"], "threads"),
    ],
)
def test_bad_arguments_are_refused_before_the_bridge(monkeypatch, kwargs, sequences, message):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match=message):
        peptidoform.convert(sequences, **kwargs)


def test_an_older_bridge_without_new_fields_still_projects():
    result = peptidoform.SequenceConversions._from_wire({"columns": {"input": ["A"], "output": [None]}})
    assert result.record_count == 0 and result.caveats == [] and result.source_formats == []
    assert result.outputs == [None]
