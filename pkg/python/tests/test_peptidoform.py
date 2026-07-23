"""Tests for the Peptidoform workflow.

The offline tier runs against a recorded payload and covers the Python layer's own behaviour.
The live tier is where the science is checked — including a closure invariant that would have
caught, in seconds, an 18 Da error a fake user nearly published.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import pymzlib
from pymzlib import _bridge, peptidoform

FIXTURE = Path(__file__).parent / "fixtures" / "peptidoform_P02768_small.json"


@pytest.fixture()
def recorded_digest(monkeypatch):
    """Serve a recorded albumin digest instead of calling the bridge."""
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: payload)
    return payload


@pytest.fixture()
def captured(monkeypatch):
    """Capture the argument list without running anything."""
    seen = {}
    payload = json.loads(FIXTURE.read_text(encoding="utf-8"))

    def fake_invoke(*args, timeout=None, **kwargs):
        seen["args"] = list(args)
        return payload

    monkeypatch.setattr(_bridge, "invoke", fake_invoke)
    return seen


# --------------------------------------------------------------------------- parsing


def test_digest_and_peptides_parse(recorded_digest):
    digest = peptidoform.fragments("P02768")
    assert digest.accession == "P02768"
    assert digest.full_name
    assert len(digest.peptides) == len(recorded_digest["peptides"])
    assert all(isinstance(p, peptidoform.Peptide) for p in digest.peptides)


def test_fragments_parse_into_typed_objects(recorded_digest):
    peptide = peptidoform.fragments("P02768").peptides[0]
    assert peptide.fragments
    first = peptide.fragments[0]
    assert isinstance(first, peptidoform.Fragment)
    assert first.product_type
    assert first.neutral_mass > 0


def test_fragment_count_sums_across_peptides(recorded_digest):
    digest = peptidoform.fragments("P02768")
    assert digest.fragment_count == sum(len(p.fragments) for p in digest.peptides)


# --------------------------------------------------------------------------- the m/z conversion


def test_mz_uses_the_proton_mass_not_hydrogen(recorded_digest):
    """0.55 mDa apart. At 500 m/z that is 1.1 ppm — on an Orbitrap, match versus miss. Libraries
    differ on this and almost never say which they used."""
    peptide = peptidoform.fragments("P02768").peptides[0]
    expected = (peptide.monoisotopic_mass + 2 * 1.00727646677) / 2
    assert peptide.mz(2) == pytest.approx(expected, abs=1e-9)

    hydrogen_answer = (peptide.monoisotopic_mass + 2 * 1.00782503) / 2
    assert abs(peptide.mz(2) - hydrogen_answer) > 1e-4, "must not silently use the hydrogen mass"


def test_mz_does_not_double_count_a_fixed_charge():
    """Trimethylation of a lysine ε-amine gives a quaternary ammonium. UniProt records the delta
    as 43.054227 — C3H7 minus an electron — so the mass already carries one charge, and only
    (z − fixed) protons may be added. Adding a full complement put a 2+ trimethylated peptide
    0.504 Th high, on the most important histone modification there is."""
    charged = peptidoform.Peptide(
        base_sequence="SAPATGGVK", full_sequence="", monoisotopic_mass=829.477794,
        length=9, one_based_start=1, one_based_end=9, missed_cleavages=0, fixed_charges=1,
    )
    assert charged.mz(2) == pytest.approx(
        (829.477794 + 1 * peptidoform.PROTON_MASS) / 2, abs=1e-9
    )
    naive = (829.477794 + 2 * peptidoform.PROTON_MASS) / 2
    assert abs(charged.mz(2) - naive) == pytest.approx(peptidoform.PROTON_MASS / 2, abs=1e-9)


def test_a_fixed_charge_peptide_is_observable_without_protonation():
    """It already carries the charge, so 1+ needs nothing added — an answer the old arithmetic
    could not produce."""
    charged = peptidoform.Peptide(
        base_sequence="K", full_sequence="", monoisotopic_mass=829.477794,
        length=1, one_based_start=1, one_based_end=1, missed_cleavages=0, fixed_charges=1,
    )
    assert charged.mz(1) == pytest.approx(829.477794, abs=1e-9)


def test_a_charge_below_the_fixed_charge_is_refused():
    """That species cannot exist, and silently returning a number for it would be worse."""
    charged = peptidoform.Peptide(
        base_sequence="K", full_sequence="", monoisotopic_mass=900.0,
        length=1, one_based_start=1, one_based_end=1, missed_cleavages=0, fixed_charges=2,
    )
    with pytest.raises(pymzlib.UsageError, match="fixed charge"):
        charged.mz(1)


@pytest.mark.parametrize("charge", [0, -1, 1.5, "2", True, None])
def test_mz_rejects_charges_that_are_not_positive_whole_numbers(recorded_digest, charge):
    peptide = peptidoform.fragments("P02768").peptides[0]
    with pytest.raises(pymzlib.UsageError):
        peptide.mz(charge)


# --------------------------------------------------------------------------- the census
#
# The point of the census is that a correct number produced by an invisible rule is still a trap.


def test_census_separates_sites_from_modifications(recorded_digest):
    census = peptidoform.fragments("P02768").modification_census
    assert census.annotated == 38
    assert census.applied == 14
    assert census.excluded == 24


def test_modification_positions_index_the_peptides_own_residues(recorded_digest):
    """mzLib's dictionary reserves slot 1 for the N-terminus, so its keys are one past the residue
    they modify. Passing that through unchanged pointed 474 of 498 modifications at the wrong
    residue — and at a residue that need not exist: a peptide MAR reported position 4 for its
    arginine."""
    checked = 0
    for peptide in proteoform_peptides(recorded_digest):
        for mod in peptide.modifications:
            residue = mod.get("one_based_residue")
            if residue is None:
                assert mod.get("terminus") in {"N", "C"}, "a non-residue position must name a terminus"
                continue
            assert 1 <= residue <= len(peptide.base_sequence)
            target = mod["id"].rsplit(" on ", 1)[-1]
            assert peptide.base_sequence[residue - 1] == target, (
                f"{peptide.base_sequence}: {mod['id']} points at "
                f"{peptide.base_sequence[residue - 1]} at residue {residue}"
            )
            checked += 1

    # Guard against the test passing vacuously: if the fixture ever loses its modified peptide,
    # the loops above run zero times and assert nothing. The whole point is to check a real
    # position, so require at least one residue-anchored modification to have been verified.
    assert checked >= 1, "fixture has no residue-anchored modification to check the position of"


def proteoform_peptides(_payload):
    """The parsed peptides from the recorded digest."""
    return [p for p in peptidoform.fragments("P02768").peptides if p.modifications]


def test_unresolved_modification_names_are_reported(recorded_digest):
    """A name UniProt annotates but cannot define vanishes silently otherwise — seven
    N6-lactoyllysine sites disappeared from histone H3.1 while the type summary still said
    'modified residue … loaded'."""
    census = peptidoform.ModificationCensus(
        sites=28, applied=114, annotated=124, by_type=[], unresolved=["N6-lactoyllysine"]
    )
    assert "N6-lactoyllysine" in census.explain()
    assert "resolved" in census.explain()


def test_the_two_exclusion_reasons_are_not_conflated():
    """Being excluded by type and failing to resolve are different failures with different fixes,
    and an explanation that merges them tells the reader the wrong thing about both."""
    census = peptidoform.ModificationCensus(
        sites=3, applied=5, annotated=9,
        by_type=[{"type": "cross-link", "count": 3, "loaded": False}],
        unresolved=["N6-lactoyllysine"],
    )
    text = census.explain()
    assert "Excluded by type" in text
    assert "Could not be resolved" in text


def test_census_explains_what_was_excluded_and_why(recorded_digest):
    text = peptidoform.fragments("P02768").modification_census.explain()
    assert "24" in text and "glycosylation site" in text
    assert "no defined chemical composition" in text


def test_census_says_so_when_nothing_was_excluded():
    census = peptidoform.ModificationCensus(sites=5, applied=9, annotated=9, by_type=[])
    assert census.excluded == 0
    assert "All 9" in census.explain()


def test_sites_and_modifications_are_not_the_same_number():
    """A histone lists K9me1/me2/me3/ac at one residue: four modifications, one site. Reporting
    the site count as a modification count made H3.1 look as though 93 annotations had been
    dropped when every one had been loaded."""
    census = peptidoform.ModificationCensus(sites=28, applied=114, annotated=124, by_type=[])
    assert census.sites < census.applied
    assert census.excluded == 10


# --------------------------------------------------------------------------- the silent cap


def test_truncation_is_visible():
    """A truncated answer and a short answer are indistinguishable from the outside."""
    assert peptidoform.Digest(
        accession="X", name="", full_name="", organism="", sequence_length=1,
        protease="trypsin|P", dissociation="ETD", terminus="Both", modifications_applied=True,
        max_modifications=4, max_isoforms=1024, peptides_at_cap=5,
        modification_census=peptidoform.ModificationCensus(0, 0, 0), peptides=[],
    ).truncated is True


def test_no_truncation_when_nothing_hit_the_cap(recorded_digest):
    assert peptidoform.fragments("P02768").truncated is False


# --------------------------------------------------------------------------- the call itself


def test_defaults_are_the_labs_opinion_and_reach_the_bridge(captured):
    peptidoform.fragments("P02768")
    args = captured["args"]
    assert args[:2] == ["peptidoform", "fragments"]
    assert args[args.index("--protease") + 1] == "trypsin|P", "the Keil rule, not mzLib's plain trypsin"
    assert args[args.index("--dissociation") + 1] == "ETD"
    assert args[args.index("--max-isoforms") + 1] == "1024"
    assert "--no-modifications" not in args


def test_every_default_is_reachable(captured):
    peptidoform.fragments(
        "P02768", protease="trypsin", dissociation="HCD", modifications=False,
        missed_cleavages=3, min_length=5, max_length=40, max_modifications=1, max_isoforms=9,
    )
    args = captured["args"]
    assert args[args.index("--protease") + 1] == "trypsin"
    assert args[args.index("--dissociation") + 1] == "HCD"
    assert args[args.index("--max-length") + 1] == "40"
    assert args[args.index("--max-isoforms") + 1] == "9"
    assert "--no-modifications" in args


@pytest.mark.parametrize("accession", ["", "   ", None, 123])
def test_blank_accession_is_refused_before_any_work(accession):
    with pytest.raises(pymzlib.UsageError):
        peptidoform.fragments(accession)


@pytest.mark.parametrize("bad", [-1, 1.5, "2", True])
def test_negative_or_non_integer_counts_are_refused(bad):
    with pytest.raises(pymzlib.UsageError):
        peptidoform.fragments("P02768", missed_cleavages=bad)
