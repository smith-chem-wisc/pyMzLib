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

FIXTURE = Path(__file__).parent / "fixtures" / "Peptidoform_P02768_small.json"


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
