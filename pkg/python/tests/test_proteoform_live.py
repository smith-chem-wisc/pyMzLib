"""Live canaries for the proteoform workflow, against the real UniProt.

These are the tests that would notice UniProt changing its XML, mzLib changing its digestion, or
the modification definitions failing to resolve. Each routes through ``external_service()`` so an
outage skips with an explanatory message rather than failing.

The closure test below is the most valuable thing in this file. It came from a fake user who,
using a different library, produced a fragment table that was silently 18 Da wrong and looked
entirely plausible — and caught it only because they had added this check on a hunch.
"""

from __future__ import annotations

import pytest

import pymzlib
from conftest import external_service

pytestmark = pytest.mark.network

#: Serum albumin: heavily annotated, and most of its annotations are the *excluded* kind.
ALBUMIN = "P02768"
#: Histone H3.1: the case where modifications are the biology, and combinatorics bite.
HISTONE = "P68431"

#: Mass of a proton, and of a hydrogen atom. The closure identity needs the latter.
HYDROGEN_MASS = 1.00782503207


def test_the_workflow_still_answers_end_to_end():
    with external_service("UniProt"):
        digest = pymzlib.proteoform.fragments(ALBUMIN, max_modifications=1)

    assert digest.accession == ALBUMIN
    assert "lbumin" in digest.full_name
    assert digest.peptides and digest.fragment_count > 0


def test_etd_produces_c_and_z_ions_not_b_and_y():
    """The dissociation type must reach mzLib, or the caller silently gets the wrong chemistry."""
    with external_service("UniProt"):
        digest = pymzlib.proteoform.fragments(ALBUMIN, max_modifications=0, min_length=20)

    kinds = {f.product_type for p in digest.peptides for f in p.fragments}
    assert any(k.startswith("c") for k in kinds), f"expected c ions from ETD, got {sorted(kinds)}"
    assert any(k.lower().startswith("z") for k in kinds), f"expected z ions from ETD, got {sorted(kinds)}"


def test_fragment_series_close_on_the_precursor_mass():
    """c_i + z•_(n-i) must equal the peptide mass plus one hydrogen, for every i.

    This is the invariant that catches a whole class of silent error: a wrong terminal group, a
    wrong ion definition, a modification applied to the wrong terminus. Any of those produce a
    fragment table with sensible spacings and monotonic series that is nonetheless wrong.
    """
    with external_service("UniProt"):
        digest = pymzlib.proteoform.fragments(ALBUMIN, max_modifications=0, min_length=15)

    checked = 0
    for peptide in digest.peptides[:5]:
        c_ions = {f.fragment_number: f.neutral_mass for f in peptide.fragments
                  if f.product_type == "c" and f.neutral_loss == 0}
        z_ions = {f.fragment_number: f.neutral_mass for f in peptide.fragments
                  if f.product_type.lower().startswith("z") and f.neutral_loss == 0}
        n = peptide.length
        for i, c_mass in c_ions.items():
            z_mass = z_ions.get(n - i)
            if z_mass is None:
                continue
            closure = c_mass + z_mass - peptide.monoisotopic_mass
            assert closure == pytest.approx(HYDROGEN_MASS, abs=5e-4), (
                f"{peptide.base_sequence}: c{i} + z{n - i} - M = {closure:.6f}, "
                f"expected {HYDROGEN_MASS:.6f}"
            )
            checked += 1

    assert checked > 10, f"expected many closure pairs to check, got {checked}"


def test_the_annotation_census_reports_what_was_excluded():
    """Albumin's annotations are mostly glycosylation sites, which have no defined mass."""
    with external_service("UniProt"):
        census = pymzlib.proteoform.fragments(ALBUMIN, max_modifications=0).modification_census

    assert census.annotated > census.applied
    assert census.excluded > 0
    assert "glycosylation" in census.explain()


def test_modifications_change_the_answer_substantially():
    """The control that shows the annotations are doing real work."""
    with external_service("UniProt"):
        with_mods = pymzlib.proteoform.fragments(ALBUMIN, max_modifications=1)
        without = pymzlib.proteoform.fragments(ALBUMIN, modifications=False)

    assert len(with_mods.peptides) > len(without.peptides)
    assert with_mods.modified_peptides
    assert not without.modified_peptides


@pytest.mark.slow
def test_modification_isoforms_are_enumerated_combinatorially():
    """Histones are where this matters: alternatives at one residue multiply across residues."""
    with external_service("UniProt"):
        one = pymzlib.proteoform.fragments(HISTONE, max_modifications=1)
        two = pymzlib.proteoform.fragments(HISTONE, max_modifications=2)

    assert len(two.peptides) > len(one.peptides) * 2, (
        "modification isoforms should multiply, not add"
    )


@pytest.mark.slow
def test_the_isoform_cap_truncates_and_says_so():
    """mzLib's default of 1024 isoforms per peptide truncates silently. It must not be silent
    here: a truncated proteoform list is indistinguishable from a short one."""
    with external_service("UniProt"):
        capped = pymzlib.proteoform.fragments(HISTONE, max_modifications=4, max_isoforms=1024)
        raised = pymzlib.proteoform.fragments(HISTONE, max_modifications=4, max_isoforms=100_000)

    assert capped.truncated, "the default cap binds on a histone at four modifications"
    assert capped.peptides_at_cap > 0
    assert not raised.truncated
    assert len(raised.peptides) > len(capped.peptides), (
        "raising the cap must recover proteoforms the default discarded"
    )


def test_an_unknown_accession_is_a_usage_error_not_an_empty_result():
    with external_service("UniProt"):
        with pytest.raises(pymzlib.UsageError):
            pymzlib.proteoform.fragments("P99999999")


def test_an_unknown_protease_names_the_alternatives():
    with external_service("UniProt"):
        with pytest.raises(pymzlib.UsageError, match="Unknown protease"):
            pymzlib.proteoform.fragments(ALBUMIN, protease="banana")
