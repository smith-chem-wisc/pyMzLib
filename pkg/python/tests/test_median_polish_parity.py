"""End-to-end parity: pyMzLib's median-polish roll-up must reproduce FlashLFQ's own protein output.

Unlike the rest of the median-polish tests, this one runs the **real bridge** — the whole point is
to prove the two halves agree with FlashLFQ end to end, not just that the Python layer shuttles
arguments correctly. It skips when no bridge is built/staged (a source checkout that has not run
``publish-bridge.ps1`` and has not set ``PYMZLIB_BRIDGE``), the same way the C# suite ignores its
real-data test when the mzLib test data is absent.

The fixtures are a curated slice of a real FlashLFQ run (the E. coli spike-in vignette,
``UseSharedPeptidesForProteinQuant = false``, no experimental design): its ``QuantifiedPeptides.tsv``
rows for 40 proteins, paired with those proteins' rows from the ``QuantifiedProteins.tsv`` FlashLFQ
wrote in the same run. Because median polish quantifies each protein independently from its own
*unique* peptides, this slice reproduces the full run's numbers exactly — so the expected values here
are genuinely FlashLFQ's, computed by FlashLFQ, not a re-derivation. The slice deliberately spans the
interesting cases: proteins quantified in both runs, proteins measured in only one run (a ``0`` in the
other), and proteins FlashLFQ could not resolve at all (``NaN`` -> ``None``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pymzlib import _bridge, flashlfq

FIXTURES = Path(__file__).parent / "fixtures"
PEPTIDES = FIXTURES / "flashlfq_vignette_peptides.tsv"
PROTEINS = FIXTURES / "flashlfq_vignette_proteins.tsv"

# The two runs the vignette quantified. With no experimental design, median_polish keys each
# protein's intensities by run base name — matching the Intensity_<run> columns FlashLFQ wrote.
RUN1 = "09-04-18_EcoliSpikeInSingleShot1x"
RUN2 = "09-04-18_EcoliSpikeInSingleShot2x"


@pytest.fixture(scope="module")
def bridge():
    """Require a real bridge that supports the verb, or skip — this test cannot be run mocked.

    Two ways to skip rather than fail: no bridge is staged/built at all, or the bridge that *is*
    staged predates the ``quant median-polish`` verb. The latter is a real local-dev hazard — an old
    staged payload lingers from a previous ``publish-bridge.ps1`` — so it is detected by probing the
    verb with a path that does not exist: a current bridge answers "file not found" (the verb ran),
    an old one answers "Unknown command" (the verb is absent).
    """
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError:
        pytest.skip(
            "No mzLib bridge is built or staged; set PYMZLIB_BRIDGE to a bridge executable to run "
            "the FlashLFQ parity test."
        )

    try:
        flashlfq.median_polish("this-path-does-not-exist.tsv")
    except _bridge.UsageError as exc:
        if "Unknown command" in str(exc):
            pytest.skip(
                "The available mzLib bridge predates the 'quant median-polish' verb; rebuild and "
                "stage it (publish-bridge.ps1) to run the FlashLFQ parity test."
            )
        # Otherwise it's the expected "file not found" — the verb exists, so carry on.


def _expected_proteins() -> dict[str, tuple[float | None, float | None]]:
    """Parse FlashLFQ's QuantifiedProteins.tsv fixture into {name: (run1, run2)}, NaN -> None."""
    lines = PROTEINS.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    i1 = header.index("Intensity_" + RUN1)
    i2 = header.index("Intensity_" + RUN2)

    def cell(value: str) -> float | None:
        return None if value == "NaN" else float(value)

    out: dict[str, tuple[float | None, float | None]] = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        out[fields[0]] = (cell(fields[i1]), cell(fields[i2]))
    return out


def test_fixture_covers_the_interesting_cases():
    """Guard the fixture itself: if a regeneration flattened it, the parity test below is hollow."""
    expected = _expected_proteins()
    assert len(expected) == 40
    both = [v for v in expected.values() if v[0] not in (None, 0.0) and v[1] not in (None, 0.0)]
    one_zero = [v for v in expected.values() if 0.0 in v]
    unresolvable = [v for v in expected.values() if None in v]
    assert both, "expected some proteins quantified in both runs"
    assert one_zero, "expected some proteins measured in only one run"
    assert unresolvable, "expected some NaN/None proteins, to exercise the unquantifiable path"


def test_median_polish_reproduces_flashlfq_quantified_proteins(bridge):
    proteins = flashlfq.median_polish(str(PEPTIDES))
    got = {p.protein_group: p for p in proteins}
    expected = _expected_proteins()

    # Same proteins, no more and no fewer — the unique-peptide slice must not conjure phantom groups.
    assert set(got) == set(expected)

    for name, (exp1, exp2) in expected.items():
        for run, exp in ((RUN1, exp1), (RUN2, exp2)):
            actual = got[name].intensity(run)
            if exp is None:
                # FlashLFQ marked this sample unquantifiable (NaN); pyMzLib must too, as None.
                assert actual is None, f"{name}/{run}: FlashLFQ has NaN, pyMzLib has {actual!r}"
            else:
                assert actual is not None, f"{name}/{run}: FlashLFQ has {exp}, pyMzLib has None"
                # Six significant figures: the algorithm is identical, so the only gap is the last
                # few ULPs of the log2/pow2 round-trip. A looser match would hide a real divergence.
                assert actual == pytest.approx(exp, rel=1e-6), f"{name}/{run}: {actual} != {exp}"
