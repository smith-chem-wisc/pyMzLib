"""Tests for :func:`pymzlib.proteins.align`, :func:`align_many` and :func:`parse_accessions`.

Offline. Every fixture is recorded from the real bridge (scripts/record_fixture.py):

* ``proteins_align_deletion.json``: mzLib's own deletion case (QRQ missing from the target).
* ``proteins_align_many.json``: that pair and mzLib's extra-N-terminus case, in one call.
* ``proteins_align_charged.json``: the extra N-terminus with ``--end-gaps charged``.
* ``proteins_align_no_map.json``: the two pairs with ``--no-residue-map``.
* ``proteins_align_dync1i2.json``: DYNC1I2 (Q13409) against its isoforms 2B and 2C, from the
  committed ``human_extra.fasta``.
* ``proteins_align_albumin.json``: human albumin (P02768) against bovine (P02769).
* ``proteins_parse_accession.json``: nine accessions, recognized and not.

What is pinned here is that pyMzLib sends only the options it was given, sends the sequences and
accessions unchanged, and projects mzLib's answer residue for residue without repairing it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pymzlib import _bridge, proteins

FIXTURES = Path(__file__).parent / "fixtures"

DELETION = ("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ")
OVERHANG = ("MLLLPEPTIDE", "PEPTIDE")
ACCESSIONS = [
    "P04637", "P04637-2", "Q13409-2", "NP_000537.3_R72P", "P04637_S70N_A80T", "P04637_2",
    "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768",
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
            seen["timeout"] = timeout
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        monkeypatch.setattr(_bridge, "require_verb", lambda verb, since: seen.setdefault("guard", (verb, since)))
        return seen

    return use


def _fasta(name: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    accession = None
    for line in (FIXTURES / "proteins" / name).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            accession = line.split("|")[1]
            sequences[accession] = ""
        elif accession:
            sequences[accession] += line.strip()
    return sequences


# ---- align ---------------------------------------------------------------------------------------


def test_align_sends_the_pair_on_stdin_and_no_option_it_was_not_given(replay):
    seen = replay("proteins_align_deletion.json")
    proteins.align(*DELETION)

    assert seen["args"] == ("proteins", "align")
    assert seen["stdin"] == "MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\n"
    assert seen["guard"] == ("proteins align", "0.6.0")
    assert seen["timeout"] is None


def test_align_projects_mzlibs_counts_and_rows(replay):
    replay("proteins_align_deletion.json")
    aln = proteins.align(*DELETION)

    assert aln.aligner_id == "global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"
    assert (aln.score, aln.identical, aln.substituted) == (80, 19, 0)   # 39 + 55 - 14, mzLib's figure
    assert (aln.source_residues_facing_gap, aln.target_residues_facing_gap) == (3, 0)
    assert (aln.source_length, aln.target_length) == (22, 19)
    assert aln.aligned_target == "MKTAYIAK---ISFVKSHFSRQ"
    assert len(aln.caveats) == 7


def test_a_residue_facing_a_gap_maps_to_none_never_a_neighbour(replay):
    replay("proteins_align_deletion.json")
    aln = proteins.align(*DELETION)

    assert [aln.map_source_to_target(p) for p in range(7, 14)] == [7, 8, None, None, None, 9, 10]
    assert aln.map_target_to_source(9) == 12
    assert len(aln.source_to_target) == 22 and len(aln.target_to_source) == 19
    assert all(aln.map_target_to_source(p) is not None for p in range(1, 20))


@pytest.mark.parametrize("position", [0, 23, -1, True, 1.0, "1"])
def test_a_position_outside_the_sequence_is_refused(replay, position):
    replay("proteins_align_deletion.json")
    aln = proteins.align(*DELETION)

    with pytest.raises(_bridge.UsageError):
        aln.map_source_to_target(position)


@pytest.mark.parametrize(
    ("kwargs", "sent"),
    [
        ({"gap_open": 5}, ("--gap-open", "5")),
        ({"gap_extend": 2}, ("--gap-extend", "2")),
        ({"max_cells": 1000}, ("--max-cells", "1000")),
        ({"end_gaps": "charged"}, ("--end-gaps", "charged")),
        ({"gap_open": 0, "end_gaps": "free"}, ("--gap-open", "0", "--end-gaps", "free")),
    ],
)
def test_each_option_given_is_sent_and_nothing_else(replay, kwargs, sent):
    seen = replay("proteins_align_deletion.json")
    proteins.align(*DELETION, **kwargs)

    assert seen["args"] == ("proteins", "align", *sent)


def test_charged_end_gaps_come_back_in_the_envelope(replay):
    replay("proteins_align_charged.json")
    aln = proteins.align(*OVERHANG, end_gaps="charged")

    assert aln.score == 39 - (11 + 4)
    assert aln.aligner_id.endswith("end-gaps=charged")


def test_align_many_keeps_input_order_and_one_alignment_per_pair(replay):
    seen = replay("proteins_align_many.json")
    alns = proteins.align_many([DELETION, list(OVERHANG)])

    assert seen["stdin"] == "MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\nMLLLPEPTIDE\tPEPTIDE\n"
    assert alns.record_count == 2
    assert (alns.matrix, alns.gap_open, alns.gap_extend, alns.end_gaps, alns.max_cells) == (
        "BLOSUM62", 11, 1, "free", 200_000_000)
    assert alns.columns["pair_index"] == [0, 1]
    assert alns.residue_map.row_count == 41 + 18
    pairs = alns.alignments
    assert [a.pair_index for a in pairs] == [0, 1]
    assert pairs[1].map_source_to_target(5) == 1                 # the free MLLL overhang
    assert [pairs[1].map_source_to_target(p) for p in range(1, 5)] == [None] * 4
    assert alns.records[1]["aligned_target"] == "----PEPTIDE"


def test_residue_map_false_sends_the_flag_and_leaves_no_map(replay):
    seen = replay("proteins_align_no_map.json")
    alns = proteins.align_many([DELETION, OVERHANG], residue_map=False)

    assert seen["args"] == ("proteins", "align", "--no-residue-map")
    assert alns.residue_map is None
    assert alns.columns["score"] == [80, 39]
    (first, _) = alns.alignments
    assert first.source_to_target is None
    with pytest.raises(_bridge.UsageError, match="no residue map"):
        first.map_source_to_target(1)


def test_real_isoforms_shift_every_position_after_the_missing_exon(replay):
    replay("proteins_align_dync1i2.json")
    seq = _fasta("human_extra.fasta")
    two_b, two_c = proteins.align_many([(seq["Q13409"], seq["Q13409-2"]), (seq["Q13409"], seq["Q13409-3"])]).alignments

    assert (two_b.identical, two_b.source_residues_facing_gap) == (632, 6)
    assert two_b.map_source_to_target(75) == 75
    assert [two_b.map_source_to_target(p) for p in range(76, 82)] == [None] * 6
    # V76 and V82 could each pair with the isoform's V76; mzLib's tie-break pairs V82.
    assert seq["Q13409"][75] == seq["Q13409"][81] == "V"
    assert two_b.map_source_to_target(82) == 76
    assert (two_c.identical, two_c.source_residues_facing_gap) == (612, 26)
    assert two_c.map_source_to_target(638) == 612


def test_real_orthologs_substituted_residues_keep_their_own_letters(replay):
    replay("proteins_align_albumin.json")
    seq = {**_fasta("human_extra.fasta"), **_fasta("contaminants.fasta")}
    aln = proteins.align(seq["P02768"], seq["P02769"])
    alns = proteins.align_many([(seq["P02768"], seq["P02769"])])
    rows = [r for r in alns.residue_map.records if r["sequence"] == "source"]

    assert (aln.identical, aln.substituted, aln.source_residues_facing_gap) == (465, 142, 2)
    assert aln.map_source_to_target(58) == 58                     # mature Cys34, conserved
    assert seq["P02768"][57] == seq["P02769"][57] == "C"
    assert sum(r["mapped_residue"] not in (None, r["residue"]) for r in rows) == aln.substituted


@pytest.mark.parametrize(
    ("pairs", "match"),
    [
        ("PEPTIDE", "not one string"),
        ([("PEPTIDE",)], "Pair 0 must be a"),
        ([("A", "B", "C")], "Pair 0 must be a"),
        (["PE"], "Pair 0 must be a"),
        ([("PEPTIDE", None)], "target must be a str"),
        ([("PEP\tTIDE", "PEPTIDE")], "tab or line break"),
        ([("PEPTIDE", "PEP\nTIDE")], "tab or line break"),
        ([("PEPTIDE", "PEPTIDE"), ("", " ")], "Pair 1 has two empty sequences"),
        ([], "At least one pair"),
    ],
)
def test_bad_pairs_are_refused_before_the_bridge(monkeypatch, pairs, match):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match=match):
        proteins.align_many(pairs)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"gap_open": True}, "gap_open must be an int"),
        ({"gap_extend": 1.0}, "gap_extend must be an int"),
        ({"max_cells": "1000"}, "max_cells must be an int"),
        ({"end_gaps": "Free"}, "end_gaps must be one of"),
        ({"end_gaps": True}, "end_gaps must be one of"),
    ],
)
def test_bad_options_are_refused_before_the_bridge(monkeypatch, kwargs, match):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match=match):
        proteins.align(*DELETION, **kwargs)


def test_an_empty_sequence_on_one_side_is_sent(replay):
    seen = replay("proteins_align_deletion.json")
    proteins.align_many([("", "PEP")])

    assert seen["stdin"] == "\tPEP\n"


def test_align_refuses_a_non_string(monkeypatch):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match="source must be a str"):
        proteins.align(["PEPTIDE"], "PEPTIDE")


# ---- parse_accessions ----------------------------------------------------------------------------


def test_parse_accessions_sends_each_accession_unchanged(replay):
    seen = replay("proteins_parse_accession.json")
    proteins.parse_accessions([" P04637", "P04637 ", "p04637", "P04637|P02768"])

    assert seen["args"] == ("proteins", "parse-accession")
    assert seen["stdin"] == " P04637\nP04637 \np04637\nP04637|P02768\n"   # not stripped, not split
    assert seen["guard"] == ("proteins parse-accession", "0.6.0")


def test_parse_accessions_projects_mzlibs_reading(replay):
    replay("proteins_parse_accession.json")
    parsed = proteins.parse_accessions(ACCESSIONS)
    rows = {r["input"]: r for r in parsed.records}

    assert parsed.record_count == 9
    assert parsed.namespace_counts == {"unrecognized": 4, "uniprot": 4, "refseq": 1}
    assert set(parsed.namespace_counts) == set(proteins.NAMESPACES)
    assert parsed.with_applied_variants_count == 2
    assert rows["Q13409-2"]["entry_accession"] == "Q13409" and rows["Q13409-2"]["isoform"] == 2
    refseq = rows["NP_000537.3_R72P"]
    assert (refseq["entry_accession"], refseq["version"], refseq["parent_accession"], refseq["applied_variants"]) == (
        "NP_000537", 3, "NP_000537.3", "R72P")
    assert rows["P04637_S70N_A80T"]["applied_variants"] == "S70N_A80T"
    assert rows["P04637_S70N_A80T"]["has_applied_variants"] is True
    for unrecognized in ("P04637_2", "DECOY_P02768", "sp|P04637|P53_HUMAN", "P04637|P02768"):
        row = rows[unrecognized]
        assert row["namespace"] == "unrecognized"
        assert row["entry_accession"] == unrecognized          # kept verbatim, never repaired
        assert row["applied_variants"] is None
    assert parsed.entry_of()["P04637-2"] == "P04637"
    assert len(parsed.caveats) == 7


@pytest.mark.parametrize(
    ("accessions", "match"),
    [
        ("P04637", "not one string"),
        ([], "At least one accession"),
        (["P04637", ""], "blank entry"),
        (["P04637", "   "], "blank entry"),
        (["P04637\nP02768"], "line break"),
        (["P04637\r"], "line break"),
        ([None], "must be a str"),
    ],
)
def test_bad_accession_lists_are_refused_before_the_bridge(monkeypatch, accessions, match):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError, match=match):
        proteins.parse_accessions(accessions)
