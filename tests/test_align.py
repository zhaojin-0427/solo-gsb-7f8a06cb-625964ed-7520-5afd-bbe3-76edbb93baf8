"""Tests for semi-global read alignment and the /align API."""

from __future__ import annotations

import random

from fastapi.testclient import TestClient

from app.aligner import (
    Alignment,
    _best_on_strand,
    align_read,
    align_reads,
    reverse_complement,
)
from app.main import app

client = TestClient(app)


def _key(a: Alignment | None):
    assert a is not None
    return (
        a.edits,
        a.cigar.count("I") + a.cigar.count("D"),
        a.start,
        a.end,
        0 if a.strand == "+" else 1,
        a.cigar,
    )


# ---------------------------------------------------------------------------
# reverse complement
# ---------------------------------------------------------------------------


def test_reverse_complement():
    assert reverse_complement("A") == "T"
    assert reverse_complement("G") == "C"
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACGT") == "ACGTT"
    assert reverse_complement(reverse_complement("AACGACGTTA")) == "AACGACGTTA"


# ---------------------------------------------------------------------------
# Algorithm: gapless / substitutions
# ---------------------------------------------------------------------------


def test_exact_forward_match():
    assert align_read("ACGTACGT", "TACG", 0) == Alignment(
        start=1, end=5, strand="-", edits=0, cigar="4M"
    )


def test_forward_wins_over_reverse_on_palindromic_read():
    # Read equals its own reverse complement; both strands match equally,
    # so the forward strand wins the tie.
    assert align_read("ACGTACGT", "ACGT", 0) == Alignment(
        start=0, end=4, strand="+", edits=0, cigar="4M"
    )


def test_substitution_edits():
    assert align_read("ACGTACGT", "TCGT", 1) == Alignment(
        start=0, end=4, strand="+", edits=1, cigar="4M"
    )
    # Two substitutions exceed max_edits=1.
    assert align_read("ACGT", "TGGT", 1) is None
    # ...but are accepted at max_edits=2.
    assert align_read("ACGT", "TGGT", 2) == Alignment(
        start=0, end=4, strand="+", edits=2, cigar="4M"
    )


def test_start_then_end_tiebreak():
    # "AA" matches at three gapless positions with equal edits; earliest
    # start wins.
    assert align_read("AAAAA", "AA", 0) == Alignment(
        start=0, end=2, strand="+", edits=0, cigar="2M"
    )


# ---------------------------------------------------------------------------
# Algorithm: insertions / deletions
# ---------------------------------------------------------------------------


def test_single_insertion():
    # read ACAGT vs reference ACGT: one inserted A at read offset 2.
    result = align_read("ACGT", "ACAGT", 1)
    assert result == Alignment(
        start=0, end=4, strand="+", edits=1, cigar="2M1I2M"
    )


def test_single_deletion():
    # ref AAGCACG, read AGACG: skipping reference base G at offset 2 of
    # interval [1,7) gives a perfect single-deletion alignment.
    result = align_read("AAGCACG", "AGACG", 1)
    assert result == Alignment(
        start=1, end=7, strand="+", edits=1, cigar="2M1D3M"
    )


def test_insertion_and_substitution():
    # ref CTTGTC, read CTAGGT: insert A at read offset 2 (one indel),
    # with the read's second G mismatching reference T -> 2 edits.
    result = align_read("CTTGTC", "CTAGGT", 2)
    assert result == Alignment(
        start=0, end=5, strand="+", edits=2, cigar="2M1I3M"
    )
    # Nothing qualifies with only one edit.
    assert align_read("CTTGTC", "CTAGGT", 1) is None


def test_edge_insertion_requires_nonempty_interval():
    # A 2-base read against a 1-base reference: the whole reference is
    # the only non-empty interval; a leading/trailing single-base
    # insertion covers it with one indel edit and zero substitutions.
    result = align_read("C", "CA", 1)
    assert result == Alignment(
        start=0, end=1, strand="+", edits=1, cigar="1M1I"
    )
    assert align_read("C", "CA", 0) is None
    # A 1-base read cannot create a non-empty interval via an insertion
    # (that would cover zero reference bases); only gapless is allowed.
    assert _best_on_strand("A", "T", "+", 1) == Alignment(
        start=0, end=1, strand="+", edits=1, cigar="1M"
    )
    # Every returned alignment covers a non-empty reference interval.
    assert align_read("A", "G", 1) == Alignment(
        start=0, end=1, strand="+", edits=1, cigar="1M"
    )
    # A 5-base read can never cover the single base with at most one
    # single-base insertion.
    assert align_read("C", "ACGTA", 5) is None


def test_cigar_tiebreak_keeps_first_indel_position():
    # ref AAAA, read AAAAA: an inserted A at any offset yields a perfect
    # alignment with identical (edits, indels, start, end); the smallest
    # cigar wins, i.e. the earliest insertion ("1I4M").
    result = align_read("AAAA", "AAAAA", 1)
    assert result == Alignment(
        start=0, end=4, strand="+", edits=1, cigar="1I4M"
    )


def test_deletion_on_minus_strand():
    # Mirror of the forward deletion case: read = rc("AGACG") aligns on
    # the minus strand over interval [1,7) with one deletion.
    read = reverse_complement("AGACG")
    result = align_read("AAGCACG", read, 1)
    assert result == Alignment(
        start=1, end=7, strand="-", edits=1, cigar="2M1D3M"
    )


# ---------------------------------------------------------------------------
# Algorithm: ordering rules
# ---------------------------------------------------------------------------


def test_fewer_edits_preferred_over_earlier_start():
    # Exact match at position 4 beats a one-substitution placement at 0.
    result = align_read("TAAAACGT", "ACGT", 1)
    assert result == Alignment(
        start=4, end=8, strand="+", edits=0, cigar="4M"
    )


def test_gapless_preferred_over_indel_at_equal_edits():
    # At position 0 a one-substitution gapless placement exists; an
    # equally good indel placement must lose on indel count.
    # ref ACGT, read ACGA -> gapless 1 substitution.
    result = align_read("TACGT", "ACGA", 2)
    assert result.cigar == "4M"
    assert (result.edits, result.start, result.end) == (1, 1, 5)


def test_cigar_lexicographic_tiebreak_between_insertions():
    # ref ACGT, read ACCGT: a perfect alignment exists via an inserted C
    # either after read offset 1 ("1M1I3M") or after offset 2
    # ("2M1I2M"); both have identical (edits, indels, start, end,
    # strand), so the lexicographically smaller cigar wins.
    result = align_read("ACGT", "ACCGT", 1)
    assert result == Alignment(
        start=0, end=4, strand="+", edits=1, cigar="1M1I3M"
    )


def test_forward_preferred_at_full_tie():
    # Reference with reverse-complement symmetry for the placement so the
    # only discriminator is the strand.
    # ref ACGT, read ACGT: forward exact at [0,4); reverse exact too.
    result = align_read("ACGT", "ACGT", 0)
    assert result.strand == "+"


def test_miss_returns_none():
    assert align_read("ACGT", "GGGG", 0) is None
    assert align_read("ACGT", "GGGG", 1) is None
    # read longer than reference by more than one inserted base: one
    # single-base insertion cannot close the gap.
    assert align_read("AC", "ACGGG", 5) is None


def test_align_reads_preserves_order():
    result = align_reads("ACGT", ["ACGT", "GGGG", "ACG"], 0)
    assert result[0] == Alignment(0, 4, "+", 0, "4M")
    assert result[1] is None
    assert result[2] == Alignment(0, 3, "+", 0, "3M")


# ---------------------------------------------------------------------------
# Exhaustive cross-check with a tiny DP oracle
# ---------------------------------------------------------------------------


def _oracle(reference: str, read: str, max_edits: int):
    """Enumerate placements directly: all intervals, strands, indel
    modes; return the best Alignment or None."""
    ref_len = len(reference)
    n = len(read)
    best = None

    def score_gapless(seq, start):
        if start + n > ref_len:
            return None
        edits = sum(seq[i] != reference[start + i] for i in range(n))
        return edits

    def score_insertion(seq, start, k):
        # insert seq[k], interval [start, start + n - 1); must be non-empty
        if n < 2 or start + n - 1 > ref_len:
            return None
        left = sum(
            seq[i] != reference[start + i] for i in range(k)
        )
        right = sum(
            seq[k + 1 + j] != reference[start + k + j]
            for j in range(n - k - 1)
        )
        edits = left + right + 1
        runs = [(k, "M"), (1, "I"), (n - k - 1, "M")]
        return edits, runs, start + n - 1

    def score_deletion(seq, start, k):
        # delete reference[start + k], k matched bases before it
        if start + n + 1 > ref_len or not (1 <= k <= n - 1):
            return None
        left = sum(seq[i] != reference[start + i] for i in range(k))
        right = sum(
            seq[k + j] != reference[start + k + 1 + j]
            for j in range(n - k)
        )
        edits = left + right + 1
        runs = [(k, "M"), (1, "D"), (n - k, "M")]
        return edits, runs, start + n + 1

    def cigar_of(runs):
        merged = []
        for count, op in runs:
            if count == 0:
                continue
            if merged and merged[-1][1] == op:
                merged[-1] = (merged[-1][0] + count, op)
            else:
                merged.append((count, op))
        return "".join(f"{c}{o}" for c, o in merged)

    for seq, strand in ((read, "+"), (reverse_complement(read), "-")):
        for start in range(ref_len):
            g = score_gapless(seq, start)
            if g is not None and g <= max_edits:
                a = Alignment(start, start + n, strand, g, f"{n}M")
                if best is None or _key(a) < _key(best):
                    best = a
            for k in range(0, n):
                ins = score_insertion(seq, start, k)
                if ins is not None and ins[0] <= max_edits:
                    edits, runs, end = ins
                    a = Alignment(start, end, strand, edits, cigar_of(runs))
                    if best is None or _key(a) < _key(best):
                        best = a
            for k in range(1, n):
                dele = score_deletion(seq, start, k)
                if dele is not None and dele[0] <= max_edits:
                    edits, runs, end = dele
                    a = Alignment(start, end, strand, edits, cigar_of(runs))
                    if best is None or _key(a) < _key(best):
                        best = a
    return best


def test_matches_exhaustive_oracle():
    rng = random.Random(2026)
    for _ in range(3000):
        ref_len = rng.randint(1, 10)
        n = rng.randint(1, 8)
        reference = "".join(rng.choice("ACGT") for _ in range(ref_len))
        read = "".join(rng.choice("ACGT") for _ in range(n))
        max_edits = rng.randint(0, 3)
        got = align_read(reference, read, max_edits)
        want = _oracle(reference, read, max_edits)
        assert got == want, (
            reference,
            read,
            max_edits,
            None if got is None else (got.start, got.end, got.strand, got.edits, got.cigar),
            None if want is None else (want.start, want.end, want.strand, want.edits, want.cigar),
        )


# ---------------------------------------------------------------------------
# API tests
# ---------------------------------------------------------------------------


def _payload(**overrides):
    payload = {
        "reference": "ACGTACGT",
        "reads": ["TACG", "AACGT", "ACGTAA", "CCCC"],
        "max_edits": 1,
    }
    payload.update(overrides)
    return payload


def test_api_ok():
    resp = client.post("/align", json=_payload())
    assert resp.status_code == 200
    assert resp.json() == {
        "alignments": [
            {"start": 1, "end": 5, "strand": "-", "edits": 0, "cigar": "4M"},
            {"start": 0, "end": 5, "strand": "-", "edits": 1, "cigar": "5M"},
            {"start": 0, "end": 6, "strand": "+", "edits": 1, "cigar": "6M"},
            None,
        ]
    }


def test_api_indel_shapes():
    resp = client.post(
        "/align",
        json={
            "reference": "AAGCACG",
            "reads": ["AACGC", "AGACG"],
            "max_edits": 1,
        },
    )
    assert resp.status_code == 200
    data = resp.json()["alignments"]
    assert data[0] == {
        "start": 0, "end": 4, "strand": "+", "edits": 1, "cigar": "2M1I2M"
    }
    assert data[1] == {
        "start": 1, "end": 7, "strand": "+", "edits": 1, "cigar": "2M1D3M"
    }


def test_api_zero_edits_requires_perfect_match():
    resp = client.post(
        "/align",
        json={"reference": "ACGT", "reads": ["ACGA", "ACGT"], "max_edits": 0},
    )
    assert resp.status_code == 200
    assert resp.json()["alignments"] == [
        None,
        {"start": 0, "end": 4, "strand": "+", "edits": 0, "cigar": "4M"},
    ]


def test_api_boundary_sizes_accepted():
    # 1000-base reference, 20 reads of length 80; all-A exact matches.
    reference = "A" * 1000
    reads = ["A" * 80 for _ in range(20)]
    resp = client.post(
        "/align",
        json={"reference": reference, "reads": reads, "max_edits": 5},
    )
    assert resp.status_code == 200
    alignments = resp.json()["alignments"]
    assert len(alignments) == 20
    assert all(
        a == {"start": 0, "end": 80, "strand": "+", "edits": 0, "cigar": "80M"}
        for a in alignments
    )


def test_api_single_read_minimums():
    resp = client.post(
        "/align",
        json={"reference": "A", "reads": ["T"], "max_edits": 0},
    )
    assert resp.status_code == 200
    # rc(T) == A matches the single-base reference on the minus strand.
    assert resp.json() == {
        "alignments": [
            {"start": 0, "end": 1, "strand": "-", "edits": 0, "cigar": "1M"}
        ]
    }


def test_api_max_edits_boundaries():
    body = {"reference": "ACGT", "reads": ["ACGT"]}
    assert client.post("/align", json={**body, "max_edits": 0}).status_code == 200
    assert client.post("/align", json={**body, "max_edits": 5}).status_code == 200


def _assert_4xx(payload):
    resp = client.post("/align", json=payload)
    assert 400 <= resp.status_code < 500, payload


def test_api_rejects_invalid():
    _assert_4xx({})
    _assert_4xx({"reads": ["A"], "max_edits": 0})
    _assert_4xx({"reference": "ACG", "max_edits": 0})
    _assert_4xx({"reference": "ACG", "reads": ["A"]})
    # reference boundaries
    _assert_4xx({"reference": "", "reads": ["A"], "max_edits": 0})
    _assert_4xx({"reference": "A" * 1001, "reads": ["A"], "max_edits": 0})
    # reads boundaries
    _assert_4xx({"reference": "ACG", "reads": [], "max_edits": 0})
    _assert_4xx(
        {"reference": "ACG", "reads": ["A"] * 21, "max_edits": 0}
    )
    _assert_4xx({"reference": "ACG", "reads": [""], "max_edits": 0})
    _assert_4xx(
        {"reference": "ACG", "reads": ["A" * 81], "max_edits": 0}
    )
    # alphabet
    _assert_4xx({"reference": "ACGN", "reads": ["A"], "max_edits": 0})
    _assert_4xx({"reference": "ACG", "reads": ["acg"], "max_edits": 0})
    _assert_4xx({"reference": "ACG", "reads": ["AC U"], "max_edits": 0})
    # max_edits boundaries / types
    _assert_4xx({"reference": "ACG", "reads": ["A"], "max_edits": -1})
    _assert_4xx({"reference": "ACG", "reads": ["A"], "max_edits": 6})
    _assert_4xx({"reference": "ACG", "reads": ["A"], "max_edits": "1"})
    _assert_4xx({"reference": "ACG", "reads": ["A"], "max_edits": 1.5})
    # wrong field types
    _assert_4xx({"reference": 123, "reads": ["A"], "max_edits": 0})
    _assert_4xx({"reference": "ACG", "reads": [123], "max_edits": 0})
    _assert_4xx({"reference": "ACG", "reads": "A", "max_edits": 0})
