"""Tests for single-base substitution calling and the /variants API."""

from __future__ import annotations

import math

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.variants import AlignedRead, call_variants

client = TestClient(app)


def read(start: int, sequence: str, quality: int | list[int]) -> AlignedRead:
    if isinstance(quality, int):
        quality = [quality] * len(sequence)
    return AlignedRead(start=start, sequence=sequence, qualities=quality)


# ---------------------------------------------------------------------------
# Algorithm tests
# ---------------------------------------------------------------------------


def test_no_reads_mismatch_no_variant():
    # A low-quality read contradicting the reference is not enough.
    variants = call_variants("ACGT", [read(0, "ATGT", 2)])
    assert variants == []


def test_single_confident_substitution():
    # One high-quality read contradicting the reference at one site.
    variants = call_variants("ACGT", [read(0, "ATGT", 40)])
    assert len(variants) == 1
    v = variants[0]
    assert (v.position, v.ref, v.alt, v.depth) == (1, "C", "T", 1)


def test_majority_wins_and_depth_counts_coverers():
    reference = "AAAA"
    variants = call_variants(
        reference,
        [
            read(0, "GAAA", 30),
            read(0, "GAAA", 30),
            read(0, "AAAA", 30),
        ],
    )
    assert len(variants) == 1
    assert (variants[0].position, variants[0].ref, variants[0].alt) == (
        0,
        "A",
        "G",
    )
    assert variants[0].depth == 3


def test_uncovered_positions_are_ignored():
    # Reference positions 3..9 are never covered.
    variants = call_variants("AAAAAAAAAA", [read(0, "GGG", [40, 40, 40])])
    assert [v.position for v in variants] == [0, 1, 2]
    assert all(v.depth == 1 for v in variants)


def test_partial_overlap_depth():
    reference = "AAAAA"
    variants = call_variants(
        reference,
        [read(0, "GGG", 40), read(2, "GGG", 40)],
    )
    by_pos = {v.position: v for v in variants}
    assert set(by_pos) == {0, 1, 2, 3, 4}
    assert by_pos[0].depth == 1
    assert by_pos[2].depth == 2  # covered by both reads
    assert by_pos[4].depth == 1


def test_results_sorted_by_position():
    # Positions 0..4 see only T -> confidently T. Positions 5..9 see T
    # and G at equal quality -> exact tie, suppressed by the cutoff.
    variants = call_variants(
        "AAAAAAAAAA",
        [read(0, "TTTTTTTTTT", 40), read(5, "GGGGG", 40)],
    )
    assert [(v.position, v.alt) for v in variants] == [
        (0, "T"),
        (1, "T"),
        (2, "T"),
        (3, "T"),
        (4, "T"),
    ]


def test_exact_tie_picks_smallest_base_then_unconfident():
    # With two equal-quality reads observing different bases, the two
    # candidates tie exactly. The fixed order A<C<G<T picks the smallest
    # one; either way the best/runner-up gap is zero, so the ln(10)
    # confidence rule suppresses the call.
    assert call_variants(
        "A",
        [read(0, "A", 40), read(0, "C", 40)],
    ) == []
    assert call_variants(
        "A",
        [read(0, "G", 40), read(0, "T", 40)],
    ) == []


def test_low_quality_discordant_reads_outvoted():
    # Two Q=2 reads observing T vs one Q=40 read observing C at an A
    # site: the strong C evidence wins confidently.
    variants = call_variants(
        "A",
        [read(0, "T", 2), read(0, "T", 2), read(0, "C", 40)],
    )
    assert len(variants) == 1
    assert variants[0].alt == "C"
    assert variants[0].depth == 3


def test_marginal_gap_below_cutoff_not_reported():
    # C@10 vs T@19 at an A site: T wins, but by less than ln(10).
    variants = call_variants(
        "A",
        [read(0, "C", 10), read(0, "T", 19)],
    )
    assert variants == []


def test_marginal_gap_above_cutoff_reported():
    # The slightly stronger T@20 pushes the same contest past ln(10).
    variants = call_variants(
        "A",
        [read(0, "C", 10), read(0, "T", 20)],
    )
    assert len(variants) == 1
    assert variants[0].alt == "T"


def test_lowest_quality_single_mismatch_is_unconfident():
    # At the minimum quality Q=2, p = 10**(-0.2) ~ 0.63: a lone
    # observation carries less than a 10-fold likelihood ratio.
    variants = call_variants("A", [read(0, "C", 2)])
    assert variants == []


def test_reference_site_reported():
    variants = call_variants("CGTA", [read(0, "CGTA", 40)])
    assert variants == []


def test_manual_likelihood_accumulation():
    # Recompute scores independently and compare with the reported call.
    reference = "A"
    reads = [read(0, "C", 10), read(0, "T", 20)]
    variants = call_variants(reference, reads)

    scores = {base: 0.0 for base in "ACGT"}
    for obs, q in [("C", 10), ("T", 20)]:
        p = 10 ** (-q / 10)
        for base in "ACGT":
            scores[base] += (
                math.log(1 - p) if base == obs else math.log(p / 3)
            )
    ranked = sorted(scores, key=lambda b: scores[b], reverse=True)
    assert ranked[0] == "T"
    gap = scores["T"] - scores["C"]
    assert gap >= math.log(10)
    assert len(variants) == 1
    assert (variants[0].ref, variants[0].alt, variants[0].depth) == (
        "A",
        "T",
        2,
    )


# ---------------------------------------------------------------------------
# API tests
# ---------------------------------------------------------------------------


def _payload(**overrides):
    payload = {
        "reference": "ACGTACGT",
        "reads": [
            {"start": 1, "sequence": "TG", "qualities": [40, 40]},
        ],
    }
    payload.update(overrides)
    return payload


def test_api_ok():
    resp = client.post("/variants", json=_payload())
    assert resp.status_code == 200
    assert resp.json() == {
        "variants": [
            {"position": 1, "ref": "C", "alt": "T", "depth": 1},
        ]
    }


def test_api_empty_variants():
    resp = client.post(
        "/variants",
        json={
            "reference": "ACGT",
            "reads": [
                {"start": 0, "sequence": "ACGT", "qualities": [40] * 4}
            ],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {"variants": []}


def test_api_depth_over_multiple_reads():
    resp = client.post(
        "/variants",
        json={
            "reference": "AAA",
            "reads": [
                {"start": 0, "sequence": "GAA", "qualities": [30] * 3},
                {"start": 0, "sequence": "GAA", "qualities": [30] * 3},
                {"start": 0, "sequence": "GAA", "qualities": [2] * 3},
            ],
        },
    )
    assert resp.status_code == 200
    data = resp.json()["variants"]
    assert len(data) == 1
    assert data[0]["position"] == 0
    assert data[0]["alt"] == "G"
    assert data[0]["depth"] == 3


@pytest.mark.parametrize(
    "payload",
    [
        # missing fields
        {},
        {"reference": "ACGT"},
        {"reads": []},
        # reference too short / too long
        {"reference": "", "reads": []},
        {
            "reference": "A" * 1001,
            "reads": [
                {"start": 0, "sequence": "A", "qualities": [40]}
            ],
        },
        # invalid reference character
        {
            "reference": "ACGN",
            "reads": [
                {"start": 0, "sequence": "A", "qualities": [40]}
            ],
        },
        # no reads / too many reads
        {"reference": "ACGT", "reads": []},
        {
            "reference": "A",
            "reads": [
                {"start": 0, "sequence": "A", "qualities": [40]}
                for _ in range(101)
            ],
        },
        # read out of range
        {
            "reference": "ACGT",
            "reads": [
                {"start": 3, "sequence": "TT", "qualities": [40, 40]}
            ],
        },
        # negative start
        {
            "reference": "ACGT",
            "reads": [
                {"start": -1, "sequence": "A", "qualities": [40]}
            ],
        },
        # sequence / qualities length mismatch
        {
            "reference": "ACGT",
            "reads": [
                {"start": 0, "sequence": "AC", "qualities": [40]}
            ],
        },
        # quality out of range
        {
            "reference": "ACGT",
            "reads": [
                {"start": 0, "sequence": "A", "qualities": [1]}
            ],
        },
        {
            "reference": "ACGT",
            "reads": [
                {"start": 0, "sequence": "A", "qualities": [41]}
            ],
        },
        # read too long / short
        {
            "reference": "A" * 81,
            "reads": [
                {
                    "start": 0,
                    "sequence": "A" * 81,
                    "qualities": [40] * 81,
                }
            ],
        },
        {
            "reference": "ACGT",
            "reads": [{"start": 0, "sequence": "", "qualities": []}],
        },
        # invalid base in read
        {
            "reference": "ACGT",
            "reads": [
                {"start": 0, "sequence": "NC", "qualities": [40, 40]}
            ],
        },
    ],
)
def test_api_rejects_invalid(payload):
    resp = client.post("/variants", json=payload)
    assert 400 <= resp.status_code < 500, payload


def test_api_rejects_wrong_types():
    resp = client.post(
        "/variants",
        json={
            "reference": 1234,
            "reads": [
                {"start": "0", "sequence": "A", "qualities": [40]}
            ],
        },
    )
    assert 400 <= resp.status_code < 500


def test_api_boundary_sizes_accepted():
    # 1000-base reference with 100 reads of length 80, all inside bounds.
    reference = "A" * 1000
    reads = [
        {
            "start": (i * 9) % (1000 - 80),
            "sequence": "A" * 80,
            "qualities": [40] * 80,
        }
        for i in range(100)
    ]
    resp = client.post(
        "/variants", json={"reference": reference, "reads": reads}
    )
    assert resp.status_code == 200
    assert resp.json() == {"variants": []}
