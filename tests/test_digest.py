"""Tests for the restriction-style digest and the /digest API."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.digester import (
    Fragment,
    Template,
    cut_positions,
    digest,
    reverse_complement_pattern,
)
from app.main import app

client = TestClient(app)


# ---------------------------------------------------------------------------
# reverse complement of degenerate patterns
# ---------------------------------------------------------------------------


def test_reverse_complement_pattern():
    assert reverse_complement_pattern("ACGT") == "ACGT"
    assert reverse_complement_pattern("CG") == "CG"
    assert reverse_complement_pattern("AATT") == "AATT"
    # N complements to N, R <-> Y.
    assert reverse_complement_pattern("ANRY") == "RYNT"
    assert reverse_complement_pattern("NNRRYY") == "RRYYNN"


# ---------------------------------------------------------------------------
# cut positions
# ---------------------------------------------------------------------------


def test_cut_positions_forward_and_reverse_complement():
    # Forward "AC" at 0/4/8 cuts at +1; reverse complement "GT" at 2/6/10
    # cuts at start + 2 - 1.
    cuts = cut_positions(
        "ACGTACGTACGT", [Template("AC", 1)], circular=False
    )
    assert cuts == [1, 3, 5, 7, 9, 11]


def test_cut_positions_deduplicates_palindromic_pattern():
    # "CG" is its own reverse complement and offset 1 sits in the middle,
    # so both orientations produce the same cuts.
    cuts = cut_positions(
        "ACGTACGTACGT", [Template("CG", 1)], circular=False
    )
    assert cuts == [2, 6, 10]


def test_cut_positions_offset_at_pattern_end():
    # Forward cuts at start + 2, reverse-complement cuts at start + 0.
    cuts = cut_positions(
        "ACGTACGTACGT", [Template("CG", 2)], circular=False
    )
    assert cuts == [1, 3, 5, 7, 9, 11]


def test_cut_positions_degenerate_bases():
    # R = A/G, Y = C/T; "RY" is its own reverse complement. Offset 0 puts
    # the reverse-complement cut on the end of the linear sequence.
    cuts = cut_positions(
        "ACGTACGTACGT", [Template("RY", 0)], circular=False
    )
    assert cuts == [0, 2, 4, 6, 8, 10, 12]


def test_cut_positions_n_matches_any_base():
    cuts = cut_positions(
        "AAAAAAAAAAAA", [Template("AN", 0)], circular=False
    )
    assert cuts == list(range(11))


def test_cut_positions_multiple_templates():
    cuts = cut_positions(
        "ACGTACGTACGT",
        [Template("CG", 1), Template("AC", 1)],
        circular=False,
    )
    assert cuts == [1, 2, 3, 5, 6, 7, 9, 10, 11]


def test_cut_positions_circular_wraps_past_end():
    sequence = "GTTAACGTTAAC"
    template = Template("CGT", 0)
    # Linear: forward "CGT" at 5, reverse complement "ACG" at 4 (cut 7).
    assert cut_positions(sequence, [template], circular=False) == [5, 7]
    # Circular adds the wrapped forward match at 11 and the wrapped
    # reverse-complement match at 10 (cut 13 % 12 == 1).
    assert cut_positions(sequence, [template], circular=True) == [
        1,
        5,
        7,
        11,
    ]


# ---------------------------------------------------------------------------
# digest: linear
# ---------------------------------------------------------------------------


def test_digest_linear_fragments():
    cuts, fragments = digest(
        "ACGTACGTACGT", "linear", [Template("CG", 1)]
    )
    assert cuts == [2, 6, 10]
    assert fragments == [
        Fragment(0, 2, "AC"),
        Fragment(2, 6, "GTAC"),
        Fragment(6, 10, "GTAC"),
        Fragment(10, 12, "GT"),
    ]


def test_digest_linear_no_cuts():
    cuts, fragments = digest(
        "AAAAAAAAAAAA", "linear", [Template("CC", 1)]
    )
    assert cuts == []
    assert fragments == [Fragment(0, 12, "AAAAAAAAAAAA")]


def test_digest_linear_cuts_on_boundaries_make_no_empty_fragments():
    cuts, fragments = digest(
        "ACGTACGTACGT", "linear", [Template("RY", 0)]
    )
    assert cuts == [0, 2, 4, 6, 8, 10, 12]
    assert fragments == [
        Fragment(0, 2, "AC"),
        Fragment(2, 4, "GT"),
        Fragment(4, 6, "AC"),
        Fragment(6, 8, "GT"),
        Fragment(8, 10, "AC"),
        Fragment(10, 12, "GT"),
    ]


# ---------------------------------------------------------------------------
# digest: circular
# ---------------------------------------------------------------------------


def test_digest_circular_last_fragment_wraps_origin():
    cuts, fragments = digest(
        "TTAACGTTAACG", "circular", [Template("CGT", 1)]
    )
    assert cuts == [5, 11]
    assert fragments == [
        Fragment(5, 11, "GTTAAC"),
        Fragment(11, 5, "GTTAAC"),
    ]


def test_digest_circular_three_cuts():
    cuts, fragments = digest(
        "ACGTACGTACGT", "circular", [Template("CG", 1)]
    )
    assert cuts == [2, 6, 10]
    assert fragments == [
        Fragment(2, 6, "GTAC"),
        Fragment(6, 10, "GTAC"),
        Fragment(10, 2, "GTAC"),
    ]


def test_digest_circular_no_cuts_returns_full_circle():
    cuts, fragments = digest(
        "AAAAAAAAAAAA", "circular", [Template("CC", 1)]
    )
    assert cuts == []
    assert fragments == [Fragment(0, 0, "AAAAAAAAAAAA")]


def test_digest_circular_single_cut_returns_full_circle_from_cut():
    cuts, fragments = digest(
        "AACAAAAAAAAA", "circular", [Template("CA", 1)]
    )
    assert cuts == [3]
    assert fragments == [Fragment(3, 3, "AAAAAAAAAAAC")]


def test_digest_circular_single_cut_at_origin():
    cuts, fragments = digest(
        "CAAAAAAAAAAA", "circular", [Template("CA", 0)]
    )
    assert cuts == [0]
    assert fragments == [Fragment(0, 0, "CAAAAAAAAAAA")]


# ---------------------------------------------------------------------------
# /digest API
# ---------------------------------------------------------------------------


def test_api_linear_digest():
    response = client.post(
        "/digest",
        json={
            "sequence": "ACGTACGTACGT",
            "topology": "linear",
            "templates": [{"pattern": "CG", "offset": 1}],
        },
    )
    assert response.status_code == 200
    assert response.json() == {
        "cuts": [2, 6, 10],
        "fragments": [
            {"start": 0, "end": 2, "sequence": "AC"},
            {"start": 2, "end": 6, "sequence": "GTAC"},
            {"start": 6, "end": 10, "sequence": "GTAC"},
            {"start": 10, "end": 12, "sequence": "GT"},
        ],
    }


def test_api_circular_digest():
    response = client.post(
        "/digest",
        json={
            "sequence": "TTAACGTTAACG",
            "topology": "circular",
            "templates": [{"pattern": "CGT", "offset": 1}],
        },
    )
    assert response.status_code == 200
    assert response.json() == {
        "cuts": [5, 11],
        "fragments": [
            {"start": 5, "end": 11, "sequence": "GTTAAC"},
            {"start": 11, "end": 5, "sequence": "GTTAAC"},
        ],
    }


def test_api_offset_equal_to_pattern_length_is_allowed():
    response = client.post(
        "/digest",
        json={
            "sequence": "ACGTACGTACGT",
            "topology": "linear",
            "templates": [{"pattern": "CG", "offset": 2}],
        },
    )
    assert response.status_code == 200
    assert response.json()["cuts"] == [1, 3, 5, 7, 9, 11]


def test_api_five_templates_allowed():
    response = client.post(
        "/digest",
        json={
            "sequence": "ACGTACGTACGT",
            "topology": "linear",
            "templates": [{"pattern": "CG", "offset": 1}] * 5,
        },
    )
    assert response.status_code == 200
    assert response.json()["cuts"] == [2, 6, 10]


_VALID = {
    "sequence": "ACGTACGTACGT",
    "topology": "linear",
    "templates": [{"pattern": "CG", "offset": 1}],
}


@pytest.mark.parametrize(
    "payload",
    [
        # sequence: too short, too long, illegal letters, lowercase
        {**_VALID, "sequence": "ACGTACGTACG"},
        {**_VALID, "sequence": "A" * 1001},
        {**_VALID, "sequence": "ACGTACGTACGN"},
        {**_VALID, "sequence": "acgtacgtacgt"},
        # topology: not linear/circular
        {**_VALID, "topology": "circle"},
        # templates: too few, too many
        {**_VALID, "templates": []},
        {**_VALID, "templates": [{"pattern": "CG", "offset": 1}] * 6},
        # pattern: too short, too long, illegal letters, lowercase
        {**_VALID, "templates": [{"pattern": "C", "offset": 0}]},
        {**_VALID, "templates": [{"pattern": "ACGTACGTACGTA", "offset": 0}]},
        {**_VALID, "templates": [{"pattern": "CX", "offset": 0}]},
        {**_VALID, "templates": [{"pattern": "cn", "offset": 0}]},
        # offset: negative, beyond pattern length, not an integer
        {**_VALID, "templates": [{"pattern": "CG", "offset": -1}]},
        {**_VALID, "templates": [{"pattern": "CG", "offset": 3}]},
        {**_VALID, "templates": [{"pattern": "CG", "offset": "1"}]},
        # missing fields
        {"topology": "linear", "templates": [{"pattern": "CG", "offset": 1}]},
        {"sequence": "ACGTACGTACGT", "templates": [{"pattern": "CG", "offset": 1}]},
        {"sequence": "ACGTACGTACGT", "topology": "linear"},
    ],
)
def test_api_invalid_payloads_return_4xx(payload):
    response = client.post("/digest", json=payload)
    assert response.status_code == 422
