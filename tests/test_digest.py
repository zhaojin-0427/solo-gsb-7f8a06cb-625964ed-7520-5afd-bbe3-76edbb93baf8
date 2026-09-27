"""Tests for in silico digestion and the /digest API."""

from __future__ import annotations

from fastapi.testclient import TestClient

from app.digester import (
    DigestTemplate,
    Fragment,
    digest,
    reverse_complement,
)
from app.main import app

client = TestClient(app)


def t(pattern: str, offset: int) -> DigestTemplate:
    return DigestTemplate(pattern=pattern, offset=offset)


# ---------------------------------------------------------------------------
# reverse complement (including degenerate symbols)
# ---------------------------------------------------------------------------


def test_reverse_complement():
    assert reverse_complement("GAATTC") == "GAATTC"
    assert reverse_complement("AATT") == "AATT"
    assert reverse_complement("AATC") == "GATT"
    assert reverse_complement("AR") == "YT"
    assert reverse_complement("NARY") == "RYTN"
    assert reverse_complement(reverse_complement("ACGNNRYY")) == "ACGNNRYY"


# ---------------------------------------------------------------------------
# Algorithm: linear topology
# ---------------------------------------------------------------------------


def test_linear_palindromic_site_two_staggered_cuts():
    # GAATTC (EcoRI-like palindrome) at origin 2, offset 1: cuts at 3
    # (forward) and 7 (reverse complement).
    cuts, fragments = digest("AAGAATTCTTAA", "linear", [t("GAATTC", 1)])
    assert cuts == [3, 7]
    assert fragments == [
        Fragment(0, 3, "AAG"),
        Fragment(3, 7, "AATT"),
        Fragment(7, 12, "CTTAA"),
    ]


def test_linear_single_interior_cut_site():
    # "AAAAATTAAAAA" (12): palindrome AATT at origin 3, offset 2 (its
    # midpoint) so both strands cut at 5.
    cuts, fragments = digest("AAAAATTAAAAA", "linear", [t("AATT", 2)])
    assert cuts == [5]
    assert fragments == [
        Fragment(0, 5, "AAAAA"),
        Fragment(5, 12, "TTAAAAA"),
    ]


def test_linear_only_reverse_complement_matches():
    # The forward pattern AATC is absent; its reverse complement GATT
    # occurs at origin 4. With offset 0 the cut is origin + 4 = 8.
    cuts, fragments = digest("GGGGGATTGGGG", "linear", [t("AATC", 0)])
    assert cuts == [8]
    assert fragments == [
        Fragment(0, 8, "GGGGGATT"),
        Fragment(8, 12, "GGGG"),
    ]


def test_linear_degenerate_pattern():
    # ARGY = [A][A/G][G][C/T] matches AAGT at origin 4; its reverse
    # complement RCYT matches nowhere. Offset 4 cuts right after the site.
    cuts, fragments = digest("CCCCAAGTCCCC", "linear", [t("ARGY", 4)])
    assert cuts == [8]
    assert fragments == [Fragment(0, 8, "CCCCAAGT"), Fragment(8, 12, "CCCC")]


def test_linear_n_matches_any_base():
    # AANA = [A][A][any][A] matches AAGA at origin 4; reverse complement
    # TTNT matches nowhere. Offset 4 cuts right after the site.
    cuts, fragments = digest("CCCCAAGACCCC", "linear", [t("AANA", 4)])
    assert cuts == [8]
    assert fragments == [Fragment(0, 8, "CCCCAAGA"), Fragment(8, 12, "CCCC")]


def test_linear_no_site_yields_whole_sequence():
    cuts, fragments = digest("CCCCCCCCCCCC", "linear", [t("AATT", 1)])
    assert cuts == []
    assert fragments == [Fragment(0, 12, "CCCCCCCCCCCC")]


def test_linear_sites_do_not_wrap_across_end():
    # On the circle AATT wraps across the end (origin 10); on the line it
    # does not match at all.
    cuts, fragments = digest("TTCCCCGGGGAA", "linear", [t("AATT", 1)])
    assert cuts == []
    assert fragments == [Fragment(0, 12, "TTCCCCGGGGAA")]


def test_linear_cut_at_start_produces_no_empty_fragment():
    # AATT at origin 0 with offset 0 gives cuts 0 and 4.
    cuts, fragments = digest("AATTAAAAAAAA", "linear", [t("AATT", 0)])
    assert cuts == [0, 4]
    assert fragments == [Fragment(0, 4, "AATT"), Fragment(4, 12, "AAAAAAAA")]
    assert all(fragment.end > fragment.start for fragment in fragments)


def test_linear_cut_at_end_produces_no_empty_fragment():
    # AATT at origin 8 with offset 4 gives cuts 8 and 12.
    cuts, fragments = digest("AAAAAAAAAATT", "linear", [t("AATT", 4)])
    assert cuts == [8, 12]
    assert fragments == [Fragment(0, 8, "AAAAAAAA"), Fragment(8, 12, "AATT")]
    assert all(fragment.end > fragment.start for fragment in fragments)


def test_linear_repeated_site():
    # "AAAATTAAAATT": palindrome AATT at origins 2 and 8; offset 2 is the
    # palindrome midpoint so both strands cut at origin + 2 -> {4, 10}.
    cuts, fragments = digest("AAAATTAAAATT", "linear", [t("AATT", 2)])
    assert cuts == [4, 10]
    assert fragments == [
        Fragment(0, 4, "AAAA"),
        Fragment(4, 10, "TTAAAA"),
        Fragment(10, 12, "TT"),
    ]


def test_multiple_templates_cuts_deduplicated():
    # "AAAAAATTCCCC": AATT (offset 1) at origin 4 cuts {5, 7}; AT
    # (self-reverse-complementary, offset 0) at origin 5 cuts {5, 7} too.
    cuts, fragments = digest(
        "AAAAAATTCCCC",
        "linear",
        [t("AATT", 1), t("AT", 0)],
    )
    assert cuts == [5, 7]
    assert fragments == [
        Fragment(0, 5, "AAAAA"),
        Fragment(5, 7, "AT"),
        Fragment(7, 12, "TCCCC"),
    ]


# ---------------------------------------------------------------------------
# Algorithm: circular topology
# ---------------------------------------------------------------------------


def test_circular_sites_wrap_across_end():
    # "TTCCCCGGGGAA": AATT reads across the boundary at origin 10;
    # offset 1 cuts at 11 (forward) and 1 (reverse complement).
    cuts, fragments = digest("TTCCCCGGGGAA", "circular", [t("AATT", 1)])
    assert cuts == [1, 11]
    assert fragments == [
        Fragment(1, 11, "TCCCCGGGGA"),
        Fragment(11, 1, "AT"),
    ]
    # The fragments in circular order reproduce the circle read from cut 1.
    assert fragments[0].sequence + fragments[1].sequence == (
        "TTCCCCGGGGAA"[1:] + "TTCCCCGGGGAA"[:1]
    )


def test_circular_two_cuts_fragments_sorted_by_start():
    # "AAAAGTACGTAC": palindrome GTAC at origins 4 and 8, offset 2 cuts
    # at 6 and 10; the second fragment wraps through position 0.
    cuts, fragments = digest("AAAAGTACGTAC", "circular", [t("GTAC", 2)])
    assert cuts == [6, 10]
    assert fragments == [
        Fragment(6, 10, "ACGT"),
        Fragment(10, 6, "ACAAAAGT"),
    ]
    assert [f.start for f in fragments] == sorted(f.start for f in fragments)


def test_circular_no_cut_returns_whole_circle_at_origin():
    sequence = "CCCCCCCCCCCC"
    cuts, fragments = digest(sequence, "circular", [t("AATT", 1)])
    assert cuts == []
    assert fragments == [Fragment(0, 0, sequence)]


def test_circular_single_cut_returns_one_full_turn():
    # "CCCCCAATTCCC": palindrome AATT at origin 5, offset 2 is its
    # midpoint so both strands cut once at 7; the fragment is the whole
    # circle read from that point.
    sequence = "CCCCCAATTCCC"
    cuts, fragments = digest(sequence, "circular", [t("AATT", 2)])
    assert cuts == [7]
    expected = sequence[7:] + sequence[:7]
    assert fragments == [Fragment(7, 7, expected)]
    assert fragments[0].sequence == "TTCCCCCCCCAA"
    assert len(fragments[0].sequence) == len(sequence)


def test_circular_cut_modulo_length():
    # A site ending past the boundary must wrap its cut modulo n: AATT
    # offset 4 (cut after the site) at wrapping origin 10 -> (10+4) % 12=2,
    # and the reverse cut (10+0) % 12 = 10.
    cuts, fragments = digest("TTCCCCGGGGAA", "circular", [t("AATT", 4)])
    assert cuts == [2, 10]
    assert fragments == [
        Fragment(2, 10, "CCCCGGGG"),
        Fragment(10, 2, "AATT"),
    ]


# ---------------------------------------------------------------------------
# API tests
# ---------------------------------------------------------------------------


def _payload(**overrides):
    payload = {
        "sequence": "AAAAAATTCCCC",
        "topology": "linear",
        "templates": [{"pattern": "AATT", "offset": 1}],
    }
    payload.update(overrides)
    return payload


def test_api_linear_ok():
    resp = client.post("/digest", json=_payload())
    assert resp.status_code == 200
    assert resp.json() == {
        "cuts": [5, 7],
        "fragments": [
            {"start": 0, "end": 5, "sequence": "AAAAA"},
            {"start": 5, "end": 7, "sequence": "AT"},
            {"start": 7, "end": 12, "sequence": "TCCCC"},
        ],
    }


def test_api_circular_ok():
    resp = client.post(
        "/digest",
        json={
            "sequence": "AAAAGTACGTAC",
            "topology": "circular",
            "templates": [{"pattern": "GTAC", "offset": 2}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "cuts": [6, 10],
        "fragments": [
            {"start": 6, "end": 10, "sequence": "ACGT"},
            {"start": 10, "end": 6, "sequence": "ACAAAAGT"},
        ],
    }


def test_api_circular_no_cut():
    resp = client.post(
        "/digest",
        json={
            "sequence": "CCCCCCCCCCCC",
            "topology": "circular",
            "templates": [{"pattern": "AATT", "offset": 0}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "cuts": [],
        "fragments": [
            {"start": 0, "end": 0, "sequence": "CCCCCCCCCCCC"}
        ],
    }


def test_api_boundary_sequence_sizes():
    assert (
        client.post(
            "/digest",
            json={
                "sequence": "A" * 12,
                "topology": "linear",
                "templates": [{"pattern": "AA", "offset": 0}],
            },
        ).status_code
        == 200
    )
    assert (
        client.post(
            "/digest",
            json={
                "sequence": "A" * 1000,
                "topology": "circular",
                "templates": [{"pattern": "A" * 12, "offset": 6}],
            },
        ).status_code
        == 200
    )


def test_api_pattern_and_template_limits():
    body = {"sequence": "A" * 12, "topology": "linear"}
    # 1 and 5 templates are fine; pattern lengths 2 and 12, offsets 0 and
    # the pattern length are fine.
    assert (
        client.post(
            "/digest",
            json={**body, "templates": [{"pattern": "AA", "offset": 0}]},
        ).status_code
        == 200
    )
    assert (
        client.post(
            "/digest",
            json={
                **body,
                "templates": [
                    {"pattern": "A" * 12, "offset": 12} for _ in range(5)
                ],
            },
        ).status_code
        == 200
    )


def _assert_4xx(payload):
    resp = client.post("/digest", json=payload)
    assert 400 <= resp.status_code < 500, payload


def test_api_rejects_invalid():
    _assert_4xx({})
    _assert_4xx({"sequence": "A" * 12})
    _assert_4xx({"topology": "linear", "templates": []})
    # sequence boundaries
    _assert_4xx(
        {
            "sequence": "A" * 11,
            "topology": "linear",
            "templates": [{"pattern": "AA", "offset": 0}],
        }
    )
    _assert_4xx(
        {
            "sequence": "A" * 1001,
            "topology": "linear",
            "templates": [{"pattern": "AA", "offset": 0}],
        }
    )
    # sequence alphabet (N/lowercase belong to patterns, not sequences)
    _assert_4xx(_payload(sequence="AAAAAAAAAAAN"))
    _assert_4xx(_payload(sequence="aaaaaaaaaaaa"))
    _assert_4xx(_payload(sequence="AAAAAAAAAA U"))
    # topology
    _assert_4xx(_payload(topology="plasmid"))
    _assert_4xx(_payload(topology=""))
    # template counts
    _assert_4xx(_payload(templates=[]))
    _assert_4xx(
        _payload(templates=[{"pattern": "AA", "offset": 0}] * 6)
    )
    # pattern lengths
    _assert_4xx(_payload(templates=[{"pattern": "A", "offset": 0}]))
    _assert_4xx(_payload(templates=[{"pattern": "A" * 13, "offset": 0}]))
    _assert_4xx(_payload(templates=[{"pattern": "", "offset": 0}]))
    # pattern alphabet
    _assert_4xx(_payload(templates=[{"pattern": "AAX", "offset": 0}]))
    _assert_4xx(_payload(templates=[{"pattern": "aa", "offset": 0}]))
    _assert_4xx(_payload(templates=[{"pattern": "A B", "offset": 0}]))
    # offset boundaries / types
    _assert_4xx(_payload(templates=[{"pattern": "AATT", "offset": -1}]))
    _assert_4xx(_payload(templates=[{"pattern": "AATT", "offset": 5}]))
    _assert_4xx(_payload(templates=[{"pattern": "AATT", "offset": "1"}]))
    _assert_4xx(_payload(templates=[{"pattern": "AATT", "offset": 1.5}]))
    _assert_4xx(_payload(templates=[{"pattern": "AATT"}]))
    # wrong field types
    _assert_4xx(_payload(sequence=123))
    _assert_4xx(_payload(templates={"pattern": "AA", "offset": 0}))
    _assert_4xx(_payload(templates=["AATT"]))
