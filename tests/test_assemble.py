"""Tests for the assembly algorithm and the HTTP API."""

from __future__ import annotations

import itertools
import random
from collections import deque

import pytest
from fastapi.testclient import TestClient

from app.assembler import assemble, shortest_common_superstring
from app.main import app

client = TestClient(app)


# ---------------------------------------------------------------------------
# Reference implementations (only used on tiny inputs)
# ---------------------------------------------------------------------------

def _overlap(a: str, b: str) -> int:
    for k in range(min(len(a), len(b)), 0, -1):
        if a.endswith(b[:k]):
            return k
    return 0


def brute_force_scs(reads: list[str]) -> str:
    """Min over all permutations merged with maximal pairwise overlaps."""
    unique = list(set(reads))
    kept = [r for r in unique if not any(r != o and r in o for o in unique)]
    best: str | None = None
    for perm in itertools.permutations(kept):
        s = perm[0]
        for i in range(1, len(perm)):
            s += perm[i][_overlap(perm[i - 1], perm[i]):]
        if best is None or (len(s), s) < (len(best), best):
            best = s
    return best or ""


def bfs_oracle(reads: list[str]) -> str:
    """Fully independent oracle: enumerate all DNA strings in (length, lex)
    order and return the first one containing every read."""
    queue: deque[str] = deque([""])
    while queue:
        s = queue.popleft()
        if all(r in s for r in reads):
            return s
        for c in "ACGT":
            queue.append(s + c)
    raise AssertionError("unreachable")


def check_solution(reads: list[str], seq: str) -> None:
    for r in reads:
        assert r in seq, f"{r} not found in {seq}"


# ---------------------------------------------------------------------------
# Algorithm tests
# ---------------------------------------------------------------------------

def test_basic_example():
    seq, pos = assemble(["ACGT", "GTAC", "TACG"])
    assert seq == "GTACGT"  # ACGT@2, GTAC@0, TACG@1
    assert pos == [2, 0, 1]


def test_duplicates_and_containment():
    seq, pos = assemble(["ACGT", "ACGT", "CG"])
    assert seq == "ACGT"
    assert pos == [0, 0, 1]


def test_single_base_reads():
    seq, pos = assemble(["A", "C", "G", "T"])
    assert seq == "ACGT"  # lexicographically smallest among shortest
    assert pos == [0, 1, 2, 3]


def test_order_invariance():
    reads = ["ACGT", "GTAC", "TACG", "CGTA", "TAC"]
    results = {assemble(list(perm))[0] for perm in itertools.permutations(reads)}
    assert results == {"ACGTACG"}


def test_positions_follow_input_order():
    reads = ["TACG", "ACGT", "GTAC"]
    seq, pos = assemble(reads)
    assert seq == "GTACGT"
    assert pos == [1, 2, 0]
    for r, p in zip(reads, pos):
        assert seq[p : p + len(r)] == r


def test_earliest_position_reported():
    seq, pos = assemble(["AA", "AA", "AA"])
    assert seq == "AA"
    assert pos == [0, 0, 0]


def test_lexicographic_tiebreak():
    # ACGTA / CGTAC / TACGT are all shortest (5); smallest wins.
    seq, _ = assemble(["AC", "CG", "GT", "TA"])
    assert seq == "ACGTA"


@pytest.mark.parametrize("seed", range(40))
def test_against_brute_force(seed: int):
    rng = random.Random(seed)
    reads = [
        "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
        for _ in range(rng.randint(2, 6))
    ]
    seq = shortest_common_superstring(reads)
    check_solution(reads, seq)
    assert seq == brute_force_scs(reads)


@pytest.mark.parametrize("seed", range(8))
def test_against_bfs_oracle(seed: int):
    rng = random.Random(500 + seed)
    reads = [
        "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3)))
        for _ in range(rng.randint(2, 4))
    ]
    assert shortest_common_superstring(reads) == bfs_oracle(reads)


@pytest.mark.parametrize("seed", range(10))
def test_order_invariance_random(seed: int):
    rng = random.Random(1000 + seed)
    reads = [
        "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 10)))
        for _ in range(6)
    ]
    expected = shortest_common_superstring(reads)
    for _ in range(5):
        shuffled = reads[:]
        rng.shuffle(shuffled)
        assert shortest_common_superstring(shuffled) == expected


def test_max_size_input():
    rng = random.Random(7)
    reads = ["".join(rng.choice("ACGT") for _ in range(80)) for _ in range(12)]
    seq, pos = assemble(reads)
    check_solution(reads, seq)
    assert len(pos) == 12
    for r, p in zip(reads, pos):
        assert seq[p : p + len(r)] == r


# ---------------------------------------------------------------------------
# API tests
# ---------------------------------------------------------------------------

def test_api_ok():
    resp = client.post("/assemble", json={"reads": ["ACGT", "GTAC", "TACG"]})
    assert resp.status_code == 200
    assert resp.json() == {
        "sequence": "GTACGT",
        "length": 6,
        "positions": [2, 0, 1],
    }


@pytest.mark.parametrize(
    "reads",
    [
        ["ACGT"],               # too few
        ["ACGT"] * 13,          # too many
        [],                     # empty
        ["", "AC"],             # too short
        ["A" * 81, "ACGT"],     # too long
        ["ACGX", "ACGT"],       # invalid character
        ["acgt", "ACGT"],       # lowercase not allowed
        ["ACGT", 123],          # not a string
    ],
)
def test_api_rejects_invalid(reads):
    resp = client.post("/assemble", json={"reads": reads})
    assert 400 <= resp.status_code < 500


def test_api_rejects_missing_field():
    resp = client.post("/assemble", json={})
    assert 400 <= resp.status_code < 500


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}
