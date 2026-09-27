"""Shortest common superstring (DNA read assembly) with deterministic tie-breaking.

Given 2..12 reads over the alphabet A/C/G/T (used as-is, no reverse
complement), find the shortest sequence containing every read as a
substring. Among all shortest candidates the lexicographically smallest
one is returned, so the result does not depend on the input order.
"""

from __future__ import annotations


def _max_overlap(a: str, b: str) -> int:
    """Length of the longest suffix of ``a`` that is a prefix of ``b``."""
    # Straightforward scan; inputs are tiny so clarity beats cleverness.
    for k in range(min(len(a), len(b)), 0, -1):
        if a.endswith(b[:k]):
            return k
    return 0


def shortest_common_superstring(reads: list[str]) -> str:
    """Return the lexicographically smallest shortest common superstring."""
    # Duplicated reads do not change the superstring; drop exact copies.
    unique = list(set(reads))
    # A read contained in another read is implied by it; remove it.
    kept = [r for r in unique if not any(r != o and r in o for o in unique)]

    n = len(kept)
    if n == 1:
        return kept[0]

    overlap = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                overlap[i][j] = _max_overlap(kept[i], kept[j])

    # Held-Karp style DP over subsets. dp[mask][i] is the best
    # (shortest, then lexicographically smallest) superstring of the
    # reads in ``mask`` whose last read is ``kept[i]`` (ending flush
    # with the end of the string). Storing the actual strings makes the
    # lexicographic tie-break exact; n <= 12 keeps this cheap.
    size = 1 << n
    dp: list[list[str | None]] = [[None] * n for _ in range(size)]
    for i in range(n):
        dp[1 << i][i] = kept[i]

    for mask in range(size):
        row = dp[mask]
        for i in range(n):
            cur = row[i]
            if cur is None:
                continue
            for j in range(n):
                if mask >> j & 1:
                    continue
                cand = cur + kept[j][overlap[i][j]:]
                slot = dp[mask | (1 << j)][j]
                if slot is None or (len(cand), cand) < (len(slot), slot):
                    dp[mask | (1 << j)][j] = cand

    full = size - 1
    return min(
        (dp[full][i] for i in range(n) if dp[full][i] is not None),
        key=lambda s: (len(s), s),
    )


def assemble(reads: list[str]) -> tuple[str, list[int]]:
    """Assemble reads and locate each read (input order, 0-based)."""
    sequence = shortest_common_superstring(reads)
    positions = [sequence.index(r) for r in reads]
    return sequence, positions
