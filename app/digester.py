"""Restriction-style in silico digestion of linear and circular DNA.

Each template carries an IUPAC-light ``pattern`` over A/C/G/T/N/R/Y
(``N`` matches any base, ``R`` matches A/G, ``Y`` matches C/T) and an
``offset`` from 0 to the pattern length. Every occurrence of the pattern
*and* of its reverse complement is a recognition site:

* a forward-strand site starting at ``origin`` cuts at ``origin + offset``;
* a reverse-complement site at ``origin`` cuts at
  ``origin + len(pattern) - offset``.

On circular sequences sites may straddle the end of the sequence: the
sequence is conceptually doubled for matching and every cut is taken
modulo the sequence length. Linear sequences never wrap, so cuts are
plain positions in ``0..len(sequence)`` (a cut at either end is kept but
produces no fragment).

Fragments are the spans between consecutive cut points, ordered by
ascending start. A circular sequence with no cuts yields one fragment
``(0, 0, sequence)``; with a single cut it yields the whole circle read
from that point, ``(cut, cut, sequence)``.
"""

from __future__ import annotations

from dataclasses import dataclass

LINEAR = "linear"
CIRCULAR = "circular"

_BASE_SETS: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "N": frozenset("ACGT"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
}

_COMPLEMENT = str.maketrans(
    {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "Y", "Y": "R"}
)


@dataclass(frozen=True)
class DigestTemplate:
    """A recognition pattern with its forward-strand cut offset."""

    pattern: str
    offset: int


@dataclass(frozen=True)
class Fragment:
    """A half-open fragment; ``end`` may precede ``start`` only by wrapping."""

    start: int
    end: int
    sequence: str


def reverse_complement(pattern: str) -> str:
    """Reverse complement of an A/C/G/T/N/R/Y pattern (R<->Y)."""

    return pattern.translate(_COMPLEMENT)[::-1]


def _match_at(haystack: str, origin: int, symbols: list[frozenset[str]]) -> bool:
    """Whether ``symbols`` match ``haystack`` starting at ``origin``."""

    for k, allowed in enumerate(symbols):
        if haystack[origin + k] not in allowed:
            return False
    return True


def digest(
    sequence: str, topology: str, templates: list[DigestTemplate]
) -> tuple[list[int], list[Fragment]]:
    """Digest ``sequence`` with ``templates``.

    ``sequence`` must use A/C/G/T and be at least as long as every
    pattern; ``topology`` is ``"linear"`` or ``"circular"``. Returns the
    de-duplicated, ascending cut positions and the fragments ordered by
    ascending start.
    """

    n = len(sequence)
    circular = topology == CIRCULAR
    # Doubling lets sites (and the final wrapping fragment) read past the
    # end of a circular sequence without per-index modular arithmetic.
    doubled = sequence + sequence
    cuts: set[int] = set()

    for template in templates:
        pattern = template.pattern
        length = len(pattern)
        symbols = [_BASE_SETS[symbol] for symbol in pattern]
        reverse_symbols = [_BASE_SETS[symbol] for symbol in reverse_complement(pattern)]
        # Every origin is possible on a circle; on a linear molecule a
        # site must fit entirely inside the sequence.
        origins = range(n) if circular else range(n - length + 1)
        for origin in origins:
            source = doubled if circular else sequence
            if _match_at(source, origin, symbols):
                cut = origin + template.offset
                cuts.add(cut % n if circular else cut)
            if _match_at(source, origin, reverse_symbols):
                cut = origin + length - template.offset
                cuts.add(cut % n if circular else cut)

    ordered = sorted(cuts)
    fragments: list[Fragment]
    if circular:
        if not ordered:
            # No recognition site: the intact circle, reported at origin 0.
            fragments = [Fragment(0, 0, sequence)]
        else:
            # One cut closes the loop onto itself (start == end); with
            # several cuts the last fragment runs across position 0.
            bounds = ordered + [ordered[0] + n]
            fragments = [
                Fragment(
                    bounds[i] % n,
                    bounds[i + 1] % n,
                    doubled[bounds[i] : bounds[i + 1]],
                )
                for i in range(len(ordered))
            ]
    else:
        # Cuts at 0/n coincide with the molecule ends and simply yield no
        # (empty) fragment; distinct interior cuts always give non-empty
        # spans.
        bounds = [0, *ordered, n]
        fragments = [
            Fragment(start, end, sequence[start:end])
            for start, end in zip(bounds, bounds[1:])
            if end > start
        ]
    return ordered, fragments
