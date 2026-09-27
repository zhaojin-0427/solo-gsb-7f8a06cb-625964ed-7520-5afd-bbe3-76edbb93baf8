"""Restriction-style digest of a DNA sequence against degenerate templates.

Each template carries a degenerate ``pattern`` (A/C/G/T plus N = any
base, R = A/G, Y = C/T) and an integer ``offset``. The pattern and its
reverse complement are located on the sequence; a forward match starting
at ``s`` cuts at ``s + offset`` and a reverse-complement match cuts at
``s + len(pattern) - offset``. On circular molecules matches may wrap
past the end of the sequence and cut positions are taken modulo the
sequence length.

Cut positions are deduplicated and sorted. Linear molecules are split at
0, every cut and the end; circular molecules are split between adjacent
cuts with the last fragment wrapping across the origin. A circle with no
cuts yields a single ``(0, 0, sequence)`` fragment and a circle with a
single cut yields the full circle starting at that cut (``start == end``
in both cases). Fragments are never empty.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Bases matched by each pattern letter (N = any, R = A/G, Y = C/T).
_MATCHES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "N": frozenset("ACGT"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
}

#: Complement of each pattern letter; degenerate letters stay degenerate.
_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "R": "Y",
    "Y": "R",
}


@dataclass(frozen=True)
class Template:
    """A digest motif: a degenerate pattern plus a cut offset."""

    pattern: str
    offset: int


@dataclass(frozen=True)
class Fragment:
    """One digest product; ``end <= start`` wraps the origin (a full
    circle when equal)."""

    start: int
    end: int
    sequence: str


def reverse_complement_pattern(pattern: str) -> str:
    """Return the reverse complement of an A/C/G/T/N/R/Y pattern."""

    return "".join(_COMPLEMENT[base] for base in reversed(pattern))


def _match_starts(sequence: str, pattern: str, circular: bool) -> list[int]:
    """Start positions where ``pattern`` matches the sequence.

    Circular molecules allow matches wrapping past the end; every
    position is then a candidate start.
    """

    n = len(sequence)
    span = len(pattern)
    limit = n if circular else n - span + 1
    starts = []
    for start in range(limit):
        if all(
            sequence[(start + i) % n] in _MATCHES[pattern[i]]
            for i in range(span)
        ):
            starts.append(start)
    return starts


def cut_positions(
    sequence: str, templates: list[Template], circular: bool
) -> list[int]:
    """Sorted, deduplicated cut positions produced by every template.

    A forward match at ``s`` cuts at ``s + offset``; a reverse-complement
    match cuts at ``s + len(pattern) - offset``. Circular cuts are taken
    modulo the sequence length; linear cuts are kept as computed (they
    may coincide with 0 or the end of the sequence).
    """

    n = len(sequence)
    cuts: set[int] = set()
    for template in templates:
        span = len(template.pattern)
        offset = template.offset
        for start in _match_starts(sequence, template.pattern, circular):
            cut = start + offset
            cuts.add(cut % n if circular else cut)
        rc_pattern = reverse_complement_pattern(template.pattern)
        for start in _match_starts(sequence, rc_pattern, circular):
            cut = start + span - offset
            cuts.add(cut % n if circular else cut)
    return sorted(cuts)


def _linear_fragments(sequence: str, cuts: list[int]) -> list[Fragment]:
    """Split at 0, every cut and the end.

    Boundaries are deduplicated before slicing, so no empty fragments
    are produced even when a cut lands on 0 or the end.
    """

    boundaries = sorted({0, len(sequence), *cuts})
    return [
        Fragment(start=a, end=b, sequence=sequence[a:b])
        for a, b in zip(boundaries, boundaries[1:])
    ]


def _circular_fragments(sequence: str, cuts: list[int]) -> list[Fragment]:
    """Split between adjacent cuts; the last fragment wraps the origin."""

    if not cuts:
        return [Fragment(start=0, end=0, sequence=sequence)]
    if len(cuts) == 1:
        (cut,) = cuts
        return [
            Fragment(
                start=cut,
                end=cut,
                sequence=sequence[cut:] + sequence[:cut],
            )
        ]
    fragments = [
        Fragment(start=a, end=b, sequence=sequence[a:b])
        for a, b in zip(cuts, cuts[1:])
    ]
    last, first = cuts[-1], cuts[0]
    fragments.append(
        Fragment(
            start=last,
            end=first,
            sequence=sequence[last:] + sequence[:first],
        )
    )
    return fragments


def digest(
    sequence: str, topology: str, templates: list[Template]
) -> tuple[list[int], list[Fragment]]:
    """Digest ``sequence`` and return ``(cuts, fragments)``.

    ``topology`` is ``"linear"`` or ``"circular"``; fragments come out
    ordered by ascending ``start``.
    """

    circular = topology == "circular"
    cuts = cut_positions(sequence, templates, circular)
    if circular:
        fragments = _circular_fragments(sequence, cuts)
    else:
        fragments = _linear_fragments(sequence, cuts)
    return cuts, fragments
