"""Single-base substitution calling from aligned, gapless forward reads.

Each read sits at a known 0-based ``start`` on the reference and contains
no insertions or deletions. For every covered reference position we sum,
independently over the four candidate bases A/C/G/T, the log-likelihood
of the observations:

* a quality value Q gives a misread probability ``p = 10**(-Q/10)``;
* observing a base equal to the candidate contributes ``ln(1 - p)``;
* observing any of the other three bases contributes ``ln(p / 3)``.

The highest-scoring candidate wins; exact ties are resolved by the fixed
order ``A < C < G < T``. A site is reported as a variant only when it is
covered, the winning candidate differs from the reference, and its score
exceeds the runner-up by at least ``ln(10)`` (a 10-fold likelihood ratio).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

BASES = "ACGT"
_BASE_INDEX = {base: i for i, base in enumerate(BASES)}

#: Minimum required log-likelihood gap between best and runner-up.
LIKELIHOOD_RATIO_CUTOFF = math.log(10.0)


@dataclass(frozen=True)
class AlignedRead:
    """A forward, gapless read placed at ``start`` on the reference."""

    start: int
    sequence: str
    qualities: list[int]


@dataclass(frozen=True)
class Variant:
    position: int
    ref: str
    alt: str
    depth: int


def _log_likelihood_table() -> tuple[list[float], list[float]]:
    """Precompute ln(1-p) and ln(p/3) for every allowed quality Q (2..40)."""

    match_ll: list[float] = [0.0] * 41
    mismatch_ll: list[float] = [0.0] * 41
    for quality in range(2, 41):
        p = 10.0 ** (-quality / 10.0)
        match_ll[quality] = math.log(1.0 - p)
        mismatch_ll[quality] = math.log(p / 3.0)
    return match_ll, mismatch_ll


def call_variants(
    reference: str, reads: list[AlignedRead]
) -> list[Variant]:
    """Call confident single-base substitutions against ``reference``.

    ``reads`` must already be validated: every read falls entirely inside
    the reference, sequences use A/C/G/T and qualities are 2..40 with one
    value per base. Returns variants ordered by ascending position.
    """

    ref_len = len(reference)
    # scores[pos][c] accumulates the log-likelihood of candidate BASES[c].
    scores = [[0.0] * len(BASES) for _ in range(ref_len)]
    depth = [0] * ref_len

    match_ll, mismatch_ll = _log_likelihood_table()

    for read in reads:
        start = read.start
        for offset, (base, quality) in enumerate(
            zip(read.sequence, read.qualities)
        ):
            pos = start + offset
            depth[pos] += 1
            row = scores[pos]
            obs = _BASE_INDEX[base]
            match_score = match_ll[quality]
            mismatch_score = mismatch_ll[quality]
            for candidate in range(len(BASES)):
                row[candidate] += (
                    match_score if candidate == obs else mismatch_score
                )

    variants: list[Variant] = []
    for pos in range(ref_len):
        if depth[pos] == 0:
            continue
        row = scores[pos]
        # Strict "greater than" keeps the first (A < C < G < T) index on
        # exact ties, which is the required tie-breaking order.
        best = 0
        for candidate in range(1, len(BASES)):
            if row[candidate] > row[best]:
                best = candidate
        runner_up = max(
            row[candidate]
            for candidate in range(len(BASES))
            if candidate != best
        )
        alt = BASES[best]
        if (
            alt != reference[pos]
            and row[best] - runner_up >= LIKELIHOOD_RATIO_CUTOFF
        ):
            variants.append(
                Variant(
                    position=pos,
                    ref=reference[pos],
                    alt=alt,
                    depth=depth[pos],
                )
            )

    return variants
