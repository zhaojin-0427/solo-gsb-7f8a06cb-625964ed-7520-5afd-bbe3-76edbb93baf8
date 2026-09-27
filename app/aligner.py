"""Local semi-global alignment of short DNA reads against a reference.

For each read, and independently for its reverse complement, every
non-empty interval of the reference is considered and the read is aligned
end-to-end against it (free end gaps are *not* allowed beyond choosing
the interval itself).

Allowed alignments:

* substitutions only (gapless, ``M`` CIGAR);
* exactly one single-base insertion (an extra read base) or deletion (a
  skipped reference base), together with any number of substitutions.

The edit count is the number of substitutions plus the number of
insertions/deletions (so at most one indel event). Only alignments at or
below the caller-supplied ``max_edits`` are kept.

Ties resolve by the global key ``(edits, indels, start, end, strand_rank,
cigar)`` with ``+`` ranking before ``-`` and the CIGAR compared as plain
ASCII text.
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Alignment:
    """One full-length placement of a read on the reference."""

    start: int
    end: int
    strand: str
    edits: int
    cigar: str


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an A/C/G/T sequence."""

    return sequence.translate(_COMPLEMENT)[::-1]


def _encode_cigar(operations: list[tuple[int, str]]) -> str:
    """Runs of ``(count, op)``; zero counts are dropped, equal adjacent
    operations merged."""

    merged: list[list] = []
    for count, op in operations:
        if count == 0:
            continue
        if merged and merged[-1][1] == op:
            merged[-1][0] += count
        else:
            merged.append([count, op])
    return "".join(f"{count}{op}" for count, op in merged)


def _indel_count(cigar: str) -> int:
    """Number of I/D operations in a CIGAR (I/D runs are never merged)."""

    return cigar.count("I") + cigar.count("D")


def _best_on_strand(
    reference: str, read: str, strand: str, max_edits: int
) -> Alignment | None:
    """Best alignment of ``read`` against ``reference`` on one strand."""

    ref_len = len(reference)
    n = len(read)
    # Best key so far, split so candidates losing on the (edits, indels,
    # start, end) prefix never pay for CIGAR construction.
    best_prefix: tuple[int, int, int, int] | None = None
    best_cigar: str | None = None
    best: Alignment | None = None

    def consider(
        start: int,
        end: int,
        edits: int,
        indels: int,
        operations: list[tuple[int, str]],
    ) -> None:
        nonlocal best_prefix, best_cigar, best
        if edits > max_edits:
            return
        prefix = (edits, indels, start, end)
        if best_prefix is not None and prefix > best_prefix:
            return
        cigar = _encode_cigar(operations)
        if (
            best_prefix is None
            or prefix < best_prefix
            or cigar < best_cigar  # type: ignore[operator]
        ):
            best_prefix = prefix
            best_cigar = cigar
            best = Alignment(
                start=start, end=end, strand=strand, edits=edits, cigar=cigar
            )

    for start in range(ref_len):
        # --- Gapless: read[k] vs reference[start + k], span n bases. ---
        mismatches = 0
        for k in range(n):
            pos = start + k
            if pos >= ref_len:
                mismatches = max_edits + 1
                break
            if read[k] != reference[pos]:
                mismatches += 1
                if mismatches > max_edits:
                    break
        else:
            consider(start, start + n, mismatches, 0, [(n, "M")])
            if (
                best_prefix is not None
                and best_prefix[0] == 0
                and best_prefix[1] == 0
            ):
                # A perfect gapless placement: starts only increase from
                # here and indels cost at least one, so nothing later can
                # tie or beat it.
                break

        if max_edits == 0:
            # Indels cost at least one edit, so with max_edits == 0 only
            # the gapless placements are eligible.
            continue

        # Mismatch counts are only interesting while below max_edits:
        # edits = left[k] + right[k] + 1 must stay <= max_edits, so any
        # entry reaching max_edits makes its candidates unusable. The
        # arrays are pre-filled with that sentinel and the scans stop as
        # soon as the running count hits it.
        sentinel = max_edits

        # --- One single-base insertion at read offset k (0..n-1). ---
        # read[:k] matches reference[start:start+k]; read[k] is the
        # inserted base; read[k+1:] matches reference[start+k:start+n-1].
        # The interval is [start, start + n - 1): n-1 reference bases.
        # It must be non-empty, so n >= 2 (a 1-base read cannot cover any
        # reference base with an insertion). k == 0 and k == n-1 are
        # leading/trailing insertions; unlike edge deletions they are not
        # dominated (they need not fit an extra reference base), kept.
        if n >= 2 and start + n - 1 <= ref_len:
            left = [sentinel] * (n + 1)
            left[0] = 0  # nothing precedes the inserted base
            count = 0
            for k in range(1, n):
                if read[k - 1] != reference[start + k - 1]:
                    count += 1
                if count >= sentinel:
                    break
                left[k] = count
            right = [sentinel] * (n + 1)
            right[n - 1] = 0  # nothing follows the inserted base
            count = 0
            for k in range(n - 2, -1, -1):
                if read[k + 1] != reference[start + k]:
                    count += 1
                if count >= sentinel:
                    break
                right[k] = count
            for k in range(0, n):
                edits = left[k] + right[k] + 1
                if edits <= max_edits:
                    consider(
                        start,
                        start + n - 1,
                        edits,
                        1,
                        [(k, "M"), (1, "I"), (n - k - 1, "M")],
                    )

        # --- One single-base deletion at reference[start+k]. ---
        # read[:k] matches reference[start:start+k]; reference[start+k]
        # is skipped; read[k:] matches reference[start+k+1:start+n+1].
        # The interval is [start, start + n + 1): n+1 reference bases.
        # k == 0 and k == n (leading/trailing deletions) are dominated by
        # the gapless placement over the same read, so they are skipped.
        if start + n + 1 <= ref_len:
            left = [sentinel] * (n + 1)
            left[0] = 0
            count = 0
            for k in range(1, n):
                if read[k - 1] != reference[start + k - 1]:
                    count += 1
                if count >= sentinel:
                    break
                left[k] = count
            right = [sentinel] * (n + 1)
            count = 0
            for k in range(n - 1, 0, -1):
                if read[k] != reference[start + k + 1]:
                    count += 1
                if count >= sentinel:
                    break
                right[k] = count
            for k in range(1, n):
                edits = left[k] + right[k] + 1
                if edits <= max_edits:
                    consider(
                        start,
                        start + n + 1,
                        edits,
                        1,
                        [(k, "M"), (1, "D"), (n - k, "M")],
                    )

    return best


def align_read(
    reference: str, read: str, max_edits: int
) -> Alignment | None:
    """Best full-length alignment of ``read`` (or its reverse complement).

    Returns ``None`` when no placement has at most ``max_edits`` edits.
    """

    candidates = [
        alignment
        for alignment in (
            _best_on_strand(reference, read, "+", max_edits),
            _best_on_strand(
                reference, reverse_complement(read), "-", max_edits
            ),
        )
        if alignment is not None
    ]
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda a: (
            a.edits,
            _indel_count(a.cigar),
            a.start,
            a.end,
            0 if a.strand == "+" else 1,
            a.cigar,
        ),
    )


def align_reads(
    reference: str, reads: list[str], max_edits: int
) -> list[Alignment | None]:
    """Align each read in input order; misses are ``None``."""

    return [align_read(reference, read, max_edits) for read in reads]
