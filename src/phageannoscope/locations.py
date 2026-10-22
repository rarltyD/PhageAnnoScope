"""Coordinate primitives. Segments always occur in translation order.

Coordinates are zero based, half open; each segment has ``start0`` and
``end0``. A minus-strand segment is read from end0-1 towards start0.
"""

from __future__ import annotations

from Bio.Seq import Seq


def _checked(segments, length=None):
    result = []
    for segment in segments:
        start, end = int(segment["start0"]), int(segment["end0"])
        if start < 0 or end <= start or (length is not None and end > length):
            raise ValueError(f"Invalid half-open segment [{start},{end})")
        result.append({"start0": start, "end0": end})
    return result


def gff_to_segments(start1, end1, strand):
    """Convert one GFF inclusive interval, retaining strand separately."""
    if strand not in (1, -1):
        raise ValueError("Strand must be +1 or -1")
    if int(start1) < 1 or int(end1) < int(start1):
        raise ValueError("GFF start/end must be ordered, positive inclusive coordinates")
    return [{"start0": int(start1) - 1, "end0": int(end1)}]


def segments_to_gff(segments):
    """Return inclusive ``(start1,end1)`` tuples in translation order."""
    return [(s["start0"] + 1, s["end0"]) for s in _checked(segments)]


def _union(segments):
    intervals = sorted((s["start0"], s["end0"]) for s in _checked(segments))
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def overlap_nt(a, b):
    """Count overlap of segment unions; gaps and repeated bases do not count."""
    left, right = _union(a), _union(b)
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def coalesce_segments(segments, strand):
    """Combine adjacent pieces, preserving the break at the circular origin."""
    result = []
    for segment in _checked(segments):
        if result and strand == 1 and result[-1]["end0"] == segment["start0"]:
            result[-1]["end0"] = segment["end0"]
        elif result and strand == -1 and result[-1]["start0"] == segment["end0"]:
            result[-1]["start0"] = segment["start0"]
        else:
            result.append(dict(segment))
    return result


def map_rotated_segments(segments, shift, length, strand):
    """Map rotated positions back using ``(new_pos + shift) % length``."""
    if strand not in (1, -1) or length <= 0:
        raise ValueError("A positive genome length and a known strand are required")
    shift %= length
    mapped = []
    for s in _checked(segments, length):
        size = s["end0"] - s["start0"]
        start = (s["start0"] + shift) % length
        if start + size <= length:
            pieces = [{"start0": start, "end0": start + size}]
        else:
            pieces = [
                {"start0": start, "end0": length},
                {"start0": 0, "end0": start + size - length},
            ]
            if strand == -1:
                pieces.reverse()
        mapped.extend(pieces)
    return coalesce_segments(mapped, strand)


def extract_segments(sequence, segments, strand):
    """Extract each piece in its supplied translation order."""
    if strand not in (1, -1):
        raise ValueError("Cannot translate an unknown strand")
    sequence = str(sequence).upper()
    pieces = [sequence[s["start0"] : s["end0"]] for s in _checked(segments, len(sequence))]
    if strand == -1:
        pieces = [str(Seq(piece).reverse_complement()) for piece in pieces]
    return "".join(pieces)


def stop_positions(segments, strand, codon_start=1):
    """Return the terminal codon's three genomic positions, or ``None``.

    This is a coordinate calculation: callers must independently establish
    that the terminal codon is a stop under the effective translation table.
    """
    if strand not in (1, -1) or codon_start not in (1, 2, 3):
        return None
    segments = _checked(segments)
    size = sum(s["end0"] - s["start0"] for s in segments) - (codon_start - 1)
    if size < 3 or size % 3:
        return None
    terminal = []
    for s in reversed(segments):
        if strand == 1:
            positions = range(s["end0"] - 1, s["start0"] - 1, -1)
        else:
            positions = range(s["start0"], s["end0"])
        for position in positions:
            terminal.append(position)
            if len(terminal) == 3:
                return list(reversed(terminal))
    return None
