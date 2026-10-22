"""Coordinate expectations use artificial sequences exclusively from fixtures."""

import json
from pathlib import Path

import pytest

from phageannoscope.locations import (
    extract_segments,
    gff_to_segments,
    map_rotated_segments,
    overlap_nt,
    segments_to_gff,
    stop_positions,
)

FIXTURES = Path(__file__).parent / "fixtures" / "synthetic"
SEQUENCES = json.loads((FIXTURES / "sequences.json").read_text())


@pytest.mark.parametrize("strand,stop", [(1, [197, 198, 199]), (-1, [101, 100, 99])])
def test_gff_inclusive_and_stop_order(strand, stop):
    segments = gff_to_segments(100, 200, strand)
    assert segments == [{"start0": 99, "end0": 200}]
    assert segments_to_gff(segments) == [(100, 200)]
    # 101 bp requires codon_start=3 to leave 99 translated bases.
    assert stop_positions(segments, strand, 3) == stop
    assert stop_positions(segments, strand, 1) is None


@pytest.mark.parametrize(
    "strand,key,segments,stop",
    [
        (1, "origin_plus", [{"start0": 9, "end0": 12}, {"start0": 0, "end0": 6}], [3, 4, 5]),
        (-1, "origin_minus", [{"start0": 0, "end0": 3}, {"start0": 6, "end0": 12}], [8, 7, 6]),
    ],
)
def test_origin_translation_order(strand, key, segments, stop):
    assert extract_segments(SEQUENCES[key], segments, strand) == SEQUENCES["complete"]
    assert stop_positions(segments, strand) == stop
    # The untouched gap must not be counted through a bounding box.
    gap = [{"start0": 6, "end0": 9}] if strand == 1 else [{"start0": 3, "end0": 6}]
    assert overlap_nt(segments, gap) == 0


@pytest.mark.parametrize(
    "strand,rotated,expected",
    [
        (1, [{"start0": 3, "end0": 12}], [{"start0": 9, "end0": 12}, {"start0": 0, "end0": 6}]),
        (-1, [{"start0": 0, "end0": 9}], [{"start0": 0, "end0": 3}, {"start0": 6, "end0": 12}]),
    ],
)
def test_rotation_back_mapping(strand, rotated, expected):
    assert map_rotated_segments(rotated, 6, 12, strand) == expected
    assert map_rotated_segments(expected, -6, 12, strand) == rotated


def test_overlap_is_union_not_bounding_box_or_double_count():
    a = [{"start0": 0, "end0": 4}, {"start0": 2, "end0": 6}, {"start0": 9, "end0": 12}]
    b = [{"start0": 3, "end0": 10}]
    assert overlap_nt(a, b) == 4


@pytest.mark.parametrize("args", [(0, 9, 1), (9, 1, -1), (1, 9, 0)])
def test_bad_gff_rejected(args):
    with pytest.raises(ValueError):
        gff_to_segments(*args)
