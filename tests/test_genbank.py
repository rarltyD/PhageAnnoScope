import json
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from phageannoscope.genbank_io import genome_from_record, normalize_cds, read_genbank

FIXTURES = Path(__file__).parent / "fixtures" / "synthetic"
SEQUENCES = json.loads((FIXTURES / "sequences.json").read_text())


def genome(key, topology="linear"):
    record = SeqRecord(Seq(SEQUENCES[key]), id="SYNTHETIC_LOCATION.2")
    record.annotations["topology"] = topology
    return genome_from_record(record)


def test_features_protein_id_missing_duplicate_and_default_code():
    data = read_genbank(FIXTURES / "synthetic.gb")
    assert len(data["features"]) == 7
    assert len(data["cds_calls"]) == 3
    calls = data["cds_calls"]
    assert len({c["feature_id"] for c in calls}) == 3
    assert calls[0]["protein_id"] == calls[1]["protein_id"]
    assert calls[2]["protein_id"] is None
    assert calls[0]["transl_table_effective"] == 1
    assert calls[0]["transl_table_qualifier_missing"] is True
    assert calls[1]["transl_table_effective"] == 11
    assert calls[1]["stop_triplet"] == [14, 13, 12]
    assert [c["recomputed_translation"] for c in calls] == ["MK", "MK", "MP"]
    assert all(c["eligible_primary"] for c in calls)
    again = read_genbank(FIXTURES / "synthetic.gb")
    assert again == data


@pytest.mark.parametrize(
    "strand,key,parts",
    [
        (1, "origin_plus", [(9, 12), (0, 6)]),
        (-1, "origin_minus", [(0, 3), (6, 12)]),
    ],
)
def test_origin_joins_both_strands(strand, key, parts):
    feature = SeqFeature(
        CompoundLocation([SimpleLocation(a, b, strand=strand) for a, b in parts], "join"),
        type="CDS",
    )
    result = normalize_cds(genome(key, "circular"), feature)
    assert result["segments"] == [{"start0": a, "end0": b} for a, b in parts]
    assert result["feature_nt_length"] == 9
    assert result["recomputed_translation"] == "MK"
    assert result["wraps_origin"] is True
    assert result["segments_contiguous"] is True
    assert result["eligible_primary"] is False
    assert "wraps_origin" in result["primary_exclusion_reasons"]


def test_actual_genbank_complement_join_parser_preserves_order():
    calls = read_genbank(FIXTURES / "circular.gb")["cds_calls"]
    assert [c["strand"] for c in calls] == [1, -1]
    assert [c["recomputed_translation"] for c in calls] == ["MK", "MK"]
    assert calls[0]["segments"] == [{"start0": 9, "end0": 12}, {"start0": 0, "end0": 6}]
    assert calls[1]["segments"] == [{"start0": 0, "end0": 3}, {"start0": 6, "end0": 12}]
    assert all(c["wraps_origin"] and not c["eligible_primary"] for c in calls)


@pytest.mark.parametrize("strand,expected5,expected3", [(1, True, False), (-1, False, True)])
def test_fuzzy_boundaries_follow_transcript_strand(strand, expected5, expected3):
    feature = SeqFeature(SimpleLocation(BeforePosition(0), 9, strand=strand), type="CDS")
    result = normalize_cds(genome("complete"), feature)
    assert result["partial_5p"] is expected5
    assert result["partial_3p"] is expected3
    assert result["eligible_primary"] is False


@pytest.mark.parametrize("key,offset", [("codon_start_2", 2), ("codon_start_3", 3)])
def test_codon_start_preserved(key, offset):
    g = genome(key)
    feature = SeqFeature(
        SimpleLocation(0, g["length_bp"], strand=1),
        type="CDS",
        qualifiers={"codon_start": [str(offset)]},
    )
    result = normalize_cds(g, feature)
    assert result["codon_start"] == offset
    assert result["feature_nt_length"] == g["length_bp"]
    assert result["translated_nt_length"] == 9
    assert result["recomputed_translation"] == "MK"


@pytest.mark.parametrize(
    "qualifiers,reason",
    [
        ({"pseudo": [""]}, "pseudo"),
        ({"exception": ["ribosomal slippage"]}, "exception"),
        ({"ribosomal_slippage": [""]}, "exception"),
        ({"transl_except": ["(pos:4..6,aa:Sec)"]}, "transl_except"),
    ],
)
def test_exceptional_features_are_retained(qualifiers, reason):
    feature = SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS", qualifiers=qualifiers)
    result = normalize_cds(genome("complete"), feature)
    assert result["eligible_primary"] is False
    assert reason in result["primary_exclusion_reasons"]
    assert json.loads(result["qualifiers_json"]) == qualifiers


def test_order_and_noncontiguous_join_excluded():
    feature = SeqFeature(
        CompoundLocation([SimpleLocation(0, 3, strand=1), SimpleLocation(6, 9, strand=1)], "order"),
        type="CDS",
    )
    result = normalize_cds(genome("complete"), feature)
    assert result["location_operator"] == "order"
    assert "order_location" in result["primary_exclusion_reasons"]
    assert "noncontiguous_join" in result["primary_exclusion_reasons"]


def test_same_geometry_sources_share_model_but_not_feature_ids():
    feature = SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS")
    a = normalize_cds(genome("complete"), feature, "genbank")
    b = normalize_cds(genome("complete"), feature, "phanotate")
    assert a["model_id"] == b["model_id"]
    assert a["feature_id"] != b["feature_id"]


def test_feature_and_call_provenance_separates_source_runs():
    feature = SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS")
    a = normalize_cds(genome("complete"), feature, "phanotate", "synthetic_run_r0")
    b = normalize_cds(genome("complete"), feature, "phanotate", "synthetic_run_rL2")
    assert a["model_id"] == b["model_id"]
    assert a["feature_id"] != b["feature_id"]
    assert a["call_id"] != b["call_id"]


def test_unknown_location_is_preserved_and_excluded():
    result = normalize_cds(genome("complete"), SeqFeature(type="CDS"))
    assert result["segments"] == []
    assert result["eligible_primary"] is False
    assert result["primary_exclusion_reason"] == "missing_location"
