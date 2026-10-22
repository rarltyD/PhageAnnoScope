"""Synthetic adapter checks; no PHANOTATE executable is substituted or faked."""

import json
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phageannoscope.callers.phanotate import parse_gff
from phageannoscope.callers.pyrodigal_gv import normalize_gene, predict
from phageannoscope.genbank_io import genome_from_record
from phageannoscope.locations import extract_segments

FIXTURES = Path(__file__).parent / "fixtures" / "synthetic"
SEQUENCES = json.loads((FIXTURES / "sequences.json").read_text())


def synthetic_genome(key, topology="linear"):
    record = SeqRecord(Seq(SEQUENCES[key]), id="SYNTHETIC_MAPPING.2")
    record.annotations["topology"] = topology
    return genome_from_record(record)


def artificial_gff(start, end, strand="+", phase="0", attributes="ID=SYNTHETIC_CALL"):
    """Hand-authored GFF fixture, explicitly independent of a predictor run."""
    return f"SYNTHETIC_MAPPING.2\tsynthetic_adapter_test\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attributes}\n"


@pytest.mark.parametrize(
    "key,strand,stop", [("complete", "+", [6, 7, 8]), ("reverse", "-", [2, 1, 0])]
)
def test_hand_gff_includes_stop_and_reverse_complement(key, strand, stop):
    call = parse_gff(artificial_gff(1, 9, strand), synthetic_genome(key), "synthetic_gff_run")[0]
    assert call["segments"] == [{"start0": 0, "end0": 9}]
    assert call["feature_nt_length"] == 9
    assert call["recomputed_translation"] == "MK"
    assert call["stop_triplet"] == stop
    assert call["terminal_stop_present"]
    assert call["eligible_primary"]
    assert call["raw_feature_id"] == "SYNTHETIC_CALL"


@pytest.mark.parametrize(
    "key,strand,start,end,expected,stop",
    [
        (
            "origin_plus",
            "+",
            4,
            12,
            [{"start0": 9, "end0": 12}, {"start0": 0, "end0": 6}],
            [3, 4, 5],
        ),
        (
            "origin_minus",
            "-",
            1,
            9,
            [{"start0": 0, "end0": 3}, {"start0": 6, "end0": 12}],
            [8, 7, 6],
        ),
    ],
)
def test_hand_gff_rotation_back_mapping(key, strand, start, end, expected, stop):
    call = parse_gff(
        artificial_gff(start, end, strand),
        synthetic_genome(key, "circular"),
        "synthetic_rotated",
        6,
    )[0]
    assert call["segments"] == expected
    assert call["recomputed_translation"] == "MK"
    assert call["stop_triplet"] == stop
    assert call["wraps_origin"]
    assert not call["partial_5p"] and not call["partial_3p"]
    assert not call["eligible_primary"]


@pytest.mark.parametrize("key,phase,end", [("codon_start_2", "1", 10), ("codon_start_3", "2", 11)])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_numeric_gff_phase_is_offset_from_transcript_start(key, phase, end, strand):
    genome = synthetic_genome(key)
    if strand == "-":
        genome = genome_from_record(
            SeqRecord(Seq(genome["sequence"]).reverse_complement(), id="SYNTHETIC_MAPPING.2")
        )
    call = parse_gff(artificial_gff(1, end, strand=strand, phase=phase), genome, "synthetic_phase")[
        0
    ]
    assert call["codon_start"] == int(phase) + 1
    assert call["recomputed_translation"] == "MK"
    assert call["translated_nt_length"] == 9


def test_unknown_phase_is_preserved_and_never_assumed_zero():
    # This internal interval has neither a recognized start nor a stop and
    # cannot be assigned a frame from the omitted GFF phase.
    call = parse_gff(
        artificial_gff(13, 18, phase="."),
        synthetic_genome("multi_feature_genome"),
        "synthetic_unknown",
    )[0]
    assert call["gff_phase_raw"] == "."
    assert call["codon_start"] is None
    assert call["translation_validity"] == "invalid_codon_start"
    assert call["effective_phase"] is None
    assert call["phase_resolution"] == "unresolved_missing_phase"
    assert not call["eligible_primary"]


@pytest.mark.parametrize(
    "key,strand,stop", [("complete", "+", [6, 7, 8]), ("reverse", "-", [2, 1, 0])]
)
def test_missing_gff_phase_resolved_only_from_complete_nucleotide_evidence(key, strand, stop):
    call = parse_gff(
        artificial_gff(1, 9, strand=strand, phase="."),
        synthetic_genome(key),
        "synthetic_missing_phase",
    )[0]
    assert call["gff_phase_raw"] == "."
    assert call["effective_phase"] == 0
    assert call["phase_resolution"] == "inferred_zero_from_complete_start_stop_and_codon_length"
    assert call["codon_start"] == 1
    assert call["recomputed_translation"] == "MK"
    assert call["stop_triplet"] == stop
    assert call["protein_validation"] == "valid"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("partial_end", ["5p", "3p"])
def test_missing_phase_edge_partial_requires_one_valid_codon_anchor(strand, partial_end):
    sequence = SEQUENCES["complete"][3:] if partial_end == "5p" else SEQUENCES["complete"][:-3]
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    genome = genome_from_record(SeqRecord(Seq(sequence), id="SYNTHETIC_MAPPING.2"))
    call = parse_gff(
        artificial_gff(1, 6, strand=strand, phase="."), genome, "synthetic_edge_phase"
    )[0]
    assert call["effective_phase"] == 0
    assert call["partial_5p"] is (partial_end == "5p")
    assert call["partial_3p"] is (partial_end == "3p")
    expected = (
        "inferred_zero_from_terminal_stop_and_5p_input_edge"
        if partial_end == "5p"
        else "inferred_zero_from_start_codon_and_3p_input_edge"
    )
    assert call["phase_resolution"] == expected
    assert not call["eligible_primary"]


def test_missing_phase_is_not_inferred_from_unanchored_edge_or_internal_incomplete_cds():
    unanchored = parse_gff(
        artificial_gff(1, 3, phase="."), synthetic_genome("partial"), "synthetic_unanchored"
    )[0]
    incomplete = parse_gff(
        artificial_gff(1, 6, phase="."), synthetic_genome("complete"), "synthetic_incomplete"
    )[0]
    assert unanchored["effective_phase"] is None
    assert incomplete["effective_phase"] is None
    assert not unanchored["eligible_primary"] and not incomplete["eligible_primary"]


@pytest.mark.parametrize("strand", [".", "?"])
def test_unknown_strand_kept_for_review(strand):
    call = parse_gff(
        artificial_gff(1, 9, strand=strand), synthetic_genome("complete"), "synthetic_unknown"
    )[0]
    assert call["strand"] is None
    assert call["segments"] == [{"start0": 0, "end0": 9}]
    assert not call["eligible_primary"]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_missing_codon_at_input_edges_is_partial(strand):
    call = parse_gff(
        artificial_gff(1, 3, strand=strand), synthetic_genome("partial"), "synthetic_edge"
    )[0]
    assert call["partial_5p"] and call["partial_3p"]
    assert not call["eligible_primary"]


def test_parser_rejects_bad_phase_and_bounds_but_accepts_successful_zero():
    genome = synthetic_genome("complete")
    assert (
        parse_gff(
            "##gff-version 3\n# SYNTHETIC successful zero-CDS output\n", genome, "synthetic_zero"
        )
        == []
    )
    with pytest.raises(ValueError, match="strand/phase"):
        parse_gff(artificial_gff(1, 9, phase="4"), genome, "synthetic_bad")
    with pytest.raises(ValueError, match="outside genome"):
        parse_gff(artificial_gff(1, 10), genome, "synthetic_bad")
    with pytest.raises(ValueError, match="identifier"):
        parse_gff(
            artificial_gff(1, 9).replace("SYNTHETIC_MAPPING.2", "SYNTHETIC_OTHER.2"),
            genome,
            "synthetic_bad",
        )


def derived_artificial_sequence():
    """Expand codons from labelled fixtures, solely for a native API test."""
    coding = (
        SEQUENCES["complete"][:3] + SEQUENCES["complete"][3:6] * 100 + SEQUENCES["complete"][-3:]
    )
    return coding + SEQUENCES["multi_feature_genome"] * 20 + str(Seq(coding).reverse_complement())


@pytest.mark.parametrize("shift", [0, 100, 401])
def test_installed_native_pyrodigal_coordinates_stop_and_partial_flags(shift):
    pyrodigal_gv = pytest.importorskip("pyrodigal_gv")
    sequence = derived_artificial_sequence()
    record = SeqRecord(Seq(sequence), id="SYNTHETIC_NATIVE.1")
    record.annotations["topology"] = "circular"
    genome = genome_from_record(record)
    rotated = sequence[shift:] + sequence[:shift]
    native = pyrodigal_gv.ViralGeneFinder(meta=True).find_genes(rotated.encode("ascii"))
    assert len(native) > 0  # The actual installed algorithm is exercised.
    for ordinal, gene in enumerate(native):
        call = normalize_gene(genome, gene, "synthetic_native", ordinal, shift)
        assert call["feature_nt_length"] == gene.end - gene.begin + 1
        assert extract_segments(sequence, call["segments"], call["strand"]) == gene.sequence()
        assert call["provided_translation"] == gene.translate(include_stop=True)
        assert call["recomputed_translation_raw"] == gene.translate(include_stop=True)
        assert call["partial_5p"] == (gene.partial_begin if gene.strand == 1 else gene.partial_end)
        assert call["partial_3p"] == (gene.partial_end if gene.strand == 1 else gene.partial_begin)
        if not call["partial_3p"]:
            assert gene.translate().endswith("*")
            assert call["terminal_stop_present"]
            assert call["stop_triplet"] is not None


def test_installed_native_parameters_recorded_and_single_mode_rejected():
    pytest.importorskip("pyrodigal_gv")
    result = predict(synthetic_genome("complete"), "synthetic_native_zero")
    assert result["calls"] == []
    assert result["actual_params"]["meta"] is True
    assert result["actual_params"]["min_gene"] == 90
    assert result["actual_params"]["min_edge_gene"] == 60
    assert result["actual_params"]["max_overlap"] == 60
    assert result["actual_params"]["min_mask"] == 50
    with pytest.raises(ValueError, match="metagenomic"):
        predict(synthetic_genome("complete"), "synthetic_bad_mode", params={"meta": False})
