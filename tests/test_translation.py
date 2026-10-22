import hashlib
import json
from pathlib import Path

import pytest

from phageannoscope.translation import translate_cds

SEQUENCES = json.loads(
    (Path(__file__).parent / "fixtures" / "synthetic" / "sequences.json").read_text()
)


@pytest.mark.parametrize(
    "key,codon_start", [("complete", 1), ("codon_start_2", 2), ("codon_start_3", 3)]
)
def test_offsets_and_terminal_stop_rule(key, codon_start):
    result = translate_cds(SEQUENCES[key], codon_start=codon_start, provided_translation="MK*")
    assert result["recomputed_translation_raw"] == "MK*"
    assert result["recomputed_translation"] == "MK"
    assert result["provided_translation"] == "MK*"
    assert result["translated_nt_length"] == 9
    assert result["translation_validity"] == "valid"
    assert result["protein_sha256"] == hashlib.sha256(b"MK").hexdigest()


def test_genbank_default_one_and_table_eleven_initiation():
    one = translate_cds(SEQUENCES["gtg"], transl_table=1)
    eleven = translate_cds(SEQUENCES["gtg"], transl_table=11)
    assert one["recomputed_translation"] == "VK"
    assert one["translation_validity"] == "invalid_start"
    assert eleven["recomputed_translation"] == "MK"
    assert eleven["translation_notes"] == ["alternative_initiation_met"]


def test_code_fifteen_tag_is_glutamine():
    alt = translate_cds(SEQUENCES["table15"], transl_table=15)
    assert alt["recomputed_translation"] == "MQ"
    assert alt["translation_validity"] == "valid"
    standard = translate_cds(SEQUENCES["table15"], transl_table=1)
    assert standard["recomputed_translation"] == "M*"
    assert standard["translation_validity"] == "internal_stop"
    assert standard["protein_sha256"] is None


def test_provided_disagreement_does_not_destroy_valid_recomputation():
    result = translate_cds(SEQUENCES["complete"], provided_translation="MM")
    assert result["provided_translation"] == "MM"
    assert result["recomputed_translation"] == "MK"
    assert result["translation_validity"] == "provided_translation_mismatch"
    assert result["protein_validation"] == "valid"
    assert result["protein_sha256"]


@pytest.mark.parametrize(
    "key,kwargs,reason",
    [
        ("invalid_symbol", {}, "invalid_nucleotide"),
        ("ambiguous", {}, "ambiguous_translation"),
        ("internal_stop", {}, "internal_stop"),
        ("complete", {"pseudo": True}, "pseudo"),
        ("complete", {"transl_except": ["synthetic exception"]}, "transl_except"),
        ("complete", {"exception": "ribosomal slippage"}, "exception"),
        ("partial", {"partial_5p": True, "partial_3p": True}, "partial"),
        ("complete", {"codon_start": 0}, "invalid_codon_start"),
        ("complete", {"transl_table": 9999}, "unsupported_table"),
    ],
)
def test_exceptional_sequences_retained_without_validated_protein(key, kwargs, reason):
    result = translate_cds(SEQUENCES[key], **kwargs)
    assert result["translation_validity"] == reason
    assert result["protein_sha256"] is None
