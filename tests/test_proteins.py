"""Synthetic protein provenance tests, independent from caller predictions."""

import hashlib

import pytest

from phageannoscope.proteins import build_catalogue, write_fasta


def origin(identifier, sequence="MKM", **extra):
    return {
        "call_id": identifier,
        "model_id": "model_" + identifier,
        "source_run_id": "run_" + identifier,
        "genome_id": "synthetic_" + identifier,
        "source": "genbank",
        "protein_id": "possibly-reused-id",
        "protein_validation": "valid",
        "translation_validity": "valid",
        "recomputed_translation": sequence,
        "provided_translation": sequence,
        "protein_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        **extra,
    }


def test_identical_proteins_deduplicate_without_losing_orfs_or_genomes():
    calls = [origin("a"), origin("b"), origin("c", "MM")]
    catalogue = build_catalogue(calls)
    assert len(catalogue["proteins"]) == 2
    assert len(catalogue["protein_origins"]) == 3
    assert catalogue == build_catalogue(reversed(calls))
    digest = hashlib.sha256(b"MKM").hexdigest()
    assert {
        r["genome_id"] for r in catalogue["protein_origins"] if r["protein_sha256"] == digest
    } == {"synthetic_a", "synthetic_b"}


def test_invalid_provided_translation_never_becomes_catalogue_protein():
    call = origin(
        "a", provided_translation="M*BAD", translation_validity="provided_translation_mismatch"
    )
    result = build_catalogue([call])
    assert result["proteins"][0]["amino_acid_sequence"] == "MKM"
    assert result["protein_origins"][0]["provided_translation_matches"] is False
    assert result["protein_origins"][0]["translation_provenance"] == "verified_recomputed"
    call["protein_validation"] = "invalid"
    call["translation_validity"] = "internal_stop"
    assert build_catalogue([call]) == {"proteins": [], "protein_origins": []}


@pytest.mark.parametrize("sequence", ["M*KM", "MKX", "MKM*", "mkm", ""])
def test_bad_sequence_cannot_be_marked_verified(sequence):
    with pytest.raises(ValueError, match="invalid amino acids"):
        build_catalogue([origin("a", sequence)])


def test_incorrect_digest_and_duplicate_origin_are_rejected():
    with pytest.raises(ValueError, match="SHA mismatch"):
        build_catalogue([origin("a", protein_sha256="incorrect")])
    with pytest.raises(ValueError, match="Duplicate protein origin"):
        build_catalogue([origin("a"), origin("a")])


def test_fasta_uses_sha_headers_and_verified_sequences(tmp_path):
    result = build_catalogue([origin("a")])
    path = write_fasta(result, tmp_path / "proteins.faa", line_width=2)
    digest = hashlib.sha256(b"MKM").hexdigest()
    assert path.read_text() == f">{digest}\nMK\nM\n"
    result["proteins"][0]["amino_acid_sequence"] = "BAD*"
    with pytest.raises(ValueError, match="unverified"):
        write_fasta(result, path)
