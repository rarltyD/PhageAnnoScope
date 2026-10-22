"""A verified amino-acid catalogue with lossless translation provenance."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable

NORMALIZATION_VERSION = "terminal-stop-only;complete-initiation-met;v1"
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def build_catalogue(calls: Iterable[dict]) -> dict[str, list[dict]]:
    """Deduplicate verified recomputations, never a dubious provided protein.

    Invalid translations remain in ``cds_calls`` for inspection.  A disagreement
    in an existing /translation does not invalidate a separately verified
    recomputation; that disagreement is explicit on the provenance row.
    """
    proteins: dict[str, dict] = {}
    origins = []
    seen: set[str] = set()
    for call in sorted(calls, key=lambda row: str(row.get("call_id"))):
        valid = call.get("protein_validation") == "valid"
        if "protein_validation" not in call:
            valid = call.get("translation_validity") in {
                "valid",
                "provided_translation_mismatch",
            } and bool(call.get("protein_sha256"))
        if not valid:
            continue
        sequence = call.get("recomputed_translation")
        if (
            not isinstance(sequence, str)
            or not sequence
            or any(aa not in AMINO_ACIDS for aa in sequence)
        ):
            raise ValueError(
                f"Verified protein contains invalid amino acids: {call.get('call_id')}"
            )
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if call.get("protein_sha256") and call["protein_sha256"] != digest:
            raise ValueError(f"Recomputed protein SHA mismatch: {call.get('call_id')}")
        call_id = str(call["call_id"])
        if call_id in seen:
            raise ValueError(f"Duplicate protein origin call_id {call_id}")
        seen.add(call_id)
        proteins[digest] = {
            "protein_sha256": digest,
            "amino_acid_sequence": sequence,
            "length_aa": len(sequence),
            "normalization_version": NORMALIZATION_VERSION,
            "validation": "valid",
            "translation_provenance": "recomputed",
        }
        provided = call.get("provided_translation")
        origins.append(
            {
                "protein_sha256": digest,
                "model_id": call.get("model_id"),
                "call_id": call["call_id"],
                "genome_id": call["genome_id"],
                "source_run_id": call["source_run_id"],
                "feature_id": call.get("feature_id"),
                "source": call.get("source"),
                "protein_id": call.get("protein_id"),
                "translation_provenance": "verified_recomputed",
                "translation_validity": call.get("translation_validity"),
                "provided_translation_present": provided is not None,
                "provided_translation_matches": (
                    provided[:-1] if provided.endswith("*") else provided
                )
                == sequence
                if isinstance(provided, str)
                else None,
            }
        )
    return {"proteins": [proteins[k] for k in sorted(proteins)], "protein_origins": origins}


def write_fasta(proteins: Iterable[dict] | dict, path: str | Path, line_width: int = 80) -> Path:
    """Write one verified sequence per SHA header, in deterministic order."""
    if line_width < 1:
        raise ValueError("line_width must be positive")
    if isinstance(proteins, dict):
        proteins = proteins["proteins"]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row in sorted(proteins, key=lambda r: r["protein_sha256"]):
            sequence = row["amino_acid_sequence"]
            digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            if (
                digest != row["protein_sha256"]
                or not sequence
                or any(aa not in AMINO_ACIDS for aa in sequence)
            ):
                raise ValueError("Protein catalogue contains an unverified sequence")
            handle.write(">" + digest + "\n")
            for start in range(0, len(sequence), line_width):
                handle.write(sequence[start : start + line_width] + "\n")
    return path
