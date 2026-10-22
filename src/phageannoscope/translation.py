"""Explicit CDS translation, with independently retained annotation and QC."""

from __future__ import annotations

import hashlib

from Bio.Data import CodonTable
from Bio.Seq import Seq

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
TRANSLATION_NORMALIZATION = "terminal-stop-only;complete-initiation-met;v1"


def translate_cds(
    nucleotide,
    *,
    codon_start=1,
    transl_table=1,
    provided_translation=None,
    partial_5p=False,
    partial_3p=False,
    pseudo=False,
    exception=None,
    transl_except=None,
):
    """Recompute without deleting internal stops or replacing invalid bases.

    ``protein_validation`` concerns the computed sequence; a disagreement
    with a provided annotation is separately visible in translation_validity.
    Partial/exceptional models are retained but not admitted to the catalogue.
    """
    nucleotide = str(nucleotide).upper()
    provided = (
        "".join(str(provided_translation).split()) if provided_translation is not None else None
    )
    result = {
        "provided_translation": provided,
        "recomputed_translation": None,
        "recomputed_translation_raw": None,
        "protein_sha256": None,
        "nucleotide_sha256": hashlib.sha256(nucleotide.encode("utf-8")).hexdigest(),
        "translation_validity": "valid",
        "protein_validation": "invalid",
        "translation_notes": [],
        "translated_nt_length": 0,
        "start_codon": None,
        "stop_codon": None,
        "terminal_stop_present": False,
        "translation_normalization": TRANSLATION_NORMALIZATION,
    }
    reasons = []
    if codon_start not in (1, 2, 3):
        result["translation_validity"] = "invalid_codon_start"
        return result
    try:
        table = CodonTable.unambiguous_dna_by_id[int(transl_table)]
    except (KeyError, TypeError, ValueError):
        result["translation_validity"] = "unsupported_table"
        return result
    coding = nucleotide[codon_start - 1 :]
    if set(coding) - IUPAC_DNA:
        result["translation_validity"] = "invalid_nucleotide"
        return result
    complete_length = len(coding) - len(coding) % 3
    result["translated_nt_length"] = complete_length
    result["start_codon"] = coding[:3] if len(coding) >= 3 else None
    result["stop_codon"] = coding[-3:] if len(coding) >= 3 and len(coding) % 3 == 0 else None
    result["terminal_stop_present"] = result["stop_codon"] in table.stop_codons
    if not coding or len(coding) % 3:
        reasons.append("invalid_length")
    raw_protein = (
        str(Seq(coding[:complete_length]).translate(table=int(transl_table)))
        if complete_length
        else ""
    )
    if not partial_5p and result["start_codon"] in table.start_codons and raw_protein:
        if raw_protein[0] != "M":
            result["translation_notes"].append("alternative_initiation_met")
        raw_protein = "M" + raw_protein[1:]
    protein = raw_protein[:-1] if raw_protein.endswith("*") else raw_protein
    result["recomputed_translation_raw"] = raw_protein
    result["recomputed_translation"] = protein
    if "*" in protein:
        reasons.append("internal_stop")
    if "X" in protein:
        reasons.append("ambiguous_translation")
    if pseudo:
        reasons.insert(0, "pseudo")
    if exception:
        reasons.insert(0, "exception")
    if transl_except:
        reasons.insert(0, "transl_except")
    if partial_5p or partial_3p:
        reasons.insert(0, "partial")
    if not partial_5p and result["start_codon"] not in table.start_codons:
        reasons.append("invalid_start")
    if not partial_3p and not result["terminal_stop_present"]:
        reasons.append("missing_terminal_stop")
    # A trailing stop in an annotation is normalized by this one explicit rule.
    comparison = provided[:-1] if provided and provided.endswith("*") else provided
    if provided is not None and comparison != protein:
        result["translation_notes"].append("provided_translation_disagrees")
    if reasons:
        result["translation_validity"] = reasons[0]
    else:
        result["protein_validation"] = "valid"
        result["translation_validity"] = (
            "provided_translation_mismatch"
            if provided is not None and comparison != protein
            else "valid"
        )
        result["protein_sha256"] = hashlib.sha256(protein.encode("ascii")).hexdigest()
    result["translation_issues"] = reasons
    return result
