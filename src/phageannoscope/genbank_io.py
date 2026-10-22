"""Lossless feature import with an explicit, conservative CDS analysis layer."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

from .locations import coalesce_segments, extract_segments, stop_positions
from .translation import translate_cds


def stable_id(prefix, value):
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, default=str
    )
    return prefix + "_" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def genome_from_record(record):
    sequence = str(record.seq).upper()
    accession_version = record.id
    accession = re.sub(r"\.\d+$", "", accession_version)
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    return {
        "genome_id": stable_id(
            "genome", {"accession_version": accession_version, "sequence_sha256": digest}
        ),
        "accession": accession,
        "accession_version": accession_version,
        "sequence_sha256": digest,
        "length_bp": len(sequence),
        "sequence": sequence,
        "topology": str(record.annotations.get("topology", "unknown")).lower(),
        "organism": record.annotations.get("organism"),
        "annotations_json": json.dumps(record.annotations, sort_keys=True, default=str),
        "source_qualifiers_json": json.dumps(
            [f.qualifiers for f in record.features if f.type == "source"],
            sort_keys=True,
            default=str,
        ),
    }


def _first(qualifiers, name, default=None):
    value = qualifiers.get(name, default)
    if isinstance(value, (list, tuple)):
        return value[0] if value else default
    return value


def _location(feature, length):
    location = feature.location
    if location is None:
        return [], None, None, False, False, ["missing_location"]
    parts = list(location.parts)
    strand = location.strand
    issues = []
    segments = []
    for part in parts:
        if part.ref or part.ref_db:
            issues.append("remote_location")
        try:
            start, end = int(part.start), int(part.end)
        except (TypeError, ValueError):
            issues.append("unknown_location_position")
            continue
        if not (0 <= start < end <= length):
            issues.append("invalid_location_bounds")
        else:
            segments.append({"start0": start, "end0": end})
    if strand not in (1, -1) or any(p.strand != strand for p in parts):
        issues.append("unknown_or_mixed_strand")
    # Fuzziness follows the transcript endpoints, not genomic left/right.
    first_boundary = parts[0].start if strand != -1 else parts[0].end
    last_boundary = parts[-1].end if strand != -1 else parts[-1].start
    partial_5p = type(first_boundary) is not ExactPosition
    partial_3p = type(last_boundary) is not ExactPosition
    internal_boundaries = []
    for index, part in enumerate(parts):
        boundaries = (part.start, part.end) if strand != -1 else (part.end, part.start)
        if index > 0:
            internal_boundaries.append(boundaries[0])
        if index < len(parts) - 1:
            internal_boundaries.append(boundaries[1])
    if any(type(boundary) is not ExactPosition for boundary in internal_boundaries):
        issues.append("internal_fuzzy_position")
    if len({(s["start0"], s["end0"]) for s in segments}) != len(segments):
        issues.append("duplicate_segments")
    intervals = sorted((s["start0"], s["end0"]) for s in segments)
    if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
        issues.append("overlapping_segments")
    return segments, strand, getattr(location, "operator", None), partial_5p, partial_3p, issues


def _feature_row(genome, feature, source, source_run_id, ordinal):
    segments, strand, operator, partial5, partial3, issues = _location(feature, genome["length_bp"])
    original_location = str(feature.location) if feature.location is not None else None
    qualifiers_json = json.dumps(
        feature.qualifiers, sort_keys=True, ensure_ascii=False, default=str
    )
    feature_id = stable_id(
        "feature",
        {
            "sequence_sha256": genome["sequence_sha256"],
            "accession_version": genome["accession_version"],
            "source": source,
            "source_run_id": source_run_id,
            "ordinal": ordinal,
            "location": original_location,
            "type": feature.type,
            "qualifiers": feature.qualifiers,
        },
    )
    return {
        "feature_id": feature_id,
        "genome_id": genome["genome_id"],
        "accession_version": genome["accession_version"],
        "source": source,
        "source_run_id": source_run_id,
        "feature_ordinal": ordinal,
        "feature_type": feature.type,
        "type": feature.type,
        "raw_feature_id": getattr(feature, "id", None),
        "strand": strand,
        "segments": segments,
        "original_location": original_location,
        "location_operator": operator,
        "qualifiers_json": qualifiers_json,
        "partial_5p": partial5,
        "partial_3p": partial3,
        "location_issues": issues,
    }


def normalize_cds(
    genome, feature, source="genbank", source_run_id=None, ordinal=0, rotation_shift=0
):
    """Normalize a Bio.SeqFeature CDS already mapped to original coordinates.

    Raw feature ordinals are record-wide, including non-CDS features. IDs never
    depend on protein_id and identical reruns produce identical identifiers.
    """
    row = _feature_row(genome, feature, source, source_run_id, ordinal)
    q = feature.qualifiers
    codon_raw = _first(q, "codon_start", "1")
    table_raw = _first(q, "transl_table")
    try:
        codon_start = int(codon_raw)
    except (ValueError, TypeError):
        codon_start = None
    try:
        table = int(table_raw) if table_raw is not None else 1
    except (ValueError, TypeError):
        table = None
    code_class = (
        "standard_stop"
        if table in (1, 11)
        else (f"ncbi_{table}" if table is not None else "unknown")
    )
    segments, strand = row["segments"], row["strand"]
    length = genome["length_bp"]
    canonical_segments = coalesce_segments(segments, strand)
    wraps = any(
        (strand == 1 and left["end0"] == length and right["start0"] == 0)
        or (strand == -1 and left["start0"] == 0 and right["end0"] == length)
        for left, right in zip(canonical_segments, canonical_segments[1:])
    )
    contiguous = len(canonical_segments) == 1 or (
        len(canonical_segments) == 2 and wraps and genome.get("topology") == "circular"
    )
    pseudo = "pseudo" in q or "pseudogene" in q
    exception = _first(q, "exception")
    if exception is None and "ribosomal_slippage" in q:
        # Programmed frameshift: same exclusion class as /exception (complex biology).
        exception = "ribosomal_slippage"
    transl_except = q.get("transl_except") or None
    nucleotide = ""
    if not row["location_issues"] and strand in (1, -1):
        nucleotide = extract_segments(genome["sequence"], segments, strand)
    translation = translate_cds(
        nucleotide,
        codon_start=codon_start,
        transl_table=table,
        provided_translation=_first(q, "translation"),
        partial_5p=row["partial_5p"],
        partial_3p=row["partial_3p"],
        pseudo=pseudo,
        exception=exception,
        transl_except=transl_except,
    )
    exclusions = list(row["location_issues"])
    if row["partial_5p"]:
        exclusions.append("partial_5p")
    if row["partial_3p"]:
        exclusions.append("partial_3p")
    if pseudo:
        exclusions.append("pseudo")
    if exception:
        exclusions.append("exception")
    if transl_except:
        exclusions.append("transl_except")
    if row["location_operator"] == "order":
        exclusions.append("order_location")
    if not contiguous:
        exclusions.append("noncontiguous_join")
    if wraps:
        exclusions.append("wraps_origin")
    if translation["protein_validation"] != "valid":
        exclusions.append(translation["translation_validity"])
    # Preserve original segmentation on calls; canonical geometry joins only
    # adjacent pieces so equivalent complete coordinate models share an ID.
    model_id = stable_id(
        "model",
        {
            "sequence_sha256": genome["sequence_sha256"],
            "strand": strand,
            "segments": canonical_segments,
            "codon_start": codon_start,
            "code_class": code_class,
        },
    )
    row.update(translation)
    row.update(
        {
            "call_id": stable_id(
                "call",
                {
                    "feature_id": row["feature_id"],
                    "source_run_id": source_run_id,
                    "rotation_shift": rotation_shift,
                },
            ),
            "model_id": model_id,
            "sequence_sha256": genome["sequence_sha256"],
            "algorithm_family": source
            if source in ("genbank", "phanotate", "pyrodigal_gv")
            else None,
            "protein_id": _first(q, "protein_id"),
            "locus_tag": _first(q, "locus_tag"),
            "codon_start": codon_start,
            "codon_start_raw": str(codon_raw),
            "transl_table_raw": str(table_raw) if table_raw is not None else None,
            "transl_table_effective": table,
            "transl_table_qualifier_missing": table_raw is None,
            "qualifier_missing": table_raw is None,
            "code_class": code_class,
            "feature_nt_length": sum(s["end0"] - s["start0"] for s in segments),
            "pseudo": pseudo,
            "exception": exception,
            "transl_except": transl_except,
            "wraps_origin": wraps,
            "segments_contiguous": contiguous,
            "rotation_shift": int(rotation_shift),
            "rotation_sensitive": False,
            "eligible_primary": not exclusions,
            "primary_exclusion_reason": exclusions[0] if exclusions else None,
            "primary_exclusion_reasons": list(dict.fromkeys(exclusions)),
            "stop_triplet": stop_positions(segments, strand, codon_start)
            if (
                translation["terminal_stop_present"]
                and not row["partial_3p"]
                and not exception
                and not transl_except
            )
            else None,
        }
    )
    return row


def record_to_tables(record, source="genbank", source_run_id=None):
    genome = genome_from_record(record)
    result = {"genomes": [genome], "features": [], "cds_calls": [], "qc_events": []}
    for ordinal, feature in enumerate(record.features):
        result["features"].append(_feature_row(genome, feature, source, source_run_id, ordinal))
        if feature.type != "CDS":
            continue
        call = normalize_cds(genome, feature, source, source_run_id, ordinal)
        result["cds_calls"].append(call)
        if not call["eligible_primary"] or call["translation_validity"] != "valid":
            result["qc_events"].append(
                {
                    "qc_id": stable_id(
                        "qc",
                        {
                            "call_id": call["call_id"],
                            "reason": call["primary_exclusion_reasons"],
                            "translation": call["translation_validity"],
                        },
                    ),
                    "severity": "warning",
                    "code": "cds_requires_review",
                    "genome_id": genome["genome_id"],
                    "source": source,
                    "source_run_id": source_run_id,
                    "call_id": call["call_id"],
                    "message": ";".join(call["primary_exclusion_reasons"])
                    or call["translation_validity"],
                    "evidence_json": json.dumps(
                        {
                            "original_location": call["original_location"],
                            "translation_validity": call["translation_validity"],
                        }
                    ),
                }
            )
    return result


def read_genbank(path, source="genbank", source_run_id=None, accession_version=None):
    result = {"genomes": [], "features": [], "cds_calls": [], "qc_events": []}
    seen = {}
    for record in SeqIO.parse(str(Path(path)), "genbank"):
        if accession_version is not None and record.id != accession_version:
            continue
        tables = record_to_tables(record, source, source_run_id)
        genome = tables["genomes"][0]
        key = genome["accession_version"]
        if key in seen:
            raise ValueError(f"Duplicate GenBank accession.version: {key}")
        seen[key] = genome["sequence_sha256"]
        for name in result:
            result[name].extend(tables[name])
    if not result["genomes"]:
        raise ValueError(f"No GenBank records found in {path}")
    return result
