"""PHANOTATE GFF import with per-call stop and phase validation.

The parser never extends an interval to manufacture a stop. Certification of
an installed PHANOTATE version's output remains a separate integration check.
"""

import importlib.metadata
import re
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

from Bio.Data import CodonTable
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    SeqFeature,
    SimpleLocation,
)

from ..genbank_io import normalize_cds
from ..locations import extract_segments, map_rotated_segments


def command(configured=None):
    if configured:
        return list(configured) if isinstance(configured, list) else [configured]
    candidate = Path(sys.executable).parent / "phanotate.py"
    found = str(candidate) if candidate.exists() else shutil.which("phanotate.py")
    if found:
        return [sys.executable, found] if found.endswith(".py") else [found]
    return None


def command_version(config=None):
    """Observe the configured executable's version, including an external WSL tool."""
    config = config or {}

    def clean(value):
        value = str(value).strip() if value is not None else ""
        return value if re.fullmatch(r"\d+(?:\.\d+)+(?:[A-Za-z0-9.+_-]*)", value) else None

    if config.get("phanotate_observed_version"):
        return clean(config["phanotate_observed_version"])
    if config.get("phanotate_version_command"):
        argv = config["phanotate_version_command"]
        if (
            not isinstance(argv, list)
            or not argv
            or not all(isinstance(item, str) for item in argv)
        ):
            return None
        try:
            result = subprocess.run(argv, capture_output=True, text=True, timeout=30, check=False)
        except (OSError, subprocess.SubprocessError):
            return None
        return clean(result.stdout) if result.returncode == 0 else None
    try:
        return clean(importlib.metadata.version("phanotate"))
    except importlib.metadata.PackageNotFoundError:
        return None


def parse_gff(text, genome, source_run_id, rotation_shift=0):
    calls = []
    for line in text.splitlines():
        if line.strip() == "##FASTA":
            break
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"Malformed PHANOTATE GFF line: {line[:150]}")
        seqid, _, kind, left, right, score, strand_text, phase, attributes = fields
        if kind != "CDS":
            continue
        if seqid != genome["accession_version"]:
            raise ValueError("GFF sequence identifier does not match submitted input")
        if strand_text not in ("+", "-", ".", "?") or phase not in ("0", "1", "2", "."):
            raise ValueError("Invalid PHANOTATE strand/phase")
        start0, end0 = int(left) - 1, int(right)
        if not 0 <= start0 < end0 <= genome["length_bp"]:
            raise ValueError("PHANOTATE coordinates outside genome")
        strand = {"+": 1, "-": -1}.get(strand_text)
        segments = map_rotated_segments(
            [{"start0": start0, "end0": end0}], rotation_shift, genome["length_bp"], strand or 1
        )
        parts = [SimpleLocation(s["start0"], s["end0"], strand=strand) for s in segments]
        # PHANOTATE can emit edge ORFs without explicit partial qualifiers.
        # Only infer a missing biological endpoint when its codon is absent
        # and the predicted endpoint is within two bases of the input edge.
        partial5 = partial3 = False
        effective_phase = int(phase) if phase != "." else None
        phase_resolution = (
            "provided_gff_phase" if effective_phase is not None else "unresolved_missing_phase"
        )
        if strand in (1, -1):
            shift = rotation_shift % genome["length_bp"]
            rotated = genome["sequence"][shift:] + genome["sequence"][:shift]
            nucleotide = extract_segments(rotated, [{"start0": start0, "end0": end0}], strand)
            table = CodonTable.unambiguous_dna_by_id[11]
            near5 = start0 <= 2 if strand == 1 else end0 >= genome["length_bp"] - 2
            near3 = end0 >= genome["length_bp"] - 2 if strand == 1 else start0 <= 2
            if effective_phase is None:
                # PHANOTATE 1.6.7's observed GFF exporter omits every phase.
                # Resolve its complete ORFs from both nucleotide endpoints;
                # edge partials need one valid endpoint to anchor the frame.
                has_start = nucleotide[:3] in table.start_codons
                has_stop = len(nucleotide) % 3 == 0 and nucleotide[-3:] in table.stop_codons
                if has_start and has_stop:
                    effective_phase = 0
                    phase_resolution = "inferred_zero_from_complete_start_stop_and_codon_length"
                elif near5 and not has_start and has_stop:
                    effective_phase = 0
                    phase_resolution = "inferred_zero_from_terminal_stop_and_5p_input_edge"
                elif near3 and has_start and not has_stop:
                    effective_phase = 0
                    phase_resolution = "inferred_zero_from_start_codon_and_3p_input_edge"
            if effective_phase is not None:
                coding = nucleotide[effective_phase:]
                partial5 = near5 and coding[:3] not in table.start_codons
                partial3 = near3 and (len(coding) % 3 != 0 or coding[-3:] not in table.stop_codons)
        if partial5:
            p = parts[0]
            parts[0] = SimpleLocation(
                BeforePosition(p.start) if strand == 1 else p.start,
                AfterPosition(p.end) if strand == -1 else p.end,
                strand=strand,
            )
        if partial3:
            p = parts[-1]
            parts[-1] = SimpleLocation(
                BeforePosition(p.start) if strand == -1 else p.start,
                AfterPosition(p.end) if strand == 1 else p.end,
                strand=strand,
            )
        loc = parts[0] if len(parts) == 1 else CompoundLocation(parts, operator="join")
        attribute_map = {
            key: unquote(value)
            for item in attributes.split(";")
            if "=" in item
            for key, value in [item.split("=", 1)]
        }
        feature = SeqFeature(
            loc,
            type="CDS",
            id=attribute_map.get("ID", "<unknown id>"),
            qualifiers={
                "transl_table": ["11"],
                "codon_start": [
                    str(effective_phase + 1) if effective_phase is not None else "unknown"
                ],
                "note": [attributes, f"score={score}"],
            },
        )
        call = normalize_cds(
            genome, feature, "phanotate", source_run_id, len(calls), rotation_shift=rotation_shift
        )
        call["raw_gff"] = line
        call["gff_phase_raw"] = phase
        call["effective_phase"] = effective_phase
        call["phase_resolution"] = phase_resolution
        call["gff_strand_raw"] = strand_text
        call["gff_attributes_raw"] = attributes
        call["edge_partial_inference"] = "missing-endpoint-codon-and-within-2bp-of-input-edge"
        call["stop_coordinate_convention"] = "1-based-inclusive;terminal-stop-validated-per-call"
        calls.append(call)
    return calls
