"""FASTA identifiers are resolved explicitly, never by substring search."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from Bio import SeqIO

from .genbank_io import stable_id

ACCESSION_VERSION = re.compile(r"^[A-Za-z][A-Za-z0-9_]*\.\d+$")
PIPE_DATABASES = {"ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd"}


def resolve_fasta_identifier(header, aliases=None):
    aliases = aliases or {}
    identifier = header.split()[0]
    explicit = aliases.get(header, aliases.get(identifier))
    if explicit is not None:
        resolved = str(explicit)
    elif "|" not in identifier:
        resolved = identifier
    else:
        pieces = identifier.split("|")
        accessions = [
            pieces[i + 1] for i, field in enumerate(pieces[:-1]) if field.lower() in PIPE_DATABASES
        ]
        if len(set(accessions)) != 1:
            raise ValueError(
                f"Unrecognized or ambiguous pipe FASTA identifier: {identifier}; provide an explicit alias"
            )
        resolved = accessions[0]
    if not ACCESSION_VERSION.fullmatch(resolved):
        raise ValueError(
            f"FASTA identifier must resolve to accession.version: {resolved}; provide an explicit alias"
        )
    return resolved


def read_fasta(path, aliases=None):
    if aliases is not None and not isinstance(aliases, dict):
        with Path(aliases).open(encoding="utf-8") as handle:
            aliases = json.load(handle)
    genomes = []
    seen_headers = set()
    seen_accessions = set()
    for record in SeqIO.parse(str(Path(path)), "fasta"):
        header = record.description
        accession_version = resolve_fasta_identifier(header, aliases)
        if record.id in seen_headers:
            raise ValueError(f"Duplicate FASTA header: {record.id}")
        if accession_version in seen_accessions:
            raise ValueError(f"Duplicate FASTA accession.version: {accession_version}")
        seen_headers.add(record.id)
        seen_accessions.add(accession_version)
        sequence = str(record.seq).upper()
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        genomes.append(
            {
                "genome_id": stable_id(
                    "genome", {"accession_version": accession_version, "sequence_sha256": digest}
                ),
                "accession": re.sub(r"\.\d+$", "", accession_version),
                "accession_version": accession_version,
                "sequence_sha256": digest,
                "length_bp": len(sequence),
                "sequence": sequence,
                "topology": "unknown",
                "fasta_header": header,
            }
        )
    if not genomes:
        raise ValueError(f"No FASTA records found in {path}")
    return genomes


def validate_fasta_genbank(fasta_genomes, genbank_genomes):
    """Require the same exact versions and sequences in the two collections."""

    def index(genomes, label):
        result = {}
        for genome in genomes:
            key = genome["accession_version"]
            if key in result:
                raise ValueError(f"Duplicate {label} accession.version: {key}")
            result[key] = genome
        return result

    fasta, genbank = index(fasta_genomes, "FASTA"), index(genbank_genomes, "GenBank")
    for key, genome in fasta.items():
        if key not in genbank:
            base = genome["accession"]
            related = [k for k, g in genbank.items() if g["accession"] == base]
            if related:
                raise ValueError(f"FASTA/GenBank version mismatch: {key} versus {related}")
            raise ValueError(f"FASTA identifier absent from GenBank: {key}")
        if genome["sequence_sha256"] != genbank[key]["sequence_sha256"]:
            raise ValueError(f"FASTA/GenBank sequence SHA mismatch: {key}")
    if set(genbank) != set(fasta):
        raise ValueError(f"GenBank records missing from FASTA: {sorted(set(genbank) - set(fasta))}")
    return True
