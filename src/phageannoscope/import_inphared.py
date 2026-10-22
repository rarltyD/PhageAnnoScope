"""Inspect a local optional INPHARED archive; import only sequence-proven matches.

Archives are never extracted. All names and member types are checked before any
annotation is accepted. The historical archive never defines cohort membership.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


def _safe_name(name):
    name = name.replace("\\", "/")
    path = PurePosixPath(name)
    if (
        not name
        or path.is_absolute()
        or ".." in path.parts
        or re.match(r"^[A-Za-z]:", name)
        or "\x00" in name
        or ":" in name
    ):
        raise ValueError(f"Unsafe archive member path: {name!r}")
    return str(path)


def _members(archive):
    """Yield (metadata, read_bytes); callers keep the generator open when reading."""
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as handle:
            for member in handle.infolist():
                mode = member.external_attr >> 16
                kind = "directory" if member.is_dir() else "file"
                if stat.S_ISLNK(mode) or (
                    stat.S_IFMT(mode) and not (stat.S_ISREG(mode) or stat.S_ISDIR(mode))
                ):
                    raise ValueError("Archive links and special files are forbidden")
                info = {
                    "name": _safe_name(member.filename),
                    "size_bytes": member.file_size,
                    "kind": kind,
                }
                yield info, lambda m=member: handle.read(m)
    else:
        with tarfile.open(archive, "r:*") as handle:
            for member in handle:
                if not (member.isfile() or member.isdir()):
                    raise ValueError("Archive links and special files are forbidden")
                info = {
                    "name": _safe_name(member.name),
                    "size_bytes": member.size,
                    "kind": "file" if member.isfile() else "directory",
                }

                def read(m=member):
                    stream = handle.extractfile(m)
                    if stream is None:
                        return b""
                    with stream:
                        return stream.read()

                yield info, read


def inventory_archive(path, out=None, *, max_members=100000, max_bytes=2_000_000_000):
    """Validate a local archive and list metadata without extracting any members."""
    path = Path(path)
    members, seen, total = [], set(), 0
    for info, _ in _members(path):
        if info["name"] in seen:
            raise ValueError("Duplicate normalized archive member name")
        seen.add(info["name"])
        total += info["size_bytes"]
        if len(members) >= max_members or total > max_bytes:
            raise ValueError("Archive exceeds configured inventory member/size budget")
        members.append(info)
    with path.open("rb") as handle:
        archive_sha = hashlib.file_digest(handle, "sha256").hexdigest()
    result = {
        "archive": str(path.resolve()),
        "archive_sha256": archive_sha,
        "archive_size_bytes": path.stat().st_size,
        "members": members,
        "member_count": len(members),
        "uncompressed_bytes": total,
        "status": "validated_inventory",
    }
    if out:
        dest = Path(out)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def _reference_rows(manifest):
    if isinstance(manifest, (str, Path)):
        with Path(manifest).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    return list(manifest)


def _record_identity(record):
    identity = str(record.id)
    if "|" in identity:
        fields = identity.split("|")
        candidates = [
            fields[i + 1]
            for i, field in enumerate(fields[:-1])
            if field in {"gb", "ref", "emb", "dbj"}
        ]
        if len(candidates) == 1:
            identity = candidates[0]
    return identity


def import_inphared(archive, manifest, out, *, max_members=100000, max_bytes=2_000_000_000):
    """Import optional GenBank annotations only after accession and exact SHA agree.

    Base accessions are useful to report potential aliases, but an unversioned
    historical record cannot silently acquire the current NCBI version. FASTA
    members are checked for sequence identity and retained as evidence only.
    """
    from Bio import SeqIO

    archive, out = Path(archive), Path(out)
    inventory = inventory_archive(archive, max_members=max_members, max_bytes=max_bytes)
    refs = _reference_rows(manifest)
    by_version = {}
    for ref in refs:
        version, sha = ref.get("accession_version"), ref.get("sequence_sha256")
        if not version or not sha:
            continue
        if version in by_version and by_version[version]["sequence_sha256"] != sha:
            raise ValueError("Reference manifest has conflicting sequence hashes")
        by_version[version] = ref
    matches, quarantine, imported = (
        [],
        [],
        {"genomes": [], "features": [], "cds_calls": [], "qc_events": []},
    )
    # Preflight completes before creating output or accepting a single member.
    out.mkdir(parents=True, exist_ok=True)
    (out / "inventory.json").write_text(
        json.dumps(inventory, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for info, read in _members(archive):
        suffix = PurePosixPath(info["name"]).suffix.lower()
        if info["kind"] != "file" or suffix not in {
            ".gb",
            ".gbk",
            ".gbf",
            ".genbank",
            ".fna",
            ".fasta",
            ".fsa",
        }:
            continue
        raw = read()
        raw_sha = hashlib.sha256(raw).hexdigest()
        fmt = "genbank" if suffix in {".gb", ".gbk", ".gbf", ".genbank"} else "fasta"
        try:
            records = list(SeqIO.parse(io.StringIO(raw.decode("utf-8")), fmt))
            if not records:
                raise ValueError("empty sequence member")
        except Exception as exc:
            quarantine.append(
                {
                    "member": info["name"],
                    "reason": "ANNOTATION_PARSE_FAILED",
                    "error_type": type(exc).__name__,
                }
            )
            continue
        identities = [_record_identity(record) for record in records]
        if len(set(identities)) != len(identities):
            quarantine.append(
                {
                    "member": info["name"],
                    "reason": "DUPLICATE_RECORD_IDENTITY",
                    "accession_versions": identities,
                }
            )
            continue
        accepted = []
        for record in records:
            # Documented NCBI pipe forms; no substring accession matching.
            identity = _record_identity(record)
            sequence_sha = hashlib.sha256(str(record.seq).upper().encode("ascii")).hexdigest()
            reference = by_version.get(identity)
            if reference is None:
                base = identity.rsplit(".", 1)[0]
                possible = sorted(v for v in by_version if v.rsplit(".", 1)[0] == base)
                quarantine.append(
                    {
                        "member": info["name"],
                        "accession_version": identity,
                        "sequence_sha256": sequence_sha,
                        "reason": "ACCESSION_VERSION_NOT_MATCHED",
                        "possible_versions": possible,
                    }
                )
                continue
            if sequence_sha != reference["sequence_sha256"]:
                quarantine.append(
                    {
                        "member": info["name"],
                        "accession_version": identity,
                        "sequence_sha256": sequence_sha,
                        "expected_sequence_sha256": reference["sequence_sha256"],
                        "reason": "SEQUENCE_SHA_MISMATCH",
                    }
                )
                continue
            matches.append(
                {
                    "member": info["name"],
                    "member_sha256": raw_sha,
                    "accession_version": identity,
                    "sequence_sha256": sequence_sha,
                    "source": "inphared_historical",
                    "format": fmt,
                    "annotation_imported": fmt == "genbank",
                }
            )
            accepted.append(identity)
        if fmt == "genbank" and accepted:
            # Preserve the exact member bytes under a content-addressed name.
            raw_path = out / "raw" / (raw_sha + ".gb")
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            if raw_path.exists() and raw_path.read_bytes() != raw:
                raise ValueError("Immutable annotation member hash collision")
            if not raw_path.exists():
                raw_path.write_bytes(raw)
            from .genbank_io import read_genbank

            parsed = read_genbank(raw_path, source="inphared", source_run_id="inphared_" + raw_sha)
            # Multi-record archives may contain unmatched records. Filter every
            # table by canonical genome identity; never import those by accident.
            genome_ids = {
                r["genome_id"]
                for r in parsed.get("genomes", [])
                if r.get("accession_version") in accepted
            }
            for table in imported:
                imported[table].extend(
                    r for r in parsed.get(table, []) if r.get("genome_id") in genome_ids
                )
    for name, value in (
        ("matches", matches),
        ("quarantine", quarantine),
        ("annotations", imported),
    ):
        (out / (name + ".json")).write_text(
            json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8"
        )
    summary = {
        "inventory_members": inventory["member_count"],
        "matched_records": len(matches),
        "quarantined_records": len(quarantine),
        "imported_cds_calls": len(imported["cds_calls"]),
        "out": str(out.resolve()),
        "annotations": imported,
        "matches": matches,
        "quarantine": quarantine,
        "limitations": "GFF/protein-only annotations require sequence-proven nucleotide linkage and are inventoried but not imported",
    }
    return summary
