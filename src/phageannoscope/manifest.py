"""Evidence-based cohort validation and reproducible, nested length sampling."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

IUPAC_DNA = set("ACGTRYSWKMBDHVN")
BINS = (0, 10000, 50000, 100000, 200000, math.inf)
BIN_NAMES = ("lt10k", "10k_50k", "50k_100k", "100k_200k", "ge200k")


def _read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _true(value):
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _int(value):
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _json(value):
    if isinstance(value, (dict, list)):
        return value
    try:
        return json.loads(value or "null")
    except (ValueError, TypeError):
        return None


def _events(path, rows):
    fields = [
        "run_id",
        "accession",
        "accession_version",
        "stage",
        "decision",
        "primary_reason_code",
        "all_reason_codes_json",
        "evidence_json",
        "evaluated_at_utc",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    empty = not path.exists() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        if empty:
            writer.writeheader()
        writer.writerows(rows)


def _raw_records(directory):
    from Bio import SeqIO

    directory = Path(directory)
    paths = sorted(
        {p for pattern in ("*.gb", "*.gbk", "*.gbf", "*.genbank") for p in directory.rglob(pattern)}
    )
    records = defaultdict(list)
    errors = []
    for path in paths:
        try:
            with path.open("rb") as handle:
                raw_sha = hashlib.file_digest(handle, "sha256").hexdigest()
            entries = []
            with path.open(encoding="utf-8") as handle:
                for ordinal, record in enumerate(SeqIO.parse(handle, "genbank")):
                    sequence = str(record.seq).upper()
                    entries.append(
                        (
                            record.id,
                            {
                                "record_ordinal": ordinal,
                                "path": str(path.resolve()),
                                "raw_sha": raw_sha,
                                "sequence_sha": _sha(sequence.encode("ascii")),
                            },
                        )
                    )
            if not entries:
                raise ValueError("no GenBank records")
            for identity, entry in entries:
                records[identity].append(entry)
        except Exception as exc:
            errors.append(
                {
                    "path": str(path.resolve()),
                    "reason": "GENBANK_PARSE_ERROR",
                    "error_type": type(exc).__name__,
                }
            )
    return records, errors


def validate_cohort(manifest, genbank_dir, out, host_evidence=None, data_root=None):
    """Write a new manifest. Unknown evidence remains a hold, never a guessed pass.

    This validates eligibility to attempt the pilot. Paired caller success still
    controls in_primary_analysis downstream, so that field is always false here.
    """
    from .genbank_io import genome_from_record
    from .hosts import _record_at

    source_rows = _read(manifest)
    verified_hosts = None
    if host_evidence:
        from .hosts import load_host_evidence

        verified_hosts = load_host_evidence(host_evidence, data_root=data_root)
    out = Path(out)
    if out.resolve() == Path(manifest).resolve():
        raise ValueError("Validation must not overwrite the frozen input manifest")
    by_version, parse_errors = _raw_records(genbank_dir)
    fetched = {}
    index = Path(genbank_dir) / "fetch_records.jsonl"
    if index.exists():
        for line in index.read_text(encoding="utf-8").splitlines():
            item = json.loads(line)
            if item.get("status") == "success":
                fetched[item["accession_version"]] = item
    by_base = defaultdict(list)
    for version in by_version:
        by_base[version.rsplit(".", 1)[0]].append(version)
    run_id = "cohort_" + uuid.uuid4().hex
    timestamp = datetime.now(timezone.utc).isoformat()
    rows, event_rows, quarantines, molecule_events = [], [], list(parse_errors), []
    seen_versions = set()
    for source in source_rows:
        row = dict(source)
        accession = row.get("accession", "")
        requested = row.get("accession_version") or accession
        evidence, reasons, exclude, quarantine = {}, [], False, False
        exclusion_codes, quarantine_codes = [], []

        def reason(code, value=None, *, excluded=False, quarantined=False):
            nonlocal exclude, quarantine
            if code not in reasons:
                reasons.append(code)
            if value is not None:
                evidence[code] = value
            exclude |= excluded
            quarantine |= quarantined
            if excluded and code not in exclusion_codes:
                exclusion_codes.append(code)
            if quarantined and code not in quarantine_codes:
                quarantine_codes.append(code)

        baltimore = row.get("baltimore_group", "")
        if baltimore not in {"Group I", "Group II"}:
            if baltimore and baltimore != "Unclassified":
                reason("OUTSIDE_BALTIMORE_I_II", baltimore, excluded=True)
            else:
                reason("MOLECULE_UNRESOLVED", baltimore)
        quality = row.get("checkv_quality", "")
        if quality not in {"Complete", "High-quality"}:
            reason(
                "BELOW_PRIMARY_QUALITY"
                if quality in {"Medium-quality", "Low-quality"}
                else "QUALITY_UNRESOLVED",
                quality,
                excluded=quality in {"Medium-quality", "Low-quality"},
            )
        if "genome-fragment" in row.get("checkv_miuvig_quality", "").lower():
            reason("MIUVIG_GENOME_FRAGMENT", row["checkv_miuvig_quality"], excluded=True)
        metadata_length = _int(row.get("metadata_length_bp"))
        if not metadata_length or metadata_length < 0:
            reason("INVALID_METADATA_LENGTH", row.get("metadata_length_bp"))
        versions = (
            [requested]
            if "." in requested and requested in by_version
            else (by_base.get(requested, []) if "." not in requested else [])
        )
        if len(versions) != 1:
            reason(
                "NCBI_RECORD_MISSING" if not versions else "NCBI_IDENTITY_AMBIGUOUS",
                versions,
                quarantined=bool(versions),
            )
            row["ncbi_fetch_status"] = "missing" if not versions else "ambiguous"
        else:
            version = versions[0]
            if version in seen_versions:
                quarantines.append(
                    {
                        "accession": accession,
                        "accession_version": version,
                        "reason": "DUPLICATE_MANIFEST_IDENTITY",
                    }
                )
                continue
            seen_versions.add(version)
            candidates = by_version[version]
            if len({(x["sequence_sha"], x["raw_sha"]) for x in candidates}) > 1:
                reason(
                    "MULTIPLE_RAW_SNAPSHOTS_REQUIRE_RESOLUTION",
                    [x["path"] for x in candidates],
                    quarantined=True,
                )
            item = sorted(candidates, key=lambda x: x["path"])[0]
            # The index holds hashes/locations only. Reuse a bounded two-batch
            # cache instead of reparsing a raw batch for each selected record.
            record = _record_at(item["path"], item["record_ordinal"])
            sequence = str(record.seq).upper()
            info = fetched.get(version, {})
            if info.get("genbank_sha256") != item["raw_sha"]:
                info = {}
            row.update(
                accession_version=version,
                genome_id=genome_from_record(record)["genome_id"],
                ncbi_fetch_status="success",
                genbank_path=item["path"],
                genbank_sha256=item["raw_sha"],
                sequence_sha256=item["sequence_sha"],
                ncbi_length_bp=len(sequence),
                topology=record.annotations.get("topology", "unknown"),
            )
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.[1-9][0-9]*", version):
                reason("NCBI_VERSION_MISSING", version, quarantined=True)
            if accession and accession != version.rsplit(".", 1)[0]:
                reason(
                    "METADATA_NCBI_ACCESSION_MISMATCH",
                    {"metadata": accession, "ncbi": version},
                    quarantined=True,
                )
            if not sequence:
                reason("EMPTY_SEQUENCE", quarantined=True)
            if metadata_length and metadata_length != len(sequence):
                reason(
                    "METADATA_NCBI_LENGTH_MISMATCH",
                    {"metadata": metadata_length, "ncbi": len(sequence)},
                    quarantined=True,
                )
            if source.get("sequence_sha256") and source["sequence_sha256"] != item["sequence_sha"]:
                reason(
                    "FROZEN_SEQUENCE_SHA_MISMATCH",
                    {"expected": source["sequence_sha256"], "actual": item["sequence_sha"]},
                    quarantined=True,
                )
            if source.get("genbank_sha256") and source["genbank_sha256"] != item["raw_sha"]:
                reason(
                    "FROZEN_GENBANK_SHA_MISMATCH",
                    {"expected": source["genbank_sha256"], "actual": item["raw_sha"]},
                    quarantined=True,
                )
            invalid = sorted(set(sequence) - IUPAC_DNA)
            if invalid:
                reason("NON_IUPAC_SEQUENCE", invalid, quarantined=True)
            if "N" in sequence:
                reason(
                    "N_CONTAINING_SEQUENCE_SENSITIVITY_ONLY",
                    {"n_count": sequence.count("N")},
                    excluded=True,
                )
            ambiguous = sorted((set(sequence) & IUPAC_DNA) - set("ACGTN"))
            if ambiguous:
                # Amendment A5: non-N ambiguity codes follow the N rule (sensitivity only).
                reason(
                    "AMBIGUOUS_IUPAC_SENSITIVITY_ONLY",
                    {c: sequence.count(c) for c in ambiguous},
                    excluded=True,
                )
            row["n_count"] = sequence.count("N")
            row["non_iupac_characters_json"] = json.dumps(invalid)
            source_qualifiers = [dict(f.qualifiers) for f in record.features if f.type == "source"]
            row["source_qualifiers_json"] = json.dumps(source_qualifiers, ensure_ascii=False)
            mol_types = sorted({m for q in source_qualifiers for m in q.get("mol_type", [])})
            row["mol_types_json"] = json.dumps(mol_types)
            if not mol_types:
                reason("MOL_TYPE_UNRESOLVED")
            elif any("dna" not in x.lower() for x in mol_types):
                reason("NON_DNA_OR_AMBIGUOUS_MOL_TYPE", mol_types, excluded=True)
            segments = [q.get("segment") for q in source_qualifiers if q.get("segment")]
            if segments or _true(row.get("segmented_genome")):
                reason("SEGMENTED_GENOME_SEPARATE_LAYER", segments, excluded=True)
            if any("proviral" in q for q in source_qualifiers):
                reason("PROVIRAL_FRAGMENT_REVIEW_REQUIRED", source_qualifiers)
            if _true(row.get("host_contamination_confirmed")):
                reason(
                    "HOST_CONTAMINATION_CONFIRMED",
                    row.get("host_contamination_evidence"),
                    excluded=True,
                )
            viral_taxids = sorted(
                {
                    v.split(":", 1)[1]
                    for q in source_qualifiers
                    for v in q.get("db_xref", [])
                    if v.startswith("taxon:")
                }
            )
            row["ncbi_taxid"] = viral_taxids[0] if len(viral_taxids) == 1 else ""
            if len(viral_taxids) != 1:
                reason("VIRUS_TAXID_UNRESOLVED", viral_taxids)
            taxonomy = record.annotations.get("taxonomy", [])
            row["ncbi_taxonomy_json"] = json.dumps(taxonomy)
            if "Viruses" not in taxonomy:
                cellular = any(
                    domain in taxonomy for domain in ("Bacteria", "Archaea", "Eukaryota")
                )
                reason(
                    "CELLULAR_SOURCE_TAXONOMY" if cellular else "VIRUS_TAXONOMY_UNRESOLVED",
                    taxonomy,
                    excluded=cellular,
                )
            tables, recoding = set(), []
            for feature in record.features:
                if feature.type == "CDS":
                    raw = feature.qualifiers.get("transl_table", ["1"])
                    tables.update(str(x) for x in raw)
                    if any(
                        x in feature.qualifiers
                        for x in ("transl_except", "exception", "ribosomal_slippage")
                    ):
                        recoding.append(str(feature.location))
            row["translation_tables_json"] = json.dumps(sorted(tables))
            if tables - {"1", "11"}:
                reason("ALTERNATIVE_TRANSLATION_CODE_REVIEW_REQUIRED", sorted(tables))
            if recoding:
                reason("RECODING_REVIEW_REQUIRED", recoding)
            # Prefer verified lookup provenance from the immutable fetch index.
            for key in (
                "host_domain",
                "host_domain_status",
                "host_taxid",
                "host_domain_evidence",
                "host_scientific_name",
            ):
                if info.get(key):
                    row[key] = info[key]
            if info.get("host_taxonomy_lineage"):
                row["host_taxonomy_lineage_json"] = json.dumps(info["host_taxonomy_lineage"])
            if verified_hosts is not None:
                host = verified_hosts.get(version)
                # An explicitly supplied evidence table is authoritative for this
                # validation pass; stale evidence from a fetch must not bypass it.
                for key in (
                    "host_domain",
                    "host_domain_status",
                    "host_taxid",
                    "host_domain_evidence",
                    "host_scientific_name",
                    "host_rank",
                    "host_evidence_level",
                ):
                    row[key] = (host or {}).get(key) or ""
                if host:
                    if host.get("sequence_sha256") != item["sequence_sha"]:
                        reason("HOST_EVIDENCE_SEQUENCE_MISMATCH", host, quarantined=True)
                        row["host_domain_status"] = "quarantine"
                    if host.get("genbank_sha256") != item["raw_sha"]:
                        reason(
                            "HOST_EVIDENCE_GENBANK_MISMATCH",
                            host.get("genbank_sha256"),
                            quarantined=True,
                        )
                        row["host_domain_status"] = "quarantine"
                    from .hosts import host_inputs_sha256

                    if host.get("host_inputs_sha256") != host_inputs_sha256(
                        source, source_qualifiers
                    ):
                        reason(
                            "HOST_EVIDENCE_INPUTS_CHANGED",
                            "Run resolve-hosts for the current host metadata",
                        )
                        row["host_domain_status"] = "stale"
                    row["host_taxonomy_rank"] = host.get("host_taxonomy_rank", "")
                    row["host_rank"] = host.get("host_taxonomy_rank", "")
                    row["host_evidence_data_root"] = host.get("host_evidence_data_root", "")
                    if host.get("host_domain_status") == "conflict":
                        reason("HOST_DOMAIN_CONFLICT", host.get("host_domains_verified"))
                    row["host_taxonomy_evidence_json"] = host.get(
                        "host_taxonomy_evidence_json", "[]"
                    )
                    row["host_evidence_table"] = str(Path(host_evidence).resolve())
                else:
                    row["host_domain_status"] = "not_resolved"
            host_domain = row.get("host_domain", "")
            host_status = row.get("host_domain_status", "")
            host_domain_evidence = row.get("host_domain_evidence", "")
            if (
                not host_domain
                or host_status not in {"verified", "resolved", "confirmed"}
                or not host_domain_evidence
            ):
                reason(
                    "HOST_DOMAIN_UNRESOLVED",
                    {
                        "domain": host_domain,
                        "status": host_status,
                        "evidence": host_domain_evidence,
                    },
                )
            elif host_domain != "Bacteria":
                reason(
                    "NON_BACTERIAL_HOST",
                    {"domain": host_domain, "evidence": host_domain_evidence},
                    excluded=True,
                )
        row["checkv_version"] = row.get("checkv_version") or "not_reported"
        row["checkv_database_version"] = row.get("checkv_database_version") or "not_reported"
        row["all_reason_codes_json"] = json.dumps(reasons)
        row["primary_reason_code"] = next(
            iter(quarantine_codes or exclusion_codes or reasons), "SEQUENCE_METADATA_ELIGIBLE"
        )
        row["validation_evidence_json"] = json.dumps(evidence, ensure_ascii=False)
        row["sequence_qc_status"] = (
            "quarantine"
            if quarantine
            else "excluded_primary"
            if exclude
            else "hold"
            if reasons
            else "passed"
        )
        row["selection_status"] = (
            "quarantine"
            if quarantine
            else "excluded_primary"
            if exclude
            else "hold_validation"
            if reasons
            else "eligible_for_pilot"
        )
        row["eligible_for_pilot"] = str(not reasons).lower()
        row["in_primary_analysis"] = "false"
        row["validation_run_id"] = run_id
        row["validated_at_utc"] = timestamp
        if baltimore not in {"Group I", "Group II"} and (
            not baltimore or baltimore == "Unclassified"
        ):
            # DNA mol_type alone does not establish Baltimore I versus II.
            # Keep this separate from other quality/exclusion decisions.
            molecule_evidence = {
                "baltimore_metadata": baltimore,
                "ncbi_fetch_status": row.get("ncbi_fetch_status"),
                "mol_types": _json(row.get("mol_types_json")),
                "ncbi_taxonomy": _json(row.get("ncbi_taxonomy_json")),
                "rule": "unresolved Baltimore group remains hold until evidence-specific review; DNA text alone insufficient",
            }
            row["molecule_resolution_status"] = (
                "review_required" if row.get("ncbi_fetch_status") == "success" else "not_downloaded"
            )
            molecule_events.append(
                {
                    "run_id": run_id,
                    "accession": accession,
                    "accession_version": row.get("accession_version"),
                    "stage": "molecule_resolution",
                    "decision": row["molecule_resolution_status"],
                    "primary_reason_code": "MOLECULE_UNRESOLVED",
                    "all_reason_codes_json": '["MOLECULE_UNRESOLVED"]',
                    "evidence_json": json.dumps(molecule_evidence),
                    "evaluated_at_utc": timestamp,
                }
            )
        rows.append(row)
        if quarantine:
            quarantines.append(
                {
                    "accession": accession,
                    "accession_version": row.get("accession_version"),
                    "reason": row["all_reason_codes_json"],
                    "evidence": row["validation_evidence_json"],
                }
            )
    duplicate_groups = defaultdict(list)
    for row in rows:
        if row.get("sequence_sha256"):
            duplicate_groups[row["sequence_sha256"]].append(row)
    for sha, group in duplicate_groups.items():
        if len(group) < 2:
            continue
        eligible = [r for r in group if _true(r["eligible_for_pilot"])]
        representative = min(eligible or group, key=lambda r: r["accession_version"])[
            "accession_version"
        ]
        for row in group:
            row["exact_sequence_duplicate_group"] = "exact_" + sha
            row["exact_sequence_representative"] = representative
            row["exact_sequence_aliases_json"] = json.dumps(
                sorted(r["accession_version"] for r in group)
            )
            if row["accession_version"] != representative and _true(row["eligible_for_pilot"]):
                row["eligible_for_pilot"] = "false"
                row["selection_status"] = "exact_sequence_alias"
                reasons = json.loads(row["all_reason_codes_json"])
                reasons.append("EXACT_SEQUENCE_DUPLICATE_ALIAS")
                row["all_reason_codes_json"] = json.dumps(reasons)
                row["primary_reason_code"] = "EXACT_SEQUENCE_DUPLICATE_ALIAS"
    for row in rows:
        event_rows.append(
            {
                "run_id": run_id,
                "accession": row.get("accession"),
                "accession_version": row.get("accession_version"),
                "stage": "ncbi_sequence_validation",
                "decision": row["selection_status"],
                "primary_reason_code": row["primary_reason_code"],
                "all_reason_codes_json": row["all_reason_codes_json"],
                "evidence_json": row["validation_evidence_json"],
                "evaluated_at_utc": timestamp,
            }
        )
    _write(
        out,
        rows,
        list(
            dict.fromkeys([*(source_rows[0] if source_rows else {}), *(k for r in rows for k in r)])
        ),
    )
    events_path = out.parent / "selection_events.csv"
    _events(events_path, event_rows)
    molecule_path = out.parent / "molecule_resolution_events.csv"
    _events(molecule_path, molecule_events)
    quarantine_path = out.with_name(out.stem + "_quarantine.csv")
    _write(
        quarantine_path,
        quarantines,
        ["accession", "accession_version", "path", "reason", "evidence", "error_type"],
    )
    summary = {
        "run_id": run_id,
        "input_rows": len(source_rows),
        "output_rows": len(rows),
        "eligible_for_pilot": sum(_true(r["eligible_for_pilot"]) for r in rows),
        "in_primary_analysis": 0,
        "counts": dict(Counter(r["selection_status"] for r in rows)),
        "out": str(out.resolve()),
        "selection_events": str(events_path.resolve()),
        "quarantine": str(quarantine_path.resolve()),
        "quarantine_rows": len(quarantines),
        "molecule_resolution_events": str(molecule_path.resolve()),
        "host_evidence": str(Path(host_evidence).resolve()) if host_evidence else None,
        "exact_duplicate_detection": "exact orientation SHA256 only; rotations and reverse complements not inferred",
    }
    out.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def _bin(row):
    length = _int(row.get("ncbi_length_bp") or row.get("length_bp"))
    if not length or length <= 0:
        raise ValueError("Eligible row lacks positive validated sequence length")
    return next(i for i in range(5) if BINS[i] <= length < BINS[i + 1])


def _rank(seed, row):
    return _sha(f"{seed},{row['accession_version']}".encode("utf-8"))


def _block(row):
    for field in ("relatedness_cluster", "taxmyphage_species", "taxmyphage_genus"):
        value = (row.get(field) or "").strip()
        if value and value.casefold() not in {
            "new_species",
            "new_genus",
            "unclassified",
            "unknown",
            "na",
            "n/a",
        }:
            return field + ":" + value
    return "unresolved:" + row["accession_version"]


def select_manifest(
    manifest, n, seed=30102026, strategy="protocol-length-strata", out=None, include=None
):
    if strategy != "protocol-length-strata":
        raise ValueError("Only protocol-length-strata is implemented")
    if not out or n <= 0:
        raise ValueError("A positive n and output path are required")
    all_rows = _read(manifest)
    eligible = [
        r
        for r in all_rows
        if _true(r.get("eligible_for_pilot"))
        or _true(r.get("in_primary_analysis"))
        or r.get("selection_status") in {"eligible_primary", "final_eligible"}
    ]
    ids = [r.get("accession_version") for r in eligible]
    if any(not x or "." not in x for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("Eligible rows need unique explicit accession.version identities")
    shas = [r.get("sequence_sha256") for r in eligible]
    if any(not x for x in shas) or len(set(shas)) != len(shas):
        raise ValueError(
            "Eligible rows need unique sequence SHA256; run cohort deduplication first"
        )
    if len(eligible) < n:
        raise ValueError(
            f"Insufficient final eligible genomes: requested {n}, available {len(eligible)}; pilot not selected"
        )
    by_id = {r["accession_version"]: r for r in eligible}
    included_rows = _read(include) if include else []
    included = set()
    for previous in included_rows:
        identity = previous.get("accession_version")
        if identity in included or identity not in by_id:
            raise ValueError("Nested include contains a duplicate or ineligible genome")
        if previous.get("sequence_sha256") != by_id[identity].get("sequence_sha256"):
            raise ValueError("Nested include sequence SHA mismatch")
        included.add(identity)
    if len(included) > n:
        raise ValueError("Nested include is larger than requested sample")
    strata = [[r for r in eligible if _bin(r) == i] for i in range(5)]
    allocation = [
        max(min(5, len(group)), sum(r["accession_version"] in included for r in group))
        for group in strata
    ]
    if sum(allocation) > n:
        raise ValueError(
            "Requested n cannot satisfy minimum per nonempty length stratum and nested include"
        )
    # Allocate remaining seats with largest remainders, respecting finite strata.
    while sum(allocation) < n:
        remaining = n - sum(allocation)
        active = [i for i in range(5) if allocation[i] < len(strata[i])]
        total_weight = sum(len(strata[i]) for i in active)
        quotas = {i: remaining * len(strata[i]) / total_weight for i in active}
        floors = {i: min(len(strata[i]) - allocation[i], math.floor(quotas[i])) for i in active}
        for i, count in floors.items():
            allocation[i] += count
        seats = n - sum(allocation)
        if seats:
            order = sorted(active, key=lambda i: (-(quotas[i] - math.floor(quotas[i])), i))
            for i in order:
                if seats and allocation[i] < len(strata[i]):
                    allocation[i] += 1
                    seats -= 1
    selected = []
    for i, group in enumerate(strata):
        chosen = [r for r in group if r["accession_version"] in included]
        block_counts = Counter(_block(r) for r in chosen)
        remaining = sorted(
            [r for r in group if r["accession_version"] not in included],
            key=lambda r: (_rank(seed, r), r["accession_version"]),
        )
        while len(chosen) < allocation[i]:
            lowest = min(block_counts[_block(r)] for r in remaining)
            candidate = next(r for r in remaining if block_counts[_block(r)] == lowest)
            remaining.remove(candidate)
            chosen.append(candidate)
            block_counts[_block(candidate)] += 1
        for row in sorted(chosen, key=lambda r: (_rank(seed, r), r["accession_version"])):
            selected.append(
                {
                    **row,
                    "selection_seed": seed,
                    "selection_strategy": strategy,
                    "selection_stratum": BIN_NAMES[i],
                    "selection_hash": _rank(seed, row),
                    "nested_included": str(row["accession_version"] in included).lower(),
                    "selection_order": len(selected) + 1,
                }
            )
    _write(out, selected)
    summary = {
        "requested_n": n,
        "selected_n": len(selected),
        "eligible_n": len(eligible),
        "seed": seed,
        "strategy": strategy,
        "nested_included_n": len(included),
        "strata": dict(zip(BIN_NAMES, allocation)),
        "hash_serialization": "UTF-8 decimal seed + comma + accession.version",
        "underfilled_minimum_strata": [
            BIN_NAMES[i] for i, g in enumerate(strata) if 0 < len(g) < 5
        ],
        "out": str(Path(out).resolve()),
    }
    Path(out).with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    timestamp = datetime.now(timezone.utc).isoformat()
    _events(
        Path(out).parent / "selection_events.csv",
        [
            {
                "run_id": "select_" + str(seed) + "_" + str(n),
                "accession": r.get("accession"),
                "accession_version": r["accession_version"],
                "stage": "protocol_sampling",
                "decision": "selected",
                "primary_reason_code": "PROTOCOL_LENGTH_STRATA",
                "all_reason_codes_json": '["PROTOCOL_LENGTH_STRATA"]',
                "evidence_json": json.dumps(
                    {"stratum": r["selection_stratum"], "hash": r["selection_hash"], "seed": seed}
                ),
                "evaluated_at_utc": timestamp,
            }
            for r in selected
        ],
    )
    return summary
