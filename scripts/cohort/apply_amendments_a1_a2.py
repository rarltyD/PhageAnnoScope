"""Apply documented protocol amendments A1/A2 (2026-10-18) to a validated manifest.

Input: manifest from ``phageannoscope validate-cohort`` run without host evidence.
Never modifies the input; writes a new manifest, an append-only event log and a
summary.  See docs/protocol_amendments.md.

A1  Host domain = NCBI Taxonomy domain of the INPHARED2 ``metadata_host`` name
    (scripts/cohort/resolve_host_genera.py).  Bacteria -> HOST_DOMAIN_UNRESOLVED removed;
    Archaea/Eukaryota -> NON_BACTERIAL_HOST (excluded); unresolved/placeholder ->
    stays on hold.
A2  RECODING_REVIEW_REQUIRED is removed when every flagged CDS of the genome
    carries only /ribosomal_slippage (programmed frameshift; the stop-codon
    system is unchanged).  Those CDS remain individually excluded as
    complex_biology by the matcher.  Genomes keep ``a2_slippage_only=true`` for
    the sensitivity analysis that drops them.  /transl_except or /exception keep
    the genome on hold.
Exact-sequence aliases are recomputed among the newly eligible rows.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO

EXCLUDED = {
    "OUTSIDE_BALTIMORE_I_II", "BELOW_PRIMARY_QUALITY", "MIUVIG_GENOME_FRAGMENT",
    "N_CONTAINING_SEQUENCE_SENSITIVITY_ONLY", "NON_DNA_OR_AMBIGUOUS_MOL_TYPE",
    "SEGMENTED_GENOME_SEPARATE_LAYER", "HOST_CONTAMINATION_CONFIRMED",
    "CELLULAR_SOURCE_TAXONOMY", "NON_BACTERIAL_HOST",
}
QUARANTINED = {
    "NCBI_IDENTITY_AMBIGUOUS", "MULTIPLE_RAW_SNAPSHOTS_REQUIRE_RESOLUTION",
    "NCBI_VERSION_MISSING", "METADATA_NCBI_ACCESSION_MISMATCH", "EMPTY_SEQUENCE",
    "METADATA_NCBI_LENGTH_MISMATCH", "FROZEN_SEQUENCE_SHA_MISMATCH",
    "FROZEN_GENBANK_SHA_MISMATCH", "NON_IUPAC_SEQUENCE",
}
RECODING_KEYS = ("transl_except", "exception", "ribosomal_slippage")


def slippage_only(path):
    flagged = []
    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                keys = [k for k in RECODING_KEYS if k in feature.qualifiers]
                if keys:
                    flagged.append(keys)
    return bool(flagged) and all(keys == ["ribosomal_slippage"] for keys in flagged), flagged


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--hosts", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    hosts = {r["metadata_host"]: r for r in csv.DictReader(open(args.hosts, encoding="utf-8"))}
    with open(args.manifest, encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames)
        rows = list(reader)
    for extra in ("host_evidence_level", "host_scientific_name", "a1_host_status",
                  "a2_slippage_only", "amendments_applied"):
        if extra not in fields:
            fields.append(extra)
    now = datetime.now(timezone.utc).isoformat()
    events = []
    for row in rows:
        reasons = json.loads(row["all_reason_codes_json"] or "[]")
        before = list(reasons)
        applied = []
        if row.get("ncbi_fetch_status") == "success":
            host = hosts.get(row.get("metadata_host", ""))
            status = host["host_domain_status"] if host else "not_queried"
            row["a1_host_status"] = f"{status}:{host['host_domain'] if host else ''}"
            if host and status == "verified":
                row.update(host_domain=host["host_domain"], host_domain_status="verified",
                           host_taxid=host["host_taxid"],
                           host_scientific_name=host["host_scientific_name"],
                           host_evidence_level="A1_metadata_host_" + host["evidence_level"],
                           host_domain_evidence=host["fetch_evidence"]
                           or host["subtree_evidence_json"])
                if "HOST_DOMAIN_UNRESOLVED" in reasons:
                    reasons.remove("HOST_DOMAIN_UNRESOLVED")
                if host["host_domain"] != "Bacteria" and "NON_BACTERIAL_HOST" not in reasons:
                    reasons.append("NON_BACTERIAL_HOST")
                applied.append("A1")
            row["a2_slippage_only"] = "false"
            if "RECODING_REVIEW_REQUIRED" in reasons:
                only, flagged = slippage_only(row["genbank_path"])
                if only:
                    reasons.remove("RECODING_REVIEW_REQUIRED")
                    row["a2_slippage_only"] = "true"
                    applied.append("A2")
        excluded = [c for c in reasons if c in EXCLUDED]
        quarantined = [c for c in reasons if c in QUARANTINED]
        if reasons != before:
            row["all_reason_codes_json"] = json.dumps(reasons)
            row["primary_reason_code"] = next(
                iter(quarantined or excluded or reasons), "SEQUENCE_METADATA_ELIGIBLE")
            row["selection_status"] = ("quarantine" if quarantined else "excluded_primary"
                                       if excluded else "hold_validation" if reasons
                                       else "eligible_for_pilot")
            row["eligible_for_pilot"] = str(not reasons).lower()
            row["sequence_qc_status"] = row["selection_status"]
        row["amendments_applied"] = ",".join(applied)
        if applied:
            events.append({"accession_version": row.get("accession_version"),
                           "stage": "protocol_amendment_A1_A2", "decision":
                           row["selection_status"], "before_json": json.dumps(before),
                           "after_json": json.dumps(reasons), "applied": ",".join(applied),
                           "evaluated_at_utc": now})
    groups = defaultdict(list)
    for row in rows:
        if row.get("sequence_sha256"):
            groups[row["sequence_sha256"]].append(row)
    for group in groups.values():
        eligible = [r for r in group if r["eligible_for_pilot"] == "true"]
        if len(eligible) < 2:
            continue
        representative = min(r["accession_version"] for r in eligible)
        for row in eligible:
            if row["accession_version"] != representative:
                row["eligible_for_pilot"] = "false"
                row["selection_status"] = "exact_sequence_alias"
                reasons = json.loads(row["all_reason_codes_json"]) + [
                    "EXACT_SEQUENCE_DUPLICATE_ALIAS"]
                row["all_reason_codes_json"] = json.dumps(reasons)
                row["primary_reason_code"] = "EXACT_SEQUENCE_DUPLICATE_ALIAS"
                row["exact_sequence_representative"] = representative
                events.append({"accession_version": row["accession_version"],
                               "stage": "exact_alias_after_amendment",
                               "decision": "exact_sequence_alias", "before_json": "",
                               "after_json": json.dumps(reasons), "applied": representative,
                               "evaluated_at_utc": now})
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with out.with_name("amendment_events.csv").open("a", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=list(events[0]))
        if h.tell() == 0:
            writer.writeheader()
        writer.writerows(events)
    primary = [r for r in rows if r.get("cohort_route") == "primary_candidate"]
    summary = {
        "input": str(Path(args.manifest).resolve()), "hosts": str(Path(args.hosts).resolve()),
        "out": str(out.resolve()), "at_utc": now,
        "primary_candidates": len(primary),
        "primary_status": dict(Counter(r["selection_status"] for r in primary)),
        "eligible_for_pilot": sum(r["eligible_for_pilot"] == "true" for r in rows),
        "eligible_slippage_only": sum(r["eligible_for_pilot"] == "true"
                                      and r.get("a2_slippage_only") == "true" for r in rows),
        "remaining_reason_counts_primary": dict(Counter(
            c for r in primary for c in json.loads(r["all_reason_codes_json"]))),
    }
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
