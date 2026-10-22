#!/usr/bin/env python3
"""Freeze an INPHARED2 TSV as an auditable, preliminary cohort manifest.

This is metadata screening, not sequence validation or a gene-calling pipeline.
Only Python's standard library is required.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

MAPPING = {
    "Accession": "accession",
    "Description": "description",
    "Genome Length (bp)": "metadata_length_bp",
    "molGC (%)": "metadata_gc_percent",
    "Host": "metadata_host",
    "Baltimore Group": "baltimore_group",
    "Lowest Taxa": "ncbi_lowest_taxa_metadata",
    "Genus": "ncbi_genus_metadata",
    "Sub-family": "ncbi_subfamily_metadata",
    "Family": "ncbi_family_metadata",
    "Order": "ncbi_order_metadata",
    "Class": "ncbi_class_metadata",
    "Phylum": "ncbi_phylum_metadata",
    "Kingdom": "ncbi_kingdom_metadata",
    "Realm": "ncbi_realm_metadata",
    "Isolation Host": "isolation_host_metadata",
    "CheckV quality": "checkv_quality",
    "CheckV MIUVIG quality": "checkv_miuvig_quality",
    "CheckV completeness": "checkv_completeness_percent",
    "CheckV completeness method": "checkv_completeness_method",
    "TaxMyPhage genus": "taxmyphage_genus",
    "TaxMyPhage species": "taxmyphage_species",
    "TaxMyPhage message": "taxmyphage_message",
    "PhageLeads lifestyle": "phageleads_lifestyle",
    "PhageLeads problematic genes": "phageleads_problematic_genes",
}
EXTRA = [
    "accession_version", "source_snapshot_id", "source_retrieved_at_utc",
    "source_url", "metadata_sha256", "source_row_number",
    "screen_decision", "cohort_route", "primary_reason_code",
    "all_reason_codes_json", "selection_status", "in_primary_analysis",
    "ncbi_fetch_status", "genbank_path", "genbank_sha256", "sequence_sha256",
    "ncbi_length_bp", "ncbi_taxid", "host_taxid", "host_domain",
    "host_domain_status", "topology", "translation_tables_json",
    "sequence_qc_status", "exact_sequence_duplicate_group", "relatedness_cluster",
    "predictor_status", "notes",
]


def screen(row):
    reasons = []
    group = row["Baltimore Group"].strip()
    quality = row["CheckV quality"].strip()
    miuvig = row["CheckV MIUVIG quality"].strip()
    if group in {"Group III", "Group IV", "Group V", "Group VI", "Group VII"}:
        reasons.append("OUTSIDE_BALTIMORE_I_II")
    elif group not in {"Group I", "Group II"}:
        reasons.append("MOLECULE_UNRESOLVED")
    if quality in {"Medium-quality", "Low-quality"}:
        reasons.append("BELOW_PRIMARY_QUALITY")
    elif quality not in {"Complete", "High-quality"}:
        reasons.append("QUALITY_UNRESOLVED")
    if miuvig == "Genome-fragment":
        reasons.append("MIUVIG_GENOME_FRAGMENT")
    try:
        if int(row["Genome Length (bp)"]) <= 0:
            reasons.append("INVALID_METADATA_LENGTH")
    except (ValueError, TypeError):
        reasons.append("INVALID_METADATA_LENGTH")
    if "OUTSIDE_BALTIMORE_I_II" in reasons:
        return "exclude_primary_metadata", "outside_primary_scope", reasons
    if any(r in reasons for r in ["BELOW_PRIMARY_QUALITY", "MIUVIG_GENOME_FRAGMENT"]):
        decisive = "BELOW_PRIMARY_QUALITY" if "BELOW_PRIMARY_QUALITY" in reasons else "MIUVIG_GENOME_FRAGMENT"
        reasons = [decisive] + [r for r in reasons if r != decisive]
        return "exclude_primary_metadata", "sensitivity_candidate", reasons
    if reasons:
        return "hold_metadata", "unresolved", reasons
    return "pass_metadata_pending_ncbi", "primary_candidate", ["PENDING_NCBI_VALIDATION"]


def write_csv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def build(input_path, output_root, snapshot_id, retrieved_at, source_url):
    payload = input_path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    with input_path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != list(MAPPING):
            raise ValueError("Source schema changed; inspect columns before adapting MAPPING")
        source_rows = list(reader)
    accessions = [r["Accession"] for r in source_rows]
    if not source_rows or len(set(accessions)) != len(accessions):
        raise ValueError("Source has zero rows or duplicate accessions; no silent deduplication")
    if any(None in r or any(v is None for v in r.values()) for r in source_rows):
        raise ValueError("Malformed TSV record")
    manifest, events = [], []
    for i, raw in enumerate(source_rows, 2):
        decision, route, reasons = screen(raw)
        row = {new: raw[old] for old, new in MAPPING.items()}
        row.update({k: "" for k in EXTRA})
        row.update(
            source_snapshot_id=snapshot_id, source_retrieved_at_utc=retrieved_at,
            source_url=source_url.split("?", 1)[0], metadata_sha256=digest, source_row_number=i,
            screen_decision=decision, cohort_route=route, primary_reason_code=reasons[0],
            all_reason_codes_json=json.dumps(reasons, separators=(",", ":")),
            selection_status="pending_ncbi_validation" if route == "primary_candidate" else decision,
            in_primary_analysis="false", ncbi_fetch_status="not_started",
            host_domain_status="not_checked", sequence_qc_status="not_started",
            predictor_status="not_started",
            notes="Metadata screening only; accession.version and host domain not verified.",
        )
        manifest.append(row)
        events.append({"accession": row["accession"], "stage": "metadata_screen",
                       "decision": decision, "cohort_route": route,
                       "primary_reason_code": reasons[0],
                       "all_reason_codes_json": row["all_reason_codes_json"],
                       "scope": "primary_cohort", "snapshot_id": snapshot_id,
                       "evaluated_at_utc": retrieved_at})
    dest = output_root / "data_manifest"
    dest.mkdir(parents=True, exist_ok=True)
    write_csv(dest / "dataset_manifest.csv", list(MAPPING.values()) + EXTRA, manifest)
    write_csv(dest / "selection_flow.csv", list(events[0]), events)
    counts = collections.Counter((r["screen_decision"], r["cohort_route"], r["primary_reason_code"]) for r in manifest)
    summary = [{"screen_decision": k[0], "cohort_route": k[1], "primary_reason_code": k[2],
                "n_records": n} for k, n in sorted(counts.items())]
    write_csv(dest / "selection_summary.csv", list(summary[0]), summary)
    stats = {"snapshot_id": snapshot_id, "retrieved_at_utc": retrieved_at,
             "source_url": source_url, "source_sha256": digest,
             "source_bytes": len(payload), "n_metadata_rows": len(manifest),
             "n_unique_accessions": len(set(accessions)),
             "n_primary_candidates": sum(r["cohort_route"] == "primary_candidate" for r in manifest),
             "n_finally_included": 0, "n_genomes_gene_called": 0,
             "accession_versions_verified": 0,
             "counts_baltimore": dict(collections.Counter(r["baltimore_group"] for r in manifest)),
             "counts_checkv_quality": dict(collections.Counter(r["checkv_quality"] for r in manifest)),
             "counts_screen_decision": dict(collections.Counter(r["screen_decision"] for r in manifest)),
             "meaning": "Real source metadata; no NCBI sequence validation or biological analysis yet."}
    (dest / "source_snapshot.json").write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stats, ensure_ascii=False, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output-root", type=Path, default=Path("."))
    p.add_argument("--snapshot-id", required=True)
    p.add_argument("--retrieved-at-utc", required=True)
    p.add_argument("--source-url", default="https://phagebank.uk/inphared2/")
    a = p.parse_args()
    build(a.input, a.output_root, a.snapshot_id, a.retrieved_at_utc, a.source_url)


if __name__ == "__main__":
    main()
