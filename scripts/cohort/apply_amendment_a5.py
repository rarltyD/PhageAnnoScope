"""Amendment A5: non-N IUPAC ambiguity codes follow the N rule (sensitivity analyses only).

Applies the rule (implemented in ``validate_cohort``) to an already validated manifest and
filters pre-drawn sample manifests WITHOUT re-selection (removed rows are listed, not replaced).

    python scripts/cohort/apply_amendment_a5.py --manifest validated_amended.csv \
        --out cohort_manifest_validated.csv --samples pilot_100.csv nested_1000.csv batch1.csv
"""

import argparse
import csv
import json
from pathlib import Path

from Bio import SeqIO

AMBIGUOUS = set("RYSWKMBDHV")


def read(path):
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def write(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--samples", nargs="*", default=[])
    args = parser.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    fields, rows = read(args.manifest)
    removed = {}
    for row in rows:
        if row.get("ncbi_fetch_status") != "success" or row["eligible_for_pilot"] != "true":
            continue
        sequence = str(next(SeqIO.parse(row["genbank_path"], "genbank")).seq).upper()
        found = {c: sequence.count(c) for c in sorted(set(sequence) & AMBIGUOUS)}
        if found:
            reasons = json.loads(row["all_reason_codes_json"]) + ["AMBIGUOUS_IUPAC_SENSITIVITY_ONLY"]
            row.update(all_reason_codes_json=json.dumps(reasons),
                       primary_reason_code="AMBIGUOUS_IUPAC_SENSITIVITY_ONLY",
                       selection_status="excluded_primary", sequence_qc_status="excluded_primary",
                       eligible_for_pilot="false",
                       amendments_applied=(row.get("amendments_applied", "") + ",A5").strip(","))
            removed[row["accession_version"]] = found
    write(out, fields, rows)
    report = {"removed_from_eligible": len(removed), "codes": removed, "samples": {}}
    for sample in map(Path, args.samples):
        sample_fields, sample_rows = read(sample)
        keep = [r for r in sample_rows if r["accession_version"] not in removed]
        write(sample.with_name(sample.stem + "_a5.csv"), sample_fields, keep)
        report["samples"][sample.stem] = {
            "before": len(sample_rows), "after": len(keep),
            "removed": sorted({r["accession_version"] for r in sample_rows} & set(removed))}
    out.with_name(out.stem + "_a5_report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps({"removed_from_eligible": len(removed),
                      "samples": {k: (v["before"], v["after"]) for k, v in report["samples"].items()}}))


if __name__ == "__main__":
    main()
