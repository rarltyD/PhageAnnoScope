#!/usr/bin/env python3
"""Check the real frozen metadata deliverables; never run bioinformatics."""
import collections
import csv
import hashlib
import json
from pathlib import Path

from freeze_manifest import MAPPING

root = Path(__file__).resolve().parents[1]
raw = root / "data/raw/metadata_2026-10-16/inphared2_metadata.tsv"
snapshot = json.loads((root / "data_manifest/source_snapshot.json").read_text())


def read_csv(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise AssertionError(f"Malformed CSV/TSV: {path.name}")
        return rows


def main():
    source = read_csv(raw, "\t")
    manifest = read_csv(root / "data_manifest/dataset_manifest.csv")
    flow = read_csv(root / "data_manifest/selection_flow.csv")
    summary = read_csv(root / "data_manifest/selection_summary.csv")
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == snapshot["source_sha256"]
    assert len(source) == len(manifest) == len(flow) == snapshot["n_metadata_rows"] == 37387
    assert len({r["accession"] for r in manifest}) == 37387
    assert {r["Accession"] for r in source} == {r["accession"] for r in flow}
    assert len({r["accession"] for r in flow}) == len(flow)
    for i, (s, m) in enumerate(zip(source, manifest), 2):
        assert all(s[old] == m[new] for old, new in MAPPING.items())
        assert int(m["source_row_number"]) == i
        assert m["metadata_sha256"] == snapshot["source_sha256"]
        assert m["accession_version"] == ""
        assert m["in_primary_analysis"] == "false"
        assert m["predictor_status"] == "not_started"
        assert m["ncbi_fetch_status"] == "not_started"
        assert m["host_domain_status"] == "not_checked"
        assert m["primary_reason_code"] in json.loads(m["all_reason_codes_json"])
        if m["screen_decision"] == "exclude_primary_metadata":
            assert m["primary_reason_code"] in {"OUTSIDE_BALTIMORE_I_II", "BELOW_PRIMARY_QUALITY", "MIUVIG_GENOME_FRAGMENT"}
    counts = collections.Counter(r["screen_decision"] for r in manifest)
    assert counts == snapshot["counts_screen_decision"]
    assert counts == {"pass_metadata_pending_ncbi": 28664,
                      "exclude_primary_metadata": 2533, "hold_metadata": 6190}
    assert sum(int(r["n_records"]) for r in summary) == len(source)
    by_acc = {r["accession"]: r for r in manifest}
    assert by_acc["AB008550"]["cohort_route"] == "primary_candidate"
    assert by_acc["AB012573"]["cohort_route"] == "sensitivity_candidate"
    assert snapshot["n_finally_included"] == snapshot["n_genomes_gene_called"] == 0
    print(json.dumps({"status": "passed", "n_source_rows": len(source),
                      "n_preserved_metadata_fields_per_row": len(MAPPING),
                      "checks": ["raw SHA-256", "CSV/TSV structure", "lossless source columns",
                                 "unique accessions", "selection event coverage", "summary totals",
                                 "documented screening examples", "no invented versions/results"]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
