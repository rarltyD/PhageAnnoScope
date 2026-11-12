"""Build run manifests for the four prespecified batches from the validated cohort manifest.

The GenBank file checksum is refreshed to the re-downloaded file (records may have been updated
at NCBI since 2026-10-16); the sequence checksum is kept, so the pipeline still refuses any genome
whose sequence differs from the one analysed in the study.  Writes manifests/batch1-4.csv,
manifests/reference_genomes_local.csv and a JSON report of changed or missing records.

    python scripts/pipeline/make_batch_manifests.py [REPORT_JSON]
"""
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

man = pd.read_csv("manifests/cohort_manifest_validated.csv.gz", low_memory=False, dtype=str)
report = {}
for b in (1, 2, 3, 4):
    acc = Path(f"manifests/accessions_batch{b}.txt").read_text().split()
    sub = man[man.accession_version.isin(acc)].drop_duplicates("accession_version")
    assert len(sub) == len(set(acc)), (b, len(sub), len(set(acc)))
    sub = sub.set_index("accession_version").loc[acc].reset_index()
    changed, missing = [], []
    for i, r in sub.iterrows():
        p = Path(r.genbank_path)
        if not p.exists():
            missing.append(r.accession_version)
            continue
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        if sha != r.genbank_sha256:
            changed.append(r.accession_version)
            sub.at[i, "genbank_sha256"] = sha
    sub = sub[~sub.accession_version.isin(missing)]
    sub.to_csv(f"manifests/batch{b}.csv", index=False)
    report[f"batch{b}"] = {"genomes": len(sub), "genbank_file_changed_since_2026-10-16": len(changed),
                           "not_downloaded": missing, "changed": changed}
ref = pd.read_csv("scripts/reproducibility/reference_genomes.csv", dtype=str)
for i, r in ref.iterrows():
    p = Path(r.genbank_path)
    if p.exists():
        ref.at[i, "genbank_sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
ref.to_csv("manifests/reference_genomes_local.csv", index=False)
report_path = Path(sys.argv[1] if len(sys.argv) > 1 else "logs/manifest_report.json")
report_path.parent.mkdir(parents=True, exist_ok=True)
report_path.write_text(json.dumps(report, indent=1))
print({k: (v["genomes"], v["genbank_file_changed_since_2026-10-16"], len(v["not_downloaded"])) for k, v in report.items()})
