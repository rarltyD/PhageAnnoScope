"""Rebuild run tables from verified per-genome checkpoints of an interrupted run.

Used for batch 3, computed on a second machine: 4 963/5 000 genomes were checkpointed,
but the run could not be resumed on the primary machine because the run identity includes
the WSL kernel string (code, versions and caller settings identical).  The
checkpoint JSON is exactly the per-genome result the pipeline would have written,
so feeding it through the same ResultWriter reproduces the run tables without
recomputation.  Genomes of the manifest without a checkpoint are written to a
'missing' manifest for a separate small run.

    python scripts/pipeline/rebuild_from_checkpoints.py --run results/batch3_checkpointed \
        --manifest manifests/batch3.csv --out results/batch3
"""

import argparse
import csv
import json
import shutil
from pathlib import Path

from phageannoscope.storage import ResultWriter, build_db, export_proteins
from phageannoscope.util import file_sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    run, out = Path(a.run), Path(a.out)
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    out.mkdir(parents=True)
    shutil.copy2(run / "run_config.json", out / "source_run_config.json")
    writer = ResultWriter(out / "parquet", 20)
    done, bad = set(), []
    for ck in sorted((run / "checkpoints").glob("*.json")):
        sha = ck.with_suffix(".sha256")
        if not sha.exists() or sha.read_text() != file_sha(ck):
            bad.append(ck.name)
            continue
        result = json.loads(ck.read_text(encoding="utf-8"))
        writer.add(result)
        done.update(g["accession_version"] for g in result["genomes"])
    counts = writer.close()
    db = out / "phageannoscope.duckdb"
    build_db(out / "parquet", db)
    export_proteins(db, out / "proteins.faa")
    with open(a.manifest, encoding="utf-8") as h:
        reader = csv.DictReader(h)
        fields, rows = reader.fieldnames, list(reader)
    missing = [r for r in rows if r["accession_version"] not in done]
    miss_path = out / "missing_manifest.csv"
    with miss_path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        w.writerows(missing)
    summary = {"source_run": str(run.resolve()), "manifest": a.manifest,
               "checkpoints_used": len(done), "checkpoints_rejected": bad,
               "missing_genomes": len(missing), "missing_manifest": str(miss_path),
               "table_counts": counts}
    (out / "rebuild_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "table_counts"}, indent=1))


if __name__ == "__main__":
    main()
