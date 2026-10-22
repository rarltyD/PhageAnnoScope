"""Combine completed batch runs into one analysis view (no recomputation, no copies).

Genome sets of the batches are disjoint by construction (nested select, new
genomes only).  Each table directory of the merged view holds symlinks to the
batch shards (renamed <run>__part-*.parquet); analysis queries use DISTINCT/GROUP
BY on identities, so proteins/models shared between batches do not double count
genome-level units.  PHROGs function tables are concatenated and de-duplicated by
protein_sha256 (identical HMM SHA enforced).

    python scripts/pipeline/merge_runs.py --out results/full_cohort --runs results/batch1 \
        results/batch2 results/batch3 results/batch4
"""

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

VARIANTS = ("primary", "e1e-3", "e1e-10", "cov0.3", "cov0.7")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--runs", nargs="+", required=True)
    a = p.parse_args()
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    (out / "parquet").mkdir(parents=True)
    hmm, genomes, info = set(), {}, {"runs": []}
    for run in map(Path, a.runs):
        tag = run.name
        for table_dir in sorted((run / "parquet").iterdir()):
            if not table_dir.is_dir():
                continue
            target = out / "parquet" / table_dir.name
            target.mkdir(exist_ok=True)
            for part in sorted(table_dir.glob("*.parquet")):
                (target / f"{tag}__{part.name}").symlink_to(part.resolve())
        g = duckdb.sql(f"SELECT DISTINCT accession_version FROM read_parquet("
                       f"'{(run / 'parquet' / 'genomes').as_posix()}/*.parquet')").df()
        for acc in g.accession_version:
            if acc in genomes:
                raise SystemExit(f"Genome {acc} in both {genomes[acc]} and {tag}")
            genomes[acc] = tag
        ri = json.loads((run / "phrogs" / "run_info.json").read_text())
        hmm.add(ri["hmm_sha256"])
        info["runs"].append({"run": str(run.resolve()), "genomes": int(len(g)),
                             "phrogs_hmm_sha256": ri["hmm_sha256"]})
    if len(hmm) != 1:
        raise SystemExit(f"Different PHROGs HMM across runs: {hmm}")
    (out / "phrogs").mkdir()
    for v in VARIANTS:
        frames = [pd.read_parquet(Path(r) / "phrogs" / f"protein_function_{v}.parquet")
                  for r in a.runs]
        merged = pd.concat(frames).drop_duplicates("protein_sha256")
        merged.to_parquet(out / "phrogs" / f"protein_function_{v}.parquet", index=False)
    info.update(total_genomes=len(genomes), hmm_sha256=hmm.pop())
    (out / "phrogs" / "run_info.json").write_text(json.dumps({"hmm_sha256": info["hmm_sha256"]}))
    (out / "merge_info.json").write_text(json.dumps(info, indent=2))
    print(json.dumps({"total_genomes": info["total_genomes"],
                      "runs": [(r["run"].split("/")[-1], r["genomes"]) for r in info["runs"]]}))


if __name__ == "__main__":
    main()
