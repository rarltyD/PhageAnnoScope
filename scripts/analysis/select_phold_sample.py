"""Stratified probability sample of proteins for Phold (protocol section 5).

Strata: support (1, 2) x length (<100, >=100 aa).  Per stratum: up to 100
proteins WITHOUT specific PHROGs function (unknown-homology / no-hit) and up to
50 WITH specific function (controls).  At most --max-per-genome proteins per
genome and stratum.  Selection order: SHA-256(seed, protein_sha256) — fixed
before Phold is run.  Inclusion probability = selected / eligible within the
stratum x function group (the per-genome cap makes it approximate; recorded).

    python scripts/analysis/select_phold_sample.py --run results/<run> --out results/<run>/phold_sample
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

SEED = 30102026
TARGETS = {"no_specific": 100, "specific": 50}


def rank(sha):
    return hashlib.sha256(f"{SEED},{sha}".encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-per-genome", type=int, default=2)
    args = parser.parse_args()
    out = Path(args.out)
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    out.mkdir(parents=True)
    loci = pd.read_parquet(Path(args.run) / "analysis" / "locus_table.parquet")
    loci = loci[loci["function_status"].isin(
        ["specific_function", "homology_unknown_function", "no_informative_hit"])]
    # One row per protein: a protein shared by loci keeps its first locus by rank.
    loci = loci.assign(
        group=lambda d: d["function_status"].map(
            lambda s: "specific" if s == "specific_function" else "no_specific"),
        length_class=lambda d: d["length_aa"].map(lambda v: "lt100" if v < 100 else "ge100"),
        order=lambda d: d["protein_sha256"].map(rank),
    ).sort_values("order").drop_duplicates("protein_sha256")
    picked = []
    for (support, length_class, group), sub in loci.groupby(["support", "length_class", "group"]):
        target = TARGETS[group]
        per_genome = {}
        chosen = []
        for row in sub.itertuples():
            if len(chosen) >= target:
                break
            if per_genome.get(row.genome_id, 0) >= args.max_per_genome:
                continue
            per_genome[row.genome_id] = per_genome.get(row.genome_id, 0) + 1
            chosen.append(row.Index)
        frame = sub.loc[chosen].assign(
            stratum=f"support{support}_{length_class}_{group}",
            eligible_in_stratum=len(sub),
            inclusion_probability=len(chosen) / len(sub) if len(sub) else None)
        picked.append(frame)
    sample = pd.concat(picked)
    proteins = pd.read_parquet(Path(args.run) / "parquet" / "proteins",
                               columns=["protein_sha256", "amino_acid_sequence"]
                               ).drop_duplicates("protein_sha256")
    sample = sample.merge(proteins, on="protein_sha256", how="left")
    assert sample["amino_acid_sequence"].notna().all()
    with (out / "phold_sample.faa").open("w") as handle:
        for row in sample.itertuples():
            handle.write(f">{row.protein_sha256}\n{row.amino_acid_sequence}\n")
    sample.drop(columns=["amino_acid_sequence"]).to_csv(out / "phold_sample.csv", index=False)
    summary = sample.groupby("stratum").agg(
        n=("protein_sha256", "size"), eligible=("eligible_in_stratum", "first"),
        inclusion_probability=("inclusion_probability", "first"),
        median_len=("length_aa", "median")).reset_index()
    summary.to_csv(out / "strata_summary.csv", index=False)
    (out / "selection.json").write_text(json.dumps({
        "seed": SEED, "targets": TARGETS, "max_per_genome": args.max_per_genome,
        "n_total": int(len(sample)), "run": str(Path(args.run).resolve())}, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
