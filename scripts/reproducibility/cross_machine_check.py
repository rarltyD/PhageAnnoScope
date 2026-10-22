"""Cross-machine reproducibility check: fresh caller runs on this machine must give
exactly the same CDS coordinates as the reference machine for the same genomes.

    .venv/bin/python scripts/reproducibility/cross_machine_check.py --run results/cross_machine_check \
        --reference scripts/reproducibility/reference_calls.parquet
Exit code 0 = identical; 1 = any difference (do NOT proceed with batches).
"""
import argparse
import json
import sys

import duckdb


def calls(path_glob):
    return duckdb.sql(f"""
        SELECT g.accession_version, c.source, c.strand,
               CAST(c.segments AS VARCHAR) AS segments, c.rotation_shift
        FROM read_parquet('{path_glob}/cds_calls/*.parquet', union_by_name=true) c
        JOIN (SELECT DISTINCT genome_id, accession_version
              FROM read_parquet('{path_glob}/genomes/*.parquet', union_by_name=true)) g
          USING (genome_id)
        WHERE c.source IN ('phanotate', 'pyrodigal_gv')""").df()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--reference", required=True)
    a = p.parse_args()
    new = calls(f"{a.run}/parquet")
    ref = duckdb.sql(f"SELECT * FROM read_parquet('{a.reference}')").df()
    genomes = set(ref.accession_version)
    new = new[new.accession_version.isin(genomes)]
    missing = sorted(genomes - set(new.accession_version))
    key = ["accession_version", "source", "strand", "segments", "rotation_shift"]
    a_set = set(map(tuple, new[key].astype(str).values))
    b_set = set(map(tuple, ref[key].astype(str).values))
    result = {"genomes": int(new.accession_version.nunique()), "calls_here": len(a_set),
              "calls_reference": len(b_set), "only_here": len(a_set - b_set),
              "only_reference": len(b_set - a_set), "missing_genomes": missing,
              "identical": a_set == b_set and not missing}
    print(json.dumps(result, indent=1))
    sys.exit(0 if result["identical"] else 1)


if __name__ == "__main__":
    main()
