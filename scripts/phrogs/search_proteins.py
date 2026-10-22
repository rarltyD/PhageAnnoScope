"""PHROGs HMM search for a run, via a study-wide protein-level cache.

    nice -n 10 python scripts/phrogs/search_proteins.py --run results/batch1 \
        --hmm data/db/phrogs_v4/phrogs_v4_m50.h3m --cpus 3

Every unique protein (protein_sha256) is searched once for the whole study:
raw hits (E<=1e-3) live in data/db/phrogs_cache/hits/*.parquet with a
members list per chunk.  A run searches only proteins not yet in the cache,
then writes <run>/phrogs/protein_function_<variant>.parquet for the
prespecified thresholds and sensitivity variants.  Cache identity (HMM SHA-256,
report E-value) is enforced.  ``--seed-from`` imports finished chunks of an
earlier run-local search (same HMM SHA) into the cache.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from phageannoscope.functional import HIT_SCHEMA, classify, load_annotations, search_chunk

SENSITIVITY = {
    "primary": (1e-5, 0.5, 0.5),
    "e1e-3": (1e-3, 0.5, 0.5),
    "e1e-10": (1e-10, 0.5, 0.5),
    "cov0.3": (1e-5, 0.3, 0.3),
    "cov0.7": (1e-5, 0.7, 0.7),
}
E_REPORT = 1e-3


def file_sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cache_init(cache, hmm):
    (cache / "hits").mkdir(parents=True, exist_ok=True)
    meta_path = cache / "cache_meta.json"
    meta = {"hmm": str(Path(hmm).resolve()), "hmm_sha256": file_sha(hmm), "e_report": E_REPORT}
    if meta_path.exists():
        old = json.loads(meta_path.read_text())
        if old["hmm_sha256"] != meta["hmm_sha256"] or old["e_report"] != E_REPORT:
            raise SystemExit("PHROGs cache identity mismatch (different HMM or E-value)")
    else:
        meta_path.write_text(json.dumps(meta, indent=2))
    return meta


def covered(cache):
    done = set()
    for path in (cache / "hits").glob("*.members.json"):
        if path.with_name(path.name.replace(".members.json", ".parquet")).exists():
            done.update(json.loads(path.read_text()))
    return done


def seed(cache, source):
    n = 0
    for members in Path(source).glob("chunk_*.members.json"):
        hits = members.with_name(members.name.replace(".members.json", ".parquet"))
        target = cache / "hits" / ("seed_" + hits.name)
        if hits.exists() and not target.exists():
            shutil.copy2(hits, target)
            shutil.copy2(members, target.with_name(target.name.replace(".parquet",
                                                                       ".members.json")))
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--hmm", required=True)
    parser.add_argument("--cache", default="data/db/phrogs_cache")
    parser.add_argument("--seed-from", nargs="*", default=[])
    parser.add_argument("--annotations", default="data/db/phrogs_v4/phrog_annot_v4.tsv")
    parser.add_argument("--cpus", type=int, default=3)
    parser.add_argument("--chunk-size", type=int, default=20000)
    args = parser.parse_args()
    import os
    limit = int(os.environ.get("PHAGEANNOSCOPE_MAX_CPUS", "4"))  # default safety limit (override for workstations)
    if args.cpus > limit:
        raise SystemExit(f"Project limit: at most {limit} CPU threads")
    run, cache = Path(args.run), Path(args.cache)
    meta = cache_init(cache, args.hmm)
    for source in args.seed_from:
        print(f"seeded {seed(cache, source)} chunks from {source}", flush=True)
    proteins_dir = run / "parquet" / "proteins"
    table = pq.read_table(proteins_dir, columns=["protein_sha256", "amino_acid_sequence",
                                                 "length_aa"])
    proteins = table.to_pandas().drop_duplicates("protein_sha256")
    done = covered(cache)
    todo = proteins[~proteins.protein_sha256.isin(done)].sort_values("protein_sha256")
    print(f"run proteins {len(proteins)}; cached {len(proteins) - len(todo)}; "
          f"to search {len(todo)}", flush=True)
    started = time.perf_counter()
    log = cache / "chunks.jsonl"
    rows = list(zip(todo.protein_sha256, todo.amino_acid_sequence))
    for i in range(0, len(rows), args.chunk_size):
        chunk = rows[i : i + args.chunk_size]
        key = hashlib.sha256(json.dumps([meta["hmm_sha256"], E_REPORT,
                                         [s for s, _ in chunk]]).encode()).hexdigest()[:16]
        path = cache / "hits" / f"chunk_{key}.parquet"
        t0 = time.perf_counter()
        hits = search_chunk(chunk, args.hmm, args.cpus, E_REPORT)
        tmp = path.with_suffix(".tmp")
        pq.write_table(pa.Table.from_pylist(hits, schema=HIT_SCHEMA), tmp)
        path.with_name(path.stem + ".members.json").write_text(json.dumps([s for s, _ in chunk]))
        tmp.replace(path)
        with log.open("a") as h:
            h.write(json.dumps({"chunk": path.name, "proteins": len(chunk), "hits": len(hits),
                                "cpus": args.cpus, "run": str(run),
                                "wall_seconds": round(time.perf_counter() - t0, 2),
                                "at": time.strftime("%Y-%m-%dT%H:%M:%S")}) + "\n")
        print(f"chunk {path.name}: {len(chunk)} proteins, {len(hits)} hits, "
              f"{time.perf_counter() - t0:.1f}s; {min(i + args.chunk_size, len(rows))}/"
              f"{len(rows)}", flush=True)
    # Summaries for this run's proteins only (hits filtered at read time).
    wanted = set(proteins.protein_sha256)
    done = covered(cache)
    hit_files = [str(p) for p in sorted((cache / "hits").glob("*.parquet"))]
    hits = ds.dataset(hit_files, schema=HIT_SCHEMA).to_table(
        filter=ds.field("protein_sha256").isin(pa.array(sorted(wanted)))).to_pandas()
    hits = hits.drop_duplicates(["protein_sha256", "phrog"])
    grouped = {k: g.to_dict("records") for k, g in hits.groupby("protein_sha256")}
    annotations = load_annotations(args.annotations)
    out = run / "phrogs"
    out.mkdir(exist_ok=True)
    outputs = {}
    for name, (e, q, p) in SENSITIVITY.items():
        recs = []
        for sha, length in zip(proteins.protein_sha256, proteins.length_aa):
            if sha in done:
                row = classify(grouped.get(sha, []), annotations, e, q, p)
            else:
                row = {"function_status": "not_searched_or_failed", "best_phrog": None,
                       "best_category": None, "best_annot": None, "best_evalue": None,
                       "best_bitscore": None, "n_passing_hits": 0}
            recs.append({"protein_sha256": sha, "length_aa": length, **row,
                         "e_max": e, "qcov_min": q, "pcov_min": p})
        target = out / f"protein_function_{name}.parquet"
        pd.DataFrame(recs).to_parquet(target, index=False)
        outputs[name] = str(target)
    info = {"run": str(run), "proteins": int(len(proteins)), "searched_now": int(len(todo)),
            "cache": str(cache.resolve()), "hmm_sha256": meta["hmm_sha256"],
            "wall_seconds": round(time.perf_counter() - started, 1), "cpus": args.cpus,
            "outputs": outputs}
    (out / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
