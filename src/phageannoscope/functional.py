"""PHROGs HMM search over the deduplicated protein catalogue (protocol section 5).

One search per exact protein sequence (``protein_sha256``); results are joined
back to models/loci via ``protein_origins``.  Proteins are processed in chunks;
each finished chunk is an immutable Parquet file, so an interrupted run resumes
without recomputing completed chunks.

Status rules (prespecified in ``configs/protocol_v1.json``; thresholds are
investigator choices, not validated biological truth):

* a hit *passes* when full-sequence E <= e_max, query coverage >= qcov_min and
  profile coverage >= pcov_min (coverages are unions of domain envelopes /
  HMM intervals of included domains);
* ``specific_function``: best passing hit belongs to a PHROG whose category is
  not ``unknown function``;
* ``homology_unknown_function``: passing hits exist, best is an unknown-function
  PHROG;
* ``ambiguous``: another passing hit with bitscore >= 0.9 x best bitscore has a
  different *known* category than the best hit;
* ``no_informative_hit``: no passing hit;
* ``not_searched_or_failed`` is assigned by callers for proteins never searched.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

UNKNOWN = "unknown function"
NEAR_TIE = 0.9

HIT_SCHEMA = pa.schema(
    [
        ("protein_sha256", pa.string()),
        ("phrog", pa.string()),
        ("evalue", pa.float64()),
        ("bitscore", pa.float64()),
        ("best_domain_ievalue", pa.float64()),
        ("query_coverage", pa.float64()),
        ("profile_coverage", pa.float64()),
        ("n_domains_included", pa.int64()),
        ("length_aa", pa.int64()),
        ("profile_length", pa.int64()),
    ]
)


def _union_length(intervals):
    total, current = 0, None
    for start, end in sorted(intervals):
        if current is None or start > current[1]:
            if current:
                total += current[1] - current[0] + 1
            current = [start, end]
        else:
            current[1] = max(current[1], end)
    if current:
        total += current[1] - current[0] + 1
    return total


def load_annotations(path):
    table = {}
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            table["phrog_" + row["phrog"]] = {
                "annot": row["annot"] if row["annot"] not in ("", "NA") else None,
                "category": row["category"] or UNKNOWN,
            }
    return table


def search_chunk(sequences, hmm_db, cpus, e_report):
    """sequences: list of (sha, aa). Returns list of hit dicts (E <= e_report)."""
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence
    from pyhmmer.plan7 import HMMFile

    alphabet = Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            TextSequence(name=sha.encode(), sequence=aa).digitize(alphabet)
            for sha, aa in sequences
        ],
    )
    lengths = {sha: len(aa) for sha, aa in sequences}
    hits = []
    with HMMFile(hmm_db) as handle:
        profiles = handle.optimized_profiles() if handle.is_pressed() else handle
        for top in pyhmmer.hmmsearch(profiles, block, cpus=cpus, E=e_report, domE=e_report):
            phrog = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
            model_len = top.query.M
            for hit in top:
                if not hit.reported:
                    continue
                name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                domains = [d for d in hit.domains if d.included] or list(hit.domains)
                q_iv = [(d.env_from, d.env_to) for d in domains]
                p_iv = [(d.alignment.hmm_from, d.alignment.hmm_to) for d in domains]
                hits.append(
                    {
                        "protein_sha256": name,
                        "phrog": phrog,
                        "evalue": float(hit.evalue),
                        "bitscore": float(hit.score),
                        "best_domain_ievalue": float(min(d.i_evalue for d in domains)),
                        "query_coverage": _union_length(q_iv) / lengths[name],
                        "profile_coverage": _union_length(p_iv) / model_len,
                        "n_domains_included": sum(d.included for d in hit.domains),
                        "length_aa": lengths[name],
                        "profile_length": model_len,
                    }
                )
    return hits


def classify(hits, annotations, e_max=1e-5, qcov_min=0.5, pcov_min=0.5):
    """Return one status row for a single protein given all its raw hits."""
    passing = sorted(
        (
            h
            for h in hits
            if h["evalue"] <= e_max
            and h["query_coverage"] >= qcov_min
            and h["profile_coverage"] >= pcov_min
        ),
        key=lambda h: (-h["bitscore"], h["evalue"], h["phrog"]),
    )
    if not passing:
        return {"function_status": "no_informative_hit", "best_phrog": None,
                "best_category": None, "best_annot": None, "best_evalue": None,
                "best_bitscore": None, "n_passing_hits": 0}
    best = passing[0]
    info = annotations.get(best["phrog"], {"annot": None, "category": UNKNOWN})
    status = "homology_unknown_function" if info["category"] == UNKNOWN else "specific_function"
    near = {
        annotations.get(h["phrog"], {}).get("category", UNKNOWN)
        for h in passing
        if h["bitscore"] >= NEAR_TIE * best["bitscore"]
    } - {UNKNOWN}
    if status == "specific_function" and len(near) > 1:
        status = "ambiguous"
    return {
        "function_status": status,
        "best_phrog": best["phrog"],
        "best_category": info["category"],
        "best_annot": info["annot"],
        "best_evalue": best["evalue"],
        "best_bitscore": best["bitscore"],
        "n_passing_hits": len(passing),
    }


def run_search(proteins_parquet, hmm_db, out_dir, *, chunk_size=20000, cpus=4,
               e_report=1e-3, limit=0):
    """Checkpointed search; returns a summary dict. Never overwrites a chunk."""
    out = Path(out_dir)
    (out / "hits").mkdir(parents=True, exist_ok=True)
    table = pq.read_table(proteins_parquet, columns=["protein_sha256", "amino_acid_sequence"])
    rows = sorted(dict(zip(table.column(0).to_pylist(), table.column(1).to_pylist())).items())
    if limit:
        rows = rows[:limit]
    db_sha = hashlib.sha256(Path(hmm_db).read_bytes()).hexdigest()
    log = out / "chunks.jsonl"
    done = 0
    for i in range(0, len(rows), chunk_size):
        chunk = rows[i : i + chunk_size]
        key = hashlib.sha256(
            json.dumps([db_sha, e_report, [s for s, _ in chunk]]).encode()
        ).hexdigest()[:16]
        path = out / "hits" / f"chunk_{i // chunk_size:05d}_{key}.parquet"
        if path.exists() and path.with_name(path.stem + ".members.json").exists():
            done += len(chunk)
            continue
        started = time.perf_counter()
        hits = search_chunk(chunk, hmm_db, cpus, e_report)
        tmp = path.with_suffix(".tmp")
        pq.write_table(pa.Table.from_pylist(hits, schema=HIT_SCHEMA), tmp)
        members = path.with_name(path.stem + ".members.json")
        members.write_text(json.dumps([s for s, _ in chunk]))
        tmp.replace(path)
        done += len(chunk)
        with log.open("a") as handle:
            handle.write(json.dumps({"chunk": path.name, "proteins": len(chunk),
                                     "hits": len(hits), "cpus": cpus,
                                     "wall_seconds": round(time.perf_counter() - started, 2),
                                     "at": time.strftime("%Y-%m-%dT%H:%M:%S")}) + "\n")
        print(f"chunk {path.name}: {len(chunk)} proteins, {len(hits)} hits, "
              f"{time.perf_counter() - started:.1f}s; {done}/{len(rows)}", flush=True)
    return {"proteins": len(rows), "hmm_db_sha256": db_sha, "out": str(out)}


def summarize(proteins_parquet, out_dir, annotation_tsv, *, e_max=1e-5, qcov_min=0.5,
              pcov_min=0.5, suffix="primary"):
    """Write protein_function_<suffix>.parquet with one row per searched protein."""
    out = Path(out_dir)
    annotations = load_annotations(annotation_tsv)
    by_protein = {}
    for path in sorted((out / "hits").glob("chunk_*.parquet")):
        for hit in pq.read_table(path).to_pylist():
            by_protein.setdefault(hit["protein_sha256"], []).append(hit)
    table = pq.read_table(proteins_parquet, columns=["protein_sha256", "length_aa"])
    chunk_files = {p.name for p in (out / "hits").glob("chunk_*.parquet")}
    if not chunk_files:
        raise ValueError("No completed search chunks")
    rows = []
    covered = _covered_shas(out)
    unique = dict(zip(table.column(0).to_pylist(), table.column(1).to_pylist()))
    for sha, length in sorted(unique.items()):
        if sha not in covered:
            row = {"function_status": "not_searched_or_failed", "best_phrog": None,
                   "best_category": None, "best_annot": None, "best_evalue": None,
                   "best_bitscore": None, "n_passing_hits": 0}
        else:
            row = classify(by_protein.get(sha, []), annotations, e_max, qcov_min, pcov_min)
        rows.append({"protein_sha256": sha, "length_aa": length, **row,
                     "e_max": e_max, "qcov_min": qcov_min, "pcov_min": pcov_min})
    target = out / f"protein_function_{suffix}.parquet"
    pq.write_table(pa.Table.from_pylist(rows), target)
    return target


def _covered_shas(out):
    """Proteins belonging to completed chunks (recorded in the sidecar list)."""
    covered = set()
    for path in (out / "hits").glob("chunk_*.members.json"):
        covered.update(json.loads(path.read_text()))
    return covered
