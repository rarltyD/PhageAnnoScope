"""Select 30 loci for manual review (protocol section 7) and dump raw evidence.

10 random primary loci + 4 each of: exact agreement (with GenBank exact),
alternative start, single-caller locus overlapping another locus, additional ORF
vs GenBank (H4 baseline_clean), circular/edge/complex context.  Selection is by
SHA-256(seed, locus_id) order, fixed before inspection.  The sheet lists every
source's calls within the locus window so a reviewer can check coordinates
before looking at function.  Decision columns are left EMPTY: a locus is
'reviewed' only after a real inspection.

    python scripts/analysis/select_manual_review.py --run results/<run>
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import duckdb
import pandas as pd

SEED = 30102026


def order(v):
    return hashlib.sha256(f"{SEED},{v}".encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    run = Path(args.run)
    out = Path(args.out or run / "analysis" / "manual_review.csv")
    if out.exists():
        raise SystemExit(f"Refusing to overwrite {out}")
    df = pd.read_parquet(run / "analysis" / "locus_table.parquet")
    df = df.assign(order=df.locus_id.map(order)).sort_values("order")
    h4 = pd.read_csv(run / "analysis" / "h4_candidates.csv")
    cats = {
        "random": df,
        "exact_agreement": df[(df.support == 2) & df.genbank_exact_match.fillna(False)],
        "single_caller_overlapping": df[(df.support == 1) & (df.overlap_any_ge1 == 1)],
        "alternative_start_disagreement": df[(df.support == 2) & df.genbank_stop_match
                                             .fillna(False) & ~df.genbank_exact_match
                                             .fillna(False)],
        "additional_orf_h4": df[df.locus_id.isin(h4[h4.h4_layer == "baseline_clean"].locus_id)],
        "complex_or_edge_context": df[df.genbank_complex_context.fillna(False)
                                      | df.genbank_partial_context.fillna(False)],
    }
    quota = {"random": 10, "exact_agreement": 4, "single_caller_overlapping": 4,
             "alternative_start_disagreement": 4, "additional_orf_h4": 4,
             "complex_or_edge_context": 4}
    picked, seen = [], set()
    for cat, sub in cats.items():
        n = 0
        for row in sub.itertuples():
            if n >= quota[cat]:
                break
            if row.locus_id in seen:
                continue
            seen.add(row.locus_id)
            picked.append((cat, row))
            n += 1
    p = lambda t: f"read_parquet('{(run / 'parquet' / t).as_posix()}/*.parquet', union_by_name=true)"  # noqa: E731
    con = duckdb.connect()
    rows = []
    for cat, r in picked:
        lo, hi = int(r.start0) - 300, int(r.end0) + 300
        calls = con.execute(f"""
            SELECT source, strand, segments, feature_nt_length, eligible_primary,
                   primary_exclusion_reasons
            FROM {p('cds_calls')}
            WHERE genome_id = ? AND coalesce(rotation_shift, 0) = 0
              AND list_extract(segments, 1).start0 < ? AND list_extract(segments, -1).end0 > ?
            ORDER BY source, list_extract(segments, 1).start0""",
                            [r.genome_id, hi, lo]).df()
        tracks = [
            f"{c.source}:{'+' if c.strand == 1 else '-'}"
            f"{','.join(f'{s['start0'] + 1}..{s['end0']}' for s in c.segments)}"
            f"{'' if c.eligible_primary else ' [excl:' + str(c.primary_exclusion_reasons) + ']'}"
            for c in calls.itertuples()]
        rows.append({
            "review_category": cat, "accession_version": r.accession_version,
            "locus_id": r.locus_id, "strand": int(r.strand),
            "representative_1based": f"{int(r.start0) + 1}..{int(r.end0)}",
            "length_aa": int(r.length_aa), "support": int(r.support),
            "phanotate": bool(r.phanotate_present), "pyrodigal_gv": bool(r.pyrodigal_gv_present),
            "genbank_exact": r.genbank_exact_match, "genbank_stop": r.genbank_stop_match,
            "tracks_window": " | ".join(tracks),
            "phrogs_status": getattr(r, "function_status", None),
            "phrogs_best": getattr(r, "best_phrog", None),
            "phrogs_annot": getattr(r, "best_annot", None),
            "phrogs_evalue": getattr(r, "best_evalue", None),
            "coordinate_check": "", "function_evidence_check": "", "reviewer": "",
            "reviewed_at": "", "notes": ""})
    pd.DataFrame(rows).to_csv(out, index=False)
    print(json.dumps(pd.Series([c for c, _ in picked]).value_counts().to_dict()), out)


if __name__ == "__main__":
    main()
