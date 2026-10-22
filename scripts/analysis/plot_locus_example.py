"""Genome-region tracks for H4 examples (figure 11): GenBank vs PHANOTATE vs Pyrodigal-gv.

    python scripts/analysis/plot_locus_example.py --run results/full_cohort --out figures \
        --example MK972686.1:39424:39895 --example ON117153.1:3326:3674

Genes are drawn as arrows (primary rotation only); the H4 locus is highlighted.
Labels: GenBank /product for GenBank genes; PHROGs annotation for caller genes
with a specific function.  Illustrative examples, not proof of gene existence.
"""

import argparse
import json
import sys
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from make_figures import BLUE, GREY, INK, INK2, ORANGE, YELLOW, save  # noqa: E402

TRACKS = [("genbank", "GenBank", GREY), ("phanotate", "PHANOTATE", ORANGE),
          ("pyrodigal_gv", "Pyrodigal-gv", BLUE)]


def short(text, n=22):
    if not text:
        return ""
    text = text.replace("hypothetical protein", "hyp.")
    return text if len(text) <= n else text[: n - 1] + "…"


def draw(ax, con, pq, fn, acc, s0, e0, pad=2500):
    lo, hi = max(0, s0 - pad), e0 + pad
    calls = con.execute(f"""
        SELECT c.source, c.strand, c.segments, json_extract_string(c.qualifiers_json,'$.product[0]') product,
               o.protein_sha256
        FROM {pq('cds_calls')} c
        JOIN (SELECT DISTINCT genome_id, accession_version FROM {pq('genomes')}) g USING (genome_id)
        LEFT JOIN (SELECT call_id, min(protein_sha256) protein_sha256 FROM {pq('protein_origins')}
                   GROUP BY call_id) o ON o.call_id = c.call_id
        WHERE g.accession_version = ? AND coalesce(c.rotation_shift,0) = 0
          AND list_extract(c.segments,1).start0 < ? AND list_extract(c.segments,-1).end0 > ?""",
                        [acc, hi, lo]).df()
    calls = calls.merge(fn, on="protein_sha256", how="left")
    for row_i, (src, name, colr) in enumerate(TRACKS):
        y = 2 - row_i
        ax.text(lo - (hi - lo) * 0.01, y, name, ha="right", va="center", fontsize=9, color=INK)
        ax.plot([lo, hi], [y, y], color="#e6e5e0", lw=1, zorder=0)
        sub = calls[calls.source == src]
        for r in sub.itertuples():
            a = min(s["start0"] for s in r.segments)
            b = max(s["end0"] for s in r.segments)
            a, b = max(a, lo), min(b, hi)
            length = b - a
            hl = min(150, length * 0.4)
            hit = a <= e0 and b >= s0 and abs(a - s0) < 50 and abs(b - e0) < 50
            fc = YELLOW if hit else colr
            if r.strand == 1:
                arrow = FancyArrow(a, y, length, 0, width=0.32, head_width=0.32, head_length=hl,
                                   length_includes_head=True, fc=fc, ec="white", lw=0.8)
            else:
                arrow = FancyArrow(b, y, -length, 0, width=0.32, head_width=0.32, head_length=hl,
                                   length_includes_head=True, fc=fc, ec="white", lw=0.8)
            ax.add_patch(arrow)
            label = (r.product if src == "genbank" else
                     (r.best_annot if r.function_status == "specific_function" else ""))
            if label and length > 300:
                ax.text((a + b) / 2, y + (0.28 if row_i % 2 == 0 else -0.33), short(label),
                        ha="center", va="center", fontsize=6.5, color=INK2)
    ax.add_patch(Rectangle((s0, -0.45), e0 - s0, 2.9, fc=YELLOW, alpha=0.18, ec="none", zorder=0))
    ax.set_xlim(lo - (hi - lo) * 0.13, hi)
    ax.set_ylim(-0.6, 2.6)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.grid(False)
    ax.set_xlabel(f"{acc}, coordinate (nt)", fontsize=9)
    hit = calls[(calls.source != "genbank") & (calls.function_status == "specific_function")]
    hit = hit[hit.segments.map(lambda s: abs(min(x['start0'] for x in s) - s0) < 50)]
    return hit.best_annot.iloc[0] if len(hit) else "", float(hit.best_evalue.iloc[0]) if len(hit) else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--example", action="append", required=True)
    a = p.parse_args()
    run = Path(a.run)
    pq = lambda t: f"read_parquet('{(run / 'parquet' / t).as_posix()}/*.parquet', union_by_name=true)"  # noqa: E731
    con = duckdb.connect()
    con.execute("SET memory_limit='2GB'; SET threads=3")
    import pandas as pd
    fn = pd.read_parquet(run / "phrogs" / "protein_function_primary.parquet",
                         columns=["protein_sha256", "function_status", "best_annot", "best_evalue"])
    fig, axes = plt.subplots(len(a.example), 1, figsize=(10, 2.6 * len(a.example)))
    axes = [axes] if len(a.example) == 1 else axes
    info = []
    for ax, ex in zip(axes, a.example):
        acc, s0, e0 = ex.split(":")
        ann, ev = draw(ax, con, pq, fn, acc, int(s0), int(e0))
        ax.set_title(f"{acc}: ORF {int(s0) + 1}..{e0} predicted by both callers, absent from GenBank "
                     f"— PHROGs: {ann} (E = {ev:.0e})", fontsize=9.5, loc="left")
        info.append({"accession": acc, "start1": int(s0) + 1, "end": int(e0), "phrogs": ann,
                     "evalue": ev})
    fig.tight_layout()
    save(fig, Path(a.out), "fig11_h4_locus_examples")
    (run / "analysis" / "extra" / "h4_examples.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, ensure_ascii=False))


if __name__ == "__main__":
    main()
