"""Publication figures (PNG 300 dpi + SVG + PDF) from <run>/analysis outputs.

    python scripts/analysis/make_figures.py --run results/<run> --out figures

Palette: reference categorical slots 1-4 (dataviz skill, documented as
CVD-validated; local validator not run because Node.js is unavailable).
All CIs: relatedness-block bootstrap (seed 30102026).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from fast_stats import cluster_bootstrap  # noqa: E402  (vectorised; same estimand)

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GREY, INK, INK2 = "#bdbcb5", "#0b0b0b", "#52514e"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "figure.dpi": 100,
})
LABEL = {"genbank": "GenBank", "phanotate": "PHANOTATE", "pyrodigal_gv": "Pyrodigal-gv"}


def save(fig, out, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"{name}.{ext}", dpi=300 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_agreement(run, out):
    m = pd.read_csv(run / "analysis" / "pairwise_agreement.csv")
    pairs = [f"{LABEL[a]}\n vs {LABEL[b]}" for a, b in zip(m.source_a, m.source_b)]
    x = np.arange(len(m))
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    w = 0.36
    b1 = ax.bar(x - w / 2 - 0.01, m.pooled_exact_jaccard, w, color=BLUE, label="Exact CDS")
    b2 = ax.bar(x + w / 2 + 0.01, m.pooled_stop_jaccard, w, color=AQUA, label="Shared stop")
    for bars in (b1, b2):
        ax.bar_label(bars, fmt="%.2f", fontsize=7, color=INK2, padding=2)
    ax.set_xticks(x, pairs, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Pooled Jaccard index")
    ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.12))
    save(fig, out, "fig1_pairwise_agreement")


def best_relations(run):
    duckdb.sql("SET memory_limit='2GB'; SET threads=3")
    p = lambda t: f"read_parquet('{(run / 'parquet' / t).as_posix()}/*.parquet', union_by_name=true)"  # noqa: E731
    order = ("CASE relation WHEN 'exact_cds' THEN 1 WHEN 'alternative_start' THEN 2 "
             "WHEN 'contained_distinct_stop' THEN 3 WHEN 'partial_overlap' THEN 4 "
             "WHEN 'antisense_overlap' THEN 5 WHEN 'no_correspondence' THEN 6 ELSE 7 END")
    q = f"""
    WITH rows AS (
      SELECT source_a AS src, source_b AS other, call_id_a AS call_id, relation,
             substantial_overlap, comparable FROM {p('matches')} WHERE call_id_a IS NOT NULL
      UNION ALL
      SELECT source_b, source_a, call_id_b, relation, substantial_overlap, comparable
      FROM {p('matches')} WHERE call_id_b IS NOT NULL),
    ranked AS (
      SELECT *, CASE WHEN relation='partial_overlap' AND substantial_overlap
                     THEN 'substantial_overlap' WHEN relation='partial_overlap'
                     THEN 'minor_overlap_only' ELSE relation END AS rel2,
             row_number() OVER (PARTITION BY src, other, call_id ORDER BY {order},
                                substantial_overlap DESC) rn
      FROM rows)
    SELECT src, other, rel2, count(*) n FROM ranked WHERE rn=1 GROUP BY ALL"""
    return duckdb.sql(q).df()


def fig_relations(run, out):
    d = best_relations(run)
    d.to_csv(run / "analysis" / "best_relation_per_call.csv", index=False)
    d = d[d.rel2 != "minor_overlap_only"].copy()
    d.loc[d.rel2.isin(["antisense_overlap"]), "rel2"] = "no_correspondence"
    cats = ["exact_cds", "alternative_start", "contained_distinct_stop", "substantial_overlap",
            "no_correspondence", "not_comparable"]
    names = ["Exact CDS", "Alternative start", "Nested, other stop", "Partial (≥80% reciprocal)",
             "No correspondence", "Not comparable"]
    colors = [BLUE, AQUA, YELLOW, "#e87ba4", ORANGE, GREY]
    piv = d.pivot_table(index=["src", "other"], columns="rel2", values="n", aggfunc="sum",
                        fill_value=0).reindex(columns=cats, fill_value=0)
    frac = piv.div(piv.sum(axis=1), axis=0)
    labels = [f"{LABEL[a]} → {LABEL[b]}" for a, b in frac.index]
    fig, ax = plt.subplots(figsize=(6.2, 2.9))
    left = np.zeros(len(frac))
    for c, name, col in zip(cats, names, colors):
        ax.barh(labels, frac[c], left=left, color=col, label=name, edgecolor="white",
                linewidth=1)
        for i, v in enumerate(frac[c]):
            if v > 0.06:
                ax.text(left[i] + v / 2, i, f"{v:.0%}", ha="center", va="center", fontsize=7,
                        color="white" if col in (BLUE, ORANGE) else INK)
        left += frac[c].values
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.set_xlabel("Share of CDS of the first source (best relation to second source)")
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=3, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.22))
    save(fig, out, "fig2_discrepancy_types")


def fig_h2(df, out, reps):
    bins = [0, 50, 100, 150, 200, 300, 500, 100000]
    names = ["<50", "50–99", "100–149", "150–199", "200–299", "300–499", "≥500"]
    d = df.assign(lbin=pd.cut(df.length_aa, bins, right=False, labels=names))
    rows = []
    for split, sub0 in (("No overlap", d[d.overlap_any_ge1 == 0]),
                        ("Antisense overlap", d[d.overlap_anti_ge1 == 1])):
        for name in names:
            sub = sub0[sub0.lbin == name]
            if len(sub) < 20:
                continue
            b = cluster_bootstrap(sub.assign(all=1), "single_support", "all", reps=reps)
            lv = b["levels"]["1"]
            rows.append((split, name, lv["proportion"], *lv["ci95"], lv["n"]))
    t = pd.DataFrame(rows, columns=["split", "bin", "p", "lo", "hi", "n"])
    t.to_csv(out / "fig3_h2_data.csv", index=False)
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    for split, col in (("No overlap", BLUE), ("Antisense overlap", ORANGE)):
        s = t[t.split == split]
        xs = [names.index(b) for b in s.bin]
        ax.errorbar(xs, s.p, yerr=[s.p - s.lo, s.hi - s.p], color=col, marker="o", ms=4,
                    lw=2, capsize=2, label=split)
    ax.set_xticks(range(len(names)), names, fontsize=7.5, rotation=20)
    ax.set_xlabel("Protein length, aa (representative model)")
    ax.set_ylabel("Supported by one caller only")
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.set_ylim(0, 1)
    ax.legend()
    save(fig, out, "fig3_h2_length_overlap")


def fig_h1(df, out, reps):
    d = df[df.function_status != "not_searched_or_failed"].assign(
        no_specific=lambda x: (x.function_status != "specific_function").astype(int),
        lclass=lambda x: np.where(x.length_aa < 100, "<100 aa", "≥100 aa"))
    rows = []
    for lc in ("<100 aa", "≥100 aa"):
        b = cluster_bootstrap(d[d.lclass == lc], "no_specific", "support", reps=reps)
        for sup, lv in b["levels"].items():
            rows.append((lc, int(sup), lv["proportion"], *lv["ci95"], lv["n"]))
    t = pd.DataFrame(rows, columns=["length", "support", "p", "lo", "hi", "n"])
    t.to_csv(out / "fig4_h1_data.csv", index=False)
    comp = d.groupby(["support", "function_status"]).size().unstack(fill_value=0)
    comp = comp.div(comp.sum(axis=1), axis=0)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1.1]})
    ax = axes[0]
    x = np.arange(2)
    for i, (sup, col) in enumerate(((1, ORANGE), (2, BLUE))):
        s = t[t.support == sup]
        bars = ax.bar(x + (i - 0.5) * 0.38, s.p, 0.36, color=col,
                      label=f"{sup} caller{'s' if sup == 2 else ''}",
                      yerr=[s.p - s.lo, s.hi - s.p], capsize=2, ecolor=INK2)
        ax.bar_label(bars, labels=[f"{v:.0%}" for v in s.p], fontsize=7, color=INK2,
                     padding=6)
    ax.set_xticks(x, ["<100 aa", "≥100 aa"])
    ax.set_ylim(0, 1.25)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.set_ylabel("No specific PHROGs function")
    ax.legend(title="Support", loc="upper right", fontsize=7, title_fontsize=7, ncol=2)
    ax = axes[1]
    cats = ["specific_function", "ambiguous", "homology_unknown_function", "no_informative_hit"]
    names = ["Specific function", "Ambiguous", "PHROG of unknown function", "No informative hit"]
    cols = [BLUE, YELLOW, AQUA, GREY]
    left = np.zeros(len(comp))
    labels = [f"Support {s}" for s in comp.index]
    for c, name, col in zip(cats, names, cols):
        v = comp.get(c, pd.Series(0, index=comp.index)).values
        ax.barh(labels, v, left=left, color=col, label=name, edgecolor="white", linewidth=1)
        for i, vi in enumerate(v):
            if vi > 0.07:
                ax.text(left[i] + vi / 2, i, f"{vi:.0%}", ha="center", va="center", fontsize=7,
                        color="white" if col == BLUE else INK)
        left += v
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.45, -0.15))
    fig.tight_layout()
    save(fig, out, "fig4_h1_function_by_support")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--reps", type=int, default=1000)
    args = parser.parse_args()
    run, out = Path(args.run), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(run / "analysis" / "locus_table.parquet",
                         columns=["length_aa", "overlap_any_ge1", "overlap_anti_ge1",
                                  "single_support", "block", "function_status", "support"])
    fig_agreement(run, out)
    fig_relations(run, out)
    fig_h2(df, out, args.reps)
    if "function_status" in df:
        fig_h1(df, out, args.reps)
    (out / "figures_info.json").write_text(json.dumps(
        {"run": str(run.resolve()), "reps": args.reps, "seed": 30102026}, indent=2))
    print(sorted(p.name for p in out.glob("*.png")))


if __name__ == "__main__":
    main()
