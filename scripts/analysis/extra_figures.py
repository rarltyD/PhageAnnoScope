"""Figures 5-8 from <run>/analysis/extra (descriptive; for the 20-min talk).

    python scripts/analysis/extra_figures.py --run results/full_cohort --out figures
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from make_figures import AQUA, BLUE, GREY, INK, INK2, LABEL, ORANGE, YELLOW, save  # noqa: E402


def fig5_genbank_products(extra, out):
    t = pd.read_csv(extra / "genbank_product_vs_phrogs.csv", index_col=0)
    t = t.loc[[i for i in ("named", "hypothetical_like") if i in t.index]]
    cats = ["specific_function", "homology_unknown_function", "no_informative_hit"]
    names = ["Specific PHROGs function", "PHROG of unknown function", "No informative hit"]
    cols = [BLUE, AQUA, GREY]
    frac = t.reindex(columns=cats, fill_value=0)
    frac = frac.div(frac.sum(axis=1), axis=0)
    labels = {"named": f"Named in GenBank\n(n={int(t.loc['named'].sum()):,})",
              "hypothetical_like": f"'Hypothetical' in GenBank\n"
                                   f"(n={int(t.loc['hypothetical_like'].sum()):,})"}
    fig, ax = plt.subplots(figsize=(6.4, 2.4))
    left = np.zeros(len(frac))
    ylab = [labels[i] for i in frac.index]
    for c, n, col in zip(cats, names, cols):
        ax.barh(ylab, frac[c], left=left, color=col, label=n, edgecolor="white", linewidth=1)
        for i, v in enumerate(frac[c]):
            if v > 0.06:
                ax.text(left[i] + v / 2, i, f"{v:.0%}", ha="center", va="center", fontsize=8,
                        color="white" if col == BLUE else INK)
        left += frac[c].values
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.grid(axis="y", visible=False)
    ax.invert_yaxis()
    ax.legend(ncol=3, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.25))
    save(fig, out, "fig5_genbank_hypothetical_vs_phrogs")


def fig6_year(extra, out):
    t = pd.read_csv(extra / "by_record_year.csv")
    order = ["≤2010", "2011–2015", "2016–2020", "2021–2023", "2024–2026"]
    t = t.set_index("year_bin").reindex([o for o in order if o in set(t.year_bin)])
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    x = np.arange(len(t))
    for col, name, colr in (("genbank_vs_pyrodigal_gv", "GenBank vs Pyrodigal-gv", BLUE),
                            ("genbank_vs_phanotate", "GenBank vs PHANOTATE", ORANGE)):
        ax.plot(x, t[col], marker="o", color=colr, lw=2, label=name)
    ax.set_xticks(x, [f"{i}\n(n={int(n)})" for i, n in zip(t.index, t.genomes)], fontsize=7)
    ax.set_ylabel("Median per-genome stop Jaccard")
    ax.set_ylim(0.7, 1.0)
    ax.legend(fontsize=8)
    save(fig, out, "fig6_agreement_by_record_year")


def fig7_hosts(extra, out):
    t = pd.read_csv(extra / "by_host_genus_top15.csv").sort_values("no_specific")
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    y = np.arange(len(t))
    ax.scatter(t.no_specific, y, color=BLUE, s=30, label="No specific function", zorder=3)
    ax.scatter(t.single_support, y, color=ORANGE, s=30, label="Single-caller loci", zorder=3)
    ax.set_yticks(y, [f"{h} (n={n})" for h, n in zip(t.metadata_host, t.genomes)], fontsize=7)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.set_xlim(0, 1)
    ax.set_xlabel("Median share per genome")
    ax.legend(fontsize=7, loc="lower right")
    save(fig, out, "fig7_by_host_genus")


def fig8_density(extra, out):
    sc = pd.read_csv(extra / "start_codon_usage.csv")
    sc = sc[sc.start_codon.isin(["ATG", "GTG", "TTG"])]
    piv = sc.pivot_table(index="source", columns="start_codon", values="share").fillna(0)
    d = pd.read_csv(extra / "gene_density_by_source.csv", index_col=0)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ax = axes[0]
    srcs = [s for s in ("genbank", "phanotate", "pyrodigal_gv") if s in d.index]
    ax.bar([LABEL[s] for s in srcs], d.loc[srcs, "cds_per_kb_median"],
           yerr=[d.loc[srcs, "cds_per_kb_median"] - d.loc[srcs, "cds_per_kb_q1"],
                 d.loc[srcs, "cds_per_kb_q3"] - d.loc[srcs, "cds_per_kb_median"]],
           color=[GREY, ORANGE, BLUE], capsize=3, ecolor=INK2)
    ax.set_ylabel("CDS per kb (median, IQR)")
    ax.tick_params(axis="x", labelsize=8)
    ax = axes[1]
    x = np.arange(len(srcs))
    for i, (codon, colr) in enumerate((("GTG", AQUA), ("TTG", YELLOW))):
        ax.bar(x + (i - 0.5) * 0.38, piv.loc[srcs, codon], 0.36, color=colr, label=codon)
    ax.set_xticks(x, [LABEL[s] for s in srcs], fontsize=8)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.set_ylabel("Share of start codons (ATG = rest)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, out, "fig8_density_start_codons")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    extra, out = Path(a.run) / "analysis" / "extra", Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for f in (fig5_genbank_products, fig6_year, fig7_hosts, fig8_density):
        f(extra, out)
    print(sorted(x.name for x in out.glob("fig[5-8]*.png")))


if __name__ == "__main__":
    main()
