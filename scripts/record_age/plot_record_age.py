"""Figure: the two components of GenBank discordance by submission period.

(a) Share of genomes deposited without any CDS; (b) share of annotated genomes whose CDS set
disagrees with both callers (stop-set Jaccard < 0.8).  Error bars: 95% cluster bootstrap over
relatedness blocks (2 000 replicates, seed 30102026).  Also writes the CIs to a CSV.

    python scripts/record_age/plot_record_age.py results/record_age/genome_table.csv.gz results/record_age
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

genome_table, out_dir = sys.argv[1], Path(sys.argv[2])
ORDER = ["≤2007", "2008–2011", "2012–2015", "2016–2019", "2020–2023", "2024–2026"]
d = pd.read_csv(genome_table)
d["unannotated"] = (d.gb_cds.fillna(0) == 0).astype(int)
d["disc_annotated"] = ((d.gb_cds > 0) & (d.gb_best < 0.8)).astype(int)


def boot(df, col, reps=2000, seed=30102026):
    rng = np.random.default_rng(seed)
    agg = df.groupby(["block", "period"])[col].agg(["sum", "count"]).unstack("period", fill_value=0)
    s, n = agg["sum"].reindex(columns=ORDER).to_numpy(float), agg["count"].reindex(columns=ORDER).to_numpy(float)
    idx = rng.integers(0, len(agg), size=(reps, len(agg)))
    b = s[idx].sum(1) / n[idx].sum(1)
    est = s.sum(0) / n.sum(0)
    lo, hi = np.nanpercentile(b, [2.5, 97.5], axis=0)
    return pd.DataFrame({"n": n.sum(0).astype(int), "events": s.sum(0).astype(int),
                         "share": est, "ci_low": lo, "ci_high": hi}, index=ORDER)


panels = [("unannotated", d, "Genomes without any GenBank CDS, %", "a"),
          ("disc_annotated", d[d.unannotated == 0], "Annotated genomes discordant\nwith both callers, %", "b")]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4))
tabs = []
for (col, df, ylab, letter), ax in zip(panels, axes):
    t = boot(df, col)
    tabs.append(t.add_prefix(col + "_"))
    x = np.arange(len(ORDER))
    ax.bar(x, t.share * 100, color="#9DB4C0" if col == "unannotated" else "#C2847A", width=.65)
    ax.errorbar(x, t.share * 100, yerr=[(t.share - t.ci_low) * 100, (t.ci_high - t.share) * 100],
                fmt="none", ecolor="black", capsize=3, lw=1)
    ax.set_xticks(x, [f"{p[:5] + p[7:] if '–' in p else p}\nn={n:,}" for p, n in zip(ORDER, t.n)],
                  fontsize=7)
    ax.set_ylabel(ylab, fontsize=9)
    ax.set_xlabel("Year of GenBank submission", fontsize=9)
    ax.set_title(letter, loc="left", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(out_dir / "record_age_components.png", dpi=300)
fig.savefig(out_dir / "record_age_components.pdf")
pd.concat(tabs, axis=1).to_csv(out_dir / "components_by_period.csv")
print(pd.concat(tabs, axis=1).round(4).to_string())
