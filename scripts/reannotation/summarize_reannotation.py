"""Summary tables and figure for the consensus re-annotation.

(a) record classes by submission period; (b) PHROGs categories of proposed additions (add_A) in
annotated records; (c) protein length of GenBank CDS without caller or PHROGs support versus
confirmed GenBank CDS.  Also a sensitivity analysis of the 10% threshold for 'major'.

    python scripts/reannotation/summarize_reannotation.py results/reannotation results/full_cohort
"""
import json
import sys
from pathlib import Path

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

res, run = Path(sys.argv[1]), Path(sys.argv[2])
ORDER = ["≤2007", "2008–2011", "2012–2015", "2016–2019", "2020–2023", "2024–2026"]
CLASSES = ["consistent", "minor", "major", "unannotated"]
COLORS = {"consistent": "#8FB996", "minor": "#E9C46A", "major": "#C2847A", "unannotated": "#9DB4C0"}
r = pd.read_csv(res / "records_reannotation.csv.gz")
g = pd.read_csv(res / "gene_proposals.csv.gz")
ann = r[r.unannotated == 0]
out = {}

# Gene-level table (Table 4).
summ = json.loads((res / "reannotation_summary.json").read_text())
gc = summ["gene_categories"]
annotated_ids = set(ann.accession_version)
addA_ann = g[(g.category_rule == "add_A") & g.accession_version.isin(annotated_ids)]
out["table4"] = {
    "gb_confirmed": gc["confirmed"], "gb_supported_one_caller": gc["supported"],
    "gb_only_function": gc["gb_only_function"], "gb_only_unsupported": gc["gb_only_unsupported"],
    "gb_only_unsupported_records": int((ann.gb_only_unsupported > 0).sum()),
    "add_A_annotated": int(len(addA_ann)), "add_A_annotated_records": int(addA_ann.accession_version.nunique()),
    "add_A_unannotated": int(gc["add_A"] - len(addA_ann)),
    "add_B": gc["add_B"], "start_corrections": summ["start_corrections"],
    "start_correction_records": int((r.start_corrections > 0).sum()),
    "consensus_genes": summ["consensus_genes"],
    "functional_relabels": summ["functional_relabels"],
    "functional_relabel_records": summ["functional_relabel_records"],
}
out["functional_relabel_by_category"] = summ["functional_relabel_by_category"]
gu = g[g.category_rule == "gb_only_unsupported"]
out["gb_only_unsupported_length"] = gu.length_aa.describe(percentiles=[.25, .5, .75]).round(1).to_dict()
out["gb_only_unsupported_lt50_share"] = round(float((gu.length_aa < 50).mean()), 4)
out["flagged_per_annotated_record"] = ann.flagged.describe(percentiles=[.5, .75, .9, .99]).round(2).to_dict()
out["minor_with_one_flagged_gene"] = int(((ann.reannotation_class == "minor") & (ann.flagged == 1)).sum())
out["major_from_part_a"] = int(((ann.reannotation_class == "major") & (ann.disc_annotated == 1)).sum())
out["record_classes"] = r.reannotation_class.value_counts().reindex(CLASSES).to_dict()
cp = pd.crosstab(r.period, r.reannotation_class).reindex(index=ORDER, columns=CLASSES).fillna(0)
out["class_by_period_n"] = cp.astype(int).to_dict(orient="index")
out["class_by_caller_share"] = pd.crosstab(r.caller, r.reannotation_class, normalize="index").round(4).to_dict(orient="index")

# Sensitivity of the 'major' threshold.
sens = {}
for t in (0.05, 0.10, 0.20):
    major = (ann.flagged_share >= t) | (ann.disc_annotated == 1)
    sens[str(t)] = {"major": int(major.sum()), "minor": int((~major & (ann.flagged > 0)).sum())}
out["major_threshold_sensitivity"] = sens

# Confirmed GenBank CDS lengths for the comparison panel.
con = duckdb.connect()
con.execute("SET memory_limit='2GB'; SET threads=2")
conf_len = con.execute(f"""
    SELECT feature_nt_length // 3 - 1 AS aa FROM read_parquet('{(run / 'parquet' / 'cds_calls' / '*.parquet').as_posix()}',
    union_by_name=true) WHERE source = 'genbank' AND eligible_primary USING SAMPLE 200000 ROWS""").df().aa
out["genbank_all_length_median_sample"] = float(conf_len.median())

fig = plt.figure(figsize=(7.4, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
ax = fig.add_subplot(gs[0, :])
share = cp.div(cp.sum(axis=1), axis=0) * 100
bottom = np.zeros(len(share))
for c in CLASSES:
    ax.bar(range(len(ORDER)), share[c], bottom=bottom, color=COLORS[c], width=.65, label=c)
    bottom += share[c].to_numpy()
ax.set_xticks(range(len(ORDER)), [f"{p}\nn={int(n):,}" for p, n in zip(ORDER, cp.sum(axis=1))], fontsize=7.5)
ax.set_ylabel("Records, %", fontsize=9)
ax.set_xlabel("Year of GenBank submission", fontsize=9)
ax.legend(frameon=False, fontsize=7.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.13))
ax.set_title("a", loc="left", fontsize=12, fontweight="bold", x=-0.07)

ax = fig.add_subplot(gs[1, 0])
cat = addA_ann.category.fillna("unknown").value_counts().head(8)[::-1]
ax.barh(range(len(cat)), cat.values, color="#6C8EBF")
ax.set_yticks(range(len(cat)), [c if len(c) < 32 else c[:30] + "…" for c in cat.index], fontsize=7)
ax.set_xlabel("Proposed additions in annotated records", fontsize=8.5)
ax.set_title("b", loc="left", fontsize=12, fontweight="bold", x=-0.75)

ax = fig.add_subplot(gs[1, 1])
bins = np.logspace(np.log10(10), np.log10(2000), 40)
ax.hist(conf_len.clip(10, 2000), bins=bins, density=True, histtype="step", color="#555555", lw=1.2,
        label="All GenBank CDS (sample)")
ax.hist(gu.length_aa.clip(10, 2000), bins=bins, density=True, histtype="stepfilled", alpha=.6,
        color="#C2847A", label="GenBank CDS without\ncaller or PHROGs support")
ax.set_xscale("log")
ax.set_xlabel("Protein length, aa", fontsize=8.5)
ax.set_ylabel("Density", fontsize=8.5)
ax.legend(frameon=False, fontsize=6.5)
ax.set_title("c", loc="left", fontsize=12, fontweight="bold", x=-0.2)
for a in fig.axes:
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=7.5)
fig.tight_layout()
fig.savefig(res / "reannotation.png", dpi=300)
fig.savefig(res / "reannotation.pdf")
(res / "summary.json").write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str))
print(json.dumps({k: out[k] for k in ("table4", "record_classes", "major_threshold_sensitivity",
                                       "gb_only_unsupported_length", "gb_only_unsupported_lt50_share",
                                       "minor_with_one_flagged_gene", "major_from_part_a",
                                       "genbank_all_length_median_sample")}, indent=1, default=str))
print(cp.astype(int).to_string())
