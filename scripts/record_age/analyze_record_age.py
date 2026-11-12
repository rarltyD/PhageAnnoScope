"""Are GenBank annotations that disagree with both gene callers associated with record age?

Unit: genome (primary cohort, n = 17 862).  Outcome: ``discordant`` = GenBank CDS stop-set
Jaccard with BOTH callers below the threshold (primary 0.8; sensitivity 0.7 and 0.9).
Exposure: submission period (earliest Direct Submission date), not the LOCUS date, which is the
last modification.  Dependence: logistic GEE (independence, cluster-robust SE) and cluster
bootstrap over the same relatedness blocks as the primary analysis (TaxMyPhage genus, else NCBI
genus, else the genome alone).  Disagreement is a flag for review, not proof of error.
Exploratory analysis added after the primary analysis (see docs/PROTOCOL.md).

    python scripts/record_age/analyze_record_age.py results/descriptive/per_genome_summary.csv \
        results/record_age/record_metadata.csv.gz manifests/cohort_manifest_validated.csv.gz results/record_age
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

per_genome, meta_csv, manifest, out = sys.argv[1:5]
out = Path(out); out.mkdir(parents=True, exist_ok=True)
SEED, REPS = 30102026, 2000
PERIODS = [(1900, 2007, "≤2007"), (2008, 2011, "2008–2011"), (2012, 2015, "2012–2015"),
           (2016, 2019, "2016–2019"), (2020, 2023, "2020–2023"), (2024, 2100, "2024–2026")]
REF_PERIOD = "2020–2023"
PLACEHOLDER_TAXA = {"", "new_genus", "new_species", "unclassified", "unknown", "na", "n/a", "nan"}


def relatedness_block(row):
    # Same rule as phageannoscope.analysis.relatedness_block.
    for field in ("taxmyphage_genus", "ncbi_genus_metadata"):
        value = str(row.get(field) or "").strip()
        if value.casefold() not in PLACEHOLDER_TAXA and not value.endswith("new_name"):
            return f"{field}:{value}"
    return "genome:" + str(row["accession_version"])


def period(y):
    for lo, hi, lab in PERIODS:
        if lo <= y <= hi:
            return lab


g = pd.read_csv(per_genome)
m = pd.read_csv(meta_csv, parse_dates=["submitted", "last_modified"])
man = pd.read_csv(manifest, usecols=["accession_version", "taxmyphage_genus", "ncbi_genus_metadata"],
                  low_memory=False).drop_duplicates("accession_version")
d = g.drop(columns=["taxmyphage_genus"]).merge(m, on="accession_version", how="left") \
     .merge(man, on="accession_version", how="left")
d["block"] = d.apply(relatedness_block, axis=1)
d["gb_best"] = d[["genbank_vs_phanotate", "genbank_vs_pyrodigal_gv"]].max(axis=1)
d["sub_year"] = d.submitted.dt.year
d["mod_year"] = d.last_modified.dt.year
d["period"] = d.sub_year.map(lambda y: period(y) if pd.notna(y) else None)
d["log10_kb"] = np.log10(d.length_bp / 1000)
d["gc_c"] = (d.gc - d.gc.mean()) * 100
d["lifestyle"] = d.phageleads_lifestyle.fillna("unknown")
d["gb_cds_per_kb"] = d.gb_cds / (d.length_bp / 1000)
caller_n = d.gene_caller.value_counts()
d["caller"] = d.gene_caller.where(d.gene_caller.map(caller_n) >= 50, "Other stated")
for t in (0.7, 0.8, 0.9):
    d[f"disc_{int(t * 100)}"] = (d.gb_best < t).astype(int)
d["discordant"] = d.disc_80
# Two components of discordance: no CDS deposited at all, and annotated but discordant.
d["unannotated"] = (d.gb_cds.fillna(0) == 0).astype(int)
d["disc_annotated"] = ((d.gb_cds > 0) & (d.gb_best < 0.8)).astype(int)

excluded = {"no_genbank_jaccard": int(d.gb_best.isna().sum()), "no_submission_date": int(d.period.isna().sum())}
a = d.dropna(subset=["gb_best", "period"]).copy()
summary = {"components": {"unannotated": int(a.unannotated.sum()),
                          "discordant_annotated": int(a.disc_annotated.sum())},
           "n_genomes_input": len(d), "n_analysed": len(a), "n_blocks": a.block.nunique(),
           "excluded": excluded, "threshold_primary": 0.8,
           "overall_discordant": {"n": int(a.discordant.sum()), "share": round(a.discordant.mean(), 4)},
           "submission_vs_modification_year": {
               "same_year_share": round((a.sub_year == a.mod_year).mean(), 4),
               "median_gap_years": float((a.mod_year - a.sub_year).median()),
               "spearman": round(a[["sub_year", "mod_year"]].corr("spearman").iloc[0, 1], 4)}}


def cluster_boot_share(df, col, by, reps=REPS, seed=SEED):
    """Cluster bootstrap (resample blocks) CI of the share of `col` within each level of `by`."""
    rng = np.random.default_rng(seed)
    agg = df.groupby(["block", by], observed=True)[col].agg(["sum", "count"]).unstack(by, fill_value=0)
    s, n = agg["sum"].to_numpy(float), agg["count"].to_numpy(float)
    idx = rng.integers(0, len(agg), size=(reps, len(agg)))
    boot = s[idx].sum(1) / np.where(n[idx].sum(1) == 0, np.nan, n[idx].sum(1))
    lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    return pd.DataFrame({"ci_low": lo, "ci_high": hi}, index=agg["sum"].columns)


order = [p[2] for p in PERIODS]
tab = a.groupby("period").agg(genomes=("discordant", "size"), discordant=("discordant", "sum"),
                              share=("discordant", "mean"), unannotated=("unannotated", "mean"),
                              disc_annotated=("disc_annotated", "mean"), median_gb_best=("gb_best", "median"),
                              median_gb_cds_per_kb=("gb_cds_per_kb", "median"),
                              share_disc_70=("disc_70", "mean"), share_disc_90=("disc_90", "mean"))
tab["median_caller_loci_per_kb"] = a.assign(x=a.loci / (a.length_bp / 1000)).groupby("period").x.median()
tab["disc_among_annotated"] = a[a.unannotated == 0].groupby("period").disc_annotated.mean()
tab = tab.join(cluster_boot_share(a, "discordant", "period"))
for col in ("unannotated", "disc_annotated"):
    tab = tab.join(cluster_boot_share(a, col, "period").add_prefix(col + "_"))
tab = tab.reindex(order)
tab.to_csv(out / "by_submission_period.csv")

year_tab = a.groupby("sub_year").agg(genomes=("discordant", "size"), share=("discordant", "mean"))
year_tab.to_csv(out / "by_submission_year.csv")

cat_tabs = {}
for col in ("seq_tech", "caller", "lifestyle"):
    t = a.groupby(col).agg(genomes=("discordant", "size"), share=("discordant", "mean"),
                           median_sub_year=("sub_year", "median"))
    t = t.join(cluster_boot_share(a, "discordant", col))
    t.to_csv(out / f"by_{col}.csv")
    cat_tabs[col] = t.round(4).to_dict(orient="index")
summary["by_category"] = cat_tabs


def gee(formula, df, label):
    model = smf.gee(formula, groups=df.block, data=df, family=sm.families.Binomial(),
                    cov_struct=sm.cov_struct.Independence())
    r = model.fit()
    ci = r.conf_int()
    rows = pd.DataFrame({"OR": np.exp(r.params), "ci_low": np.exp(ci[0]), "ci_high": np.exp(ci[1]),
                         "p": r.pvalues}).drop(index="Intercept")
    rows.insert(0, "model", label)
    return rows, {"formula": formula, "n_obs": int(r.nobs), "n_clusters": int(df.block.nunique()),
                  "n_events": int(df[formula.split(" ~")[0]].sum())}


P = f"C(period, Treatment('{REF_PERIOD}'))"
ADJ = "C(seq_tech, Treatment('Illumina/Ion/MGI')) + log10_kb + gc_c + C(lifestyle)"
models, fits = [], {}
specs = [
    ("M1_period", f"discordant ~ {P}", a),
    ("M2_period_adjusted", f"discordant ~ {P} + {ADJ}", a),
    ("M3_plus_caller", f"discordant ~ {P} + {ADJ} + C(caller, Treatment('Not stated'))", a),
    ("U_unannotated_adjusted", f"unannotated ~ {P} + {ADJ}", a),
    ("D_discordant_annotated_adjusted", f"disc_annotated ~ {P} + {ADJ}", a[a.unannotated == 0]),
    ("S70_period_adjusted", f"disc_70 ~ {P} + {ADJ}", a),
    ("S90_period_adjusted", f"disc_90 ~ {P} + {ADJ}", a),
    ("S_trend_per_5y_adjusted", f"discordant ~ I((sub_year - 2020) / 5) + {ADJ}", a),
]
for label, f, df in specs:
    rows, info = gee(f, df, label)
    models.append(rows); fits[label] = info
models = pd.concat(models)
models.to_csv(out / "gee_models.csv")
summary["gee"] = fits
summary["gee_period_or"] = {
    lab: models[models.model == lab].filter(like="period", axis=0)
         .rename(index=lambda s: s.split("T.")[-1].rstrip("]")).round(4).drop(columns="model")
         .to_dict(orient="index")
    for lab in ("M1_period", "M2_period_adjusted", "M3_plus_caller", "U_unannotated_adjusted",
                "D_discordant_annotated_adjusted", "S70_period_adjusted", "S90_period_adjusted")}
summary["gee_trend_per_5y"] = models[models.model == "S_trend_per_5y_adjusted"].filter(
    like="sub_year", axis=0).round(4).drop(columns="model").to_dict(orient="records")[0]

# Direction of discordance: is GenBank under- or over-calling relative to the callers?
a["gb_minus_callers_per_kb"] = a.gb_cds_per_kb - a.loci / (a.length_bp / 1000)
summary["direction"] = a.groupby("discordant").gb_minus_callers_per_kb.describe()[
    ["count", "25%", "50%", "75%"]].round(3).to_dict(orient="index")

# Discordant records (oldest first), for manual review.
cols = ["unannotated", "disc_annotated", "accession_version", "submitted", "last_modified", "gb_best",
        "genbank_vs_phanotate", "genbank_vs_pyrodigal_gv", "gb_cds", "loci", "length_bp", "seq_tech", "caller",
        "block"]
a[a.discordant == 1].sort_values("submitted")[cols].to_csv(out / "discordant_records.csv.gz", index=False)
a[cols + ["discordant", "sub_year", "period"]].to_csv(out / "genome_table.csv.gz", index=False)

# Figure: share by submission year (genome counts as bars) with period CIs.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True, height_ratios=[3, 1.3])
yt = year_tab[year_tab.index >= 1995]
ax1.plot(yt.index, yt.share * 100, "o", color="#4C72B0", ms=3.5, alpha=.6, label="per submission year")
mid = {"≤2007": 2003, "2008–2011": 2009.5, "2012–2015": 2013.5, "2016–2019": 2017.5,
       "2020–2023": 2021.5, "2024–2026": 2025}
for p, r in tab.iterrows():
    x = mid[p]
    ax1.errorbar(x, r.share * 100, yerr=[[(r.share - r.ci_low) * 100], [(r.ci_high - r.share) * 100]],
                 fmt="s", color="#C44E52", ms=6, capsize=3)
ax1.plot([], [], "s", color="#C44E52", label="period (cluster bootstrap 95% CI)")
ax1.set_ylabel("Discordant GenBank\nannotations, %")
ax1.legend(frameon=False, fontsize=8)
ax2.bar(yt.index, yt.genomes, color="#8C8C8C", width=.8)
ax2.set_ylabel("Genomes")
ax2.set_xlabel("Year of GenBank submission")
for ax in (ax1, ax2):
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(out / "fig_discordance_by_submission_year.png", dpi=300)
fig.savefig(out / "fig_discordance_by_submission_year.pdf")

summary["by_submission_period"] = tab.round(4).to_dict(orient="index")
(out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str))
print(json.dumps({k: summary[k] for k in ("n_analysed", "n_blocks", "excluded", "overall_discordant",
                                          "submission_vs_modification_year")}, indent=1, default=str))
print(tab.round(3).to_string())
print(models.round(3).to_string())
