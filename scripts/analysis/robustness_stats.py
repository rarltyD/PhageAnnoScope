"""Robustness of H1 (international-level checks; descriptive + GEE, no new compute).

    python scripts/analysis/robustness_stats.py --run results/full_cohort

(a) H1 OR under the five prespecified PHROGs threshold variants;
(b) H1 OR within protein-length classes (addresses confounding by length);
(c) H1 OR within the largest host genera (forest plot; heterogeneity).
GEE: logistic, independence, cluster-robust SE by relatedness block; within
strata the length adjustment is log-linear.  Writes <run>/analysis/robustness/*.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from make_figures import BLUE, INK2, ORANGE, save  # noqa: E402

from lowmem_gee import gee_logit_lowmem as gee_logit  # noqa: E402

VARIANTS = {"primary": "E≤1e-5, cov≥0.5", "e1e-3": "E≤1e-3", "e1e-10": "E≤1e-10",
            "cov0.3": "cov≥0.3", "cov0.7": "cov≥0.7"}
FORMULA = ("no_specific ~ support2 + bs(np.log(length_aa), df=4) + gc + np.log(length_bp)"
           " + overlap_any_ge1")  # same length adjustment as the primary H1 model


def odds(df):
    try:
        t = gee_logit(df, FORMULA)["terms"]["support2"]
        return {"n": int(len(df)), "or": t["odds_ratio"], "lo": t["or_ci95"][0],
                "hi": t["or_ci95"][1], "p": t["p"]}
    except Exception as error:  # report, do not hide
        return {"n": int(len(df)), "error": f"{type(error).__name__}: {error}"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--manifest", default="manifests/cohort_manifest_validated.csv.gz")
    p.add_argument("--fig-out", default=None)
    a = p.parse_args()
    run = Path(a.run)
    out = run / "analysis" / "robustness"
    out.mkdir(parents=True, exist_ok=True)
    base = pd.read_parquet(run / "analysis" / "locus_table.parquet",
                           columns=["locus_id", "protein_sha256", "support", "length_aa", "gc",
                                    "length_bp", "overlap_any_ge1", "block", "accession_version"])
    base = base[base.length_aa > 0].assign(support2=lambda d: (d.support == 2).astype(int))
    res = {"threshold_variants": {}, "length_classes": {}, "host_genera": {}}
    for v, label in VARIANTS.items():
        fn = pd.read_parquet(run / "phrogs" / f"protein_function_{v}.parquet",
                             columns=["protein_sha256", "function_status"])
        d = base.merge(fn, on="protein_sha256", how="inner")
        d = d[d.function_status != "not_searched_or_failed"]
        d = d.assign(no_specific=(d.function_status != "specific_function").astype(int))
        res["threshold_variants"][label] = {
            **odds(d), "prop_support1": float(d[d.support == 1].no_specific.mean()),
            "prop_support2": float(d[d.support == 2].no_specific.mean())}
        if v == "primary":
            primary = d
    bins = [(0, 50, "<50 aa"), (50, 100, "50–99 aa"), (100, 200, "100–199 aa"),
            (200, 10**6, "≥200 aa")]
    for lo, hi, label in bins:
        res["length_classes"][label] = odds(primary[(primary.length_aa >= lo)
                                                    & (primary.length_aa < hi)])
    meta = pd.read_csv(a.manifest, dtype=str, keep_default_na=False,
                       usecols=["accession_version", "metadata_host"])
    primary = primary.merge(meta, on="accession_version", how="left")
    top = primary.groupby("metadata_host").accession_version.nunique().sort_values(
        ascending=False)
    top = [h for h in top.index if h not in ("", "Unspecified")][:10]
    for h in top:
        res["host_genera"][h] = odds(primary[primary.metadata_host == h])
        res["host_genera"][h]["genomes"] = int(
            primary[primary.metadata_host == h].accession_version.nunique())
    (out / "robustness.json").write_text(json.dumps(res, indent=2, ensure_ascii=False,
                                                    default=float))
    rows = [(sec, k, v.get("or"), v.get("lo"), v.get("hi"), v.get("n"))
            for sec, d in res.items() for k, v in d.items()]
    pd.DataFrame(rows, columns=["section", "stratum", "or", "lo", "hi", "n"]).to_csv(
        out / "robustness_table.csv", index=False)
    # Forest plot (all three sections).
    fig_rows = [(f"{k}", v) for sec in ("threshold_variants", "length_classes", "host_genera")
                for k, v in res[sec].items() if "or" in v]
    heads = {"threshold_variants": "PHROGs thresholds", "length_classes": "Protein length",
             "host_genera": "Host genus"}
    fig, ax = plt.subplots(figsize=(7.4, 0.36 * (len(fig_rows) + 4) + 0.9))
    y, ticks, labels = 0, [], []
    for sec in ("threshold_variants", "length_classes", "host_genera"):
        labels.append(f"$\\bf{{{heads[sec].replace(' ', '~')}}}$")
        ticks.append(y)
        y -= 1
        for k, v in res[sec].items():
            if "or" not in v:
                continue
            col = ORANGE if sec == "threshold_variants" and k.startswith("E≤1e-5") else BLUE
            ax.errorbar(v["or"], y, xerr=[[v["or"] - v["lo"]], [v["hi"] - v["or"]]], fmt="o",
                        color=col, ms=4, capsize=2, lw=1.5)
            ticks.append(y)
            labels.append(f"{k} (n={v['n']:,})" if sec != "host_genera"
                          else f"{k} ({v['genomes']} genomes)")
            y -= 1
        y -= 0.5
    ax.axvline(1, color=INK2, lw=1, ls="--")
    ax.set_yticks(ticks, labels, fontsize=10)
    ax.set_xscale("log")
    ax.set_xticks([0.2, 0.3, 0.5, 0.7, 1, 1.5], ["0.2", "0.3", "0.5", "0.7", "1", "1.5"])
    ax.set_xlabel("Adjusted OR: no specific function, 2 vs 1 caller (95% CI)", fontsize=11)
    ax.grid(axis="y", visible=False)
    fig_out = Path(a.fig_out) if a.fig_out else out
    fig_out.mkdir(parents=True, exist_ok=True)
    save(fig, fig_out, "fig9_h1_robustness_forest")
    print(pd.read_csv(out / "robustness_table.csv").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
