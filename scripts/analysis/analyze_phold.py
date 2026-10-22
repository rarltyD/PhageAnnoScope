"""H3: Phold functional gain on the stratified probability sample (protocol section 5/6).

    python scripts/analysis/analyze_phold.py --run results/<run> --phold results/<run>/phold

gain      = Phold function != 'unknown function' (default Phold Foldseek E <= 1e-3)
            among sampled proteins WITHOUT specific PHROGs-HMM function.
confirmed = gain AND qCov >= 0.7 AND tCov >= 0.7.
controls  = sampled proteins WITH specific PHROGs-HMM function: category agreement.
Per-stratum Wilson 95% CI (per-genome cap 2 limits clustering; noted).  The
design-weighted estimate is for the first nested sample (981 genomes) sampling frame only.
Phold and PHROGs share reference information: H3 is not independent validation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [c - h, c + h]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--phold", required=True)
    args = parser.parse_args()
    run, ph = Path(args.run), Path(args.phold)
    sample = pd.read_csv(run / "phold_sample" / "phold_sample.csv")
    pred = pd.read_csv(ph / "compare" / "phold_per_cds_predictions.tsv", sep="\t")
    pred = pred.rename(columns={"cds_id": "protein_sha256", "function": "phold_function",
                                "phrog": "phold_phrog", "product": "phold_product",
                                "evalue": "phold_evalue"})
    df = sample.merge(pred, on="protein_sha256", how="left")
    df["phold_processed"] = df["phold_function"].notna()
    df["phold_informative"] = df["phold_processed"] & (df["phold_function"] != "unknown function")
    df["phold_confirmed"] = df["phold_informative"] & (df["qCov"] >= 0.7) & (df["tCov"] >= 0.7)
    df.to_csv(ph / "h3_protein_table.csv", index=False)
    out = {"n_sampled": int(len(df)), "n_processed": int(df.phold_processed.sum()),
           "technical_failures": int((~df.phold_processed).sum()), "strata": {}}
    for stratum, sub in df.groupby("stratum"):
        sub = sub[sub.phold_processed]
        k, kc, n = int(sub.phold_informative.sum()), int(sub.phold_confirmed.sum()), len(sub)
        row = {"n": n, "eligible_in_frame": int(sub.eligible_in_stratum.iloc[0]),
               "informative": k, "informative_prop": k / n if n else None,
               "informative_ci95": wilson(k, n), "confirmed": kc,
               "confirmed_prop": kc / n if n else None, "confirmed_ci95": wilson(kc, n)}
        if sub.group.iloc[0] == "specific":
            agree = (sub.phold_function == sub.best_category).sum()
            row.update(control_category_agreement=int(agree),
                       control_agreement_prop=agree / n if n else None,
                       control_agreement_ci95=wilson(int(agree), n))
        out["strata"][stratum] = row
    gain = df[(df.group == "no_specific") & df.phold_processed]
    out["gain_unweighted_all_no_specific"] = {
        "n": int(len(gain)), "informative": int(gain.phold_informative.sum()),
        "prop": float(gain.phold_informative.mean()),
        "ci95": wilson(int(gain.phold_informative.sum()), len(gain)),
        "confirmed_prop": float(gain.phold_confirmed.mean())}
    # Design-weighted (Hajek) estimate for the sampling frame of proteins without a specific function.
    w = 1 / gain.inclusion_probability
    out["gain_weighted_sample_frame"] = {
        "informative_prop": float((w * gain.phold_informative).sum() / w.sum()),
        "confirmed_prop": float((w * gain.phold_confirmed).sum() / w.sum()),
        "note": "Hajek estimator, weights = 1/inclusion probability per stratum; frame = "
                "unique representative proteins of primary loci of the first nested sample (981 genomes) without specific "
                "PHROGs function; not generalized beyond this frame; no CI (approximate "
                "inclusion probabilities due to per-genome cap)."}
    out["gain_by_function_status"] = {
        s: {"n": int(len(g)), "informative_prop": float(g.phold_informative.mean())}
        for s, g in gain.groupby("function_status")}
    out["top_new_categories"] = gain[gain.phold_informative].phold_function.value_counts() \
        .head(10).to_dict()
    out["caveat"] = ("Phold reference DB includes PHROG-derived annotations; agreement and gain "
                     "are not independent validation. Functional assignment != experimental proof.")
    (ph / "h3_summary.json").write_text(json.dumps(out, indent=2, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
