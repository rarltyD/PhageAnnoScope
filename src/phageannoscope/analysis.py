"""Locus-level analysis table and H1/H2 statistics (protocol section 6).

Unit: primary-eligible stop-locus in the union of successfully completed
callers (PHANOTATE, Pyrodigal-gv).  ``support`` is the number of eligible
caller families supporting the locus (1 or 2).  Function status is taken for
the prespecified representative model's protein.  Dependence: all loci of all
genomes in one relatedness block (TaxMyPhage genus; otherwise NCBI genus
metadata; otherwise the genome alone) form one cluster for GEE robust SEs and
cluster bootstrap.  GenBank is a baseline, never a vote.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

PLACEHOLDER_TAXA = {"", "new_genus", "new_species", "unclassified", "unknown", "na", "n/a"}


def relatedness_block(row) -> str:
    for field in ("taxmyphage_genus", "ncbi_genus_metadata"):
        value = str(row.get(field) or "").strip()
        if value.casefold() not in PLACEHOLDER_TAXA and not value.endswith("new_name"):
            return f"{field}:{value}"
    return "genome:" + str(row["accession_version"])


def _parquet(run: Path, table: str) -> str:
    return (run / "parquet" / table / "*.parquet").as_posix()


def locus_table(run, manifest, function_parquet=None) -> pd.DataFrame:
    """One row per primary-eligible stop-locus with covariates."""
    run = Path(run)
    con = duckdb.connect()
    loci = con.execute(
        f"""
        SELECT l.locus_id, l.genome_id, l.representative_model_id, l.caller_support_eligible
               AS support, l.strand, l.genbank_present, l.genbank_exact_match,
               l.genbank_stop_match, l.additional_stop_model, l.phanotate_present,
               l.pyrodigal_gv_present, l.genbank_comparable, l.genbank_complex_context,
               l.genbank_partial_context, l.genbank_antisense_context,
               l.unannotated_same_strand_region, m.segments, m.feature_nt_length
        FROM read_parquet('{_parquet(run, "stop_loci")}', union_by_name=true) l
        JOIN (SELECT DISTINCT model_id, segments, feature_nt_length
              FROM read_parquet('{_parquet(run, "models")}', union_by_name=true)) m
          ON m.model_id = l.representative_model_id
        WHERE l.eligible_primary
        """
    ).df()
    proteins = con.execute(
        f"""
        SELECT model_id, min(protein_sha256) AS protein_sha256, count(DISTINCT protein_sha256) n
        FROM read_parquet('{_parquet(run, "protein_origins")}', union_by_name=true)
        GROUP BY model_id
        """
    ).df()
    genbank_n = con.execute(
        f"""SELECT genome_id, count(*) FILTER (WHERE source = 'genbank') AS genbank_cds_n
            FROM read_parquet('{_parquet(run, "cds_calls")}', union_by_name=true)
            GROUP BY genome_id"""
    ).df()
    genomes = con.execute(
        f"""SELECT genome_id, any_value(accession_version) AS accession_version,
                   any_value(length_bp) AS length_bp,
                   any_value(length(regexp_replace(upper(sequence), '[^GC]', '', 'g'))
                             / greatest(1, length(sequence))) AS gc
            FROM read_parquet('{_parquet(run, "genomes")}', union_by_name=true)
            GROUP BY genome_id"""
    ).df()
    genomes = genomes.merge(genbank_n, on="genome_id", how="left")
    genomes["genbank_cds_n"] = genomes["genbank_cds_n"].fillna(0).astype(int)
    genomes["genbank_baseline"] = genomes["genbank_cds_n"] > 0
    meta = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    meta = meta[meta["accession_version"].isin(genomes["accession_version"])
                & (meta["eligible_for_pilot"] == "true")]
    meta["block"] = meta.apply(relatedness_block, axis=1)
    genomes = genomes.merge(
        meta[["accession_version", "block", "baltimore_group", "a2_slippage_only",
              "taxmyphage_genus", "topology"]], on="accession_version", how="left")
    genomes = genomes[genomes["accession_version"].isin(meta["accession_version"])]
    df = loci.merge(proteins, left_on="representative_model_id", right_on="model_id",
                    how="left").merge(genomes, on="genome_id", how="left")
    df["length_aa"] = df["feature_nt_length"] // 3 - 1
    df["start0"] = df["segments"].map(lambda s: min(x["start0"] for x in s))
    df["end0"] = df["segments"].map(lambda s: max(x["end0"] for x in s))
    df = _overlaps(df)
    df["single_support"] = (df["support"] == 1).astype(int)
    df["short"] = (df["length_aa"] < 100).astype(int)
    if function_parquet:
        fn = pd.read_parquet(function_parquet)[
            ["protein_sha256", "function_status", "best_phrog", "best_category", "best_annot",
             "best_evalue"]]
        df = df.merge(fn, on="protein_sha256", how="left")
        df["function_status"] = df["function_status"].fillna("not_searched_or_failed")
    return df.drop(columns=["segments", "model_id"])


def _overlaps(df: pd.DataFrame) -> pd.DataFrame:
    """Representative-model overlaps between different loci (>=1 and >=30 nt)."""
    same1 = np.zeros(len(df), bool)
    same30 = np.zeros(len(df), bool)
    anti1 = np.zeros(len(df), bool)
    for _, idx in df.groupby("genome_id").groups.items():
        sub = df.loc[idx].sort_values("start0")
        rows = list(zip(sub.index, sub["start0"], sub["end0"], sub["strand"]))
        for i, (ia, sa, ea, stra) in enumerate(rows):
            for ib, sb, eb, strb in rows[i + 1 :]:
                if sb >= ea:
                    break
                ov = min(ea, eb) - sb
                if ov <= 0:
                    continue
                if stra == strb:
                    for k in (ia, ib):
                        same1[df.index.get_loc(k)] = True
                        if ov >= 30:
                            same30[df.index.get_loc(k)] = True
                else:
                    for k in (ia, ib):
                        anti1[df.index.get_loc(k)] = True
    df["overlap_same_ge1"] = same1.astype(int)
    df["overlap_same_ge30"] = same30.astype(int)
    df["overlap_anti_ge1"] = anti1.astype(int)
    df["overlap_any_ge1"] = (same1 | anti1).astype(int)
    return df


def cluster_bootstrap(df, value, by, cluster="block", reps=2000, seed=30102026):
    """Proportion of ``value`` within each level of ``by``; percentile CI over clusters."""
    rng = np.random.default_rng(seed)
    clusters = df[cluster].unique()
    grouped = {c: g for c, g in df.groupby(cluster)}
    levels = sorted(df[by].unique())
    point = {lv: df.loc[df[by] == lv, value].mean() for lv in levels}
    sums = {c: g.groupby(by)[value].agg(["sum", "count"]) for c, g in grouped.items()}
    draws = {lv: [] for lv in levels}
    diffs = []
    for _ in range(reps):
        pick = rng.choice(clusters, size=len(clusters), replace=True)
        total = pd.concat([sums[c] for c in pick]).groupby(level=0).sum()
        props = {lv: (total.loc[lv, "sum"] / total.loc[lv, "count"]) if lv in total.index
                 else np.nan for lv in levels}
        for lv in levels:
            draws[lv].append(props[lv])
        if len(levels) == 2:
            diffs.append(props[levels[1]] - props[levels[0]])
    out = {"levels": {}, "n_clusters": int(len(clusters)), "reps": reps}
    for lv in levels:
        arr = np.array(draws[lv], float)
        out["levels"][str(lv)] = {
            "n": int((df[by] == lv).sum()), "proportion": float(point[lv]),
            "ci95": [float(np.nanpercentile(arr, 2.5)), float(np.nanpercentile(arr, 97.5))]}
    if diffs:
        arr = np.array(diffs, float)
        out["difference_pp"] = {
            "estimate": float(100 * (point[levels[1]] - point[levels[0]])),
            "ci95": [float(100 * np.nanpercentile(arr, 2.5)),
                     float(100 * np.nanpercentile(arr, 97.5))],
            "contrast": f"{levels[1]} minus {levels[0]}"}
    return out


def gee_logit(df, formula, cluster="block"):
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    data = df.sort_values(cluster).reset_index(drop=True)
    model = smf.gee(formula, groups=data[cluster], data=data,
                    family=sm.families.Binomial(), cov_struct=sm.cov_struct.Independence())
    fit = model.fit(cov_type="robust")
    ci = fit.conf_int()
    terms = {}
    for name in fit.params.index:
        terms[name] = {
            "coef": float(fit.params[name]), "se_robust": float(fit.bse[name]),
            "odds_ratio": float(math.exp(fit.params[name])),
            "or_ci95": [float(math.exp(ci.loc[name, 0])), float(math.exp(ci.loc[name, 1]))],
            "p": float(fit.pvalues[name])}
    return {"formula": formula, "n_obs": int(fit.nobs), "n_clusters": int(data[cluster].nunique()),
            "cov_struct": "independence", "cov_type": "robust (cluster sandwich)",
            "terms": terms}


def holm(pvalues: dict) -> dict:
    items = sorted(pvalues.items(), key=lambda kv: kv[1])
    m, running, adjusted = len(items), 0.0, {}
    for i, (key, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        adjusted[key] = running
    return adjusted


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=float) + "\n", encoding="utf-8")
