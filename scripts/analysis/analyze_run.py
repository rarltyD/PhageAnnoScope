"""Descriptive tables + H1/H2 GEE + H4 candidate list for one pipeline run.

    python scripts/analysis/analyze_run.py --run results/<run> \
        --manifest manifests/cohort_manifest_validated.csv.gz

Outputs to <run>/analysis/.  Proportion CIs: cluster (relatedness-block)
bootstrap, 2000 reps, seed 30102026.  GEE: logistic, independence working
correlation, cluster-robust SE; Holm over the two primary tests (H1, H2).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

from fast_stats import cluster_bootstrap  # vectorised; same estimand
from lowmem_gee import gee_logit_lowmem as gee_logit  # identical estimates, low memory
from phageannoscope.analysis import holm, locus_table, write_json

COVARIATES = "bs(log_len, df=4) + gc + log_genome"
SLIM_COLUMNS = ["genome_id", "accession_version", "block", "support", "support2", "single_support",
                "short", "length_aa", "log_len", "gc", "length_bp", "log_genome",
                "overlap_any_ge1", "overlap_same_ge1", "overlap_anti_ge1", "overlap_same_ge30",
                "function_status", "a2_slippage_only", "genbank_baseline", "genbank_cds_n",
                "phanotate_present", "pyrodigal_gv_present", "genbank_stop_match",
                "additional_stop_model"]


def prepare(df):
    df = df.copy()
    df["log_len"] = (df["length_aa"].clip(lower=1)).map(float).pipe(lambda s: s.apply(
        lambda v: __import__("math").log(v)))
    df["log_genome"] = df["length_bp"].map(float).apply(lambda v: __import__("math").log(v))
    df["support2"] = (df["support"] == 2).astype(int)
    return df


def h1(df, reps):
    data = df[df["function_status"].isin(
        ["specific_function", "homology_unknown_function", "no_informative_hit", "ambiguous"])]
    data = data.assign(no_specific=(data["function_status"] != "specific_function").astype(int))
    out = {"n_loci": int(len(data)),
           "excluded_not_searched": int((df["function_status"] == "not_searched_or_failed").sum()),
           "status_by_support": data.groupby(["support", "function_status"]).size()
           .unstack(fill_value=0).to_dict(orient="index"),
           "bootstrap_no_specific_by_support": cluster_bootstrap(
               data, "no_specific", "support", reps=reps)}
    formula = f"no_specific ~ support2 + {COVARIATES} + overlap_any_ge1"
    out["gee"] = gee_logit(data, formula)
    sens = {}
    variants = {
        "exclude_slippage_genomes": data[data["a2_slippage_only"] != "true"],
        "length_ge_30aa": data[data["length_aa"] >= 30],
        "exclude_no_genbank_cds_genomes": data[data["genbank_baseline"]],
        "short_lt100aa_only": data[data["length_aa"] < 100],
        "long_ge100aa_only": data[data["length_aa"] >= 100],
        "one_genome_per_block": data[data["accession_version"].isin(
            data.sort_values("accession_version").groupby("block")["accession_version"].first())],
    }
    for name, sub in variants.items():
        try:
            sens[name] = {"n": int(len(sub)), "support2_or": gee_logit(sub, formula)
                          ["terms"]["support2"]}
        except Exception as error:  # report, never hide
            sens[name] = {"n": int(len(sub)), "error": f"{type(error).__name__}: {error}"}
    hit = data[data["function_status"] != "specific_function"].assign(
        any_homology=lambda d: (d["function_status"] != "no_informative_hit").astype(int))
    sens["unknown_homology_vs_nohit_by_support"] = cluster_bootstrap(
        hit, "any_homology", "support", reps=reps)
    out["sensitivity"] = sens
    return out


def h2(df, reps):
    data = df
    out = {"n_loci": int(len(data)),
           "bootstrap_single_support_by_short": cluster_bootstrap(
               data, "single_support", "short", reps=reps),
           "bootstrap_single_support_by_overlap_same": cluster_bootstrap(
               data, "single_support", "overlap_same_ge1", reps=reps)}
    formula = f"single_support ~ {COVARIATES} + overlap_same_ge1 + overlap_anti_ge1"
    out["gee"] = gee_logit(data, formula)
    linear = f"single_support ~ log_len + gc + log_genome + overlap_same_ge1 + overlap_anti_ge1"
    out["gee_linear_loglen"] = gee_logit(data, linear)
    sens = {}
    for name, sub in {
        "exclude_slippage_genomes": data[data["a2_slippage_only"] != "true"],
        "length_ge_30aa": data[data["length_aa"] >= 30],
        "overlap_ge30nt": data.assign(overlap_same_ge1=data["overlap_same_ge30"]),
    }.items():
        try:
            sens[name] = {"n": int(len(sub)), "terms": {
                k: v for k, v in gee_logit(sub, linear)["terms"].items()
                if k in ("log_len", "overlap_same_ge1", "overlap_anti_ge1")}}
        except Exception as error:
            sens[name] = {"n": int(len(sub)), "error": f"{type(error).__name__}: {error}"}
    out["sensitivity"] = sens
    return out


def descriptive(run, df):
    con = duckdb.connect()
    con.execute("SET memory_limit='2GB'; SET threads=3")
    p = lambda t: f"read_parquet('{(Path(run) / 'parquet' / t).as_posix()}/*.parquet', union_by_name=true)"  # noqa: E731
    metrics = con.execute(f"""
        SELECT source_a, source_b, count(*) genomes,
               sum(exact_intersection) exact_inter, sum(exact_union) exact_union,
               sum(stop_intersection) stop_inter, sum(stop_union) stop_union,
               sum(exact_a_count) a_n, sum(exact_b_count) b_n,
               sum(exact_cds_count) exact_cds, sum(alternative_start_count) alt_start,
               sum(contained_distinct_stop_count) contained, sum(partial_overlap_count) partial_overlap,
               sum(antisense_overlap_count) antisense,
               sum(no_correspondence_count) no_corr, sum(not_comparable_count) not_comparable
        FROM {p('genome_metrics')}
        WHERE metric_scope='source_pair'
        GROUP BY 1,2 ORDER BY 1,2""").df()
    if len(metrics):
        metrics["pooled_exact_jaccard"] = metrics.exact_inter / metrics.exact_union
        metrics["pooled_stop_jaccard"] = metrics.stop_inter / metrics.stop_union
        metrics["pooled_exact_dice"] = 2 * metrics.exact_inter / (metrics.a_n + metrics.b_n)
    calls = con.execute(f"""
        SELECT source, count(*) calls, count(DISTINCT genome_id) genomes
        FROM {p('cds_calls')} WHERE coalesce(rotation_shift,0)=0 GROUP BY 1""").df()
    loci = {
        "n_loci_primary": int(len(df)),
        "support_counts": df["support"].value_counts().to_dict(),
        "single_support_by_caller": {
            "phanotate_only": int(((df.support == 1) & df.phanotate_present).sum()),
            "pyrodigal_gv_only": int(((df.support == 1) & df.pyrodigal_gv_present).sum())},
        "median_length_aa_by_support": df.groupby("support")["length_aa"].median().to_dict(),
        "genbank_stop_match_by_support": df.groupby("support")["genbank_stop_match"].mean()
        .to_dict(),
        "additional_vs_genbank_by_support": df.groupby("support")["additional_stop_model"].sum()
        .to_dict(),
        "n_genomes": int(df["genome_id"].nunique()), "n_blocks": int(df["block"].nunique()),
    }
    return metrics, calls, loci


def h4(df, lt_path=None):
    """Split: (a) baseline genomes, clean context; (b) flagged context; (c) no GenBank CDS."""
    if lt_path is not None:  # low-memory mode: per-genome counts from slim df, rows filtered
        sup2_n = df[df["support"] == 2].groupby("genome_id", observed=True).size()
        gb_n = df.groupby("genome_id", observed=True)["genbank_cds_n"].first()
        df = pd.read_parquet(lt_path, filters=[("additional_stop_model", "==", True),
                                               ("function_status", "==", "specific_function")])
        df = df.assign(_sup2=df["genome_id"].map(sup2_n), _gbn=df["genome_id"].map(gb_n))
    base = (df["additional_stop_model"]) & (df["function_status"] == "specific_function")
    clean = ~(df["genbank_complex_context"].fillna(False) | df["genbank_partial_context"]
              .fillna(False))
    # Amendment A6: baseline completeness = GenBank CDS / loci supported by both callers.
    if "_sup2" in df:
        ratio = df["_gbn"] / df["_sup2"]
    else:
        sup2 = df[df["support"] == 2].groupby("genome_id").size()
        ratio = df["genome_id"].map(df.groupby("genome_id")["genbank_cds_n"].first() / sup2)
    df = df.assign(genbank_completeness_ratio=ratio)
    dense = ratio >= 0.8
    layers = {
        "baseline_clean": df[base & df["genbank_baseline"] & clean & dense],
        "baseline_sparse_annotation": df[base & df["genbank_baseline"] & clean & ~dense],
        "baseline_complex_or_partial_context": df[base & df["genbank_baseline"] & ~clean],
        "no_genbank_cds_genome": df[base & ~df["genbank_baseline"]],
    }
    cols = ["accession_version", "genbank_cds_n", "locus_id", "strand", "start0", "end0", "length_aa", "support",
            "phanotate_present", "pyrodigal_gv_present", "best_phrog", "best_category",
            "best_annot", "best_evalue", "block"]
    return pd.concat([v.assign(h4_layer=k) for k, v in layers.items()]).sort_values(
        ["h4_layer", "best_evalue"])[cols + ["h4_layer", "genbank_antisense_context",
                                              "genbank_completeness_ratio"]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--function", default=None)
    parser.add_argument("--reps", type=int, default=2000)
    parser.add_argument("--reuse-locus-table", action="store_true",
                        help="Low-memory mode: read the existing analysis/locus_table.parquet "
                             "(modelling columns only) instead of rebuilding it")
    args = parser.parse_args()
    run = Path(args.run)
    out = run / "analysis"
    out.mkdir(exist_ok=True)
    fn = args.function or run / "phrogs" / "protein_function_primary.parquet"
    lt_path = out / "locus_table.parquet"
    if args.reuse_locus_table and lt_path.exists():
        df = pd.read_parquet(lt_path, columns=SLIM_COLUMNS)
        for col in ("block", "accession_version", "genome_id", "function_status",
                    "a2_slippage_only"):
            df[col] = df[col].astype("category")
        for col in ("support", "support2", "single_support", "short", "overlap_any_ge1",
                    "overlap_same_ge1", "overlap_anti_ge1", "overlap_same_ge30"):
            df[col] = df[col].astype("int8")
        for col in ("log_len", "gc", "log_genome"):
            df[col] = df[col].astype("float32")
        import gc as _gc
        _gc.collect()
    else:
        df = prepare(locus_table(run, args.manifest, fn if Path(fn).exists() else None))
        df.to_parquet(lt_path, index=False)
    metrics, calls, loci = descriptive(run, df)
    metrics.to_csv(out / "pairwise_agreement.csv", index=False)
    calls.to_csv(out / "calls_per_source.csv", index=False)
    write_json(out / "loci_summary.json", loci)
    result = {"h2": h2(df, args.reps)}
    if "function_status" in df:
        result["h1"] = h1(df, args.reps)
        h4(df, lt_path if args.reuse_locus_table else None).to_csv(
            out / "h4_candidates.csv", index=False)
        p = {"H1_support2": result["h1"]["gee"]["terms"]["support2"]["p"],
             "H2_log_len": result["h2"]["gee_linear_loglen"]["terms"]["log_len"]["p"]}
        result["holm_primary"] = holm(p)
        result["holm_note"] = ("H2 primary p uses the linear log-length term of the "
                               "prespecified covariate set; spline model reported for shape.")
    result["caveats"] = [
        "Associations among detected ORFs; not error rates against a gold standard.",
        "Relatedness blocks approximate (TaxMyPhage genus / NCBI genus metadata); no ANI.",
        "PHROGs thresholds are prespecified investigator choices.",
    ]
    write_json(out / "hypotheses_h1_h2.json", result)
    print(json.dumps({"loci": loci, "out": str(out)}, indent=1, default=str))


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    main()
