"""Additional descriptive statistics for the talk (no recomputation; existing tables only).

    python scripts/analysis/extra_stats.py --run results/full_cohort \
        --manifest manifests/cohort_manifest_validated.csv.gz

Outputs <run>/analysis/extra/*.csv + extra_summary.json.  Descriptive only:
genome-level summaries are medians over genomes; no hypothesis tests here.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

HYPO = re.compile(r"hypothetical|uncharacteri[sz]ed|unknown|putative protein$|^protein$|duf\d+",
                  re.I)


def product_class(product):
    if not isinstance(product, str) or not product:
        return "no_product"
    return "hypothetical_like" if HYPO.search(product) else "named"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--manifest", required=True)
    a = p.parse_args()
    run = Path(a.run)
    out = run / "analysis" / "extra"
    out.mkdir(parents=True, exist_ok=True)
    pq = lambda t: f"read_parquet('{(run / 'parquet' / t).as_posix()}/*.parquet', union_by_name=true)"  # noqa: E731
    con = duckdb.connect()
    con.execute("SET threads=3; SET memory_limit='2GB'")
    fn = pd.read_parquet(run / "phrogs" / "protein_function_primary.parquet",
                         columns=["protein_sha256", "function_status", "best_category"])
    con.register("fn", fn)
    loci = pd.read_parquet(run / "analysis" / "locus_table.parquet",
                           columns=["locus_id", "protein_sha256", "accession_version",
                                    "function_status", "best_category", "single_support", "gc",
                                    "length_bp"])
    meta = pd.read_csv(a.manifest, dtype=str, keep_default_na=False)
    meta = meta[meta.accession_version.isin(loci.accession_version.unique())]
    summary = {}

    # 1. GenBank product names vs PHROGs status of the same GenBank protein.
    gb = con.execute(f"""
        SELECT c.call_id, json_extract_string(c.qualifiers_json, '$.product[0]') AS product,
               o.protein_sha256
        FROM {pq('cds_calls')} c
        LEFT JOIN (SELECT call_id, min(protein_sha256) protein_sha256
                   FROM {pq('protein_origins')} GROUP BY call_id) o USING (call_id)
        WHERE c.source = 'genbank'""").df()
    gb["product_class"] = gb["product"].map(product_class)
    gb = gb.merge(fn, on="protein_sha256", how="left")
    gb["function_status"] = gb["function_status"].fillna("no_valid_protein_or_not_searched")
    t1 = gb.groupby(["product_class", "function_status"]).size().unstack(fill_value=0)
    t1.to_csv(out / "genbank_product_vs_phrogs.csv")
    hyp = gb[gb.product_class == "hypothetical_like"]
    named = gb[gb.product_class == "named"]
    summary["genbank_cds"] = {
        "n": int(len(gb)),
        "hypothetical_like_share": float((gb.product_class == "hypothetical_like").mean()),
        "hypothetical_with_specific_phrogs": float((hyp.function_status == "specific_function")
                                                   .mean()),
        "named_without_informative_phrogs_hit": float((named.function_status
                                                       == "no_informative_hit").mean()),
        "top_hypothetical_to_category": hyp[hyp.function_status == "specific_function"]
        .best_category.value_counts().head(8).to_dict()}

    # 2. Start codon usage and gene density per source (primary rotation only).
    sc = con.execute(f"""
        SELECT source, upper(start_codon) AS start_codon, count(*) n
        FROM {pq('cds_calls')} WHERE coalesce(rotation_shift, 0) = 0
          AND NOT coalesce(partial_5p, false) GROUP BY ALL""").df()
    sc["share"] = sc.n / sc.groupby("source").n.transform("sum")
    sc.sort_values(["source", "n"], ascending=[True, False]).to_csv(
        out / "start_codon_usage.csv", index=False)
    summary["start_codons"] = {
        s: g.set_index("start_codon").share.round(4).head(4).to_dict()
        for s, g in sc.sort_values("n", ascending=False).groupby("source")}
    dens = con.execute(f"""
        SELECT c.source, g.genome_id, count(*) * 1000.0 / any_value(g.length_bp) AS cds_per_kb,
               median(c.feature_nt_length / 3 - 1) AS median_aa
        FROM {pq('cds_calls')} c JOIN (SELECT DISTINCT genome_id, length_bp
                                       FROM {pq('genomes')}) g USING (genome_id)
        WHERE coalesce(c.rotation_shift, 0) = 0 GROUP BY ALL""").df()
    d = dens.groupby("source").agg(genomes=("genome_id", "nunique"),
                                   cds_per_kb_median=("cds_per_kb", "median"),
                                   cds_per_kb_q1=("cds_per_kb", lambda s: s.quantile(.25)),
                                   cds_per_kb_q3=("cds_per_kb", lambda s: s.quantile(.75)),
                                   protein_len_median=("median_aa", "median"))
    d.to_csv(out / "gene_density_by_source.csv")
    summary["gene_density"] = d.round(3).to_dict(orient="index")

    # 3. Genome-level agreement and function by lifestyle / host genus / family / record year.
    gm = con.execute(f"""
        SELECT g.accession_version, m.source_a, m.source_b, m.stop_jaccard
        FROM {pq('genome_metrics')} m
        JOIN (SELECT DISTINCT genome_id, accession_version FROM {pq('genomes')}) g USING (genome_id)
        WHERE m.metric_scope = 'source_pair'""").df()
    gm["pair"] = gm.source_a + "_vs_" + gm.source_b
    jac = gm.pivot_table(index="accession_version", columns="pair", values="stop_jaccard")
    years = con.execute(f"""SELECT DISTINCT accession_version,
        try_cast(right(json_extract_string(annotations_json, '$.date'), 4) AS INTEGER) AS year
        FROM {pq('genomes')}""").df()
    per_genome = loci.assign(
        no_specific=(loci.function_status != "specific_function").astype(float)
    ).groupby("accession_version").agg(
        loci=("locus_id", "size"), single_support=("single_support", "mean"),
        no_specific=("no_specific", "mean"), gc=("gc", "first"), length_bp=("length_bp", "first"))
    per_genome = per_genome.join(jac).reset_index().merge(years, on="accession_version",
                                                          how="left")
    per_genome = per_genome.merge(
        meta[["accession_version", "phageleads_lifestyle", "metadata_host", "ncbi_family_metadata",
              "taxmyphage_genus"]], on="accession_version", how="left")
    per_genome["year_bin"] = pd.cut(per_genome.year, [0, 2010, 2015, 2020, 2023, 2100],
                                    labels=["≤2010", "2011–2015", "2016–2020", "2021–2023",
                                            "2024–2026"])
    per_genome.to_csv(out / "per_genome_summary.csv", index=False)
    cols = ["single_support", "no_specific", "genbank_vs_pyrodigal_gv", "genbank_vs_phanotate",
            "phanotate_vs_pyrodigal_gv"]

    def by(col, top=None, min_n=20):
        g = per_genome.groupby(col, observed=True)
        tab = g[cols].median().join(g.size().rename("genomes"))
        tab = tab[tab.genomes >= min_n].sort_values("genomes", ascending=False)
        return tab.head(top) if top else tab

    for name, col, top in (("lifestyle", "phageleads_lifestyle", None),
                           ("host_genus_top15", "metadata_host", 15),
                           ("family_top12", "ncbi_family_metadata", 12),
                           ("record_year", "year_bin", None)):
        tab = by(col, top)
        tab.to_csv(out / f"by_{name}.csv")
        summary[f"by_{name}"] = tab.round(3).reset_index().to_dict(orient="records")

    # 4. Protein "dark matter" over unique proteins of primary representative loci.
    rep = loci.drop_duplicates("protein_sha256")
    summary["unique_representative_proteins"] = {
        "n": int(len(rep)),
        "status_share": rep.function_status.value_counts(normalize=True).round(4).to_dict(),
        "specific_category_counts": rep[rep.function_status == "specific_function"]
        .best_category.value_counts().to_dict()}

    # 5. Alternative starts PHANOTATE vs Pyrodigal-gv: length difference in codons.
    alt = con.execute(f"""
        SELECT abs(ca.feature_nt_length - cb.feature_nt_length) / 3 AS diff_codons,
               ca.feature_nt_length > cb.feature_nt_length AS a_longer, m.source_a
        FROM {pq('matches')} m
        JOIN {pq('cds_calls')} ca ON ca.call_id = m.call_id_a
        JOIN {pq('cds_calls')} cb ON cb.call_id = m.call_id_b
        WHERE m.relation = 'alternative_start' AND m.source_a = 'phanotate'
          AND m.source_b = 'pyrodigal_gv'""").df()
    summary["alt_start_phanotate_vs_pyrodigal"] = {
        "n_pairs": int(len(alt)),
        "median_diff_codons": float(alt.diff_codons.median()) if len(alt) else None,
        "phanotate_longer_share": float(alt.a_longer.mean()) if len(alt) else None}
    alt.diff_codons.describe().to_csv(out / "alt_start_diff_codons.csv")

    (out / "extra_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                                       default=str))
    print(json.dumps({k: v for k, v in summary.items() if not k.startswith("by_")},
                     indent=1, ensure_ascii=False, default=str)[:4000])


if __name__ == "__main__":
    main()
