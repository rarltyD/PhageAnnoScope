"""Functional breakdown of GenBank CDS: named vs 'hypothetical', by PHROGs category (fig10).

    python scripts/analysis/genbank_function_breakdown.py --run results/full_cohort --out figures

For every GenBank CDS: product class (named / hypothetical-like, same regex as
extra_stats.py) and the PHROGs-HMM result of the same GenBank protein.  Writes
<run>/analysis/extra/genbank_function_breakdown.csv, top annotations of
re-annotatable 'hypothetical' proteins and figure 10.
"""

import argparse
import json
import sys
from pathlib import Path

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from extra_stats import product_class  # noqa: E402
from make_figures import AQUA, BLUE, INK2, save  # noqa: E402

CAT_LABEL = {"head and packaging": "Head and packaging", "tail": "Tail",
             "connector": "Connector", "DNA, RNA and nucleotide metabolism":
             "DNA, RNA and nucleotide metabolism", "transcription regulation":
             "Transcription regulation", "lysis": "Lysis",
             "integration and excision": "Integration and excision",
             "moron, auxiliary metabolic gene and host takeover": "Moron, AMG, host takeover",
             "other": "Other"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    run, out = Path(a.run), Path(a.out)
    pq = lambda t: f"read_parquet('{(run / 'parquet' / t).as_posix()}/*.parquet', union_by_name=true)"  # noqa: E731
    con = duckdb.connect()
    con.execute("SET memory_limit='2GB'; SET threads=3")
    fn = pd.read_parquet(run / "phrogs" / "protein_function_primary.parquet",
                         columns=["protein_sha256", "function_status", "best_category",
                                  "best_annot"])
    con.register("fn", fn)
    gb = con.execute(f"""
        SELECT json_extract_string(c.qualifiers_json, '$.product[0]') AS product,
               f.function_status, f.best_category, f.best_annot
        FROM {pq('cds_calls')} c
        LEFT JOIN (SELECT call_id, min(protein_sha256) protein_sha256
                   FROM {pq('protein_origins')} GROUP BY call_id) o ON o.call_id = c.call_id
        LEFT JOIN fn f ON f.protein_sha256 = o.protein_sha256
        WHERE c.source = 'genbank'""").df()
    gb["product_class"] = gb["product"].map(product_class)
    named = gb[gb.product_class == "named"]
    hyp = gb[gb.product_class == "hypothetical_like"]
    spec_named = named[named.function_status == "specific_function"]
    spec_hyp = hyp[hyp.function_status == "specific_function"]
    tab = pd.DataFrame({
        "named_with_specific_phrogs": spec_named.best_category.value_counts(),
        "hypothetical_reannotatable": spec_hyp.best_category.value_counts()}).fillna(0).astype(int)
    tab["named_share"] = tab.named_with_specific_phrogs / tab.named_with_specific_phrogs.sum()
    tab["hyp_share"] = tab.hypothetical_reannotatable / tab.hypothetical_reannotatable.sum()
    tab = tab.sort_values("hypothetical_reannotatable", ascending=False)
    extra = run / "analysis" / "extra"
    tab.to_csv(extra / "genbank_function_breakdown.csv")
    top = spec_hyp.best_annot.value_counts().head(15)
    top.to_csv(extra / "hypothetical_reannotatable_top_annotations.csv")
    summary = {"genbank_cds": int(len(gb)), "named": int(len(named)),
               "hypothetical_like": int(len(hyp)), "no_product": int((gb.product_class
                                                                       == "no_product").sum()),
               "named_specific": int(len(spec_named)), "hyp_specific": int(len(spec_hyp)),
               "top_reannotatable": top.to_dict()}
    (extra / "genbank_function_breakdown.json").write_text(json.dumps(summary, indent=2,
                                                                      ensure_ascii=False))
    # Figure 10: two panels, Russian category labels.
    cats = list(tab.index)
    y = np.arange(len(cats))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, col, share, colr, title in (
            (axes[0], "named_with_specific_phrogs", "named_share", BLUE,
             "Named in GenBank, specific PHROGs function\n(n = " + f"{len(spec_named):,}" + ")"),
            (axes[1], "hypothetical_reannotatable", "hyp_share", AQUA,
             "'Hypothetical' in GenBank, specific PHROGs function\n(n = " + f"{len(spec_hyp):,}"+ ")")):
        ax.barh(y, tab[share], color=colr, height=0.7)
        for yi, (n, s) in enumerate(zip(tab[col], tab[share])):
            ax.text(s + 0.005, yi, f"{s:.0%} ({n:,})", va="center",
                    fontsize=8, color=INK2)
        ax.set_title(title, fontsize=10)
        ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
        ax.set_xlim(0, max(tab[share].max() * 1.35, 0.1))
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(y, [CAT_LABEL.get(c, c) for c in cats], fontsize=9)
    axes[0].invert_yaxis()
    fig.tight_layout()
    save(fig, out, "fig10_genbank_functional_categories")
    print(tab.to_string())
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
