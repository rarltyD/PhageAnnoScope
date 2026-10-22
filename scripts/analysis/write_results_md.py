"""Write a numbers-only results digest (Markdown) from the final analysis outputs.

    python scripts/analysis/write_results_md.py --analysis results/statistics \
        --robustness results/robustness --descriptive results/descriptive \
        --phold results/phold/h3_summary.json --out docs/RESULTS.md
"""

import argparse
import json
from pathlib import Path

import pandas as pd


SENS = {"exclude_slippage_genomes": "Excluding genomes with ribosomal slippage",
        "length_ge_30aa": "ORFs ≥ 30 aa", "exclude_no_genbank_cds_genomes": "Excluding genomes without GenBank CDS",
        "short_lt100aa_only": "ORFs < 100 aa only", "long_ge100aa_only": "ORFs ≥ 100 aa only",
        "one_genome_per_block": "One genome per relatedness block"}


def pct(x, d=1):
    return f"{100 * x:.{d}f}%"


def ci(v, d=3, scale=1):
    return f"{v[0] * scale:.{d}f}–{v[1] * scale:.{d}f}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--analysis", required=True)
    p.add_argument("--robustness", required=True)
    p.add_argument("--descriptive", required=True)
    p.add_argument("--phold", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    an = Path(a.analysis)
    L = json.loads((an / "loci_summary.json").read_text())
    H = json.loads((an / "hypotheses_h1_h2.json").read_text())
    merge = json.loads((an / "merge_info.json").read_text())
    pw = pd.read_csv(an / "pairwise_agreement.csv")
    h4 = pd.read_csv(an / "h4_candidates.csv.gz")
    out = ["# Results (full cohort)", "",
           "Generated from the analysis outputs in `results/`; numbers only, no interpretation.", "",
           "## Cohort", "",
           f"- Genomes analysed: {merge['total_genomes']:,}; primary analysis: {L['n_genomes']:,}; "
           f"relatedness blocks: {L['n_blocks']:,}",
           f"- Primary stop-loci: {L['n_loci_primary']:,}; supported by both callers: "
           f"{L['support_counts']['2']:,}; by one caller: {L['support_counts']['1']:,} "
           f"(PHANOTATE only {L['single_support_by_caller']['phanotate_only']:,}, "
           f"Pyrodigal-gv only {L['single_support_by_caller']['pyrodigal_gv_only']:,})",
           f"- Median protein length (aa): one caller {L['median_length_aa_by_support']['1']:.0f}, "
           f"both callers {L['median_length_aa_by_support']['2']:.0f}", "",
           "## Pairwise agreement (pooled Jaccard)", "", "| Pair | Exact CDS | Shared stop |",
           "|---|---:|---:|"]
    for r in pw.itertuples():
        out.append(f"| {r.source_a} – {r.source_b} | {r.pooled_exact_jaccard:.3f} | "
                   f"{r.pooled_stop_jaccard:.3f} |")
    h1 = H["h1"]
    b = h1["bootstrap_no_specific_by_support"]
    t = h1["gee"]["terms"]["support2"]
    out += ["", "## H1: caller support and functional uncertainty", "",
            f"- No specific PHROGs function: one caller {pct(b['levels']['1']['proportion'])} "
            f"(95% CI {ci(b['levels']['1']['ci95'], 1, 100)}%), both callers "
            f"{pct(b['levels']['2']['proportion'])} ({ci(b['levels']['2']['ci95'], 1, 100)}%); "
            f"difference {b['difference_pp']['estimate']:.1f} pp ({ci(b['difference_pp']['ci95'], 1)})",
            f"- Adjusted OR (both vs one caller) = {t['odds_ratio']:.3f} (95% CI {ci(t['or_ci95'])}), "
            f"p = {t['p']:.1e}; n = {h1['gee']['n_obs']:,} loci, {h1['gee']['n_clusters']:,} clusters",
            "", "| Sensitivity analysis | OR | 95% CI |", "|---|---:|---|"] + [
                f"| {SENS.get(k, k)} | {v['support2_or']['odds_ratio']:.3f} | {ci(v['support2_or']['or_ci95'])} |"
                for k, v in h1["sensitivity"].items() if "support2_or" in v]
    h2 = H["h2"]
    s = h2["bootstrap_single_support_by_short"]["levels"]
    out += ["", "## H2: length and overlaps", "",
            f"- Single-caller support: proteins < 100 aa {pct(s['1']['proportion'])}, "
            f"≥ 100 aa {pct(s['0']['proportion'])}"]
    labels = {"log_len": "log(protein length), per unit", "gc": "Genome GC fraction, per unit",
              "log_genome": "log(genome length), per unit",
              "overlap_same_ge1": "Same-strand overlap with another locus (≥ 1 nt)",
              "overlap_anti_ge1": "Antisense overlap with another locus (≥ 1 nt)"}
    out += ["", "Logistic GEE for single-caller support (cluster-robust 95% CI):", "",
            "| Covariate | OR | 95% CI |", "|---|---:|---|"]
    for k, v in h2["gee_linear_loglen"]["terms"].items():
        if k != "Intercept":
            out.append(f"| {labels.get(k, k)} | {v['odds_ratio']:.3f} | {ci(v['or_ci95'])} |")
    out.append("")
    ov = h2["sensitivity"]["overlap_ge30nt"]["terms"]["overlap_same_ge1"]
    out.append(f"- Same-strand overlap ≥ 30 nt: OR {ov['odds_ratio']:.2f}")
    out += ["", "## H3: Phold (probability sample of 600 proteins)", ""]
    h3 = json.loads(Path(a.phold).read_text())
    g = h3["gain_unweighted_all_no_specific"]
    out.append(f"- Informative assignment among proteins without a specific function: {pct(g['prop'])} "
               f"(95% CI {ci(g['ci95'], 1, 100)}%), n = {g['n']}; design-weighted "
               f"{pct(h3['gain_weighted_sample_frame']['informative_prop'])}")
    out += ["", "| Caller support | Protein length | PHROGs function | n | Phold informative | Category agreement |",
            "|---|---|---|---:|---:|---:|"]
    for k, v in h3["strata"].items():
        support, length, group = k.split("_", 2)
        agree = pct(v["control_agreement_prop"], 0) if "control_agreement_prop" in v else "—"
        out.append(f"| {support[-1]} caller{'s' if support[-1] == '2' else ''} | "
                   f"{'< 100 aa' if length == 'lt100' else '≥ 100 aa'} | "
                   f"{'none (target)' if group == 'no_specific' else 'specific (control)'} | {v['n']} | "
                   f"{pct(v['informative_prop'], 0)} | {agree} |")
    out += ["", "## H4: ORFs absent from GenBank with a specific PHROGs function", ""]
    layer_names = {"baseline_clean": "Main layer: complete GenBank annotation, no complex context",
                   "baseline_sparse_annotation": "Sparse GenBank annotation (completeness ratio < 0.8)",
                   "baseline_complex_or_partial_context": "Next to frameshift / partial GenBank CDS",
                   "no_genbank_cds_genome": "Genomes without any GenBank CDS"}
    out += ["| Layer | ORFs | Genomes |", "|---|---:|---:|"]
    for layer in layer_names:
        sub = h4[h4.h4_layer == layer]
        out.append(f"| {layer_names[layer]} | {len(sub):,} | {sub.accession_version.nunique():,} |")
    out.append("")
    clean = h4[h4.h4_layer == "baseline_clean"]
    out.append("- Main layer by category: " + ", ".join(
        f"{k} {v:,}" for k, v in clean.best_category.value_counts().head(6).items()))
    rb = pd.read_csv(Path(a.robustness) / "robustness_table.csv").rename(columns={"or": "odds"})
    out += ["", "## H1 robustness (adjusted OR, 95% CI)", "", "| Section | Stratum | OR | 95% CI | n |",
            "|---|---|---:|---|---:|"]
    for r in rb.dropna(subset=["odds"]).itertuples():
        out.append(f"| {r.section} | {r.stratum} | {r.odds:.3f} | {r.lo:.3f}–{r.hi:.3f} | {int(r.n):,} |")
    e = json.loads((Path(a.descriptive) / "extra_summary.json").read_text())
    gb = json.loads((Path(a.descriptive) / "genbank_function_breakdown.json").read_text())
    out += ["", "## Descriptive", "",
            f"- GenBank CDS: {gb['genbank_cds']:,}; named {gb['named']:,}; 'hypothetical' "
            f"{gb['hypothetical_like']:,} ({pct(e['genbank_cds']['hypothetical_like_share'])}); "
            f"'hypothetical' with a specific PHROGs function {gb['hyp_specific']:,} "
            f"({pct(e['genbank_cds']['hypothetical_with_specific_phrogs'])})",
            f"- Unique representative proteins: {e['unique_representative_proteins']['n']:,}",
            "- CDS per kb (median): " + ", ".join(
                f"{k} {v['cds_per_kb_median']}" for k, v in e["gene_density"].items()),
            f"- Alternative starts PHANOTATE vs Pyrodigal-gv: median "
            f"{e['alt_start_phanotate_vs_pyrodigal']['median_diff_codons']:.0f} codons "
            f"(n = {e['alt_start_phanotate_vs_pyrodigal']['n_pairs']:,} pairs)"]
    Path(a.out).write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"written {a.out}")


if __name__ == "__main__":
    main()
