"""Consensus re-annotation of GenBank phage records from two gene callers and PHROGs.

Rules were fixed before the results were inspected.  Genes are keyed by genome, strand and the
coordinates of the stop codon (the stop-locus definition of the study); only primary-eligible
calls are assessed (partial, joined, slippage and origin-spanning models are left as they are).

Gene level
  confirmed          GenBank CDS whose stop is predicted by both callers
  supported          GenBank CDS whose stop is predicted by one caller
  gb_only_function   GenBank CDS predicted by neither caller, protein has a specific PHROGs function
  gb_only_unsupported GenBank CDS predicted by neither caller, no specific function  -> review
  add_A              both callers, absent from GenBank, specific PHROGs function     -> add
  add_B              both callers, absent from GenBank, no specific function, >=100 aa (reported only)
  start_correction   shared stop, callers agree on one start, GenBank start differs  (reported only)
Record level
  unannotated        no CDS in the GenBank record
  major              (gb_only_unsupported + add_A) >= 10% of GenBank CDS, or discordant with both callers
                     (record-age analysis)
  minor              at least one gb_only_unsupported or add_A gene
  consistent         none
Outputs: per-record table, per-gene proposals, consensus GFF3 for all records, summary JSON.
Exploratory analysis added after the primary analysis (see docs/PROTOCOL.md).

    python scripts/reannotation/reannotate.py results/full_cohort results/record_age/genome_table.csv.gz \
        results/reannotation
"""
import gzip
import json
import sys
from pathlib import Path

import duckdb
import pandas as pd

run, part_a, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
pq = lambda t: (run / "parquet" / t / "*.parquet").as_posix()
fn_path = (run / "phrogs" / "protein_function_primary.parquet").as_posix()
con = duckdb.connect()
con.execute("SET memory_limit='4GB'; SET threads=3; SET preserve_insertion_order=false")

genomes_a = pd.read_csv(part_a, usecols=["accession_version", "submitted", "period", "unannotated",
                                          "disc_annotated", "block", "caller", "seq_tech"])
con.register("cohort", genomes_a)

# One row per primary-eligible CDS call with its stop key, start and protein function.
con.execute(f"""
CREATE TEMP TABLE calls AS
SELECT c.genome_id, c.accession_version, c.source, c.call_id, c.strand,
       c.stop_triplet::VARCHAR AS stop_key,
       list_min(list_transform(c.segments, s -> s.start0)) AS lo,
       list_max(list_transform(c.segments, s -> s.end0)) AS hi,
       c.feature_nt_length // 3 - 1 AS length_aa, c.protein_sha256,
       f.function_status, f.best_phrog, f.best_category, f.best_annot, f.best_evalue
FROM read_parquet('{pq("cds_calls")}', union_by_name=true) c
JOIN cohort USING (accession_version)
LEFT JOIN read_parquet('{fn_path}') f USING (protein_sha256)
WHERE c.eligible_primary
""")
con.execute("""
CREATE TEMP TABLE calls2 AS
SELECT *, CASE WHEN strand = 1 THEN lo ELSE hi END AS start_coord FROM calls
""")
# Per stop key: which sources predict it, their starts, and the best function among the models.
con.execute("""
CREATE TEMP TABLE keys AS
SELECT genome_id, any_value(accession_version) AS accession_version, strand, stop_key,
       bool_or(source = 'genbank') AS gb, bool_or(source = 'phanotate') AS pha,
       bool_or(source = 'pyrodigal_gv') AS pyr,
       min(start_coord) FILTER (WHERE source = 'genbank') AS start_gb,
       min(start_coord) FILTER (WHERE source = 'phanotate') AS start_pha,
       min(start_coord) FILTER (WHERE source = 'pyrodigal_gv') AS start_pyr,
       max(length_aa) AS length_aa,
       min(lo) AS lo, max(hi) AS hi,
       bool_or(function_status = 'specific_function') FILTER (WHERE source = 'genbank') AS gb_specific,
       bool_or(function_status = 'specific_function') FILTER (WHERE source <> 'genbank') AS caller_specific,
       arg_min(best_annot, best_evalue) FILTER (WHERE function_status = 'specific_function') AS annot,
       arg_min(best_category, best_evalue) FILTER (WHERE function_status = 'specific_function') AS category,
       arg_min(best_phrog, best_evalue) FILTER (WHERE function_status = 'specific_function') AS phrog
FROM calls2
GROUP BY genome_id, strand, stop_key
""")
# Caller support counts only at primary-eligible stop-loci (the study's unit), so caller totals
# match the H1-H4 analysis exactly.
con.execute(f"""
CREATE TEMP TABLE keys_el AS
SELECT k.* REPLACE (k.pha AND el.genome_id IS NOT NULL AS pha, k.pyr AND el.genome_id IS NOT NULL AS pyr)
FROM keys k
LEFT JOIN (SELECT DISTINCT genome_id, strand, stop_triplet::VARCHAR AS stop_key
           FROM read_parquet('{pq("stop_loci")}', union_by_name=true) WHERE eligible_primary) el
  USING (genome_id, strand, stop_key)
""")
con.execute("DELETE FROM keys_el WHERE NOT gb AND NOT pha AND NOT pyr")
con.execute("""
CREATE TEMP TABLE genes AS
SELECT *, (pha::INT + pyr::INT) AS n_callers,
  CASE
    WHEN gb AND pha AND pyr THEN 'confirmed'
    WHEN gb AND (pha OR pyr) THEN 'supported'
    WHEN gb AND coalesce(gb_specific, false) THEN 'gb_only_function'
    WHEN gb THEN 'gb_only_unsupported'
    WHEN pha AND pyr AND coalesce(caller_specific, false) THEN 'add_A'
    WHEN pha AND pyr AND length_aa >= 100 THEN 'add_B'
    WHEN pha AND pyr THEN 'two_caller_other'
    WHEN coalesce(caller_specific, false) THEN 'one_caller_function'
    ELSE 'one_caller_other'
  END AS category_rule,
  (gb AND pha AND pyr AND start_pha = start_pyr AND start_gb <> start_pha) AS start_correction
FROM keys_el
""")

genes = con.execute("SELECT * FROM genes").df()
counts = genes.pivot_table(index="accession_version", columns="category_rule", values="stop_key",
                           aggfunc="size", fill_value=0)
for c in ("confirmed", "supported", "gb_only_function", "gb_only_unsupported", "add_A", "add_B",
          "two_caller_other", "one_caller_function", "one_caller_other"):
    if c not in counts:
        counts[c] = 0
rec = genomes_a.set_index("accession_version").join(counts, how="left").fillna(
    {c: 0 for c in counts.columns})
rec["start_corrections"] = genes[genes.start_correction == True].groupby("accession_version").size()
rec["start_corrections"] = rec.start_corrections.fillna(0).astype(int)
rec["gb_cds_assessed"] = rec[["confirmed", "supported", "gb_only_function", "gb_only_unsupported"]].sum(axis=1)
rec["flagged"] = rec.gb_only_unsupported + rec.add_A
rec["flagged_share"] = rec.flagged / rec.gb_cds_assessed.where(rec.gb_cds_assessed > 0)


def record_class(r):
    if r.unannotated == 1:
        return "unannotated"
    if r.flagged_share >= 0.10 or r.disc_annotated == 1:
        return "major"
    if r.flagged > 0:
        return "minor"
    return "consistent"


rec["reannotation_class"] = rec.apply(record_class, axis=1)
rec.reset_index().to_csv(out / "records_reannotation.csv.gz", index=False)

# Per-gene proposals (everything except plain confirmed/supported genes and one-caller noise).
prop = genes[genes.category_rule.isin(["gb_only_unsupported", "gb_only_function", "add_A", "add_B"])
             | (genes.start_correction == True)].copy()
prop["start1"] = prop.lo + 1
prop["end1"] = prop.hi
prop[["accession_version", "category_rule", "start_correction", "strand", "start1", "end1",
      "length_aa", "start_gb", "start_pha", "start_pyr", "phrog", "category", "annot"]] \
    .sort_values(["accession_version", "start1"]).to_csv(out / "gene_proposals.csv.gz", index=False)

# Consensus GFF3: two-caller genes, one-caller genes with a specific function, GenBank-only genes
# with a specific function.  Coordinates are those of the longest model at the stop.
keep = genes[genes.category_rule.isin(["confirmed", "supported", "gb_only_function", "add_A", "add_B",
                                        "two_caller_other", "one_caller_function"])]
with gzip.open(out / "consensus_annotation.gff3.gz", "wt") as g:
    g.write("##gff-version 3\n")
    for r in keep.sort_values(["accession_version", "lo"]).itertuples():
        ev = {"confirmed": "GenBank+2 callers", "supported": "GenBank+1 caller",
              "gb_only_function": "GenBank+PHROGs", "add_A": "2 callers+PHROGs (absent from GenBank)",
              "add_B": "2 callers (absent from GenBank)", "two_caller_other": "2 callers (absent from GenBank)",
              "one_caller_function": "1 caller+PHROGs"}[r.category_rule]
        product = r.annot if isinstance(r.annot, str) and r.annot else "hypothetical protein"
        attrs = (f"ID={r.accession_version}_{r.lo + 1}_{'+' if r.strand == 1 else '-'};"
                 f"product={product.replace(';', ',')};evidence={ev}")
        if isinstance(r.phrog, str):
            attrs += f";phrog={r.phrog};phrog_category={str(r.category).replace(';', ',')}"
        g.write(f"{r.accession_version}\tPhageAnnoScope\tCDS\t{r.lo + 1}\t{r.hi}\t.\t"
                f"{'+' if r.strand == 1 else '-'}\t0\t{attrs}\n")

summ = {
    "genomes": int(len(rec)),
    "gene_categories": genes.category_rule.value_counts().to_dict(),
    "start_corrections": int((genes.start_correction == True).sum()),
    "record_classes": rec.reannotation_class.value_counts().to_dict(),
    "records_with_add_A": int((rec.add_A > 0).sum()),
    "records_with_gb_only_unsupported": int((rec.gb_only_unsupported > 0).sum()),
    "consensus_genes": int(len(keep)),
    "add_A_by_category": genes[genes.category_rule == "add_A"].category.value_counts().head(12).to_dict(),
    "gb_only_unsupported_length_aa_median": float(genes[genes.category_rule == "gb_only_unsupported"].length_aa.median()),
    "class_by_period": pd.crosstab(rec.period, rec.reannotation_class, normalize="index").round(4).to_dict(orient="index"),
}
(out / "reannotation_summary.json").write_text(json.dumps(summ, indent=2, ensure_ascii=False, default=str))
print(json.dumps({k: summ[k] for k in ("genomes", "gene_categories", "start_corrections", "record_classes",
                                        "records_with_add_A", "records_with_gb_only_unsupported",
                                        "consensus_genes")}, indent=1, default=str))
