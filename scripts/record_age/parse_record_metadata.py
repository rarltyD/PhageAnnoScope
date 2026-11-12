"""Parse truncated GenBank headers and feature tables into one per-record metadata table.

Outputs one row per accession_version: submission date (earliest Direct Submission),
last modification date (LOCUS), sequencing technology, assembly method, annotation
pipeline from structured comments, and the gene-calling tool named in CDS
/inference and /note qualifiers.

    python scripts/record_age/parse_record_metadata.py data/raw/ncbi/headers data/raw/ncbi/features \
        results/record_age/record_metadata.csv.gz
"""
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import pandas as pd

hdr_dir, ft_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])

# Gene-calling tool families, searched in CDS /inference and /note qualifiers.
TOOLS = {
    "PGAP (GeneMarkS-2+)": r"genemarks-2\+",
    "GeneMark": r"genemark(?!s-2\+)",
    "Glimmer": r"glimmer",
    "Prodigal": r"prodigal",
    "PHANOTATE": r"phanotate",
    "RAST": r"\brast\b",
    "DNA Master": r"dna ?master",
}
COMMENT_TOOLS = {
    "PGAP (GeneMarkS-2+)": r"prokaryotic genome annotation pipeline|\bpgap\b",
    "Pharokka": r"pharokka",
    "Prokka": r"prokka",
    "RAST": r"\brast\b|myrast",
    "DNA Master": r"dna ?master",
    "PHANOTATE": r"phanotate",
    "GeneMark": r"genemark",
    "Glimmer": r"glimmer",
    "Prodigal": r"prodigal",
}


def date(s):
    return datetime.strptime(s.title(), "%d-%b-%Y").date() if s else None


def parse_headers():
    rows = []
    for f in sorted(hdr_dir.glob("*.gb")):
        for rec in f.read_text().split("\n//"):
            m = re.search(r"^VERSION\s+(\S+)", rec, re.M)
            if not m:
                continue
            locus = re.search(r"^LOCUS.*?(\d{2}-[A-Z]{3}-\d{4})\s*$", rec, re.M)
            subs = [date(d) for d in re.findall(r"Submitted \((\d{2}-[A-Z]{3}-\d{4})\)", rec)]
            comment = " ".join(re.findall(r"^COMMENT\s+(.*?)(?=^\S)", rec, re.M | re.S))
            comment = re.sub(r"\s+", " ", comment)

            def sc(key):
                v = re.search(rf"{key}\s*::\s*(.*?)(?= \S[\w /-]*::| ##|$)", comment)
                return v.group(1).strip() if v else None

            rows.append({
                "accession_version": m.group(1),
                "submitted": min(subs) if subs else None,
                "last_modified": date(locus.group(1)) if locus else None,
                "sequencing_technology": sc("Sequencing Technology"),
                "assembly_method": sc("Assembly Method"),
                "annotation_pipeline": sc("Annotation Pipeline"),
                "comment_tools": ";".join(t for t, p in COMMENT_TOOLS.items() if re.search(p, comment, re.I)),
                "sequence_updated": bool(re.search(r"sequence updated", comment, re.I)),
            })
    return pd.DataFrame(rows)


def parse_features():
    rows = []
    for f in sorted(ft_dir.glob("*.ft")):
        for block in f.read_text().split(">Feature ")[1:]:
            acc = re.search(r"\|([A-Z]+_?\d+\.\d+)\|", block.split("\n", 1)[0])
            if not acc:
                continue
            n_cds = len(re.findall(r"^\S+\t\S+\tCDS\s*$", block, re.M))
            quals = " ".join(re.findall(r"^\t\t\t(?:inference|note)\t(.*)$", block, re.M))
            hits = Counter({t: len(re.findall(p, quals, re.I)) for t, p in TOOLS.items()})
            hits = +hits
            rows.append({"accession_version": acc.group(1), "gb_cds": n_cds,
                         "feature_tools": ";".join(f"{t}:{n}" for t, n in hits.most_common()),
                         "feature_tool_top": hits.most_common(1)[0][0] if hits else None})
    return pd.DataFrame(rows)


def gene_caller(r):
    """Feature-level evidence first (dominant tool), then the comment, else not stated."""
    if isinstance(r.annotation_pipeline, str) and re.search(r"prokaryotic genome annotation", r.annotation_pipeline, re.I):
        return "PGAP (GeneMarkS-2+)"
    if isinstance(r.feature_tool_top, str):
        return r.feature_tool_top
    if isinstance(r.comment_tools, str) and r.comment_tools:
        return r.comment_tools.split(";")[0]
    return "Not stated"


def seq_tech(s):
    if not isinstance(s, str):
        return "Not stated"
    s = s.lower()
    long_read = re.search(r"nanopore|minion|ont\b|gridion|pacbio|sequel", s)
    short = re.search(r"illumina|miseq|hiseq|nextseq|novaseq|ion ?torrent|bgi|dnbseq|mgi", s)
    if long_read and short:
        return "Hybrid short+long"
    if long_read:
        return "Long-read"
    if re.search(r"454|roche|gs flx", s):
        return "454"
    if short:
        return "Illumina/Ion/MGI"
    if "sanger" in s:
        return "Sanger"
    return "Other"


h, ft = parse_headers(), parse_features()
d = h.merge(ft, on="accession_version", how="outer")
d["gene_caller"] = d.apply(gene_caller, axis=1)
d["seq_tech"] = d.sequencing_technology.map(seq_tech)
d.to_csv(out, index=False)
print(len(h), "headers;", len(ft), "feature tables;", len(d), "rows")
print(d.gene_caller.value_counts().to_string())
print(d.seq_tech.value_counts().to_string())
