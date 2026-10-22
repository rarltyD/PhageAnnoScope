#!/usr/bin/env python3
"""Download the official INPHARED2 TSV export into a NEW snapshot directory."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
import urllib.parse
import urllib.request

FIELDS = ["accession", "description", "genome_length_bp", "mol_gc_percent", "host",
          "baltimore_group", "lowest_taxa", "genus", "subfamily", "family",
          "taxonomic_order", "taxonomic_class", "phylum", "kingdom", "realm",
          "isolation_host", "checkv_quality", "checkv_miuvig_quality",
          "checkv_completeness", "checkv_completeness_method", "taxmyphage_genus",
          "taxmyphage_species", "taxmyphage_message", "phageleads_lifestyle",
          "phageleads_problematic_genes"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--max-bytes", type=int, default=64 * 1024 * 1024)
    p.add_argument("--timeout", type=int, default=60)
    a = p.parse_args()
    if a.out.exists():
        raise SystemExit("Output already exists. Use a NEW directory to preserve frozen snapshots.")
    a.out.mkdir(parents=True)
    url = "https://phagebank.uk/inphared2/?" + urllib.parse.urlencode(
        {"export": "inphared2_tsv", "download_field": FIELDS}, doseq=True)
    started = datetime.now(timezone.utc).isoformat()
    part = a.out / "inphared2_metadata.tsv.part"
    req = urllib.request.Request(url, headers={"User-Agent": "PhageAnnoScope/metadata-research"})
    try:
        with urllib.request.urlopen(req, timeout=a.timeout) as r, part.open("wb") as f:
            headers, total = dict(r.headers), 0
            if "tab-separated-values" not in r.headers.get("Content-Type", ""):
                raise ValueError("Expected TSV; received a different content type")
            while True:
                chunk = r.read(1024 * 1024)
                if not chunk:
                    break
                total += len(chunk)
                if total > a.max_bytes:
                    raise ValueError("Export exceeds explicit byte limit")
                f.write(chunk)
            expected = r.headers.get("Content-Length")
            if expected and total != int(expected):
                raise ValueError("Truncated transfer")
        payload = part.read_bytes()
        rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8-sig")), delimiter="\t"))
        if not rows or "Accession" not in rows[0]:
            raise ValueError("Missing source records or Accession column")
        if len({r["Accession"] for r in rows}) != len(rows):
            raise ValueError("Duplicate source accessions")
        meta = {"url": url, "started_at_utc": started,
                "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                "response_headers": headers, "n_records": len(rows),
                "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}
        part.rename(a.out / "inphared2_metadata.tsv")
        (a.out / "download_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
        print(json.dumps({k: meta[k] for k in ["n_records", "sha256", "bytes"]}))
    except Exception as exc:
        (a.out / "download_error.json").write_text(json.dumps({"url": url, "error": str(exc), "started_at_utc": started}) + "\n")
        raise


if __name__ == "__main__":
    main()
