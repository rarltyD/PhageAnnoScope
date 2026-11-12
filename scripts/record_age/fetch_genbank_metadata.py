"""Download GenBank record headers or feature tables for the cohort (NCBI E-utilities).

headers   truncated flat files (header, comments and references; 1 bp of sequence)
features  5-column feature tables (CDS with /inference and /note qualifiers)

One file per batch of accessions; existing files are skipped, so the command can be repeated to
resume.  No e-mail address is sent; the request rate stays below the anonymous NCBI limit.

    python scripts/record_age/fetch_genbank_metadata.py --mode headers \
        --accessions results/descriptive/per_genome_summary.csv --out data/raw/ncbi/headers
    python scripts/record_age/fetch_genbank_metadata.py --mode features \
        --accessions results/descriptive/per_genome_summary.csv --out data/raw/ncbi/features
"""
import argparse
import csv
import time
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MODES = {  # efetch parameters, batch size, file suffix, record marker
    "headers": ({"rettype": "gb", "seq_start": 1, "seq_stop": 1}, 200, ".gb", "\n//"),
    "features": ({"rettype": "ft"}, 100, ".ft", ">Feature"),
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=MODES, required=True)
    p.add_argument("--accessions", required=True,
                   help="CSV with an accession_version column, or a text file with one accession per line")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    params, size, suffix, marker = MODES[a.mode]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = Path(a.accessions)
    if path.suffix == ".csv":
        acc = [r["accession_version"] for r in csv.DictReader(path.open())]
    else:
        acc = path.read_text().split()
    for i in range(0, len(acc), size):
        f = out / f"batch_{i // size:04d}{suffix}"
        if f.exists() and f.stat().st_size > 0:
            continue
        ids = acc[i:i + size]
        data = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(ids), "retmode": "text",
                                       "tool": "PhageAnnoScope", **params}).encode()
        for attempt in range(6):
            try:
                txt = urllib.request.urlopen(URL, data=data, timeout=180).read().decode()
                if txt.count(marker) < len(ids) * 0.95:
                    raise ValueError(f"only {txt.count(marker)} records returned")
                tmp = f.with_suffix(".tmp")
                tmp.write_text(txt)
                tmp.rename(f)
                break
            except Exception as e:
                print(f"batch {i // size} attempt {attempt}: {e}", flush=True)
                time.sleep(5 * (attempt + 1))
        print(f"done {min(i + size, len(acc))}/{len(acc)}", flush=True)
        time.sleep(0.5)


if __name__ == "__main__":
    main()
