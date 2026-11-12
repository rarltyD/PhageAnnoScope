"""Download full GenBank records for the cohort and store one file per accession.

Writes <out>/<accession>.gbk (the record bytes as returned by efetch, gbwithparts), matching the
genbank_path layout of manifests/cohort_manifest_validated.csv.gz; the files are byte-identical to
those of the study unless NCBI has updated a record since (see make_batch_manifests.py).
Resumable: existing files are skipped.  No e-mail is sent to NCBI; the request rate stays below the
anonymous limit.

    python scripts/pipeline/fetch_genbank_records.py manifests/accessions_batch*.txt data/raw/ncbi/primary
"""
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

acc_files, out = sys.argv[1:-1], Path(sys.argv[-1])
out.mkdir(parents=True, exist_ok=True)
acc = []
for f in acc_files:
    acc += [a.strip() for a in Path(f).read_text().split() if a.strip()]
acc = list(dict.fromkeys(acc))
todo = [a for a in acc if not (out / (a.rsplit(".", 1)[0] + ".gbk")).exists()]
print(f"{len(acc)} accessions, {len(todo)} to download", flush=True)
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
B = 50
for i in range(0, len(todo), B):
    ids = todo[i:i + B]
    data = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(ids), "rettype": "gbwithparts",
                                   "retmode": "text", "tool": "PhageAnnoScope"}).encode()
    for attempt in range(8):
        try:
            txt = urllib.request.urlopen(URL, data=data, timeout=600).read().decode()
            recs = [r + "//\n" for r in re.split(r"^//\n", txt, flags=re.M) if r.strip()]
            got = {}
            for r in recs:
                m = re.search(r"^VERSION\s+(\S+)", r, re.M)
                if m:
                    got[m.group(1)] = r.lstrip("\n")
            missing = [a for a in ids if a not in got]
            if len(missing) > len(ids) * 0.05:
                raise ValueError(f"{len(missing)} records missing")
            for a in ids:
                if a in got:
                    p = out / (a.rsplit(".", 1)[0] + ".gbk")
                    p.with_suffix(".tmp").write_text(got[a]); p.with_suffix(".tmp").rename(p)
            if missing:
                print("missing:", ",".join(missing), flush=True)
            break
        except Exception as e:
            print(f"chunk {i // B} attempt {attempt}: {e}", flush=True)
            time.sleep(10 * (attempt + 1))
    else:
        print(f"chunk {i // B} FAILED", flush=True)
    if (i // B) % 20 == 0:
        print(f"done {min(i + B, len(todo))}/{len(todo)}", flush=True)
    time.sleep(0.5)
print("fetch finished", flush=True)
