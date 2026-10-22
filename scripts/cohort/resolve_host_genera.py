"""Protocol amendment A1 (2026-10-18): host domain from INPHARED2 metadata_host.

Each unique ``metadata_host`` value of primary candidates is resolved once with
the project's audited resolver ``hosts.resolve_host_taxonomy`` (exact
[Scientific Name] search + taxonomy efetch; immutable raw XML + SHA-256).

Exact names that are homonyms across kingdoms (e.g. a bacterial and an animal
genus) return >1 taxon and stay unresolved in the strict resolver.  For those
only, a second, explicitly recorded query is made within each prokaryotic
domain subtree: the domain is accepted only if exactly one taxon is found in
exactly one of Bacteria/Archaea.  Placeholders ("Unspecified") stay unresolved.

Output: one row per metadata_host value (CSV), raw evidence under --out.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

from phageannoscope.acquire import _request
from phageannoscope.hosts import _search_taxids, resolve_host_taxonomy


def subtree_query(name, domain, raw_dir):
    term = f'"{name}"[Scientific Name] AND "{domain}"[Subtree]'
    data, meta = _request(
        "esearch.fcgi",
        {"db": "taxonomy", "term": term, "retmode": "xml", "retmax": 2},
        email=os.environ.get("NCBI_EMAIL"),
        api_key=os.environ.get("NCBI_API_KEY"),
    )
    if data is None:
        return None, None, term, None
    digest = hashlib.sha256(data).hexdigest()
    path = raw_dir / f"{digest}.xml"
    if not path.exists():
        path.write_bytes(data)
    count, ids = _search_taxids(data)
    return count, ids, term, str(path)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    out = Path(args.out)
    (out / "subtree_raw").mkdir(parents=True, exist_ok=True)
    with open(args.manifest, encoding="utf-8") as handle:
        names = sorted(
            {
                r["metadata_host"]
                for r in csv.DictReader(handle)
                if r.get("cohort_route") == "primary_candidate"
            }
        )
    rows = []
    for i, name in enumerate(names):
        res = resolve_host_taxonomy(host_name=name, out=out, data_root=out)
        row = {
            "metadata_host": name,
            "lookup_status": res.get("lookup_status"),
            "host_domain": res.get("host_domain", ""),
            "host_domain_status": res.get("host_domain_status"),
            "host_taxid": res.get("host_taxid", ""),
            "host_scientific_name": res.get("host_scientific_name", ""),
            "host_taxonomy_rank": res.get("host_taxonomy_rank", ""),
            "evidence_level": res.get("host_evidence_level", ""),
            "fetch_evidence": res.get("fetch_evidence", ""),
            "fetch_sha256": res.get("fetch_sha256", ""),
            "subtree_evidence_json": "",
        }
        if res.get("lookup_status") == "ambiguous":
            found = {}
            evidence = []
            for domain in ("Bacteria", "Archaea"):
                count, ids, term, path = subtree_query(name, domain, out / "subtree_raw")
                evidence.append({"term": term, "count": count, "ids": ids, "raw": path})
                if count:
                    found[domain] = (count, ids)
            row["subtree_evidence_json"] = json.dumps(evidence)
            if len(found) == 1:
                domain, (count, ids) = next(iter(found.items()))
                if count == 1:
                    row.update(
                        host_domain=domain,
                        host_domain_status="verified",
                        host_taxid=ids[0],
                        lookup_status="verified_prokaryote_subtree",
                        evidence_level="taxonomy_exact_name_unique_within_prokaryote_subtree",
                    )
        rows.append(row)
        if i % 25 == 0:
            print(f"{i + 1}/{len(names)} {name}: {row['lookup_status']} {row['host_domain']}",
                  flush=True)
    target = out / "metadata_host_domains.csv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    counts = {}
    for row in rows:
        key = f"{row['lookup_status']}|{row['host_domain']}"
        counts[key] = counts.get(key, 0) + 1
    print(json.dumps({"names": len(rows), "counts": counts, "out": str(target)}, indent=1))


if __name__ == "__main__":
    main()
