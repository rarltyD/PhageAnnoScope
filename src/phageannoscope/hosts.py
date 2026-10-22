"""Independent, portable host-domain evidence; host relations remain metadata claims.

Only explicit taxids or whole exact scientific names are queried. A verified genus
establishes its taxonomic domain, never a species, strain, or experimental relation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import uuid
from collections import Counter
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from xml.etree import ElementTree as ET

VERSION = "host-evidence-v1"
PLACEHOLDERS = frozenset(
    {
        "",
        "unspecified",
        "unknown",
        "unclassified",
        "none",
        "null",
        "na",
        "n/a",
        "n.a.",
        "not known",
        "not specified",
        "not provided",
        "not available",
        "no data",
        "-",
        "?",
        "uncultured",
        "environmental sample",
    }
)
DOMAINS = ("Bacteria", "Archaea", "Eukaryota")


class _InvalidTaxonomyResponse(ValueError):
    pass


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _immutable(path, data):
    from .acquire import _immutable as save

    save(Path(path), data)


def _relative(path, root):
    path, root = Path(path).resolve(), Path(root).resolve()
    try:
        return path.relative_to(root).as_posix()
    except ValueError as exc:
        raise ValueError(f"Evidence path is outside data root: {path}") from exc


def _resolve(path, root):
    candidate = Path(path)
    if candidate.is_absolute():
        raise ValueError("Portable evidence paths must be relative to data root")
    candidate = (Path(root) / candidate).resolve()
    _relative(candidate, root)
    return candidate


def host_candidates(source_row, source_qualifiers):
    """Keep exact values and every source, including skipped placeholders."""
    candidates = []

    def add(field, value, kind="name", **extra):
        if value is None or str(value) == "":
            return
        value = str(value)
        candidates.append(
            {
                "source_field": field,
                "source_value": value,
                "query_kind": kind,
                "query_value": value.strip(),
                "relation_provenance": "explicit metadata host relation; not experimentally validated",
                **extra,
            }
        )

    add("host_taxid", source_row.get("host_taxid"), "taxid")
    for ordinal, qualifiers in enumerate(source_qualifiers):
        for field in ("host", "lab_host"):
            values = qualifiers.get(field, [])
            if isinstance(values, str):
                values = [values]
            for value in values:
                add("source/" + field, value, source_ordinal=ordinal)
    for field in ("isolation_host_metadata", "metadata_host"):
        add(field, source_row.get(field))
    return candidates


def host_inputs_sha256(source_row, source_qualifiers):
    """Bind evidence to all exact source host values, independent of genome fetch."""
    return _sha(_json(host_candidates(source_row, source_qualifiers)).encode())


def _response_xml(data, expected_root):
    """HTTP 200 alone is not a valid E-utilities biological response."""
    try:
        root = ET.fromstring(data)
    except ET.ParseError as error:
        raise _InvalidTaxonomyResponse("Malformed taxonomy XML") from error
    if root.tag != expected_root or any(
        node.tag.rsplit("}", 1)[-1].casefold() == "error" for node in root.iter()
    ):
        raise _InvalidTaxonomyResponse("Invalid taxonomy service response")
    return root


def _search_taxids(data):
    root = _response_xml(data, "eSearchResult")
    count_text = root.findtext("Count")
    id_list = root.find("IdList")
    if count_text is None or not count_text.isdecimal() or id_list is None:
        raise _InvalidTaxonomyResponse("Taxonomy search lacks Count/IdList")
    count = int(count_text)
    ids = [(node.text or "").strip() for node in id_list.findall("Id")]
    if (
        any(not tid.isdecimal() or int(tid) < 1 for tid in ids)
        or (count == 0 and ids)
        or (count == 1 and len(ids) != 1)
        or (count > 1 and not ids)
    ):
        raise _InvalidTaxonomyResponse("Taxonomy search count and returned identities disagree")
    return count, ids


def _taxonomy_values(data, taxid, name=""):
    taxa = _response_xml(data, "TaxaSet").findall("./Taxon")
    if len(taxa) != 1 or not all(
        taxa[0].findtext(field) for field in ("TaxId", "ScientificName")
    ):
        raise _InvalidTaxonomyResponse("Taxonomy fetch lacks a single complete Taxon")
    if taxa[0].findtext("TaxId") != taxid:
        return {"lookup_status": "identity_mismatch", "reason": "HOST_TAXON_ID_MISMATCH"}
    taxon = taxa[0]
    scientific_name = taxon.findtext("ScientificName") or ""
    if name and scientific_name.casefold() != name.strip().casefold():
        return {"lookup_status": "name_mismatch", "reason": "HOST_SCIENTIFIC_NAME_MISMATCH"}
    lineage = [t.findtext("ScientificName") for t in taxon.findall("./LineageEx/Taxon")]
    lineage.append(scientific_name)
    domains = [domain for domain in DOMAINS if domain in lineage]
    rank = taxon.findtext("Rank") or ""
    domain = domains[0] if len(domains) == 1 else ""
    return {
        "host_taxid": taxid,
        "host_scientific_name": scientific_name,
        "host_taxonomy_rank": rank,
        "host_taxonomy_lineage": lineage,
        "host_taxonomy_lineage_text": taxon.findtext("Lineage") or "",
        "host_domain": domain,
        "host_domain_status": "verified" if domain else "unresolved",
        "lookup_status": "verified" if domain else "domain_unresolved",
        "reason": "" if domain else "HOST_DOMAIN_UNRESOLVED",
        "host_evidence_level": (
            "taxonomy_exact_genus_domain_only"
            if name and rank == "genus"
            else "taxonomy_exact_scientific_name_domain_only"
            if name
            else "taxonomy_explicit_taxid_domain_only"
        ),
    }


def _raw(data, result, label, out, root):
    digest = _sha(data)
    path = out / "taxonomy" / "raw" / (digest + ".xml")
    _immutable(path, data)
    result[label + "_evidence"] = _relative(path, root)
    result[label + "_sha256"] = digest
    result[label + "_historical_path"] = str(path.resolve())


def _verify_candidate(result, root):
    payloads = {}
    for label in ("search", "fetch"):
        path = result.get(label + "_evidence")
        if path:
            raw = _resolve(path, root).read_bytes()
            if _sha(raw) != result.get(label + "_sha256"):
                raise ValueError(f"Taxonomy SHA mismatch: {path}")
            payloads[label] = raw
    # Validate service response structure for cached negative/mismatch outcomes
    # too. Older cache entries may have mistaken a service ERROR for ambiguity.
    if "search" in payloads:
        _search_taxids(payloads["search"])
    if "fetch" in payloads:
        _taxonomy_values(payloads["fetch"], result.get("host_taxid", ""))
    query_tid = result.get("requested_host_taxid", "")
    name = result.get("host_name", "")
    if result.get("query_kind"):
        requested = query_tid if result["query_kind"] == "taxid" else name
        if (
            result.get("query_value") != requested
            or str(result.get("source_value", "")).strip() != requested
        ):
            raise ValueError("Host candidate source value differs from exact query")
    if result.get("host_domain_status") == "verified" and result.get("lookup_status") != "verified":
        raise ValueError("Verified host domain lacks a verified taxonomy lookup")
    if result.get("lookup_status") == "verified":
        tid = result.get("host_taxid", "")
        if query_tid and tid != query_tid:
            raise ValueError("Requested/returned host taxid mismatch")
        if not query_tid:
            if result.get("search_query", {}).get("term") != f'"{name}"[Scientific Name]':
                raise ValueError("Host name evidence query differs from the exact name")
            if "search" not in payloads:
                raise ValueError("Exact-name host evidence lacks search response")
            count, ids = _search_taxids(payloads["search"])
            if count != 1 or ids != [tid]:
                raise ValueError("Host name search was not a unique exact taxid")
        if "fetch" not in payloads:
            raise ValueError("Verified host evidence lacks taxonomy response")
        if result.get("fetch_query", {}).get("id") != tid:
            raise ValueError("Taxonomy fetch query identity mismatch")
        verified = _taxonomy_values(payloads["fetch"], tid, name if not query_tid else "")
        for key in (
            "host_taxid",
            "host_scientific_name",
            "host_taxonomy_rank",
            "host_taxonomy_lineage",
            "host_domain",
            "host_domain_status",
        ):
            if result.get(key) != verified.get(key):
                raise ValueError(f"Taxonomy result mismatch: {key}")
    return result


@lru_cache(maxsize=2)
def _record_batch(path, mtime_ns, size):
    """Bound biological objects to two raw batches while reusing adjacent records."""
    from Bio import SeqIO

    with Path(path).open(encoding="utf-8") as handle:
        return tuple(SeqIO.parse(handle, "genbank"))


def _record_at(path, ordinal):
    stat = Path(path).stat()
    records = _record_batch(str(path), stat.st_mtime_ns, stat.st_size)
    if not 0 <= ordinal < len(records):
        raise ValueError("Host evidence GenBank record ordinal mismatch")
    return records[ordinal]


def resolve_host_taxonomy(
    *, host_taxid=None, host_name=None, out, refresh=False, data_root=None, email=None, api_key=None
):
    """Resolve one full query with independent versioned cache and immutable attempts."""
    from .acquire import _append_json, _request

    out = Path(out).resolve()
    root = Path(data_root).resolve() if data_root else out
    query = {
        "requested_host_taxid": str(host_taxid or "").strip(),
        "host_name": str(host_name or "").strip(),
        "resolver_version": VERSION,
    }
    key = _sha(_json(query).encode())
    index = out / "taxonomy" / "queries" / (key + ".jsonl")
    if not refresh and index.exists():
        latest = json.loads(index.read_text(encoding="utf-8").splitlines()[-1])
        attempt_bytes = _resolve(latest["attempt_path"], out).read_bytes()
        if _sha(attempt_bytes) != latest["attempt_sha256"]:
            raise ValueError("Taxonomy cache attempt SHA mismatch")
        cached = json.loads(attempt_bytes)
        cached_root = (out / latest["data_root_relative_to_cache"]).resolve()
        for label in ("search", "fetch"):
            if cached.get(label + "_evidence"):
                cached[label + "_evidence"] = _relative(
                    _resolve(cached[label + "_evidence"], cached_root), root
                )
        if cached.get("fetch_evidence"):
            cached["host_domain_evidence"] = cached["fetch_evidence"]
        if cached.get("lookup_status") not in {"request_failed", "invalid_response"}:
            try:
                return {**_verify_candidate(cached, root), "cache_hit": True}
            except _InvalidTaxonomyResponse:
                # Raw bytes stay immutable, but a service failure is retryable.
                pass
    result = {
        **query,
        "query_key": key,
        "looked_up_at_utc": _utc(),
        "host_domain": "",
        "host_domain_status": "unresolved",
        "cache_hit": False,
    }
    tid, name = query["requested_host_taxid"], query["host_name"]

    def finish(status, reason):
        result.setdefault("lookup_status", status)
        result.setdefault("reason", reason)
        attempt = out / "taxonomy" / "attempts" / (uuid.uuid4().hex + ".json")
        attempt_bytes = (_json(result) + "\n").encode()
        _immutable(attempt, attempt_bytes)
        _append_json(
            index,
            {
                "attempt_path": _relative(attempt, out),
                "attempt_sha256": _sha(attempt_bytes),
                "data_root_relative_to_cache": os.path.relpath(root, out),
                "lookup_status": result["lookup_status"],
                "at_utc": _utc(),
            },
        )
        return result

    if (not tid and name.casefold() in PLACEHOLDERS) or (
        tid.casefold() in PLACEHOLDERS and not name
    ):
        return finish("placeholder", "HOST_PLACEHOLDER_SKIPPED")
    if tid and (not tid.isdecimal() or int(tid) < 1):
        return finish("invalid_query", "INVALID_HOST_TAXID")
    try:
        if not tid and name:
            params = {
                "db": "taxonomy",
                "term": f'"{name}"[Scientific Name]',
                "retmode": "xml",
                "retmax": 2,
            }
            result["search_query"] = params
            data, request = _request("esearch.fcgi", params, email=email, api_key=api_key)
            result["search_request"] = request
            if data is None:
                return finish("request_failed", "HOST_TAXONOMY_REQUEST_FAILED")
            _raw(data, result, "search", out, root)
            count, ids = _search_taxids(data)
            if count != 1:
                return finish(
                    "negative" if count == 0 else "ambiguous", "HOST_NAME_NOT_UNIQUE_EXACT_TAXON"
                )
            tid = ids[0]
        if not tid:
            return finish("no_query", "NO_EXPLICIT_HOST_TAXON_QUERY")
        params = {"db": "taxonomy", "id": tid, "retmode": "xml"}
        result["fetch_query"] = params
        data, request = _request("efetch.fcgi", params, email=email, api_key=api_key)
        result["fetch_request"] = request
        if data is None:
            return finish("request_failed", "HOST_TAXONOMY_REQUEST_FAILED")
        _raw(data, result, "fetch", out, root)
        result.update(
            _taxonomy_values(data, tid, name if not query["requested_host_taxid"] else "")
        )
        result.update(
            host_domain_evidence=result["fetch_evidence"],
            host_taxonomy_sha256=result["fetch_sha256"],
        )
        return finish("domain_unresolved", "HOST_DOMAIN_UNRESOLVED")
    except (ET.ParseError, UnicodeError, ValueError):
        return finish("invalid_response", "HOST_TAXONOMY_INVALID_RESPONSE")


def aggregate_host_evidence(evidence):
    """A negative candidate does not veto an independently verified one."""
    verified = [x for x in evidence if x.get("host_domain_status") == "verified"]
    domains = sorted({x["host_domain"] for x in verified})
    result = {
        "host_domain": "",
        "host_domain_status": "unresolved",
        "host_taxonomy_evidence": evidence,
        "host_domains_verified": domains,
        "host_reason": "HOST_UNRESOLVED",
    }
    if len(domains) > 1:
        result.update(host_domain_status="conflict", host_reason="HOST_DOMAIN_CONFLICT")
    elif len(domains) == 1:
        result.update({key: value for key, value in verified[0].items() if key.startswith("host_")})
        result.update(host_domain_status="verified", host_reason="")
    elif not evidence:
        result["host_reason"] = "NO_EXPLICIT_HOST_TAXON_QUERY"
    return result


def resolve_record_hosts(
    source_row, source_qualifiers, *, out, refresh=False, data_root=None, email=None, api_key=None
):
    candidates = host_candidates(source_row, source_qualifiers)
    evidence = []
    for candidate in candidates:
        query = {
            "host_taxid" if candidate["query_kind"] == "taxid" else "host_name": candidate[
                "query_value"
            ]
        }
        resolved = resolve_host_taxonomy(
            **query, out=out, refresh=refresh, data_root=data_root, email=email, api_key=api_key
        )
        evidence.append({**resolved, **candidate})
    result = aggregate_host_evidence(evidence)
    result["host_inputs_sha256"] = _sha(_json(candidates).encode())
    return result


def resolve_hosts(manifest, genbank_dir, out, refresh=False, data_root=None, include_records=True):
    """Update host evidence from existing GenBank only; never downloads sequences."""
    from .acquire import _append_json, _targets
    from .manifest import _raw_records

    out, genbank_dir = Path(out).resolve(), Path(genbank_dir).resolve()
    if out == Path(manifest).resolve():
        raise ValueError("Host resolution must not overwrite its input manifest")
    root = (
        Path(data_root).resolve()
        if data_root
        else Path(os.path.commonpath([out.parent, genbank_dir]))
    )
    by_version, parse_errors = _raw_records(genbank_dir)
    by_base = {}
    for identity in by_version:
        by_base.setdefault(identity.rsplit(".", 1)[0], []).append(identity)
    run_id = "hosts_" + uuid.uuid4().hex
    rows = []
    counts = Counter()
    fields = [
        "run_id",
        "resolver_version",
        "requested_accession",
        "accession_version",
        "host_domain",
        "host_domain_status",
        "host_reason",
        "evaluated_at_utc",
        "genbank_path",
        "genbank_historical_path",
        "genbank_sha256",
        "sequence_sha256",
        "record_ordinal",
        "host_taxonomy_evidence_json",
        "host_domains_verified",
        "host_name",
        "host_taxid",
        "host_scientific_name",
        "host_taxonomy_rank",
        "host_taxonomy_lineage",
        "host_taxonomy_lineage_text",
        "host_evidence_level",
        "host_domain_evidence",
        "host_taxonomy_sha256",
        "host_inputs_sha256",
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = out.with_name(out.name + ".tmp")
    row_events = Path(str(out) + ".rows.jsonl")
    handle = temporary.open("w", encoding="utf-8", newline="")
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    emitted = 0

    def emit(row):
        nonlocal emitted
        writer.writerow(row)
        _append_json(row_events, row)
        emitted += 1
        counts[row["host_domain_status"]] += 1
        if include_records:
            rows.append(row)
        if emitted % 100 == 0:
            from .util import write_json

            handle.flush()
            write_json(
                str(out) + ".progress.json",
                {
                    "run_id": run_id,
                    "processed_n": emitted,
                    "counts": dict(counts),
                    "at_utc": _utc(),
                },
            )

    for requested, source in _targets(manifest):
        versions = (
            [requested]
            if requested in by_version
            else by_base.get(requested, [])
            if "." not in requested
            else []
        )
        row = {
            "run_id": run_id,
            "resolver_version": VERSION,
            "requested_accession": requested,
            "accession_version": "",
            "host_domain": "",
            "host_domain_status": "unresolved",
            "host_reason": "GENBANK_NOT_AVAILABLE",
            "evaluated_at_utc": _utc(),
        }
        if len(versions) != 1:
            if len(versions) > 1:
                row["host_reason"] = "GENBANK_IDENTITY_AMBIGUOUS"
            emit(row)
            continue
        version = versions[0]
        candidates = by_version[version]
        row["accession_version"] = version
        if len({(x["sequence_sha"], x["raw_sha"]) for x in candidates}) != 1:
            row["host_reason"] = "MULTIPLE_RAW_SNAPSHOTS_REQUIRE_RESOLUTION"
            emit(row)
            continue
        item = candidates[0]
        record = _record_at(item["path"], item["record_ordinal"])
        if record.id != version or "." not in version:
            row["host_reason"] = "GENBANK_EXPLICIT_VERSION_REQUIRED"
            emit(row)
            continue
        row.update(
            genbank_path=_relative(item["path"], root),
            genbank_historical_path=item["path"],
            genbank_sha256=item["raw_sha"],
            sequence_sha256=item["sequence_sha"],
            record_ordinal=item["record_ordinal"],
        )
        evidence = resolve_record_hosts(
            source,
            [dict(f.qualifiers) for f in record.features if f.type == "source"],
            out=genbank_dir,
            refresh=refresh,
            data_root=root,
            email=os.environ.get("NCBI_EMAIL"),
            api_key=os.environ.get("NCBI_API_KEY"),
        )
        row.update({k: _json(v) if isinstance(v, (list, dict)) else v for k, v in evidence.items()})
        # Explicit JSON name keeps CSV consumers from interpreting repr strings.
        row["host_taxonomy_evidence_json"] = row.pop("host_taxonomy_evidence")
        emit(row)
    handle.close()
    temporary.replace(out)
    from .util import file_sha

    meta = {
        "resolver_version": VERSION,
        "data_root_relative_to_output": os.path.relpath(root, out.parent),
        "historical_data_root": str(root),
        "csv_sha256": file_sha(out),
        "run_id": run_id,
        "evaluated_at_utc": _utc(),
        "parse_errors": parse_errors,
    }
    Path(str(out) + ".meta.json").write_text(_json(meta) + "\n", encoding="utf-8")
    events = Path(str(out) + ".events.jsonl")
    _append_json(
        events,
        {
            **meta,
            "refresh": refresh,
            "records": rows,
            "record_count": emitted,
            "row_events": str(row_events),
        },
    )
    return {
        "run_id": run_id,
        "requested": emitted,
        "records": rows,
        "counts": dict(counts),
        "out": str(out),
        "events": str(events),
        "data_root": str(root),
    }


def load_host_evidence(path, data_root=None):
    """Verify portable raw hashes, taxonomy assertions, and exact input identities."""
    path = Path(path).resolve()
    from .util import file_sha

    meta = json.loads(Path(str(path) + ".meta.json").read_text(encoding="utf-8"))
    if meta.get("resolver_version") != VERSION or meta.get("csv_sha256") != file_sha(path):
        raise ValueError("Host evidence CSV/version checksum mismatch")
    root = (
        Path(data_root).resolve()
        if data_root
        else (path.parent / meta["data_root_relative_to_output"]).resolve()
    )
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    mapping = {}
    checked_genbank = {}
    for row in rows:
        version = row.get("accession_version")
        if not version:
            continue
        if version in mapping:
            raise ValueError("Duplicate accession.version in host evidence")
        if row.get("genbank_path"):
            gb = _resolve(row["genbank_path"], root)
            if gb not in checked_genbank:
                with gb.open("rb") as handle:
                    checked_genbank[gb] = hashlib.file_digest(handle, "sha256").hexdigest()
            if checked_genbank[gb] != row["genbank_sha256"]:
                raise ValueError("Host evidence GenBank SHA mismatch")
            ordinal = int(row["record_ordinal"])
            record = _record_at(gb, ordinal)
            if (
                record.id != version
                or _sha(str(record.seq).upper().encode("ascii")) != row["sequence_sha256"]
            ):
                raise ValueError("Host evidence record identity/sequence SHA mismatch")
        evidence = json.loads(row.get("host_taxonomy_evidence_json") or "[]")
        for candidate in evidence:
            _verify_candidate(candidate, root)
        if evidence:
            aggregate = aggregate_host_evidence(evidence)
            for key in ("host_domain", "host_domain_status"):
                if row.get(key) != aggregate[key]:
                    raise ValueError("Host evidence aggregate contradicts candidate evidence")
            row.update(aggregate)
        row["host_evidence_data_root"] = str(root)
        row["host_evidence_verified"] = True
        mapping[version] = row
    return mapping
