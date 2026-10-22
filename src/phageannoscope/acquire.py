"""Bounded NCBI retrieval with content-addressed raw responses and explicit identities.

No accession version is invented. The saved GenBank bytes are the response bytes,
including batch responses; the index identifies each record within that snapshot.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import os
import re
import time
import uuid
from collections import Counter
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ACCESSION = re.compile(r"^[A-Za-z][A-Za-z0-9_]*(?:\.[1-9][0-9]*)?$")
EVENT_FIELDS = [
    "run_id",
    "stage",
    "requested_accession",
    "returned_accession_version",
    "status",
    "started_at_utc",
    "ended_at_utc",
    "request_url",
    "http_status",
    "attempts",
    "genbank_path",
    "genbank_sha256",
    "sequence_sha256",
    "message",
]
_last_request = 0.0
LOG = logging.getLogger(__name__)


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _append_json(path: Path, row: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def _immutable(path: Path, data: bytes):
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as handle:
            handle.write(data)
    except FileExistsError:
        if path.read_bytes() != data:
            raise ValueError(f"Immutable raw path has different content: {path}")


def _rate_limit(api_key=None):
    # Conservative limit also applies with a key. Every retry/search shares it.
    global _last_request
    delay = 0.34 - (time.monotonic() - _last_request)
    if delay > 0:
        time.sleep(delay)
    _last_request = time.monotonic()


def _request(endpoint, params, *, email=None, api_key=None, max_attempts=4):
    """Return bytes and safe request metadata; never expose credential values."""
    safe = dict(params, tool="PhageAnnoScope")
    actual = dict(safe)
    if email:
        actual["email"] = email
    if api_key:
        actual["api_key"] = api_key
    url = EUTILS + endpoint + "?" + urlencode(actual)
    meta = {
        "request_url": EUTILS + endpoint + "?" + urlencode(safe),
        "started_at_utc": _utc(),
        "http_status": None,
        "attempts": 0,
    }
    for attempt in range(max_attempts):
        _rate_limit(api_key)
        meta["attempts"] = attempt + 1
        try:
            with urlopen(
                Request(url, headers={"User-Agent": "PhageAnnoScope/0.1"}), timeout=60
            ) as response:
                data = response.read()
                meta["http_status"] = response.status
            meta["ended_at_utc"] = _utc()
            return data, meta
        except HTTPError as exc:
            meta["http_status"] = exc.code
            if exc.code not in {429, 500, 502, 503, 504} or attempt + 1 == max_attempts:
                break
            value = exc.headers.get("Retry-After", "") if exc.headers else ""
            try:
                delay = min(30.0, max(0, float(value)))
            except ValueError:
                try:
                    delay = min(
                        30.0,
                        max(
                            0,
                            (
                                parsedate_to_datetime(value) - datetime.now(timezone.utc)
                            ).total_seconds(),
                        ),
                    )
                except (TypeError, ValueError, OverflowError):
                    delay = min(2**attempt, 15)
            time.sleep(delay)
        except (URLError, TimeoutError, OSError):
            if attempt + 1 == max_attempts:
                break
            time.sleep(min(2**attempt, 15))
    meta["ended_at_utc"] = _utc()
    meta["error"] = "NCBI request failed after bounded retries"
    return None, meta


def _parse_genbank(data):
    from Bio import SeqIO

    return list(SeqIO.parse(io.StringIO(data.decode("utf-8")), "genbank"))


def _identity(record):
    rid = str(record.id)
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.[1-9][0-9]*", rid):
        raise ValueError("NCBI record has no explicit accession.version")
    sequence = str(record.seq).upper()
    if not sequence:
        raise ValueError("NCBI record has no sequence")
    return {
        "accession_version": rid,
        "accession": rid.rsplit(".", 1)[0],
        "sequence_sha256": _sha(sequence.encode("ascii")),
        "length_bp": len(sequence),
        "topology": record.annotations.get("topology", "unknown"),
        "organism": record.annotations.get("organism"),
        "source_qualifiers": [dict(f.qualifiers) for f in record.features if f.type == "source"],
    }


def _targets(manifest=None, accessions=None):
    rows = []
    if manifest:
        with Path(manifest).open(encoding="utf-8-sig", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    if accessions:
        if isinstance(accessions, (str, Path)) and Path(accessions).exists():
            values = Path(accessions).read_text(encoding="utf-8-sig").splitlines()
        elif isinstance(accessions, str):
            values = accessions.split(",")
        else:
            values = list(accessions)
        wanted = [v.strip() for v in values if v.strip() and not v.lstrip().startswith("#")]
        if rows:
            by_id = {
                identity: r
                for r in rows
                for identity in (r.get("accession_version"), r.get("accession"))
                if identity
            }
            rows = [by_id.get(v, {"accession": v}) for v in wanted]
        else:
            rows = [{"accession": v} for v in wanted]
    if not rows:
        raise ValueError("A nonempty frozen manifest or explicit accession list is required")
    result = []
    seen = set()
    for row in rows:
        requested = (row.get("accession_version") or row.get("accession") or "").strip()
        if not ACCESSION.fullmatch(requested):
            raise ValueError(f"Invalid accession: {requested!r}")
        if requested not in seen:
            seen.add(requested)
            result.append((requested, row))
    return result


def resolve_host_taxonomy(
    *, host_taxid=None, host_name=None, out, email=None, api_key=None, refresh=False, data_root=None
):
    """Compatibility entry point for the independent versioned host resolver."""
    from .hosts import resolve_host_taxonomy as resolve

    return resolve(
        host_taxid=host_taxid,
        host_name=host_name,
        out=out,
        email=email,
        api_key=api_key,
        refresh=refresh,
        data_root=data_root,
    )


def fetch_ncbi(
    manifest=None,
    accessions=None,
    *,
    out,
    events,
    resume=True,
    limit=None,
    email=None,
    api_key=None,
    resolve_hosts=True,
    batch_size=20,
    progress_path=None,
):
    out, events = Path(out), Path(events)
    out.mkdir(parents=True, exist_ok=True)
    events.parent.mkdir(parents=True, exist_ok=True)
    email = email or os.environ.get("NCBI_EMAIL")
    api_key = api_key or os.environ.get("NCBI_API_KEY")
    targets = _targets(manifest, accessions)
    if not 1 <= batch_size <= 100:
        raise ValueError("batch_size must be between 1 and 100")
    if limit is not None:
        if limit < 0:
            raise ValueError("limit must be nonnegative")
        targets = targets[:limit]
    index = out / "fetch_records.jsonl"
    cached = {}
    checked_raw = {}
    if resume and index.exists():
        for line in index.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            if row.get("status") == "success":
                path = Path(row["genbank_path"])
                if path not in checked_raw:
                    checked_raw[path] = _sha(path.read_bytes()) if path.exists() else None
                if checked_raw[path] == row["genbank_sha256"]:
                    cached[row["requested_accession"]] = row
    run_id = "fetch_" + uuid.uuid4().hex
    output = []
    wall_start = time.perf_counter()

    def record_event(row):
        row.setdefault("run_id", run_id)
        row.setdefault("stage", "ncbi_fetch")
        empty = not events.exists() or events.stat().st_size == 0
        with events.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=EVENT_FIELDS, extrasaction="ignore")
            if empty:
                writer.writeheader()
            writer.writerow(row)
        output.append(row)

    pending = []
    for requested, source_row in targets:
        if requested in cached:
            record_event({**cached[requested], "status": "cache_hit", "run_id": run_id})
        else:
            pending.append((requested, source_row))
    for offset in range(0, len(pending), batch_size):
        batch = pending[offset : offset + batch_size]
        data, request = _request(
            "efetch.fcgi",
            {
                "db": "nuccore",
                "id": ",".join(x[0] for x in batch),
                "rettype": "gbwithparts",
                "retmode": "text",
            },
            email=email,
            api_key=api_key,
        )
        if data is None:
            for requested, _ in batch:
                row = {**request, "requested_accession": requested, "status": "request_failed"}
                _append_json(index, row)
                record_event(row)
            continue
        digest = _sha(data)
        path = out / "raw" / (digest + ".gb")
        _immutable(path, data)
        request.update(genbank_path=str(path.resolve()), genbank_sha256=digest)
        _immutable(
            out / "requests" / (digest + "_" + run_id + "_" + str(offset) + ".json"),
            (json.dumps(request, ensure_ascii=False, indent=2) + "\n").encode(),
        )
        try:
            records = _parse_genbank(data)
            identities = [_identity(r) for r in records]
        except Exception as exc:
            records, identities = [], []
            request["message"] = f"GenBank parse failed: {type(exc).__name__}"
        assigned = set()
        for requested, source_row in batch:
            matches = [
                i
                for i, ident in enumerate(identities)
                if ident["accession_version"] == requested
                or ("." not in requested and ident["accession"] == requested)
            ]
            row = {**request, "requested_accession": requested}
            if len(matches) != 1:
                different_versions = [
                    x["accession_version"]
                    for x in identities
                    if x["accession"] == requested.rsplit(".", 1)[0]
                ]
                row["status"] = (
                    "ambiguous"
                    if len(matches) > 1
                    else "version_mismatch"
                    if "." in requested and different_versions
                    else "missing"
                )
                row["identity_resolution_status"] = (
                    "unresolved_not_returned" if not matches else "ambiguous"
                )
                row["returned_candidates"] = different_versions
                row["message"] = (
                    "Exact requested identity was not uniquely returned; no substitute accepted"
                )
            else:
                i = matches[0]
                assigned.add(i)
                row.update(
                    identities[i],
                    status="success",
                    returned_accession_version=identities[i]["accession_version"],
                )
                if resolve_hosts:
                    from .hosts import resolve_record_hosts

                    row.update(
                        resolve_record_hosts(
                            source_row,
                            row["source_qualifiers"],
                            out=out,
                            email=email,
                            api_key=api_key,
                        )
                    )
            _append_json(index, row)
            record_event(row)
        if len(assigned) != len(records):
            _append_json(
                out / "quarantine.jsonl",
                {
                    **request,
                    "run_id": run_id,
                    "stage": "ncbi_fetch",
                    "status": "unexpected_or_duplicate_records",
                    "returned_ids": [x.get("accession_version") for x in identities],
                    "requested_ids": [x[0] for x in batch],
                },
            )
        progress = {
            "run_id": run_id,
            "planned_n": len(targets),
            "attempted_n": len(output),
            "counts": dict(Counter(r["status"] for r in output)),
            "wall_seconds": time.perf_counter() - wall_start,
            "batch_size": batch_size,
            "last_batch_n": len(batch),
            "at_utc": _utc(),
        }
        if progress_path:
            from .util import write_json

            write_json(progress_path, progress)
        LOG.info("NCBI retrieval %d/%d; %s", len(output), len(targets), progress["counts"])
    summary = {
        "run_id": run_id,
        "requested": len(targets),
        "counts": dict(Counter(r["status"] for r in output)),
        "records": output,
        "events": str(events.resolve()),
        "index": str(index.resolve()),
        "wall_seconds": time.perf_counter() - wall_start,
        "batch_size": batch_size,
    }
    return summary
