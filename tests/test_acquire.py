"""Network-free tests. All sequences below are explicitly synthetic test data."""

import csv
import hashlib
import io
import json
from urllib.error import HTTPError

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phageannoscope import acquire


def synthetic_gb(*ids):
    records = [
        SeqRecord(
            Seq("ATGAAATAA"),
            id=identity,
            name=identity.split(".")[0],
            description="SYNTHETIC TEST ONLY",
            annotations={"molecule_type": "DNA"},
        )
        for identity in ids
    ]
    handle = io.StringIO()
    SeqIO.write(records, handle, "genbank")
    return handle.getvalue().encode()


def response(data):
    return data, {
        "request_url": "https://example.invalid/efetch",
        "started_at_utc": "t0",
        "ended_at_utc": "t1",
        "http_status": 200,
        "attempts": 1,
    }


def test_fetch_exact_versions_raw_bytes_and_resume(tmp_path, monkeypatch):
    raw = synthetic_gb("TEST001.4", "TEST002.8")
    calls = []
    monkeypatch.setattr(
        acquire, "_request", lambda *a, **kw: calls.append((a, kw)) or response(raw)
    )
    result = acquire.fetch_ncbi(
        accessions=["TEST001", "TEST002.8"],
        out=tmp_path / "ncbi",
        events=tmp_path / "events.csv",
        resolve_hosts=False,
    )
    assert result["counts"] == {"success": 2}
    assert result["records"][0]["accession_version"] == "TEST001.4"
    assert result["records"][1]["accession_version"] == "TEST002.8"
    saved = next((tmp_path / "ncbi" / "raw").glob("*.gb"))
    assert saved.read_bytes() == raw
    assert saved.stem == hashlib.sha256(raw).hexdigest()
    again = acquire.fetch_ncbi(
        accessions=["TEST001", "TEST002.8"],
        out=tmp_path / "ncbi",
        events=tmp_path / "events.csv",
        resolve_hosts=False,
    )
    assert again["counts"] == {"cache_hit": 2}
    assert len(calls) == 1
    with (tmp_path / "events.csv").open() as handle:
        assert len(list(csv.DictReader(handle))) == 4


def test_no_version_substitution_and_unexpected_records_quarantined(tmp_path, monkeypatch):
    monkeypatch.setattr(
        acquire, "_request", lambda *a, **kw: response(synthetic_gb("TEST001.2", "UNREQUESTED.9"))
    )
    result = acquire.fetch_ncbi(
        accessions=["TEST001.1"],
        out=tmp_path / "ncbi",
        events=tmp_path / "events.csv",
        resolve_hosts=False,
    )
    assert result["records"][0]["status"] != "success"
    assert (tmp_path / "ncbi" / "quarantine.jsonl").exists()


def test_ambiguous_duplicate_record_never_success(tmp_path, monkeypatch):
    monkeypatch.setattr(
        acquire, "_request", lambda *a, **kw: response(synthetic_gb("TEST001.1", "TEST001.1"))
    )
    result = acquire.fetch_ncbi(
        accessions=["TEST001"],
        out=tmp_path / "ncbi",
        events=tmp_path / "events.csv",
        resolve_hosts=False,
    )
    assert result["counts"] == {"ambiguous": 1}


def test_repeated_fresh_request_preserves_raw_with_new_request_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(acquire, "_request", lambda *a, **kw: response(synthetic_gb("TEST001.3")))
    for _ in range(2):
        result = acquire.fetch_ncbi(
            accessions=["TEST001"],
            out=tmp_path / "ncbi",
            events=tmp_path / "events.csv",
            resume=False,
            resolve_hosts=False,
        )
        assert result["counts"] == {"success": 1}
    assert len(list((tmp_path / "ncbi" / "raw").glob("*.gb"))) == 1
    assert len(list((tmp_path / "ncbi" / "requests").glob("*.json"))) == 2


def test_bounded_retry_and_credentials_redacted(monkeypatch):
    calls = []
    monkeypatch.setattr(acquire, "_rate_limit", lambda *a: None)
    monkeypatch.setattr(acquire.time, "sleep", lambda *a: None)

    def fail(request, **kwargs):
        calls.append(request.full_url)
        raise HTTPError(
            request.full_url, 429, "api_key=SECRET email=private", {"Retry-After": "0"}, None
        )

    monkeypatch.setattr(acquire, "urlopen", fail)
    data, meta = acquire._request(
        "efetch.fcgi",
        {"db": "nuccore", "id": "TEST001"},
        email="private",
        api_key="SECRET",
        max_attempts=3,
    )
    assert data is None and len(calls) == 3
    assert "SECRET" not in json.dumps(meta) and "private" not in json.dumps(meta)
    assert meta["attempts"] == 3


def test_taxonomy_requires_exact_unique_scientific_name(tmp_path, monkeypatch):
    xml = b"<eSearchResult><Count>2</Count><IdList><Id>1</Id><Id>2</Id></IdList></eSearchResult>"
    monkeypatch.setattr(acquire, "_request", lambda *a, **kw: response(xml))
    result = acquire.resolve_host_taxonomy(host_name="Synthetic ambiguous host", out=tmp_path)
    assert result["host_domain"] == ""
    assert result["reason"] == "HOST_NAME_NOT_UNIQUE_EXACT_TAXON"


def test_taxonomy_domain_from_preserved_lineage(tmp_path, monkeypatch):
    search = b"<eSearchResult><Count>1</Count><IdList><Id>42</Id></IdList></eSearchResult>"
    taxonomy = b"<TaxaSet><Taxon><TaxId>42</TaxId><ScientificName>Synthetic test host</ScientificName><LineageEx><Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName></Taxon></LineageEx></Taxon></TaxaSet>"
    responses = iter([search, taxonomy])
    monkeypatch.setattr(acquire, "_request", lambda *a, **kw: response(next(responses)))
    result = acquire.resolve_host_taxonomy(host_name="Synthetic test host", out=tmp_path)
    assert result["host_domain"] == "Bacteria" and result["host_domain_status"] == "verified"
    assert result["host_taxid"] == "42"
    assert hashlib.sha256(taxonomy).hexdigest() == result["host_taxonomy_sha256"]


def test_invalid_accession_rejected_before_network(tmp_path):
    with pytest.raises(ValueError, match="Invalid accession"):
        acquire.fetch_ncbi(
            accessions=["A&api_key=bad"], out=tmp_path, events=tmp_path / "events.csv"
        )
