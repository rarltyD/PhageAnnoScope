"""Network-free host regression tests; non-fixture taxa/sequences are synthetic."""

import csv
import hashlib
import shutil
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from phageannoscope import acquire, hosts

REAL = Path(__file__).parent / "fixtures" / "real" / "host_taxonomy"
SEARCH = "754a9d5a62c3f39cc1a3d7c1feed820ca6da9b6630e62a268f7037aade33f8a2.xml"
TAXONOMY = "9f0b9c5c5028969e93146e577cfaab4eb45bcd837bedda5842e422dd369e7b80.xml"


def response(data):
    return data, {
        "request_url": "https://example.invalid/synthetic-test",
        "http_status": 200,
        "started_at_utc": "2026-10-16T00:00:00+00:00",
        "ended_at_utc": "2026-10-16T00:00:01+00:00",
        "attempts": 1,
    }


def search(tid="42", count=1):
    ids = f"<Id>{tid}</Id>" if count else ""
    return f"<eSearchResult><Count>{count}</Count><IdList>{ids}</IdList></eSearchResult>".encode()


def taxonomy(name="Synthetic genus", tid="42", domain="Bacteria", rank="genus"):
    return (
        f"<TaxaSet><Taxon><TaxId>{tid}</TaxId><ScientificName>{name}</ScientificName>"
        f"<Rank>{rank}</Rank><Lineage>{domain}</Lineage><LineageEx><Taxon>"
        f"<TaxId>2</TaxId><ScientificName>{domain}</ScientificName></Taxon></LineageEx>"
        "</Taxon></TaxaSet>"
    ).encode()


def responder(monkeypatch, payloads):
    calls = []
    iterator = iter(payloads)

    def request(endpoint, params, **kwargs):
        assert params["db"] == "taxonomy", "Host refresh must never fetch nuccore"
        calls.append((endpoint, params))
        return response(next(iterator))

    monkeypatch.setattr(acquire, "_request", request)
    return calls


def prepare(tmp_path, rows=None, *, source=None):
    raw = tmp_path / "data" / "ncbi"
    raw.mkdir(parents=True)
    record = SeqRecord(
        Seq("ATGAAATAA"),
        id="SYNTHETIC001.7",
        name="SYNTHETIC001",
        description="SYNTHETIC TEST ONLY",
        annotations={"molecule_type": "DNA"},
    )
    record.features = [SeqFeature(SimpleLocation(0, 9), type="source", qualifiers=source or {})]
    SeqIO.write(record, raw / "record.gb", "genbank")
    manifest = tmp_path / "manifest.csv"
    rows = (
        rows
        if rows is not None
        else [{"accession": "SYNTHETIC001", "metadata_host": "Synthetic genus"}]
    )
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)
    return manifest, raw, tmp_path / "data_manifest" / "hosts.csv"


def test_real_escherichia_exact_genus_domain_only(tmp_path, monkeypatch):
    calls = responder(monkeypatch, [(REAL / SEARCH).read_bytes(), (REAL / TAXONOMY).read_bytes()])
    result = hosts.resolve_record_hosts(
        {"isolation_host_metadata": "Unspecified", "metadata_host": "Escherichia"}, [], out=tmp_path
    )
    assert result["host_domain_status"] == "verified"
    assert result["host_taxid"] == "561" and result["host_domain"] == "Bacteria"
    assert result["host_taxonomy_rank"] == "genus"
    assert result["host_evidence_level"] == "taxonomy_exact_genus_domain_only"
    evidence = result["host_taxonomy_evidence"]
    assert [c["source_field"] for c in evidence] == ["isolation_host_metadata", "metadata_host"]
    assert evidence[0]["lookup_status"] == "placeholder" and len(calls) == 2
    assert calls[0][1]["term"] == '"Escherichia"[Scientific Name]'
    assert (tmp_path / result["host_domain_evidence"]).read_bytes() == (
        REAL / TAXONOMY
    ).read_bytes()
    assert "not experimentally validated" in evidence[1]["relation_provenance"]


@pytest.mark.parametrize("domain", ["Bacteria", "Archaea", "Eukaryota"])
def test_domains_from_lineage_not_name(tmp_path, monkeypatch, domain):
    responder(monkeypatch, [search(), taxonomy(domain=domain)])
    result = hosts.resolve_record_hosts({"metadata_host": "Synthetic genus"}, [], out=tmp_path)
    assert result["host_domain"] == domain


@pytest.mark.parametrize("placeholder", ["Unspecified", "UNKNOWN", "unclassified", "N/A", "-", ""])
def test_placeholders_never_query_ncbi(tmp_path, monkeypatch, placeholder):
    monkeypatch.setattr(
        acquire, "_request", lambda *a, **k: pytest.fail("placeholder network call")
    )
    result = hosts.resolve_record_hosts({"metadata_host": placeholder}, [], out=tmp_path)
    assert result["host_domain_status"] == "unresolved"


@pytest.mark.parametrize("count,status", [(0, "negative"), (2, "ambiguous")])
def test_incomplete_ambiguous_name_is_not_split(tmp_path, monkeypatch, count, status):
    calls = responder(monkeypatch, [search(count=count)])
    result = hosts.resolve_record_hosts(
        {"metadata_host": "Synthetic genus strain strange words"}, [], out=tmp_path
    )
    assert result["host_domain_status"] == "unresolved"
    assert result["host_taxonomy_evidence"][0]["lookup_status"] == status
    assert calls[0][1]["term"] == '"Synthetic genus strain strange words"[Scientific Name]'
    assert len(calls) == 1


def test_taxid_response_must_match_requested_identity(tmp_path, monkeypatch):
    responder(monkeypatch, [taxonomy(tid="43")])
    result = hosts.resolve_record_hosts({"host_taxid": "42"}, [], out=tmp_path)
    assert result["host_domain_status"] == "unresolved"
    assert result["host_taxonomy_evidence"][0]["reason"] == "HOST_TAXON_ID_MISMATCH"


def test_exact_scientific_name_must_match(tmp_path, monkeypatch):
    responder(monkeypatch, [search(), taxonomy(name="Synthetic different")])
    result = hosts.resolve_record_hosts({"metadata_host": "Synthetic genus"}, [], out=tmp_path)
    assert result["host_domain_status"] == "unresolved"
    assert result["host_taxonomy_evidence"][0]["reason"] == "HOST_SCIENTIFIC_NAME_MISMATCH"


def test_preserve_conflicting_domains_from_all_sources(tmp_path, monkeypatch):
    responder(
        monkeypatch,
        [
            taxonomy(domain="Bacteria"),
            search("43"),
            taxonomy(name="Synthetic archaeon", tid="43", domain="Archaea"),
        ],
    )
    result = hosts.resolve_record_hosts(
        {"host_taxid": "42", "metadata_host": "Synthetic archaeon"}, [], out=tmp_path
    )
    assert result["host_domain_status"] == "conflict" and result["host_domain"] == ""
    assert result["host_domains_verified"] == ["Archaea", "Bacteria"]
    assert len(result["host_taxonomy_evidence"]) == 2


def test_unresolved_source_does_not_veto_other_verified_source(tmp_path, monkeypatch):
    responder(monkeypatch, [search(count=0), search(), taxonomy()])
    result = hosts.resolve_record_hosts(
        {"metadata_host": "Synthetic genus"}, [{"host": ["Synthetic absent"]}], out=tmp_path
    )
    assert result["host_domain_status"] == "verified"
    assert len(result["host_taxonomy_evidence"]) == 2


def test_missing_host_fields_no_organism_guess(tmp_path, monkeypatch):
    monkeypatch.setattr(acquire, "_request", lambda *a, **k: pytest.fail("invented host query"))
    result = hosts.resolve_record_hosts(
        {"organism": "Escherichia phage synthetic"}, [], out=tmp_path
    )
    assert result["host_domain_status"] == "unresolved" and not result["host_taxonomy_evidence"]


def test_refresh_negative_and_failure_automatically_retried(tmp_path, monkeypatch):
    calls = responder(monkeypatch, [search(count=0), search(), taxonomy()])
    first = hosts.resolve_host_taxonomy(host_name="Synthetic genus", out=tmp_path)
    assert first["lookup_status"] == "negative"
    cached = hosts.resolve_host_taxonomy(host_name="Synthetic genus", out=tmp_path)
    assert cached["cache_hit"] and len(calls) == 1
    refreshed = hosts.resolve_host_taxonomy(host_name="Synthetic genus", out=tmp_path, refresh=True)
    assert refreshed["host_domain_status"] == "verified" and len(calls) == 3
    assert len(list((tmp_path / "taxonomy" / "attempts").glob("*.json"))) == 2

    responses = iter([(None, {"http_status": 503, "attempts": 4}), response(search(count=0))])
    monkeypatch.setattr(acquire, "_request", lambda *a, **k: next(responses))
    failure = hosts.resolve_host_taxonomy(host_name="Synthetic network failure", out=tmp_path)
    assert failure["lookup_status"] == "request_failed"
    retry = hosts.resolve_host_taxonomy(host_name="Synthetic network failure", out=tmp_path)
    assert retry["lookup_status"] == "negative" and not retry["cache_hit"]


@pytest.mark.parametrize(
    "payload",
    [
        b"<eSearchResult><ERROR>temporary backend failure</ERROR></eSearchResult>",
        b"<ERROR>temporary backend failure</ERROR>",
        b"<eSearchResult><IdList/></eSearchResult>",
        b"<eSearchResult><Count>0</Count></eSearchResult>",
        b"<eSearchResult><Count>1</Count><IdList/></eSearchResult>",
    ],
)
def test_http_200_search_error_or_missing_structure_is_retried(tmp_path, monkeypatch, payload):
    calls = responder(monkeypatch, [payload, search(), taxonomy()])
    first = hosts.resolve_host_taxonomy(host_name="Synthetic genus", out=tmp_path)
    assert first["lookup_status"] == "invalid_response"
    assert first["search_request"]["http_status"] == 200
    assert (tmp_path / first["search_evidence"]).read_bytes() == payload
    retry = hosts.resolve_host_taxonomy(host_name="Synthetic genus", out=tmp_path)
    assert retry["host_domain_status"] == "verified" and not retry["cache_hit"]
    assert len(calls) == 3
    assert len(list((tmp_path / "taxonomy" / "attempts").glob("*.json"))) == 2


@pytest.mark.parametrize(
    "payload",
    [
        b"<TaxaSet><ERROR>temporary backend failure</ERROR></TaxaSet>",
        b"<ERROR>temporary backend failure</ERROR>",
        b"<TaxaSet/>",
        b"<TaxaSet><Taxon><TaxId>42</TaxId></Taxon></TaxaSet>",
    ],
)
def test_http_200_fetch_error_or_missing_structure_is_retried(tmp_path, monkeypatch, payload):
    calls = responder(monkeypatch, [payload, taxonomy()])
    first = hosts.resolve_host_taxonomy(host_taxid="42", out=tmp_path)
    assert first["lookup_status"] == "invalid_response"
    assert first["fetch_request"]["http_status"] == 200
    assert (tmp_path / first["fetch_evidence"]).read_bytes() == payload
    retry = hosts.resolve_host_taxonomy(host_taxid="42", out=tmp_path)
    assert retry["host_domain_status"] == "verified" and not retry["cache_hit"]
    assert len(calls) == 2


def test_standalone_refresh_never_downloads_genbank_and_binds_sequence(tmp_path, monkeypatch):
    manifest, raw, out = prepare(tmp_path)
    calls = responder(monkeypatch, [search(), taxonomy()])
    before = (raw / "record.gb").read_bytes()
    result = hosts.resolve_hosts(manifest, raw, out)
    assert result["counts"] == {"verified": 1} and len(calls) == 2
    assert (raw / "record.gb").read_bytes() == before
    loaded = hosts.load_host_evidence(out)
    evidence = loaded["SYNTHETIC001.7"]
    assert evidence["sequence_sha256"] == hashlib.sha256(b"ATGAAATAA").hexdigest()
    assert evidence["host_evidence_verified"]
    assert not Path(evidence["genbank_path"]).is_absolute()


def test_corrupt_taxonomy_sha_is_rejected(tmp_path, monkeypatch):
    manifest, raw, out = prepare(tmp_path)
    responder(monkeypatch, [search(), taxonomy()])
    hosts.resolve_hosts(manifest, raw, out)
    evidence = hosts.load_host_evidence(out)["SYNTHETIC001.7"]
    path = tmp_path / evidence["host_domain_evidence"]
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="Taxonomy SHA"):
        hosts.load_host_evidence(out)


def test_portable_project_root_and_shared_lookup_cache(tmp_path, monkeypatch):
    original, relocated = tmp_path / "original", tmp_path / "relocated"
    manifest, raw, out = prepare(original)
    responder(monkeypatch, [search(), taxonomy()])
    # fetch wrapper may create a cache whose default data root is the NCBI folder.
    acquire.resolve_host_taxonomy(host_name="Synthetic genus", out=raw)
    monkeypatch.setattr(acquire, "_request", lambda *a, **k: pytest.fail("cache not reused"))
    hosts.resolve_hosts(manifest, raw, out)
    shutil.copytree(original, relocated)
    relocated_out = relocated / out.relative_to(original)
    loaded = hosts.load_host_evidence(relocated_out)
    assert loaded["SYNTHETIC001.7"]["host_domain"] == "Bacteria"
    assert loaded["SYNTHETIC001.7"]["host_evidence_data_root"] == str(relocated.resolve())
    # Explicit root is also supported after transferring the CSV to another location.
    copied = tmp_path / "copied.csv"
    shutil.copyfile(relocated_out, copied)
    shutil.copyfile(Path(str(relocated_out) + ".meta.json"), Path(str(copied) + ".meta.json"))
    assert (
        hosts.load_host_evidence(copied, data_root=relocated)["SYNTHETIC001.7"]["host_domain"]
        == "Bacteria"
    )


def test_input_identity_or_sequence_tampering_is_rejected(tmp_path, monkeypatch):
    manifest, raw, out = prepare(tmp_path)
    responder(monkeypatch, [search(), taxonomy()])
    hosts.resolve_hosts(manifest, raw, out)
    gb = raw / "record.gb"
    gb.write_text(gb.read_text().replace("atgaaataa", "atgccctaa"), encoding="utf-8")
    with pytest.raises(ValueError, match="GenBank SHA"):
        hosts.load_host_evidence(out)


def test_changed_metadata_is_re_resolved_without_sequence_fetch(tmp_path, monkeypatch):
    manifest, raw, out = prepare(
        tmp_path, [{"accession": "SYNTHETIC001", "metadata_host": "Unspecified"}]
    )
    monkeypatch.setattr(acquire, "_request", lambda *a, **k: pytest.fail("placeholder queried"))
    hosts.resolve_hosts(manifest, raw, out)
    manifest.write_text("accession,metadata_host\nSYNTHETIC001,Synthetic genus\n", encoding="utf-8")
    calls = responder(monkeypatch, [search(), taxonomy()])
    result = hosts.resolve_hosts(manifest, raw, out)
    assert result["counts"] == {"verified": 1} and len(calls) == 2
    assert len(Path(str(out) + ".events.jsonl").read_text().splitlines()) == 2


def test_missing_genbank_does_not_trigger_taxonomy_or_invent_version(tmp_path, monkeypatch):
    manifest, raw, out = prepare(
        tmp_path, [{"accession": "NOTFETCHED", "metadata_host": "Synthetic genus"}]
    )
    monkeypatch.setattr(
        acquire, "_request", lambda *a, **k: pytest.fail("missing sequence queried")
    )
    result = hosts.resolve_hosts(manifest, raw, out)
    assert result["records"][0]["accession_version"] == ""
    assert hosts.load_host_evidence(out) == {}


def test_evidence_csv_tampering_rejected(tmp_path, monkeypatch):
    manifest, raw, out = prepare(tmp_path)
    responder(monkeypatch, [search(), taxonomy()])
    hosts.resolve_hosts(manifest, raw, out)
    out.write_text(out.read_text().replace("Bacteria", "Archaea"), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        hosts.load_host_evidence(out)


def test_streaming_cli_mode_counts_rows_without_retaining_report_payload(tmp_path, monkeypatch):
    manifest, raw, out = prepare(tmp_path)
    responder(monkeypatch, [search(), taxonomy()])
    result = hosts.resolve_hosts(manifest, raw, out, include_records=False)
    assert result["requested"] == 1 and result["counts"] == {"verified": 1}
    assert result["records"] == []
    assert hosts.load_host_evidence(out)["SYNTHETIC001.7"]["host_domain"] == "Bacteria"
    assert len(Path(str(out) + ".rows.jsonl").read_text().splitlines()) == 1
