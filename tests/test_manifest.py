"""All generated sequences and manifests here are synthetic test data only."""

import csv
import hashlib
import json

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from phageannoscope.genbank_io import genome_from_record
from phageannoscope.manifest import select_manifest, validate_cohort


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def genome(
    path, version="TEST001.7", sequence="ATGAAATAA", *, mol_type="genomic DNA", segment=None
):
    record = SeqRecord(
        Seq(sequence),
        id=version,
        name=version.split(".")[0],
        description="SYNTHETIC TEST ONLY",
        annotations={
            "molecule_type": "DNA",
            "topology": "linear",
            "taxonomy": ["Viruses", "Synthetic test clade"],
        },
    )
    qualifiers = {
        "mol_type": [mol_type],
        "db_xref": ["taxon:999999999"],
        "host": ["Synthetic test host"],
    }
    if segment:
        qualifiers["segment"] = [segment]
    record.features = [
        SeqFeature(SimpleLocation(0, len(sequence)), type="source", qualifiers=qualifiers)
    ]
    SeqIO.write(record, path, "genbank")


def metadata(**changes):
    row = {
        "accession": "TEST001",
        "accession_version": "",
        "metadata_length_bp": 9,
        "baltimore_group": "Group I",
        "checkv_quality": "Complete",
        "checkv_miuvig_quality": "Complete",
        "host_domain": "Bacteria",
        "host_domain_status": "verified",
        "host_domain_evidence": "synthetic_test_taxonomy_evidence",
    }
    return {**row, **changes}


def test_unknown_host_is_hold_and_frozen_input_preserved(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "synthetic.gb")
    manifest = tmp_path / "frozen.csv"
    write_csv(manifest, [metadata(host_domain="", host_domain_status="not_checked")])
    frozen = manifest.read_bytes()
    result = validate_cohort(manifest, raw, tmp_path / "validated.csv")
    row = read_csv(tmp_path / "validated.csv")[0]
    assert result["eligible_for_pilot"] == 0
    assert "HOST_DOMAIN_UNRESOLVED" in json.loads(row["all_reason_codes_json"])
    assert row["accession_version"] == "TEST001.7"
    assert (
        row["genome_id"]
        == genome_from_record(SeqIO.read(raw / "synthetic.gb", "genbank"))["genome_id"]
    )
    assert manifest.read_bytes() == frozen
    assert row["in_primary_analysis"] == "false"


def test_all_reasons_retained_and_quarantine_separate(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "synthetic.gb", sequence="ATGANATAA", segment="1")
    write_csv(
        tmp_path / "frozen.csv",
        [metadata(metadata_length_bp=10, checkv_quality="Low-quality", host_domain="Archaea")],
    )
    result = validate_cohort(tmp_path / "frozen.csv", raw, tmp_path / "validated.csv")
    row = read_csv(tmp_path / "validated.csv")[0]
    reasons = json.loads(row["all_reason_codes_json"])
    assert {
        "BELOW_PRIMARY_QUALITY",
        "METADATA_NCBI_LENGTH_MISMATCH",
        "N_CONTAINING_SEQUENCE_SENSITIVITY_ONLY",
        "SEGMENTED_GENOME_SEPARATE_LAYER",
        "NON_BACTERIAL_HOST",
    } <= set(reasons)
    assert row["selection_status"] == "quarantine" and result["quarantine_rows"] == 1
    validate_cohort(tmp_path / "frozen.csv", raw, tmp_path / "validated.csv")
    assert len(read_csv(tmp_path / "selection_events.csv")) == 2


def test_sequence_dedup_has_one_representative_aliases_preserved(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "a.gb", "AAA001.2")
    genome(raw / "b.gb", "BBB001.4")
    write_csv(tmp_path / "frozen.csv", [metadata(accession="BBB001"), metadata(accession="AAA001")])
    summary = validate_cohort(tmp_path / "frozen.csv", raw, tmp_path / "validated.csv")
    rows = read_csv(tmp_path / "validated.csv")
    assert summary["eligible_for_pilot"] == 1
    assert {r["exact_sequence_representative"] for r in rows} == {"AAA001.2"}
    assert (
        next(r for r in rows if r["accession"] == "BBB001")["selection_status"]
        == "exact_sequence_alias"
    )


def test_validation_reuses_a_bounded_raw_batch_cache(tmp_path, monkeypatch):
    from phageannoscope.hosts import _record_batch

    raw = tmp_path / "raw"
    raw.mkdir()
    batch = raw / "batch.gb"
    records = []
    rows = []
    for number, sequence in enumerate(("ATGAAATAA", "ATGCCCTAA", "ATGGGGTAA"), 1):
        temporary = tmp_path / f"record{number}.gb"
        accession = f"TEST{number:03}"
        genome(temporary, f"{accession}.1", sequence)
        records.append(SeqIO.read(temporary, "genbank"))
        rows.append(metadata(accession=accession))
    SeqIO.write(records, batch, "genbank")
    manifest = tmp_path / "frozen.csv"
    write_csv(manifest, rows)
    original_parse = SeqIO.parse
    parsed_batches = []

    def counted_parse(handle, format, *args, **kwargs):
        if getattr(handle, "name", None) == str(batch.resolve()):
            parsed_batches.append(handle.name)
        return original_parse(handle, format, *args, **kwargs)

    _record_batch.cache_clear()
    monkeypatch.setattr(SeqIO, "parse", counted_parse)
    result = validate_cohort(manifest, raw, tmp_path / "validated.csv")
    assert result["eligible_for_pilot"] == 3
    assert len(parsed_batches) == 2  # One index scan, then one shared biological batch.
    assert _record_batch.cache_info().maxsize == 2
    _record_batch.cache_clear()


@pytest.mark.parametrize(
    "mol_type,sequence,expected",
    [
        ("genomic RNA", "ATGAAATAA", "NON_DNA_OR_AMBIGUOUS_MOL_TYPE"),
        ("genomic DNA", "ATGXAATAA", "NON_IUPAC_SEQUENCE"),
    ],
)
def test_rna_and_invalid_alphabet_never_eligible(tmp_path, mol_type, sequence, expected):
    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "synthetic.gb", mol_type=mol_type, sequence=sequence)
    write_csv(tmp_path / "frozen.csv", [metadata()])
    result = validate_cohort(tmp_path / "frozen.csv", raw, tmp_path / "validated.csv")
    row = read_csv(tmp_path / "validated.csv")[0]
    assert result["eligible_for_pilot"] == 0
    assert expected in json.loads(row["all_reason_codes_json"])


def selection_rows():
    rows = []
    for stratum, length in enumerate([9000, 10000, 50000, 100000, 200000]):
        for j in range(30):
            version = f"SYNTH{stratum}{j:03}.9"
            rows.append(
                {
                    "accession_version": version,
                    "sequence_sha256": hashlib.sha256(version.encode()).hexdigest(),
                    "ncbi_length_bp": length,
                    "eligible_for_pilot": "true",
                    "taxmyphage_species": f"test_species_{j % 20}",
                }
            )
    return rows


def test_deterministic_strata_and_nested_sample(tmp_path):
    rows = selection_rows()
    write_csv(tmp_path / "validated.csv", rows)
    summary = select_manifest(tmp_path / "validated.csv", 100, 30102026, out=tmp_path / "pilot.csv")
    assert summary["selected_n"] == 100
    assert list(summary["strata"].values()) == [20] * 5
    pilot = read_csv(tmp_path / "pilot.csv")
    write_csv(tmp_path / "reordered.csv", rows[::-1])
    select_manifest(tmp_path / "reordered.csv", 100, 30102026, out=tmp_path / "same.csv")
    assert (tmp_path / "same.csv").read_bytes() == (tmp_path / "pilot.csv").read_bytes()
    nested = select_manifest(
        tmp_path / "validated.csv",
        125,
        30102026,
        out=tmp_path / "stress.csv",
        include=tmp_path / "pilot.csv",
    )
    assert nested["nested_included_n"] == 100
    assert {r["accession_version"] for r in pilot} <= {
        r["accession_version"] for r in read_csv(tmp_path / "stress.csv")
    }


def test_preliminary_rows_cannot_be_pilot(tmp_path):
    rows = selection_rows()
    for row in rows:
        row["eligible_for_pilot"] = "false"
        row["selection_status"] = "pending_ncbi_validation"
    write_csv(tmp_path / "frozen.csv", rows)
    with pytest.raises(ValueError, match="Insufficient final eligible"):
        select_manifest(tmp_path / "frozen.csv", 100, out=tmp_path / "pilot.csv")
    assert not (tmp_path / "pilot.csv").exists()


def test_nested_sequence_mismatch_refused(tmp_path):
    rows = selection_rows()
    write_csv(tmp_path / "validated.csv", rows)
    write_csv(tmp_path / "include.csv", [{**rows[0], "sequence_sha256": "different"}])
    with pytest.raises(ValueError, match="SHA mismatch"):
        select_manifest(
            tmp_path / "validated.csv",
            100,
            out=tmp_path / "pilot.csv",
            include=tmp_path / "include.csv",
        )


def test_primary_reason_tracks_quarantine_or_exclusion_not_earlier_hold(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "synthetic.gb")
    manifest = tmp_path / "frozen.csv"
    write_csv(manifest, [metadata(baltimore_group="Unclassified", checkv_quality="Low-quality")])
    out = tmp_path / "validated.csv"
    validate_cohort(manifest, raw, out)
    row = read_csv(out)[0]
    assert row["selection_status"] == "excluded_primary"
    assert row["primary_reason_code"] == "BELOW_PRIMARY_QUALITY"
    assert "MOLECULE_UNRESOLVED" in json.loads(row["all_reason_codes_json"])
    events = read_csv(tmp_path / "molecule_resolution_events.csv")
    assert events[0]["decision"] == "review_required"
    assert json.loads(events[0]["evidence_json"])["mol_types"] == ["genomic DNA"]
    assert row["eligible_for_pilot"] == "false"
    write_csv(manifest, [metadata(baltimore_group="Unclassified", genbank_sha256="wrong")])
    validate_cohort(manifest, raw, out)
    assert read_csv(out)[0]["primary_reason_code"] == "FROZEN_GENBANK_SHA_MISMATCH"


def test_new_host_evidence_multiline_provenance_and_stale_inputs_hold(tmp_path, monkeypatch):
    from phageannoscope import acquire, hosts

    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "synthetic.gb")
    genome(raw / "second.gb", "TEST002.3", "ATGCCCTAA")
    manifest = tmp_path / "frozen.csv"
    write_csv(
        manifest,
        [
            metadata(metadata_host="Escherichia"),
            metadata(accession="TEST002", metadata_host="Escherichia"),
        ],
    )
    fixtures = __import__("pathlib").Path(__file__).parent / "fixtures/real/host_taxonomy"
    search = (
        fixtures / "754a9d5a62c3f39cc1a3d7c1feed820ca6da9b6630e62a268f7037aade33f8a2.xml"
    ).read_bytes()
    taxonomy = (
        fixtures / "9f0b9c5c5028969e93146e577cfaab4eb45bcd837bedda5842e422dd369e7b80.xml"
    ).read_bytes()
    missing = b"<eSearchResult><Count>0</Count><IdList/></eSearchResult>"

    def request(endpoint, params, **kwargs):
        assert params["db"] == "taxonomy"
        data = (
            taxonomy
            if endpoint == "efetch.fcgi"
            else search
            if "Escherichia" in params["term"]
            else missing
        )
        return data, {"http_status": 200, "request_url": "https://example.invalid/test"}

    monkeypatch.setattr(acquire, "_request", request)
    evidence = tmp_path / "hosts.csv"
    hosts.resolve_hosts(manifest, raw, evidence)
    out = tmp_path / "validated.csv"
    result = validate_cohort(manifest, raw, out, host_evidence=evidence)
    assert result["eligible_for_pilot"] == 2
    rows = read_csv(out)
    assert all(row["host_taxid"] == "561" and row["host_rank"] == "genus" for row in rows)
    assert {row["host_evidence_table"] for row in rows} == {str(evidence.resolve())}
    assert result["host_evidence"] == str(evidence.resolve())
    assert json.loads(out.with_suffix(".summary.json").read_text())["host_evidence"] == str(
        evidence.resolve()
    )
    write_csv(manifest, [metadata(metadata_host="Synthetic changed relation")])
    result = validate_cohort(manifest, raw, out, host_evidence=evidence)
    assert result["eligible_for_pilot"] == 0
    assert "HOST_EVIDENCE_INPUTS_CHANGED" in json.loads(read_csv(out)[0]["all_reason_codes_json"])


def test_non_n_ambiguity_codes_follow_n_rule(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    genome(raw / "synthetic.gb", sequence="ATGARYTAA")
    write_csv(tmp_path / "frozen.csv", [metadata(metadata_length_bp=9)])
    validate_cohort(tmp_path / "frozen.csv", raw, tmp_path / "validated.csv")
    row = read_csv(tmp_path / "validated.csv")[0]
    assert "AMBIGUOUS_IUPAC_SENSITIVITY_ONLY" in json.loads(row["all_reason_codes_json"])
    assert row["selection_status"] == "excluded_primary"
    assert json.loads(row["validation_evidence_json"])["AMBIGUOUS_IUPAC_SENSITIVITY_ONLY"] == {
        "R": 1, "Y": 1}
