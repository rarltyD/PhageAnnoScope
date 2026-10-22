"""Typed output and actual frozen NCBI + installed-caller integration."""

import csv
import json
import os
import shutil
from pathlib import Path

import duckdb
import pyarrow as pa
import pytest

from phageannoscope.callers.phanotate import command
from phageannoscope.callers.runner import run_caller
from phageannoscope.config import load_config
from phageannoscope.genbank_io import read_genbank
from phageannoscope.pipeline import (
    apply_primary_code_policy,
    mark_rotation_sensitivity,
    run_pipeline,
)
from phageannoscope.storage import ResultWriter, arrow_table, build_db
from phageannoscope.util import file_sha, versions
from phageannoscope.validation import validate_results

REAL = Path(__file__).parent / "fixtures" / "real" / "ncbi_2026-10-16.gb"
EXTERNAL_CONFIG = os.environ.get("PHAGEANNOSCOPE_TEST_CONFIG")


def test_empty_schema_is_typed_and_has_real_nulls():
    table = arrow_table(
        "cds_calls", [{"call_id": "synthetic", "strand": None, "stop_triplet": None}]
    )
    assert table.schema.field("strand").type == pa.int64()
    assert table.schema.field("stop_triplet").type == pa.list_(pa.int64())
    assert table["strand"].to_pylist() == [None]


def test_writer_deduplicates_on_disk_and_removes_stale_shards(tmp_path):
    root = tmp_path / "parquet"
    writer = ResultWriter(root, batch_genomes=1)
    writer.add({"genomes": [{"genome_id": "synthetic1"}]})
    writer.add({"genomes": [{"genome_id": "synthetic2"}]})
    writer.close()
    assert len(list((root / "genomes").glob("*.parquet"))) == 2
    writer = ResultWriter(root)
    for _ in range(2):
        writer.add(
            {
                "genomes": [{"genome_id": "synthetic1"}],
                "proteins": [
                    {"protein_sha256": "same", "amino_acid_sequence": "M", "length_aa": 1}
                ],
            }
        )
    counts = writer.close()
    assert counts["genomes"] == counts["proteins"] == 1
    assert len(list((root / "genomes").glob("*.parquet"))) == 1
    db_path = tmp_path / "x.duckdb"
    build_db(root, db_path)
    with duckdb.connect(str(db_path), read_only=True) as db:
        assert db.execute("SELECT count(*) FROM genomes").fetchone()[0] == 1
        assert db.execute("SELECT count(*) FROM cds_calls").fetchone()[0] == 0
    moved = tmp_path / "moved"
    shutil.copytree(root, moved)
    build_db(moved, tmp_path / "moved.duckdb")
    with duckdb.connect(str(tmp_path / "moved.duckdb"), read_only=True) as db:
        assert db.execute("SELECT genome_id FROM genomes").fetchone()[0] == "synthetic1"


def test_rotation_exclusions_preserve_reasons_and_require_real_two_tracks():
    calls = [
        {
            "source": "pyrodigal_gv",
            "model_id": "same",
            "rotation_shift": shift,
            "eligible_primary": False,
            "primary_exclusion_reason": "partial",
            "primary_exclusion_reasons": ["partial"],
        }
        for shift in (0, 50)
    ]
    runs = [
        {
            "source": "pyrodigal_gv",
            "source_run_id": str(shift),
            "rotation_shift": shift,
            "status": "success",
        }
        for shift in (0, 50)
    ]
    mark_rotation_sensitivity(calls, runs, "circular")
    assert not calls[0]["rotation_sensitive"]
    assert calls[1]["eligible_primary"] is False
    assert "partial" in calls[1]["primary_exclusion_reasons"]
    assert "rotation_track" in calls[1]["primary_exclusion_reasons"]
    runs[1]["rotation_shift"] = 0
    mark_rotation_sensitivity(calls, runs, "circular")
    assert all("rotation_check_incomplete" in c["primary_exclusion_reasons"] for c in calls)


def test_alternative_code_holds_whole_primary_pair_without_changing_execution_status():
    calls = [
        {"source": "phanotate", "transl_table_effective": 11, "eligible_primary": True},
        {"source": "pyrodigal_gv", "transl_table_effective": 15, "eligible_primary": True},
    ]
    runs = [{"source": source, "status": "success"} for source in ("phanotate", "pyrodigal_gv")]
    assert apply_primary_code_policy(calls, runs) is False
    assert all(c["eligible_primary"] is False for c in calls)
    assert all(
        r["status"] == "success" and r["comparison_status"] == "not_comparable" for r in runs
    )


def test_real_snapshot_integrity_and_known_coordinates():
    provenance = json.loads(REAL.with_suffix(".provenance.json").read_text())
    assert file_sha(REAL) == provenance["sha256"]
    tables = read_genbank(REAL)
    assert {g["accession_version"]: g["length_bp"] for g in tables["genomes"]} == {
        "J02459.1": 48502,
        "J02482.1": 5386,
        "V01146.1": 39937,
    }
    phix = next(g for g in tables["genomes"] if g["accession_version"] == "J02482.1")
    assert phix["topology"] == "circular"
    assert any(
        c["wraps_origin"] for c in tables["cds_calls"] if c["genome_id"] == phix["genome_id"]
    )


@pytest.mark.integration
def test_real_pyrodigal_cache_preserves_time_and_detects_corruption(tmp_path):
    genome = next(g for g in read_genbank(REAL)["genomes"] if g["accession_version"] == "J02482.1")
    first, calls = run_caller(genome, "pyrodigal_gv", load_config(), tmp_path)
    assert first["status"] == "success"
    assert len(calls) > 0
    second, cached = run_caller(genome, "pyrodigal_gv", load_config(), tmp_path)
    assert second["status"] == "cache_hit"
    assert second["wall_seconds"] == first["wall_seconds"]
    assert cached == calls
    saved_calls = next(tmp_path.rglob("calls.json"))
    saved_calls.write_text("[]")
    third, restored = run_caller(genome, "pyrodigal_gv", load_config(), tmp_path)
    assert third["status"] == "success"
    assert len(restored) == len(calls)


@pytest.mark.integration
def test_real_end_to_end_roundtrip(tmp_path, monkeypatch):
    # This test intentionally documents a missing independent caller.
    monkeypatch.setattr("phageannoscope.callers.runner.command", lambda configured=None: None)
    parsed = read_genbank(REAL)
    genome = next(g for g in parsed["genomes"] if g["accession_version"] == "J02482.1")
    manifest = tmp_path / "diagnostic.csv"
    row = {
        "accession_version": genome["accession_version"],
        "sequence_sha256": genome["sequence_sha256"],
        "genbank_path": str(REAL.resolve()),
        "genbank_sha256": file_sha(REAL),
        "eligible_for_pilot": "false",
    }
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    config = load_config()
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="Unvalidated"):
        run_pipeline(manifest, config, output, raw_root=tmp_path / "raw")
    summary = run_pipeline(manifest, config, output, diagnostic=True, raw_root=tmp_path / "raw")
    assert summary["imported_n"] == 1
    assert summary["caller_completed_n"] == summary["eligible_primary_n"] == 0
    assert validate_results(output)["status"] == "passed"
    with duckdb.connect(str(output / "phageannoscope.duckdb"), read_only=True) as db:
        assert (
            db.execute("SELECT count(*) FROM stop_loci WHERE eligible_primary").fetchone()[0] == 0
        )
        assert (
            db.execute("SELECT count(*) FROM source_runs WHERE source='pyrodigal_gv'").fetchone()[0]
            == 2
        )
        assert (
            db.execute(
                "SELECT count(*) FROM cds_calls WHERE rotation_shift>0 AND eligible_primary"
            ).fetchone()[0]
            == 0
        )
    repeated = run_pipeline(manifest, config, output, diagnostic=True, raw_root=tmp_path / "raw")
    assert repeated["table_counts"] == summary["table_counts"]
    assert repeated["source_status_counts"]["cache_hit"] == 2
    assert validate_results(output)["status"] == "passed"


@pytest.mark.integration
@pytest.mark.skipif(
    not EXTERNAL_CONFIG and (not versions()["phanotate"] or not command()),
    reason="Set PHAGEANNOSCOPE_TEST_CONFIG for external PHANOTATE or install native PHANOTATE",
)
def test_frozen_real_genome_with_both_callers(tmp_path):
    genome = next(g for g in read_genbank(REAL)["genomes"] if g["accession_version"] == "J02482.1")
    for source in ("phanotate", "pyrodigal_gv"):
        run, calls = run_caller(genome, source, load_config(EXTERNAL_CONFIG), tmp_path)
        assert run["status"] == "success"
        assert calls
        if source == "phanotate":
            assert all(c["terminal_stop_present"] for c in calls if not c["partial_3p"])
    manifest = tmp_path / "paired.csv"
    row = {
        "accession_version": genome["accession_version"],
        "sequence_sha256": genome["sequence_sha256"],
        "genbank_path": str(REAL.resolve()),
        "genbank_sha256": file_sha(REAL),
        "eligible_for_pilot": "true",
    }
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    config = load_config(EXTERNAL_CONFIG)
    output = tmp_path / "paired-run"
    first = run_pipeline(
        manifest, config, output, diagnostic=True, raw_root=tmp_path / "paired-raw"
    )
    assert first["caller_completed_n"] == 1
    assert first["eligible_primary_n"] == 0  # Diagnostic never enters the research N.
    assert validate_results(output)["status"] == "passed"
    second = run_pipeline(
        manifest, config, output, diagnostic=True, raw_root=tmp_path / "paired-raw"
    )
    assert second["checkpoint_hits"] == 1
    assert second["table_counts"] == first["table_counts"]
    assert second["source_status_counts"]["cache_hit"] == 4  # Both rotations of both callers.
    assert second["fresh_caller_cpu_seconds"] is None
    assert validate_results(output)["status"] == "passed"
    checkpoint = next((output / "checkpoints").glob("*.json"))
    checkpoint.write_text("{}", encoding="utf-8")
    repaired = run_pipeline(
        manifest, config, output, diagnostic=True, raw_root=tmp_path / "paired-raw"
    )
    assert repaired["checkpoint_hits"] == 0
    assert repaired["table_counts"] == first["table_counts"]
    assert validate_results(output)["status"] == "passed"
