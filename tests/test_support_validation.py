"""Independent support checks, including the unchanged real pre-review artifact."""

from pathlib import Path

import duckdb

from phageannoscope.matching import build_loci
from phageannoscope.storage import arrow_table, sql_quote
from phageannoscope.validation import audit_support_gate


def run(source, status="success"):
    return {
        "source_run_id": source,
        "genome_id": "synthetic",
        "source": source,
        "algorithm_family": source,
        "input_sha256": "synthetic-sha",
        "status": status,
        "rotation_shift": 0,
    }


def call(source, excluded=False):
    return {
        "call_id": source,
        "model_id": "synthetic-model",
        "source_run_id": source,
        "genome_id": "synthetic",
        "sequence_sha256": "synthetic-sha",
        "source": source,
        "strand": 1,
        "segments": [{"start0": 0, "end0": 9}],
        "codon_start": 1,
        "code_class": "standard_stop",
        "transl_table_effective": 11,
        "translated_nt_length": 9,
        "stop_triplet": [6, 7, 8],
        "eligible_primary": not excluded,
        "rotation_sensitive": excluded,
        "primary_exclusion_reason": "rotation_sensitive" if excluded else None,
        "rotation_shift": 0,
    }


def connection(calls, runs, mutate=None):
    tables = {"cds_calls": calls, "source_runs": runs, **build_loci(calls, runs)}
    for locus in tables["stop_loci"]:
        locus["geometry_eligible_primary"] = locus["eligible_primary"]
        locus["eligible_primary"] = False  # Cohort overlay cannot conceal a technical bug.
    if mutate:
        mutate(tables)
    db = duckdb.connect()
    for name, rows in tables.items():
        db.register(name, arrow_table(name, rows))
    return db


def test_real_original_geometry_error_is_detected_without_mutating_artifact():
    original = Path(__file__).resolve().parents[1] / "results/diagnostic/parquet"
    if not original.exists():
        import pytest

        pytest.skip("Historical local run not distributed; unit regression remains active")
    with duckdb.connect() as db:
        for name in ("cds_calls", "source_runs", "stop_loci", "locus_models"):
            glob = (original / name / "part-*.parquet").as_posix()
            db.execute(
                f"CREATE VIEW {name} AS SELECT * FROM read_parquet({sql_quote(glob)},union_by_name=true)"
            )
        report = audit_support_gate(db)
    assert report["n_admission_inconsistencies"] == 1
    assert (
        report["admission_cases"][0]["locus_id"]
        == "locus_d4a6db6bdb2eacaf3c1bc28c81946f326e7f0d993d9f366debc10420129b3403"
    )


def test_new_gate_matches_all_defining_calls():
    with connection(
        [call("phanotate", True), call("pyrodigal_gv")], [run("phanotate"), run("pyrodigal_gv")]
    ) as db:
        assert audit_support_gate(db)["status"] == "passed"
        row = db.execute(
            "SELECT caller_support_raw,caller_support_eligible,support_assessable_primary FROM stop_loci"
        ).fetchone()
        assert row == (2, 1, False)


def test_forged_admission_and_mask_cannot_pass():
    def mutate(tables):
        locus = tables["stop_loci"][0]
        locus.update(
            geometry_eligible_primary=True,
            support_assessable_primary=True,
            caller_support_eligible=2,
            caller_support_mask_eligible=3,
        )

    with connection(
        [call("phanotate", True), call("pyrodigal_gv")],
        [run("phanotate"), run("pyrodigal_gv")],
        mutate,
    ) as db:
        report = audit_support_gate(db)
    assert report["status"] == "failed" and report["n_admission_inconsistencies"] == 1
    assert any("counts/masks" in error for error in report["errors"])


def test_representative_must_belong_to_admissible_call_domain():
    def mutate(tables):
        tables["stop_loci"][0]["representative_call_id"] = "phanotate"

    with connection(
        [call("phanotate", True), call("pyrodigal_gv")],
        [run("phanotate"), run("pyrodigal_gv")],
        mutate,
    ) as db:
        assert any("Representative" in e for e in audit_support_gate(db)["errors"])


def test_successful_zero_peer_is_assessable_but_failed_peer_is_not():
    for status, expected in [("success", True), ("timeout", False)]:
        with connection([call("phanotate")], [run("phanotate"), run("pyrodigal_gv", status)]) as db:
            assert audit_support_gate(db)["status"] == "passed"
            assert (
                db.execute("SELECT support_assessable_primary FROM stop_loci").fetchone()[0]
                is expected
            )
