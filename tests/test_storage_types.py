"""Null/failed batches must preserve the types of later biological results."""

from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from phageannoscope.genbank_io import read_genbank
from phageannoscope.matching import build_loci, build_models, compare_calls, genome_metrics
from phageannoscope.proteins import build_catalogue
from phageannoscope.storage import SCHEMAS, ResultWriter, arrow_table, build_db


def test_failed_then_successful_shards_keep_numeric_bool_and_list_types(tmp_path):
    fields = {
        "genome_metrics": {
            "exact_a_count": 7,
            "stop_union": 9,
            "length_nt_median": 123.5,
            "same_strand_overlap_pairs_ge1nt": 3,
            "antisense_overlapping_models_ge30nt": 2,
            "short_protein_fraction": 0.25,
            "not_comparable_count": 1,
        },
        "source_runs": {
            "selected_translation_table": 11,
            "cds_count": 7,
            "prediction_cpu_seconds": 0.125,
            "command_argv": ["python", "worker.py"],
        },
        "cds_calls": {
            "native_begin": 1,
            "native_partial_begin": False,
            "automatic_code_selection": True,
            "location_issues": [],
            "translation_notes": ["alternative_initiation_met"],
            "transl_except": ["(pos:3,aa:Sec)"],
        },
        "models": {"translation_tables": [1, 11], "genome_ids": ["g"], "has_partial_source": False},
        "stop_loci": {
            "caller_support_mask": 3,
            "phanotate_completed": True,
            "additional_stop_model": False,
            "caller_families": ["phanotate", "pyrodigal_gv"],
        },
        "locus_models": {"call_ids": ["c"], "source_run_ids": ["r"], "is_representative": True},
        "protein_origins": {
            "provided_translation_present": True,
            "provided_translation_matches": False,
        },
        "features": {"feature_ordinal": 0, "partial_5p": False, "location_issues": ["test_issue"]},
        "matches": {"substantial_overlap": True, "protein_translation_mismatch": False},
    }
    writer = ResultWriter(tmp_path / "parquet", batch_genomes=1)
    failed = {name: [{key: None for key in values}] for name, values in fields.items()}
    successful = {name: [dict(values)] for name, values in fields.items()}
    for name in ("models", "locus_models"):
        failed[name][0]["model_id"] = "failed"
        successful[name][0]["model_id"] = "successful"
    writer.add(failed)
    writer.add(successful)
    writer.close()
    build_db(tmp_path / "parquet", tmp_path / "results.duckdb")
    with duckdb.connect(str(tmp_path / "results.duckdb"), read_only=True) as connection:
        for name, values in fields.items():
            first = pq.read_table(tmp_path / "parquet" / name / "part-00000.parquet")
            second = pq.read_table(tmp_path / "parquet" / name / "part-00001.parquet")
            for field, value in values.items():
                assert first.schema.field(field).type == second.schema.field(field).type
                stored = connection.execute(f'SELECT "{field}" FROM "{name}"').fetchall()
                assert stored == [(None,), (value,)]
            types = dict(
                (row[0], row[1]) for row in connection.execute(f'DESCRIBE "{name}"').fetchall()
            )
            for field, value in values.items():
                if isinstance(value, bool):
                    assert types[field] == "BOOLEAN"
                elif isinstance(value, int):
                    assert types[field] == "BIGINT"
                elif isinstance(value, float):
                    assert types[field] == "DOUBLE"
                elif isinstance(value, list):
                    assert types[field] == (
                        "BIGINT[]" if value and isinstance(value[0], int) else "VARCHAR[]"
                    )


def test_unknown_all_null_extension_can_later_be_numeric(tmp_path):
    first = arrow_table("genome_metrics", [{"future_metric": None}])
    second = arrow_table(
        "genome_metrics", [{"future_metric": 2.5, "future_count": 3, "future_flag": True}]
    )
    pq.write_table(first, tmp_path / "first.parquet")
    pq.write_table(second, tmp_path / "second.parquet")
    with duckdb.connect() as connection:
        rows = connection.execute(
            "SELECT future_metric FROM read_parquet(?,union_by_name=true)",
            [[str(tmp_path / "first.parquet"), str(tmp_path / "second.parquet")]],
        ).fetchall()
    assert rows == [(None,), (2.5,)]
    assert pa.types.is_null(first.schema.field("future_metric").type)


def test_declared_open_params_dict_is_json_and_known_numeric_rejects_strings():
    table = arrow_table(
        "source_runs",
        [{"actual_params": {"meta": True, "min_gene": 90}, "selected_translation_table": 11}],
    )
    assert table.schema.field("selected_translation_table").type == pa.int64()
    assert table.column("actual_params").to_pylist() == ['{"meta": true, "min_gene": 90}']
    with pytest.raises(pa.ArrowInvalid):
        arrow_table("source_runs", [{"selected_translation_table": "11"}])


def test_current_normalized_biological_outputs_have_declared_schema_fields():
    # A real frozen input exercises nullable/raw biological fields that ordinary
    # valid CDS would miss. No predictions or research results are fabricated.
    tables = read_genbank(
        Path(__file__).parent / "fixtures" / "real" / "ncbi_2026-10-16.gb", source_run_id="baseline"
    )
    genome = next(g for g in tables["genomes"] if g["accession_version"] == "J02482.1")
    calls = [c for c in tables["cds_calls"] if c["genome_id"] == genome["genome_id"]]
    runs = [
        {
            "source_run_id": "baseline",
            "genome_id": genome["genome_id"],
            "source": "genbank",
            "input_sha256": genome["sequence_sha256"],
            "status": "success",
        },
        {
            "source_run_id": "failed_caller",
            "genome_id": genome["genome_id"],
            "source": "phanotate",
            "input_sha256": genome["sequence_sha256"],
            "status": "failed",
        },
    ]
    matches = compare_calls(calls, runs)
    loci = build_loci(calls, runs)
    tables.update(build_models(calls))
    tables.update(build_catalogue(calls))
    tables.update(loci)
    tables["matches"] = matches
    tables["genome_metrics"] = genome_metrics(calls, matches, loci, runs)
    tables["source_runs"] = runs
    for name, rows in tables.items():
        observed = {field for row in rows for field in row}
        assert observed <= set(SCHEMAS[name]), (
            f"Undeclared {name} fields: {observed - set(SCHEMAS[name])}"
        )
        # Declared types must accept the actual values, including all-null lists.
        arrow_table(name, rows)
