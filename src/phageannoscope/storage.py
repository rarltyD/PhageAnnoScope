"""Bounded Parquet batches and relocatable DuckDB rebuilds."""

import json
import sqlite3
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from . import SCHEMA_VERSION

S = pa.string()
INT = pa.int64()
B = pa.bool_()
F = pa.float64()
SEG = pa.list_(pa.struct([("start0", INT), ("end0", INT)]))
STRINGS = pa.list_(S)
INTEGERS = pa.list_(INT)

# Every current normalized output field has a declared type. In particular, a
# failed first batch cannot turn a later numeric/boolean/list column into VARCHAR.
SCHEMAS = {
    "genomes": {
        "primary_exclusion_reason": S,
        "host_evidence_data_root": S,
        "genome_id": S,
        "accession": S,
        "accession_version": S,
        "sequence_sha256": S,
        "length_bp": INT,
        "topology": S,
        "host_domain": S,
        "host_taxid": INT,
        "host_domain_evidence": S,
        "host_scientific_name": S,
        "host_taxonomy_rank": S,
        "host_taxonomy_lineage": S,
        "host_evidence_level": S,
        "host_taxonomy_evidence_json": S,
        "host_inputs_sha256": S,
        "host_evidence_verified": B,
        "in_primary_analysis": B,
        "selection_status": S,
        "sequence": S,
        "organism": S,
        "annotations_json": S,
        "source_qualifiers_json": S,
        "fasta_header": S,
        "host_domain_status": S,
        "checkv_quality": S,
        "baltimore_group": S,
        "relatedness_cluster": S,
        "exact_sequence_duplicate_group": S,
    },
    "source_runs": {
        "source_run_id": S,
        "genome_id": S,
        "source": S,
        "algorithm_family": S,
        "status": S,
        "rotation_shift": INT,
        "input_sha256": S,
        "tool_version": S,
        "params_json": S,
        "started_at_utc": S,
        "ended_at_utc": S,
        "exit_code": INT,
        "wall_seconds": F,
        "cpu_seconds": F,
        "peak_rss_bytes": INT,
        "cache_key": S,
        "cache_origin_status": S,
        "error": S,
        "cache_read_at_utc": S,
        "adapter_sha256": S,
        "raw_genbank_sha256": S,
        "raw_path": S,
        "raw_stdout": S,
        "raw_stderr": S,
        "command_argv": STRINGS,
        "input_fasta_sha256": S,
        "calls_sha256": S,
        "cds_count": INT,
        "selected_translation_table": INT,
        "selected_model": S,
        "actual_params": S,
        "worker_cpu_seconds": F,
        "prediction_cpu_seconds": F,
        "resource_method": S,
        "cpu_unavailable_reason": S,
        "peak_rss_unavailable_reason": S,
        "resource_monitor_status": S,
        "launcher_cpu_unavailable_reason": S,
        "launcher_peak_rss_unavailable_reason": S,
        "comparison_status": S,
        "comparison_reason": S,
        "sequence_sha256": S,
    },
    "features": {
        "feature_id": S,
        "genome_id": S,
        "source_run_id": S,
        "type": S,
        "original_location": S,
        "qualifiers_json": S,
        "accession_version": S,
        "source": S,
        "feature_ordinal": INT,
        "feature_type": S,
        "raw_feature_id": S,
        "strand": INT,
        "segments": SEG,
        "location_operator": S,
        "partial_5p": B,
        "partial_3p": B,
        "location_issues": STRINGS,
    },
    "cds_calls": {
        "call_id": S,
        "feature_id": S,
        "model_id": S,
        "source_run_id": S,
        "genome_id": S,
        "source": S,
        "sequence_sha256": S,
        "strand": INT,
        "segments": SEG,
        "codon_start": INT,
        "transl_table_effective": INT,
        "stop_triplet": INTEGERS,
        "feature_nt_length": INT,
        "translated_nt_length": INT,
        "protein_sha256": S,
        "provided_translation": S,
        "recomputed_translation": S,
        "eligible_primary": B,
        "rotation_shift": INT,
        "rotation_sensitive": B,
        "accession_version": S,
        "feature_ordinal": INT,
        "feature_type": S,
        "type": S,
        "raw_feature_id": S,
        "original_location": S,
        "location_operator": S,
        "qualifiers_json": S,
        "partial_5p": B,
        "partial_3p": B,
        "location_issues": STRINGS,
        "algorithm_family": S,
        "protein_id": S,
        "locus_tag": S,
        "codon_start_raw": S,
        "transl_table_raw": S,
        "transl_table_qualifier_missing": B,
        "qualifier_missing": B,
        "code_class": S,
        "pseudo": B,
        "exception": S,
        "transl_except": STRINGS,
        "wraps_origin": B,
        "segments_contiguous": B,
        "primary_exclusion_reason": S,
        "primary_exclusion_reasons": STRINGS,
        "recomputed_translation_raw": S,
        "nucleotide_sha256": S,
        "translation_validity": S,
        "protein_validation": S,
        "translation_notes": STRINGS,
        "translation_issues": STRINGS,
        "start_codon": S,
        "stop_codon": S,
        "terminal_stop_present": B,
        "translation_normalization": S,
        "raw_gff": S,
        "gff_phase_raw": S,
        "gff_strand_raw": S,
        "gff_attributes_raw": S,
        "edge_partial_inference": S,
        "stop_coordinate_convention": S,
        "native_begin": INT,
        "native_end": INT,
        "native_partial_begin": B,
        "native_partial_end": B,
        "native_translation": S,
        "native_nucleotide": S,
        "selected_model": S,
        "automatic_code_selection": B,
        "code_selection_uncertainty": S,
    },
    "models": {
        "model_id": S,
        "genome_id": S,
        "sequence_sha256": S,
        "segments": SEG,
        "strand": INT,
        "codon_start": INT,
        "code_class": S,
        "genome_ids": STRINGS,
        "feature_nt_length": INT,
        "source_call_count": INT,
        "translation_tables": INTEGERS,
        "has_partial_source": B,
        "has_primary_eligible_source": B,
    },
    "model_sources": {
        "call_id": S,
        "model_id": S,
        "source_run_id": S,
        "source": S,
        "genome_id": S,
        "feature_id": S,
    },
    "model_segments": {"model_id": S, "segment_order": INT, "start0": INT, "end0": INT},
    "matches": {
        "match_id": S,
        "genome_id": S,
        "call_id_a": S,
        "call_id_b": S,
        "model_id_a": S,
        "model_id_b": S,
        "source_a": S,
        "source_b": S,
        "relation": S,
        "nested": B,
        "overlap_nt": INT,
        "coverage_a": F,
        "coverage_b": F,
        "relative_frame": INT,
        "comparable": B,
        "reason": S,
        "sequence_sha256": S,
        "source_run_id_a": S,
        "source_run_id_b": S,
        "substantial_overlap": B,
        "protein_translation_mismatch": B,
    },
    "stop_loci": {
        "locus_id": S,
        "genome_id": S,
        "sequence_sha256": S,
        "representative_model_id": S,
        "caller_support": INT,
        "caller_support_raw": INT,
        "caller_support_eligible": INT,
        "support_assessable_primary": B,
        "eligible_primary": B,
        "geometry_eligible_primary": B,
        "cohort_eligible_primary": B,
        "stop_triplet": INTEGERS,
        "strand": INT,
        "code_class": S,
        "representative_call_id": S,
        "representative_rule": S,
        "caller_support_mask": INT,
        "caller_support_mask_raw": INT,
        "caller_support_mask_eligible": INT,
        "caller_families": STRINGS,
        "caller_families_eligible": STRINGS,
        "source_presence": STRINGS,
        "phanotate_present": B,
        "pyrodigal_gv_present": B,
        "phanotate_completed": B,
        "pyrodigal_gv_completed": B,
        "genbank_present": B,
        "genbank_exact_match": B,
        "genbank_stop_match": B,
        "genbank_comparable": B,
        "universe": S,
        "additional_stop_model": B,
        "unannotated_same_strand_region": B,
        "genbank_antisense_context": B,
        "genbank_partial_context": B,
        "genbank_complex_context": B,
        "primary_exclusion_reason": S,
    },
    "locus_models": {
        "locus_id": S,
        "model_id": S,
        "call_ids": STRINGS,
        "source_run_ids": STRINGS,
        "sources": STRINGS,
        "genome_id": S,
        "is_representative": B,
    },
    "proteins": {
        "protein_sha256": S,
        "amino_acid_sequence": S,
        "length_aa": INT,
        "normalization_version": S,
        "validation": S,
        "translation_provenance": S,
    },
    "protein_origins": {
        "protein_sha256": S,
        "model_id": S,
        "call_id": S,
        "genome_id": S,
        "source_run_id": S,
        "feature_id": S,
        "source": S,
        "protein_id": S,
        "translation_provenance": S,
        "translation_validity": S,
        "provided_translation_present": B,
        "provided_translation_matches": B,
    },
    "selection_events": {
        "genome_id": S,
        "accession": S,
        "stage": S,
        "decision": S,
        "primary_reason": S,
        "evidence": S,
        "run_id": S,
        "accession_version": S,
        "primary_reason_code": S,
        "all_reason_codes_json": S,
        "evidence_json": S,
        "evaluated_at_utc": S,
    },
    "genome_metrics": {
        "genome_id": S,
        "source_a": S,
        "source_b": S,
        "exact_jaccard": F,
        "exact_dice": F,
        "stop_jaccard": F,
        "stop_dice": F,
        "sequence_sha256": S,
        "metric_scope": S,
        "source_run_id_a": S,
        "source_run_id_b": S,
        "status": S,
        "reason": S,
        "rotation_shift": INT,
        "denominator_definition": S,
        "cds_count": INT,
        "observed_call_count": INT,
        "unique_model_count": INT,
        "eligible_call_count": INT,
        "stop_locus_count": INT,
        "primary_locus_count": INT,
        "length_defined_model_count": INT,
        "length_nt_min": INT,
        "length_nt_median": F,
        "length_nt_max": INT,
        "protein_length_defined_model_count": INT,
        "length_aa_min": INT,
        "length_aa_median": F,
        "length_aa_max": INT,
        "short_protein_threshold_aa": INT,
        "short_protein_model_count": INT,
        "short_protein_fraction": F,
        "partial_call_count": INT,
        "wraps_origin_call_count": INT,
        "rotation_sensitive_call_count": INT,
        "complex_biology_call_count": INT,
        "primary_excluded_call_count": INT,
        "overlap_denominator_model_count": INT,
        "same_strand_overlap_pairs_ge1nt": INT,
        "same_strand_overlap_pairs_ge30nt": INT,
        "same_strand_overlapping_models_ge1nt": INT,
        "same_strand_overlapping_models_ge30nt": INT,
        "antisense_overlap_pairs_ge1nt": INT,
        "antisense_overlap_pairs_ge30nt": INT,
        "antisense_overlapping_models_ge1nt": INT,
        "antisense_overlapping_models_ge30nt": INT,
        "exact_a_count": INT,
        "exact_b_count": INT,
        "exact_intersection": INT,
        "exact_union": INT,
        "exact_denominator_status": S,
        "stop_a_count": INT,
        "stop_b_count": INT,
        "stop_intersection": INT,
        "stop_union": INT,
        "stop_denominator_status": S,
        "exact_cds_count": INT,
        "alternative_start_count": INT,
        "contained_distinct_stop_count": INT,
        "partial_overlap_count": INT,
        "antisense_overlap_count": INT,
        "no_correspondence_count": INT,
        "not_comparable_count": INT,
        "substantial_overlap_count": INT,
    },
    "qc_events": {
        "genome_id": S,
        "source_run_id": S,
        "call_id": S,
        "severity": S,
        "code": S,
        "message": S,
        "qc_id": S,
        "source": S,
        "evidence_json": S,
    },
}

# Open-ended caller parameters are intentionally serialized as JSON, rather than
# allowing a batch-dependent Arrow struct layout or coercing scalar measurements.
JSON_COLUMNS = {"source_runs": frozenset({"actual_params"})}


def arrow_table(name, rows):
    """Keep unknown scalar fields; store open-ended dicts as JSON, never numeric strings."""
    rows = [dict(row) for row in rows]
    fixed = SCHEMAS[name]
    columns = list(fixed)
    columns.extend(sorted({key for row in rows for key in row} - set(columns)))
    arrays = []
    for key in columns:
        values = [r.get(key) for r in rows]
        kind = fixed.get(key)
        if key in JSON_COLUMNS.get(name, ()):
            values = [
                json.dumps(v, sort_keys=True) if v is not None and not isinstance(v, str) else v
                for v in values
            ]
        if kind is None and any(isinstance(v, dict) for v in values if v is not None):
            values = [json.dumps(v, sort_keys=True) if v is not None else None for v in values]
            kind = S
        if kind is None and all(v is None for v in values):
            # Unknown extension columns have no justified type yet. Null is
            # promoted by a later non-null shard during DuckDB union_by_name.
            kind = pa.null()
        arrays.append(pa.array(values, type=kind))
    table = pa.Table.from_arrays(arrays, names=columns)
    return table.replace_schema_metadata({b"phageannoscope_schema": SCHEMA_VERSION.encode()})


class ResultWriter:
    """One shard per table per bounded batch, not one file per CDS/genome."""

    def __init__(self, root, batch_genomes=20):
        self.root = Path(root)
        # Regenerate only our named shards. Interrupted/retried runs may have fewer batches.
        for name in SCHEMAS:
            for shard in (self.root / name).glob("part-*.parquet"):
                shard.unlink()
        self.batch_genomes = batch_genomes
        self.buffer = {name: [] for name in SCHEMAS}
        self.counts = {name: 0 for name in SCHEMAS}
        self.batch = 0
        self.pending = 0
        # Dedup keys live on disk; millions of proteins/models must not become Python sets.
        self.root.mkdir(parents=True, exist_ok=True)
        self.dedup = sqlite3.connect(str(self.root / "dedup.sqlite3"))
        self.dedup.execute(
            "CREATE TABLE IF NOT EXISTS seen (name TEXT, identity TEXT, PRIMARY KEY(name,identity))"
        )
        self.dedup.execute("DELETE FROM seen")
        self.unique = {
            "genomes": ("genome_id",),
            "models": ("model_id",),
            "model_segments": ("model_id", "segment_order"),
            "proteins": ("protein_sha256",),
            "locus_models": ("locus_id", "model_id"),
        }

    def add(self, tables):
        for name, rows in tables.items():
            if name not in self.buffer:
                continue
            for row in rows:
                if name in self.unique:
                    key = json.dumps([row.get(k) for k in self.unique[name]])
                    cursor = self.dedup.execute(
                        "INSERT OR IGNORE INTO seen VALUES (?,?)", (name, key)
                    )
                    if not cursor.rowcount:
                        continue
                self.buffer[name].append(row)
                self.counts[name] += 1
        self.pending += 1
        if self.pending >= self.batch_genomes:
            self.flush()

    def flush(self):
        for name, rows in self.buffer.items():
            if rows or self.batch == 0:
                target = self.root / name / f"part-{self.batch:05d}.parquet"
                target.parent.mkdir(parents=True, exist_ok=True)
                pq.write_table(arrow_table(name, rows), target, compression="zstd")
            rows.clear()
        self.batch += 1
        self.pending = 0
        self.dedup.commit()

    def close(self):
        if self.pending or self.batch == 0:
            self.flush()
        self.dedup.close()
        return self.counts


def sql_quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def build_db(data_root, out):
    root, out = Path(data_root).resolve(), Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with duckdb.connect(str(out)) as db:
        for name in SCHEMAS:
            files = sorted((root / name).glob("*.parquet"))
            if not files:
                raise ValueError(f"Missing Parquet table: {name}")
            glob = (root / name / "*.parquet").as_posix()
            db.execute(
                f'CREATE OR REPLACE VIEW "{name}" AS SELECT * FROM '
                f"read_parquet({sql_quote(glob)}, union_by_name=true)"
            )
        db.execute("CREATE OR REPLACE TABLE schema_info AS SELECT ? AS version", [SCHEMA_VERSION])
    return str(out)


def export_proteins(db_path, out):
    with (
        duckdb.connect(str(db_path), read_only=True) as db,
        Path(out).open("w", encoding="ascii") as handle,
    ):
        cursor = db.execute(
            "SELECT protein_sha256, amino_acid_sequence FROM proteins ORDER BY protein_sha256"
        )
        while batch := cursor.fetchmany(1000):
            for sha, seq in batch:
                handle.write(f">{sha}\n")
                for offset in range(0, len(seq), 80):
                    handle.write(seq[offset : offset + 80] + "\n")
