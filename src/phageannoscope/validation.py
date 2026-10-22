"""Independent output checks using database constraints and sequence hashes."""

from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from .storage import SCHEMAS
from .util import digest, utcnow, write_json

PRIMARY_KEYS = {
    "genomes": ["genome_id"],
    "source_runs": ["source_run_id"],
    "features": ["feature_id"],
    "cds_calls": ["call_id"],
    "models": ["model_id"],
    "model_sources": ["call_id"],
    "model_segments": ["model_id", "segment_order"],
    "stop_loci": ["locus_id"],
    "locus_models": ["locus_id", "model_id"],
    "proteins": ["protein_sha256"],
    "protein_origins": ["protein_sha256", "call_id"],
    "matches": ["match_id"],
}


def audit_support_gate(db):
    """Recompute support from source links, independently of the locus builder."""
    columns = {row[0] for row in db.execute("DESCRIBE stop_loci").fetchall()}
    required = {
        "caller_support_raw",
        "caller_support_eligible",
        "caller_support_mask_raw",
        "caller_support_mask_eligible",
        "support_assessable_primary",
    }
    errors = []
    if missing := required - columns:
        errors.append("Support fields missing: " + ", ".join(sorted(missing)))
    broken_links = db.execute("""SELECT count(*) FROM locus_models lm
        JOIN stop_loci l USING(locus_id)
        CROSS JOIN UNNEST(lm.call_ids) AS u(cid)
        JOIN cds_calls d ON d.call_id = u.cid
        WHERE lm.model_id IS DISTINCT FROM d.model_id
           OR l.genome_id IS DISTINCT FROM d.genome_id
           OR l.sequence_sha256 IS DISTINCT FROM d.sequence_sha256
           OR l.stop_triplet IS DISTINCT FROM d.stop_triplet
           OR l.strand IS DISTINCT FROM d.strand
           OR l.code_class IS DISTINCT FROM d.code_class""").fetchone()[0]
    if broken_links:
        errors.append(f"Locus links disagree with defining call geometry/identity: {broken_links}")
    raw = "l.caller_support_raw" if "caller_support_raw" in columns else "l.caller_support"
    admitted = (
        "l.geometry_eligible_primary"
        if "geometry_eligible_primary" in columns
        else "l.eligible_primary"
    )
    # DISTINCT avoids multiplying a call when malformed redundant links exist;
    # those links are separately checked by the PK/FK checks.
    db.execute("""CREATE OR REPLACE TEMP VIEW support_call_evidence AS
        SELECT DISTINCT lm.locus_id, d.call_id, d.model_id,
            CASE WHEN replace(lower(coalesce(r.algorithm_family,r.source)),'-','_')
                    IN ('pyrodigal_gv','pyrodigal','prodigal','prokka')
                 THEN 'pyrodigal_gv'
                 ELSE replace(lower(coalesce(r.algorithm_family,r.source)),'-','_') END AS family,
            d.source IN ('phanotate','pyrodigal_gv') AND coalesce(d.rotation_shift,0)=0
                AND r.source IN ('phanotate','pyrodigal_gv')
                AND r.status IN ('success','cache_hit') AS raw_contribution,
            coalesce(d.eligible_primary,false) AND NOT coalesce(d.wraps_origin,false)
                AND NOT coalesce(d.rotation_sensitive,false)
                AND r.comparison_status IS DISTINCT FROM 'not_comparable' AS admissible
        FROM locus_models lm CROSS JOIN UNNEST(lm.call_ids) AS u(cid)
        JOIN cds_calls d ON d.call_id = u.cid
        JOIN source_runs r ON d.source_run_id=r.source_run_id""")
    db.execute("""CREATE OR REPLACE TEMP VIEW support_observed AS
        SELECT locus_id,
            count(DISTINCT CASE WHEN raw_contribution AND family IN ('phanotate','pyrodigal_gv') THEN family END) AS raw_count,
            count(DISTINCT CASE WHEN raw_contribution AND admissible AND family IN ('phanotate','pyrodigal_gv') THEN family END) AS eligible_count,
            coalesce(bit_or(CASE WHEN raw_contribution THEN CASE family WHEN 'phanotate' THEN 1 WHEN 'pyrodigal_gv' THEN 2 ELSE 0 END ELSE 0 END),0) AS raw_mask,
            coalesce(bit_or(CASE WHEN raw_contribution AND admissible THEN CASE family WHEN 'phanotate' THEN 1 WHEN 'pyrodigal_gv' THEN 2 ELSE 0 END ELSE 0 END),0) AS eligible_mask,
            count(*) FILTER (WHERE raw_contribution AND family IN ('phanotate','pyrodigal_gv') AND NOT admissible) AS blocked_count
        FROM support_call_evidence GROUP BY locus_id""")
    db.execute("""CREATE OR REPLACE TEMP VIEW support_expected AS
        SELECT l.locus_id, coalesce(o.raw_count,0) AS raw_count,
            coalesce(o.eligible_count,0) AS eligible_count,
            coalesce(o.raw_mask,0) AS raw_mask, coalesce(o.eligible_mask,0) AS eligible_mask,
            coalesce(o.blocked_count,0) AS blocked_count,
            coalesce(o.raw_count,0)>0 AND coalesce(o.eligible_count,0)>0
            AND coalesce(o.blocked_count,0)=0 AND l.code_class='standard_stop'
            AND EXISTS (SELECT 1 FROM source_runs r WHERE r.genome_id=l.genome_id AND r.source='phanotate'
                AND coalesce(r.rotation_shift,0)=0 AND r.status IN ('success','cache_hit')
                AND r.comparison_status IS DISTINCT FROM 'not_comparable')
            AND EXISTS (SELECT 1 FROM source_runs r WHERE r.genome_id=l.genome_id AND r.source='pyrodigal_gv'
                AND coalesce(r.rotation_shift,0)=0 AND r.status IN ('success','cache_hit')
                AND r.comparison_status IS DISTINCT FROM 'not_comparable') AS assessable
        FROM stop_loci l LEFT JOIN support_observed o USING(locus_id)""")
    bad_admitted = db.execute(f"""SELECT l.locus_id,{raw},e.eligible_count,e.blocked_count
        FROM stop_loci l JOIN support_expected e USING(locus_id)
        WHERE ({admitted} OR l.eligible_primary) AND NOT coalesce(e.assessable,false)
        ORDER BY l.locus_id""").fetchall()
    if bad_admitted:
        errors.append(f"Primary support admission contradicts defining calls: {len(bad_admitted)}")
    if not missing:
        n = db.execute("""SELECT count(*) FROM stop_loci l JOIN support_expected e USING(locus_id)
            WHERE l.caller_support_raw IS DISTINCT FROM e.raw_count
               OR l.caller_support_eligible IS DISTINCT FROM e.eligible_count
               OR l.caller_support_mask_raw IS DISTINCT FROM e.raw_mask
               OR l.caller_support_mask_eligible IS DISTINCT FROM e.eligible_mask
               OR l.caller_support IS DISTINCT FROM e.raw_count
               OR l.caller_support_mask IS DISTINCT FROM e.raw_mask
               OR l.support_assessable_primary IS DISTINCT FROM e.assessable""").fetchone()[0]
        if n:
            errors.append(f"Support counts/masks/assessability disagree with source calls: {n}")
    n = db.execute("""SELECT count(*) FROM stop_loci l
        WHERE l.representative_call_id IS NOT NULL AND NOT EXISTS
          (SELECT 1 FROM support_call_evidence e WHERE e.locus_id=l.locus_id
            AND e.call_id=l.representative_call_id AND e.model_id=l.representative_model_id
            AND e.raw_contribution AND e.admissible AND e.family IN ('phanotate','pyrodigal_gv'))""").fetchone()[
        0
    ]
    if n:
        errors.append(f"Representative lies outside eligible support domain: {n}")
    return {
        "status": "failed" if errors else "passed",
        "errors": errors,
        "n_admission_inconsistencies": len(bad_admitted),
        "admission_cases": [
            {
                "locus_id": row[0],
                "caller_support_raw": row[1],
                "caller_support_eligible": row[2],
                "excluded_calls": row[3],
            }
            for row in bad_admitted[:20]
        ],
    }


def validate_results(run, out=None):
    run = Path(run)
    errors, counts, checks = [], {}, []
    with duckdb.connect(str(run / "phageannoscope.duckdb"), read_only=True) as db:
        for table in SCHEMAS:
            counts[table] = db.execute(f'SELECT count(*) FROM "{table}"').fetchone()[0]
            for path in (run / "parquet" / table).glob("*.parquet"):
                if pq.read_schema(path).metadata.get(b"phageannoscope_schema") != b"1.0":
                    errors.append(f"Schema version mismatch: {path.name}")
        checks.append("schema version and all table counts")
        for table, columns in PRIMARY_KEYS.items():
            grouped = ",".join(f'"{c}"' for c in columns)
            nulls = " OR ".join(f'"{c}" IS NULL' for c in columns)
            if db.execute(f'SELECT count(*) FROM "{table}" WHERE {nulls}').fetchone()[0]:
                errors.append(f"Null primary key in {table}")
            if db.execute(
                f'SELECT count(*) FROM (SELECT {grouped} FROM "{table}" GROUP BY {grouped} HAVING count(*)>1)'
            ).fetchone()[0]:
                errors.append(f"Duplicate primary key in {table}")
        checks.append("primary keys and duplicate associations")
        references = [
            ("source_runs", "genome_id", "genomes", "genome_id"),
            ("features", "genome_id", "genomes", "genome_id"),
            ("features", "source_run_id", "source_runs", "source_run_id"),
            ("cds_calls", "source_run_id", "source_runs", "source_run_id"),
            ("cds_calls", "genome_id", "genomes", "genome_id"),
            ("cds_calls", "feature_id", "features", "feature_id"),
            ("cds_calls", "model_id", "models", "model_id"),
            ("model_sources", "call_id", "cds_calls", "call_id"),
            ("model_sources", "model_id", "models", "model_id"),
            ("model_sources", "source_run_id", "source_runs", "source_run_id"),
            ("models", "genome_id", "genomes", "genome_id"),
            ("model_segments", "model_id", "models", "model_id"),
            ("locus_models", "model_id", "models", "model_id"),
            ("locus_models", "locus_id", "stop_loci", "locus_id"),
            ("stop_loci", "representative_model_id", "models", "model_id"),
            ("stop_loci", "genome_id", "genomes", "genome_id"),
            ("protein_origins", "protein_sha256", "proteins", "protein_sha256"),
            ("protein_origins", "call_id", "cds_calls", "call_id"),
            ("protein_origins", "feature_id", "features", "feature_id"),
            ("protein_origins", "model_id", "models", "model_id"),
            ("protein_origins", "genome_id", "genomes", "genome_id"),
            ("protein_origins", "source_run_id", "source_runs", "source_run_id"),
            ("matches", "call_id_a", "cds_calls", "call_id"),
            ("matches", "call_id_b", "cds_calls", "call_id"),
        ]
        for child, foreign, parent, key in references:
            n = db.execute(
                f'SELECT count(*) FROM "{child}" c LEFT JOIN "{parent}" p ON c."{foreign}"=p."{key}" WHERE c."{foreign}" IS NOT NULL AND p."{key}" IS NULL'
            ).fetchone()[0]
            if n:
                errors.append(f"Broken FK {child}.{foreign} -> {parent}.{key}: {n}")
        checks.append("genome, model, run, protein and locus foreign keys")
        invalid = db.execute(
            "SELECT count(*) FROM model_segments s JOIN models m USING(model_id) JOIN genomes g USING(genome_id) WHERE s.start0<0 OR s.end0<=s.start0 OR s.end0>g.length_bp"
        ).fetchone()[0]
        if invalid:
            errors.append(f"Invalid segment boundaries: {invalid}")
        cursor = db.execute("SELECT protein_sha256,amino_acid_sequence,length_aa FROM proteins")
        while rows := cursor.fetchmany(1000):
            for sha, seq, length in rows:
                if digest(seq.encode("ascii")) != sha or len(seq) != length or "*" in seq:
                    errors.append(f"Invalid normalized protein: {sha}")
        checks.append("segment boundaries and normalized protein hashes")
        bad_support = db.execute(
            "SELECT count(*) FROM stop_loci l WHERE eligible_primary AND NOT EXISTS (SELECT 1 FROM source_runs r WHERE r.genome_id=l.genome_id AND r.source='phanotate' AND r.rotation_shift=0 AND r.status IN ('success','cache_hit')) OR eligible_primary AND NOT EXISTS (SELECT 1 FROM source_runs r WHERE r.genome_id=l.genome_id AND r.source='pyrodigal_gv' AND r.rotation_shift=0 AND r.status IN ('success','cache_hit'))"
        ).fetchone()[0]
        if bad_support:
            errors.append(f"Loci falsely eligible without both callers: {bad_support}")
        checks.append("primary eligibility requires both successful callers")
        support_report = audit_support_gate(db)
        errors.extend(support_report["errors"])
        checks.append(
            "support counts, masks, assessability and representative from all defining calls"
        )
    summary_path = run / "run_summary.json"
    if summary_path.exists():
        import json

        expected = json.loads(summary_path.read_text(encoding="utf-8"))
        if expected["table_counts"] != counts:
            errors.append("Run summary counts disagree with Parquet/DuckDB")
        if counts["genomes"] != expected["imported_n"]:
            errors.append("Imported genome IDs were lost")
    report = {
        "status": "passed" if not errors else "failed",
        "checked_at_utc": utcnow(),
        "table_counts": counts,
        "checks": checks,
        "errors": errors,
        "support_gate": support_report,
    }
    write_json(out or run / "validation_report.json", report)
    return report
