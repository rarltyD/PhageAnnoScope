"""Genome-sized jobs, deterministic checkpoints and bounded result assembly."""

import concurrent.futures
import csv
import json
import logging
import time
from pathlib import Path

from .callers.runner import run_caller
from .config import resolve_external_tools
from .genbank_io import read_genbank
from .matching import build_loci, build_models, compare_calls, genome_metrics
from .proteins import build_catalogue
from .resources import BackgroundResourceMonitor
from .storage import ResultWriter, build_db, export_proteins
from .util import digest, file_sha, hardware, utcnow, versions, write_json

LOG = logging.getLogger(__name__)


def apply_primary_code_policy(calls, runs):
    """Protocol v1 holds the entire primary pair for alternative code/recoding."""
    reasons = []
    if any(c.get("transl_table_effective") not in (1, 11) for c in calls):
        reasons.append("nonstandard_or_invalid_translation_code")
    if any(r.get("selected_translation_table") not in (None, 1, 11) for r in runs):
        reasons.append("caller_selected_alternative_code")
    # Amendment A2: a programmed frameshift (/ribosomal_slippage, recorded as
    # exception="ribosomal_slippage" by genbank_io) excludes that CDS only; it does
    # not change the genome's stop-codon system, so it does not hold the genome.
    if any(
        c.get("source") == "genbank"
        and (
            (c.get("exception") and c.get("exception") != "ribosomal_slippage")
            or c.get("transl_except")
        )
        for c in calls
    ):
        reasons.append("baseline_recoding_review_required")
    if reasons:
        for run in runs:
            run.update(comparison_status="not_comparable", comparison_reason=";".join(reasons))
        for call in calls:
            prior = list(call.get("primary_exclusion_reasons") or [])
            call["primary_exclusion_reasons"] = list(dict.fromkeys(prior + reasons))
            call["primary_exclusion_reason"] = ";".join(call["primary_exclusion_reasons"])
            call["eligible_primary"] = False
    return not reasons


def mark_rotation_sensitivity(calls, runs, topology):
    if topology != "circular":
        return
    successful = {r["source_run_id"] for r in runs if r["status"] in ("success", "cache_hit")}
    for source in ("phanotate", "pyrodigal_gv"):
        tracks = [
            [c for c in calls if c["source"] == source and bool(c.get("rotation_shift")) == rotated]
            for rotated in (False, True)
        ]
        source_runs = [r for r in runs if r["source"] == source]
        successful_shifts = {
            r["rotation_shift"] for r in source_runs if r["source_run_id"] in successful
        }
        both_ok = len(source_runs) == 2 and len(successful_shifts) == 2 and 0 in successful_shifts
        keys = [{c["model_id"] for c in track} for track in tracks]
        for index, track in enumerate(tracks):
            for call in track:
                call["rotation_sensitive"] = both_ok and call["model_id"] not in keys[1 - index]
                reason = (
                    "rotation_check_incomplete"
                    if not both_ok
                    else "rotation_sensitive"
                    if call["rotation_sensitive"]
                    else "rotation_track"
                    if index
                    else None
                )
                if reason:
                    call["eligible_primary"] = False
                    prior = list(call.get("primary_exclusion_reasons") or [])
                    if (
                        call.get("primary_exclusion_reason")
                        and call["primary_exclusion_reason"] not in prior
                    ):
                        prior.append(call["primary_exclusion_reason"])
                    if reason not in prior:
                        prior.append(reason)
                    call["primary_exclusion_reasons"] = prior
                    call["primary_exclusion_reason"] = ";".join(prior)


def process_genome(row, config, raw_root):
    path = Path(row["genbank_path"])
    raw_sha = file_sha(path)
    if row.get("genbank_sha256") and raw_sha != row["genbank_sha256"]:
        raise ValueError("Raw GenBank checksum mismatch")
    baseline_id = "run_genbank_" + digest(
        {"accession_version": row["accession_version"], "raw_sha": raw_sha}
    )
    imported = read_genbank(
        path, source_run_id=baseline_id, accession_version=row["accession_version"]
    )
    genomes = [g for g in imported["genomes"] if g["accession_version"] == row["accession_version"]]
    if len(genomes) != 1:
        raise ValueError("Requested accession.version is not unique in GenBank input")
    genome = genomes[0]
    if row.get("sequence_sha256") and row["sequence_sha256"] != genome["sequence_sha256"]:
        raise ValueError("Manifest sequence checksum mismatch")
    for key in (
        "host_domain",
        "host_domain_evidence",
        "host_domain_status",
        "host_taxonomy_rank",
        "host_evidence_level",
        "host_taxonomy_evidence_json",
        "host_evidence_data_root",
        "checkv_quality",
        "baltimore_group",
        "relatedness_cluster",
        "exact_sequence_duplicate_group",
        "selection_status",
    ):
        genome[key] = row.get(key) or None
    genome["host_taxid"] = int(row["host_taxid"]) if row.get("host_taxid") else None
    genome["in_primary_analysis"] = False
    calls = [c for c in imported["cds_calls"] if c["genome_id"] == genome["genome_id"]]
    features = [f for f in imported["features"] if f["genome_id"] == genome["genome_id"]]
    runs = [
        {
            "source_run_id": baseline_id,
            "genome_id": genome["genome_id"],
            "source": "genbank",
            "algorithm_family": "genbank",
            "status": "success",
            "rotation_shift": 0,
            "input_sha256": genome["sequence_sha256"],
            "tool_version": versions()["biopython"],
            "raw_genbank_sha256": raw_sha,
            "raw_path": str(path),
            "cds_count": len(calls),
        }
    ]
    codes = sorted(
        {c["transl_table_effective"] for c in calls if c.get("transl_table_effective") is not None}
    )
    for source in ("phanotate", "pyrodigal_gv"):
        for shift in (0, genome["length_bp"] // 2) if genome["topology"] == "circular" else (0,):
            source_run, predictions = run_caller(genome, source, config, raw_root, shift, codes)
            runs.append(source_run)
            calls.extend(predictions)
            fields = (
                "feature_id",
                "genome_id",
                "accession_version",
                "source",
                "source_run_id",
                "feature_ordinal",
                "feature_type",
                "type",
                "raw_feature_id",
                "strand",
                "segments",
                "original_location",
                "location_operator",
                "qualifiers_json",
                "partial_5p",
                "partial_3p",
                "location_issues",
            )
            features.extend({key: call.get(key) for key in fields} for call in predictions)
    mark_rotation_sensitivity(calls, runs, genome["topology"])
    compatible = apply_primary_code_policy(calls, runs)
    matches = compare_calls(calls, runs, config["partial_overlap_reciprocal_primary"])
    loci = build_loci(calls, runs)
    models = build_models(calls)
    proteins = build_catalogue(calls)
    qc = [q for q in imported["qc_events"] if q.get("genome_id") == genome["genome_id"]]
    qc.extend(
        {
            "genome_id": genome["genome_id"],
            "source_run_id": r["source_run_id"],
            "severity": "warning",
            "code": r["status"],
            "message": r.get("error", r["status"]),
        }
        for r in runs
        if r["status"] not in ("success", "cache_hit")
    )
    complete = all(r["status"] in ("success", "cache_hit") for r in runs)
    cohort_ok = (
        row.get("eligible_for_pilot", "").lower() == "true"
        or row.get("in_primary_analysis", "").lower() == "true"
    )
    genome["in_primary_analysis"] = (
        complete and cohort_ok and compatible and not config.get("diagnostic", False)
    )
    genome["primary_exclusion_reason"] = (
        ";".join(
            filter(
                None,
                (
                    "diagnostic_only" if config.get("diagnostic", False) else None,
                    "cohort_not_eligible" if not cohort_ok else None,
                    "caller_not_successful" if not complete else None,
                    "genetic_code_or_recoding_review" if not compatible else None,
                ),
            )
        )
        or None
    )
    for locus in loci["stop_loci"]:
        locus["geometry_eligible_primary"] = locus["eligible_primary"]
        locus["cohort_eligible_primary"] = genome["in_primary_analysis"]
        if not genome["in_primary_analysis"]:
            locus["eligible_primary"] = False
            prior = locus.get("primary_exclusion_reason")
            locus["primary_exclusion_reason"] = ";".join(
                filter(None, (prior, "genome_not_in_primary_cohort"))
            )
    return {
        "genomes": [genome],
        "source_runs": runs,
        "features": features,
        "cds_calls": calls,
        "matches": matches,
        **loci,
        **models,
        **proteins,
        "qc_events": qc,
        "genome_metrics": genome_metrics(calls, matches, loci, runs),
        "selection_events": [
            {
                "genome_id": genome["genome_id"],
                "accession": genome["accession"],
                "stage": "caller_qc",
                "decision": "eligible_primary" if genome["in_primary_analysis"] else "hold_primary",
                "primary_reason": genome["primary_exclusion_reason"] or "caller_success",
                "evidence": json.dumps(
                    {r["source"] + str(r["rotation_shift"]): r["status"] for r in runs}
                ),
            }
        ],
    }


def run_pipeline(manifest, config, out, workers=1, diagnostic=False, raw_root="data/raw/callers"):
    config = resolve_external_tools(config)
    config["diagnostic"] = bool(diagnostic)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "checkpoints"
    work.mkdir(exist_ok=True)
    if workers < 1 or workers > config["max_workers"]:
        raise ValueError("workers must be between 1 and configured max_workers")
    rows = []
    with Path(manifest).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if not row.get("genbank_path") or not row.get("accession_version"):
                raise ValueError(
                    "Run manifest requires resolved accession_version and genbank_path"
                )
            if not diagnostic and not (
                row.get("eligible_for_pilot", "").lower() == "true"
                or row.get("in_primary_analysis", "").lower() == "true"
            ):
                raise ValueError(
                    "Unvalidated cohort row; use --diagnostic for a separate technical run"
                )
            rows.append(row)
    if not rows or len({r["accession_version"] for r in rows}) != len(rows):
        raise ValueError("Run manifest must contain unique accession.version rows")
    code_hash = digest(
        {
            str(p.relative_to(Path(__file__).parent)): file_sha(p)
            for p in Path(__file__).parent.rglob("*.py")
        }
    )
    identity = digest(
        {
            "manifest_sha256": file_sha(manifest),
            "config": config,
            "code_sha256": code_hash,
            "versions": versions(),
            "diagnostic": diagnostic,
        }
    )
    run_meta = out / "run_config.json"
    if (
        run_meta.exists()
        and json.loads(run_meta.read_text(encoding="utf-8"))["identity"] != identity
    ):
        raise ValueError(
            "Output directory belongs to a different run/config/code; choose a new --out"
        )
    write_json(
        run_meta,
        {
            "identity": identity,
            "config": config,
            "versions": versions(),
            "code_sha256": code_hash,
            "manifest": str(Path(manifest).resolve()),
            "manifest_sha256": file_sha(manifest),
            "diagnostic": diagnostic,
        },
    )
    started, wall_start, cpu_start = utcnow(), time.perf_counter(), time.process_time()
    resource_monitor = BackgroundResourceMonitor(disk_path=out.resolve())
    resource_monitor.start()
    writer = ResultWriter(out / "parquet", config["storage_batch_genomes"])
    completed = 0
    failed = []
    checkpoint_hits = 0

    def job(row):
        key = digest({"row": row, "raw_sha": file_sha(row["genbank_path"]), "identity": identity})
        checkpoint = work / f"{key}.json"
        checksum = work / f"{key}.sha256"
        if (
            config.get("use_caller_cache", True)
            and checkpoint.exists()
            and checksum.exists()
            and checksum.read_text() == file_sha(checkpoint)
        ):
            result = json.loads(checkpoint.read_text(encoding="utf-8"))
            # Retry incomplete callers, including dependencies installed after a failed attempt.
            if all(r["status"] in ("success", "cache_hit") for r in result["source_runs"]):
                for source_run in result["source_runs"]:
                    if source_run["source"] != "genbank":
                        source_run["cache_origin_status"] = source_run.get(
                            "cache_origin_status", source_run["status"]
                        )
                        source_run["status"] = "cache_hit"
                        source_run["cache_read_at_utc"] = utcnow()
                return result, True
        result = process_genome(row, config, raw_root)
        write_json(checkpoint, result)
        checksum.write_text(file_sha(checkpoint), encoding="ascii")
        return result, False

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        # Bounded in-flight queue: avoid retaining a complete collection of CDS in Futures.
        iterator = iter(rows)
        pending = {}
        for _ in range(min(workers, len(rows))):
            row = next(iterator)
            pending[pool.submit(job, row)] = row
        while pending:
            done, _ = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                row = pending.pop(future)
                try:
                    result, hit = future.result()
                except Exception as error:
                    LOG.exception("Genome failed: %s", row["accession_version"])
                    failed.append(
                        {
                            "accession_version": row["accession_version"],
                            "error": f"{type(error).__name__}: {error}",
                        }
                    )
                else:
                    # Storage failures are fatal; never describe a partially written genome as skipped.
                    writer.add(result)
                    checkpoint_hits += int(hit)
                    completed += 1
                    LOG.info("Processed %s (%d/%d)", row["accession_version"], completed, len(rows))
                next_row = next(iterator, None)
                if next_row is not None:
                    pending[pool.submit(job, next_row)] = next_row
    counts = writer.close()
    db = out / "phageannoscope.duckdb"
    build_db(out / "parquet", db)
    export_proteins(db, out / "proteins.faa")
    import duckdb

    with duckdb.connect(str(db), read_only=True) as connection:
        status_counts = dict(
            connection.execute(
                "SELECT status, count(*) FROM source_runs GROUP BY status"
            ).fetchall()
        )
        primary_n = connection.execute(
            "SELECT count(*) FROM genomes WHERE in_primary_analysis"
        ).fetchone()[0]
        caller_complete_n = connection.execute(
            "SELECT count(*) FROM (SELECT genome_id FROM source_runs WHERE source IN ('phanotate','pyrodigal_gv') GROUP BY genome_id HAVING count(DISTINCT source)=2 AND bool_and(status IN ('success','cache_hit')))"
        ).fetchone()[0]
        fresh_cpu, fresh_runs, observed_cpu_runs = connection.execute(
            "SELECT sum(cpu_seconds), count(*), count(cpu_seconds) FROM source_runs "
            "WHERE source!='genbank' AND status='success'"
        ).fetchone()
    resources = resource_monitor.stop()
    fresh_cpu_reason = None
    if not fresh_runs:
        fresh_cpu_reason = "no_fresh_successful_callers"
    elif observed_cpu_runs != fresh_runs:
        fresh_cpu = None
        fresh_cpu_reason = "fresh_caller_cpu_measurements_incomplete"
    orchestrator_cpu = time.process_time() - cpu_start
    summary = {
        "run_id": out.name,
        "identity": identity,
        "started_at_utc": started,
        "ended_at_utc": utcnow(),
        "diagnostic": diagnostic,
        "planned_n": len(rows),
        "imported_n": completed,
        "caller_completed_n": caller_complete_n,
        "eligible_primary_n": primary_n,
        "failed_imports": failed,
        "checkpoint_hits": checkpoint_hits,
        "source_status_counts": status_counts,
        "table_counts": counts,
        "wall_seconds": time.perf_counter() - wall_start,
        "orchestrator_cpu_seconds": orchestrator_cpu,
        "fresh_caller_cpu_seconds": fresh_cpu,
        "fresh_caller_cpu_unavailable_reason": fresh_cpu_reason,
        "total_cpu_observed_seconds": None
        if fresh_runs and fresh_cpu is None
        else orchestrator_cpu + (fresh_cpu if fresh_cpu is not None else 0),
        "total_cpu_unavailable_reason": fresh_cpu_reason
        if fresh_runs and fresh_cpu is None
        else None,
        "peak_rss_process_tree_bytes": resources["peak_rss_bytes"],
        "peak_rss_unavailable_reason": resources["peak_rss_unavailable_reason"],
        "resource_monitor_status": resources["resource_monitor_status"],
        "disk_free_start_bytes": resources["disk_free_start_bytes"],
        "disk_free_min_bytes": resources["disk_free_min_bytes"],
        "disk_peak_delta_bytes": resources["disk_peak_delta_bytes"],
        "disk_unavailable_reason": resources["disk_unavailable_reason"],
        "disk_measurement_scope": resources["disk_measurement_scope"],
        "resource_method": "Windows process-tree RSS sampled 50ms excludes Linux VM; native WSL RSS separate. Total CPU adds fresh native WSL child CPU and sampled Windows caller CPU (lower bound) to orchestrator CPU; cache historical CPU excluded.",
        "workers": workers,
        "hardware": hardware(),
        "output_bytes": sum(p.stat().st_size for p in out.rglob("*") if p.is_file()),
        "pilot_status": "diagnostic_only"
        if diagnostic
        else ("completed" if primary_n == len(rows) else "incomplete"),
    }
    write_json(out / "run_summary.json", summary)
    with Path(manifest).open(encoding="utf-8-sig") as source:
        (out / "input_manifest.csv").write_text(source.read(), encoding="utf-8")
    return summary
