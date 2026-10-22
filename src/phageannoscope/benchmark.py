"""Export observed caller timings; no synthetic speed estimates."""

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import duckdb

from .manifest import BIN_NAMES, BINS, _true
from .util import utcnow, write_json


def export_benchmark(run, out=None):
    run = Path(run)
    target = Path(out or run / "benchmark.csv")
    with duckdb.connect(str(run / "phageannoscope.duckdb"), read_only=True) as db:
        # All statuses and rotations are exported. Cache times are historical, not fresh work.
        result = db.execute(
            "SELECT r.*, g.accession_version, g.length_bp FROM source_runs r JOIN genomes g USING(genome_id) ORDER BY genome_id,source,rotation_shift"
        )
        names = [c[0] for c in result.description]
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(names + ["fresh_successful_measurement", "timeout_censored"])
            while rows := result.fetchmany(1000):
                status_index = names.index("status")
                writer.writerows(
                    list(row)
                    + [
                        str(row[status_index] == "success").lower(),
                        str(row[status_index] == "timeout").lower(),
                    ]
                    for row in rows
                )
    return {"out": str(target)}


def _percentile(values, fraction=0.9):
    values = sorted(values)
    if not values:
        return None
    position = (len(values) - 1) * fraction
    left = int(position)
    right = min(left + 1, len(values) - 1)
    return values[left] + (values[right] - values[left]) * (position - left)


def resource_report(run, population_manifest, out, available_hours=None, reserve_fraction=0.25):
    """Cost a metadata/sequence-eligible frame using fresh observations by length.

    Timeouts/failures remain in measured costs. Cached and unexecuted runs do not
    supply fresh observations. A forecast requires >=5 fresh successful r0 calls
    in each populated source/length cell. Circular calls are budgeted twice.
    """
    if available_hours is not None and (not math.isfinite(available_hours) or available_hours <= 0):
        raise ValueError("available_hours must be finite and positive")
    if not 0 <= reserve_fraction < 1:
        raise ValueError("reserve_fraction must be in [0,1)")
    run = Path(run)
    summary = json.loads((run / "run_summary.json").read_text(encoding="utf-8"))
    population = {name: {"n": 0, "bp": 0, "circular_n": 0, "circular_bp": 0} for name in BIN_NAMES}
    with Path(population_manifest).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if not _true(row.get("eligible_for_pilot")):
                continue
            length = int(row.get("ncbi_length_bp") or row.get("length_bp") or 0)
            if length <= 0:
                raise ValueError("Forecast population lacks a validated positive length")
            name = BIN_NAMES[next(i for i in range(5) if BINS[i] <= length < BINS[i + 1])]
            population[name]["n"] += 1
            population[name]["bp"] += length
            if row.get("topology") == "circular":
                population[name]["circular_n"] += 1
                population[name]["circular_bp"] += length
    groups = defaultdict(list)
    with duckdb.connect(str(run / "phageannoscope.duckdb"), read_only=True) as db:
        cursor = db.execute("""SELECT r.source,r.status,r.rotation_shift,r.wall_seconds,
            r.cpu_seconds,r.peak_rss_bytes,g.length_bp FROM source_runs r JOIN genomes g USING(genome_id)
            WHERE r.source IN ('phanotate','pyrodigal_gv')""")
        for source, status, rotation, wall, cpu, rss, length in cursor.fetchall():
            name = BIN_NAMES[next(i for i in range(5) if BINS[i] <= length < BINS[i + 1])]
            groups[source, name].append(
                {
                    "status": status,
                    "rotation": rotation,
                    "wall": wall,
                    "cpu": cpu,
                    "rss": rss,
                    "bp": length,
                }
            )
    cells, warnings = [], []
    total_cost = 0.0
    fresh_wall_total = 0.0
    for source in ("phanotate", "pyrodigal_gv"):
        for name in BIN_NAMES:
            rows = groups[source, name]
            fresh = [r for r in rows if r["status"] in {"success", "failed", "timeout"}]
            success = [r for r in fresh if r["status"] == "success" and not r["rotation"]]
            cost_rows = [
                r
                for r in fresh
                if r["wall"] is not None and math.isfinite(r["wall"]) and r["wall"] >= 0
            ]
            ratios = [r["wall"] / r["bp"] for r in cost_rows]
            fresh_wall_total += sum(r["wall"] for r in cost_rows)
            rate = max(_percentile(ratios) or 0, sum(ratios) / len(ratios)) if ratios else None
            if population[name]["n"] and (
                len(success) < 5 or rate is None or len(cost_rows) != len(fresh)
            ):
                warnings.append(f"insufficient_fresh_r0_observations:{source}:{name}")
            if population[name]["circular_n"] and sum(bool(r["rotation"]) for r in fresh) < 5:
                warnings.append(f"limited_rotation_cost_observations:{source}:{name}")
            cost = (
                rate * (population[name]["bp"] + population[name]["circular_bp"])
                if rate is not None
                else None
            )
            total_cost += cost or 0
            cpus, rss_values = (
                [r["cpu"] for r in fresh if r["cpu"] is not None],
                [r["rss"] for r in fresh if r["rss"] is not None],
            )
            cells.append(
                {
                    "source": source,
                    "stratum": name,
                    "status_counts": dict(Counter(r["status"] for r in rows)),
                    "fresh_attempt_n": len(fresh),
                    "fresh_r0_success_n": len(success),
                    "cpu_observed_n": len(cpus),
                    "peak_rss_observed_n": len(rss_values),
                    "cpu_seconds_sum_observed": sum(cpus) if cpus else None,
                    "max_native_or_process_rss_bytes": max(rss_values) if rss_values else None,
                    "conservative_seconds_per_bp": rate,
                    "population": population[name],
                    "estimated_source_work_seconds": cost,
                }
            )
    population_n = sum(v["n"] for v in population.values())
    blocking = [w for w in warnings if w.startswith("insufficient_")]
    if population_n == 0:
        blocking.append("no_final_sequence_metadata_eligible_population")
    elapsed = summary["wall_seconds"]
    effective = (
        min(summary["workers"], fresh_wall_total / elapsed)
        if elapsed > 0 and fresh_wall_total > 0
        else None
    )
    estimate = total_cost / effective if effective and not blocking else None
    calendar = estimate / 3600 / (1 - reserve_fraction) if estimate is not None else None
    report = {
        "at_utc": utcnow(),
        "run": str(run.resolve()),
        "population_n": population_n,
        "method": "per-source length strata; maximum of p90 and mean fresh attempt wall/bp including failed/timeouts; extra full call for circular genomes; observed effective parallelism",
        "cells": cells,
        "observed_effective_parallelism": effective,
        "forecast_compute_wall_seconds": estimate,
        "required_calendar_hours_with_reserve": calendar,
        "reserve_fraction": reserve_fraction,
        "available_hours_assumption": available_hours,
        "fits_supplied_time_budget": calendar <= available_hours
        if calendar is not None and available_hours is not None
        else None,
        "blocking_forecast_reasons": blocking,
        "warnings": warnings,
        "observed_run_peak_rss_bytes": summary.get("peak_rss_process_tree_bytes"),
        "observed_disk_peak_delta_bytes": summary.get("disk_peak_delta_bytes"),
        "disk_measurement_scope": summary.get("disk_measurement_scope"),
        "run_output_bytes": summary.get("output_bytes"),
        "limitations": [
            "Forecast is a resource estimate, not a guarantee or biological inclusion decision",
            "Windows process tree excludes native WSL VM memory; per-call native measurements are separate",
            "Rotation/reverse-complement duplicates and unresolved relatedness require resolution before final biological N/inference",
        ],
        "main_run_launched": False,
    }
    write_json(out, report)
    return report
