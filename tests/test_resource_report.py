"""Hand-calculated synthetic timing observations, never biological benchmarks."""

import csv
import json

import duckdb
import pytest

from phageannoscope.benchmark import resource_report


def example(tmp_path, status="success", timeout=False):
    run = tmp_path / "run"
    run.mkdir()
    (run / "run_summary.json").write_text(
        json.dumps(
            {
                "wall_seconds": 50,
                "workers": 2,
                "output_bytes": 1000,
                "peak_rss_process_tree_bytes": None,
            }
        )
    )
    with duckdb.connect(str(run / "phageannoscope.duckdb")) as db:
        db.execute("CREATE TABLE genomes(genome_id VARCHAR,length_bp BIGINT)")
        db.execute("INSERT INTO genomes VALUES ('synthetic',9000)")
        db.execute(
            "CREATE TABLE source_runs(genome_id VARCHAR,source VARCHAR,status VARCHAR,rotation_shift BIGINT,wall_seconds DOUBLE,cpu_seconds DOUBLE,peak_rss_bytes BIGINT)"
        )
        for source in ("phanotate", "pyrodigal_gv"):
            db.executemany(
                "INSERT INTO source_runs VALUES (?,?,?,?,?,?,?)",
                [("synthetic", source, status, 0, 10, None, None)] * 5,
            )
        if timeout:
            db.execute(
                "INSERT INTO source_runs VALUES ('synthetic','phanotate','timeout',0,300,NULL,NULL)"
            )
    manifest = tmp_path / "population.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["ncbi_length_bp", "eligible_for_pilot", "topology"]
        )
        writer.writeheader()
        writer.writerows(
            [{"ncbi_length_bp": 9000, "eligible_for_pilot": "true", "topology": "linear"}] * 2
        )
    return run, manifest


def test_stratified_forecast_with_observed_parallelism_and_calendar_reserve(tmp_path):
    run, population = example(tmp_path)
    report = resource_report(run, population, tmp_path / "report.json", available_hours=1)
    assert report["observed_effective_parallelism"] == 2
    assert report["forecast_compute_wall_seconds"] == 20
    assert report["required_calendar_hours_with_reserve"] == pytest.approx(20 / 3600 / 0.75)
    assert report["fits_supplied_time_budget"] is True
    assert report["main_run_launched"] is False
    observed = [c for c in report["cells"] if c["fresh_attempt_n"]]
    assert all(c["cpu_seconds_sum_observed"] is None for c in observed)
    assert all(c["max_native_or_process_rss_bytes"] is None for c in observed)


def test_historical_cache_times_cannot_supply_forecast(tmp_path):
    run, population = example(tmp_path, status="cache_hit")
    report = resource_report(run, population, tmp_path / "report.json")
    assert report["forecast_compute_wall_seconds"] is None
    assert report["fits_supplied_time_budget"] is None
    assert len(report["blocking_forecast_reasons"]) == 2
    assert all(c["fresh_attempt_n"] == 0 for c in report["cells"])


def test_timeout_is_in_cost_and_not_discarded_as_fast_success(tmp_path):
    run, population = example(tmp_path, timeout=True)
    report = resource_report(run, population, tmp_path / "report.json")
    cell = next(
        c for c in report["cells"] if c["source"] == "phanotate" and c["stratum"] == "lt10k"
    )
    assert cell["status_counts"] == {"success": 5, "timeout": 1}
    assert cell["fresh_r0_success_n"] == 5 and cell["fresh_attempt_n"] == 6
    assert report["forecast_compute_wall_seconds"] > 20


@pytest.mark.parametrize("hours", [0, -1, float("nan"), float("inf")])
def test_invalid_time_budget_rejected(tmp_path, hours):
    run, population = example(tmp_path)
    with pytest.raises(ValueError, match="available_hours"):
        resource_report(run, population, tmp_path / "report.json", available_hours=hours)
