"""Unavailable telemetry must not alter subprocess or biological result status."""

import csv
import sys
import threading
from types import SimpleNamespace

import duckdb
import psutil
import pytest

from phageannoscope import resources
from phageannoscope.callers.runner import monitor
from phageannoscope.pipeline import run_pipeline
from phageannoscope.resources import BackgroundResourceMonitor, ProcessTreeSampler
from phageannoscope.util import hardware

PROBE_ERRORS = [
    OSError("/proc unavailable"),
    NotImplementedError("platform probe unavailable"),
    PermissionError("restricted platform probe"),
    psutil.NoSuchProcess(123),
    psutil.AccessDenied(123),
]


def failing(error):
    def fail(*args, **kwargs):
        raise error

    return fail


class MeasuredProcess:
    def children(self, recursive=False):
        return []

    def memory_info(self):
        return SimpleNamespace(rss=0)

    def cpu_times(self):
        return SimpleNamespace(user=0.0, system=0.0)


@pytest.mark.parametrize("error", PROBE_ERRORS)
@pytest.mark.parametrize("exit_code", [0, 7])
def test_missing_process_probe_preserves_exit_status_and_logs(
    tmp_path, monkeypatch, error, exit_code
):
    monkeypatch.setattr(resources.psutil, "Process", failing(error))
    stdout, stderr = tmp_path / "stdout", tmp_path / "stderr"
    result = monitor(
        [
            sys.executable,
            "-c",
            "import sys; print('kept stdout'); print('kept stderr', file=sys.stderr); "
            f"raise SystemExit({exit_code})",
        ],
        stdout,
        stderr,
        10,
    )
    assert result["status"] == ("success" if exit_code == 0 else "failed")
    assert result["exit_code"] == exit_code
    assert "kept stdout" in stdout.read_text()
    assert "kept stderr" in stderr.read_text()
    assert result["cpu_seconds"] is result["peak_rss_bytes"] is None
    assert type(error).__name__ in result["cpu_unavailable_reason"]
    assert type(error).__name__ in result["peak_rss_unavailable_reason"]


@pytest.mark.parametrize("error", PROBE_ERRORS)
def test_missing_process_probe_does_not_disable_timeout_kill(tmp_path, monkeypatch, error):
    monkeypatch.setattr(resources.psutil, "Process", failing(error))
    result = monitor(
        [sys.executable, "-c", "import time; print('started', flush=True); time.sleep(20)"],
        tmp_path / "stdout",
        tmp_path / "stderr",
        0.15,
    )
    assert result["status"] == "timeout"
    assert result["exit_code"] is not None
    assert result["wall_seconds"] < 5
    assert result["cpu_seconds"] is result["peak_rss_bytes"] is None


@pytest.mark.parametrize("error", PROBE_ERRORS)
def test_tree_probe_failure_does_not_disable_timeout(tmp_path, monkeypatch, error):
    process = MeasuredProcess()
    monkeypatch.setattr(process, "children", failing(error))
    monkeypatch.setattr(resources.psutil, "Process", lambda pid=None: process)
    result = monitor(
        [sys.executable, "-c", "import time; time.sleep(20)"],
        tmp_path / "stdout",
        tmp_path / "stderr",
        0.1,
    )
    assert result["status"] == "timeout"
    assert result["wall_seconds"] < 5
    assert result["cpu_seconds"] is result["peak_rss_bytes"] is None
    assert "process_tree" in result["peak_rss_unavailable_reason"]


def test_real_measured_zero_is_not_missing(monkeypatch):
    monkeypatch.setattr(resources.psutil, "Process", lambda pid=None: MeasuredProcess())
    sampler = ProcessTreeSampler()
    sampler.sample()
    result = sampler.result()
    assert result["cpu_seconds"] == 0.0
    assert result["peak_rss_bytes"] == 0
    assert result["cpu_unavailable_reason"] is result["peak_rss_unavailable_reason"] is None
    assert result["resource_monitor_status"] == "available"


@pytest.mark.parametrize("error", PROBE_ERRORS)
def test_later_probe_failure_invalidates_stale_values(monkeypatch, error):
    process = MeasuredProcess()
    monkeypatch.setattr(resources.psutil, "Process", lambda pid=None: process)
    sampler = ProcessTreeSampler()
    sampler.sample()
    monkeypatch.setattr(process, "memory_info", failing(error))
    monkeypatch.setattr(process, "cpu_times", failing(error))
    sampler.sample()
    result = sampler.result()
    assert result["cpu_seconds"] is result["peak_rss_bytes"] is None
    assert type(error).__name__ in result["cpu_unavailable_reason"]
    # Successful later probes cannot make an interrupted observation complete again.
    monkeypatch.setattr(process, "memory_info", lambda: SimpleNamespace(rss=12))
    monkeypatch.setattr(process, "cpu_times", lambda: SimpleNamespace(user=3.0, system=2.0))
    sampler.sample()
    assert sampler.result()["peak_rss_bytes"] is None


def test_cpu_and_memory_probes_fail_independently(monkeypatch):
    process = MeasuredProcess()
    monkeypatch.setattr(process, "cpu_times", failing(PermissionError("CPU denied")))
    monkeypatch.setattr(resources.psutil, "Process", lambda pid=None: process)
    sampler = ProcessTreeSampler()
    sampler.sample()
    result = sampler.result()
    assert result["cpu_seconds"] is None
    assert result["peak_rss_bytes"] == 0
    assert result["resource_monitor_status"] == "partial"


@pytest.mark.parametrize("error", PROBE_ERRORS)
def test_hardware_probe_errors_are_explicit_nulls(monkeypatch, error):
    monkeypatch.setattr(resources.psutil, "cpu_count", failing(error))
    monkeypatch.setattr(resources.psutil, "virtual_memory", failing(error))
    result = hardware()
    for field in ("cpu_logical", "cpu_physical", "ram_total_bytes"):
        assert result[field] is None
        assert type(error).__name__ in result[f"{field}_unavailable_reason"]


def test_unavailable_physical_cpu_does_not_erase_other_hardware(monkeypatch):
    monkeypatch.setattr(resources.psutil, "cpu_count", lambda logical=True: 4 if logical else None)
    monkeypatch.setattr(resources.psutil, "virtual_memory", lambda: SimpleNamespace(total=1024))
    result = hardware()
    assert result["cpu_logical"] == 4
    assert result["ram_total_bytes"] == 1024
    assert result["cpu_physical"] is None
    assert "no value" in result["cpu_physical_unavailable_reason"]


@pytest.mark.parametrize("stage", ["create", "start", "body"])
def test_background_monitor_guards_entire_thread_lifecycle(monkeypatch, stage):
    monkeypatch.setattr(resources.psutil, "Process", lambda pid=None: MeasuredProcess())
    background = BackgroundResourceMonitor()
    if stage == "body":
        monkeypatch.setattr(background.sampler, "sample", failing(OSError("/proc unavailable")))
    else:
        factory = (
            failing(NotImplementedError("thread unavailable"))
            if stage == "create"
            else lambda **kwargs: SimpleNamespace(start=failing(PermissionError("thread denied")))
        )
        monkeypatch.setattr(
            resources, "threading", SimpleNamespace(Event=threading.Event, Thread=factory)
        )
    background.start()
    result = background.stop()
    assert result["cpu_seconds"] is result["peak_rss_bytes"] is None
    assert "monitor_thread" in result["peak_rss_unavailable_reason"]


@pytest.mark.parametrize(
    "failure", ["process_init", "thread_create", "thread_start", "thread_body"]
)
def test_pipeline_writes_tables_despite_monitor_failure(tmp_path, monkeypatch, failure):
    from phageannoscope import pipeline

    if failure == "process_init":
        monkeypatch.setattr(resources.psutil, "Process", failing(OSError("/proc unavailable")))
    elif failure == "thread_body":
        monkeypatch.setattr(ProcessTreeSampler, "sample", failing(NotImplementedError("no RSS")))
    else:
        factory = (
            failing(OSError("thread creation denied"))
            if failure == "thread_create"
            else lambda **kwargs: SimpleNamespace(start=failing(OSError("thread start denied")))
        )
        monkeypatch.setattr(
            resources, "threading", SimpleNamespace(Event=threading.Event, Thread=factory)
        )
    raw = tmp_path / "synthetic-input.gb"
    raw.write_text("synthetic test input; biological parser is replaced in this unit test")
    manifest = tmp_path / "synthetic-manifest.csv"
    row = {
        "accession_version": "SYNTHETIC.1",
        "genbank_path": str(raw),
        "eligible_for_pilot": "true",
    }
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    monkeypatch.setattr(pipeline, "resolve_external_tools", lambda config: config)
    monkeypatch.setattr(
        pipeline,
        "process_genome",
        lambda *args: {
            "genomes": [{"genome_id": "synthetic", "in_primary_analysis": False}],
            "source_runs": [
                {
                    "genome_id": "synthetic",
                    "source_run_id": "unmeasured",
                    "source": "phanotate",
                    "status": "success",
                    "cpu_seconds": None,
                },
                {
                    "genome_id": "synthetic",
                    "source_run_id": "measured",
                    "source": "pyrodigal_gv",
                    "status": "success",
                    "cpu_seconds": 2.0,
                },
            ],
        },
    )
    output = tmp_path / "synthetic-run"
    summary = run_pipeline(
        manifest,
        {"max_workers": 1, "storage_batch_genomes": 1},
        output,
        diagnostic=True,
    )
    assert summary["imported_n"] == 1
    assert summary["peak_rss_process_tree_bytes"] is None
    assert summary["peak_rss_unavailable_reason"]
    assert summary["fresh_caller_cpu_seconds"] is None
    assert (
        summary["fresh_caller_cpu_unavailable_reason"] == "fresh_caller_cpu_measurements_incomplete"
    )
    assert summary["total_cpu_observed_seconds"] is None
    assert summary["total_cpu_unavailable_reason"] == "fresh_caller_cpu_measurements_incomplete"
    assert (output / "run_summary.json").is_file()
    assert list((output / "parquet" / "genomes").glob("*.parquet"))
    with duckdb.connect(str(output / "phageannoscope.duckdb"), read_only=True) as db:
        assert db.execute("SELECT count(*) FROM genomes").fetchone()[0] == 1
