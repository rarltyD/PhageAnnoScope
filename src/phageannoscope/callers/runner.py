"""Timeouts, immutable raw logs, cache provenance and measured resources."""

import json
import subprocess
import sys
import time
import uuid
from pathlib import Path

from ..resources import ProcessTreeSampler, probe
from ..util import digest, file_sha, utcnow, versions, write_json
from .phanotate import command, command_version, parse_gff


def adapter_sha(source):
    root = Path(__file__).resolve().parents[1]
    names = [
        "genbank_io.py",
        "locations.py",
        "translation.py",
        "callers/runner.py",
        "callers/worker.py",
        "resources.py",
        f"callers/{source}.py",
    ]
    hashes = {name: file_sha(root / name) for name in names}
    bridge = root.parents[1] / "scripts" / "phanotate_wsl.py"
    if source == "phanotate" and bridge.exists():
        hashes["scripts/phanotate_wsl.py"] = file_sha(bridge)
    return digest(hashes)


def cache_key(genome, source, version, params, rotation_shift=0, adapter_hash=None):
    return digest(
        {
            "genome_id": genome["genome_id"],
            "input_sha256": genome["sequence_sha256"],
            "source": source,
            "version": version,
            "params": params,
            "rotation_shift": rotation_shift,
            "adapter_sha256": adapter_hash or adapter_sha(source),
        }
    )


def monitor(argv, stdout, stderr, timeout):
    started = time.perf_counter()
    status = "success"
    with Path(stdout).open("wb") as out, Path(stderr).open("wb") as err:
        child = subprocess.Popen(argv, stdout=out, stderr=err, shell=False)
        sampler = ProcessTreeSampler(child.pid)
        while child.poll() is None:
            sampler.sample()
            if time.perf_counter() - started > timeout:
                status = "timeout"
                sampler.kill_descendants()
                try:
                    child.kill()
                except ProcessLookupError:
                    # The child can exit between poll() and kill(). wait() still reaps it.
                    pass
                break
            time.sleep(0.02)
        code = child.wait()
    if status != "timeout" and code:
        status = "failed"
    return {
        "status": status,
        "exit_code": code,
        "wall_seconds": time.perf_counter() - started,
        **sampler.result(),
        "resource_method": "psutil process tree, sampled 20ms; CPU/RSS lower bound",
    }


def merge_native_resources(run, resource_path):
    """Keep launcher/native observations separate; telemetry cannot override success."""
    for metric, launcher_field, reason_field, launcher_reason in (
        (
            "cpu_seconds",
            "launcher_cpu_sampled_seconds",
            "cpu_unavailable_reason",
            "launcher_cpu_unavailable_reason",
        ),
        (
            "peak_rss_bytes",
            "launcher_peak_rss_sampled_bytes",
            "peak_rss_unavailable_reason",
            "launcher_peak_rss_unavailable_reason",
        ),
    ):
        run[launcher_field] = run.get(metric)
        run[launcher_reason] = run.get(reason_field)
        run[metric] = None
        run[reason_field] = "native_resource_sidecar_unavailable"

    def read_sidecar():
        native = json.loads(resource_path.read_text(encoding="utf-8"))
        if not isinstance(native, dict):
            raise ValueError("native resource sidecar must be an object")
        if native.get("exit_code") != run["exit_code"]:
            raise ValueError("native resource sidecar exit code differs from subprocess")
        return native

    native, reason = probe(read_sidecar, "native_resource_sidecar")
    if reason:
        run.update(
            cpu_unavailable_reason=reason,
            peak_rss_unavailable_reason=reason,
            resource_monitor_status="unavailable",
            resource_method="native WSL measurements unavailable; launcher observations separate",
        )
        return
    for field in ("prediction_wall_seconds", "resource_method", "native_command_argv"):
        if field in native:
            run[field] = native[field]
    for metric, reason_field in (
        ("cpu_seconds", "cpu_unavailable_reason"),
        ("peak_rss_bytes", "peak_rss_unavailable_reason"),
    ):
        run[metric] = native.get(metric)
        run[reason_field] = (
            native.get(reason_field) or "native_resource_measurement_unavailable"
            if run[metric] is None
            else None
        )
    observed = sum(run[metric] is not None for metric in ("cpu_seconds", "peak_rss_bytes"))
    run["resource_monitor_status"] = (
        "available" if observed == 2 else "partial" if observed else "unavailable"
    )
    # The bridge's native deadline uses this exit code; launcher monitoring sees a failure.
    if native.get("status") == "timeout" and run["exit_code"] == 124:
        run["status"] = "timeout"


def run_caller(genome, source, config, raw_root, rotation_shift=0, supported_codes=None):
    installed = versions()
    version = (
        command_version(config)
        if source == "phanotate"
        else installed.get(source.replace("_", "-"))
    )
    params = (
        config["pyrodigal_gv"]
        if source == "pyrodigal_gv"
        else {"format": "gff", "translation_table": 11}
    )
    phanotate = command(config.get("phanotate_command")) if source == "phanotate" else None
    key = cache_key(
        genome,
        source,
        version,
        {
            **params,
            "command": phanotate,
            "dependencies": installed,
            "external_environment": config.get("phanotate_observed_environment"),
            "supported_codes": supported_codes,
        },
        rotation_shift,
    )
    raw = Path(raw_root).resolve() / key[:2] / key
    raw.mkdir(parents=True, exist_ok=True)
    run = {
        "source_run_id": f"run_{key}",
        "genome_id": genome["genome_id"],
        "source": source,
        "algorithm_family": source,
        "tool_version": version,
        "params_json": json.dumps(params, sort_keys=True),
        "input_sha256": genome["sequence_sha256"],
        "rotation_shift": rotation_shift,
        "cache_key": key,
        "adapter_sha256": adapter_sha(source),
        "started_at_utc": utcnow(),
        "status": "not_run",
        "exit_code": None,
        "wall_seconds": None,
        "cpu_seconds": None,
        "peak_rss_bytes": None,
        "cpu_unavailable_reason": "subprocess_not_started",
        "peak_rss_unavailable_reason": "subprocess_not_started",
        "resource_monitor_status": "not_started",
        "raw_stdout": str(raw / "stdout.txt"),
        "raw_stderr": str(raw / "stderr.txt"),
    }
    if source == "phanotate":
        run["external_environment_json"] = json.dumps(
            config.get("phanotate_observed_environment"), sort_keys=True
        )
    metadata = raw / "source_run.json"
    calls_path = raw / "calls.json"
    if config.get("use_caller_cache", True) and metadata.exists() and calls_path.exists():
        prior = json.loads(metadata.read_text(encoding="utf-8"))
        if prior["status"] == "success" and prior.get("calls_sha256") == file_sha(calls_path):
            prior.update(
                status="cache_hit", cache_origin_status="success", cache_read_at_utc=utcnow()
            )
            return prior, json.loads(calls_path.read_text(encoding="utf-8"))
    attempt = raw / "attempts" / uuid.uuid4().hex
    attempt.mkdir(parents=True)
    run.update(
        attempt_id=attempt.name,
        raw_stdout=str(attempt / "stdout.txt"),
        raw_stderr=str(attempt / "stderr.txt"),
    )
    reason = None
    if genome["length_bp"] > config["max_genome_length_bp"]:
        reason = "genome_exceeds_configured_resource_limit"
    elif not version or (source == "phanotate" and phanotate is None):
        reason = "caller_not_installed"
    elif (
        source == "phanotate" and supported_codes and any(c not in (1, 11) for c in supported_codes)
    ):
        run["status"] = "unsupported_code"
        reason = "PHANOTATE adapter supports standard stops only"
    if reason:
        run.update(error=reason, ended_at_utc=utcnow())
        write_json(attempt / "source_run.json", run)
        write_json(metadata, run)
        return run, []
    seq = genome["sequence"]
    rotated = seq[rotation_shift:] + seq[:rotation_shift]
    fasta = attempt / "input.fna"
    fasta.write_text(f">{genome['accession_version']}\n{rotated}\n", encoding="ascii")
    if source == "phanotate":
        argv = phanotate + [str(fasta), "-f", "gff"]
        if config.get("phanotate_wsl_bridge"):
            argv.extend(["--bridge-timeout", str(max(0.1, config["caller_timeout_seconds"] - 10))])
    else:
        request = attempt / "request.json"
        write_json(
            request,
            {
                "genome": genome,
                "source_run_id": run["source_run_id"],
                "rotation_shift": rotation_shift,
                "params": params,
            },
        )
        argv = [
            sys.executable,
            "-m",
            "phageannoscope.callers.worker",
            str(request),
            str(attempt / "worker_result.json"),
        ]
    run["command_argv"] = argv
    run["input_fasta_sha256"] = file_sha(fasta)
    try:
        run.update(
            monitor(argv, run["raw_stdout"], run["raw_stderr"], config["caller_timeout_seconds"])
        )
        resource_path = fasta.with_suffix(".resources.json")
        if source == "phanotate" and config.get("phanotate_wsl_bridge"):
            merge_native_resources(run, resource_path)
        calls = []
        if run["status"] == "success":
            if source == "phanotate":
                calls = parse_gff(
                    Path(run["raw_stdout"]).read_text(),
                    genome,
                    run["source_run_id"],
                    rotation_shift,
                )
            else:
                result = json.loads((attempt / "worker_result.json").read_text(encoding="utf-8"))
                calls = result.pop("calls")
                run.update(result)
                run["prediction_cpu_seconds"] = result["worker_cpu_seconds"]
            run["cds_count"] = len(calls)
            write_json(calls_path, calls)
            write_json(attempt / "calls.json", calls)
            run["calls_sha256"] = file_sha(calls_path)
        else:
            run["error"] = Path(run["raw_stderr"]).read_text(encoding="utf-8", errors="replace")[
                -4000:
            ]
    except (OSError, ValueError, KeyError) as error:
        run.update(status="failed", error=f"{type(error).__name__}: {error}")
        calls = []
    run["ended_at_utc"] = utcnow()
    write_json(attempt / "source_run.json", run)
    write_json(metadata, run)
    return run, calls
