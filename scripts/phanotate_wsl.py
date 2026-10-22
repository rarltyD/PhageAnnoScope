"""Run the project's real Linux PHANOTATE installation from Windows.

The Linux process owns its deadline and child process group. Resource data
describe that native predictor, independently of the Windows WSL launcher.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
import platform
import signal
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_DISTRIBUTION = "MetaPool-Ubuntu-24.04"
DEFAULT_TIMEOUT = 290.0


def _decode_output(value) -> str:
    # WSL startup failures are UTF-16 on Windows; native Linux stdout is UTF-8.
    if isinstance(value, str):
        return value
    if value.startswith((b"\xff\xfe", b"\xfe\xff")):
        return value.decode("utf-16", errors="replace")
    return value.decode("utf-16-le" if b"\x00" in value else "utf-8", errors="replace")


def _wsl_path(path: Path, distribution: str) -> str:
    result = subprocess.run(
        [
            "wsl.exe",
            "--distribution",
            distribution,
            "--exec",
            "wslpath",
            "-a",
            "-u",
            str(path.resolve()),
        ],
        capture_output=True,
        timeout=30,
        check=False,
    )
    converted = _decode_output(result.stdout).strip()
    if result.returncode or not converted.startswith("/") or "\n" in converted or "\r" in converted:
        diagnostic = _decode_output(result.stderr).strip() or converted
        raise RuntimeError(f"WSL path conversion failed for {path}: {diagnostic}")
    return converted


def _windows(args) -> int:
    distribution = os.environ.get("PHAGEANNOSCOPE_WSL_DISTRIBUTION", DEFAULT_DISTRIBUTION)
    root = Path(__file__).resolve().parents[1]
    linux_root = _wsl_path(root, distribution)
    linux_script = _wsl_path(Path(__file__), distribution)
    argv = [
        "wsl.exe",
        "--distribution",
        distribution,
        "--exec",
        linux_root + "/.venv-phanotate-wsl/bin/python",
        linux_script,
        "--native",
    ]
    if args.version:
        argv.append("--version")
    elif args.environment:
        argv.append("--environment")
    else:
        argv.extend(
            [
                _wsl_path(Path(args.fasta), distribution),
                "-f",
                args.format,
                "--bridge-timeout",
                str(args.bridge_timeout),
            ]
        )
    # The native process enforces the deadline. Killing this launcher does not
    # reliably terminate a Linux child, so the Windows wrapper owns no timer.
    return subprocess.run(argv, check=False).returncode


def _native(args) -> int:
    if args.version:
        print(importlib.metadata.version("phanotate"))
        return 0
    if args.environment:
        packages = {}
        for distribution in ("phanotate", "fastpath", "genbank", "setuptools"):
            try:
                packages[distribution] = importlib.metadata.version(distribution)
            except importlib.metadata.PackageNotFoundError:
                packages[distribution] = None
        print(
            json.dumps(
                {
                    "platform": platform.platform(),
                    "python": platform.python_version(),
                    "python_implementation": platform.python_implementation(),
                    "packages": packages,
                },
                sort_keys=True,
            )
        )
        return 0

    # Native telemetry is independent of the predictor's exit status.
    def usage():
        try:
            import resource

            value = resource.getrusage(resource.RUSAGE_CHILDREN)
            cpu = float(value.ru_utime + value.ru_stime)
            rss = int(value.ru_maxrss * 1024)
            if not math.isfinite(cpu) or cpu < 0 or rss < 0:
                raise ValueError("invalid native resource measurement")
            return (cpu, rss), None
        except Exception as error:
            return None, f"native_getrusage: {type(error).__name__}: {error}"

    fasta = Path(args.fasta).resolve()
    root = Path(__file__).resolve().parents[1]
    argv = [
        sys.executable,
        str(root / ".venv-phanotate-wsl" / "bin" / "phanotate.py"),
        str(fasta),
        "-f",
        args.format,
    ]
    before, before_reason = usage()
    started = time.perf_counter()
    code, status = 127, "failed"
    try:
        child = subprocess.Popen(argv, start_new_session=True)
        try:
            code = child.wait(timeout=args.bridge_timeout)
            status = "success" if code == 0 else "failed"
        except subprocess.TimeoutExpired:
            # The predictor and any descendants share this native session.
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
            code, status = 124, "timeout"
    except OSError as error:
        print(f"PHANOTATE native launch failed: {error}", file=sys.stderr)
    elapsed = time.perf_counter() - started
    after, after_reason = usage()
    cpu = after[0] - before[0] if before is not None and after is not None else None
    cpu_reason = (after_reason or before_reason) if cpu is None else None
    if cpu is not None and cpu < 0:
        cpu, cpu_reason = None, "native_getrusage: CPU counter decreased"
    metrics = {
        "status": status,
        "exit_code": code,
        "prediction_wall_seconds": elapsed,
        "cpu_seconds": cpu,
        "cpu_unavailable_reason": cpu_reason,
        "peak_rss_bytes": after[1] if after is not None else None,
        "peak_rss_unavailable_reason": after_reason,
        "resource_method": "Linux resource.getrusage(RUSAGE_CHILDREN); ru_maxrss KiB converted to bytes",
        "native_command_argv": argv,
    }
    destination = fasta.with_suffix(".resources.json")
    temporary = destination.with_name(destination.name + f".{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
        temporary.replace(destination)
    except OSError as error:
        print(f"PHANOTATE resource sidecar unavailable: {error}", file=sys.stderr)
    return code


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", nargs="?")
    parser.add_argument("-f", "--format", default="gff")
    parser.add_argument("--bridge-timeout", type=float, default=DEFAULT_TIMEOUT)
    inspection = parser.add_mutually_exclusive_group()
    inspection.add_argument("--version", action="store_true")
    inspection.add_argument("--environment", action="store_true")
    parser.add_argument("--native", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if not args.version and not args.environment and not args.fasta:
        parser.error("a FASTA path is required")
    if not 0 < args.bridge_timeout < float("inf"):
        parser.error("--bridge-timeout must be finite and positive")
    try:
        if args.native:
            if not sys.platform.startswith("linux"):
                parser.error("--native requires Linux")
            return _native(args)
        return _windows(args)
    except (
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        importlib.metadata.PackageNotFoundError,
    ) as error:
        print(f"PHANOTATE WSL bridge failed: {error}", file=sys.stderr)
        return 127


if __name__ == "__main__":
    raise SystemExit(main())
