"""Bridge orchestration unit tests; mocks never stand in for research outputs."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from phageannoscope.callers.phanotate import command_version

BRIDGE_FILE = Path(__file__).resolve().parents[1] / "scripts" / "phanotate_wsl.py"
SPEC = importlib.util.spec_from_file_location("phanotate_wsl_test_module", BRIDGE_FILE)
bridge = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bridge)


def arguments(path=None, **kwargs):
    return SimpleNamespace(
        fasta=str(path) if path is not None else None,
        version=False,
        environment=False,
        bridge_timeout=0.5,
        format="gff",
        **kwargs,
    )


def test_wslpath_uses_argv_and_preserves_spaces_and_unicode(monkeypatch, tmp_path):
    path = tmp_path / "données with spaces.fna"
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        return SimpleNamespace(returncode=0, stdout="/mnt/c/données with spaces.fna\n", stderr="")

    monkeypatch.setattr(bridge.subprocess, "run", run)
    assert bridge._wsl_path(path, "Synthetic Test Distribution") == "/mnt/c/données with spaces.fna"
    argv, kwargs = calls[0]
    assert argv == [
        "wsl.exe",
        "--distribution",
        "Synthetic Test Distribution",
        "--exec",
        "wslpath",
        "-a",
        "-u",
        str(path.resolve()),
    ]
    assert "shell" not in kwargs
    assert kwargs["timeout"] == 30


@pytest.mark.parametrize(
    "returncode,stdout", [(1, "/mnt/c/example"), (0, "relative"), (0, "/mnt/c/a\n/mnt/c/b")]
)
def test_bad_path_conversion_is_rejected(monkeypatch, tmp_path, returncode, stdout):
    monkeypatch.setattr(
        bridge.subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(
            returncode=returncode, stdout=stdout, stderr="synthetic diagnostic"
        ),
    )
    with pytest.raises(RuntimeError, match="path conversion failed"):
        bridge._wsl_path(tmp_path, "Synthetic")


def test_wsl_startup_error_on_stdout_is_decoded_from_utf16(monkeypatch, tmp_path):
    diagnostic = "Accès refusé. Wsl/E_ACCESSDENIED"
    monkeypatch.setattr(
        bridge.subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(
            returncode=1, stdout=diagnostic.encode("utf-16-le"), stderr=b""
        ),
    )
    with pytest.raises(RuntimeError, match="Wsl/E_ACCESSDENIED"):
        bridge._wsl_path(tmp_path, "Synthetic")


def test_windows_command_uses_real_project_python_native_script_and_deadline(monkeypatch, tmp_path):
    monkeypatch.setenv("PHAGEANNOSCOPE_WSL_DISTRIBUTION", "Synthetic Test Distribution")
    monkeypatch.setattr(bridge, "_wsl_path", lambda path, distro: "/mnt/c/" + path.name)
    calls = []
    monkeypatch.setattr(
        bridge.subprocess,
        "run",
        lambda argv, **kwargs: calls.append((argv, kwargs)) or SimpleNamespace(returncode=124),
    )
    args = arguments(tmp_path / "input with spaces.fna")
    assert bridge._windows(args) == 124
    argv, kwargs = calls[0]
    assert argv[:4] == ["wsl.exe", "--distribution", "Synthetic Test Distribution", "--exec"]
    assert argv[4].endswith("/.venv-phanotate-wsl/bin/python")
    assert argv[5].endswith("/phanotate_wsl.py")
    assert argv[6:] == [
        "--native",
        "/mnt/c/input with spaces.fna",
        "-f",
        "gff",
        "--bridge-timeout",
        "0.5",
    ]
    assert "shell" not in kwargs and "timeout" not in kwargs


def test_windows_version_reaches_native_interpreter(monkeypatch):
    monkeypatch.setattr(bridge, "_wsl_path", lambda path, distro: "/mnt/c/" + path.name)
    calls = []
    monkeypatch.setattr(
        bridge.subprocess,
        "run",
        lambda argv, **kwargs: calls.append(argv) or SimpleNamespace(returncode=0),
    )
    args = arguments()
    args.version = True
    assert bridge._windows(args) == 0
    assert calls[0][-2:] == ["--native", "--version"]


def test_windows_environment_reaches_native_interpreter_without_fasta(monkeypatch):
    monkeypatch.setattr(bridge, "_wsl_path", lambda path, distro: "/mnt/c/" + path.name)
    calls = []
    monkeypatch.setattr(
        bridge.subprocess,
        "run",
        lambda argv, **kwargs: calls.append(argv) or SimpleNamespace(returncode=0),
    )
    args = arguments()
    args.environment = True
    assert bridge._windows(args) == 0
    assert calls[0][-2:] == ["--native", "--environment"]


def test_native_environment_reports_observed_packages_and_missing_values(monkeypatch, capsys):
    observed = {"phanotate": "1.6.7", "fastpath": "1.0.3"}
    queried = []

    def version(name):
        queried.append(name)
        if name not in observed:
            raise bridge.importlib.metadata.PackageNotFoundError(name)
        return observed[name]

    monkeypatch.setattr(bridge.importlib.metadata, "version", version)
    monkeypatch.setattr(bridge.platform, "platform", lambda: "Synthetic observed platform")
    monkeypatch.setattr(bridge.platform, "python_version", lambda: "3.12.3")
    monkeypatch.setattr(bridge.platform, "python_implementation", lambda: "CPython")
    args = arguments()
    args.environment = True
    assert bridge._native(args) == 0
    result = json.loads(capsys.readouterr().out)
    assert queried == ["phanotate", "fastpath", "genbank", "setuptools"]
    assert result == {
        "platform": "Synthetic observed platform",
        "python": "3.12.3",
        "python_implementation": "CPython",
        "packages": {
            "phanotate": "1.6.7",
            "fastpath": "1.0.3",
            "genbank": None,
            "setuptools": None,
        },
    }


@pytest.mark.parametrize(
    "child_code,timeout,status,expected_code",
    [(0, False, "success", 0), (7, False, "failed", 7), (0, True, "timeout", 124)],
)
def test_native_deadline_resource_measurements_and_sidecar(
    monkeypatch, tmp_path, child_code, timeout, status, expected_code
):
    # Only subprocess orchestration is mocked; this is not a predictor fixture.
    usages = iter(
        [
            SimpleNamespace(ru_utime=1.0, ru_stime=0.5, ru_maxrss=0),
            SimpleNamespace(ru_utime=2.25, ru_stime=0.75, ru_maxrss=4321),
        ]
    )
    monkeypatch.setitem(
        sys.modules,
        "resource",
        SimpleNamespace(RUSAGE_CHILDREN=-1, getrusage=lambda who: next(usages)),
    )
    times = iter([10.0, 12.5])
    monkeypatch.setattr(bridge.time, "perf_counter", lambda: next(times))
    waited = []

    class Child:
        pid = 12345

        def wait(self, timeout=None):
            waited.append(timeout)
            if timeout is not None and timeout:
                if len(waited) == 1 and should_timeout:
                    raise subprocess.TimeoutExpired("synthetic orchestration", timeout)
            return child_code

    should_timeout = timeout
    popped, killed = [], []
    monkeypatch.setattr(
        bridge.subprocess, "Popen", lambda argv, **kwargs: popped.append((argv, kwargs)) or Child()
    )
    monkeypatch.setattr(
        bridge.os, "killpg", lambda pid, sig: killed.append((pid, sig)), raising=False
    )
    monkeypatch.setattr(bridge.signal, "SIGKILL", 9, raising=False)
    path = tmp_path / "input.fna"
    assert bridge._native(arguments(path)) == expected_code
    resources = json.loads(path.with_suffix(".resources.json").read_text())
    assert resources["status"] == status
    assert resources["exit_code"] == expected_code
    assert resources["prediction_wall_seconds"] == 2.5
    assert resources["cpu_seconds"] == 1.5
    assert resources["peak_rss_bytes"] == 4321 * 1024
    assert "RUSAGE_CHILDREN" in resources["resource_method"]
    assert popped[0][0][1].endswith(".venv-phanotate-wsl\\bin\\phanotate.py") or popped[0][0][
        1
    ].endswith(".venv-phanotate-wsl/bin/phanotate.py")
    assert popped[0][1] == {"start_new_session": True}
    assert killed == ([(12345, 9)] if timeout else [])
    assert waited == ([0.5, None] if timeout else [0.5])


def test_version_probe_preference_failure_and_clean_output(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "phageannoscope.callers.phanotate.subprocess.run",
        lambda argv, **kwargs: (
            calls.append((argv, kwargs)) or SimpleNamespace(returncode=0, stdout="1.6.7\n")
        ),
    )
    assert (
        command_version(
            {"phanotate_observed_version": "1.6.8", "phanotate_version_command": ["ignored"]}
        )
        == "1.6.8"
    )
    assert calls == []
    assert (
        command_version({"phanotate_version_command": ["synthetic-probe", "--version"]}) == "1.6.7"
    )
    assert calls[0][1]["timeout"] == 30
    assert calls[0][1]["capture_output"] is True
    monkeypatch.setattr(
        "phageannoscope.callers.phanotate.subprocess.run",
        lambda *a, **k: SimpleNamespace(returncode=0, stdout="banner\n1.6.7\n"),
    )
    assert command_version({"phanotate_version_command": ["synthetic-probe"]}) is None
    monkeypatch.setattr(
        "phageannoscope.callers.phanotate.subprocess.run",
        lambda *a, **k: SimpleNamespace(returncode=2, stdout="1.6.7"),
    )
    assert command_version({"phanotate_version_command": ["synthetic-probe"]}) is None


def test_version_probe_timeout_returns_unavailable(monkeypatch):
    def fail(*args, **kwargs):
        raise subprocess.TimeoutExpired("synthetic-probe", 30)

    monkeypatch.setattr("phageannoscope.callers.phanotate.subprocess.run", fail)
    assert command_version({"phanotate_version_command": ["synthetic-probe"]}) is None


def test_no_assumed_version_when_distribution_is_absent(monkeypatch):
    def absent(*args):
        raise bridge.importlib.metadata.PackageNotFoundError("phanotate")

    monkeypatch.setattr("phageannoscope.callers.phanotate.importlib.metadata.version", absent)
    assert command_version() is None
