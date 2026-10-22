"""Test native orchestration telemetry failures, without simulated gene calls."""

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "native_bridge", Path(__file__).resolve().parents[1] / "scripts/phanotate_wsl.py"
)
bridge = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bridge)


def arguments(path):
    return SimpleNamespace(
        version=False, environment=False, fasta=str(path), format="gff", bridge_timeout=1
    )


@pytest.mark.parametrize("exit_code", [0, 7])
@pytest.mark.parametrize(
    "error",
    [
        PermissionError("resource denied"),
        NotImplementedError("no rusage"),
        OSError("resource unavailable"),
    ],
)
def test_native_resource_failure_does_not_change_predictor_exit(
    monkeypatch, tmp_path, exit_code, error
):
    def unavailable(who):
        raise error

    monkeypatch.setitem(
        sys.modules, "resource", SimpleNamespace(RUSAGE_CHILDREN=-1, getrusage=unavailable)
    )
    monkeypatch.setattr(
        bridge.subprocess, "Popen", lambda *a, **k: SimpleNamespace(wait=lambda **kwargs: exit_code)
    )
    path = tmp_path / "input.fna"
    assert bridge._native(arguments(path)) == exit_code
    measured = json.loads(path.with_suffix(".resources.json").read_text())
    assert measured["status"] == ("success" if exit_code == 0 else "failed")
    assert measured["cpu_seconds"] is measured["peak_rss_bytes"] is None
    assert type(error).__name__ in measured["cpu_unavailable_reason"]
    assert type(error).__name__ in measured["peak_rss_unavailable_reason"]


def test_native_sidecar_write_failure_preserves_success(monkeypatch, tmp_path, capsys):
    monkeypatch.setitem(
        sys.modules,
        "resource",
        SimpleNamespace(
            RUSAGE_CHILDREN=-1,
            getrusage=lambda who: SimpleNamespace(ru_utime=0, ru_stime=0, ru_maxrss=0),
        ),
    )
    monkeypatch.setattr(
        bridge.subprocess, "Popen", lambda *a, **k: SimpleNamespace(wait=lambda **kwargs: 0)
    )

    def denied(*args, **kwargs):
        raise PermissionError("sidecar denied")

    monkeypatch.setattr(Path, "write_text", denied)
    assert bridge._native(arguments(tmp_path / "input.fna")) == 0
    assert "resource sidecar unavailable" in capsys.readouterr().err
