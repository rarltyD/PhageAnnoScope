"""Process/cache tests use no simulated biological research results."""

import sys

from phageannoscope.callers.runner import cache_key, monitor, run_caller
from phageannoscope.config import load_config


def test_process_error_and_timeout_are_distinct(tmp_path):
    error = monitor(
        [sys.executable, "-c", "raise SystemExit(7)"], tmp_path / "o", tmp_path / "e", 10
    )
    assert error["status"] == "failed"
    assert error["exit_code"] == 7
    timeout = monitor(
        [sys.executable, "-c", "import time; time.sleep(10)"], tmp_path / "o2", tmp_path / "e2", 0.1
    )
    assert timeout["status"] == "timeout"
    assert timeout["wall_seconds"] < 5


def test_cache_identity_includes_every_prediction_input():
    genome = {"genome_id": "synthetic", "sequence_sha256": "a"}
    base = cache_key(genome, "pyrodigal_gv", "0.3.2", {"meta": True}, adapter_hash="A")
    alternatives = [
        cache_key(
            {**genome, "sequence_sha256": "b"},
            "pyrodigal_gv",
            "0.3.2",
            {"meta": True},
            adapter_hash="A",
        ),
        cache_key(genome, "pyrodigal_gv", "changed", {"meta": True}, adapter_hash="A"),
        cache_key(genome, "pyrodigal_gv", "0.3.2", {"meta": False}, adapter_hash="A"),
        cache_key(genome, "pyrodigal_gv", "0.3.2", {"meta": True}, adapter_hash="B"),
        cache_key(
            genome, "pyrodigal_gv", "0.3.2", {"meta": True}, rotation_shift=10, adapter_hash="A"
        ),
    ]
    assert all(value != base for value in alternatives)


def test_missing_caller_is_not_successful_zero(tmp_path, monkeypatch):
    monkeypatch.setattr("phageannoscope.callers.runner.versions", lambda: {"phanotate": None})
    genome = {
        "genome_id": "synthetic",
        "sequence_sha256": "a",
        "length_bp": 9,
        "sequence": "ATGAAATAA",
        "accession_version": "SYNTHETIC.1",
    }
    result, calls = run_caller(genome, "phanotate", load_config(), tmp_path)
    assert result["status"] == "not_run"
    assert result["wall_seconds"] is None
    assert calls == []


def test_changed_code_policy_cannot_reuse_cache_key(tmp_path, monkeypatch):
    monkeypatch.setattr("phageannoscope.callers.runner.versions", lambda: {"phanotate": None})
    genome = {"genome_id": "synthetic", "sequence_sha256": "a", "length_bp": 9}
    a, _ = run_caller(genome, "phanotate", load_config(), tmp_path, supported_codes=[11])
    b, _ = run_caller(genome, "phanotate", load_config(), tmp_path, supported_codes=[15])
    assert a["cache_key"] != b["cache_key"]
