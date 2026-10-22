"""Run settings are separate from the frozen scientific protocol."""

import json
import subprocess
from pathlib import Path

DEFAULTS = {
    "caller_timeout_seconds": 300,
    "max_genome_length_bp": 2_000_000,
    "max_workers": 4,
    "partial_overlap_reciprocal_primary": 0.8,
    "pyrodigal_gv": {"meta": True},
    "phanotate_command": None,
    "storage_batch_genomes": 20,
}


def load_config(path=None):
    config = dict(DEFAULTS)
    if path:
        config.update(json.loads(Path(path).read_text(encoding="utf-8")))
    if config["caller_timeout_seconds"] <= 0 or config["max_workers"] < 1:
        raise ValueError("Timeout and max_workers must be positive")
    if not 0 < config["partial_overlap_reciprocal_primary"] <= 1:
        raise ValueError("Reciprocal coverage threshold must be in (0,1]")
    # Native parameters are mandatory for the protocol's primary comparison.
    if config.get("native_caller_settings", True) and config["pyrodigal_gv"] != {"meta": True}:
        raise ValueError("Primary protocol requires native ViralGeneFinder(meta=True) settings")
    return config


def resolve_external_tools(config):
    """Probe configured external environments once per run, before cache identity."""
    from .callers.phanotate import command_version

    config = dict(config)
    config.pop("phanotate_observed_version", None)
    config.pop("phanotate_observed_environment", None)
    if config.get("phanotate_version_command"):
        config["phanotate_observed_version"] = command_version(config)
    if config.get("phanotate_environment_command"):
        try:
            result = subprocess.run(
                config["phanotate_environment_command"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                timeout=45,
                check=True,
            )
            config["phanotate_observed_environment"] = json.loads(result.stdout)
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            config["phanotate_environment_probe_error"] = type(error).__name__
    return config
