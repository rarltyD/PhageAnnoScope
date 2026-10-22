"""Small deterministic provenance helpers."""

import hashlib
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def digest(value):
    if not isinstance(value, bytes):
        value = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(value).hexdigest()


def file_sha(path):
    sha = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8"
    )
    temp.replace(path)


def versions():
    found = {}
    for name in (
        "phageannoscope",
        "biopython",
        "pyarrow",
        "duckdb",
        "pyrodigal-gv",
        "pyrodigal",
        "phanotate",
        "psutil",
    ):
        try:
            found[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            found[name] = None
    return found


def hardware():
    import psutil

    from .resources import probe

    result = {
        "os": platform.platform(),
        "python": sys.version,
        "executable": sys.executable,
        "machine": platform.machine(),
    }
    for name, producer in {
        "cpu_logical": psutil.cpu_count,
        "cpu_physical": lambda: psutil.cpu_count(logical=False),
        "ram_total_bytes": lambda: psutil.virtual_memory().total,
    }.items():
        result[name], result[f"{name}_unavailable_reason"] = probe(producer, name)
    return result
