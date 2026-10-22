"""Synthetic resource-pilot replacements; no scientific results are generated."""

import csv
import hashlib
import json
from collections import Counter

import duckdb
import pytest

from phageannoscope.manifest import BIN_NAMES, select_manifest
from phageannoscope.replacement import replace_selection
from phageannoscope.util import file_sha


def write(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)


def read(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def cohort(tmp_path, per_bin=8):
    rows = []
    for i, length in enumerate((9000, 10000, 50000, 100000, 200000)):
        for j in range(per_bin):
            version = f"SYNTH{i}{j}.1"
            rows.append(
                {
                    "accession_version": version,
                    "sequence_sha256": hashlib.sha256(version.encode()).hexdigest(),
                    "ncbi_length_bp": length,
                    "eligible_for_pilot": "true",
                }
            )
    path = tmp_path / "cohort.csv"
    write(path, rows)
    selected = tmp_path / "pilot.csv"
    select_manifest(path, 25, out=selected)
    return path, selected


def evaluated(tmp_path, selected, failed, name="run"):
    out = tmp_path / name
    out.mkdir()
    (out / "run_config.json").write_text(json.dumps({"manifest_sha256": file_sha(selected)}))
    with duckdb.connect(str(out / "phageannoscope.duckdb")) as db:
        db.execute(
            "CREATE TABLE genomes(genome_id VARCHAR, accession_version VARCHAR, in_primary_analysis BOOLEAN)"
        )
        db.execute(
            "CREATE TABLE source_runs(genome_id VARCHAR,source VARCHAR,status VARCHAR,comparison_reason VARCHAR)"
        )
        for row in read(selected):
            identity = row["accession_version"]
            db.execute(
                "INSERT INTO genomes VALUES (?,?,?)", [identity, identity, identity not in failed]
            )
            db.execute(
                "INSERT INTO source_runs VALUES (?,?,?,?)",
                [identity, "phanotate", "timeout" if identity in failed else "success", None],
            )
    return out


def test_replacement_keeps_successes_and_original_stratum_counts(tmp_path):
    source, selected = cohort(tmp_path)
    old = read(selected)
    failed = {old[0]["accession_version"], old[8]["accession_version"]}
    run = evaluated(tmp_path, selected, failed)
    out = tmp_path / "replacement.csv"
    result = replace_selection(source, selected, run, out)
    rows = read(out)
    assert result["selected_n"] == 25 and result["retained_successful_n"] == 23
    assert not failed & {r["accession_version"] for r in rows}
    assert result["strata"] == dict.fromkeys(BIN_NAMES, 5)
    assert result["unique_attempted_n"] == 27
    events = read(tmp_path / "replacement_events.csv")
    assert len(events) == 2 and all(e["decision"] == "replaced" for e in events)
    second = tmp_path / "same.csv"
    replace_selection(source, selected, run, second)
    assert out.read_bytes() == second.read_bytes()


def test_second_round_does_not_reuse_any_prior_failure(tmp_path):
    source, selected = cohort(tmp_path, per_bin=10)
    failed = {read(selected)[0]["accession_version"]}
    first_run = evaluated(tmp_path, selected, failed)
    first = tmp_path / "round1.csv"
    replace_selection(source, selected, first_run, first)
    new_ids = {r["accession_version"] for r in read(first)} - {
        r["accession_version"] for r in read(selected)
    }
    second_run = evaluated(tmp_path, first, new_ids, "run2")
    second = tmp_path / "round2.csv"
    result = replace_selection(source, first, second_run, second)
    assert not (failed | new_ids) & {r["accession_version"] for r in read(second)}
    assert result["unique_attempted_n"] == 27


def test_exhausted_reserve_reports_deficit_without_false_n(tmp_path):
    source, selected = cohort(tmp_path, per_bin=5)
    failed = {read(selected)[0]["accession_version"]}
    result = replace_selection(
        source, selected, evaluated(tmp_path, selected, failed), tmp_path / "deficit.csv"
    )
    assert result["status"] == "incomplete" and result["selected_n"] == 24
    assert result["deficit_by_stratum"] == {"lt10k": 1}
    assert (
        sum(Counter(r["selection_stratum"] for r in read(tmp_path / "deficit.csv")).values()) == 24
    )


def test_changed_manifest_or_seed_cannot_reinterpret_run(tmp_path):
    source, selected = cohort(tmp_path)
    run = evaluated(tmp_path, selected, set())
    with pytest.raises(ValueError, match="seed"):
        replace_selection(source, selected, run, tmp_path / "wrong.csv", seed=999)
    selected.write_text(selected.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="match"):
        replace_selection(source, selected, run, tmp_path / "wrong.csv")
