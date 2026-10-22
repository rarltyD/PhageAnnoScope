"""Predeclared same-stratum reserves with complete failed-attempt provenance."""

import json
from collections import Counter
from pathlib import Path

import duckdb

from .manifest import BIN_NAMES, _bin, _block, _events, _rank, _read, _true, _write
from .util import file_sha, utcnow, write_json


def replace_selection(manifest, previous, run, out, seed=30102026, attempted=None):
    """Keep eligible successes and replace failures in their original length strata.

    This changes a resource pilot's attempted set, not the study population.
    All prior failures remain in the journal and in the previous run directories.
    """
    run, out = Path(run), Path(out)
    if out.resolve() in {Path(manifest).resolve(), Path(previous).resolve()}:
        raise ValueError("Replacement must create a new manifest")
    config = json.loads((run / "run_config.json").read_text(encoding="utf-8"))
    if file_sha(previous) != config["manifest_sha256"]:
        raise ValueError("Previous manifest does not match the evaluated run")
    prior = _read(previous)
    if not prior or len({r["accession_version"] for r in prior}) != len(prior):
        raise ValueError("Previous manifest must contain unique accession.version rows")
    if any(r.get("selection_seed") and int(r["selection_seed"]) != seed for r in prior):
        raise ValueError("Replacement seed differs from the original pilot")
    eligible = [r for r in _read(manifest) if _true(r.get("eligible_for_pilot"))]
    by_id = {r["accession_version"]: r for r in eligible}
    for row in prior:
        identity = row["accession_version"]
        if identity not in by_id or row["sequence_sha256"] != by_id[identity]["sequence_sha256"]:
            raise ValueError("Previous accession.version/SHA is no longer eligible in this cohort")
    old_summary = Path(previous).with_suffix(".summary.json")
    history = set()
    if old_summary.exists():
        history.update(
            json.loads(old_summary.read_text(encoding="utf-8")).get("attempted_accessions", [])
        )
    for path in attempted or []:
        history.update(r["accession_version"] for r in _read(path))
    history.update(r["accession_version"] for r in prior)
    with duckdb.connect(str(run / "phageannoscope.duckdb"), read_only=True) as db:
        accepted = {
            r[0]
            for r in db.execute(
                "SELECT accession_version FROM genomes WHERE in_primary_analysis"
            ).fetchall()
        }
        failure_reasons = dict(
            db.execute("""SELECT g.accession_version,
            string_agg(DISTINCT r.source || ':' || r.status || coalesce(':' || r.comparison_reason,''),';' ORDER BY r.source || ':' || r.status || coalesce(':' || r.comparison_reason,''))
            FROM genomes g LEFT JOIN source_runs r USING(genome_id)
            WHERE NOT coalesce(g.in_primary_analysis,false) GROUP BY g.accession_version""").fetchall()
        )
    if not accepted <= {r["accession_version"] for r in prior}:
        raise ValueError("Run contains genomes outside the previous selection")
    chosen = [
        dict(by_id[r["accession_version"]]) for r in prior if r["accession_version"] in accepted
    ]
    failures = [r for r in prior if r["accession_version"] not in accepted]
    pool = sorted(
        (r for r in eligible if r["accession_version"] not in history),
        key=lambda r: (_rank(seed, r), r["accession_version"]),
    )
    timestamp = utcnow()
    events, deficit = [], Counter()
    for failed in sorted(failures, key=lambda r: (_bin(r), _rank(seed, r), r["accession_version"])):
        stratum = _bin(failed)
        blocks = Counter(_block(r) for r in chosen if _bin(r) == stratum)
        candidates = [r for r in pool if _bin(r) == stratum]
        replacement = (
            min(
                candidates,
                key=lambda r: (blocks[_block(r)], _rank(seed, r), r["accession_version"]),
            )
            if candidates
            else None
        )
        if replacement is None:
            deficit[BIN_NAMES[stratum]] += 1
        else:
            pool.remove(replacement)
            chosen.append(dict(replacement))
            history.add(replacement["accession_version"])
        events.append(
            {
                "run_id": "replace_" + run.name,
                "accession": failed.get("accession"),
                "accession_version": failed["accession_version"],
                "stage": "pilot_replacement",
                "decision": "replaced" if replacement else "reserve_exhausted",
                "primary_reason_code": "CALLER_QC_NOT_ELIGIBLE",
                "all_reason_codes_json": json.dumps(["CALLER_QC_NOT_ELIGIBLE"]),
                "evidence_json": json.dumps(
                    {
                        "failed_run": str(run.resolve()),
                        "failure": failure_reasons.get(
                            failed["accession_version"], "import_failed_or_not_evaluated"
                        ),
                        "stratum": BIN_NAMES[stratum],
                        "replacement": replacement["accession_version"] if replacement else None,
                        "seed": seed,
                    }
                ),
                "evaluated_at_utc": timestamp,
            }
        )
    chosen.sort(key=lambda r: (_bin(r), _rank(seed, r), r["accession_version"]))
    for index, row in enumerate(chosen, 1):
        row.update(
            selection_seed=seed,
            selection_strategy="protocol-length-strata",
            selection_stratum=BIN_NAMES[_bin(row)],
            selection_hash=_rank(seed, row),
            selection_order=index,
            replacement_retained=str(row["accession_version"] in accepted).lower(),
        )
    _write(out, chosen)
    event_path = out.parent / "replacement_events.csv"
    _events(event_path, events)
    summary = {
        "status": "incomplete" if deficit else "selected",
        "requested_n": len(prior),
        "selected_n": len(chosen),
        "retained_successful_n": len(accepted),
        "previous_failed_n": len(failures),
        "deficit_by_stratum": dict(deficit),
        "strata": dict(Counter(BIN_NAMES[_bin(r)] for r in chosen)),
        "seed": seed,
        "attempted_accessions": sorted(history),
        "unique_attempted_n": len(history),
        "previous_manifest_sha256": file_sha(previous),
        "cohort_sha256": file_sha(manifest),
        "out": str(out.resolve()),
        "events": str(event_path.resolve()),
    }
    write_json(out.with_suffix(".summary.json"), summary)
    return summary
