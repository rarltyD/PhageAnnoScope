"""Coordinate agreement and stop loci, without treating a baseline as truth.

All public functions accept and return ordinary row dictionaries.  Coordinates
are 0-based half-open and segment order is translation order.  Failed callers
are never converted into successful callers that predicted zero genes.
"""

from __future__ import annotations

import hashlib
import json
from bisect import bisect_left
from collections import defaultdict
from itertools import combinations
from statistics import median
from typing import Any, Iterable

from .locations import coalesce_segments, stop_positions

SUCCESS = frozenset({"success", "cache_hit"})
CALLER_FAMILIES = ("phanotate", "pyrodigal_gv")


def _stable(prefix: str, value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return prefix + hashlib.sha256(raw.encode()).hexdigest()


def _family(value: Any) -> str:
    value = str(value or "unknown").lower().replace("-", "_")
    # Prokka/Prodigal/Pyrodigal do not create additional independent votes.
    return "pyrodigal_gv" if value in {"pyrodigal_gv", "pyrodigal", "prodigal", "prokka"} else value


def _is_r0(row: dict) -> bool:
    shift = row.get("rotation_shift", row.get("rotation", 0))
    return shift in (None, 0, "0", "r0", "original")


def _sha(row: dict) -> str:
    return str(
        row.get("sequence_sha256")
        or row.get("input_sha256")
        or row.get("input_sequence_sha256")
        or ""
    )


def _code(row: dict) -> str | None:
    if row.get("code_class"):
        return str(row["code_class"])
    table = row.get(
        "transl_table_effective",
        row.get("translation_table", row.get("selected_translation_table")),
    )
    if table is None:
        return None
    return "standard_stop" if int(table) in (1, 11) else "table_" + str(table)


def _geometry(row: dict) -> dict:
    segments = [
        {"start0": int(s["start0"]), "end0": int(s["end0"])} for s in row.get("segments", [])
    ]
    if all(0 <= s["start0"] < s["end0"] for s in segments):
        segments = coalesce_segments(segments, row.get("strand"))
    return {
        "sequence_sha256": _sha(row),
        "strand": row.get("strand"),
        "segments": segments,
        "codon_start": row.get("codon_start", 1),
        "code_class": _code(row),
    }


def _model_id(row: dict) -> str:
    return str(row.get("model_id") or _stable("model_", _geometry(row)))


def _union(row: dict) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted((int(s["start0"]), int(s["end0"])) for s in row.get("segments", [])):
        if start < 0 or end <= start:
            continue
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def _length(row: dict) -> int:
    return sum(end - start for start, end in _union(row))


def _overlap(a: dict, b: dict) -> int:
    ua, ub = _union(a), _union(b)
    ia = ib = total = 0
    while ia < len(ua) and ib < len(ub):
        total += max(0, min(ua[ia][1], ub[ib][1]) - max(ua[ia][0], ub[ib][0]))
        if ua[ia][1] <= ub[ib][1]:
            ia += 1
        else:
            ib += 1
    return total


def _invalid_reason(row: dict) -> str | None:
    if row.get("strand") not in (1, -1):
        return "unknown_strand"
    if row.get("codon_start", 1) not in (1, 2, 3):
        return "invalid_codon_start"
    if not row.get("segments") or any(
        int(s["start0"]) < 0 or int(s["end0"]) <= int(s["start0"]) for s in row["segments"]
    ):
        return "invalid_coordinates"
    if row.get("pseudo") or row.get("exception") or row.get("transl_except"):
        return "complex_biology"
    if row.get("location_operator") == "order":
        return "unordered_location"
    if row.get("location_issues"):
        return str(row["location_issues"][0])
    validity = row.get("translation_validity")
    if validity in {
        "invalid_codon_start",
        "unsupported_table",
        "invalid_nucleotide",
        "invalid_length",
        "internal_stop",
        "invalid_start",
        "ambiguous_translation",
    }:
        return str(validity)
    if _code(row) in (None, "unknown"):
        return "unknown_genetic_code"
    return None


def _simple_path(row: dict) -> bool:
    if row.get("segments_contiguous") is False:
        return False
    segments = _geometry(row)["segments"]
    if len(segments) == 1:
        return True
    if len(segments) != 2 or not row.get("wraps_origin"):
        return False
    if row.get("strand") == 1:
        return segments[1]["start0"] == 0
    return segments[0]["start0"] == 0


def _stop_key(row: dict) -> tuple | None:
    stop = row.get("stop_triplet")
    if (
        _invalid_reason(row)
        or row.get("partial_5p")
        or row.get("partial_3p")
        or not _simple_path(row)
        or not stop
        or len(stop) != 3
        or len(set(stop)) != 3
    ):
        return None
    translated = row.get("translated_nt_length", _length(row) - int(row.get("codon_start", 1)) + 1)
    if translated % 3:
        return None
    if list(stop) != stop_positions(row["segments"], row["strand"], row.get("codon_start", 1)):
        return None
    return (_sha(row), row["strand"], tuple(int(x) for x in stop), _code(row))


def _first_position(row: dict) -> int:
    segment = row["segments"][0]
    first = segment["start0"] if row["strand"] == 1 else segment["end0"] - 1
    return int(first) + (int(row.get("codon_start", 1)) - 1) * row["strand"]


def _frame(a: dict, b: dict) -> int | None:
    """Reading offset at a shared nucleotide, including joined locations."""
    if a.get("strand") != b.get("strand") or a.get("strand") not in (1, -1):
        return None

    def offset(row: dict, position: int) -> int:
        previous = 0
        for segment in row["segments"]:
            if segment["start0"] <= position < segment["end0"]:
                local = (
                    position - segment["start0"]
                    if row["strand"] == 1
                    else segment["end0"] - 1 - position
                )
                return (previous + local - int(row.get("codon_start", 1)) + 1) % 3
            previous += segment["end0"] - segment["start0"]
        raise ValueError("Position outside CDS")

    for sa in _union(a):
        for sb in _union(b):
            position = max(sa[0], sb[0])
            if position < min(sa[1], sb[1]):
                return (offset(a, position) - offset(b, position)) % 3
    return None


class _IntervalIndex:
    """Start-sorted interval index with prefix maxima, one entry per segment."""

    def __init__(self, rows: list[dict]):
        self.entries = sorted(
            (start, end, index) for index, row in enumerate(rows) for start, end in _union(row)
        )
        self.starts = [s for s, _, _ in self.entries]
        self.max_ends: list[int] = []
        maximum = -1
        for _, end, _ in self.entries:
            maximum = max(maximum, end)
            self.max_ends.append(maximum)

    def query(self, row: dict) -> set[int]:
        found: set[int] = set()
        for start, end in _union(row):
            i = bisect_left(self.starts, end) - 1
            while i >= 0 and self.max_ends[i] > start:
                if self.entries[i][1] > start:
                    found.add(self.entries[i][2])
                i -= 1
        return found


def _prepare(
    calls: Iterable[dict], source_runs: Iterable[dict], *, include_rotated: bool = False
) -> tuple[list[dict], dict[tuple[str, str], list[dict]]]:
    calls = [dict(c) for c in calls]
    by_run: dict[str, list[dict]] = defaultdict(list)
    by_genome: dict[str, set[str]] = defaultdict(set)
    for call in calls:
        by_run[str(call.get("source_run_id"))].append(call)
        by_genome[str(call["genome_id"])].add(_sha(call))
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for original in source_runs:
        run = dict(original)
        if not include_rotated and not _is_r0(run):
            continue
        run_id = str(run.get("source_run_id") or _stable("source_run_", run))
        run["source_run_id"] = run_id
        genome = str(run["genome_id"])
        shas = {_sha(c) for c in by_run[run_id]}
        sequence_sha = _sha(run) or (next(iter(shas)) if len(shas) == 1 else "")
        if not sequence_sha and len(by_genome[genome]) == 1:
            sequence_sha = next(iter(by_genome[genome]))
        run["sequence_sha256"] = sequence_sha
        groups[(genome, sequence_sha)].append(run)
    for runs in groups.values():
        runs.sort(key=lambda r: (str(r.get("source", "")), r["source_run_id"]))
    return calls, groups


def _pair_reason(a: dict | None, b: dict | None, ra: dict, rb: dict) -> str | None:
    for run in (ra, rb):
        if run.get("comparison_status") == "not_comparable":
            return run.get("comparison_reason") or "genome_primary_code_review"
        if run.get("status") not in SUCCESS:
            return f"source_{run.get('source')}_{run.get('status', 'not_run')}"
    for call in (a, b):
        if call is not None and (reason := _invalid_reason(call)):
            return reason
    ca, cb = _code(a or ra), _code(b or rb)
    if ca and cb and ca != cb:
        return "incompatible_genetic_codes"
    return None


def _match_row(
    a: dict | None, b: dict | None, ra: dict, rb: dict, threshold: float, reason: str | None = None
) -> dict:
    overlap = _overlap(a, b) if a and b else 0
    la, lb = _length(a) if a else 0, _length(b) if b else 0
    ca, cb = (overlap / la if la else None), (overlap / lb if lb else None)
    frame = (
        _frame(a, b)
        if a and b and overlap and not _invalid_reason(a) and not _invalid_reason(b)
        else None
    )
    reason = reason or _pair_reason(a, b, ra, rb)
    nested = bool(a and b and overlap and overlap == min(la, lb))
    substantial = bool(
        a
        and b
        and la
        and lb
        and a.get("strand") == b.get("strand")
        and ca >= threshold
        and cb >= threshold
    )
    exact = bool(a and b and _geometry(a) == _geometry(b))
    same_stop = bool(
        a and b and _stop_key(a) is not None and _stop_key(a) == _stop_key(b) and frame == 0
    )
    if reason:
        relation = "not_comparable"
    elif a is None or b is None:
        relation = "no_correspondence"
    elif a["strand"] != b["strand"]:
        relation = "antisense_overlap"
    elif exact and (not _stop_key(a) or not _stop_key(b) or _stop_key(a) == _stop_key(b)):
        relation = "exact_cds"
    elif same_stop and _first_position(a) != _first_position(b):
        relation = "alternative_start"
    elif nested:
        relation = "contained_distinct_stop"
    else:
        relation = "partial_overlap"
    value = {
        "genome_id": ra["genome_id"],
        "sequence_sha256": _sha(ra),
        "source_a": ra.get("source"),
        "source_b": rb.get("source"),
        "source_run_id_a": ra["source_run_id"],
        "source_run_id_b": rb["source_run_id"],
        "call_id_a": a.get("call_id") if a else None,
        "call_id_b": b.get("call_id") if b else None,
        "model_id_a": _model_id(a) if a else None,
        "model_id_b": _model_id(b) if b else None,
        "relation": relation,
        "nested": nested,
        "overlap_nt": overlap,
        "coverage_a": ca,
        "coverage_b": cb,
        "relative_frame": frame,
        "substantial_overlap": substantial,
        "comparable": reason is None,
        "reason": reason,
        "protein_translation_mismatch": bool(
            a
            and b
            and a.get("recomputed_translation")
            and b.get("recomputed_translation")
            and a["recomputed_translation"] != b["recomputed_translation"]
        ),
    }
    value["match_id"] = _stable("match_", value)
    return value


def compare_calls(
    calls: Iterable[dict], source_runs: Iterable[dict], reciprocal_threshold: float = 0.8
) -> list[dict]:
    """Return every r0 cross-source intersecting pair and directional absences.

    A small same-strand intersection or antisense intersection is retained AND
    receives a separate ``no_correspondence`` row when no suitable agreement
    exists.  A failed/missing run receives ``not_comparable`` instead.
    """
    if not 0 <= reciprocal_threshold <= 1:
        raise ValueError("reciprocal_threshold must be between zero and one")
    calls, groups = _prepare(calls, source_runs)
    by_run: dict[str, list[dict]] = defaultdict(list)
    for call in calls:
        if _is_r0(call):
            by_run[str(call.get("source_run_id"))].append(call)
    output: list[dict] = []
    for (genome, sha), runs in sorted(groups.items()):
        for ra, rb in combinations(runs, 2):
            if ra.get("source") == rb.get("source"):
                continue
            aa = sorted(
                (
                    c
                    for c in by_run[ra["source_run_id"]]
                    if str(c["genome_id"]) == genome and _sha(c) == sha
                ),
                key=lambda c: str(c.get("call_id")),
            )
            bb = sorted(
                (
                    c
                    for c in by_run[rb["source_run_id"]]
                    if str(c["genome_id"]) == genome and _sha(c) == sha
                ),
                key=lambda c: str(c.get("call_id")),
            )
            index = _IntervalIndex(bb)
            appropriate_a: set[int] = set()
            appropriate_b: set[int] = set()
            incomparable_a: dict[int, str] = {}
            incomparable_b: dict[int, str] = {}
            for ia, a in enumerate(aa):
                for ib in sorted(index.query(a)):
                    row = _match_row(a, bb[ib], ra, rb, reciprocal_threshold)
                    output.append(row)
                    if row["comparable"] and (
                        row["relation"] in {"exact_cds", "alternative_start"}
                        or row["substantial_overlap"]
                    ):
                        appropriate_a.add(ia)
                        appropriate_b.add(ib)
                    if not row["comparable"]:
                        incomparable_a[ia] = row["reason"]
                        incomparable_b[ib] = row["reason"]
            for rows, other, own_run, other_run, appropriate, incomparable, reverse in (
                (aa, bb, ra, rb, appropriate_a, incomparable_a, False),
                (bb, aa, rb, ra, appropriate_b, incomparable_b, True),
            ):
                other_codes = {_code(c) for c in other if _code(c)}
                for i, call in enumerate(rows):
                    if i in appropriate:
                        continue
                    reason = _pair_reason(call, None, own_run, other_run)
                    if not reason and other_codes and _code(call) not in other_codes:
                        reason = "incompatible_genetic_codes"
                    reason = reason or incomparable.get(i)
                    row = (
                        _match_row(None, call, ra, rb, reciprocal_threshold, reason)
                        if reverse
                        else _match_row(call, None, ra, rb, reciprocal_threshold, reason)
                    )
                    output.append(row)
            if not aa and not bb:
                row = _match_row(None, None, ra, rb, reciprocal_threshold)
                if row["comparable"]:
                    row["reason"] = "both_sources_zero_calls"
                    row["match_id"] = _stable(
                        "match_", {k: v for k, v in row.items() if k != "match_id"}
                    )
                output.append(row)
    return sorted(
        output,
        key=lambda r: (
            str(r["genome_id"]),
            r["source_run_id_a"],
            r["source_run_id_b"],
            str(r["call_id_a"]),
            str(r["call_id_b"]),
            r["relation"],
        ),
    )


def build_models(calls: Iterable[dict]) -> dict[str, list[dict]]:
    """Deduplicate coordinate models while keeping every original call link."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for call in calls:
        groups[_model_id(call)].append(call)
    models, sources, segments = [], [], []
    seen_calls: set[str] = set()
    for model_id, rows in sorted(groups.items()):
        rows.sort(key=lambda r: str(r.get("call_id")))
        canonical = _geometry(rows[0])
        if any(_geometry(row) != canonical for row in rows[1:]):
            raise ValueError(f"Conflicting geometry for model_id {model_id}")
        models.append(
            {
                "model_id": model_id,
                **canonical,
                "genome_id": rows[0]["genome_id"],
                "genome_ids": sorted({str(r["genome_id"]) for r in rows}),
                "feature_nt_length": _length(rows[0]),
                "source_call_count": len(rows),
                "translation_tables": sorted(
                    {
                        int(r["transl_table_effective"])
                        for r in rows
                        if r.get("transl_table_effective") is not None
                    }
                ),
                "has_partial_source": any(r.get("partial_5p") or r.get("partial_3p") for r in rows),
                "has_primary_eligible_source": any(r.get("eligible_primary", False) for r in rows),
            }
        )
        for order, segment in enumerate(canonical["segments"]):
            segments.append({"model_id": model_id, "segment_order": order, **segment})
        for row in rows:
            call_id = str(row["call_id"])
            if call_id in seen_calls:
                raise ValueError(f"Duplicate call_id {call_id}")
            seen_calls.add(call_id)
            sources.append(
                {
                    "call_id": row["call_id"],
                    "model_id": model_id,
                    "source_run_id": row["source_run_id"],
                    "source": row.get("source"),
                    "genome_id": row["genome_id"],
                    "feature_id": row.get("feature_id"),
                }
            )
    return {"models": models, "model_sources": sources, "model_segments": segments}


def build_loci(calls: Iterable[dict], source_runs: Iterable[dict]) -> dict[str, list[dict]]:
    """Group complete compatible CDS by their exact ordered terminal stop.

    Raw support describes observed r0 calls, not evidence that an ORF is a gene.
    Primary support is assessable only when both current independent callers
    completed (including successful zero calls) and every defining caller call
    is admissible. An excluded contribution must not become biological absence.
    """
    calls, groups = _prepare(calls, source_runs)
    run_map = {r["source_run_id"]: r for runs in groups.values() for r in runs}
    by_genome: dict[tuple[str, str], list[dict]] = defaultdict(list)
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for call in calls:
        if not _is_r0(call):
            continue
        genome_key = (str(call["genome_id"]), _sha(call))
        by_genome[genome_key].append(call)
        if (key := _stop_key(call)) is not None:
            grouped[(genome_key, key)].append(call)
    baseline_context = {}
    for genome_key, genome_calls in by_genome.items():
        baseline = [r for r in genome_calls if r.get("source") == "genbank"]
        baseline_context[genome_key] = (
            baseline,
            _IntervalIndex(baseline),
            {json.dumps(_geometry(r), sort_keys=True) for r in baseline},
            {_stop_key(r) for r in baseline if _stop_key(r) is not None},
        )
    loci, links = [], []
    for (genome_key, key), rows in sorted(grouped.items(), key=lambda item: str(item[0])):
        genome, sha = genome_key
        locus_id = _stable("locus_", [genome, *key])
        runs = groups.get(genome_key, [])
        current_by_family: dict[str, list[dict]] = defaultdict(list)
        for run in runs:
            # A Prokka import shares the Prodigal family, but is not a current
            # Pyrodigal-gv run and must not substitute for its execution.
            source = str(run.get("source", "")).lower().replace("-", "_")
            family = _family(run.get("algorithm_family", source))
            if source in CALLER_FAMILIES and family in CALLER_FAMILIES:
                current_by_family[family].append(run)
        completed = {
            f: bool(current_by_family[f])
            and all(
                r.get("status") in SUCCESS and r.get("comparison_status") != "not_comparable"
                for r in current_by_family[f]
            )
            for f in CALLER_FAMILIES
        }

        def current_caller(row: dict) -> bool:
            run = run_map.get(str(row.get("source_run_id")))
            return bool(
                run
                and str(run.get("source", "")).lower().replace("-", "_") in CALLER_FAMILIES
                and _family(run.get("algorithm_family", run.get("source"))) in CALLER_FAMILIES
                and run.get("status") in SUCCESS
            )

        caller_rows = [r for r in rows if current_caller(r)]
        families = sorted(
            {
                _family(run_map[str(r["source_run_id"])].get("algorithm_family", r.get("source")))
                for r in caller_rows
            }
            & set(CALLER_FAMILIES)
        )

        def exclusion_reasons(row: dict) -> list[str]:
            """Use the same admissible domain for representatives and support."""
            run = run_map[str(row["source_run_id"])]
            reasons = []
            if run.get("comparison_status") == "not_comparable":
                reasons.extend(
                    str(run.get("comparison_reason") or "genome_primary_code_review").split(";")
                )
            if row.get("wraps_origin"):
                reasons.append("wraps_origin")
            if row.get("rotation_sensitive"):
                reasons.append("rotation_sensitive")
            if not row.get("eligible_primary", False):
                recorded = list(row.get("primary_exclusion_reasons") or [])
                recorded.extend(str(row.get("primary_exclusion_reason") or "").split(";"))
                reasons.extend(reason for reason in recorded if reason)
                if not reasons:
                    reasons.append("caller_model_not_eligible")
            return list(dict.fromkeys(reasons))

        eligible = []
        excluded_contributions = []
        for row in caller_rows:
            row_reasons = exclusion_reasons(row)
            if row_reasons:
                family = _family(
                    run_map[str(row["source_run_id"])].get("algorithm_family", row.get("source"))
                )
                excluded_contributions.extend(
                    f"support_source_excluded:{family}:{reason}" for reason in row_reasons
                )
            else:
                eligible.append(row)
        eligible_families = sorted(
            {
                _family(run_map[str(r["source_run_id"])].get("algorithm_family", r.get("source")))
                for r in eligible
            }
            & set(CALLER_FAMILIES)
        )
        eligible.sort(
            key=lambda r: (
                -int(r.get("translated_nt_length", _length(r))),
                _model_id(r),
                str(r.get("call_id")),
            )
        )
        baseline, baseline_index, baseline_geometries, baseline_stops = baseline_context[genome_key]
        baseline_same_stop = key in baseline_stops
        baseline_exact = any(
            json.dumps(_geometry(c), sort_keys=True) in baseline_geometries for c in caller_rows
        )
        context_indices: set[int] = set()
        for caller in caller_rows:
            context_indices.update(baseline_index.query(caller))
        context = [baseline[i] for i in sorted(context_indices)]
        same_strand = [b for b in context if b.get("strand") == key[1]]
        antisense = any(b.get("strand") != key[1] for b in context)
        baseline_runs = [r for r in runs if r.get("source") == "genbank"]
        baseline_comparable = bool(baseline_runs) and all(
            r.get("status") in SUCCESS and r.get("comparison_status") != "not_comparable"
            for r in baseline_runs
        )
        reasons = []
        if not families:
            reasons.append(
                "genbank_only"
                if any(r.get("source") == "genbank" for r in rows)
                else "no_current_caller"
            )
        for family in CALLER_FAMILIES:
            if not completed[family]:
                gated = [
                    r
                    for r in current_by_family[family]
                    if r.get("comparison_status") == "not_comparable"
                ]
                if gated:
                    reasons.extend(
                        r.get("comparison_reason") or "genome_primary_code_review" for r in gated
                    )
                else:
                    reasons.append(f"{family}_not_successful")
        if not eligible:
            reasons.append("no_eligible_caller_model")
        reasons.extend(sorted(set(excluded_contributions)))
        if key[3] != "standard_stop":
            reasons.append("nonstandard_genetic_code")
        raw_mask = sum(1 << CALLER_FAMILIES.index(f) for f in families)
        eligible_mask = sum(1 << CALLER_FAMILIES.index(f) for f in eligible_families)
        locus = {
            "locus_id": locus_id,
            "genome_id": genome,
            "sequence_sha256": sha,
            "strand": key[1],
            "stop_triplet": list(key[2]),
            "code_class": key[3],
            "representative_model_id": _model_id(eligible[0]) if eligible else None,
            "representative_call_id": eligible[0].get("call_id") if eligible else None,
            "representative_rule": "longest_eligible_current_caller_then_model_id",
            # Compatibility fields retain observed raw support. Primary analyses
            # must first require support_assessable_primary and eligible_primary.
            "caller_support": len(families),
            "caller_support_raw": len(families),
            "caller_support_eligible": len(eligible_families),
            "caller_support_mask": raw_mask,
            "caller_support_mask_raw": raw_mask,
            "caller_support_mask_eligible": eligible_mask,
            "caller_families": families,
            "caller_families_eligible": eligible_families,
            "support_assessable_primary": not reasons,
            "source_presence": sorted({str(r.get("source")) for r in rows}),
            "phanotate_present": "phanotate" in families,
            "pyrodigal_gv_present": "pyrodigal_gv" in families,
            "phanotate_completed": completed["phanotate"],
            "pyrodigal_gv_completed": completed["pyrodigal_gv"],
            "genbank_present": any(r.get("source") == "genbank" for r in rows),
            "genbank_exact_match": baseline_exact,
            "genbank_stop_match": bool(baseline_same_stop),
            "genbank_comparable": baseline_comparable,
            "universe": "caller_union" if families else "baseline_only",
            "additional_stop_model": (not baseline_exact and not baseline_same_stop)
            if families and baseline_comparable
            else None,
            "unannotated_same_strand_region": not bool(same_strand)
            if families and baseline_comparable
            else None,
            "genbank_antisense_context": antisense,
            "genbank_partial_context": any(
                b.get("partial_5p") or b.get("partial_3p") for b in same_strand
            ),
            "genbank_complex_context": any(
                b.get("pseudo") or b.get("exception") or b.get("transl_except") for b in same_strand
            ),
            "eligible_primary": not reasons,
            "primary_exclusion_reason": ";".join(dict.fromkeys(reasons)) or None,
        }
        loci.append(locus)
        model_rows: dict[str, list[dict]] = defaultdict(list)
        for row in rows:
            model_rows[_model_id(row)].append(row)
        for model_id, origins in sorted(model_rows.items()):
            links.append(
                {
                    "locus_id": locus_id,
                    "model_id": model_id,
                    "call_ids": sorted({str(r["call_id"]) for r in origins}),
                    "source_run_ids": sorted({str(r["source_run_id"]) for r in origins}),
                    "sources": sorted({str(r.get("source")) for r in origins}),
                    "genome_id": genome,
                    "is_representative": bool(eligible and model_id == _model_id(eligible[0])),
                }
            )
    return {"stop_loci": loci, "locus_models": links}


def _track_metrics(rows: list[dict], success: bool) -> dict:
    """Within-source overlap, excluding alternative starts at the same stop."""
    by_model = {}
    for row in sorted(rows, key=lambda r: str(r.get("call_id"))):
        by_model.setdefault(_model_id(row), row)
    models = list(by_model.values())
    lengths = [_length(r) for r in models if _length(r)]
    aa_lengths = []
    for row in models:
        if isinstance(row.get("recomputed_translation"), str):
            aa_lengths.append(len(row["recomputed_translation"]))
        elif row.get("translated_nt_length") is not None:
            aa_lengths.append(
                int(row["translated_nt_length"]) // 3
                - int(bool(row.get("terminal_stop_present") or row.get("stop_triplet")))
            )
    result = {
        "length_defined_model_count": len(lengths),
        "length_nt_min": min(lengths) if lengths else None,
        "length_nt_median": float(median(lengths)) if lengths else None,
        "length_nt_max": max(lengths) if lengths else None,
        "protein_length_defined_model_count": len(aa_lengths),
        "length_aa_min": min(aa_lengths) if aa_lengths else None,
        "length_aa_median": float(median(aa_lengths)) if aa_lengths else None,
        "length_aa_max": max(aa_lengths) if aa_lengths else None,
        "short_protein_threshold_aa": 100,
        "short_protein_model_count": sum(n < 100 for n in aa_lengths),
        "short_protein_fraction": sum(n < 100 for n in aa_lengths) / len(aa_lengths)
        if aa_lengths
        else None,
        "partial_call_count": sum(bool(r.get("partial_5p") or r.get("partial_3p")) for r in rows),
        "wraps_origin_call_count": sum(bool(r.get("wraps_origin")) for r in rows),
        "rotation_sensitive_call_count": sum(bool(r.get("rotation_sensitive")) for r in rows),
        "complex_biology_call_count": sum(
            bool(r.get("pseudo") or r.get("exception") or r.get("transl_except")) for r in rows
        ),
        "primary_excluded_call_count": sum(not bool(r.get("eligible_primary")) for r in rows),
    }
    index = _IntervalIndex(models)
    pairs = {
        ("same_strand", 1): 0,
        ("same_strand", 30): 0,
        ("antisense", 1): 0,
        ("antisense", 30): 0,
    }
    overlapping = {key: set() for key in pairs}
    for ia, a in enumerate(models):
        for ib in sorted(index.query(a)):
            if ib <= ia:
                continue
            b = models[ib]
            if a.get("strand") not in (1, -1) or b.get("strand") not in (1, -1):
                continue
            if _stop_key(a) is not None and _stop_key(a) == _stop_key(b):
                continue
            category = "same_strand" if a["strand"] == b["strand"] else "antisense"
            overlap = _overlap(a, b)
            for threshold in (1, 30):
                if overlap >= threshold:
                    pairs[(category, threshold)] += 1
                    overlapping[(category, threshold)].update((ia, ib))
    for (category, threshold), count in pairs.items():
        result[f"{category}_overlap_pairs_ge{threshold}nt"] = count
        result[f"{category}_overlapping_models_ge{threshold}nt"] = len(
            overlapping[(category, threshold)]
        )
    result["overlap_denominator_model_count"] = sum(
        r.get("strand") in (1, -1) and bool(_union(r)) for r in models
    )
    if not success:
        return {key: None for key in result}
    return result


def genome_metrics(
    calls: Iterable[dict],
    matches: Iterable[dict],
    loci: dict | Iterable[dict],
    source_runs: Iterable[dict],
) -> list[dict]:
    """Unique-model agreement with explicit denominators and run statuses."""
    calls, groups = _prepare(calls, source_runs, include_rotated=True)
    matches = list(matches)
    loci = loci["stop_loci"] if isinstance(loci, dict) else list(loci)
    output = []
    for (genome, sha), runs in sorted(groups.items()):

        def rows_for(run: dict) -> list[dict]:
            return [
                c
                for c in calls
                if str(c["genome_id"]) == genome
                and _sha(c) == sha
                and c.get("source_run_id") == run["source_run_id"]
            ]

        for run in runs:
            rows = rows_for(run)
            success = run.get("status") in SUCCESS
            output.append(
                {
                    "genome_id": genome,
                    "sequence_sha256": sha,
                    "metric_scope": "source",
                    "source_a": run.get("source"),
                    "source_b": None,
                    "source_run_id_a": run["source_run_id"],
                    "source_run_id_b": None,
                    "status": run.get("status", "not_run"),
                    "rotation_shift": run.get("rotation_shift", 0),
                    "cds_count": len(rows) if success else None,
                    "observed_call_count": len(rows),
                    "unique_model_count": len({_model_id(c) for c in rows}) if success else None,
                    "eligible_call_count": sum(bool(c.get("eligible_primary")) for c in rows)
                    if success
                    else None,
                    "stop_locus_count": len({_stop_key(c) for c in rows if _stop_key(c)})
                    if success
                    else None,
                    "primary_locus_count": sum(
                        locus["genome_id"] == genome and locus.get("eligible_primary", False)
                        for locus in loci
                    ),
                    "denominator_definition": "one source/rotation track; failed counts null; models deduplicated by canonical geometry; within-track overlaps exclude shared-stop alternative starts; short threshold strictly <100 aa",
                    **_track_metrics(rows, success),
                }
            )
        for ra, rb in combinations(runs, 2):
            if ra.get("source") == rb.get("source") or not _is_r0(ra) or not _is_r0(rb):
                continue
            aa, bb = rows_for(ra), rows_for(rb)
            reason = _pair_reason(None, None, ra, rb)
            code_a = {_code(c) for c in aa if not _invalid_reason(c)}
            code_b = {_code(c) for c in bb if not _invalid_reason(c)}
            if not reason and code_a and code_b and code_a != code_b:
                reason = "incompatible_genetic_codes"
            row = {
                "genome_id": genome,
                "sequence_sha256": sha,
                "metric_scope": "source_pair",
                "source_a": ra.get("source"),
                "source_b": rb.get("source"),
                "source_run_id_a": ra["source_run_id"],
                "source_run_id_b": rb["source_run_id"],
                "status": "not_comparable" if reason else "success",
                "reason": reason,
                "denominator_definition": "unique comparable r0 coordinate models for exact; full simple complete stop keys for stop; baseline agreement is not accuracy",
            }
            for name, key_function in (
                ("exact", lambda c: json.dumps(_geometry(c), sort_keys=True)),
                ("stop", _stop_key),
            ):
                sets = [
                    {
                        key_function(c)
                        for c in rows
                        if not _invalid_reason(c) and key_function(c) is not None
                    }
                    for rows in (aa, bb)
                ]
                a_set, b_set = sets
                intersection, union = len(a_set & b_set), len(a_set | b_set)
                denominator = len(a_set) + len(b_set)
                row.update(
                    {
                        f"{name}_a_count": len(a_set) if not reason else None,
                        f"{name}_b_count": len(b_set) if not reason else None,
                        f"{name}_intersection": intersection if not reason else None,
                        f"{name}_union": union if not reason else None,
                        f"{name}_jaccard": intersection / union if union and not reason else None,
                        f"{name}_dice": 2 * intersection / denominator
                        if denominator and not reason
                        else None,
                        f"{name}_denominator_status": "not_comparable"
                        if reason
                        else ("zero_denominator" if union == 0 else "defined"),
                    }
                )
            pair_matches = [
                m
                for m in matches
                if m.get("source_run_id_a") == ra["source_run_id"]
                and m.get("source_run_id_b") == rb["source_run_id"]
            ]
            for relation in (
                "exact_cds",
                "alternative_start",
                "contained_distinct_stop",
                "partial_overlap",
                "antisense_overlap",
                "no_correspondence",
                "not_comparable",
            ):
                row[relation + "_count"] = sum(m.get("relation") == relation for m in pair_matches)
            row["substantial_overlap_count"] = sum(
                m.get("comparable", False) and m.get("substantial_overlap", False)
                for m in pair_matches
            )
            output.append(row)
    return output
