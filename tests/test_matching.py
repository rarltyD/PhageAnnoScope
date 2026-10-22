"""Manually calculated synthetic coordinates; no research data are generated."""

import hashlib
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from phageannoscope.matching import build_loci, build_models, compare_calls, genome_metrics

SHA = hashlib.sha256(b"explicitly-synthetic-coordinate-test").hexdigest()


def run(source, status="success", **extra):
    return {
        "source_run_id": source + "_run",
        "genome_id": "synthetic.1",
        "source": source,
        "algorithm_family": source,
        "sequence_sha256": SHA,
        "status": status,
        "rotation_shift": 0,
        **extra,
    }


def call(identifier, source, start, end, strand=1, **extra):
    stop = [end - 3, end - 2, end - 1] if strand == 1 else [start + 2, start + 1, start]
    return {
        "call_id": identifier,
        "feature_id": "feature_" + identifier,
        "source_run_id": source + "_run",
        "genome_id": "synthetic.1",
        "source": source,
        "algorithm_family": source,
        "sequence_sha256": SHA,
        "segments": [{"start0": start, "end0": end}],
        "strand": strand,
        "codon_start": 1,
        "transl_table_effective": 11,
        "code_class": "standard_stop",
        "feature_nt_length": end - start,
        "translated_nt_length": end - start,
        "stop_triplet": stop,
        "partial_5p": False,
        "partial_3p": False,
        "pseudo": False,
        "exception": None,
        "transl_except": None,
        "wraps_origin": False,
        "rotation_shift": 0,
        "rotation_sensitive": False,
        "eligible_primary": True,
        "translation_validity": "valid",
        **extra,
    }


def pair_rows(calls, threshold=0.8):
    return compare_calls(calls, [run("genbank"), run("phanotate")], threshold)


def overlap_rows(rows):
    return [r for r in rows if r["call_id_a"] is not None and r["call_id_b"] is not None]


@pytest.mark.parametrize("strand", [1, -1])
def test_exact_and_alternative_start_on_both_strands(strand):
    a = call("a", "genbank", 30, 150, strand)
    b = call("b", "phanotate", 30, 150, strand, transl_table_effective=1)
    assert overlap_rows(pair_rows([a, b]))[0]["relation"] == "exact_cds"
    b = call("b", "phanotate", 60 if strand == 1 else 30, 150 if strand == 1 else 120, strand)
    result = overlap_rows(pair_rows([a, b]))[0]
    assert result["relation"] == "alternative_start"
    assert result["relative_frame"] == 0
    assert result["nested"] is True


def test_every_split_merge_overlap_is_retained():
    a = call("a", "genbank", 0, 120)
    b = call("b", "phanotate", 3, 63)
    c = call("c", "phanotate", 63, 123)
    result = overlap_rows(pair_rows([a, b, c]))
    assert {(r["call_id_a"], r["call_id_b"]) for r in result} == {("a", "b"), ("a", "c")}
    assert {r["relation"] for r in result} == {"contained_distinct_stop", "partial_overlap"}


@pytest.mark.parametrize("offset,expected", [(24, True), (25, False)])
def test_reciprocal_threshold_includes_point_eight(offset, expected):
    a = call("a", "genbank", 0, 120)
    b = call("b", "phanotate", offset, 120 + offset)
    rows = pair_rows([a, b])
    match = overlap_rows(rows)[0]
    assert match["overlap_nt"] == 120 - offset
    assert match["coverage_a"] == (120 - offset) / 120
    assert match["coverage_b"] == (120 - offset) / 120
    assert match["substantial_overlap"] is expected
    assert sum(r["relation"] == "no_correspondence" for r in rows) == (0 if expected else 2)


def test_unequal_coverage_and_antisense_are_context():
    a = call("a", "genbank", 0, 120)
    b = call("b", "phanotate", 60, 240)
    result = overlap_rows(pair_rows([a, b]))[0]
    assert result["coverage_a"] == 0.5
    assert result["coverage_b"] == pytest.approx(1 / 3)
    b["strand"] = -1
    result = pair_rows([a, b])
    assert overlap_rows(result)[0]["relation"] == "antisense_overlap"
    assert sum(r["relation"] == "no_correspondence" for r in result) == 2


def test_origin_segments_use_union_and_translation_order():
    a = call(
        "a",
        "genbank",
        0,
        90,
        wraps_origin=True,
        segments=[{"start0": 960, "end0": 1000}, {"start0": 0, "end0": 50}],
        stop_triplet=[47, 48, 49],
    )
    b = call("b", "phanotate", 300, 600)
    assert not overlap_rows(pair_rows([a, b]))
    c = call("c", "phanotate", 0, 30)
    match = overlap_rows(pair_rows([a, c]))[0]
    assert match["overlap_nt"] == 30
    assert match["coverage_a"] == pytest.approx(1 / 3)


def test_different_code_and_complex_features_are_not_comparable():
    a = call("a", "genbank", 0, 120)
    b = call("b", "phanotate", 0, 120, code_class="table_15", transl_table_effective=15)
    assert {r["relation"] for r in pair_rows([a, b])} == {"not_comparable"}
    b = call("b", "phanotate", 0, 120, transl_except="(pos:7..9,aa:Sec)")
    assert {r["relation"] for r in pair_rows([a, b])} == {"not_comparable"}


@pytest.mark.parametrize("status", ["failed", "timeout", "not_run", "unsupported_code"])
def test_failed_caller_is_not_zero_predictions(status):
    calls = [call("a", "phanotate", 0, 120)]
    runs = [run("phanotate"), run("pyrodigal_gv", status)]
    matches = compare_calls(calls, runs)
    assert len(matches) == 1
    assert matches[0]["relation"] == "not_comparable"
    loci = build_loci(calls, runs)
    assert loci["stop_loci"][0]["caller_support"] == 1
    assert loci["stop_loci"][0]["caller_support_raw"] == 1
    assert loci["stop_loci"][0]["caller_support_eligible"] == 1
    assert loci["stop_loci"][0]["support_assessable_primary"] is False
    assert loci["stop_loci"][0]["eligible_primary"] is False
    metrics = genome_metrics(calls, matches, loci, runs)
    assert (
        next(
            r for r in metrics if r["metric_scope"] == "source" and r["source_a"] == "pyrodigal_gv"
        )["cds_count"]
        is None
    )
    assert next(r for r in metrics if r["metric_scope"] == "source_pair")["stop_jaccard"] is None


def test_successful_zero_calls_is_valid_single_support_and_empty_denominator():
    calls = [call("a", "phanotate", 0, 120)]
    runs = [run("phanotate"), run("pyrodigal_gv", "cache_hit")]
    matches = compare_calls(calls, runs)
    assert matches[0]["relation"] == "no_correspondence"
    loci = build_loci(calls, runs)
    assert loci["stop_loci"][0]["eligible_primary"] is True
    assert loci["stop_loci"][0]["caller_support_raw"] == 1
    assert loci["stop_loci"][0]["caller_support_eligible"] == 1
    assert loci["stop_loci"][0]["support_assessable_primary"] is True
    empty = genome_metrics([], compare_calls([], runs), build_loci([], runs), runs)
    pair = next(r for r in empty if r["metric_scope"] == "source_pair")
    assert pair["status"] == "success"
    assert pair["exact_a_count"] == pair["exact_b_count"] == 0
    assert pair["exact_jaccard"] is None
    assert pair["stop_dice"] is None
    assert pair["stop_denominator_status"] == "zero_denominator"


def test_no_overlap_transitivity_and_prokka_not_extra_vote():
    calls = [
        call("a", "phanotate", 0, 120),
        call("b", "pyrodigal_gv", 60, 180),
        call("c", "phanotate", 120, 240),
        call("d", "prokka", 0, 120, algorithm_family="prodigal"),
    ]
    runs = [run(s) for s in ("phanotate", "pyrodigal_gv", "prokka")]
    loci = build_loci(calls, runs)["stop_loci"]
    assert len(loci) == 3
    assert all(locus["caller_support"] == 1 for locus in loci)


def test_alternative_starts_representative_and_genbank_additional_flags():
    calls = [
        call("a", "phanotate", 30, 150),
        call("b", "pyrodigal_gv", 0, 150),
        call("baseline", "genbank", 60, 180),
    ]
    runs = [run(s) for s in ("phanotate", "pyrodigal_gv", "genbank")]
    result = build_loci(calls, runs)
    locus = next(row for row in result["stop_loci"] if row["caller_support"])
    assert locus["caller_support"] == 2
    assert locus["representative_call_id"] == "b"
    assert locus["additional_stop_model"] is True
    assert locus["unannotated_same_strand_region"] is False
    assert len([r for r in result["locus_models"] if r["locus_id"] == locus["locus_id"]]) == 2
    baseline = next(row for row in result["stop_loci"] if not row["caller_support"])
    assert baseline["universe"] == "baseline_only"
    assert baseline["eligible_primary"] is False
    assert baseline["support_assessable_primary"] is False
    assert baseline["caller_support_raw"] == baseline["caller_support_eligible"] == 0


def test_duplicates_and_rotations_do_not_add_support():
    calls = [
        call("a", "phanotate", 0, 120),
        call("b", "phanotate", 0, 120),
        call("c", "pyrodigal_gv", 0, 120, rotation_shift=500, source_run_id="rotated"),
    ]
    runs = [
        run("phanotate"),
        run("pyrodigal_gv"),
        run("pyrodigal_gv", rotation_shift=500, source_run_id="rotated"),
    ]
    locus = build_loci(calls, runs)["stop_loci"][0]
    assert locus["caller_support"] == 1
    assert locus["representative_call_id"] == "a"
    assert all(r["call_id_b"] != "c" for r in compare_calls(calls, runs))


def test_partial_and_no_stop_are_not_stop_loci():
    calls = [
        call("a", "phanotate", 0, 120, partial_3p=True),
        call("b", "pyrodigal_gv", 0, 120, stop_triplet=None),
    ]
    runs = [run("phanotate"), run("pyrodigal_gv")]
    assert build_loci(calls, runs)["stop_loci"] == []
    assert overlap_rows(compare_calls(calls, runs))[0]["relation"] == "exact_cds"


def test_models_and_metrics_are_deterministic_without_losing_sources():
    calls = [call("a", "phanotate", 0, 120), call("b", "pyrodigal_gv", 0, 120)]
    runs = [run("phanotate"), run("pyrodigal_gv")]
    models = build_models(calls)
    assert len(models["models"]) == len(models["model_segments"]) == 1
    assert len(models["model_sources"]) == 2
    assert models == build_models(reversed(calls))
    matches = compare_calls(calls, runs)
    loci = build_loci(calls, runs)
    assert matches == compare_calls(reversed(calls), reversed(runs))
    assert loci == build_loci(reversed(calls), reversed(runs))
    pair = next(
        r for r in genome_metrics(calls, matches, loci, runs) if r["metric_scope"] == "source_pair"
    )
    assert pair["exact_jaccard"] == pair["stop_dice"] == 1


def test_duplicate_call_id_and_model_geometry_conflict_rejected():
    a = call("a", "phanotate", 0, 120, model_id="same")
    b = call("b", "pyrodigal_gv", 3, 120, model_id="same")
    with pytest.raises(ValueError, match="Conflicting geometry"):
        build_models([a, b])
    with pytest.raises(ValueError, match="Duplicate call_id"):
        build_models([a, a])


@pytest.mark.parametrize(
    "strand,segments",
    [
        (1, [{"start0": 0, "end0": 60}, {"start0": 60, "end0": 120}]),
        (-1, [{"start0": 60, "end0": 120}, {"start0": 0, "end0": 60}]),
    ],
)
def test_equivalent_contiguous_segmentation_has_one_canonical_model(strand, segments):
    a = call("a", "phanotate", 0, 120, strand)
    b = call(
        "b",
        "pyrodigal_gv",
        0,
        120,
        strand,
        segments=segments,
        location_operator="join",
        segments_contiguous=True,
    )
    models = build_models([a, b])
    assert len(models["models"]) == 1
    assert models["models"][0]["segments"] == [{"start0": 0, "end0": 120}]
    runs = [run("phanotate"), run("pyrodigal_gv")]
    assert len(build_loci([a, b], runs)["locus_models"]) == 1
    assert overlap_rows(compare_calls([a, b], runs))[0]["relation"] == "exact_cds"


def test_invalid_codon_start_is_preserved_as_model_without_comparison():
    a = call("a", "phanotate", 0, 120, codon_start=None, translation_validity="invalid_codon_start")
    b = call("b", "pyrodigal_gv", 0, 120)
    assert len(build_models([a, b])["models"]) == 2
    runs = [run("phanotate"), run("pyrodigal_gv")]
    assert {r["relation"] for r in compare_calls([a, b], runs)} == {"not_comparable"}


def test_declared_stop_must_match_terminal_path_and_phase():
    a = call("a", "phanotate", 0, 120, stop_triplet=[114, 115, 116])
    assert build_loci([a], [run("phanotate"), run("pyrodigal_gv")])["stop_loci"] == []


def test_within_track_overlap_excludes_artificial_alternative_start_overlap():
    calls = [
        call("a", "phanotate", 0, 120),
        call("alt", "phanotate", 30, 120),
        call("b", "phanotate", 90, 210),
        call("anti", "phanotate", 0, 30, -1),
        call("long", "phanotate", 303, 606),
        call("a_duplicate", "phanotate", 0, 120),
    ]
    runs = [run("phanotate"), run("pyrodigal_gv")]
    metrics = genome_metrics(calls, compare_calls(calls, runs), build_loci(calls, runs), runs)
    track = next(
        r for r in metrics if r["metric_scope"] == "source" and r["source_a"] == "phanotate"
    )
    assert track["cds_count"] == 6
    assert track["unique_model_count"] == 5
    assert track["length_nt_median"] == 120.0
    assert track["short_protein_model_count"] == 4
    assert track["short_protein_fraction"] == 0.8
    assert track["same_strand_overlap_pairs_ge1nt"] == 2
    assert track["same_strand_overlap_pairs_ge30nt"] == 2
    assert track["same_strand_overlapping_models_ge1nt"] == 3
    assert track["antisense_overlap_pairs_ge1nt"] == 1
    assert track["antisense_overlap_pairs_ge30nt"] == 1


def test_rotated_tracks_have_separate_metrics_but_no_primary_pair_comparison():
    calls = [
        call("a", "phanotate", 0, 120),
        call("rot", "phanotate", 0, 120, rotation_shift=500, source_run_id="rot_run"),
    ]
    runs = [
        run("phanotate"),
        run("pyrodigal_gv"),
        run("phanotate", rotation_shift=500, source_run_id="rot_run"),
    ]
    metrics = genome_metrics(calls, compare_calls(calls, runs), build_loci(calls, runs), runs)
    tracks = [r for r in metrics if r["metric_scope"] == "source" and r["source_a"] == "phanotate"]
    assert {(t["rotation_shift"], t["cds_count"]) for t in tracks} == {(0, 1), (500, 1)}
    pairs = [r for r in metrics if r["metric_scope"] == "source_pair"]
    assert len(pairs) == 1
    assert pairs[0]["source_run_id_a"] != "rot_run"


@pytest.mark.parametrize("comparison_reason", ["alternative_code_or_recoding", None])
def test_successful_runs_can_be_held_from_primary_comparison(comparison_reason):
    calls = [
        call("a", "phanotate", 0, 120),
        call("b", "pyrodigal_gv", 30, 120),
        call("gb", "genbank", 0, 120),
    ]
    runs = [
        run(source, comparison_status="not_comparable", comparison_reason=comparison_reason)
        for source in ("phanotate", "pyrodigal_gv", "genbank")
    ]
    expected_reason = comparison_reason or "genome_primary_code_review"
    matches = compare_calls(calls, runs)
    assert matches
    assert {row["relation"] for row in matches} == {"not_comparable"}
    assert {row["reason"] for row in matches} == {expected_reason}
    loci = build_loci(calls, runs)
    assert loci["stop_loci"]
    assert not any(row["eligible_primary"] for row in loci["stop_loci"])
    assert all(expected_reason in row["primary_exclusion_reason"] for row in loci["stop_loci"])
    assert not any(row["genbank_comparable"] for row in loci["stop_loci"])
    metrics = genome_metrics(calls, matches, loci, runs)
    pairs = [row for row in metrics if row["metric_scope"] == "source_pair"]
    assert all(row["status"] == "not_comparable" for row in pairs)
    assert all(row["reason"] == expected_reason for row in pairs)
    assert all(row["exact_jaccard"] is None and row["stop_dice"] is None for row in pairs)
    assert all(row["cds_count"] == 1 for row in metrics if row["metric_scope"] == "source")


def test_successful_empty_caller_with_comparison_hold_is_not_absence():
    calls = [call("a", "phanotate", 0, 120)]
    runs = [
        run("phanotate"),
        run(
            "pyrodigal_gv",
            comparison_status="not_comparable",
            comparison_reason="alternative_code_or_recoding",
        ),
    ]
    matches = compare_calls(calls, runs)
    assert len(matches) == 1
    assert matches[0]["relation"] == "not_comparable"
    assert matches[0]["reason"] == "alternative_code_or_recoding"
    assert not build_loci(calls, runs)["stop_loci"][0]["eligible_primary"]


@pytest.mark.parametrize("unstable_source", ["phanotate", "pyrodigal_gv"])
def test_excluded_support_is_uncertainty_not_single_caller_biology(unstable_source):
    calls = [
        call(
            source,
            source,
            0,
            120,
            eligible_primary=source != unstable_source,
            rotation_sensitive=source == unstable_source,
        )
        for source in ("phanotate", "pyrodigal_gv")
    ]
    result = build_loci(calls, [run("phanotate"), run("pyrodigal_gv")])
    locus = result["stop_loci"][0]
    stable_source = "pyrodigal_gv" if unstable_source == "phanotate" else "phanotate"
    assert locus["caller_support"] == locus["caller_support_raw"] == 2
    assert locus["caller_support_mask"] == locus["caller_support_mask_raw"] == 3
    assert locus["caller_support_eligible"] == 1
    assert locus["caller_families_eligible"] == [stable_source]
    assert locus["caller_support_mask_eligible"] == (2 if stable_source == "pyrodigal_gv" else 1)
    assert locus["representative_call_id"] == stable_source
    assert locus["support_assessable_primary"] is False
    assert locus["eligible_primary"] is False
    assert (
        f"support_source_excluded:{unstable_source}:rotation_sensitive"
        in locus["primary_exclusion_reason"]
    )
    assert result["locus_models"][0]["call_ids"] == ["phanotate", "pyrodigal_gv"]


def test_both_rotation_sensitive_sources_preserve_raw_support_without_representative():
    calls = [
        call(source, source, 0, 120, eligible_primary=False, rotation_sensitive=True)
        for source in ("phanotate", "pyrodigal_gv")
    ]
    locus = build_loci(calls, [run("phanotate"), run("pyrodigal_gv")])["stop_loci"][0]
    assert locus["caller_support_raw"] == 2
    assert locus["caller_support_eligible"] == locus["caller_support_mask_eligible"] == 0
    assert locus["representative_call_id"] is None
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is False
    for source in ("phanotate", "pyrodigal_gv"):
        assert (
            f"support_source_excluded:{source}:rotation_sensitive"
            in locus["primary_exclusion_reason"]
        )


def test_stable_alternative_starts_share_primary_support_domain():
    calls = [
        call("short", "phanotate", 30, 150),
        call("long", "pyrodigal_gv", 0, 150),
        # An excluded baseline model does not define current caller support.
        call("baseline", "genbank", 0, 150, eligible_primary=False, rotation_sensitive=True),
    ]
    result = build_loci(calls, [run(s) for s in ("phanotate", "pyrodigal_gv", "genbank")])
    locus = result["stop_loci"][0]
    assert locus["caller_support_raw"] == locus["caller_support_eligible"] == 2
    assert locus["caller_support_mask_raw"] == locus["caller_support_mask_eligible"] == 3
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is True
    assert locus["representative_call_id"] == "long"
    assert len(result["locus_models"]) == 2


def test_excluded_alternative_start_is_not_hidden_by_stable_model_from_same_family():
    calls = [
        call("stable_phano", "phanotate", 30, 150),
        call("unstable_phano", "phanotate", 0, 150, rotation_sensitive=True),
        call("stable_pyro", "pyrodigal_gv", 30, 150),
    ]
    result = build_loci(calls, [run("phanotate"), run("pyrodigal_gv")])
    locus = result["stop_loci"][0]
    assert locus["caller_support_raw"] == locus["caller_support_eligible"] == 2
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is False
    assert locus["representative_call_id"] != "unstable_phano"
    assert len(result["locus_models"]) == 2
    assert (
        "support_source_excluded:phanotate:rotation_sensitive" in locus["primary_exclusion_reason"]
    )


def test_origin_model_is_explicitly_excluded_even_with_stable_same_stop_caller():
    calls = [
        call(
            "origin",
            "phanotate",
            960,
            1050,
            segments=[{"start0": 960, "end0": 1000}, {"start0": 0, "end0": 50}],
            wraps_origin=True,
            stop_triplet=[47, 48, 49],
        ),
        call("stable", "pyrodigal_gv", 2, 50),
    ]
    result = build_loci(calls, [run("phanotate"), run("pyrodigal_gv")])
    locus = result["stop_loci"][0]
    assert locus["caller_support_raw"] == 2
    assert locus["caller_support_eligible"] == 1
    assert locus["representative_call_id"] == "stable"
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is False
    assert "support_source_excluded:phanotate:wraps_origin" in locus["primary_exclusion_reason"]
    assert len(result["locus_models"]) == 2


@pytest.mark.parametrize("reason", ["rotation_check_incomplete", "manual_technical_review", None])
def test_other_excluded_defining_contribution_has_explicit_reason(reason):
    calls = [
        call("held", "phanotate", 0, 120, eligible_primary=False, primary_exclusion_reason=reason),
        call("stable", "pyrodigal_gv", 0, 120),
    ]
    locus = build_loci(calls, [run("phanotate"), run("pyrodigal_gv")])["stop_loci"][0]
    assert locus["caller_support_raw"] == 2
    assert locus["caller_support_eligible"] == 1
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is False
    expected = reason or "caller_model_not_eligible"
    assert f"support_source_excluded:phanotate:{expected}" in locus["primary_exclusion_reason"]


def test_run_comparison_hold_also_limits_support_and_representative_domain():
    calls = [call("held", "phanotate", 0, 150), call("stable", "pyrodigal_gv", 30, 150)]
    runs = [
        run("phanotate", comparison_status="not_comparable", comparison_reason="technical_hold"),
        run("pyrodigal_gv"),
    ]
    locus = build_loci(calls, runs)["stop_loci"][0]
    assert locus["caller_support_raw"] == 2
    assert locus["caller_support_eligible"] == 1
    assert locus["caller_families_eligible"] == ["pyrodigal_gv"]
    assert locus["representative_call_id"] == "stable"
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is False
    assert "support_source_excluded:phanotate:technical_hold" in locus["primary_exclusion_reason"]


def test_real_phix174_original_parquet_regression_preserves_raw_links():
    """Read the frozen diagnostic as evidence; never regenerate or overwrite it."""
    root = Path(__file__).resolve().parents[1] / "results" / "diagnostic" / "parquet"
    if not root.is_dir():
        pytest.skip("Original diagnostic Parquet is not bundled in this source distribution")
    target = "locus_d4a6db6bdb2eacaf3c1bc28c81946f326e7f0d993d9f366debc10420129b3403"
    old = next(r for r in pq.read_table(root / "stop_loci").to_pylist() if r["locus_id"] == target)
    assert old["caller_support"] == 2
    assert old["geometry_eligible_primary"] is True  # Historical defect must remain observable.
    calls = pq.read_table(root / "cds_calls").to_pylist()
    defining = [
        r
        for r in calls
        if r["genome_id"] == old["genome_id"]
        and r["segments"] == [{"start0": 2930, "end0": 3917}]
        and r["rotation_shift"] == 0
    ]
    assert {r["accession_version"] for r in defining} == {"J02482.1"}
    assert {r["strand"] for r in defining} == {1}
    phano = next(r for r in defining if r["source"] == "phanotate")
    pyro = next(r for r in defining if r["source"] == "pyrodigal_gv")
    assert phano["rotation_sensitive"] is True and phano["eligible_primary"] is False
    assert pyro["rotation_sensitive"] is False and pyro["eligible_primary"] is True
    rebuilt = build_loci(calls, pq.read_table(root / "source_runs").to_pylist())
    locus = next(r for r in rebuilt["stop_loci"] if r["locus_id"] == target)
    assert locus["caller_support"] == locus["caller_support_raw"] == 2
    assert locus["caller_support_eligible"] == 1
    assert locus["caller_support_mask_raw"] == 3
    assert locus["caller_support_mask_eligible"] == 2
    assert locus["caller_families_eligible"] == ["pyrodigal_gv"]
    assert locus["support_assessable_primary"] is locus["eligible_primary"] is False
    assert locus["representative_call_id"] == pyro["call_id"]
    assert (
        "support_source_excluded:phanotate:rotation_sensitive" in locus["primary_exclusion_reason"]
    )
    links = [r for r in rebuilt["locus_models"] if r["locus_id"] == target]
    assert {identifier for link in links for identifier in link["call_ids"]} == {
        r["call_id"] for r in defining
    }


def test_ribosomal_slippage_holds_cds_not_genome():
    from phageannoscope.pipeline import apply_primary_code_policy

    base = {"source": "genbank", "transl_table_effective": 11, "eligible_primary": True}
    calls = [dict(base, exception="ribosomal_slippage", eligible_primary=False), dict(base)]
    runs = [{"source": "genbank"}]
    assert apply_primary_code_policy(calls, runs) is True
    assert calls[1]["eligible_primary"] is True
    calls = [dict(base, exception="ribosomal slippage"), dict(base)]
    assert apply_primary_code_policy(calls, [{"source": "genbank"}]) is False
    assert calls[1]["eligible_primary"] is False
