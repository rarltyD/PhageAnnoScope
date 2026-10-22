from phageannoscope.functional import _union_length, classify

ANN = {
    "phrog_1": {"annot": "integrase", "category": "integration and excision"},
    "phrog_2": {"annot": "terminase", "category": "head and packaging"},
    "phrog_3": {"annot": None, "category": "unknown function"},
}


def hit(phrog, e=1e-20, bits=100.0, q=0.9, p=0.9):
    return {"phrog": phrog, "evalue": e, "bitscore": bits, "query_coverage": q,
            "profile_coverage": p}


def test_union_length_merges_overlaps_inclusive():
    assert _union_length([(1, 10), (5, 20), (30, 30)]) == 21


def test_no_hits_is_no_informative_hit():
    assert classify([], ANN)["function_status"] == "no_informative_hit"


def test_thresholds_are_inclusive_and_all_required():
    assert classify([hit("phrog_1", e=1e-5, q=0.5, p=0.5)], ANN)["function_status"] == (
        "specific_function")
    for bad in (hit("phrog_1", e=2e-5), hit("phrog_1", q=0.49), hit("phrog_1", p=0.49)):
        assert classify([bad], ANN)["function_status"] == "no_informative_hit"


def test_unknown_best_hit_is_homology_unknown():
    row = classify([hit("phrog_3", bits=200), hit("phrog_1", bits=50)], ANN)
    assert row["function_status"] == "homology_unknown_function"
    assert row["best_phrog"] == "phrog_3"


def test_near_tie_with_conflicting_known_categories_is_ambiguous():
    row = classify([hit("phrog_1", bits=100), hit("phrog_2", bits=95)], ANN)
    assert row["function_status"] == "ambiguous"
    row = classify([hit("phrog_1", bits=100), hit("phrog_2", bits=80)], ANN)
    assert row["function_status"] == "specific_function"


def test_unannotated_phrog_id_defaults_to_unknown():
    assert classify([hit("phrog_999")], ANN)["function_status"] == "homology_unknown_function"
