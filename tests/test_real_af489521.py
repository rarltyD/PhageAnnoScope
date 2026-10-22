"""Real-data regression: AF489521 first reverse-strand CDS shares its stop.

GenBank complement(210..752); Pyrodigal-gv 210..833; PHANOTATE 210..773 (earlier
smoke test).  On the minus strand the stop codon is at 210..212, so every pair
must be classified ``alternative_start`` -- never exact or partial overlap.
"""

import os
from pathlib import Path

import pytest

from phageannoscope.config import load_config
from phageannoscope.pipeline import process_genome
from phageannoscope.util import file_sha

GBK = Path(__file__).resolve().parents[1] / "data/raw/ncbi/pilot_100/AF489521.gbk"
CONFIG = os.environ.get("PHAGEANNOSCOPE_TEST_CONFIG")


@pytest.mark.integration
@pytest.mark.skipif(not GBK.exists() or not CONFIG, reason="needs local data and caller config")
def test_af489521_shared_stop_alternative_starts(tmp_path):
    row = {"accession_version": "AF489521.1", "genbank_path": str(GBK),
           "genbank_sha256": file_sha(GBK), "eligible_for_pilot": "true"}
    result = process_genome(row, load_config(CONFIG), tmp_path)
    runs = {r["source"]: r for r in result["source_runs"] if r.get("rotation_shift", 0) == 0}
    assert {s: runs[s]["status"] for s in runs} == {
        "genbank": "success", "phanotate": "success", "pyrodigal_gv": "success"}
    calls = {c["call_id"]: c for c in result["cds_calls"]}

    def first_reverse(source):
        rows = [c for c in calls.values() if c["source"] == source and c["strand"] == -1
                and c.get("rotation_shift", 0) == 0]
        return min(rows, key=lambda c: c["segments"][0]["start0"])

    first = {s: first_reverse(s) for s in ("genbank", "phanotate", "pyrodigal_gv")}
    assert {s: c["segments"][0]["start0"] for s, c in first.items()} == {
        "genbank": 209, "phanotate": 209, "pyrodigal_gv": 209}
    ids = {c["call_id"] for c in first.values()}
    pairs = [m for m in result["matches"] if m["call_id_a"] in ids and m["call_id_b"] in ids]
    assert len(pairs) == 3
    assert {m["relation"] for m in pairs} == {"alternative_start"}
    assert all(m["relative_frame"] == 0 for m in pairs)
