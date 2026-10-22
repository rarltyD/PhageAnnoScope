#!/usr/bin/env bash
# Run gene calling + matching and the PHROGs search for one or more prespecified batches.
# Resumable: finished steps are skipped (markers in logs/batches/); a stopped run restarts from
# per-genome checkpoints when the same command is repeated.
#
#   scripts/pipeline/run_batches.sh manifests/batch1.csv results/batch1 [manifests/batch2.csv results/batch2 ...]
#
# Environment (defaults suit a 4-thread machine):
#   CONFIG=configs/run_linux.json WORKERS=3 CORES=0-3 PHROGS_CPUS=3
#   workstation: CONFIG=configs/run_linux_workstation.json WORKERS=14 CORES=0-15 PHROGS_CPUS=14 PHAGEANNOSCOPE_MAX_CPUS=16
set -uo pipefail
[ $(( $# % 2 )) -eq 0 ] && [ $# -gt 0 ] || { echo "usage: $0 MANIFEST OUT [MANIFEST OUT ...]"; exit 2; }
CONFIG=${CONFIG:-configs/run_linux.json}; WORKERS=${WORKERS:-3}; CORES=${CORES:-0-3}; PHROGS_CPUS=${PHROGS_CPUS:-3}
HMM=${HMM:-data/db/phrogs_v4/phrogs_v4_m50.h3m}; PY=.venv/bin/python; L=logs/batches; mkdir -p $L
while [ $# -gt 0 ]; do
  man=$1; out=$2; shift 2; name=$(basename "$out")
  if [ ! -f $L/$name.callers.done ]; then
    echo "[$(date +%F' '%T)] callers $name"
    $PY scripts/pipeline/resource_guard.py --log $L/$name.guard.jsonl -- taskset -c $CORES nice -n 10 \
      .venv/bin/phageannoscope run --manifest "$man" --config $CONFIG --workers $WORKERS --out "$out" \
      --raw-root data/raw/callers > $L/$name.callers.out 2> $L/$name.callers.err
    # Exit code 1 only means that some genomes were held by the protocol (e.g. alternative code).
    grep -q '"caller_completed_n"' "$out/run_summary.json" 2>/dev/null && touch $L/$name.callers.done \
      || { echo "callers for $name did not finish; rerun to resume"; exit 3; }
  fi
  if [ ! -f $L/$name.phrogs.done ]; then
    echo "[$(date +%F' '%T)] PHROGs $name"
    $PY scripts/pipeline/resource_guard.py --log $L/$name.phrogs.guard.jsonl -- taskset -c $CORES nice -n 10 \
      $PY scripts/phrogs/search_proteins.py --run "$out" --hmm $HMM --cpus $PHROGS_CPUS > $L/$name.phrogs.log 2>&1 \
      && touch $L/$name.phrogs.done || { echo "PHROGs for $name failed"; exit 3; }
  fi
done
echo "[$(date +%F' '%T)] all batches done"
