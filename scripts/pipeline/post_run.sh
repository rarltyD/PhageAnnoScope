#!/usr/bin/env bash
# Statistics and figures for a finished (or merged) run.
#   scripts/pipeline/post_run.sh results/full_cohort figures/
set -euo pipefail
RUN=${1:?run directory}; FIG=${2:-figures}
MAN=${MANIFEST:-manifests/cohort_manifest_validated.csv.gz}
PY="taskset -c ${CORES:-0-3} nice -n 10 .venv/bin/python"
$PY scripts/analysis/analyze_run.py --run "$RUN" --manifest "$MAN" --reps 2000 ${REUSE:+--reuse-locus-table}
$PY scripts/analysis/make_figures.py --run "$RUN" --out "$FIG" --reps 1000
$PY scripts/analysis/extra_stats.py --run "$RUN" --manifest "$MAN"
(cd scripts/analysis && $PY extra_figures.py --run "../../$RUN" --out "../../$FIG")
$PY scripts/analysis/robustness_stats.py --run "$RUN" --manifest "$MAN" --fig-out "$FIG"
$PY scripts/analysis/genbank_function_breakdown.py --run "$RUN" --out "$FIG"
echo "post_run done: $RUN"
