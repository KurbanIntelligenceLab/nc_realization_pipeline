#!/usr/bin/env bash
# Run the whole analysis from the released inputs in data/, including the exploratory analyses,
# and compare the results with the numbers reported in the paper. About 4 to 30 minutes on a
# laptop (227 s on an Apple M2 Max with 8 workers); most of it is the hull.
#
#   scripts/reproduce.sh                 # uses python3 and 4 worker processes
#   PYTHON=.venv/bin/python3 WORKERS=8 scripts/reproduce.sh
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
CFG="${CONFIG:-config/study.json}"
WORKERS="${WORKERS:-4}"

if [ -d work ] && [ -n "$(ls -A work 2>/dev/null)" ]; then
  echo "work/ is not empty; move or delete it first so that no earlier output is reused" >&2
  exit 1
fi
run() { echo "== $*"; "$PY" -m ncr --config "$CFG" "$@"; }
run check
run hull --workers "$WORKERS"
run reconstruction
run cohorts
run power
run premise
run coverage-known
run match
run estimate
run robustness
run figures
# the exploratory analyses (Supplementary Notes 6 and 7 and the novelty checks), on this run
echo "== exploratory analyses"
"$PY" scripts/novelty_checks.py --derived work
"$PY" scripts/literature_sensitivity.py --derived work
"$PY" scripts/audit_estimates.py --derived work
"$PY" scripts/compare_results.py work
