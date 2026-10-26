#!/usr/bin/env bash
# Copy the derived outputs of a run (work/) into data_derived/, the released derived data.
set -euo pipefail
cd "$(dirname "$0")/.."
test -f work/results/summary.json || { echo "run scripts/reproduce.sh first" >&2; exit 1; }
rm -rf data_derived/cohorts data_derived/matches data_derived/checks data_derived/results data_derived/figures
mkdir -p data_derived
cp -R work/cohorts work/matches work/checks work/results work/figures data_derived/
cp work/manifest.json data_derived/manifest.json
echo "data_derived/ updated from work/"
