# Derived data

Everything the analysis produced, from which every number in the paper can be checked.

| Path | Content |
|---|---|
| `cohorts/materials_project.csv`, `cohorts/oqmd.csv` | the predictions: one row per composition, with entry release and date, features and prototype |
| `cohorts/threshold__oqmd.csv` | the threshold population: displayed decomposition energy (`x_mev`), side, block and recomputed value |
| `matches/*.csv` | for each prediction, the earliest matching experimental record at each identity tier, overall and per source |
| `checks/reconstruction.json` | recomputed against displayed stability labels for each release |
| `checks/power.json` | window counts, minimum detectable effect and balance (check b) for each blocking scheme and window |
| `checks/premise.json` | cross-database label disagreement by distance to the threshold (check a) |
| `checks/coverage_known.json` | share of each release's experimental compositions found in the outcome record |
| `checks/threshold_blocks__oqmd.json` | the blocking scheme and its electronegativity quintile edges, fixed before outcomes |
| `results/summary.json` | every estimate in the paper: pre-existence and realization by cohort and release, capture-recapture, threshold design |
| `results/robustness.json` | every robustness and sensitivity analysis (Supplementary Table 3) |
| `results/figure_data.json`, `results/figure_checks.json`, `figures/` | the data behind Figs 1 to 3, their presentation checks, and the figures |
| `manifest.json` | SHA-256 of every input and the software versions of the run |
| `plan_freeze.txt` | checksums of the analysis plan and inputs, recorded before any prediction was matched to an outcome |

`tests/test_reference_results.py` recomputes every estimate from `cohorts/` and `matches/` and
checks it against `results/summary.json`.
