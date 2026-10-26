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
| `study_manifest.json` | the manifest of the original study run (October 2026, Python 3.11.2) |

`tests/test_reference_results.py` recomputes every estimate from `cohorts/` and `matches/` and
checks it against `results/summary.json`.

## Provenance

`plan_freeze.txt` and `study_manifest.json` come from the original study run. The plan was fixed on 25 October 2026 at 18:52 UTC, before any prediction was matched to an outcome. `manifest.json` and everything else in this folder come from the release 2.4.0 rerun from the inputs in `data/`, which reproduces every number in the paper.

Of the files whose checksums were frozen, the rerun reproduces `cohorts/oqmd.csv`, `cohorts/threshold__oqmd.csv` and `checks/threshold_blocks__oqmd.json` byte for byte; the experimental records match their frozen checksum in `data/CHECKSUMS.txt`. `cohorts/materials_project.csv` differs in one value: the decomposition energy of BeFeRh2 is -0.000371 eV/atom in the frozen file and -0.000372 eV/atom in the rerun, a last-digit floating-point difference between numerical-library builds. It changes no label and no estimate. `tests/test_released_data.py` checks these identities.

`study_manifest.json` is kept exactly as the original run wrote it. Its package list includes the Anthropic SDK, installed for an earlier design with a language-model labeling step that was dropped before the analysis plan was fixed; no language model is used in the analysis. Its file paths are those of the machine that ran the study.
