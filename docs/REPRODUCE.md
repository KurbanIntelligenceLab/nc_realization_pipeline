# Reproducing the study

All commands run from the repository root, in the environment set up as in the README.

## Level 1: recompute the estimates from the released derived data (seconds)

```bash
python -m pytest -q tests/test_reference_results.py tests/test_released_data.py
```

This does two things:
1. It copies `data_derived/cohorts/` and `data_derived/matches/` to a temporary directory, recomputes every estimate, and compares them with `data_derived/results/summary.json`.
2. It checks that `data_derived/` reproduces the headline numbers of the paper (`tests/paper_numbers.json`) and that the inputs in `data/` match the checksums recorded in the study.

## Level 2: run the whole analysis from the released inputs (4 to 30 minutes)

```bash
WORKERS=4 scripts/reproduce.sh
```

An independent rerun of release 2.4.15 took 227 s on an Apple M2 Max with `WORKERS=8`; machines with fewer cores take longer, up to about 30 minutes.

The script runs every step into a fresh `work/` directory and ends with `scripts/compare_results.py`, which compares the run with the paper:
- the counts must match exactly;
- the effect estimate must match to 1e-9.

| Step | Output in `work/` | What to expect |
|---|---|---|
| `check` | `manifest.json` | 279,210, 310,219 and 154,807 snapshot rows; 624,213 experimental records |
| `hull` | `hull/*.csv` | the longest step; parallel with `--workers` |
| `reconstruction` | `checks/reconstruction.json` | 2,860 and 2,282 displayed stable new compositions for OQMD v1.0 and v1.1, and 10,966 for the Materials Project |
| `cohorts` | `cohorts/*.csv` | 3,552 OQMD and 10,966 Materials Project predictions; 9,703 entries in the threshold population |
| `power` | `checks/power.json` | largest blocked standardized difference 0.029 at ±25 meV per atom; 1,408 and 3,738 entries in usable blocks |
| `premise` | `checks/premise.json` | label disagreement 26.9% (v1.0) and 29.4% (v1.1) within 25 meV per atom |
| `coverage-known` | `checks/coverage_known.json` | COD contains 52.7% of the OQMD v1.0 experimental compositions |
| `match` | `matches/*.csv` | 150 OQMD and 2,074 Materials Project predictions with a matching record |
| `estimate` | `results/summary.json` | the estimates in the paper |
| `robustness` | `results/robustness.json` | Supplementary Table 3 |
| `figures` | `figures/`, `results/figure_data.json` | Figures 1 to 3 and Supplementary Figures 1 and 2, drawn from the exported results only |
| exploratory scripts | `checks/novelty_checks.json`, `checks/literature_sensitivity.json`, `checks/audit_estimates.json` | the novelty checks, the literature-audit sensitivity (Supplementary Note 6) and the record-level audit estimates (Supplementary Note 7); identical to the files in `data_derived/checks/` |

To copy a run into the released derived data, run `scripts/export_derived.sh`.

## Level 3: rebuild the inputs from the original sources (several hours)

`scripts/rebuild_inputs.sh` rebuilds every file in `data/` from the public sources:
1. It downloads the OQMD releases v1.0, v1.1 and v1.8 and loads them into a local MySQL server, then exports and converts them.
2. It downloads the Materials Project build of 28 October 2022 from the AWS Open Data registry.
3. It exports the COD data table from COD's public server.
4. It merges the outcome records.

It needs a local MySQL server, the AWS command-line client and about 100 GB of free disk. The OQMD dumps are checked against `data/CHECKSUMS_RAW.txt`.

COD changes daily. An export made on a later day than the study's (25 October 2026) contains more records, so the matches can differ slightly from the paper's. All other inputs are fixed releases.

## Configuration

`config/study.json` holds every parameter of the analysis:
- the releases and their dates;
- the data cutoff (the COD export date, 25 October 2026);
- the horizons;
- the threshold windows and the blocking scheme;
- the premise-check pairs;
- the seed for the one-entry-per-chemical-system sensitivity analysis.

Paths are relative to the configuration file.
