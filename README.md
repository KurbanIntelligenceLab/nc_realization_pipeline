# nc_realization_pipeline

Code and data for **"Hypothetical labels in computational materials databases are not a reliable proxy for compositional novelty"** by Hasan Kurban and Mustafa Kurban.

The study follows every composition that three dated public database releases displayed as stable and new:
- OQMD v1.0 (April 2014);
- OQMD v1.1 (October 2016);
- the Materials Project build of 28 October 2022.

Each prediction is matched deterministically against dated records of experimentally determined crystal structures from three sources:
- the Crystallography Open Database (COD);
- the ICSD-derived records of the Materials Project;
- the ICSD-derived records of OQMD.

Near the stability threshold, compounds just inside and just outside the stable label are compared to estimate the effect of the label itself.

## Main results

Prior-record counts are lower bounds. Realization counts can err in either direction, because a missed earlier report can make a known composition appear realized.

| | OQMD v1.0 (2014) | OQMD v1.1 (2016) | Materials Project (2022) |
|---|---|---|---|
| Compositions labeled hypothetical and stable | 2,860 | 692 | 10,966 |
| Already had a recorded structure | 103 (3.6%) | 43 (6.2%) | 2,057 (18.8%) |
| Acquired one later | 4 within 10 years (0.14%) | 0 within 5 years | 11 within 2 years (0.10%) |

These counts are database-matched: each date is that of the earliest matching record in the analyzed sources. A literature audit (`data/audit/literature_audit.csv`; paper, Supplementary Note 6) found earlier experimental reports for two of the four OQMD v1.0 compositions recorded within ten years: EuNbO3 had been synthesized by 1998, and InCl3 had an experimentally characterized structure by 1954. After incorporating the earlier structural report for InCl3, three of the four retain post-release dates for their earliest located structure reports; under the broader criterion of any experimental report, two remain. Under the outcome definition, an experimentally determined crystal structure, the audit reclassifies none of the four threshold events. `scripts/literature_sensitivity.py` reproduces the sensitivity analysis that counts any earlier experimental report (0.10 percentage points, 95% CI -0.10 to 0.30).

**Threshold window.** Within ±25 meV/atom of the OQMD stability threshold, 3 of 1,418 compounds labeled stable and 1 of 3,783 labeled unstable acquired a recorded structure within five years. The comparison is too imprecise to determine whether the label changed later structure reporting: an approximate normal-theory power calculation puts the smallest difference detectable with 80% power at 0.24 percentage points, which describes the design rather than the effects the data exclude (the 95% interval extends to 0.38). An exploratory bound, added after the outcomes were seen, puts the stable-side recording probability below 0.63%.

**Record-level audit** (`data/audit/prior_record_audit.csv`; `scripts/audit_estimates.py`, output in `data_derived/checks/audit_estimates.json`). One author checked a stratified random sample of 60 prior-record matches against the original papers. Under the paper's composition rule, 53 were confirmed (49 with exact composition); 4 records were not experimental structures (3 first-principles structures imported through OQMD's ICSD library, 1 with no experimental source found), and 3 sources report a different composition. Weighted by stratum, an estimated 1,880 (95% interval 1,380-2,060) of the 2,203 matched entries are confirmed, and Materials Project pre-existence falls from 18.8% to 15.9% (11.4-17.5%).

**Exploratory checks on the already-known compositions** (`scripts/novelty_checks.py`, output in `data_derived/checks/novelty_checks.json`). Matching is at the level of composition. Among already-known predictions with a COD record dated by the release, 25.0% (OQMD v1.0), 13.0% (OQMD v1.1) and 27.4% (Materials Project) had such a record in the predicted space group (a diagnostic, not a structural match); restricted to exact composition matches, pre-existence is 3.4%, 5.6% and 16.2% rather than 3.6%, 6.2% and 18.8%; the others have recorded structures of the same composition in another space group. For 1,900 of the 2,057 already-known Materials Project compositions (92.4%), the earliest record predates 2011, when the project started, so these cannot be later realizations of its predictions.

## Quick start

Requires Python 3.11.

```bash
git clone https://github.com/KurbanIntelligenceLab/nc_realization_pipeline.git
cd nc_realization_pipeline
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt -r requirements-dev.txt
python -m pytest -q
```

`requirements.txt` pins the direct dependencies. `requirements-lock.txt` lists every package of the environment that reproduced the paper with release 2.4.0 (Python 3.11.17, macOS arm64); install it instead to match that environment exactly.

The tests run in about a minute. They cover three things:
- an end-to-end run on a synthetic world with known truth;
- the converters and the statistics;
- the published results: every estimate is recomputed from `data_derived/`, and the inputs in `data/` are checked against the checksums recorded in the study.

## Reproducing the study

There are three levels, each documented in [`docs/REPRODUCE.md`](docs/REPRODUCE.md).

| Level | What it does | Time | Command |
|---|---|---|---|
| 1 | Recompute every estimate from the released cohorts and matches | seconds | `python -m pytest -q tests/test_reference_results.py` |
| 2 | Run the whole analysis from the released inputs and compare with the paper | 15 to 30 minutes | `scripts/reproduce.sh` |
| 3 | Rebuild the inputs from the original public sources | several hours, about 100 GB | `scripts/rebuild_inputs.sh` |

## Repository layout

| Path | Content |
|---|---|
| `ncr/` | the pipeline: converters, convex hulls, cohorts, matching, checks, estimates, robustness and figures |
| `config/study.json` | every parameter of the study |
| `data/` | converted inputs (compressed) with their checksums; see [`data/README.md`](data/README.md) |
| `data_derived/` | all derived data and results; see [`data_derived/README.md`](data_derived/README.md) |
| `scripts/` | reproduction, comparison and rebuild scripts |
| `tests/` | the test suite and the headline numbers of the paper (`tests/paper_numbers.json`) |

## Pipeline commands

All commands are run as `python -m ncr --config config/study.json <command>`, in this order:

| Command | Purpose |
|---|---|
| `check` | record the checksum of every input |
| `hull` | recompute the decomposition energy of every hypothetical entry in each release |
| `reconstruction` | compare recomputed and displayed stability labels |
| `cohorts` | build the predictions and the threshold population |
| `power`, `premise`, `coverage-known` | checks that use no outcome |
| `match` | match every prediction against the experimental records |
| `estimate`, `robustness` | all estimates and sensitivity analyses |
| `figures` | Figures 1 to 3 and Supplementary Figures 1 and 2 from the exported results |

The converters (`convert-oqmd`, `convert-mp`, `export-cod`, `convert-cod`, `convert-oqmd-icsd`, `convert-mp-icsd`, `merge-experiments`) build the inputs from the original sources (level 3).

## Repository history

Commit `43f659e` (release 2.2.2) predates the final design. It contains a language-model labeling step and a literature-search step that were dropped before the analysis plan was fixed, together with working files from the AI-assisted development disclosed in the paper. Releases 2.4.0 and later hold the code and data for the paper.

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite the data sources: OQMD, the Materials Project and COD.

## License

The code is released under the MIT License ([`LICENSE`](LICENSE)). The data in `data/` and `data_derived/` derive from OQMD and the Materials Project (CC BY 4.0) and from COD; reuse them under the terms of those sources, with attribution.
