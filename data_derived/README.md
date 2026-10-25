# Derived data

- cohorts/: predictions (one row per composition) and the threshold population, with entry dates, displayed labels and energies.
- dbmatch/: the earliest matching experimental structure record per prediction, overall and per source (COD, MP-ICSD, OQMD-ICSD).
- checks/: reconstruction (x0), power (x10), premise check, coverage of known compounds, and the fixed block edges.
- results/, results_for_claude.json, robustness_for_claude.json: every estimate in the paper and Supplementary Tables S2 and S3.
- manifest.json, plan_freeze.txt: SHA-256 checksums of all inputs and the plan freeze with its addendum.
- config_used.json: the configuration of the run.

Values derive from OQMD, the Materials Project and COD; cite those sources (see the paper) when reusing them.
