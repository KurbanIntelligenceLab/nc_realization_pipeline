# Changelog

## 2.4.14 (29 October 2026)

Figure revisions for submission. No data, estimate or Source Data value changes. Every plotted value was verified identical to 2.4.13, and all ten presentation checks pass.
- **Fig. 1:** the outcome text becomes an aligned table with each release's shares. The after-release axis is labeled "(expanded scale)". Records for which the literature audit found an earlier experimental report are drawn as open circles: `render_all` and `fig1` take the audit, which the figure step reads beside the workdir, as for the Source Data. The before-release points are slightly smaller and lighter.
- **Fig. 2:** all three panels share one row per release. Panel c shows horizons as offset points, with a horizon and count column. Release labels carry the cohort size once. Panel a gives the combined "Either source" coverage the dominant mark.
- **Fig. 3:**
  - Panel b draws the pooled window disagreement, derived exactly from the bins with a Wilson interval, and limits the criterion line to the window it applies to.
  - Panel d groups the windows and the placebo thresholds under headings that state the unit once, marks "(no events)", and draws the detectable-effect lines on the primary row only.
- **Supplementary Fig. 1:** follow-up is drawn as in Fig. 1, horizons are labeled, release labels match Figs 1 and 2, and a bracket marks the period in which the record rests on COD alone.
- **Supplementary Fig. 2:**
  - The primary window comes first, as in Fig. 3d.
  - The axis is tightened to the data while keeping the full limit band; its lower limit uses signed values.
  - Each panel states the largest remaining difference; the legend marks the first-specified scheme and says "Prespecified".
- **Tests:** `tests/test_figure_audit_marks.py`.

## 2.4.13 (27 October 2026)

Documentation only: README and CITATION carry the paper's final title, "Hypothetical labels in computational materials databases are not a reliable proxy for compositional novelty". No code, data or estimate changes.

## 2.4.12 (27 October 2026)

The record-level audit of prior-record matches, and the paper's revised title. No prespecified result changes.
- **`data/audit/prior_record_audit.csv`:** the 60 audited matches (stratified, seed 20261026), with the composition, year and experimental status each original paper reports, the evidence tier, and the auditor's notes.
- **`scripts/audit_estimates.py`** (output `data_derived/checks/audit_estimates.json`): recomputes each verdict under the paper's composition rule with the pipeline's own distance, weights confirmation by stratum, and gives Jeffreys-based 95% intervals by simulation.
  - Confirmed: 53 of 60 under the rule, 49 exactly.
  - Matched entries confirmed: an estimated 1,880 of 2,203.
- **README and CITATION:** the title is "More than a thousand inorganic compositions labeled hypothetical and stable already had recorded experimental structures", because the audit-adjusted count is below two thousand.

## 2.4.11 (26 October 2026)

Reviewer-driven corrections. No prespecified result changes.
- **`scripts/novelty_checks.py`:** a third exploratory check, pre-existence restricted to exact composition matches (tolerance 1e-6), using the same matching code: 3.4%, 5.6% and 16.2% (98 of 2,860; 39 of 692; 1,779 of 10,966), against 3.6%, 6.2% and 18.8% under the paper's rule. The tolerance counts it recomputes reproduce the paper's exactly. The output in `data_derived/checks/novelty_checks.json` now includes it.
- **Figure 3 and Source Data Fig. 3:** the detectable-effect line is described as an approximate normal-theory 80%-power calculation; "reliably detectable" is removed. The legend reads "Detectable at 80% power, approx."; all ten figure presentation checks pass, including the one that rejects text past the panel edge.
- **README:** states the power calculation as approximate, with the interval extending to 0.38; states the exploratory 0.63% bound as a bound on the stable-side recording probability, not on the effect; and adds the exact-composition rates.

## 2.4.10 (26 October 2026)

Wording only. The audit and README now describe post-release cases as the "earliest located report" rather than a "first report", because a literature search cannot establish that no earlier report exists. The unverified description of the Li2AgGe paper as a single-crystal determination is removed. No classification or result changes.

## 2.4.9 (26 October 2026)

Audit entries resolved. No database-defined result changes, and the threshold sensitivity analysis is unchanged.
- **MgSc2Ru:** confirmed as a first report after release. It is the 2015 crystal-structure report of Sc2RuMg (Kersting, Rodewald and Poettgen, Z. Kristallogr. 230, 151-155, 2015), tetragonal rather than the predicted Heusler structure.
- **InCl3:** already known under the outcome definition, from a 1954 powder-diffraction structure assignment (Templeton and Carter, J. Phys. Chem. 58, 940-944, 1954), which replaces the secondary citation.
- **README:** the ten-year realizations are now stated under both definitions.

## 2.4.8 (26 October 2026)

The literature audit of individual outcomes, kept separate from the database-defined results, which are unchanged.
- **`data/audit/literature_audit.csv`:** the four OQMD v1.0 compositions recorded within ten years and the four events of the primary threshold window, with what each source establishes and the classification under the outcome definition and under a broader one that counts any earlier experimental report.
- **Source Data:** `ncr/source_data.py` marks the audited rows of Source Data Fig. 1 and adds a Literature_audit sheet to Figs 1 and 3, plus a note to Fig. 2. Counts and formulas are unchanged. The audit is read from `data/audit/` beside the run's workdir, so synthetic worlds and the frozen study configuration are unaffected.
- **`scripts/literature_sensitivity.py`** (output `data_derived/checks/literature_sensitivity.json`): reruns the blocked estimator and the exact test with the events recoded under the broader definition. It refuses to run if any threshold event lacks an audit entry.
- **Tests:** `tests/test_literature_audit.py`.

## 2.4.7 (26 October 2026)

Documentation only. The README reports the literature check of the four ten-year OQMD v1.0 realizations: two (EuNbO3, InCl3) had been reported before the release in sources outside the open record. The threshold summary now describes the comparison as too imprecise to determine whether the label changed later structure reporting.

## 2.4.6 (26 October 2026)

Documentation only. The README states the bounds as the paper does: prior-record counts are lower bounds, while realization counts can err in either direction. It also describes the space-group check as a diagnostic, not a structural match.

## 2.4.5 (26 October 2026)

Two exploratory checks and the paper's final title; no estimate changes.
- **`scripts/novelty_checks.py`** (output `data_derived/checks/novelty_checks.json`): for the predictions already known at release, whether a matching COD record dated by the release has the predicted space group (25.0%, 13.0% and 27.4% for OQMD v1.0, OQMD v1.1 and the Materials Project), and the share of already-known Materials Project compositions whose earliest record predates the project's start in 2011 (92.4%). Both were added after the outcomes were measured and are reported as exploratory.
- **README and CITATION:** the title is "Thousands of inorganic compositions labeled hypothetical and stable already had recorded experimental structures", and the results table names the cohorts accordingly.

## 2.4.4 (26 October 2026)

Documentation only: `README.md` and `CITATION.cff` carry the paper's revised title, "Thousands of inorganic compounds released as new and stable by computational databases were already known". No code, data or estimate changes.

## 2.4.3 (26 October 2026)

Figure 1 and its Source Data only; no estimate changes.
- **Accurate outcome labels:** Figure 1 now labels its counts "recorded within N yr" and "not recorded within N yr", where N is the release's longest horizon of complete follow-up. The third count had read "no record", but six Materials Project predictions were first recorded in 2025 or 2026, after their two-year follow-up, so they do have a record. No analysis counts them, and no estimate changes.
- **Every point visible:** points that share a date (year-only dates placed on 1 July) are now spaced evenly across the band without random offsets, so the nine Materials Project predictions first recorded in 2023 no longer overlap.
- **Source Data:** the Figure 1 outcome columns are renamed to match, and the README explains the third count.

## 2.4.2 (26 October 2026)

Documentation only.
- `data_derived/README.md` explains the packages and paths in the original study manifest, which is kept unchanged.
- `README.md` notes what release 2.2.2 in the history contains.

## 2.4.1 (26 October 2026)

Documentation and provenance only; no code path or estimate changes.
- **Provenance restored:** `data_derived/plan_freeze.txt` (the plan frozen before matching) and the original study's manifest (`data_derived/study_manifest.json`), which the 2.4.0 update had replaced.
- **Frozen-file identity test:** the released OQMD cohort, threshold cohort and blocking scheme equal their frozen checksums byte for byte. The released experimental records equal theirs. The one-value difference in the Materials Project cohort is documented.
- **Exact environment:** `requirements-lock.txt` lists every package of the environment that reproduced the paper.
- **README:** the threshold result is stated as in the paper (detectable difference 0.24 percentage points; the 0.6 bound is exploratory).

## 2.4.0 (26 October 2026)

Changes after the first public release. Every estimate in the paper is unchanged; the coverage fix changes what a rerun of the coverage check computes, so that it reproduces the published 52.7%.
- **Fix: the known-compound coverage check now uses COD records only,** as the paper defines it (COD share of each release's experimental compositions). It had read the merged record, whose ICSD-derived OQMD source contains most of these compositions, so a rerun would not have reproduced the reported 52.7%. It now also reports coverage by chemistry class (Supplementary Table 2c), and `compare_results.py` checks the coverage counts.
- **Fig. 1 redesigned:** a single panel on a split time axis aligned at each release (decades before, years after), with median-age markers and the outcome of every prediction; the hull schematic moves to Fig. 3a.
- **Fig. 2 redesigned:** (a) capture-recapture coverage of the record with intervals, (b) already known by source, (c) recorded later at every complete horizon with simultaneous intervals and counts; a check reads every plotted point back.
- **Fig. 3 linked across panels:** the near-hull window is named in the schematic (a) and shaded in the premise panel (b) as in (c); the criterion line is labeled; panel d marks the smallest reliably detectable effect (exported as `fig3.mde`).
- **Source Data workbooks (Nature Portfolio requirement):** `python -m ncr figures` now also writes one Excel workbook per figure to `work/figures/source_data/`, one sheet per panel. Every plotted record is included. Shares, Wilson and Bonferroni-adjusted intervals, coverage bounds and counts are Excel formulas, so they can be audited in Excel. Empty counts leave blank cells, never #DIV/0!. `openpyxl==3.1.5` is added to the requirements, and a new test checks the workbooks on the synthetic world.
- **Figure lint in the pipeline:** a new figure check fails on overlapping or clipped text, hyphens used as minus signs, and axes whose tick labels mix decimal places; it runs on the study data and on the synthetic test world. Ticks are chosen from the data range at nice steps with consistent decimals. Fig. 3c's two histogram sides now differ in greyscale (contrast ratio 1.09 -> about 1.6). Fig. 1 no longer clips the oldest record.
- **Cross-figure consistency:** release names with years in every figure (as in Fig. 1); Fig. 2c names each release once and lists horizons; Fig. 2a uses the source swatches of Fig. 2b, so a filled black circle means "recorded later" throughout; Fig. 3 gives units on every energy label, writes "percentage points", and names the unstable side the same way in a and c; unused palette removed.
- **Supplementary figures:** Supp. Fig. 1 keys every element, marks that the outcome sources continue before 2010, and shows the 1-yr reporting-lag buffer; Supp. Fig. 2 shows all three windows (10, 25, 50 meV/atom), as the text cites, against the labeled pre-specified limit, with true minus signs. New check: the figure draws every exported window, read back against the export.

## 2.3.0 (25 October 2026)

First public release.
- **Package:** reorganized into one package (`python -m ncr`) with descriptive command names.
- **Data:** the converted inputs (`data/`) and all derived data (`data_derived/`) are included.
- **Reproduction:** scripts reproduce the study at three levels (`docs/REPRODUCE.md`).
- **Results unchanged:** every estimate is identical to those of the study run (613 values compared).
- **Figure 2 redesigned:** three panels: already known at release; already-known predictions by record source; recorded later.
- **Figure 3 redesigned:** four panels: label disagreement between databases, balance by blocking scheme, outcomes by bin with entries per bin, and the label effect under every sensitivity analysis. Figures are drawn by a separate `figures` step after `robustness`.
- **All figures restyled:** one restrained house style (Arial-type 6–6.5 pt, 0.5 pt rules, muted colorblind-safe palette, no in-plot annotations). Figure 1 is drawn from data: a convex-hull schematic, the dated releases against the years each outcome source covers, and the outcome classes. Each color keeps one meaning across figures (release, record source, side of the threshold); Fig. 3c shows each event as a dot above the entries per bin; the balance and effect panels are monochrome. Figures follow Nature's published specifications (183 mm width, 5-7 pt text, 8 pt bold panel letters, strokes of 0.25-1 pt, black text only, colors keyed inside the figure), and each is checked in code.
- **Figures simplified:** each main figure carries one message in at most three panels: Fig. 1 (labeling and outcome classes), Fig. 2 (one square per OQMD v1.0 prediction; already known and recorded later by release), Fig. 3 (premise check, the two sides of the threshold, the effect with its windows and placebos). The release timeline, the balance plot and the events-by-bin plot move to Supplementary Figs 1 to 3. Fig. 1b now shows every matched prediction at the date of its earliest record, Fig. 2 shows the sources of already-known predictions and cumulative time-to-record curves, and Fig. 3b shows every event above the entries per bin. Checks count every plotted point back against the estimates. The ledger export asserts that its counts equal the summary.

## 2.2.0

- Added the `robustness` command for every sensitivity analysis in the paper.

## 2.1.2

- The Wilson interval is now exactly 0 (or 1) at zero (or all) events; floating-point rounding had left about 1e-18.
- Figure intervals: Wilson intervals per bin, and Bonferroni intervals across horizons in Fig. 2b.

None of these changes alters an estimate.

## 2.1.1

- The exact binomial density test is computed in log space; it overflowed for more than about 1,030 entries.

## 2.1.0

- Analysis plan fixed before any prediction was matched to an outcome (`data_derived/plan_freeze.txt`).
