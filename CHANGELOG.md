# Changelog

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
