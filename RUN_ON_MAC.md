# Running the study pipeline on a Mac

Everything runs locally in Terminal. The pipeline downloads nothing except two optional things you switch on yourself:
- literature searches through OpenAlex;
- language-model calls through the Anthropic API.

Every step writes files under your work folder, and you can stop and resume at any point.

---

## 1. One-time setup (about 10 minutes)

Open **Terminal** (Applications → Utilities → Terminal) and check your Python version:

    python3 --version

You need **Python 3.12**. If you have something else, install 3.12 from python.org (macOS installer) or with Homebrew (`brew install python@3.12`), then use `python3.12` in place of `python3` in the next command.

Unzip the bundle and create an isolated environment:

    cd ~/Downloads
    unzip nc_realization_pipeline.zip
    cd nc_realization_pipeline
    python3 -m venv .venv
    source .venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements-pipeline.txt

Run the end-to-end test on synthetic data. It uses a mock language model and costs nothing:

    python3 tests/test_pipeline.py

It takes one to two minutes and must end with `ALL PIPELINE CHECKS PASSED` (25 checks). If it does not, send me the full output and stop.

**Each time you open a new Terminal window**, reactivate the environment:

    cd ~/Downloads/nc_realization_pipeline && source .venv/bin/activate

---

## 2. Set up your study folder

    python3 ncr_cli.py init --dir ~/ncr_study

Then:
1. Edit `~/ncr_study/config.json` in any text editor.
2. Put your data files where the config points.
3. Read `DATA_FORMATS.md` for every column.

Every `SET_ME` must be filled before the step that needs it, and the pipeline stops with a clear message if one is missing.

In every command below, add `--config ~/ncr_study/config.json`. To save typing, define a shortcut once per Terminal window:

    alias ncr="python3 $PWD/ncr_cli.py --config ~/ncr_study/config.json"

**Cutoff date.** `cutoff_date` is the date of your COD export, the latest outcome record that can have been observed. It may not lie in the future.

---

## 3. Milestone 1: power and premise checks

This milestone needs no language model and no experts, so it costs nothing. **Do this first.** It decides whether contribution 2 is viable (REDFLAGS C5 and M2).

    ncr check                         # validates inputs, records SHA-256 of every file
    ncr hull --all --workers 8        # recompute decomposition energies (8 = your CPU cores; see: sysctl -n hw.ncpu)
    ncr x0                            # snapshot reconstruction: OQMD 2015 must give 3,231 and 2,290
    ncr cohorts                       # cohorts, threshold populations, phase-1 samples
    ncr x10                           # window census and minimum detectable effect
    ncr premise                       # cross-database label disagreement near the hull

**Run time.** `hull` is the slow step: it can take hours for a full database. For long runs, keep the Mac awake for as long as this Terminal window is open, then run the steps as usual:

    caffeinate -i -w $$ &

(`caffeinate` starts its command directly, without the shell, so it cannot run the `ncr` shortcut itself.)

If it stops, run it again; finished releases are skipped.

**Send me these four files:**
- `~/ncr_study/work/checks/x0.json`
- `~/ncr_study/work/checks/x10.json`
- `~/ncr_study/work/checks/premise.json`
- `~/ncr_study/work/manifest.json`

They contain only counts, hashes and rates.

---

## 3b. OQMD from MySQL dumps

See the conversation log for loading the dumps. Once `qmdb_v1_0` and `qmdb_v1_1` are loaded and exported to `.tsv` files, convert each one with `ncr convert-oqmd` (`DATA_FORMATS.md`, section 1).

**Displayed labels are strict.** An entry whose release shows no stability value is not treated as presented stable. It is excluded from cohorts and from the threshold population, and the counts are logged.

**One entry per composition.** For composition-level sources, each composition is represented by its lowest-energy entry in each release.

## 4. Milestone 2: the OQMD pilot (E1)

This milestone uses the language model and your experts. Replace `oqmd_2015` with the cohort name that `ncr list` shows.

### Before extraction

1. **Register the analysis plan** (REDFLAGS M1). Registration must come before any outcome extraction.
2. **Fill in the language-model settings** in the config: `llm.model`, plus `retrieval.backend` and `retrieval.openalex_email`.
3. **Set your API key** in the Terminal window. It is never stored by the pipeline:

       export ANTHROPIC_API_KEY="sk-ant-..."

### Run the pilot

    ncr list                                   # shows cohort names and progress
    ncr dbmatch  --cohort oqmd_2015            # deterministic matches against your ICSD/COD export
    ncr retrieve --cohort oqmd_2015            # literature search (cached)
    ncr plan     --cohort oqmd_2015            # prints the number of model calls BEFORE you spend anything
    ncr extract  --cohort oqmd_2015            # language-model extraction (cached, resumable)
    ncr label    --cohort oqmd_2015
    ncr freeze   --cohort oqmd_2015            # hashes and locks the labels; required before sampling
    ncr sample   --cohort oqmd_2015            # expert sample with recorded seed

### Experts

The sheets are in `~/ncr_study/work/experts/oqmd_2015/`:
- Send `expert_A.csv` to one expert and `expert_B.csv` to the other.
- **Do not send `sample_key.csv`.**
- Follow the protocol in the runbook (Appendix C) and in `DATA_FORMATS.md`, section 6.
- Put the completed sheets back in the same folder with the same names.

### After the experts return their sheets

    ncr adjudicate --cohort oqmd_2015              # kappa, and a list of disagreements
    # the third expert fills adjudication.csv in the same folder, then:
    ncr adjudicate --cohort oqmd_2015 --finalize
    ncr audit-sample --cohort oqmd_2015            # X3: 50 intermetallic + 50 other negatives for the third expert
    # the third expert fills residual_audit_sheet.csv, then:
    ncr estimate

**Send me** `~/ncr_study/work/results_for_claude.json`. It holds aggregates only, with no formulas, identifiers, documents or ICSD content. I will:
- apply the pilot's go/adjust/stop gate;
- set the expert allocation for the other cohorts from the pilot's error rates;
- reconcile the manuscript.

---

## 5. Milestone 3: the remaining cohorts and checks

Repeat the Milestone 2 sequence for each remaining cohort (`ncr list` shows them):
- the other base cohorts;
- the `threshold__...` populations (E6 and the X9 attention check);
- `relabel__...` after running `ncr relabel-pop` (E7);
- `..._phase1` for large cohorts.

Then run the two validation checks:

    ncr control-cohort                                    # X1 control panel
    ncr retrieve --cohort control_panel && ncr extract --cohort control_panel && ncr label --cohort control_panel
    ncr rerun-cohort --cohort oqmd_2015 --n 300           # X5 rerun
    ncr extract --cohort oqmd_2015 --rerun && ncr label --cohort oqmd_2015 --rerun

Finally, run `ncr estimate` and send me `results_for_claude.json` again. Figures are written to `work/figures/`.

---

## 6. What never leaves your Mac

- Your ICSD export, and anything derived row by row from it.
- The API key.
- The `sample_key.csv` files.

Expert sheets stay with your team.

`results_for_claude.json` is designed to be shareable. Open it before sending if you want to confirm what it contains.

---

## 7. Implementation choices for the authors to confirm

These were not fixed by the manuscript, so the code makes a documented choice in each case.

1. **Month-only dates** (`YYYY-MM`) are set to the 15th. Year-only dates use 1 July, as in Table S1.
2. **Literature structure proxy.** Without a refined structure, a literature record counts as structure-level (T2) on the model side when its reported space group or prototype equals the prediction's. This affects only the model labels, which the expert sample corrects; experts still decide T1 and T2.
3. **When adjudication is needed.** Two experts' records go to adjudication only when they imply different outcomes: different pre-existence, or a report falling in a different horizon window. Records that differ only within the same window are resolved by taking the earlier date.
4. **Premise check (E6a).** Compounds are matched across databases by reduced formula, optionally with space group. The lowest-energy entry per composition represents each database.
5. **Balance covariates (check b).** These are the number of elements, the mean electronegativity, the largest atomic number and indicators for each chemistry class.
6. **Generation settings.** The current Anthropic Python SDK (2.3.0) has no temperature argument, so the pipeline sends none unless you set `llm.temperature` (passed as an extra field). The manuscript's "temperature zero" was changed to "fixed, recorded generation settings". Run-to-run variation is measured by X5, and the frozen outputs make every estimate reproducible.
7. **Two-phase cohorts.** Ratio intervals use a delta-method phase-1 term. Simulated coverage was 0.950 in the test.
8. **Relabeling scheme application** is smoke-tested only. It needs real Materials Project entries with calculation parameters. Check the group counts it prints against your expectations before you trust them.

---

## 8. Troubleshooting

- **`pip install` fails on pymatgen.** Install Apple's command-line tools (`xcode-select --install`) and retry. Make sure you are on Python 3.12.
- **`config value ... is still SET_ME`.** Fill that field in `config.json`.
- **`labels ... changed after freezing` or `is frozen`.** This is intended: labels cannot change once the expert sample depends on them. If they must change before any expert has started, delete the cohort's `labels/<name>.csv.sha256` and its `experts/<name>/` folder, and record why in your lab notes.
- **OpenAlex errors.** The pipeline retries automatically. Rerunning `retrieve` continues where it stopped.
- **Any other error.** Rerun the same command with `NCR_DEBUG=1` in front (for example `NCR_DEBUG=1 ncr hull`) to get the full traceback. Send me the last 30 lines and the command you ran.
