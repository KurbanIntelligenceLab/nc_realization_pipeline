# nc_realization_pipeline

Code and derived data for "Computational predictions of stable inorganic materials and their later experimental structures" (Hasan Kurban and Mustafa Kurban).

The pipeline follows every composition that three dated releases displayed as stable and new: OQMD v1.0 (April 2014), OQMD v1.1 (October 2016) and the Materials Project build of 28 October 2022. It matches each one deterministically against dated experimental crystal-structure records from three sources:

- the Crystallography Open Database (COD);
- the ICSD-derived records of the Materials Project;
- the ICSD-derived records of OQMD.

No language model and no human annotation take part in the measurement.

## Contents

| Path | What it is |
|---|---|
| `ncr_cli.py`, `ncr/` | the pipeline (version in `ncr/__init__.py`) |
| `nc_realization_metrics.py`, `nc_realization_cohorts.py`, `nc_realization_figures.py` | estimators, identity tiers and figure script |
| `tests/test_pipeline.py` | end-to-end tests on a synthetic world (`python3 tests/test_pipeline.py`) |
| `RUN_ON_MAC.md`, `DATA_FORMATS.md` | how to rebuild every input and run every step |
| `config.example.json` | configuration template |
| `data_derived/` | the ledger and all aggregate results of the study (see its README) |

## Reproducing the study

1. Install the pinned dependencies from `requirements-pipeline.txt` (Python 3.11).
2. Follow `RUN_ON_MAC.md` to download and convert the public inputs. Checksums are recorded in `data_derived/manifest.json`.
3. Run the steps in order: `convert-*`, `hull`, `cohorts`, `x10`, `premise`, `export-cod`/`convert-cod`, `convert-oqmd-icsd`, `convert-mp-icsd`, `merge-experiments`, `dbmatch --cohort all`, `estimate-census` and `robustness`.

The raw database files are not redistributed here because they are public at their sources: the OQMD MySQL releases, the Materials Project open-data archive, and COD.

## License

MIT; see `LICENSE`.
