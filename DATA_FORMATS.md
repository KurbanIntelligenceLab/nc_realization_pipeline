# Input file formats

All tables are CSV files with a header row, encoded as UTF-8. Dates use `YYYY-MM-DD`; `YYYY` alone is accepted and set to 1 July. Paths inside files may be absolute, or relative to the file that contains them.

## 1. Release snapshots

Use one file per database release. List each one under `releases` in `config.json`.

| Column | Required | Meaning |
|---|---|---|
| `entry_id` | yes | The database's identifier for the entry. It must be unique within the file and stable across releases of the same source. |
| `formula` | yes | The composition as computed, for example `Li2FeO3`. |
| `energy_per_atom` | yes | The energy used for the convex hull, in eV/atom. Use the release's own corrected energy, so that the hull matches what the release displayed. |
| `is_hypothetical` | yes | `1` if the source presented the entry as hypothetical (for example, the Materials Project "theoretical" flag or an OQMD prototype decoration), else `0`. |
| `displayed_stable` | recommended | `1` if the release displayed the entry as stable. The threshold design's side is this displayed label (Table S1). |
| `displayed_dhd` | optional | The displayed energy above the hull, in eV/atom, kept for comparison. |
| `structure_file` | optional | A CIF of the computed structure, used for identity tiers T1 and T2. |
| `spacegroup_number` | optional | 1 to 230. Used by the model-side structure proxy and by `premise` in its `composition+spacegroup` matching mode. |
| `prototype` | optional | Prototype name. The X0 count of Heusler decorations looks for "heusler" in this column. |
| `origin_icsd_id` | optional | If non-empty, the entry is of experimental origin and is excluded from cohorts (Table S1, "Cohort hygiene"). |
| `deprecated` | optional | `1` if the entry was later deprecated. Deprecated entries are kept and flagged. |

**The snapshot must include every computed entry of that release**, including elements, known compounds and unstable entries, because the hull is recomputed from it. If your column names differ, map them in the config: `"columns": {"energy_per_atom": "e_corrected", ...}`.

**Entry dumps from pymatgen.** If a release comes as pymatgen `ComputedEntry` or `ComputedStructureEntry` dictionaries (a JSON list, optionally `.gz` or `.bz2`), convert it:

    python3 ncr_cli.py convert-entries --in mp_entries.json.gz --out data/snapshots/mp_X.csv \
        --hypothetical-ids mp_theoretical_ids.txt --cif-dir data/cifs/mp_X

A generic dump does not record which entries are hypothetical. Provide that list from the source itself; the code never guesses it.

**OQMD MySQL dumps.** Export two tab-separated files from the loaded database (commands in `RUN_ON_MAC.md`, section 3b) and convert them:

    python3 ncr_cli.py convert-oqmd --entries qmdb_v1_0_entries.tsv --structures qmdb_v1_0_structures.tsv \
        --out data/snapshots/oqmd_v1_0.csv

The converter applies these rules, fixed from diagnostics of v1.0 and v1.1:
- **Fit:** the `standard` fit, the only fit in either release. When an entry has several rows, the lowest formation energy is kept: every row is a point on the hull, so an entry is stable exactly when its lowest row is.
- **Duplicates:** entries marked as duplicates are dropped.
- **Experimental origin:** paths under `/libraries/icsd/`.
- **Test entries:** `/test/` paths are kept in the hull but never enter a cohort.
- **Heusler type:** taken from the library path.
- **References:** elemental references are added at zero formation energy.
- **Source-known compositions:** `source_known_composition` marks hypothetical entries whose composition the release already lists from ICSD, including ICSD entries without a formation energy. These are excluded from cohorts by default, because the source did not present them as new.

It also writes a `_conversion_report.json` with every count.

## 2. Experimental reports (`experiments.path`)

This is your export from ICSD or COD. Use one structure database per analysis series.

| Column | Required | Meaning |
|---|---|---|
| `exp_id` | yes | The database identifier. |
| `database` | yes | `ICSD` or `COD`. |
| `formula` | yes | The reported composition. |
| `first_report_date` | yes | The publication date of the report (first online publication if known). |
| `experimental` | yes | `1` for experimental structures. Rows with `0` (theoretical structures) are ignored, as Table S1 requires. |
| `structure_file` | optional | A CIF, which enables tiers T1 and T2 against predicted structures. |
| `ss_end_a`, `ss_end_b` | optional | End members of a reported continuous solid solution (the T2 segment rule). |
| `high_pressure`, `thin_film_only` | optional | `1` flags the report. |

**ICSD is licensed.** Keep this file on your machine. The pipeline sends nothing from it anywhere, and `results_for_claude.json` contains aggregates only.

## 3. Literature corpus (optional; `retrieval.backend = "local"`)

Use a JSON Lines file with one document per line:

    {"doc_id": "...", "doi": "...", "title": "...", "text": "full text or abstract",
     "publication_date": "YYYY-MM-DD", "received_date": "YYYY-MM-DD"}

Use this for licensed full text. With `"backend": "openalex"`, the pipeline queries the public OpenAlex index instead, using titles and abstracts only, and needs your e-mail in the config.

## 4. Control panel (`control_panel.path`, experiment X1)

| Column | Meaning |
|---|---|
| `formula` | The compound. |
| `panel_kind` | `positive` (a documented case of prediction followed by synthesis), `old` (a long-known compound) or `decoy` (more than 300 meV/atom above the hull, so no report is expected). |
| `entry_date` | The prediction date for positives; any fixed date for the other kinds. |
| `known_first_report_date` | The documented report date (positives and old compounds). |
| `granularity`, `structure_file`, `spacegroup_number` | Optional. |

## 5. Relabeling design (`relabeling`, experiment E7)

For each pair of consecutive Materials Project builds:
- `release_r` and `release_r1` are release IDs from `releases`.
- `raw_entries_r1` is a JSON list of `ComputedStructureEntry` dictionaries for build r+1. They must contain the calculation parameters the correction schemes need (run type, Hubbard values, POTCAR specification). Entries without them are dropped by pymatgen and appear as `flip_unclassifiable`.
- `scheme_r` and `scheme_r1` are the pymatgen compatibility class names for each build. Confirm both against the builds' documentation.
- `corrected_species` lists the species whose correction changed between the two schemes. It is author-verified and is not inferred by the code.

## 6. Files the pipeline writes for people

### Expert sheets: `work/experts/<cohort>/expert_A.csv` and `expert_B.csv`

Give one to each expert.

**Each expert fills:**
- `first_report_date`: the earliest matching report (`YYYY-MM-DD`, or `YYYY`);
- `tier`: `T1`, `T2`, `T3` or `none`;
- `first_report_reference`;
- `first_report_received_date` (if known);
- `high_pressure` and `thin_film_only` (`yes` or `no`);
- `structure_decision_without_refined_structure` (`yes` or `no`);
- `first_mention_date`: the earliest publication of any kind that mentions the compound;
- `minutes_spent`;
- `sources_searched`.

**Never share `sample_key.csv` with the experts**, because it links the sheets to model labels.

### Adjudication: `adjudication.csv`

The third expert fills `final_has_report` (`yes` or `no`), `final_first_report_date` and `final_first_mention_date` for each disagreement.

### Residual-miss audit (X3): `residual_audit_sheet.csv`

The third expert fills `found_missed_report` (`yes` or `no`), `report_date` and `reference`.
