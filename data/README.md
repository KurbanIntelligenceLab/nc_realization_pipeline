# Inputs

These are the inputs of the analysis, converted from the public sources and compressed with
`gzip -n`. `CHECKSUMS.txt` lists the SHA-256 of each uncompressed file, which matches the
checksums recorded during the study.

| File | Content | Source and license |
|---|---|---|
| `snapshots/oqmd_v1_0.csv.gz` | every entry of OQMD v1.0 (April 2014): formula, formation energy, displayed stability, origin, prototype | OQMD, CC BY 4.0 |
| `snapshots/oqmd_v1_1.csv.gz` | the same for OQMD v1.1 (October 2016) | OQMD, CC BY 4.0 |
| `snapshots/mp_2022_10_28.csv.gz` | every material of the Materials Project build of 28 October 2022 (GGA/GGA+U hull) | Materials Project, CC BY 4.0 |
| `experiments/experimental_records.csv.gz` | 624,213 dated experimental-structure records: COD (export of 25 October 2026), Materials Project materials with ICSD identifiers and dated references, OQMD v1.8 ICSD-derived entries | COD; Materials Project and OQMD as above |

From the ICSD-derived records, only identifiers, compositions and publication years are
included.

## Snapshot columns

| Column | Meaning |
|---|---|
| `entry_id` | identifier in the source (prefixed `oqmd-` or `mp-`); `ref-<El>` rows are elemental references at zero formation energy |
| `formula` | composition |
| `energy_per_atom` | formation energy (eV/atom) |
| `is_hypothetical` | 1 if the source presented the entry as hypothetical |
| `displayed_stable`, `displayed_dhd` | stability label and decomposition energy (eV/atom) as the release displayed them; empty if not displayed |
| `spacegroup_number`, `prototype`, `library` | structure information |
| `origin_icsd_id` | set for entries of experimental origin |
| `source_known_composition` | 1 if the release lists the same composition as experimental |

## Experimental-record columns

| Column | Meaning |
|---|---|
| `exp_id`, `database` | record identifier and source (`COD`, `MP-ICSD`, `OQMD-ICSD`) |
| `formula` | composition of the recorded structure |
| `first_report_date` | publication year (dates with only a year are read as 1 July) |
| `experimental` | 1 for experimental records (theoretical and retracted COD entries were removed) |
| `high_pressure`, `disordered` | COD flags (pressure above 1 GPa; disorder) |

To rebuild these files from the original sources, see `docs/REPRODUCE.md`, level 3.
