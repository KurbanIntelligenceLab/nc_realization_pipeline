#!/usr/bin/env bash
# Rebuild the released inputs in data/ from the original public sources. This is optional:
# data/ already holds them. It needs about 100 GB of free disk, a local MySQL server, the
# AWS command-line client and several hours (mostly downloading and loading OQMD).
#
#   RAW=/path/with/space scripts/rebuild_inputs.sh
#
# COD changes daily, so an export made on a different day than the study's (25 October 2026)
# can give slightly different matches. The SHA-256 of every file the study used is listed in
# data/CHECKSUMS.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
RAW="${RAW:-raw}"
mkdir -p "$RAW/oqmd" "$RAW/mp" "$RAW/cod" data/snapshots data/experiments
ncr() { "$PY" -m ncr "$@"; }

# 1. OQMD releases v1.0 (April 2014), v1.1 (October 2016) and v1.8 (February 2026)
for f in qmdb__v1_0__042014.sql.gz qmdb__v1_1__102016.sql.gz qmdb__v1_8__022026.sql.gz; do
  [ -f "$RAW/oqmd/$f" ] || curl -L --fail --retry 10 -C - -o "$RAW/oqmd/$f" "https://static.oqmd.org/static/downloads/$f"
done
shasum -a 256 -c data/CHECKSUMS_RAW.txt --ignore-missing
for v in v1_0__042014 v1_1__102016 v1_8__022026; do
  db="qmdb_${v%%__*}"
  mysql -u root -e "DROP DATABASE IF EXISTS $db; CREATE DATABASE $db"
  gunzip -c "$RAW/oqmd/qmdb__${v}.sql.gz" | mysql -u root "$db"
done
for db in qmdb_v1_0 qmdb_v1_1; do
  mysql -u root --quick -B "$db" -e "SELECT e.id AS entry_id, e.composition_id, e.path, e.label, e.duplicate_of_id, e.prototype_id, e.ntypes, e.natoms, e.delta_e AS entry_delta_e, e.stability AS entry_stability, fe.id AS fe_id, fe.fit_id, fe.delta_e AS fe_delta_e, fe.stability AS fe_stability FROM entries e LEFT JOIN formation_energies fe ON fe.entry_id = e.id" > "$RAW/oqmd/${db}_entries.tsv"
  mysql -u root --quick -B "$db" -e "SELECT entry_id, id AS structure_id, label, spacegroup_id FROM structures WHERE entry_id IS NOT NULL" > "$RAW/oqmd/${db}_structures.tsv"
  ncr convert-oqmd --entries "$RAW/oqmd/${db}_entries.tsv" --structures "$RAW/oqmd/${db}_structures.tsv" --out "$RAW/oqmd/${db#qmdb_}.csv"
done
mysql -u root --quick -B qmdb_v1_8 -e "SELECT e.id AS entry_id, e.composition_id, p.year FROM entries e LEFT JOIN publications p ON p.id = e.reference_id WHERE e.path LIKE '%/icsd/%'" > "$RAW/oqmd/qmdb_v1_8_icsd.tsv"
ncr convert-oqmd-icsd --in "$RAW/oqmd/qmdb_v1_8_icsd.tsv" --out "$RAW/oqmd/oqmd_icsd.csv"

# 2. Materials Project build of 28 October 2022 (AWS Open Data, no account needed)
for c in thermo provenance; do
  aws s3 sync --no-sign-request "s3://materialsproject-build/collections/2022-10-28/$c/" "$RAW/mp/2022-10-28/$c/" --only-show-errors
done
ncr convert-mp --build-dir "$RAW/mp/2022-10-28" --out "$RAW/mp/mp_2022_10_28.csv"
ncr convert-mp-icsd --build-dir "$RAW/mp/2022-10-28" --out "$RAW/mp/mp_icsd.csv"

# 3. COD data table from its public server
ncr export-cod --out "$RAW/cod/cod_data.tsv"
ncr convert-cod --in "$RAW/cod/cod_data.tsv" --out "$RAW/cod/cod_experiments.csv"

# 4. Merge the outcome records and write the compressed inputs
ncr merge-experiments --in "$RAW/cod/cod_experiments.csv" "$RAW/mp/mp_icsd.csv" "$RAW/oqmd/oqmd_icsd.csv" \
    --out "$RAW/experimental_records.csv"
gzip -n -c "$RAW/oqmd/v1_0.csv" > data/snapshots/oqmd_v1_0.csv.gz
gzip -n -c "$RAW/oqmd/v1_1.csv" > data/snapshots/oqmd_v1_1.csv.gz
gzip -n -c "$RAW/mp/mp_2022_10_28.csv" > data/snapshots/mp_2022_10_28.csv.gz
gzip -n -c "$RAW/experimental_records.csv" > data/experiments/experimental_records.csv.gz
echo "inputs rebuilt in data/; compare with data/CHECKSUMS.txt"
