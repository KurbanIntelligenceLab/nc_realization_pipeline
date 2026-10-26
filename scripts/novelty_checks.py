"""Two exploratory checks on the predictions that were already known at release (pre-existence).
They were added after the outcomes were measured and are reported as exploratory.

1. Structural identity. A composition-level match shows that the composition had a recorded
   experimental structure, not that the predicted structure had. For every already-known
   prediction with a COD record dated on or before its release, this check asks whether any such
   record, matched by the paper's composition rule (tier T3), reports the space group of the
   predicted structure. COD is the only source that gives space groups; the ICSD-derived records
   do not, so predictions known only from them are not assessed.

2. Timing of the Materials Project cohort. The 2022 build is the first archived build, so it does
   not date when each prediction was first published. A record dated before 2011, before the
   project started in October 2011 (Jain et al., APL Materials 1, 011002, 2013), cannot be a
   realization of a Materials Project prediction.

Usage (from the repository root):
    python3 scripts/novelty_checks.py [--derived data_derived] [--records data/experiments/experimental_records.csv.gz]
Writes data_derived/checks/novelty_checks.json and prints a summary.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from ncr import matching  # noqa: E402
from ncr.util import parse_date  # noqa: E402


def _decimal_year(d):
    import datetime as dt
    y0 = dt.date(d.year, 1, 1)
    return d.year + (d - y0).days / (dt.date(d.year + 1, 1, 1) - y0).days


def _read_records(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        text = fh.read()
    rows = list(csv.DictReader(io.StringIO(text)))
    sg = {r["exp_id"]: r.get("sg_number", "") for r in rows}
    tmp = os.path.join(os.path.dirname(os.path.abspath(path)), ".records_for_checks.csv")
    with open(tmp, "w") as out:
        out.write(text)
    try:
        exps = matching.read_experiments(tmp)
    finally:
        os.remove(tmp)
    return exps, sg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data_derived")
    ap.add_argument("--records", default="data/experiments/experimental_records.csv.gz")
    a = ap.parse_args(argv)
    exps, sg_of = _read_records(a.records)
    cod_by_sys = matching.build_exp_index([e for e in exps if e["database"] == "COD"])
    out = {"note": "exploratory; added after the outcomes were measured", "structure": {}, "timing": {}}
    for stem in ("oqmd", "materials_project"):
        cohort = {r["pred_id"]: r for r in csv.DictReader(open(os.path.join(a.derived, "cohorts", f"{stem}.csv")))}
        dbm = {r["pred_id"]: r for r in csv.DictReader(open(os.path.join(a.derived, "dbmatch", f"{stem}.csv")))}
        by_release = defaultdict(list)
        for pid, r in cohort.items():
            by_release[r["entry_release"]].append(r)
        for rel, rows in sorted(by_release.items()):
            known = known_cod = same_sg = 0
            pre2011 = 0
            for r in rows:
                e, _ = parse_date(r["entry_date"])
                first, _ = parse_date(dbm[r["pred_id"]]["db_first_T3"]) if dbm[r["pred_id"]]["db_first_T3"] else (None, None)
                if first is None or first > e:
                    continue
                known += 1
                if first.year < 2011:
                    pre2011 += 1
                cod_first = dbm[r["pred_id"]].get("db_first_T3__COD", "")
                if not cod_first or parse_date(cod_first)[0] > e:
                    continue
                known_cod += 1
                pred_sg = str(r["spacegroup_number"]).strip()
                for x in cod_by_sys.get(r["chemsys"], []):
                    if x["date"] and x["date"] <= e and str(sg_of.get(x["exp_id"], "")).strip() == pred_sg \
                            and matching.match_tier(r["formula"], "", x) is not None:
                        same_sg += 1
                        break
            out["structure"][rel] = {"already_known": known, "with_COD_record_by_release": known_cod,
                                     "COD_record_in_predicted_space_group": same_sg,
                                     "share": same_sg / known_cod if known_cod else None}
            if stem == "materials_project":
                out["timing"][rel] = {"already_known": known, "earliest_record_before_2011": pre2011,
                                      "share": pre2011 / known if known else None}
    os.makedirs(os.path.join(a.derived, "checks"), exist_ok=True)
    with open(os.path.join(a.derived, "checks", "novelty_checks.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    for rel, v in out["structure"].items():
        print(f"{rel}: {v['COD_record_in_predicted_space_group']} of {v['with_COD_record_by_release']} already-known "
              f"predictions with a COD record by release had one in the predicted space group "
              f"({100 * v['share']:.1f}%)" if v["share"] is not None else f"{rel}: none")
    for rel, v in out["timing"].items():
        print(f"{rel}: {v['earliest_record_before_2011']} of {v['already_known']} already-known compositions "
              f"have an earliest record dated before 2011 ({100 * v['share']:.1f}%)")
    return out


if __name__ == "__main__":
    main()
