"""Record-level audit of the prior-record matches (exploratory; Supplementary Note 7).

A stratified random sample of 60 already-known predictions (10 per release x earliest-record source,
seed 20261026) was checked by one author against the original papers: the composition each paper
reports, its publication year, and whether its structure was experimentally determined
(data/audit/prior_record_audit.csv). A match is confirmed under the paper's rule when the reported
composition has the predicted element set and lies within the composition tolerance (0.02, the same
distance the pipeline uses) and the structure is experimentally determined; it is confirmed exactly when
the auditor found the exact composition. Per-stratum confirmation rates are weighted by stratum size;
uncertainty uses a Jeffreys Beta posterior per stratum, combined by simulation (95% equal-tailed).

Usage (from the repository root):
    python3 scripts/audit_estimates.py [--derived data_derived] [--audit data/audit/prior_record_audit.csv]
Writes data_derived/checks/audit_estimates.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from ncr.identity import composition_distance  # noqa: E402

TOL = 0.02
RELEASES = {"oqmd_v1_0": "oqmd", "oqmd_v1_1": "oqmd", "mp_2022_10_28": "materials_project"}


def within_rule(pred, paper):
    from pymatgen.core import Composition
    if not paper:
        return False
    a, b = Composition(pred), Composition(paper)
    return set(a.elements) == set(b.elements) and composition_distance(a, b) <= TOL


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data_derived")
    ap.add_argument("--audit", default="data/audit/prior_record_audit.csv")
    ap.add_argument("--draws", type=int, default=200000)
    a = ap.parse_args(argv)
    rows = list(csv.DictReader(open(a.audit, newline="")))
    for r in rows:
        exp = r["experimentally_determined"] == "Yes"
        r["confirmed_rule"] = exp and within_rule(r["predicted_composition"], r["composition_in_paper_parseable"])
        r["confirmed_exact"] = exp and r["composition_exact_match"] == "Yes"
    # stratum sizes from the released ledger: already-known predictions by release and earliest-record source
    import datetime as dt
    size = defaultdict(int)
    known = {}
    for stem in ("oqmd", "materials_project"):
        dbm = {x["pred_id"]: x for x in csv.DictReader(open(os.path.join(a.derived, "matches", f"{stem}.csv")))}
        for x in csv.DictReader(open(os.path.join(a.derived, "cohorts", f"{stem}.csv"))):
            m = dbm[x["pred_id"]]
            if m["db_first_T3"] and dt.date.fromisoformat(m["db_first_T3"]) <= dt.date.fromisoformat(x["entry_date"]):
                cod = m.get("db_first_T3__COD", "")
                src = "COD" if cod and cod == m["db_first_T3"] else "ICSD-derived"
                size[(x["entry_release"], src)] += 1
                known.setdefault(x["entry_release"], [0, 0])
                known[x["entry_release"]][0] += 1
            known.setdefault(x["entry_release"], [0, 0])
            known[x["entry_release"]][1] += 1
    rng = np.random.default_rng(20261026)
    out = {"note": "exploratory; one-author record-level audit; added after the outcomes were measured",
           "n_audited": len(rows),
           "confirmed_rule": sum(r["confirmed_rule"] for r in rows),
           "confirmed_exact": sum(r["confirmed_exact"] for r in rows),
           "experimentally_determined": sum(r["experimentally_determined"] == "Yes" for r in rows),
           "composition_exact": sum(r["composition_exact_match"] == "Yes" for r in rows),
           "strata": {}, "releases": {}}
    by = defaultdict(list)
    for r in rows:
        rel, src = [s.strip() for s in r["stratum"].split("|")]
        by[(rel, src)].append(r)
    for (rel, src), rs in sorted(by.items()):
        out["strata"][f"{rel} | {src}"] = {"size": size[(rel, src)], "audited": len(rs),
                                          "confirmed_rule": sum(r["confirmed_rule"] for r in rs),
                                          "confirmed_exact": sum(r["confirmed_exact"] for r in rs)}
    total = {"rule": np.zeros(a.draws), "exact": np.zeros(a.draws)}
    for rel in RELEASES:
        k_known, n_pred = known[rel]
        res = {"predictions": n_pred, "already_known": k_known}
        for crit in ("rule", "exact"):
            draws = np.zeros(a.draws)
            point = 0.0
            for (r2, src), rs in by.items():
                if r2 != rel:
                    continue
                w = size[(rel, src)] / k_known
                k = sum(r[f"confirmed_{crit}"] for r in rs)
                n = len(rs)
                point += w * k / n
                draws += w * rng.beta(k + 0.5, n - k + 0.5, a.draws)
            lo, hi = np.percentile(draws, [2.5, 97.5])
            res[f"confirmed_share_{crit}"] = [point, lo, hi]
            res[f"adjusted_preexistence_{crit}"] = [point * k_known / n_pred, lo * k_known / n_pred, hi * k_known / n_pred]
            res[f"adjusted_count_{crit}"] = [point * k_known, lo * k_known, hi * k_known]
            total[crit] += draws * k_known
        out["releases"][rel] = res
    for crit in ("rule", "exact"):
        pt = sum(out["releases"][r][f"adjusted_count_{crit}"][0] for r in RELEASES)
        lo, hi = np.percentile(total[crit], [2.5, 97.5])
        out[f"adjusted_entries_{crit}"] = [pt, lo, hi]
    os.makedirs(os.path.join(a.derived, "checks"), exist_ok=True)
    with open(os.path.join(a.derived, "checks", "audit_estimates.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"audited {out['n_audited']}: confirmed under the rule {out['confirmed_rule']}, exactly {out['confirmed_exact']}; "
          f"experimental {out['experimentally_determined']}; exact composition {out['composition_exact']}")
    for s, v in out["strata"].items():
        print(f"  {s:30s} size {v['size']:5d}: rule {v['confirmed_rule']}/{v['audited']}, exact {v['confirmed_exact']}/{v['audited']}")
    for rel, v in out["releases"].items():
        for crit in ("rule", "exact"):
            p, lo, hi = v[f"adjusted_preexistence_{crit}"]
            c = v[f"confirmed_share_{crit}"]
            print(f"  {rel:14s} {crit:5s}: confirmed {100*c[0]:.1f}% ({100*c[1]:.1f}-{100*c[2]:.1f}); "
                  f"pre-existence {100*v['already_known']/v['predictions']:.1f}% -> {100*p:.1f}% ({100*lo:.1f}-{100*hi:.1f})")
    for crit in ("rule", "exact"):
        pt, lo, hi = out[f"adjusted_entries_{crit}"]
        print(f"  adjusted entries ({crit}): {pt:.0f} ({lo:.0f}-{hi:.0f}) of 2,203")
    return out


if __name__ == "__main__":
    main()
