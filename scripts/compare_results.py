#!/usr/bin/env python3
"""Compare a run's results with the headline numbers reported in the paper.

Usage: python3 scripts/compare_results.py [WORKDIR]   (default: work)
Exit status 0 if every number matches, 1 otherwise."""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def collect(work: str) -> dict:
    S = load(os.path.join(work, "results", "summary.json"))
    R = load(os.path.join(work, "results", "robustness.json"))
    P = load(os.path.join(work, "checks", "premise.json"))
    K = load(os.path.join(work, "checks", "coverage_known.json"))
    W = load(os.path.join(work, "checks", "power.json"))
    B, T = S["base"], S["threshold"]["threshold__oqmd"]
    sch = W["OQMD"]["block_scheme_candidates_25meV"]["schemes"]["class_en_nel"]
    mp_r = lambda v: v["materials_project"]["materials_project::mp_2022_10_28"]  # noqa: E731
    return {
        "pre_existence": {k: {"k": B[k]["pre_existence_all"]["k"], "n": B[k]["pre_existence_all"]["n"]}
                          for k in ("oqmd::oqmd_v1_0", "oqmd::oqmd_v1_1", "materials_project")},
        "realized": {"oqmd::oqmd_v1_0": {h: B["oqmd::oqmd_v1_0"]["horizons"][h]["realized_share"]["k"] for h in ("3", "5", "10")},
                     "oqmd::oqmd_v1_1": {h: B["oqmd::oqmd_v1_1"]["horizons"][h]["realized_share"]["k"] for h in ("3", "5")},
                     "materials_project": {"2": B["materials_project"]["horizons"]["2"]["realized_share"]["k"]}},
        "capture_recapture": {k: {q: S["capture_recapture"][k][q] for q in ("n_COD", "n_ICSD_derived", "both")}
                              for k in ("oqmd::oqmd_v1_0", "materials_project")},
        "threshold_25": {**{q: T["windows"]["25"][q] for q in ("N_A", "N_B", "events_A", "events_B")},
                         "effect_estimate": T["windows"]["25"]["effect"]["estimate"]},
        "density_test": {q: T["density_test"][q] for q in ("n_left", "n_right")},
        "premise": {k: {q: v[q] for q in ("within_25_n", "within_25_rate", "check_a_passes")} for k, v in P.items()},
        "cod_coverage_of_known": {k: {q: K[k][q] for q in ("icsd_origin_compositions", "found_in_reports")}
                                  for k in ("oqmd_v1_0", "mp_2022_10_28")},
        "power_primary_blocks": {"max_abs_blocked_smd_25": sch["max_abs_blocked_smd"],
                                 "N_A_usable": sch["N_A_usable"], "N_B_usable": sch["N_B_usable"]},
        "robustness": {"composition_tolerance_halved_mp_pre": mp_r(R["rematch"]["composition_tolerance_halved"])["pre_existence"]["k"],
                       "composition_tolerance_doubled_mp_pre": mp_r(R["rematch"]["composition_tolerance_doubled"])["pre_existence"]["k"],
                       "composition_tolerance_doubled_oqmd_v1_0_realized_10y":
                           R["rematch"]["composition_tolerance_doubled"]["oqmd"]["oqmd::oqmd_v1_0"]["realized_10y"]["k"]},
    }


def diff(a, b, path=""):
    out = []
    if isinstance(b, dict):
        for k, v in b.items():
            if k.startswith("_"):
                continue
            out += diff(a.get(k) if isinstance(a, dict) else None, v, f"{path}/{k}")
    elif isinstance(b, float):
        if not isinstance(a, (int, float)) or abs(a - b) > 1e-9 * max(1.0, abs(b)):
            out.append(f"{path}: got {a}, paper {b}")
    elif a != b:
        out.append(f"{path}: got {a}, paper {b}")
    return out


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "work")
    got = collect(work)
    ref = load(os.path.join(ROOT, "tests", "paper_numbers.json"))
    bad = diff(got, ref)
    n = sum(1 for _ in json.dumps(ref).split(":")) - 1
    if bad:
        print(f"{len(bad)} of about {n} headline numbers differ from the paper:")
        print("\n".join("  " + b for b in bad))
        sys.exit(1)
    print(f"all headline numbers match the paper ({work})")


if __name__ == "__main__":
    main()
