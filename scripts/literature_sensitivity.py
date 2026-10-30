"""Literature-informed sensitivity analysis of the threshold design (exploratory; Supplementary Note 6).

The primary analysis counts a structure as recorded when its earliest matching record in the analyzed
sources (COD and the ICSD-derived records) falls within five years after the release. The literature
audit (data/audit/literature_audit.csv) checked all four events in the +/-25 meV/atom window. Under the
paper's outcome definition, an experimentally determined crystal structure, none is reclassified. This
script recodes as non-events the events for which any earlier experimental report was found, as the
design treats any compound with an earlier record, and reruns the blocked estimator and the exact test.
For both, it also reports the upper limit of the Wilson 95% interval for the stable-side recording
probability (exploratory; quoted in Results as the upper 95% bound).

Usage (from the repository root):
    python3 scripts/literature_sensitivity.py [--derived data_derived] [--audit data/audit/literature_audit.csv]
Writes data_derived/checks/literature_sensitivity.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from ncr import stats as S  # noqa: E402
from ncr.analysis import Outcomes  # noqa: E402


def _exact_one_sided(k, n, p):
    """P(at least k of n events on the stable side), conditional on n events, under no label effect."""
    return sum(math.comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(k, n + 1))


def _summary(y, A, B, block):
    e = S.blocked_difference(y, A, B, block)
    nA, nB, kA, kB = int(A.sum()), int(B.sum()), int(y[A].sum()), int(y[B].sum())
    return {"events_stable": kA, "n_stable": nA, "events_unstable": kB, "n_unstable": nB,
            "estimate_pp": 100 * e["estimate"], "ci_low_pp": 100 * e["ci_low"], "ci_high_pp": 100 * e["ci_high"],
            "exact_one_sided_p": _exact_one_sided(kA, kA + kB, nA / (nA + nB)),
            # upper limit of the Wilson 95% interval for the stable-side recording probability, kA of nA,
            # in percent; the paper reports it rounded outward (0.6202 -> 0.63)
            "stable_side_upper_95_pct": 100 * S.share(kA, nA)["hi"]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data_derived")
    ap.add_argument("--audit", default="data/audit/literature_audit.csv")
    ap.add_argument("--config", default=None)
    a = ap.parse_args(argv)
    cfg = json.load(open(a.config or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config", "study.json")))
    o = Outcomes(cfg, os.path.join(a.derived, "cohorts", "threshold__oqmd.csv"),
                 os.path.join(a.derived, "matches", "threshold__oqmd.csv"))
    h = float(cfg.get("threshold", {}).get("horizon_years", 5))
    x = np.array([float(r["x_mev"]) for r in o.R])
    side = np.array([r["side"] for r in o.R])
    block = np.array([r["block"] for r in o.R])
    y = o.within(o.m, o.e, o.e + h)
    inw = np.abs(x) <= 25
    A, B = inw & (side == "A"), inw & (side == "B")
    ids = [r["pred_id"] for r in o.R]
    events = [ids[i] for i in np.where(inw & (y == 1))[0]]
    audit = {r["pred_id"]: r for r in csv.DictReader(open(a.audit, newline=""))}
    unaudited = [e for e in events if e not in audit]
    if unaudited:
        raise SystemExit(f"events without an audit entry: {unaudited}")
    recode = [e for e in events if audit[e]["classification_if_any_earlier_experimental_report_counts"].startswith("non-event")]
    y2 = y.copy()
    for e in recode:
        y2[ids.index(e)] = 0.0
    out = {"note": "exploratory; added after the outcomes were measured",
           "events_in_primary_window": events, "recoded_as_non_events": recode,
           "database_defined": _summary(y, A, B, block),
           "any_earlier_experimental_report_counts": _summary(y2, A, B, block)}
    os.makedirs(os.path.join(a.derived, "checks"), exist_ok=True)
    with open(os.path.join(a.derived, "checks", "literature_sensitivity.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    for k in ("database_defined", "any_earlier_experimental_report_counts"):
        v = out[k]
        print(f"{k}: {v['events_stable']}/{v['n_stable']} vs {v['events_unstable']}/{v['n_unstable']}; "
              f"estimate {v['estimate_pp']:.2f} pp (95% CI {v['ci_low_pp']:.3f} to {v['ci_high_pp']:.3f}); "
              f"exact one-sided p = {v['exact_one_sided_p']:.3f}")
    return out


if __name__ == "__main__":
    main()
