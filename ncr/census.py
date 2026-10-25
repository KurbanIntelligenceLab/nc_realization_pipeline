"""Deterministic (annotation-free) analysis.

Outcomes come only from matching each prediction to dated experimental-structure records
(COD, MP-ICSD, OQMD-ICSD) in `dbmatch`. Every prediction is "labeled", so shares are exact
for the cohort; intervals are Wilson 95% intervals (superpopulation view) and are reported as
such. The threshold and relabeling designs reuse the same estimator with a complete census,
so only the assignment (randomization) variance remains. Capture-recapture between COD and
the ICSD-derived sources estimates coverage for reports dated up to the last year both
cover (2016), under an independence assumption that is stated with every result.
"""
from __future__ import annotations

import glob
import math
import os

import numpy as np

from . import config as C
from .estimate import Bundle, relabel, threshold, _clip
from .util import decimal_year, log, parse_date, read_csv, versions, write_json

import sys
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import nc_realization_metrics as NCM   # noqa: E402

ICSD_DERIVED = ("MP-ICSD", "OQMD-ICSD")


def _y(s):
    d, _ = parse_date(s)
    return decimal_year(d) if d else math.nan


class CensusBundle(Bundle):
    """A Bundle in which every unit is observed and the outcome is the matched record."""

    def __init__(self, cfg, cohort_file, dbmatch_file, sources=None, row_filter=None):
        self.cfg = cfg
        self.name = os.path.splitext(os.path.basename(cohort_file))[0]
        rows = read_csv(cohort_file)
        if row_filter is not None:
            rows = [r for r in rows if row_filter(r)]
        db = {x["pred_id"]: x for x in read_csv(dbmatch_file)}
        self.R = rows
        self.labs = [{"pred_id": r["pred_id"]} for r in rows]
        self.N = len(rows)
        self.kappa = None
        self.strata = np.array(["census"] * self.N, dtype=object)
        self.idx = np.arange(self.N)
        self.e = np.array([_y(r["entry_date"]) for r in rows])

        def first(r, srcs):
            x = db.get(r["pred_id"], {})
            tier = "T2" if r.get("granularity", "structure") == "structure" else "T3"
            if srcs is None:
                return _y(x.get(f"db_first_{tier}", ""))
            ds = [_y(x.get(f"db_first_{tier}__{s}", "")) for s in srcs]
            ds = [d for d in ds if np.isfinite(d)]
            return min(ds) if ds else math.nan
        self.first = lambda srcs=None: np.array([first(r, srcs) for r in rows])
        self.m = self.first(sources)
        self.m_rec = self.m
        self.m_men = np.full(self.N, np.nan)
        self.h, self.h_rec, self.h_men, self.e_s = self.m, self.m_rec, self.m_men, self.e


def _wilson(k, n):
    if n == 0:
        return {"est": None, "lo": None, "hi": None, "k": k, "n": n}
    lo, hi = NCM.wilson_ci(int(k), int(n))
    return {"est": k / n, "lo": lo, "hi": hi, "k": int(k), "n": int(n)}


def census_rates(cfg, b: CensusBundle, source: str) -> dict:
    cutoff = decimal_year(C.cutoff(cfg))
    ml = source in cfg.get("ml_sources", [])
    horizons = (cfg.get("horizons_by_source", {}).get(source)
                or (cfg.get("ml_horizons_years", [2]) if ml else cfg.get("horizons_years", [3, 5, 10])))
    pre = b.pre(b.m, b.e)
    res = {"source": source, "N": b.N, "pre_existence_all": _wilson(pre.sum(), b.N), "horizons": {}}
    for t in horizons:
        at = (b.e + t <= cutoff - 1.0)
        n = int(at.sum())
        if n < 2:
            res["horizons"][str(t)] = {"note": "fewer than two entries with complete follow-up"}
            continue
        Y = b.within(b.m, b.e, b.e + t)
        kp, ky = pre[at].sum(), Y[at].sum()
        rec = {"n_at_risk": n, "pre_existence_within_R": _wilson(kp, n), "realized_share": _wilson(ky, n),
               "not_realized_share": _wilson(n - kp - ky, n), "realized_among_not_known": _wilson(ky, n - kp)}
        at2 = (b.e + t <= cutoff - 2.0)
        if at2.sum() >= 2:
            rec["sensitivity_lag_2yr"] = _wilson(Y[at2].sum(), int(at2.sum()))
        for name, srcs in (("COD_only", ("COD",)), ("ICSD_derived_only", ICSD_DERIVED)):
            mm = b.first(srcs)
            rec[f"realized_share_{name}"] = _wilson(b.within(mm, b.e, b.e + t)[at].sum(), n)
            rec[f"pre_existence_{name}"] = _wilson(b.pre(mm, b.e)[at].sum(), n)
        res["horizons"][str(t)] = rec
    return res


def capture_recapture(b: CensusBundle, end_year: float = 2017.0) -> dict:
    """Chapman estimator for compounds with a record dated before end_year in either COD
    or the ICSD-derived sources. Assumes the two sources capture independently; positive
    dependence (likely: well-studied compounds are in both) biases the estimated total
    downward, so the coverage estimates are upper bounds on true coverage."""
    a = np.isfinite(b.first(("COD",))) & (b.first(("COD",)) < end_year)
    c = np.isfinite(b.first(ICSD_DERIVED)) & (b.first(ICSD_DERIVED) < end_year)
    nA, nB, m = int(a.sum()), int(c.sum()), int((a & c).sum())
    if nA == 0 or nB == 0:
        return {"note": "a source has no records in the window", "n_COD": nA, "n_ICSD_derived": nB, "both": m}
    N = (nA + 1) * (nB + 1) / (m + 1) - 1
    var = (nA + 1) * (nB + 1) * (nA - m) * (nB - m) / ((m + 1) ** 2 * (m + 2))
    return {"window": f"records dated before {int(end_year)}", "n_COD": nA, "n_ICSD_derived": nB, "both": m,
            "union": int((a | c).sum()), "chapman_total": N, "chapman_se": math.sqrt(var),
            "coverage_COD": nA / N, "coverage_ICSD_derived": nB / N, "coverage_union": int((a | c).sum()) / N,
            "assumption": "independent capture; positive dependence makes these coverages upper bounds"}


def step_estimate_census(cfg: dict) -> dict:
    wd = cfg["workdir"]
    res = {"mode": "deterministic (no annotation)", "base": {}, "threshold": {}, "relabel": {},
           "capture_recapture": {}, "status": {}}
    for cf in sorted(glob.glob(os.path.join(wd, "cohorts", "*.csv"))):
        stem = os.path.splitext(os.path.basename(cf))[0]
        dbp = os.path.join(wd, "dbmatch", stem + ".csv")
        if not os.path.exists(dbp):
            res["status"][stem] = "no dbmatch output"
            continue
        b = CensusBundle(cfg, cf, dbp)
        if stem.startswith("threshold__"):
            t = threshold(cfg, b)
            for w in t["windows"].values():
                w.pop("attention_check", None)
            t["note"] = "attention check (X9) not measurable without literature retrieval; omitted"
            res["threshold"][stem] = t
        elif stem.startswith("relabel__"):
            res["relabel"][stem] = relabel(cfg, b)
        else:
            src = b.R[0]["source"] if b.R else stem
            res["base"][stem] = census_rates(cfg, b, src)
            res["capture_recapture"][stem] = capture_recapture(b)
            rels = sorted({r.get("entry_release", "") for r in b.R})
            if len(rels) > 1:                       # one cohort per entry release (Table 1)
                for rel in rels:
                    sb = CensusBundle(cfg, cf, dbp, row_filter=lambda r, rel=rel: r.get("entry_release") == rel)
                    if sb.N >= 2:
                        res["base"][f"{stem}::{rel}"] = census_rates(cfg, sb, f"{src} {rel}")
                        res["capture_recapture"][f"{stem}::{rel}"] = capture_recapture(sb)
        res["status"][stem] = "estimated"
    for name in ("x0", "x10", "premise", "coverage_known"):
        p = os.path.join(wd, "checks", f"{name}.json")
        if os.path.exists(p):
            from .util import read_json
            res[name] = read_json(p)
    figs = {"fig2a": [], "fig2b": []}
    per_release = [k for k in res["base"] if "::" in k]
    pooled_with_split = {k.split("::")[0] for k in per_release}
    for stem, r in res["base"].items():
        if stem in pooled_with_split:
            continue                                   # plot each release, not the pooled cohort
        hs = [h for h, v in r["horizons"].items() if "realized_share" in v and v["realized_share"]["est"] is not None]
        if not hs:
            continue
        v = r["horizons"][max(hs, key=float)]
        figs["fig2a"].append({"source": r["source"], "horizon": float(max(hs, key=float)),
                              "already_known": _clip(v["pre_existence_within_R"]), "realized": _clip(v["realized_share"]),
                              "not_realized": _clip(v["not_realized_share"])})
        for h in hs:
            rr = r["horizons"][h]["realized_among_not_known"]
            if rr["est"] is not None and rr["n"]:
                lo, hi = NCM.wilson_ci(rr["k"], rr["n"], 0.05 / len(hs))      # simultaneous across horizons
                figs["fig2b"].append({"source": r["source"], "t": float(h), "est": rr["est"], "lo": lo, "hi": hi})
    for stem, t in res["threshold"].items():
        bins = []
        for bn in t["fig3a"]["bins"]:                     # census bins have no sampling error; show Wilson intervals
            k = int(round(bn["est"] * bn["n"]))
            lo, hi = NCM.wilson_ci(k, bn["n"])
            bins.append({"center": bn["center"], "est": k / bn["n"], "lo": lo, "hi": hi, "n": bn["n"]})
        figs["fig3a"] = {"window": t["fig3a"]["window"], "bins": bins}
        break
    res["figures"] = figs
    if figs.get("fig2a") and figs.get("fig3a"):
        try:
            import nc_realization_figures as NCF
            res["figure_checks"] = NCF.render_all(figs, os.path.join(wd, "figures"), synthetic=cfg.get("_synthetic", False))
        except Exception as exc:  # noqa: BLE001
            res["figure_checks"] = {"error": str(exc)}
    res["versions"] = versions()
    write_json(os.path.join(wd, "results", "census_results.json"), res)
    for_claude = {k: res[k] for k in ("mode", "base", "threshold", "relabel", "capture_recapture", "status",
                                      "x0", "x10", "premise", "coverage_known", "figure_checks", "versions") if k in res}
    for t in for_claude.get("threshold", {}).values():
        for w in t.get("windows", {}).values():
            w.pop("blocks_used", None)
    write_json(os.path.join(wd, "results_for_claude.json"), for_claude)
    log(f"estimate-census: wrote {os.path.join(wd, 'results_for_claude.json')}")
    return res
