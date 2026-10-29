"""Estimates from the matched records: pre-existence and realization shares for every
cohort and entry release, capture-recapture coverage, and the threshold design.

Every prediction's outcome is observed, so shares are exact for the cohort and are reported
with Wilson 95% intervals (Bonferroni across horizons where several are reported together).
"""
from __future__ import annotations

import glob
import math
import os

import numpy as np

from . import checks as K
from . import config as C
from . import figures as FIG
from . import stats as S
from .util import decimal_year, log, parse_date, read_csv, read_json, versions, write_json

ICSD_DERIVED = ("MP-ICSD", "OQMD-ICSD")


def _year(s):
    d, _ = parse_date(s)
    return decimal_year(d) if d else math.nan


class Outcomes:
    """A cohort (or the threshold population) with the earliest matching record of each unit."""

    def __init__(self, cfg, cohort_file, match_file, sources=None, row_filter=None):
        self.cfg = cfg
        self.name = os.path.splitext(os.path.basename(cohort_file))[0]
        rows = read_csv(cohort_file)
        if row_filter is not None:
            rows = [r for r in rows if row_filter(r)]
        self._db = {x["pred_id"]: x for x in read_csv(match_file)}
        self.R = rows
        self.N = len(rows)
        self.e = np.array([_year(r["entry_date"]) for r in rows])
        self.m = self.first(sources)

    def first(self, srcs=None) -> np.ndarray:
        """Earliest matching record (decimal year) at the primary tier, overall or from srcs."""
        out = []
        for r in self.R:
            x = self._db.get(r["pred_id"], {})
            tier = "T2" if r.get("granularity", "structure") == "structure" else "T3"
            if srcs is None:
                out.append(_year(x.get(f"db_first_{tier}", "")))
                continue
            ds = [_year(x.get(f"db_first_{tier}__{s}", "")) for s in srcs]
            ds = [d for d in ds if np.isfinite(d)]
            out.append(min(ds) if ds else math.nan)
        return np.array(out)

    @staticmethod
    def pre(first, e):
        return np.where(np.isfinite(first) & (first <= e), 1.0, 0.0)

    @staticmethod
    def within(first, lo, hi):
        return np.where(np.isfinite(first) & (first > lo) & (first <= hi), 1.0, 0.0)


def horizons_for(cfg, source: str) -> list:
    return cfg.get("horizons_by_source", {}).get(source) or cfg.get("horizons_years", [3, 5, 10])


def cohort_rates(cfg, o: Outcomes, source: str) -> dict:
    cutoff = decimal_year(C.cutoff(cfg))
    pre = o.pre(o.m, o.e)
    res = {"source": source, "N": o.N, "pre_existence_all": S.share(pre.sum(), o.N), "horizons": {}}
    for t in horizons_for(cfg, source):
        at = (o.e + t <= cutoff - 1.0)
        n = int(at.sum())
        if n < 2:
            res["horizons"][str(t)] = {"note": "fewer than two entries with complete follow-up"}
            continue
        Y = o.within(o.m, o.e, o.e + t)
        kp, ky = pre[at].sum(), Y[at].sum()
        rec = {"n_at_risk": n, "pre_existence_within_R": S.share(kp, n), "realized_share": S.share(ky, n),
               "not_realized_share": S.share(n - kp - ky, n), "realized_among_not_known": S.share(ky, n - kp)}
        at2 = (o.e + t <= cutoff - 2.0)
        if at2.sum() >= 2:
            rec["sensitivity_lag_2yr"] = S.share(Y[at2].sum(), int(at2.sum()))
        for name, srcs in (("COD_only", ("COD",)), ("ICSD_derived_only", ICSD_DERIVED)):
            mm = o.first(srcs)
            rec[f"realized_share_{name}"] = S.share(o.within(mm, o.e, o.e + t)[at].sum(), n)
            rec[f"pre_existence_{name}"] = S.share(o.pre(mm, o.e)[at].sum(), n)
        res["horizons"][str(t)] = rec
    return res


def capture_recapture(o: Outcomes, end_year: float = 2017.0) -> dict:
    """Chapman estimator for compositions with a record dated before end_year in COD or in the
    ICSD-derived sources. Positive dependence between the sources (well-studied compounds
    appear in both) makes the estimated total too small, so the coverages are upper bounds."""
    cod, icsd = o.first(("COD",)), o.first(ICSD_DERIVED)
    a = np.isfinite(cod) & (cod < end_year)
    c = np.isfinite(icsd) & (icsd < end_year)
    n1, n2, m = int(a.sum()), int(c.sum()), int((a & c).sum())
    if n1 == 0 or n2 == 0:
        return {"note": "a source has no records in the window", "n_COD": n1, "n_ICSD_derived": n2, "both": m}
    ch = S.chapman(n1, n2, m)
    return {"window": f"records dated before {int(end_year)}", "n_COD": n1, "n_ICSD_derived": n2, "both": m,
            "union": int((a | c).sum()), "chapman_total": ch["total"], "chapman_se": ch["se"],
            "coverage_COD": n1 / ch["total"], "coverage_ICSD_derived": n2 / ch["total"],
            "coverage_union": int((a | c).sum()) / ch["total"],
            "assumption": "independent capture; positive dependence makes these coverages upper bounds"}


def threshold_analysis(cfg, o: Outcomes) -> dict:
    tc = cfg.get("threshold", {})
    h = float(tc.get("horizon_years", 5))
    x = np.array([float(r["x_mev"]) for r in o.R])
    side = np.array([r["side"] for r in o.R])
    block = np.array([r["block"] for r in o.R])
    dis = np.array([r.get("side_disagrees", "0") == "1" for r in o.R])
    xr = np.array([float(r["x_recomputed_mev"]) if r.get("x_recomputed_mev", "") not in ("", None) else np.nan
                   for r in o.R])
    y = o.within(o.m, o.e, o.e + h)
    out = {"population": o.N, "horizon_years": h, "windows": {}}
    for w in tc.get("windows_mev", [10, 25, 50]):
        inw = np.abs(x) <= w
        A, B = inw & (side == "A"), inw & (side == "B")
        rec = {"N_A": int(A.sum()), "N_B": int(B.sum()), "events_A": int(y[A].sum()), "events_B": int(y[B].sum())}
        try:
            eff = S.blocked_difference(y, A, B, block)
            p = float(y[inw].mean())
            rec.update({"effect": {k: eff[k] for k in ("estimate", "se", "ci_low", "ci_high")},
                        "blocks_used": eff["blocks_used"], "blocks_dropped": eff["blocks_dropped"],
                        "outcome_rate_in_window": p,
                        "mde": S.minimum_detectable_effect(p, int(A.sum()), int(B.sum()))})
            if tc.get("effect_of_interest") is not None:
                rec["mde_exceeds_effect_of_interest"] = rec["mde"] > float(tc["effect_of_interest"])
            keep = ~dis
            e2 = S.blocked_difference(y, A & keep, B & keep, block)
            rec["sensitivity_excluding_label_disagreements"] = {k: e2[k] for k in ("estimate", "ci_low", "ci_high")}
            inr = np.isfinite(xr) & (np.abs(xr) <= w)
            if inr.any():
                e3 = S.blocked_difference(y, inr & (side == "A"), inr & (side == "B"), block)
                rec["sensitivity_recomputed_running_variable"] = {k: e3[k] for k in ("estimate", "ci_low", "ci_high")}
        except ValueError as exc:
            rec["note"] = str(exc)
        out["windows"][str(w)] = rec
    out["balance_primary_window"] = K.balance(o.R, np.abs(x) <= 25, side, block)
    nl, nr = int(((x > -5) & (x <= 0)).sum()), int(((x > 0) & (x <= 5)).sum())
    out["density_test"] = {"n_left": nl, "n_right": nr, "p_value": S.binomial_two_sided(nl, nl + nr)}
    pl = {}
    for c in (-15.0, 15.0):
        left, right = (x > c - 10) & (x <= c), (x > c) & (x <= c + 10)
        try:
            e = S.blocked_difference(y, left, right, block)
            pl[str(c)] = {k: e[k] for k in ("estimate", "ci_low", "ci_high")}
        except ValueError as exc:
            pl[str(c)] = {"note": str(exc)}
    out["placebo"] = pl
    bins = []
    for lo in range(-50, 50, 5):
        dom = (x > lo) & (x <= lo + 5)
        if dom.sum() == 0:
            continue
        sh = S.share(int(y[dom].sum()), int(dom.sum()))
        bins.append({"center": lo + 2.5, "est": sh["est"], "lo": sh["lo"], "hi": sh["hi"], "n": sh["n"], "k": sh["k"]})
    out["bins_5mev"] = bins
    return out


def _release_label(cfg, release_id):
    for r in cfg.get("releases", []):
        if r["release_id"] == release_id:
            return r.get("label", release_id)
    return release_id


def figure_data(cfg, res: dict, robustness: dict | None = None, checks: dict | None = None,
                source_years: dict | None = None) -> dict:
    """Fig. 2: for each entry release, the already-known share at the longest complete horizon
    (Wilson 95%); realization among predictions not already known, with intervals simultaneous
    across the horizons of each release (Bonferroni); and the already-known predictions split by
    record source. Fig. 3: shares in 5 meV/atom bins."""
    fig2a, fig2b, fig2c = [], [], []
    for key, r in res["base"].items():
        if "::" not in key and any(k.startswith(key + "::") for k in res["base"]):
            continue                                       # plot each release, not the pooled cohort
        hs = [h for h, v in r["horizons"].items() if "realized_share" in v]
        if not hs:
            continue
        if "::" in key:
            label = _release_label(cfg, key.split("::")[1])
        else:                                                # a cohort from a single release
            rels = C.releases(cfg, r["source"])
            label = rels[0].get("label", rels[0]["release_id"]) if len(rels) == 1 else r["source"]
        v = r["horizons"][max(hs, key=float)]
        fig2a.append({"source": label, "horizon": float(max(hs, key=float)),
                      **{part: {q: v[name][q] for q in ("est", "lo", "hi", "k", "n")}
                         for part, name in (("already_known", "pre_existence_within_R"),
                                            ("realized", "realized_share"), ("not_realized", "not_realized_share"))}})
        hmax = max(hs, key=float)                       # realized share at the longest horizon, simultaneous
        rz = r["horizons"][hmax]["realized_share"]       # across this release's horizons (as in the text)
        lo, hi = S.wilson_ci(rz["k"], rz["n"], 0.05 / len(hs))
        fig2a[-1]["realized_simultaneous"] = {"est": rz["est"], "lo": lo, "hi": hi, "k": rz["k"], "n": rz["n"]}
        for h in hs:
            rr = r["horizons"][h]["realized_among_not_known"]
            if rr["n"]:
                lo, hi = S.wilson_ci(rr["k"], rr["n"], 0.05 / len(hs))
                fig2b.append({"source": label, "t": float(h), "est": rr["est"], "lo": lo, "hi": hi,
                              "k": rr["k"], "n": rr["n"]})
        v0 = r["horizons"][min(hs, key=float)]                   # already-known split by record source
        tot, cod, icsd = (v0[q]["k"] for q in ("pre_existence_within_R", "pre_existence_COD_only",
                                                "pre_existence_ICSD_derived_only"))
        both = cod + icsd - tot                                  # a prediction is known if either source records it
        fig2c.append({"source": label, "COD only": cod - both, "both": both, "ICSD-derived only": icsd - both,
                      "total": tot})
    order = [r.get("label", r["release_id"]) for r in C.releases(cfg)]
    fig2a.sort(key=lambda d: order.index(d["source"]) if d["source"] in order else len(order), reverse=True)
    out = {"fig2a": fig2a, "fig2b": fig2b, "fig2c": fig2c}
    rels = []
    for r in C.releases(cfg):
        key = next((k for k in res["base"] if k.endswith("::" + r["release_id"])), None)
        if key is None:
            key = next((k for k, v in res["base"].items() if "::" not in k and v["source"] == r["source"]
                        and len(C.releases(cfg, r["source"])) == 1), None)
        if key is None:
            continue
        hs = [float(h) for h, v in res["base"][key]["horizons"].items() if "realized_share" in v]
        rels.append({"label": r.get("label", r["release_id"]).split(" (")[0], "date": decimal_year(r["date"]),
                     "n": res["base"][key]["N"], "t_max": max(hs) if hs else 0.0})
    names = {"COD": "COD", "MP-ICSD": "Materials Project, ICSD-derived", "OQMD-ICSD": "OQMD, ICSD-derived"}
    srcs = [{"name": names.get(k, k), "last_year": v} for k, v in (source_years or {}).items()]
    out["fig1"] = {"releases": rels, "sources": srcs, "cutoff": decimal_year(C.cutoff(cfg))} if rels else None
    out["fig2"] = {"coverage": _coverage_rows(cfg, res), "realized_by_horizon": _realized_rows(cfg, res)}
    if out["fig1"] is not None and cfg.get("workdir"):
        out["fig1"]["ledger"] = _ledger(cfg, res)
    for t in res["threshold"].values():
        out["fig3"] = {"window": 25, "bins": t["bins_5mev"], **_fig3_parts(cfg, t, checks or {}, robustness)}
        break
    return out


def _fig3_parts(cfg, t: dict, checks: dict, robustness: dict | None) -> dict:
    """Panels a, b and d of Fig. 3: label disagreement between databases by distance to the
    threshold (check a), balance by blocking scheme (check b), and the label effect under the
    primary specification and every sensitivity analysis."""
    labels = {r["release_id"]: r.get("label", r["release_id"]) for r in cfg.get("releases", [])}
    premise = []
    for key, v in (checks.get("premise") or {}).items():
        if isinstance(v, dict) and "bins" in v:
            a = key.split(" vs ")[0]
            premise.append({"pair": f"{labels.get(a, a).split(' (')[0]} vs Materials Project", "bins": v["bins"]})
    covs = ["n_elements", "mean_electronegativity", "max_atomic_number"]
    balance = {"covariates": covs, "series": []}
    ov = t["balance_primary_window"]["overall"] or {}
    balance["series"].append({"label": "No blocks", "smd": {c: ov.get(c) for c in covs}})
    pw = next(iter((checks.get("power") or {}).values()), {})
    names = {"class": "Chemistry class", "class_en": "+ electronegativity",
             "class_en_nel": "+ number of elements (primary)"}
    for sc, lab in names.items():
        sm = pw.get("block_scheme_candidates_25meV", {}).get("schemes", {}).get(sc, {}).get("blocked_smd")
        if sm:
            balance["series"].append({"label": lab, "smd": {c: sm.get(c) for c in covs}})
    rows = []

    def add(label, eff, group, note=""):
        if eff and "estimate" in eff:
            rows.append({"label": label, "group": group, "est": eff["estimate"], "lo": eff["ci_low"],
                         "hi": eff["ci_high"], "note": note})
    w = t["windows"]
    add("Primary (\u00b125 meV/atom)", w.get("25", {}).get("effect"), "primary")
    e10 = w.get("10", {})
    add("Window \u00b110 meV/atom", e10.get("effect"), "window",
        "no events" if e10.get("events_A", 1) + e10.get("events_B", 1) == 0 else "")
    add("Window \u00b150 meV/atom", w.get("50", {}).get("effect"), "window")
    add("Recomputed running variable", w.get("25", {}).get("sensitivity_recomputed_running_variable"), "specification")
    add("Excluding label disagreements", w.get("25", {}).get("sensitivity_excluding_label_disagreements"), "specification")
    if robustness:
        stem = next((k for k in robustness.get("reanalysis", {}) if k.startswith("threshold__")), None)
        rr = robustness["reanalysis"].get(stem, {}) if stem else {}
        for key, lab in (("blocks_without_number_of_elements", "Blocks without number of elements"),
                         ("one_entry_per_chemical_system", "One entry per chemical system"),
                         ("labels_agree_with_materials_project", "Labels agree with Materials Project"),
                         ("heusler_type_only", "Heusler-type prototypes only"),
                         ("other_prototypes_only", "Other prototypes only"),
                         ("COD_only", "COD records only"), ("ICSD_derived_only", "ICSD-derived records only")):
            add(lab, rr.get(key, {}).get("effect"), "specification")
        rm = robustness.get("rematch", {}).get("year_only_dates_31_december", {}).get(stem, {})
        add("Year-only dates on 31 December", rm.get("effect"), "specification")
        add("Local-linear regression", rr.get("local_linear_25meV"), "specification")
    for c, lab in (("-15.0", "Placebo threshold at \u221215 meV/atom"), ("15.0", "Placebo threshold at +15 meV/atom")):
        add(lab, t["placebo"].get(c), "placebo")
    w25 = w.get("25", {})
    sides = []
    for lab, k, nn in (("Labeled stable", "events_A", "N_A"), ("Labeled unstable", "events_B", "N_B")):
        if w25.get(nn):
            lo, hi = S.wilson_ci(w25[k], w25[nn])
            sides.append({"label": lab, "k": w25[k], "n": w25[nn], "est": w25[k] / w25[nn], "lo": lo, "hi": hi})
    windows = {}
    pw_all = next(iter((checks.get("power") or {}).values()), {})
    for w_ in ("10", "25", "50"):
        ovw = ((pw_all.get(w_) or {}).get("balance") or {}).get("overall") or {}
        byw = (pw_all.get("block_scheme_balance_by_window") or {}).get(w_) or {}
        ser = [{"label": "No blocks", "smd": {c: ovw.get(c) for c in covs}}]
        for sc, lab in names.items():
            if sc in byw:
                ser.append({"label": lab, "smd": {c: byw[sc]["blocked_smd"].get(c) for c in covs}})
        if ovw:
            windows[w_] = ser
    balance["windows"] = windows
    return {"premise": premise, "balance": balance, "forest": rows, "sides": sides, "mde": w25.get("mde"),
            "events": {k: w25.get(k) for k in ("events_A", "N_A", "events_B", "N_B")}}


def step_estimate(cfg: dict) -> dict:
    wd = cfg["workdir"]
    res = {"base": {}, "threshold": {}, "capture_recapture": {}, "status": {}}
    for cf in sorted(glob.glob(os.path.join(wd, "cohorts", "*.csv"))):
        stem = os.path.splitext(os.path.basename(cf))[0]
        mp = os.path.join(wd, "matches", stem + ".csv")
        if not os.path.exists(mp):
            res["status"][stem] = "no match file; run 'match' first"
            continue
        o = Outcomes(cfg, cf, mp)
        if stem.startswith("threshold__"):
            res["threshold"][stem] = threshold_analysis(cfg, o)
        else:
            src = o.R[0]["source"] if o.R else stem
            res["base"][stem] = cohort_rates(cfg, o, src)
            res["capture_recapture"][stem] = capture_recapture(o)
            rels = sorted({r.get("entry_release", "") for r in o.R})
            if len(rels) > 1:
                for rel in rels:
                    so = Outcomes(cfg, cf, mp, row_filter=lambda r, rel=rel: r.get("entry_release") == rel)
                    if so.N >= 2:
                        res["base"][f"{stem}::{rel}"] = cohort_rates(cfg, so, f"{src} {rel}")
                        res["capture_recapture"][f"{stem}::{rel}"] = capture_recapture(so)
        res["status"][stem] = "estimated"
    res["checks"] = {}
    for name in ("reconstruction", "power", "premise", "coverage_known"):
        p = os.path.join(wd, "checks", f"{name}.json")
        if os.path.exists(p):
            res["checks"][name] = read_json(p)
    res["versions"] = versions()
    write_json(os.path.join(wd, "results", "summary.json"), res)
    log(f"estimate: wrote {os.path.join(wd, 'results', 'summary.json')}")
    return res


def _label_for(cfg, res, key):
    if "::" in key:
        return _release_label(cfg, key.split("::")[1])
    rels = C.releases(cfg, res["base"][key]["source"]) if key in res.get("base", {}) else []
    return rels[0].get("label", rels[0]["release_id"]) if len(rels) == 1 else key


def _coverage_rows(cfg, res) -> list:
    """Capture-recapture coverage of each source and of their union, for records dated before
    2017, with an interval from the Chapman standard error; the estimated total is never below
    the observed union."""
    rows = []
    for key, c in res.get("capture_recapture", {}).items():
        if "chapman_total" not in c:
            continue
        if "::" not in key and len(C.releases(cfg, res["base"][key]["source"])) > 1:
            continue                                     # pooled cohort; its releases appear separately
        nh, se, uni = c["chapman_total"], c["chapman_se"], c["union"]
        lo_n, hi_n = max(nh - 1.96 * se, uni), nh + 1.96 * se
        def share(k):
            return {"est": k / nh, "lo": k / hi_n, "hi": min(1.0, k / lo_n), "k": k}
        rows.append({"source": _label_for(cfg, res, key), "COD": share(c["n_COD"]),
                     "ICSD-derived": share(c["n_ICSD_derived"]), "either": share(uni),
                     "chapman_total": nh, "chapman_se": se, "union": uni})
    return rows


def _realized_rows(cfg, res) -> list:
    """Share of all predictions with a first record within t years, for every complete horizon,
    with Wilson intervals simultaneous across the horizons of each release (Bonferroni)."""
    rows = []
    for key, r in res["base"].items():
        if "::" not in key and len(C.releases(cfg, r["source"])) > 1:
            continue
        hs = sorted((h for h, v in r["horizons"].items() if "realized_share" in v), key=float)
        for h in hs:
            v = r["horizons"][h]["realized_share"]
            lo, hi = S.wilson_ci(v["k"], v["n"], 0.05 / len(hs))
            rows.append({"source": _label_for(cfg, res, key), "t": float(h), "k": v["k"], "n": v["n"],
                         "est": v["est"], "lo": lo, "hi": hi})
    return rows


def _ledger(cfg, res) -> dict:
    """Per-release points for Figs 1 and 2: the date of the earliest matching record of every
    already-known prediction, and the time to the first record of every prediction recorded
    later within the release's longest complete horizon. Counts are checked against the
    summary, so the figures cannot drift from the estimates."""
    import csv
    wd = cfg["workdir"]
    labels = {r["release_id"]: r.get("label", r["release_id"]) for r in C.releases(cfg)}
    out = {}
    for key, r in res["base"].items():
        if "::" not in key and len(C.releases(cfg, r["source"])) > 1:
            continue                                     # pooled cohort; its releases appear separately
        stem = key.split("::")[0]
        cp, mp = (os.path.join(wd, d, f"{stem}.csv") for d in ("cohorts", "matches"))
        if not (os.path.exists(cp) and os.path.exists(mp)):
            continue
        hs = [h for h, v in r["horizons"].items() if "realized_share" in v]
        if not hs:
            continue
        hmax = max(hs, key=float)
        rel = key.split("::")[1] if "::" in key else None
        with open(cp, newline="") as f:
            coh = {row["pred_id"]: row for row in csv.DictReader(f)}
        with open(mp, newline="") as f:
            rows = [(coh[m["pred_id"]], m["db_first_T3"]) for m in csv.DictReader(f) if m["pred_id"] in coh]
        if rel is None:
            rel = rows[0][0]["entry_release"] if rows else None
        rows = [(c, d) for c, d in rows if c["entry_release"] == rel]
        known, later = [], []
        for c, d in rows:
            if not d:
                continue
            e, x = decimal_year(parse_date(c["entry_date"])[0]), decimal_year(parse_date(d)[0])
            if d <= c["entry_date"]:
                known.append(round(x, 4))
            elif x <= e + float(hmax):
                later.append(round(x - e, 4))
        v = r["horizons"][hmax]
        assert len(rows) == r["N"], f"{key}: ledger size {len(rows)} != {r['N']}"
        assert len(known) == v["pre_existence_within_R"]["k"], f"{key}: already known differs"
        assert len(later) == v["realized_share"]["k"], f"{key}: recorded later differs"
        e0 = decimal_year(parse_date(rows[0][0]["entry_date"])[0])
        ages = sorted(e0 - x for x in known)
        q = (lambda a, p: float(np.quantile(a, p))) if ages else (lambda a, p: None)
        out[_release_label(cfg, rel) if rel in labels else key] = {
            "date": e0, "n": r["N"], "t_max": float(hmax), "known_years": sorted(known), "later_offsets": sorted(later),
            "age_median": q(ages, 0.5), "age_q1": q(ages, 0.25), "age_q3": q(ages, 0.75),
            "earliest_record_year": int(min(known)) if known else None}
    return out


def _source_last_years(cfg, min_records: int = 10) -> dict:
    """Last publication year with at least min_records dated records, for each outcome source."""
    from collections import Counter
    from .matching import read_experiments
    p = cfg.get("experiments", {}).get("path")
    if not p or not os.path.exists(C.resolve(cfg, p)):
        return {}
    counts = Counter((e["database"], e["date"].year) for e in read_experiments(C.resolve(cfg, p)) if e["date"])
    out = {}
    for (db, y), n in counts.items():
        if n >= min_records:
            out[db] = max(out.get(db, y), y)
    order = ["COD", "MP-ICSD", "OQMD-ICSD"]
    return {k: out[k] for k in sorted(out, key=lambda k: order.index(k) if k in order else len(order))}


def step_figures(cfg: dict) -> dict:
    """Figures 1 to 3 and Supplementary Figs 1 and 2 from the exported results only
    (results/summary.json, robustness.json and the checks, plus the last dated year of each outcome
    source), with the data behind them, one Source Data workbook per figure and the presentation
    checks."""
    wd = cfg["workdir"]
    res = read_json(os.path.join(wd, "results", "summary.json"))
    rp = os.path.join(wd, "results", "robustness.json")
    rob = read_json(rp) if os.path.exists(rp) else None
    figs = figure_data(cfg, res, rob, res.get("checks"), _source_last_years(cfg))
    write_json(os.path.join(wd, "results", "figure_data.json"), figs)
    from . import source_data
    # The literature audit belongs to the study's own data (beside its workdir), so a synthetic
    # world or any other run without the file gets none; the study configuration is unchanged.
    audit = source_data.read_audit(os.path.join(os.path.dirname(os.path.abspath(wd)),
                                                "data", "audit", "literature_audit.csv"))
    checks = FIG.render_all(figs, os.path.join(wd, "figures"), audit)
    write_json(os.path.join(wd, "results", "figure_checks.json"), checks)
    sd = source_data.write_all(figs, os.path.join(wd, "figures", "source_data"), audit)
    log(f"figures: wrote {os.path.join(wd, 'figures')} and {len(sd)} Source Data workbooks; checks {checks}")
    return checks
