"""Robustness and sensitivity analyses specified in Methods and Supplementary Note 5.

Each variant either rematches every cohort against a modified outcome record (composition
tolerance halved or doubled; high-pressure reports removed; year-only dates set to 1 January
or 31 December) or reanalyses the primary matches (blocks without the number of elements; one
entry per chemical system; labels that agree with the Materials Project; Heusler-type and other
prototypes; each outcome source alone; local-linear regression). Results go to
work/robustness_for_claude.json (aggregates only)."""
from __future__ import annotations

import copy
import datetime as dt
import glob
import math
import os

import numpy as np

from . import config as C
from .census import ICSD_DERIVED, CensusBundle, census_rates
from .cohorts import build_exp_index, read_experiments, step_dbmatch
from .estimate import block_keys, en_quintile_edges, threshold
from .snapshots import read_snapshot
from .util import log, read_csv, write_json


def _cohort_files(cfg):
    return sorted(glob.glob(os.path.join(cfg["workdir"], "cohorts", "*.csv")))


def _summ(cfg, cf, dbp, sources=None, row_filter=None, rows_edit=None):
    """Primary numbers for one cohort file under one variant."""
    stem = os.path.splitext(os.path.basename(cf))[0]
    out = {}
    if stem.startswith("threshold__"):
        b = CensusBundle(cfg, cf, dbp, sources, row_filter)
        if rows_edit:
            rows_edit(b)
        if b.N < 10:
            return {"note": "too few entries"}
        t = threshold(cfg, b)["windows"].get("25", {})
        x = np.array([float(r["x_mev"]) for r in b.R])
        side = np.array([r["side"] for r in b.R])
        Y = b.within(b.m, b.e, b.e + float(cfg.get("threshold", {}).get("horizon_years", 5)))
        inw = np.abs(x) <= 25
        out = {"events_A": int(Y[inw & (side == "A")].sum()), "N_A": int((inw & (side == "A")).sum()),
               "events_B": int(Y[inw & (side == "B")].sum()), "N_B": int((inw & (side == "B")).sum())}
        if "effect" in t:
            out["effect"] = {k: t["effect"][k] for k in ("estimate", "ci_low", "ci_high")}
        else:
            out["note"] = t.get("note", "not estimable")
        return out
    rels = sorted({r.get("entry_release", "") for r in read_csv(cf)})
    for rel in rels:
        flt = (lambda r, rel=rel: r.get("entry_release") == rel and (row_filter is None or row_filter(r)))
        b = CensusBundle(cfg, cf, dbp, sources, flt)
        if b.N < 2:
            continue
        r = census_rates(cfg, b, b.R[0]["source"])
        rec = {"N": b.N, "pre_existence": r["pre_existence_all"]}
        for h, v in r["horizons"].items():
            if "realized_share" in v:
                rec[f"realized_{h}y"] = v["realized_share"]
        out[f"{stem}::{rel}"] = rec
    return out


def _local_linear(cfg, cf, dbp, bw=25.0):
    """Sharp local-linear estimate of the jump at the threshold (triangular kernel, HC1 errors)."""
    b = CensusBundle(cfg, cf, dbp)
    x = np.array([float(r["x_mev"]) for r in b.R])
    Y = b.within(b.m, b.e, b.e + float(cfg.get("threshold", {}).get("horizon_years", 5)))
    m = np.abs(x) <= bw
    x, Y = x[m], Y[m]
    D = (x <= 0).astype(float)
    X = np.column_stack([np.ones_like(x), x, D, D * x])
    w = 1 - np.abs(x) / bw
    W = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * W[:, None], Y * W, rcond=None)
    res = Y - X @ beta
    XtWX = (X * w[:, None]).T @ X
    inv = np.linalg.pinv(XtWX)
    meat = (X * (w * res)[:, None]).T @ (X * (w * res)[:, None])
    n, k = X.shape
    V = inv @ meat @ inv * n / max(n - k, 1)
    se = float(math.sqrt(max(V[2, 2], 0.0)))
    return {"bandwidth_mev": bw, "n": int(n), "events": int(Y.sum()), "estimate": float(beta[2]), "se": se,
            "ci_low": float(beta[2] - 1.96 * se), "ci_high": float(beta[2] + 1.96 * se)}


def step_robustness(cfg: dict) -> dict:
    wd = cfg["workdir"]
    cfs = _cohort_files(cfg)
    base_exps = read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
    tols0 = dict(cfg.get("identity_tolerances") or {})
    d0 = float(tols0.get("comp_delta", 0.02))
    res = {"primary": {}, "rematch": {}, "reanalysis": {}}
    for cf in cfs:
        stem = os.path.splitext(os.path.basename(cf))[0]
        res["primary"][stem] = _summ(cfg, cf, os.path.join(wd, "dbmatch", stem + ".csv"))

    def year_date(e, day):
        if e["date"] is not None and e["date"].month == 7 and e["date"].day == 1:
            e = dict(e)
            e["date"] = dt.date(e["date"].year, *day)
        return e
    variants = {
        "composition_tolerance_halved": (base_exps, dict(tols0, comp_delta=d0 / 2)),
        "composition_tolerance_doubled": (base_exps, dict(tols0, comp_delta=d0 * 2)),
        "without_high_pressure_reports": ([e for e in base_exps if not e["high_pressure"]], tols0),
        "year_only_dates_1_january": ([year_date(e, (1, 1)) for e in base_exps], tols0),
        "year_only_dates_31_december": ([year_date(e, (12, 31)) for e in base_exps], tols0),
    }
    for name, (exps, tols) in variants.items():
        log(f"robustness: rematching for {name}")
        index = build_exp_index(exps)
        vd = os.path.join(wd, "robust", name)
        os.makedirs(vd, exist_ok=True)
        res["rematch"][name] = {}
        for cf in cfs:
            stem = os.path.splitext(os.path.basename(cf))[0]
            dbp = step_dbmatch(cfg, cf, out_path=os.path.join(vd, stem + ".csv"), exps=exps, index=index, tols=tols)
            res["rematch"][name][stem] = _summ(cfg, cf, dbp)
    # Reanalyses on the primary matches
    thr = [cf for cf in cfs if os.path.basename(cf).startswith("threshold__")]
    for cf in thr:
        stem = os.path.splitext(os.path.basename(cf))[0]
        dbp = os.path.join(wd, "dbmatch", stem + ".csv")
        rows = read_csv(cf)
        edges = en_quintile_edges(rows)

        def blocks_class_en(b):
            for r, k in zip(b.R, block_keys(b.R, "class_en", edges)):
                r["block"] = k
        rng = np.random.default_rng(int(cfg.get("seeds", {}).get("expert_sample", 1)))
        keep = {}
        for r in rows:
            keep.setdefault(r["chemsys"], []).append(r["pred_id"])
        one_per = {v[int(rng.integers(len(v)))] for v in keep.values()}
        mp_label = _mp_labels(cfg)
        agree = (lambda r: r["reduced_formula"] in mp_label and mp_label[r["reduced_formula"]] == (r["side"] == "A"))
        rr = {"blocks_without_number_of_elements": _summ(cfg, cf, dbp, rows_edit=blocks_class_en),
              "one_entry_per_chemical_system": _summ(cfg, cf, dbp, row_filter=lambda r: r["pred_id"] in one_per),
              "labels_agree_with_materials_project": _summ(cfg, cf, dbp, row_filter=agree) if mp_label else {"note": "no Materials Project release"},
              "heusler_type_only": _summ(cfg, cf, dbp, row_filter=lambda r: "heusler" in (r.get("prototype") or "").lower()),
              "other_prototypes_only": _summ(cfg, cf, dbp, row_filter=lambda r: "heusler" not in (r.get("prototype") or "").lower()),
              "COD_only": _summ(cfg, cf, dbp, sources=("COD",)),
              "ICSD_derived_only": _summ(cfg, cf, dbp, sources=ICSD_DERIVED),
              "local_linear_25meV": _local_linear(cfg, cf, dbp)}
        res["reanalysis"][stem] = rr
    for cf in cfs:
        stem = os.path.splitext(os.path.basename(cf))[0]
        if stem.startswith("threshold__"):
            continue
        dbp = os.path.join(wd, "dbmatch", stem + ".csv")
        res["reanalysis"][stem] = {
            "heusler_type_only": _summ(cfg, cf, dbp, row_filter=lambda r: "heusler" in (r.get("prototype") or "").lower()),
            "other_prototypes_only": _summ(cfg, cf, dbp, row_filter=lambda r: "heusler" not in (r.get("prototype") or "").lower())}
    write_json(os.path.join(wd, "robustness_for_claude.json"), res)
    log(f"robustness: wrote {os.path.join(wd, 'robustness_for_claude.json')}")
    return res


def _mp_labels(cfg):
    """Displayed stable label of each composition in the Materials Project release (lowest-energy entry)."""
    from pymatgen.core import Composition
    rel = [r for r in C.releases(cfg) if r["source"] == "Materials Project"]
    if not rel:
        return {}
    best = {}
    for s in read_snapshot(C.resolve(cfg, rel[0]["path"]), rel[0].get("columns")):
        if s["entry_id"].startswith("ref-") or s["displayed_stable"] is None:
            continue
        k = Composition(s["formula"]).reduced_formula
        if k not in best or s["energy_per_atom"] < best[k][0]:
            best[k] = (s["energy_per_atom"], bool(s["displayed_stable"]))
    return {k: v[1] for k, v in best.items()}
