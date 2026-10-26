"""Robustness and sensitivity analyses (Methods, Supplementary Note 5).

Rematch variants rebuild the matches against a modified outcome record: the composition
tolerance halved and doubled, high-pressure reports removed, and year-only dates set to
1 January or 31 December. Reanalyses reuse the primary matches: blocks without the number of
elements, one entry per chemical system, labels that agree with the Materials Project,
Heusler-type and other prototypes, each outcome source alone, and a local-linear regression.
"""
from __future__ import annotations

import datetime as dt
import glob
import os

import numpy as np

from . import config as C
from . import stats as S
from .analysis import ICSD_DERIVED, Outcomes, cohort_rates, threshold_analysis
from .checks import block_keys, en_quintile_edges
from .matching import build_exp_index, read_experiments, step_match
from .snapshots import read_snapshot
from .util import log, read_csv, write_json


def _cohort_files(cfg):
    return sorted(glob.glob(os.path.join(cfg["workdir"], "cohorts", "*.csv")))


def _summary(cfg, cf, mp, sources=None, row_filter=None, rows_edit=None) -> dict:
    stem = os.path.splitext(os.path.basename(cf))[0]
    if stem.startswith("threshold__"):
        o = Outcomes(cfg, cf, mp, sources, row_filter)
        if rows_edit:
            rows_edit(o)
        if o.N < 10:
            return {"note": "too few entries"}
        w = threshold_analysis(cfg, o)["windows"].get("25", {})
        out = {k: w.get(k) for k in ("events_A", "N_A", "events_B", "N_B")}
        if "effect" in w:
            out["effect"] = {k: w["effect"][k] for k in ("estimate", "ci_low", "ci_high")}
        else:
            out["note"] = w.get("note", "not estimable")
        return out
    out = {}
    for rel in sorted({r.get("entry_release", "") for r in read_csv(cf)}):
        flt = (lambda r, rel=rel: r.get("entry_release") == rel and (row_filter is None or row_filter(r)))
        o = Outcomes(cfg, cf, mp, sources, flt)
        if o.N < 2:
            continue
        r = cohort_rates(cfg, o, o.R[0]["source"])
        rec = {"N": o.N, "pre_existence": r["pre_existence_all"]}
        for h, v in r["horizons"].items():
            if "realized_share" in v:
                rec[f"realized_{h}y"] = v["realized_share"]
        out[f"{stem}::{rel}"] = rec
    return out


def _mp_labels(cfg) -> dict:
    """Displayed stable label of each composition in the Materials Project release
    (lowest-energy entry per composition)."""
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


def step_robustness(cfg: dict) -> dict:
    wd = cfg["workdir"]
    cfs = _cohort_files(cfg)
    base_exps = read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
    tols0 = dict(cfg.get("identity_tolerances") or {})
    d0 = float(tols0.get("comp_delta", 0.02))
    res = {"primary": {}, "rematch": {}, "reanalysis": {}}
    for cf in cfs:
        stem = os.path.splitext(os.path.basename(cf))[0]
        res["primary"][stem] = _summary(cfg, cf, os.path.join(wd, "matches", stem + ".csv"))

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
            mp = step_match(cfg, cf, out_path=os.path.join(vd, stem + ".csv"), exps=exps, index=index, tols=tols)
            res["rematch"][name][stem] = _summary(cfg, cf, mp)
    seed = int(cfg.get("seeds", {}).get("one_entry_per_chemical_system", 20261211))
    for cf in cfs:
        stem = os.path.splitext(os.path.basename(cf))[0]
        mp = os.path.join(wd, "matches", stem + ".csv")
        heusler = lambda r: "heusler" in (r.get("prototype") or "").lower()  # noqa: E731
        if not stem.startswith("threshold__"):
            res["reanalysis"][stem] = {"heusler_type_only": _summary(cfg, cf, mp, row_filter=heusler),
                                       "other_prototypes_only": _summary(cfg, cf, mp, row_filter=lambda r: not heusler(r))}
            continue
        rows = read_csv(cf)
        edges = en_quintile_edges(rows)

        def blocks_class_en(o):
            for r, k in zip(o.R, block_keys(o.R, "class_en", edges)):
                r["block"] = k
        rng = np.random.default_rng(seed)
        by_sys = {}
        for r in rows:
            by_sys.setdefault(r["chemsys"], []).append(r["pred_id"])
        one_per = {v[int(rng.integers(len(v)))] for v in by_sys.values()}
        mp_label = _mp_labels(cfg)
        agree = (lambda r: r["reduced_formula"] in mp_label and mp_label[r["reduced_formula"]] == (r["side"] == "A"))
        o = Outcomes(cfg, cf, mp)
        x = np.array([float(r["x_mev"]) for r in o.R])
        y = o.within(o.m, o.e, o.e + float(cfg.get("threshold", {}).get("horizon_years", 5)))
        res["reanalysis"][stem] = {
            "blocks_without_number_of_elements": _summary(cfg, cf, mp, rows_edit=blocks_class_en),
            "one_entry_per_chemical_system": _summary(cfg, cf, mp, row_filter=lambda r: r["pred_id"] in one_per),
            "labels_agree_with_materials_project": (_summary(cfg, cf, mp, row_filter=agree) if mp_label
                                                    else {"note": "no Materials Project release"}),
            "heusler_type_only": _summary(cfg, cf, mp, row_filter=heusler),
            "other_prototypes_only": _summary(cfg, cf, mp, row_filter=lambda r: not heusler(r)),
            "COD_only": _summary(cfg, cf, mp, sources=("COD",)),
            "ICSD_derived_only": _summary(cfg, cf, mp, sources=ICSD_DERIVED),
            "local_linear_25meV": S.local_linear_jump(x, y, 25.0)}
    write_json(os.path.join(wd, "results", "robustness.json"), res)
    log(f"robustness: wrote {os.path.join(wd, 'results', 'robustness.json')}")
    return res
