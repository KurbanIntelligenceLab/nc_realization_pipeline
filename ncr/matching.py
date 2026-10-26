"""Deterministic matching of every prediction against dated experimental-structure records.
For each prediction the earliest matching record is kept at each identity tier, overall and
for each source separately (COD, MP-ICSD, OQMD-ICSD)."""
from __future__ import annotations

import os
from collections import defaultdict

from . import config as C
from . import identity
from .util import date_str, log, parse_date, read_csv, truthy, write_csv


def _comp(formula):
    from pymatgen.core import Composition
    return Composition(formula)


def read_experiments(path: str) -> list[dict]:
    rows = read_csv(path)
    out = []
    for i, r in enumerate(rows):
        for k in ("exp_id", "database", "formula", "first_report_date", "experimental"):
            if r.get(k, "") == "":
                raise ValueError(f"{path} row {i + 2}: '{k}' is empty")
        if not truthy(r["experimental"]):
            continue                          # theoretical records are excluded (Supplementary Table 1)
        d, _ = parse_date(r["first_report_date"])
        sf = (r.get("structure_file") or "").strip()
        if sf and not os.path.isabs(sf):
            sf = os.path.join(os.path.dirname(os.path.abspath(path)), sf)
        out.append({"exp_id": r["exp_id"], "database": r["database"], "formula": r["formula"],
                    "date": d, "structure_file": sf,
                    "ss_end_a": (r.get("ss_end_a") or "").strip(), "ss_end_b": (r.get("ss_end_b") or "").strip(),
                    "high_pressure": truthy(r.get("high_pressure", "")),
                    "thin_film_only": truthy(r.get("thin_film_only", ""))})
    return out


_struct_cache: dict = {}


def _structure(path):
    if not path:
        return None
    if path not in _struct_cache:
        from pymatgen.core import Structure
        _struct_cache[path] = Structure.from_file(path)
    return _struct_cache[path]


def match_tier(pred_formula, pred_struct_file, exp, tols=None):
    """Tier of one experimental report for one prediction (Table 2 rules)."""
    pred = _structure(pred_struct_file) or pred_formula
    ex = _structure(exp["structure_file"]) or exp["formula"]
    tier = identity.outcome_tier(pred, ex, tols)
    if tier is None and exp["ss_end_a"] and exp["ss_end_b"]:
        if identity.in_solid_solution_range(pred_formula, exp["ss_end_a"], exp["ss_end_b"],
                                       (tols or {}).get("comp_delta", identity.DEFAULT_TOLS["comp_delta"])):
            tier = "T2" if pred_struct_file else "T3"
    return tier


def build_exp_index(exps: list[dict]) -> dict:
    by_sys = defaultdict(list)
    for e in exps:
        by_sys["-".join(sorted(str(x) for x in _comp(e["formula"]).elements))].append(e)
    return by_sys


def match_rows(rows: list[dict], by_sys: dict, sources: list[str], tols=None) -> list[dict]:
    out = []
    for r in rows:
        best = {"T1": None, "T2": None, "T3": None}
        by_src = {}
        flags = set()
        for e in by_sys.get(r["chemsys"], []):
            tier = match_tier(r["formula"], r.get("structure_file", ""), e, tols)
            if tier is None:
                continue
            for t in ("T1", "T2", "T3"):
                if identity.meets(tier, t) and e["date"] and (best[t] is None or e["date"] < best[t][0]):
                    best[t] = (e["date"], e["exp_id"])
                if identity.meets(tier, t) and e["date"]:
                    k = (t, e["database"])
                    if k not in by_src or e["date"] < by_src[k]:
                        by_src[k] = e["date"]
            if e["high_pressure"]:
                flags.add("high_pressure")
            if e["thin_film_only"]:
                flags.add("thin_film_only")
        out.append({"pred_id": r["pred_id"],
                    **{f"db_first_{t}": date_str(best[t][0]) if best[t] else "" for t in best},
                    **{f"db_id_{t}": best[t][1] if best[t] else "" for t in best},
                    **{f"db_first_{t}__{db}": date_str(by_src.get((t, db))) for t in ("T2", "T3") for db in sources},
                    "db_flags": " ".join(sorted(flags))})
    return out


def step_match(cfg: dict, cohort_file: str, out_path: str | None = None, exps: list[dict] | None = None,
                 index: dict | None = None, tols: dict | None = None) -> str:
    if exps is None:
        exps = read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
    by_sys = index if index is not None else build_exp_index(exps)
    sources = sorted({e["database"] for e in exps})
    out = match_rows(read_csv(cohort_file), by_sys, sources, tols if tols is not None else cfg.get("identity_tolerances"))
    path = out_path or cohort_file.replace(os.sep + "cohorts" + os.sep, os.sep + "matches" + os.sep)
    write_csv(path, out)
    log(f"match: {os.path.basename(cohort_file)}: {sum(1 for o in out if o['db_first_T3'])} predictions with a matching record")
    return path
