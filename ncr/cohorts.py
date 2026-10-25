"""Hull recomputation per release, cohort construction, the threshold population,
the relabeling population, and deterministic matching against structure-database
reports (ICSD or COD exports)."""
from __future__ import annotations

import datetime as dt
import math
import os
import sys
from collections import defaultdict

from . import config as C
from .snapshots import compute_dhd, read_snapshot
from .util import date_str, decimal_year, log, parse_date, read_csv, truthy, write_csv

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(HERE))          # bundle root: nc_realization_cohorts.py
sys.path.insert(0, HERE)
import nc_realization_cohorts as NCC                # noqa: E402

WINDOW_MAX_MEV = 50.0


def _comp(formula):
    from pymatgen.core import Composition
    return Composition(formula)


def compound_features(formula: str) -> dict:
    c = _comp(formula)
    els = sorted(str(e) for e in c.elements)
    xs = [e.X for e in c.elements if e.X == e.X]
    return {"reduced_formula": c.reduced_formula, "elements": " ".join(els),
            "chemsys": "-".join(els), "n_elements": len(els),
            "chemistry_class": NCC.chemistry_class(formula),
            "mean_electronegativity": round(sum(xs) / len(xs), 4) if xs else "",
            "max_atomic_number": max(e.Z for e in c.elements)}


# ----------------------------------------------------------------------------
# Step: hull
# ----------------------------------------------------------------------------
def step_hull(cfg: dict, workers: int = 1, all_entries: bool = False) -> None:
    outdir = os.path.join(cfg["workdir"], "hull")
    for r in C.releases(cfg):
        out = os.path.join(outdir, f"{r['release_id']}.csv")
        if os.path.exists(out):
            log(f"hull: {r['release_id']} exists, skipping (delete the file to recompute)")
            continue
        rows = read_snapshot(C.resolve(cfg, r["path"]), r.get("columns"))
        ids = [x["entry_id"] for x in rows]
        if len(set(ids)) != len(ids):
            raise ValueError(f"{r['release_id']}: entry_id values are not unique")
        tau = float(C.need(r, "tau_ev_per_atom"))
        filt = None if all_entries else (lambda x: x["is_hypothetical"] and not x["origin_icsd_id"])
        dhd = compute_dhd(rows, filt, workers)
        res, n_dis = [], 0
        for x in rows:
            if x["entry_id"] not in dhd:
                continue
            d, err = dhd[x["entry_id"]]
            rec_stable = (d <= tau) if not math.isnan(d) else None
            disp = x["displayed_stable"]
            dis = (disp is not None and rec_stable is not None and disp != rec_stable)
            n_dis += int(dis)
            res.append({"entry_id": x["entry_id"], "formula": x["formula"],
                        "is_hypothetical": int(x["is_hypothetical"]),
                        "dhd_recomputed": "" if math.isnan(d) else round(d, 6),
                        "stable_recomputed": "" if rec_stable is None else int(rec_stable),
                        "displayed_stable": "" if disp is None else int(disp),
                        "displayed_dhd": "" if x["displayed_dhd"] is None else x["displayed_dhd"],
                        "label_disagrees": int(dis), "error": err})
        write_csv(out, res)
        log(f"hull: {r['release_id']}: {len(res)} entries, {n_dis} displayed/recomputed label disagreements")


def _hull_table(cfg, release_id):
    p = os.path.join(cfg["workdir"], "hull", f"{release_id}.csv")
    if not os.path.exists(p):
        raise FileNotFoundError(f"run 'hull' first ({p} missing)")
    return {r["entry_id"]: r for r in read_csv(p)}


def _stable(r_hull: dict, side_from: str):
    """Stable label of one entry. 'displayed' is strict: an entry whose release did not
    display a label is not treated as presented stable (None). 'recomputed' uses our hull."""
    if side_from == "displayed":
        return None if r_hull["displayed_stable"] == "" else r_hull["displayed_stable"] == "1"
    if r_hull["stable_recomputed"] == "":
        return None
    return r_hull["stable_recomputed"] == "1"


# ----------------------------------------------------------------------------
# Step: x0 (snapshot reconstruction)
# ----------------------------------------------------------------------------
def step_x0(cfg: dict) -> dict:
    """Snapshot reconstruction: stable hypothetical entries and new compositions under the
    displayed and the recomputed labels, compared with any expected counts in the config
    (relative tolerance x0_tolerance, default 1%)."""
    out = {}
    tol = float(cfg.get("x0_tolerance", 0.01))
    for r in C.releases(cfg):
        h = _hull_table(cfg, r["release_id"])
        snap = {x["entry_id"]: x for x in read_snapshot(C.resolve(cfg, r["path"]), r.get("columns"))}
        rec = {}
        for mode in ("displayed", "recomputed"):
            st = [eid for eid, x in h.items() if x["is_hypothetical"] == "1" and _stable(x, mode)]
            new = [eid for eid in st if not snap[eid]["source_known_composition"] and not snap[eid]["origin_icsd_id"]]
            comps = {_comp(snap[eid]["formula"]).reduced_formula for eid in new}
            heus = {_comp(snap[eid]["formula"]).reduced_formula for eid in new
                    if "heusler" in (snap[eid].get("prototype") or "").lower()}
            rec[mode] = {"stable_hypothetical_entries": len(st), "stable_new_compositions": len(comps),
                         "stable_new_heusler_type_compositions": len(heus)}
        rec["recomputed_label_disagreements"] = sum(int(x["label_disagrees"]) for x in h.values())
        conf = {"disp_S_rec_S": 0, "disp_S_rec_U": 0, "disp_U_rec_S": 0, "disp_U_rec_U": 0}
        cd, cr = set(), set()
        for eid, x in h.items():
            if x["is_hypothetical"] != "1" or x["displayed_stable"] == "" or x["stable_recomputed"] == "":
                continue
            ds, rs = x["displayed_stable"] == "1", x["stable_recomputed"] == "1"
            conf[f"disp_{'S' if ds else 'U'}_rec_{'S' if rs else 'U'}"] += 1
            if not snap[eid]["source_known_composition"] and not snap[eid]["origin_icsd_id"] and (ds or rs):
                red = _comp(snap[eid]["formula"]).reduced_formula
                if ds:
                    cd.add(red)
                if rs:
                    cr.add(red)
        rec["labeled_confusion"] = conf
        rec["recomputed_recall_of_displayed_stable"] = (conf["disp_S_rec_S"] / (conf["disp_S_rec_S"] + conf["disp_S_rec_U"])
                                                         if conf["disp_S_rec_S"] + conf["disp_S_rec_U"] else None)
        rec["recomputed_precision_for_displayed_stable"] = (conf["disp_S_rec_S"] / (conf["disp_S_rec_S"] + conf["disp_U_rec_S"])
                                                            if conf["disp_S_rec_S"] + conf["disp_U_rec_S"] else None)
        rec["new_composition_overlap"] = {"displayed": len(cd), "recomputed": len(cr), "both": len(cd & cr),
                                          "jaccard": len(cd & cr) / len(cd | cr) if (cd | cr) else None}
        rec["hull_errors"] = sum(1 for x in h.values() if x["error"])
        exp = r.get("expected_counts") or {}
        cmp_ = {}
        for k, v in exp.items():
            got = rec["displayed"].get(k)
            if got is not None and v:
                cmp_[k] = {"expected": v, "displayed": got, "recomputed": rec["recomputed"].get(k),
                           "rel_diff_displayed": (got - v) / v,
                           "rel_diff_recomputed": (rec["recomputed"].get(k, 0) - v) / v}
        rec["expected"] = cmp_
        rec["within_tolerance"] = (all(abs(c["rel_diff_displayed"]) <= tol for c in cmp_.values()) if cmp_ else None)
        n_lab = sum(1 for x in h.values() if x["is_hypothetical"] == "1" and x["displayed_stable"] != "" and x["stable_recomputed"] != "")
        rec["recomputed_agrees_with_displayed_share"] = (1 - rec["recomputed_label_disagreements"] / n_lab) if n_lab else None
        out[r["release_id"]] = rec
        log(f"x0: {r['release_id']}: {rec}")
    return out


# ----------------------------------------------------------------------------
# Step: cohorts (prediction cohorts, threshold population, two-phase sample)
# ----------------------------------------------------------------------------
def step_cohorts(cfg: dict) -> None:
    import numpy as np
    cutoff = C.cutoff(cfg)
    h_thr = float(cfg.get("threshold", {}).get("horizon_years", 5))
    cdir = os.path.join(cfg["workdir"], "cohorts")
    for src in C.sources(cfg):
        rels = [r for r in C.releases(cfg, src) if r.get("use_for_cohorts", True)]
        if not rels:
            log(f"cohorts: {src}: no release is used for cohorts (use_for_cohorts false); skipped")
            continue
        gran = rels[0].get("granularity", "structure")
        first_stable, first_seen = {}, {}
        for r in rels:
            h = _hull_table(cfg, r["release_id"])
            snap = {x["entry_id"]: x for x in read_snapshot(C.resolve(cfg, r["path"]), r.get("columns"))}
            tau = float(C.need(r, "tau_ev_per_atom"))
            side_from = r.get("side_from", "displayed")
            per_key = {}
            for eid, x in h.items():
                s = snap[eid]
                if not s["is_hypothetical"] or s["origin_icsd_id"]:
                    continue                                   # cohort hygiene (Table S1)
                if s["source_known_composition"] and r.get("exclude_source_known_compositions", True):
                    continue                                   # the source itself listed this composition as known
                key = eid if gran == "structure" else _comp(s["formula"]).reduced_formula
                base = {"key": key, "entry_id": eid, "formula": s["formula"], "release": r["release_id"],
                        "date": r["date"], "tau": tau, "structure_file": s["structure_file"],
                        "spacegroup_number": s["spacegroup_number"], "prototype": s["prototype"],
                        "deprecated": int(s["deprecated"]),
                        "dhd": float(x["dhd_recomputed"]) if x["dhd_recomputed"] != "" else math.nan,
                        "dhd_displayed": float(x["displayed_dhd"]) if x.get("displayed_dhd", "") != "" else math.nan,
                        "energy": s["energy_per_atom"],
                        "side_displayed": _stable(x, "displayed"), "side_recomputed": _stable(x, "recomputed"),
                        "stable_rule": _stable(x, side_from)}
                # For composition-level keys the composition is represented by its lowest-energy
                # entry in this release (the only one that can be on the hull).
                cur = per_key.get(key)
                if cur is None or base["energy"] < cur["energy"]:
                    per_key[key] = base
            for key, base in per_key.items():
                if key not in first_seen:
                    first_seen[key] = base
                if base["stable_rule"] and key not in first_stable:
                    first_stable[key] = base
        rows = []
        log(f"cohorts: {src}: entries without a displayed label are not treated as presented stable "
            f"({sum(1 for b in first_seen.values() if b['side_displayed'] is None)} keys first seen without one)")
        for key, b in first_stable.items():
            f = compound_features(b["formula"])
            rows.append({"pred_id": f"{src}::{key}", "cohort": src, "source": src, "granularity": gran,
                         "entry_release": b["release"], "entry_date": date_str(b["date"]),
                         "formula": b["formula"], **f, "structure_file": b["structure_file"],
                         "spacegroup_number": b["spacegroup_number"], "prototype": b["prototype"],
                         "dhd_entry": "" if math.isnan(b["dhd"]) else b["dhd"], "deprecated": b["deprecated"]})
        write_csv(os.path.join(cdir, f"{_slug(src)}.csv"), rows)
        log(f"cohorts: {src}: {len(rows)} predictions released as stable")
        # Two-phase sample for very large cohorts (Table S1, "Large cohorts")
        n1 = cfg.get("two_phase", {}).get(src)
        if n1 and len(rows) > int(n1):
            rng = np.random.default_rng(int(C.need(cfg, "seeds.two_phase")))
            idx = sorted(rng.choice(len(rows), size=int(n1), replace=False).tolist())
            write_csv(os.path.join(cdir, f"{_slug(src)}__phase1.csv"), [rows[i] for i in idx])
            log(f"cohorts: {src}: phase-1 random sample of {n1} of {len(rows)} written")
        # Threshold population: all hypothetical entries at first appearance
        if src in cfg.get("threshold", {}).get("sources", []):
            thr, n_unknown, n_rv_fallback = [], 0, 0
            side_rule = rels[0].get("side_from", "displayed")
            rv_rule = cfg.get("threshold", {}).get("running_variable", "displayed")
            for key, b in first_seen.items():
                if side_rule == "displayed" and b["side_displayed"] is None:
                    n_unknown += 1
                    continue
                if decimal_year(b["date"]) + h_thr > decimal_year(cutoff) - 1.0:
                    continue
                # Running variable: the decomposition energy the release itself displayed, so that
                # side and running variable come from the same published quantity (a sharp design);
                # our recomputed value is kept for the sensitivity analysis.
                if rv_rule == "displayed" and not math.isnan(b["dhd_displayed"]):
                    dval = b["dhd_displayed"]
                else:
                    dval = b["dhd"]
                    n_rv_fallback += int(rv_rule == "displayed")
                if math.isnan(dval):
                    continue
                x = (dval - b["tau"]) * 1000.0
                if abs(x) > WINDOW_MAX_MEV:
                    continue
                side = b["side_displayed"] if side_rule == "displayed" else b["side_recomputed"]
                f = compound_features(b["formula"])
                thr.append({"pred_id": f"{src}::{key}", "cohort": f"threshold::{src}", "source": src,
                            "granularity": gran, "entry_release": b["release"], "entry_date": date_str(b["date"]),
                            "formula": b["formula"], **f, "structure_file": b["structure_file"],
                            "spacegroup_number": b["spacegroup_number"], "prototype": b["prototype"],
                            "x_mev": round(x, 3),
                            "x_recomputed_mev": "" if math.isnan(b["dhd"]) else round((b["dhd"] - b["tau"]) * 1000.0, 3),
                            "side": "A" if side else "B",
                            "side_recomputed": "A" if b["side_recomputed"] else "B",
                            "side_disagrees": int(bool(side) != bool(b["side_recomputed"])),
                            "block": f"{src}|{f['chemistry_class']}"})
            scheme = cfg.get("threshold", {}).get("block_scheme", "class")
            if scheme != "class" and thr:
                from .estimate import block_keys, en_quintile_edges
                edges = en_quintile_edges(thr)
                for r, k in zip(thr, block_keys(thr, scheme, edges)):
                    r["block"] = k
                from .util import write_json as _wj
                _wj(os.path.join(cfg["workdir"], "checks", f"threshold_blocks__{_slug(src)}.json"),
                    {"scheme": scheme, "electronegativity_quintile_edges": edges})
            write_csv(os.path.join(cdir, f"threshold__{_slug(src)}.csv"), thr)
            log(f"cohorts: {src}: threshold population of {len(thr)} entries within ±{WINDOW_MAX_MEV:.0f} meV/atom "
                f"({n_unknown} excluded because no label was displayed; running variable '{rv_rule}', "
                f"{n_rv_fallback} entries fell back to the recomputed value)")


def _slug(s: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_").lower()


# ----------------------------------------------------------------------------
# Step: relabel (correction-driven flips between consecutive builds)
# ----------------------------------------------------------------------------
def _apply_scheme(raw_dicts, scheme_name):
    from pymatgen.entries import compatibility as compat
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
    scheme = getattr(compat, scheme_name)()
    ents = [(ComputedStructureEntry if "structure" in d else ComputedEntry).from_dict(d) for d in raw_dicts]
    kept = scheme.process_entries(ents, clean=True, inplace=False, on_error="ignore")
    return {str(e.entry_id): (e.composition.formula, e.energy) for e in kept}


def step_relabel(cfg: dict) -> None:
    import json
    for spec in cfg.get("relabeling", []):
        name = spec["name"]
        r_id, r1_id = spec["release_r"], spec["release_r1"]
        rel = {r["release_id"]: r for r in C.releases(cfg)}
        hr, hr1 = _hull_table(cfg, r_id), _hull_table(cfg, r1_id)
        snap1 = {x["entry_id"]: x for x in read_snapshot(C.resolve(cfg, rel[r1_id]["path"]), rel[r1_id].get("columns"))}
        with open(C.resolve(cfg, spec["raw_entries_r1"]), encoding="utf-8") as fh:
            raw = json.load(fh)
        raw = raw["entries"] if isinstance(raw, dict) else raw
        e_j = _apply_scheme(raw, C.need(spec, "scheme_r"))
        e_k = _apply_scheme(raw, C.need(spec, "scheme_r1"))
        common_ids = set(e_j) & set(e_k)
        by_sys = defaultdict(list)
        for i in common_ids:
            by_sys["-".join(sorted(str(x) for x in _comp(e_k[i][0]).elements))].append(i)
        species = list(C.need(spec, "corrected_species"))
        cutoff = C.cutoff(cfg)
        h = float(cfg.get("threshold", {}).get("horizon_years", 5))
        if decimal_year(rel[r1_id]["date"]) + h > decimal_year(cutoff) - 1.0:
            raise ValueError(f"{name}: release {r1_id} + {h} yr is later than the cutoff minus one year")
        rows = []
        for eid in sorted(set(hr) & set(hr1)):
            a, b = hr[eid], hr1[eid]
            if a["is_hypothetical"] != "1" or snap1.get(eid, {}).get("origin_icsd_id"):
                continue
            s0, s1 = _stable(a, "displayed"), _stable(b, "displayed")
            if s0 is None or s1 is None or a["dhd_recomputed"] == "":
                continue
            if (not s0) and s1:
                els = sorted(str(x) for x in _comp(b["formula"]).elements)
                pool = []
                for k in range(1, len(els) + 1):
                    from itertools import combinations
                    for sub in combinations(els, k):
                        pool.extend(by_sys.get("-".join(sub), []))
                if eid not in common_ids:
                    group = "flip_unclassifiable"
                else:
                    ej = {i: e_j[i] for i in pool}
                    ek = {i: e_k[i] for i in pool}
                    group = "F" if NCC.flip_is_correction_driven(ej, ek, eid) else "flip_other"
            elif (not s0) and (not s1):
                group = "M"
            else:
                continue
            f = compound_features(b["formula"])
            dbin = math.floor(float(a["dhd_recomputed"]) * 1000.0 / 5.0)
            pres = "".join("1" if sp in f["elements"].split() else "0" for sp in species)
            rows.append({"pred_id": f"{rel[r1_id]['source']}::{eid}", "cohort": f"relabel::{name}",
                         "source": rel[r1_id]["source"], "granularity": "structure",
                         "entry_release": r1_id, "entry_date": date_str(rel[r1_id]["date"]),
                         "release_r_date": date_str(rel[r_id]["date"]),
                         "formula": b["formula"], **f, "group": group,
                         "structure_file": snap1.get(eid, {}).get("structure_file", ""),
                         "spacegroup_number": snap1.get(eid, {}).get("spacegroup_number", ""),
                         "cem_key": f"{f['n_elements']}|{f['chemistry_class']}|{dbin}|{pres}"})
        write_csv(os.path.join(cfg["workdir"], "cohorts", f"relabel__{_slug(name)}.csv"), rows)
        cnt = defaultdict(int)
        for r in rows:
            cnt[r["group"]] += 1
        log(f"relabel: {name}: {dict(cnt)}")


# ----------------------------------------------------------------------------
# Experimental reports and deterministic database matching
# ----------------------------------------------------------------------------
def read_experiments(path: str) -> list[dict]:
    rows = read_csv(path)
    out = []
    for i, r in enumerate(rows):
        for k in ("exp_id", "database", "formula", "first_report_date", "experimental"):
            if r.get(k, "") == "":
                raise ValueError(f"{path} row {i + 2}: '{k}' is empty")
        if not truthy(r["experimental"]):
            continue                          # theoretical structures are excluded (Table S1)
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
    tier = NCC.outcome_tier(pred, ex, tols)
    if tier is None and exp["ss_end_a"] and exp["ss_end_b"]:
        if NCC.in_solid_solution_range(pred_formula, exp["ss_end_a"], exp["ss_end_b"],
                                       (tols or {}).get("comp_delta", NCC.DEFAULT_TOLS["comp_delta"])):
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
                if NCC.meets(tier, t) and e["date"] and (best[t] is None or e["date"] < best[t][0]):
                    best[t] = (e["date"], e["exp_id"])
                if NCC.meets(tier, t) and e["date"]:
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


def step_dbmatch(cfg: dict, cohort_file: str, out_path: str | None = None, exps: list[dict] | None = None,
                 index: dict | None = None, tols: dict | None = None) -> str:
    if exps is None:
        exps = read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
    by_sys = index if index is not None else build_exp_index(exps)
    sources = sorted({e["database"] for e in exps})
    out = match_rows(read_csv(cohort_file), by_sys, sources, tols if tols is not None else cfg.get("identity_tolerances"))
    path = out_path or cohort_file.replace(os.sep + "cohorts" + os.sep, os.sep + "dbmatch" + os.sep)
    write_csv(path, out)
    log(f"dbmatch: {os.path.basename(cohort_file)}: {sum(1 for o in out if o['db_first_T3'])} entries with a database report")
    return path
