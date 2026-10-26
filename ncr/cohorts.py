"""Convex-hull recomputation for every release, prediction cohorts and the threshold
population. A prediction is a composition that a release displayed as stable and presented
as hypothetical, excluding compositions the release lists as experimental; it enters its
cohort at the first release that displayed it as stable."""
from __future__ import annotations

import math
import os

from . import config as C
from . import identity
from .snapshots import compute_dhd, read_snapshot
from .util import date_str, decimal_year, log, read_csv, write_csv, write_json

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
            "chemistry_class": identity.chemistry_class(formula),
            "mean_electronegativity": round(sum(xs) / len(xs), 4) if xs else "",
            "max_atomic_number": max(e.Z for e in c.elements)}


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


def step_cohorts(cfg: dict) -> None:
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
                    continue                                   # experimental origin (Supplementary Table 1)
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
                from .checks import block_keys, en_quintile_edges
                edges = en_quintile_edges(thr)
                for r, k in zip(thr, block_keys(thr, scheme, edges)):
                    r["block"] = k
                write_json(os.path.join(cfg["workdir"], "checks", f"threshold_blocks__{_slug(src)}.json"),
                    {"scheme": scheme, "electronegativity_quintile_edges": edges})
            write_csv(os.path.join(cdir, f"threshold__{_slug(src)}.csv"), thr)
            log(f"cohorts: {src}: threshold population of {len(thr)} entries within ±{WINDOW_MAX_MEV:.0f} meV/atom "
                f"({n_unknown} excluded because no label was displayed; running variable '{rv_rule}', "
                f"{n_rv_fallback} entries fell back to the recomputed value)")


def _slug(s: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_").lower()
