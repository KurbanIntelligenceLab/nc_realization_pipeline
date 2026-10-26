"""Checks that use no prediction outcome: reconstruction of each release, power and balance
of the threshold design (check b), cross-database label disagreement (check a), and the
coverage of the outcome record on compounds the databases already list as experimental."""
from __future__ import annotations

import math
import os

import numpy as np

from . import config as C
from . import stats as S
from .cohorts import _comp, _hull_table, _slug, _stable
from .matching import read_experiments
from .snapshots import read_snapshot
from .util import log, read_csv


def step_reconstruction(cfg: dict) -> dict:
    """Snapshot reconstruction: stable hypothetical entries and new compositions under the
    displayed and the recomputed labels, compared with any expected counts in the config
    (relative tolerance reconstruction_tolerance, default 1%)."""
    out = {}
    tol = float(cfg.get("reconstruction_tolerance", 0.01))
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
        log(f"reconstruction: {r['release_id']}: displayed {rec['displayed']}, recomputed {rec['recomputed']}")
    return out


def balance(rows, inw, side, block) -> dict:
    """Check (b), run before any outcome: standardized differences of covariates fixed before
    release between the two sides, overall in the window and for the blocked estimator."""
    covs = {"n_elements": np.array([float(r["n_elements"]) for r in rows]),
            "mean_electronegativity": np.array([float(r["mean_electronegativity"] or "nan") for r in rows]),
            "max_atomic_number": np.array([float(r["max_atomic_number"]) for r in rows])}
    for cls in sorted({r["chemistry_class"] for r in rows}):
        covs[f"class_{cls}"] = np.array([1.0 if r["chemistry_class"] == cls else 0.0 for r in rows])

    def smds(mask):
        A, B = mask & (side == "A"), mask & (side == "B")
        if A.sum() < 2 or B.sum() < 2:
            return None
        return {k: S.standardized_difference(v[A & np.isfinite(v)], v[B & np.isfinite(v)]) for k, v in covs.items()}
    overall = smds(inw)
    # Balance that matters for the blocked estimator: block-size-weighted within-block mean
    # differences over usable blocks (two or more entries per side), standardized by the pooled
    # standard deviation in the window.
    blocked = {}
    usable = [bl for bl in np.unique(block[inw])
              if ((inw & (block == bl) & (side == "A")).sum() >= 2 and (inw & (block == bl) & (side == "B")).sum() >= 2)]
    nb = {bl: int((inw & (block == bl)).sum()) for bl in usable}
    tot = sum(nb.values())
    for k, v in covs.items():
        if k.startswith("class_") or not tot:
            continue
        fin = inw & np.isfinite(v)
        sd = math.sqrt((np.var(v[fin & (side == "A")], ddof=1) + np.var(v[fin & (side == "B")], ddof=1)) / 2) \
            if (fin & (side == "A")).sum() > 1 and (fin & (side == "B")).sum() > 1 else 0.0
        diff = 0.0
        for bl in usable:
            mA = v[fin & (block == bl) & (side == "A")]
            mB = v[fin & (block == bl) & (side == "B")]
            if mA.size and mB.size:
                diff += nb[bl] / tot * (mA.mean() - mB.mean())
        blocked[k] = 0.0 if sd == 0 else float(diff / sd)
    within = {}
    for bl in np.unique(block[inw]):
        r = smds(inw & (block == bl))
        if r is not None:
            within[str(bl)] = {k: v for k, v in r.items() if not k.startswith("class_")}
    mx = max((abs(v) for v in (overall or {}).values()), default=None)
    mxb = max((abs(v) for v in blocked.values()), default=None)
    return {"overall": overall, "max_abs_smd_overall": mx, "check_b_passes": (mx is not None and mx < 0.1),
            "blocked_smd": blocked, "max_abs_blocked_smd": mxb,
            "check_b_passes_blocked": (mxb is not None and mxb < 0.1), "within_blocks": within}


def block_keys(rows, scheme: str, en_edges=None):
    """Block labels for the threshold design. Schemes: 'class' (source x chemistry class, as
    first specified); 'class_en' adds quintiles of mean electronegativity; 'class_en_nel' (the
    primary scheme) also splits by number of elements. Quintile edges come from the whole
    threshold population and use no outcome."""
    out = []
    for r in rows:
        k = f"{r['source']}|{r['chemistry_class']}"          # base key, independent of any applied scheme
        if scheme in ("class_en", "class_en_nel"):
            x = float(r["mean_electronegativity"]) if r.get("mean_electronegativity") not in ("", None) else float("nan")
            q = int(np.searchsorted(en_edges, x, side="right")) if np.isfinite(x) else -1
            k += f"|en{q}"
        if scheme == "class_en_nel":
            k += f"|n{r['n_elements']}"
        out.append(k)
    return np.array(out, dtype=object)


def en_quintile_edges(rows):
    v = np.array([float(r["mean_electronegativity"]) for r in rows if r.get("mean_electronegativity") not in ("", None)])
    v = v[np.isfinite(v)]
    return np.quantile(v, [0.2, 0.4, 0.6, 0.8]).tolist() if v.size else []


def step_power(cfg: dict) -> dict:
    """Window counts, minimum detectable effect at an assumed outcome rate, and check (b)
    for each blocking scheme and window. Uses no outcome."""
    tc = cfg.get("threshold", {})
    p = float(C.need(cfg, "threshold.assumed_outcome_rate"))
    out = {}
    for src in tc.get("sources", []):
        path = os.path.join(cfg["workdir"], "cohorts", f"threshold__{_slug(src)}.csv")
        rows = read_csv(path)
        x = np.array([float(r["x_mev"]) for r in rows])
        side = np.array([r["side"] for r in rows])
        block = np.array([r["block"] for r in rows])
        rec = {}
        for w in tc.get("windows_mev", [10, 25, 50]):
            inw = np.abs(x) <= w
            nA, nB = int((inw & (side == "A")).sum()), int((inw & (side == "B")).sum())
            blocks = {}
            for bl in np.unique(block[inw]):
                m = inw & (block == bl)
                blocks[str(bl)] = {"A": int((m & (side == "A")).sum()), "B": int((m & (side == "B")).sum())}
            usable = sum(v["A"] + v["B"] for v in blocks.values() if v["A"] >= 2 and v["B"] >= 2)
            mde = S.minimum_detectable_effect(p, nA, nB) if nA >= 2 and nB >= 2 else None
            rec[str(w)] = {"N_A": nA, "N_B": nB, "entries_in_usable_blocks": usable, "blocks": blocks,
                           "mde_at_assumed_rate": mde,
                           "balance": balance(rows, inw, side, block)}
            if mde is not None and tc.get("effect_of_interest") is not None:
                rec[str(w)]["mde_exceeds_effect_of_interest"] = mde > float(tc["effect_of_interest"])
        edges = en_quintile_edges(rows)
        by_window = {}
        for wv in tc.get("windows_mev", [10, 25, 50]):
            cw = {}
            inw_w = np.abs(x) <= wv
            for scheme in ("class", "class_en", "class_en_nel"):
                bk = block_keys(rows, scheme, edges)
                bal = balance(rows, inw_w, side, bk)
                cw[scheme] = {"max_abs_blocked_smd": bal["max_abs_blocked_smd"], "passes": bal["check_b_passes_blocked"],
                              "blocked_smd": bal["blocked_smd"]}
            by_window[str(wv)] = cw
        cand = {}
        inw = np.abs(x) <= 25
        for scheme in ("class", "class_en", "class_en_nel"):
            bk = block_keys(rows, scheme, edges)
            bal = balance(rows, inw, side, bk)
            uA = uB = 0
            for bl in np.unique(bk[inw]):
                m = inw & (bk == bl)
                a, b_ = int((m & (side == "A")).sum()), int((m & (side == "B")).sum())
                if a >= 2 and b_ >= 2:
                    uA += a
                    uB += b_
            cand[scheme] = {"blocked_smd": bal["blocked_smd"], "max_abs_blocked_smd": bal["max_abs_blocked_smd"],
                            "passes": bal["check_b_passes_blocked"], "N_A_usable": uA, "N_B_usable": uB,
                            "usable_blocks": sum(1 for bl in np.unique(bk[inw])
                                                 if ((inw & (bk == bl) & (side == "A")).sum() >= 2
                                                     and (inw & (bk == bl) & (side == "B")).sum() >= 2)),
                            "mde_at_assumed_rate": S.minimum_detectable_effect(p, uA, uB) if uA >= 2 and uB >= 2 else None}
        rec["block_scheme_candidates_25meV"] = {"electronegativity_quintile_edges": edges, "schemes": cand}
        rec["block_scheme_balance_by_window"] = by_window
        out[src] = rec
        log(f"power: {src}: 25 meV window N_A={rec.get('25', {}).get('N_A')}, N_B={rec.get('25', {}).get('N_B')}, "
            f"MDE at the assumed rate {rec.get('25', {}).get('mde_at_assumed_rate')}")
        for sc, c in cand.items():
            log(f"power: {src}: blocks '{sc}': max |blocked SMD| {c['max_abs_blocked_smd']}, passes {c['passes']}, "
                f"usable A/B {c['N_A_usable']}/{c['N_B_usable']}, MDE {c['mde_at_assumed_rate']}")
    return out


def step_premise(cfg: dict) -> dict:
    """Check (a): cross-database label disagreement by distance to the threshold.
    Works from the snapshots: each database's composition is represented by its
    lowest-energy entry; label and distance are the values each database displayed
    (falling back to the recomputed hull where a value was not displayed)."""
    from pymatgen.core import Composition
    out = {}
    rel = {r["release_id"]: r for r in C.releases(cfg)}

    def table(r, mode):
        snap = read_snapshot(C.resolve(cfg, r["path"]), r.get("columns"))
        hp = os.path.join(cfg["workdir"], "hull", f"{r['release_id']}.csv")
        hull = {x["entry_id"]: x for x in read_csv(hp)} if os.path.exists(hp) else {}
        tau = float(C.need(r, "tau_ev_per_atom"))
        best = {}
        for s in snap:
            if s["entry_id"].startswith("ref-") or s["deprecated"]:
                continue
            comp = Composition(s["formula"])
            if len(comp.elements) < 2:
                continue
            k = comp.reduced_formula
            if mode == "composition+spacegroup":
                if not s["spacegroup_number"]:
                    continue
                k = f"{k}|{s['spacegroup_number']}"
            if k in best and s["energy_per_atom"] >= best[k]["e"]:
                continue
            lab, dhd = s["displayed_stable"], s["displayed_dhd"]
            h = hull.get(s["entry_id"])
            if lab is None and h and h["stable_recomputed"] != "":
                lab = h["stable_recomputed"] == "1"
            if dhd is None and h and h["dhd_recomputed"] != "":
                dhd = float(h["dhd_recomputed"])
            if lab is None or dhd is None:
                continue
            best[k] = {"e": s["energy_per_atom"], "stable": bool(lab), "x": (dhd - tau) * 1000.0}
        return best

    for pair in C.need(cfg, "premise.pairs"):
        key = f"{pair['a']} vs {pair['b']}"
        ta, tb = table(rel[pair["a"]], pair.get("match")), table(rel[pair["b"]], pair.get("match"))
        common = sorted(set(ta) & set(tb))
        if not common:
            out[key] = {"note": "no matched compounds"}
            continue
        dist = np.array([ta[k]["x"] for k in common])
        sa = np.array([ta[k]["stable"] for k in common])
        sb = np.array([tb[k]["stable"] for k in common])
        edges = cfg.get("premise", {}).get("edges_mev", [0, 10, 25, 50, 100, 200, 1e9])
        rows = S.label_disagreement_by_bin(dist, sa, sb, edges)
        nmask, far = np.abs(dist) < 25, np.abs(dist) >= 100
        rate = lambda m: float((sa[m] != sb[m]).mean()) if m.any() else None
        lo25 = S.wilson_ci(int((sa[nmask] != sb[nmask]).sum()), int(nmask.sum()))[0] if nmask.any() else None
        out[key] = {"matched": len(common), "stable_in_a": int(sa.sum()), "stable_in_b": int(sb.sum()),
                    "bins": rows, "within_25_n": int(nmask.sum()), "within_25_rate": rate(nmask),
                    "within_25_wilson_low": lo25, "beyond_100_rate": rate(far),
                    "check_a_passes": bool(lo25 is not None and lo25 > 0.20 and rate(far) is not None
                                           and rate(far) < rate(nmask))}
        log(f"premise: {key}: matched {len(common)}; within 25 meV disagreement {rate(nmask)}; "
            f"beyond 100 {rate(far)}; check (a) passes = {out[key]['check_a_passes']}")
    return out


def step_coverage_known(cfg: dict) -> dict:
    """Share of each release's ICSD-origin compositions that also appear in the experimental
    records (for example COD). It reads only compounds the source already listed as
    experimental, which are excluded from every cohort, so no prediction outcome is touched.
    Matching is by exact reduced formula, a conservative approximation of the T3 rule."""
    from pymatgen.core import Composition
    import warnings
    warnings.filterwarnings("ignore")
    exps = read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
    # The paper reports this coverage for COD, the source that is independent of both databases:
    # the ICSD-derived records are themselves drawn from the releases' experimental entries.
    exps = [e for e in exps if e["database"] == "COD"]
    exp_red, exp_first = set(), {}
    for e in exps:
        try:
            red = Composition(e["formula"]).reduced_formula
        except Exception:  # noqa: BLE001
            continue
        exp_red.add(red)
        if e["date"] and (red not in exp_first or e["date"] < exp_first[red]):
            exp_first[red] = e["date"]
    out = {"experimental_reports_compositions": len(exp_red)}
    for r in C.releases(cfg):
        snap = read_snapshot(C.resolve(cfg, r["path"]), r.get("columns"))
        comps = {}
        for x in snap:
            if not x["origin_icsd_id"] or x["entry_id"].startswith("ref-"):
                continue
            c = Composition(x["formula"])
            if len(c.elements) < 2:
                continue
            comps[c.reduced_formula] = len(c.elements)
        if not comps:
            continue
        found = {k for k in comps if k in exp_red}
        from .identity import chemistry_class
        by_cls = {}
        for k in comps:
            cls = chemistry_class(k)
            d = by_cls.setdefault(cls, {"n": 0, "found": 0})
            d["n"] += 1
            d["found"] += int(k in exp_red)
        for d in by_cls.values():
            d["share"] = d["found"] / d["n"]
        by_n = {}
        for n in sorted(set(comps.values())):
            ks = [k for k, v in comps.items() if v == n]
            by_n[str(n)] = {"n": len(ks), "found": sum(1 for k in ks if k in exp_red),
                            "share": sum(1 for k in ks if k in exp_red) / len(ks)}
        rel_date = r["date"]
        before = sum(1 for k in found if exp_first[k] <= rel_date)
        out[r["release_id"]] = {"icsd_origin_compositions": len(comps), "found_in_reports": len(found),
                                "share_found": len(found) / len(comps), "by_number_of_elements": by_n,
                                "by_chemistry_class": by_cls,
                                "found_with_report_dated_on_or_before_release": before}
        log(f"coverage-known: {r['release_id']}: {len(found)}/{len(comps)} ICSD-origin compositions "
            f"({len(found) / len(comps):.1%}) appear in COD")
    return out
