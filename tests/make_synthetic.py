"""Synthetic world for the end-to-end test: snapshots, experimental-database export,
literature corpus and control panel, with a known truth table. Nothing here is data."""
from __future__ import annotations

import datetime as dt
import json
import os
import sys

import numpy as np
from pymatgen.core import Composition

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from ncr.snapshots import compute_dhd          # noqa: E402
from ncr.util import write_csv, write_json     # noqa: E402

CATIONS = ["Li", "Na", "K", "Mg", "Ca", "Sr", "Ba", "Al", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga"]
ANIONS = ["O", "S", "Se", "F", "Cl", "N", "P", "Sb"]
STOICH = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (3, 4), (2, 5)]


def _date(y: float) -> str:
    year = int(y)
    d = dt.date(year, 1, 1) + dt.timedelta(days=int((y - year) * 365))
    return d.isoformat()


def make(outdir: str, seed: int = 7, effect: float = 0.15) -> dict:
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    snaps = {"srca_r1": [], "srca_r2": [], "srcb": []}
    systems = [(c, a) for c in CATIONS for a in ANIONS] + [(c1, c2) for i, c1 in enumerate(CATIONS[:10])
                                                          for c2 in CATIONS[i + 1:10]]
    rng.shuffle(systems)
    truth = []
    elements = sorted({e for s in systems for e in s})
    for snap in snaps:
        for e in elements:
            snaps[snap].append({"entry_id": f"el-{e}", "formula": e, "energy_per_atom": 0.0, "is_hypothetical": 0})
    for k, (A, X) in enumerate(systems):
        # one known ground state per system sets the hull scale
        gs = STOICH[k % 3]
        e_gs = -rng.uniform(0.3, 1.2)
        known = {"entry_id": f"k{k}", "formula": f"{A}{gs[0]}{X}{gs[1]}", "energy_per_atom": round(e_gs, 5),
                 "is_hypothetical": 0, "origin_icsd_id": f"icsd-{k}"}
        for snap in snaps:
            snaps[snap].append(dict(known))
        for j, (p, q) in enumerate(STOICH[3:9]):
            f = f"{A}{p}{X}{q}"
            x = q / (p + q)
            xg = gs[1] / sum(gs)
            hull = e_gs * (x / xg if x <= xg else (1 - x) / (1 - xg))
            base_e = hull + rng.normal(0.0, 0.03)
            eid = f"h{k}_{j}"
            sg = int(rng.integers(1, 231))
            r1 = {"entry_id": eid, "formula": f, "energy_per_atom": round(base_e, 5), "is_hypothetical": 1,
                  "spacegroup_number": sg}
            snaps["srca_r1"].append(dict(r1))
            snaps["srcb"].append({"entry_id": f"b{k}_{j}", "formula": f, "energy_per_atom": round(base_e + rng.normal(0, 0.02), 5),
                                  "is_hypothetical": 1})
            if rng.random() < 0.85:
                snaps["srca_r2"].append(dict(r1, energy_per_atom=round(base_e + rng.normal(0, 0.005), 5)))
    meta = {"srca_r1": ("SrcA", "2016-06-01", 0.0, "structure"), "srca_r2": ("SrcA", "2018-06-01", 0.0, "structure"),
            "srcb": ("SrcB", "2015-12-11", 0.025, "composition")}
    # displayed labels: recomputed with a few deliberate disagreements
    for snap, rows in snaps.items():
        tau = meta[snap][2]
        dhd = compute_dhd(rows, lambda r: r["is_hypothetical"] == 1)
        for r in rows:
            if r["entry_id"] in dhd:
                d = dhd[r["entry_id"]][0]
                st = d <= tau
                if rng.random() < 0.02:
                    st = not st
                r["displayed_stable"] = int(st)
                r["_dhd"] = d
        write_csv(os.path.join(outdir, "snapshots", f"{snap}.csv"),
                  [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows],
                  ["entry_id", "formula", "energy_per_atom", "is_hypothetical", "displayed_stable",
                   "spacegroup_number", "origin_icsd_id"])
    # Truth per compound (formula): first report date and first mention date
    docs, exps = [], []
    first_seen = {}
    for snap in ("srcb", "srca_r1", "srca_r2"):
        src, rd, tau, gran = meta[snap]
        e0 = int(rd[:4]) + 0.5
        for r in snaps[snap]:
            if r.get("is_hypothetical") != 1 or "_dhd" not in r:
                continue
            f = r["formula"]
            if f in first_seen:
                continue
            x = (r["_dhd"] - tau) * 1000
            stable = r["displayed_stable"] == 1
            first_seen[f] = (e0, x, stable, r.get("spacegroup_number"))
    sg_by_formula = {r["formula"]: r.get("spacegroup_number") for r in snaps["srca_r1"] if r.get("is_hypothetical") == 1}
    for f, (e0, x, stable, sg) in first_seen.items():
        sg = sg or sg_by_formula.get(f)
        u = rng.random()
        rep = None
        if u < 0.10:
            rep = e0 - rng.uniform(0.5, 20)
        elif u < 0.10 + 0.06 + (effect if stable and abs(x) <= 50 else 0.0) * 0.9:
            rep = e0 + rng.uniform(0.05, 9.5)
        men = rep if rep is not None else (e0 + rng.uniform(0.1, 6) if rng.random() < 0.3 else None)
        truth.append({"formula": f, "first_report": rep, "first_mention": men})
        if rep is not None:
            sgtxt = f" The structure was refined in space group {sg}." if sg else ""
            rf = Composition(f).reduced_formula
            docs.append({"doc_id": f"d{len(docs)}", "doi": f"10.0000/syn.{len(docs)}", "title": f"Synthesis of {rf}",
                         "text": f"We synthesized {rf} by solid-state reaction at 1100 K and characterized it by "
                                 f"powder X-ray diffraction.{sgtxt}", "publication_date": _date(rep)})
            if rng.random() < 0.6:
                exps.append({"exp_id": f"x{len(exps)}", "database": "ICSD", "formula": f,
                             "first_report_date": _date(rep), "experimental": 1})
        if men is not None and (rep is None or men < rep):
            rf = Composition(f).reduced_formula
            docs.append({"doc_id": f"d{len(docs)}", "doi": f"10.0000/calc.{len(docs)}", "title": f"Stability of {rf}",
                         "text": f"First-principles calculations suggest that {rf} could be stable at ambient pressure.",
                         "publication_date": _date(men)})
    # theoretical record that must be ignored
    exps.append({"exp_id": "xtheo", "database": "ICSD", "formula": list(first_seen)[0],
                 "first_report_date": "2000-01-01", "experimental": 0})
    with open(os.path.join(outdir, "corpus.jsonl"), "w", encoding="utf-8") as fh:
        for d in docs:
            fh.write(json.dumps(d) + "\n")
    write_csv(os.path.join(outdir, "experiments.csv"), exps)
    tf = [t for t in truth if t["first_report"] is not None]
    panel = [{"formula": t["formula"], "panel_kind": "positive", "entry_date": _date(t["first_report"] - 1),
              "known_first_report_date": _date(t["first_report"])} for t in tf[:15]]
    panel += [{"formula": t["formula"], "panel_kind": "decoy", "entry_date": "2016-01-01"}
              for t in truth if t["first_report"] is None and t["first_mention"] is None][:15]
    write_csv(os.path.join(outdir, "control_panel.csv"), panel)
    write_json(os.path.join(outdir, "truth.json"), {"truth": truth, "meta": meta, "effect": effect})
    cfg = {
        "workdir": "work", "cutoff_date": "2026-12-31", "allow_future_cutoff_for_testing": True,
        "releases": [{"source": meta[s][0], "release_id": s, "release_date": meta[s][1],
                      "path": f"snapshots/{s}.csv", "granularity": meta[s][3], "tau_ev_per_atom": meta[s][2],
                      "side_from": "displayed"} for s in ("srcb", "srca_r1", "srca_r2")],
        "experiments": {"path": "experiments.csv"},
        "retrieval": {"backend": "local", "local_corpus": "corpus.jsonl", "k": 20},
        "llm": {"model": "mock"},
        "horizons_years": [3, 5, 10], "ml_sources": [], "two_phase": {},
        "sampling": {"n_per_cohort": 160, "min_per_stratum": 30,
                     "assumed_error_rates": {"report": 0.2, "none": 0.05, "abstain": 0.3}},
        "threshold": {"sources": ["SrcB", "SrcA"], "figure_source": "SrcB", "windows_mev": [10, 25, 50],
                      "horizon_years": 5, "assumed_outcome_rate": 0.1, "effect_of_interest": 0.05},
        "premise": {"pairs": [{"a": "srcb", "b": "srca_r1", "match": "composition"}]},
        "control_panel": {"path": "control_panel.csv"},
        "seeds": {"expert_sample": 1, "two_phase": 2, "residual_audit": 3, "rerun": 4},
    }
    write_json(os.path.join(outdir, "config.json"), cfg)
    return cfg
