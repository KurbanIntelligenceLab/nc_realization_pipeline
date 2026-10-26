"""A synthetic world with known truth for the end-to-end tests: three release snapshots,
two experimental-record sources with partial overlap, and a stable-label effect near the
threshold. Nothing here is data."""
from __future__ import annotations

import datetime as dt
import json
import os

import numpy as np

from ncr.snapshots import compute_dhd
from ncr.util import write_csv

CATIONS = ["Li", "Na", "K", "Mg", "Ca", "Sr", "Ba", "Al", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga"]
ANIONS = ["O", "S", "Se", "F", "Cl", "N", "P", "Sb"]
STOICH = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (3, 4), (2, 5)]


def _date(y: float) -> str:
    year = int(y)
    return (dt.date(year, 1, 1) + dt.timedelta(days=int((y - year) * 365))).isoformat()


def make(outdir: str, seed: int = 7, effect: float = 0.15) -> dict:
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    snaps = {"srca_r1": [], "srca_r2": [], "srcb": []}
    systems = [(c, a) for c in CATIONS for a in ANIONS] + [(c1, c2) for i, c1 in enumerate(CATIONS[:10])
                                                          for c2 in CATIONS[i + 1:10]]
    rng.shuffle(systems)
    elements = sorted({e for s in systems for e in s})
    for snap in snaps:
        for e in elements:
            snaps[snap].append({"entry_id": f"ref-{e}", "formula": e, "energy_per_atom": 0.0, "is_hypothetical": 0})
    for k, (A, X) in enumerate(systems):
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
            r1 = {"entry_id": f"h{k}_{j}", "formula": f, "energy_per_atom": round(base_e, 5), "is_hypothetical": 1,
                  "prototype": "L2_1 (heusler)" if j % 2 else "B1"}
            snaps["srca_r1"].append(dict(r1))
            snaps["srcb"].append({"entry_id": f"b{k}_{j}", "formula": f, "is_hypothetical": 1,
                                  "energy_per_atom": round(base_e + rng.normal(0, 0.02), 5), "prototype": r1["prototype"]})
            if rng.random() < 0.85:
                snaps["srca_r2"].append(dict(r1, energy_per_atom=round(base_e + rng.normal(0, 0.005), 5)))
    meta = {"srca_r1": ("SrcA", "2016-06-01", "composition"), "srca_r2": ("SrcA", "2018-06-01", "composition"),
            "srcb": ("SrcB", "2015-12-11", "composition")}
    for snap, rows in snaps.items():
        dhd = compute_dhd(rows, lambda r: r["is_hypothetical"] == 1)
        for r in rows:
            if r["entry_id"] in dhd:
                d = dhd[r["entry_id"]][0]
                r["displayed_stable"] = int(d <= 0.0)
                r["displayed_dhd"] = round(d, 6)
        write_csv(os.path.join(outdir, "snapshots", f"{snap}.csv"), rows,
                  ["entry_id", "formula", "energy_per_atom", "is_hypothetical", "displayed_stable", "displayed_dhd",
                   "prototype", "origin_icsd_id"])
    first_seen = {}
    for snap in ("srcb", "srca_r1", "srca_r2"):
        e0 = int(meta[snap][1][:4]) + 0.5
        for r in snaps[snap]:
            if r.get("is_hypothetical") == 1 and r["formula"] not in first_seen:
                first_seen[r["formula"]] = (e0, r["displayed_dhd"] * 1000, r["displayed_stable"] == 1)
    cod, icsd, truth = [], [], {}
    for f, (e0, x, stable) in first_seen.items():
        u = rng.random()
        rep = None
        if u < 0.10:
            rep = e0 - rng.uniform(0.5, 20)
        elif u < 0.16 + (effect if stable and abs(x) <= 50 else 0.0) * 0.9:
            rep = e0 + rng.uniform(0.05, 9.5)
        truth[f] = rep
        if rep is None:
            continue
        v = rng.random()
        if v < 0.7:
            cod.append({"exp_id": f"c{len(cod)}", "database": "COD", "formula": f,
                        "first_report_date": str(int(rep)), "experimental": 1, "high_pressure": int(rng.random() < 0.05)})
        if v > 0.4:
            icsd.append({"exp_id": f"m{len(icsd)}", "database": "MP-ICSD", "formula": f,
                         "first_report_date": _date(rep), "experimental": 1})
    cod.append({"exp_id": "ctheo", "database": "COD", "formula": list(first_seen)[0],
                "first_report_date": "2000", "experimental": 0})
    write_csv(os.path.join(outdir, "cod.csv"), cod)
    write_csv(os.path.join(outdir, "mp_icsd.csv"), icsd)
    cfg = {"workdir": "work", "cutoff_date": "2026-12-31", "allow_future_cutoff": True,
           "releases": [{"source": meta[s][0], "release_id": s, "release_date": meta[s][1],
                         "label": s, "path": f"snapshots/{s}.csv", "granularity": meta[s][2],
                         "tau_ev_per_atom": 0.0, "side_from": "displayed"} for s in ("srcb", "srca_r1", "srca_r2")],
           "experiments": {"path": "experiments.csv"},
           "threshold": {"sources": ["SrcB"], "windows_mev": [10, 25, 50], "horizon_years": 5,
                         "assumed_outcome_rate": 0.1, "effect_of_interest": 0.05, "block_scheme": "class_en_nel"},
           "premise": {"pairs": [{"a": "srcb", "b": "srca_r1", "match": "composition"}]},
           "seeds": {"one_entry_per_chemical_system": 1}}
    with open(os.path.join(outdir, "config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2)
    with open(os.path.join(outdir, "truth.json"), "w") as fh:
        json.dump({"first_report": truth, "effect": effect}, fh)
    return cfg
