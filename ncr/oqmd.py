"""Converter from OQMD MySQL exports (two tab-separated files) to a release snapshot.

Rules, set from the contents of qmdb v1.0 (April 2014) and v1.1 (October 2016):
  energy        formation_energies.delta_e for the fit 'standard' (the only fit in both)
  displayed     formation_energies.stability <= 0 means displayed as stable; NULL = unknown
  several rows  per entry (v1.0: 48,484 entries, 32,861 with values differing by more than
                1 meV/atom): the row with the lowest formation energy is kept, because every
                row is a point on the convex hull, so an entry is stable exactly when its
                lowest row is. The reconstruction check compares the recomputed hull with
                OQMD's displayed stability
  duplicates    entries with duplicate_of_id set to another entry are dropped
  origin        path under /libraries/icsd/ = experimental origin; '/test/' paths = test
                entries (kept in the hull, never in cohorts); other compounds = hypothetical
  prototype     taken from the path for heusler and half_heusler libraries (v1.0 leaves
                prototype_id NULL there), otherwise prototype_id
  known         source_known_composition = 1 when any ICSD-origin entry of the same reduced
                composition exists in the release, including ICSD entries without a
                formation energy and duplicates (the source listed the compound as known)
  references    an elemental reference at formation energy 0 is added for every element,
                matching the definition of the formation energy
"""
from __future__ import annotations

import csv
import sys

from .util import log, write_csv, write_json

STRUCTURE_LABEL_PRIORITY = ["standard", "fine_relax", "relaxed", "final", "static", "relaxation",
                            "coarse_relax", "initialize", "input"]
csv.field_size_limit(sys.maxsize)


def _null(v):
    return v is None or v == "" or v == "NULL"


def _tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        yield from csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)


def convert_oqmd(entries_tsv: str, structures_tsv: str | None, out_csv: str, fit: str = "standard") -> dict:
    from pymatgen.core import Composition
    rep = {"fit": fit, "rows_read": 0, "entries": 0, "entries_with_several_fe_rows": 0,
           "entries_with_conflicting_delta_e_gt_1meV": 0, "dropped_duplicates": 0, "dropped_no_formula": 0,
           "dropped_no_formation_energy": 0, "dropped_unparseable_formula": 0}
    best, latest, seen_de, icsd_comps = {}, {}, {}, set()
    for r in _tsv(entries_tsv):
        rep["rows_read"] += 1
        eid = r["entry_id"]
        if "/libraries/icsd/" in r.get("path", "") and not _null(r.get("composition_id")):
            try:
                icsd_comps.add(Composition(r["composition_id"].replace(" ", "")).reduced_formula)
            except Exception:  # noqa: BLE001
                pass
        if r.get("fit_id") != fit or _null(r.get("fe_delta_e")):
            best.setdefault(eid, None)
            continue
        de = float(r["fe_delta_e"])
        seen_de.setdefault(eid, []).append(de)
        cur = best.get(eid)
        if cur is None or de < float(cur["fe_delta_e"]) or (de == float(cur["fe_delta_e"]) and int(r["fe_id"]) > int(cur["fe_id"])):
            best[eid] = r
        lt = latest.get(eid)
        if lt is None or int(r["fe_id"]) > int(lt["fe_id"]):
            latest[eid] = r
    rep["entries_where_lowest_row_is_not_latest_row"] = sum(
        1 for e, r in best.items() if r is not None and latest.get(e) is not None and r["fe_id"] != latest[e]["fe_id"])
    for eid, des in seen_de.items():
        if len(des) > 1:
            rep["entries_with_several_fe_rows"] += 1
            if max(des) - min(des) > 0.001:
                rep["entries_with_conflicting_delta_e_gt_1meV"] += 1
    sg = {}
    if structures_tsv:
        rank = {lab: i for i, lab in enumerate(STRUCTURE_LABEL_PRIORITY)}
        for r in _tsv(structures_tsv):
            lab, e, s = r.get("label", ""), r.get("entry_id"), r.get("spacegroup_id")
            if _null(e) or _null(s) or lab not in rank:
                continue
            if e not in sg or rank[lab] < sg[e][0]:
                sg[e] = (rank[lab], s)
    rows, elements = [], set()
    counts = {"hypothetical": 0, "icsd_origin": 0, "elements": 0, "test": 0}
    for eid, r in best.items():
        rep["entries"] += 1
        if r is None:
            rep["dropped_no_formation_energy"] += 1
            continue
        dup = r.get("duplicate_of_id")
        if not _null(dup) and dup != eid:
            rep["dropped_duplicates"] += 1
            continue
        comp = r.get("composition_id")
        if _null(comp):
            rep["dropped_no_formula"] += 1
            continue
        try:
            c = Composition(comp.replace(" ", ""))
        except Exception:  # noqa: BLE001
            rep["dropped_unparseable_formula"] += 1
            continue
        path = r.get("path", "")
        n_el = len(c.elements)
        elements.update(str(e) for e in c.elements)
        if "/libraries/icsd/" in path:
            kind, origin = "icsd_origin", (r.get("label") if not _null(r.get("label")) else "icsd")
        elif "/test/" in path:
            kind, origin = "test", ""
        elif n_el == 1:
            kind, origin = "elements", ""
        else:
            kind, origin = "hypothetical", ""
        counts[kind] += 1
        if "/ternaries/half_heusler/" in path:
            proto = "C1_b (half_heusler)"
        elif "/ternaries/heusler/" in path:
            proto = "L2_1 (heusler)"
        else:
            proto = "" if _null(r.get("prototype_id")) else r["prototype_id"]
        stab = None if _null(r.get("fe_stability")) else float(r["fe_stability"])
        rows.append({"entry_id": f"oqmd-{eid}", "formula": c.formula.replace(" ", ""),
                     "energy_per_atom": float(r["fe_delta_e"]), "is_hypothetical": int(kind == "hypothetical"),
                     "displayed_stable": "" if stab is None else int(stab <= 0.0),
                     "displayed_dhd": "" if stab is None else stab,
                     "spacegroup_number": sg.get(eid, (0, ""))[1], "prototype": proto,
                     "origin_icsd_id": origin, "library": "/".join(path.split("/")[3:6]),
                     "_reduced": c.reduced_formula})
    for r in rows:
        r["source_known_composition"] = int(r["is_hypothetical"] == 1 and r["_reduced"] in icsd_comps)
    for el in sorted(elements):
        rows.append({"entry_id": f"ref-{el}", "formula": el, "energy_per_atom": 0.0, "is_hypothetical": 0,
                     "displayed_stable": "", "displayed_dhd": "", "spacegroup_number": "", "prototype": "reference",
                     "origin_icsd_id": "", "library": "reference", "source_known_composition": 0, "_reduced": el})
    stable_new = [r for r in rows if r["is_hypothetical"] and r["displayed_stable"] == 1 and not r["source_known_composition"]]
    rep.update(counts)
    rep["icsd_compositions_known_to_source"] = len(icsd_comps)
    rep["elemental_references_added"] = len(elements)
    rep["rows_written"] = len(rows)
    rep["displayed_stable_hypothetical_entries"] = sum(1 for r in rows if r["is_hypothetical"] and r["displayed_stable"] == 1)
    rep["displayed_stable_new_compositions"] = len({r["_reduced"] for r in stable_new})
    rep["displayed_stable_new_heusler_type_compositions"] = len({r["_reduced"] for r in stable_new if "heusler" in r["prototype"]})
    rep["hypothetical_with_unknown_displayed_label"] = sum(1 for r in rows if r["is_hypothetical"] and r["displayed_stable"] == "")
    fields = ["entry_id", "formula", "energy_per_atom", "is_hypothetical", "displayed_stable", "displayed_dhd",
              "spacegroup_number", "prototype", "origin_icsd_id", "source_known_composition", "library"]
    write_csv(out_csv, rows, fields)
    write_json(out_csv.replace(".csv", "_conversion_report.json"), rep)
    log(f"convert-oqmd: {rep}")
    return rep
