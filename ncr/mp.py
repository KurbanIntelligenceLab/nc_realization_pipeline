"""Converter from a Materials Project build (AWS Open Data, collections/<date>/thermo and
provenance) to a release snapshot.

Rules, fixed from the field listing of the 2022-10-28 build:
  thermo type   chosen by the record's own `thermo_type` field (folder names are not trusted;
                the report counts records by folder and by field)
  energy        formation_energy_per_atom of the chosen thermo type
  displayed     is_stable; the displayed decomposition energy is decomposition_enthalpy
                (negative on the hull), falling back to energy_above_hull
  hypothetical  provenance.theoretical is True; otherwise experimental origin
  known         source_known_composition = 1 when a non-theoretical material of the same
                reduced composition exists in the build
  references    an elemental reference at formation energy 0 is added for every element
The report also gives the year distribution of provenance.created_at, which decides
whether materials can be dated before the build.
"""
from __future__ import annotations

import collections
import glob
import gzip
import json
import os

from .util import log, write_csv, write_json


def _date(v):
    if isinstance(v, dict) and "$date" in v:
        v = v["$date"]
    return str(v)[:10] if v else ""


def _iter(files):
    for f in files:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for ln in fh:
                if ln.strip():
                    yield f, json.loads(ln)


def convert_mp(build_dir: str, out_csv: str, thermo_type: str = "GGA_GGA+U") -> dict:
    from pymatgen.core import Composition
    rep = {"build_dir": build_dir, "thermo_type": thermo_type}
    prov, created_years, created_years_theo = {}, collections.Counter(), collections.Counter()
    for _, d in _iter(glob.glob(os.path.join(build_dir, "provenance", "**", "*.jsonl.gz"), recursive=True)):
        mid = d.get("material_id")
        if not mid:
            continue
        prov[mid] = {"theoretical": bool(d.get("theoretical")), "created_at": _date(d.get("created_at")),
                     "icsd": bool((d.get("database_IDs") or {}).get("icsd"))}
        y = prov[mid]["created_at"][:4]
        created_years[y] += 1
        if prov[mid]["theoretical"]:
            created_years_theo[y] += 1
    rep["provenance_records"] = len(prov)
    rep["provenance_created_at_years_all"] = dict(sorted(created_years.items()))
    rep["provenance_created_at_years_theoretical"] = dict(sorted(created_years_theo.items()))
    by_folder = collections.Counter()
    chosen, dup = {}, 0
    for f, d in _iter(glob.glob(os.path.join(build_dir, "thermo", "**", "*.jsonl.gz"), recursive=True)):
        folder = next((p for p in f.split(os.sep) if p.startswith("thermo_type=")), "?")
        by_folder[f"{folder} | field={d.get('thermo_type')}"] += 1
        if d.get("thermo_type") != thermo_type:
            continue
        mid = d["material_id"]
        if mid in chosen:
            dup += 1
            continue
        chosen[mid] = {"formula": "".join(f"{el}{amt:g}" for el, amt in (d.get("composition_reduced") or {}).items()),
                       "fe": d.get("formation_energy_per_atom"), "stable": d.get("is_stable"),
                       "dhd": d.get("decomposition_enthalpy", d.get("energy_above_hull")),
                       "ehull": d.get("energy_above_hull"), "deprecated": bool(d.get("deprecated")),
                       "sg": (d.get("symmetry") or {}).get("number", "")}
    rep["thermo_records_by_folder_and_field"] = dict(sorted(by_folder.items()))
    rep["chosen_materials"] = len(chosen)
    rep["duplicate_records_skipped"] = dup
    if not chosen:
        raise ValueError(f"no thermo records with thermo_type == {thermo_type!r}; see the counts in the report")
    known, rows, elements = set(), [], set()
    no_prov = 0
    for mid, c in chosen.items():
        p = prov.get(mid)
        if p is None:
            no_prov += 1
        if p is not None and not p["theoretical"]:
            known.add(Composition(c["formula"]).reduced_formula)
    for mid, c in chosen.items():
        if c["fe"] is None:
            continue
        comp = Composition(c["formula"])
        elements.update(str(e) for e in comp.elements)
        p = prov.get(mid)
        hyp = bool(p and p["theoretical"]) and len(comp.elements) > 1
        rows.append({"entry_id": mid, "formula": c["formula"], "energy_per_atom": float(c["fe"]),
                     "is_hypothetical": int(hyp),
                     "displayed_stable": "" if c["stable"] is None else int(bool(c["stable"])),
                     "displayed_dhd": "" if c["dhd"] is None else float(c["dhd"]),
                     "spacegroup_number": c["sg"], "prototype": "",
                     "origin_icsd_id": "" if (p is None or p["theoretical"]) else ("icsd" if p["icsd"] else "experimental"),
                     "deprecated": int(c["deprecated"]),
                     "source_known_composition": int(hyp and comp.reduced_formula in known),
                     "created_at": p["created_at"] if p else ""})
    for el in sorted(elements):
        rows.append({"entry_id": f"ref-{el}", "formula": el, "energy_per_atom": 0.0, "is_hypothetical": 0,
                     "displayed_stable": "", "displayed_dhd": "", "spacegroup_number": "", "prototype": "reference",
                     "origin_icsd_id": "", "deprecated": 0, "source_known_composition": 0, "created_at": ""})
    hyp_rows = [r for r in rows if r["is_hypothetical"]]
    rep.update({"materials_without_provenance": no_prov, "rows_written": len(rows),
                "hypothetical": len(hyp_rows), "experimental_origin": sum(1 for r in rows if r["origin_icsd_id"]),
                "deprecated": sum(r["deprecated"] for r in rows),
                "displayed_stable_hypothetical_entries": sum(1 for r in hyp_rows if r["displayed_stable"] == 1),
                "displayed_stable_new_compositions": len({Composition(r["formula"]).reduced_formula for r in hyp_rows
                                                          if r["displayed_stable"] == 1 and not r["source_known_composition"]}),
                "elemental_references_added": len(elements)})
    fields = ["entry_id", "formula", "energy_per_atom", "is_hypothetical", "displayed_stable", "displayed_dhd",
              "spacegroup_number", "prototype", "origin_icsd_id", "deprecated", "source_known_composition", "created_at"]
    write_csv(out_csv, rows, fields)
    write_json(out_csv.replace(".csv", "_conversion_report.json"), rep)
    log(f"convert-mp: {rep}")
    return rep
