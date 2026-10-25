"""Release snapshots: reading, converting and recomputing the decomposition energy.

A snapshot is one CSV per database release (format in DATA_FORMATS.md). Required
columns: entry_id, formula, energy_per_atom, is_hypothetical. Optional columns:
displayed_stable, displayed_dhd, structure_file, spacegroup_number, prototype,
origin_icsd_id, deprecated. A config "columns" mapping renames your own column
names to these.
"""
from __future__ import annotations

import bz2
import gzip
import json
import math
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

from .util import log, read_csv, truthy, write_csv

REQUIRED = ["entry_id", "formula", "energy_per_atom", "is_hypothetical"]
OPTIONAL = ["displayed_stable", "displayed_dhd", "structure_file", "spacegroup_number",
            "prototype", "origin_icsd_id", "deprecated", "source_known_composition"]


def read_snapshot(path: str, columns: dict | None = None) -> list[dict]:
    """Read one release snapshot and normalize column names and types."""
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"{path} has no rows")
    columns = columns or {}
    out = []
    for i, r in enumerate(rows):
        g = {k: r.get(columns.get(k, k), "") for k in REQUIRED + OPTIONAL}
        for k in REQUIRED:
            if g[k] in ("", None):
                raise ValueError(f"{path} row {i + 2}: required column '{k}' is empty "
                                 f"(map your column names with 'columns' in the config)")
        g["energy_per_atom"] = float(g["energy_per_atom"])
        g["is_hypothetical"] = truthy(g["is_hypothetical"])
        g["displayed_stable"] = None if g["displayed_stable"] in ("", None) else truthy(g["displayed_stable"])
        g["displayed_dhd"] = None if g["displayed_dhd"] in ("", None) else float(g["displayed_dhd"])
        g["deprecated"] = truthy(g["deprecated"]) if g["deprecated"] not in ("", None) else False
        g["origin_icsd_id"] = (g["origin_icsd_id"] or "").strip()
        g["source_known_composition"] = truthy(g["source_known_composition"]) if g["source_known_composition"] not in ("", None) else False
        sf = (g["structure_file"] or "").strip()
        g["structure_file"] = (sf if (not sf or os.path.isabs(sf))
                               else os.path.join(os.path.dirname(os.path.abspath(path)), sf))
        out.append(g)
    return out


# ----------------------------------------------------------------------------
# Converter from pymatgen ComputedEntry / ComputedStructureEntry JSON dumps
# ----------------------------------------------------------------------------
def _open_any(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    if path.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def convert_entries(in_path: str, out_csv: str, hypothetical_ids: str | None = None,
                    all_hypothetical: bool = False, cif_dir: str | None = None,
                    id_field: str = "entry_id") -> int:
    """Write a snapshot CSV from a JSON file holding a list of pymatgen entry
    dictionaries (or {"entries": [...]}). Energies are the entries' corrected
    energies per atom. Which entries are hypothetical is NOT in a generic entry
    dump, so it must be given: a file of hypothetical IDs (one per line) or
    --all-hypothetical when every entry in the file is a hypothetical prediction."""
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
    if not hypothetical_ids and not all_hypothetical:
        raise ValueError("say which entries are hypothetical: --hypothetical-ids FILE or --all-hypothetical")
    hyp = None
    if hypothetical_ids:
        with open(hypothetical_ids, encoding="utf-8") as fh:
            hyp = {ln.strip() for ln in fh if ln.strip()}
    with _open_any(in_path) as fh:
        data = json.load(fh)
    items = data["entries"] if isinstance(data, dict) and "entries" in data else data
    if cif_dir:
        os.makedirs(cif_dir, exist_ok=True)
    rows = []
    for d in items:
        cls = ComputedStructureEntry if "structure" in d else ComputedEntry
        e = cls.from_dict(d)
        eid = str(d.get(id_field) or e.entry_id or d.get("data", {}).get(id_field) or "")
        if not eid:
            raise ValueError(f"an entry has no '{id_field}'; pass --id-field")
        row = {"entry_id": eid, "formula": e.composition.formula.replace(" ", ""),
               "energy_per_atom": e.energy_per_atom,
               "is_hypothetical": 1 if (all_hypothetical or eid in hyp) else 0}
        if cif_dir and cls is ComputedStructureEntry:
            p = os.path.join(cif_dir, f"{eid.replace('/', '_')}.cif")
            e.structure.to(filename=p)
            row["structure_file"] = os.path.abspath(p)
        rows.append(row)
    write_csv(out_csv, rows, REQUIRED + ["structure_file"])
    return len(rows)


# ----------------------------------------------------------------------------
# Decomposition energy over a whole snapshot
# ----------------------------------------------------------------------------
def _chemsys(formula: str) -> str:
    from pymatgen.core import Composition
    return "-".join(sorted(str(e) for e in Composition(formula).elements))


def _hull_for_system(args):
    """Worker: decomposition energies for the target entries of one chemical system."""
    system, pool, targets = args
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Composition
    ents = []
    for eid, f, epa in pool:
        c = Composition(f)
        ents.append(PDEntry(c, epa * c.num_atoms, name=eid))
    by_id = {e.name: e for e in ents}
    out = {}
    try:
        pd = PhaseDiagram(ents)
    except Exception as exc:  # noqa: BLE001  (missing elemental references, etc.)
        return {t: (math.nan, f"phase diagram failed: {exc}") for t in targets}
    for t in targets:
        try:
            # Energy above the hull for anything not on it (including higher-energy polymorphs of a
            # stable composition); the (negative) phase-separation energy only for hull members.
            # get_phase_separation_energy alone scores a higher polymorph against the hull without
            # its own composition, which wrongly makes it look stable (fixed in pipeline 1.5).
            eh = pd.get_e_above_hull(by_id[t])
            if eh is not None and eh > 1e-8:
                out[t] = (float(eh), "")
            else:
                out[t] = (float(pd.get_phase_separation_energy(by_id[t])), "")
        except Exception as exc:  # noqa: BLE001
            out[t] = (math.nan, f"{exc}")
    return out


def compute_dhd(rows: list[dict], target_filter=None, workers: int = 1) -> dict:
    """Decomposition energy (eV/atom) of every target entry against all other
    entries of the same snapshot: negative on the hull, energy above the hull
    otherwise (Methods). Entries with experimental origin stay in the hull,
    because the hull of a release includes everything the release computed.
    Returns {entry_id: (dhd, error_message)}."""
    by_sys = defaultdict(list)
    for r in rows:
        by_sys[_chemsys(r["formula"])].append((r["entry_id"], r["formula"], r["energy_per_atom"]))
    targets_by_sys = defaultdict(list)
    for r in rows:
        if target_filter is None or target_filter(r):
            targets_by_sys[_chemsys(r["formula"])].append(r["entry_id"])
    jobs = []
    for system, targets in targets_by_sys.items():
        els = system.split("-")
        pool = []
        for k in range(1, len(els) + 1):
            for sub in combinations(els, k):
                pool.extend(by_sys.get("-".join(sub), []))
        jobs.append((system, pool, targets))
    log(f"decomposition energies: {sum(len(t) for _, _, t in jobs)} targets in {len(jobs)} chemical systems")
    out = {}
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for res in ex.map(_hull_for_system, jobs, chunksize=8):
                out.update(res)
    else:
        for i, job in enumerate(jobs):
            out.update(_hull_for_system(job))
            if (i + 1) % 500 == 0:
                log(f"  {i + 1}/{len(jobs)} systems done")
    return out
