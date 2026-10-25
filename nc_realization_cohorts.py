#!/usr/bin/env python3
"""
nc_realization_cohorts.py

Cohort construction and outcome-identity rules for the realization audit.
This file codifies the three outcome tiers (Table 2 of the manuscript), the
decomposition-energy running variable of the threshold design, and the cohort
entry / pre-existence / realization logic. It does not download data: every
loader raises if the expected snapshot file is missing (no silent fallbacks).

Tiers (a predicted entry P versus an experimental report E):
  T1  structure match: E is ordered and P matches E with pymatgen StructureMatcher
      (ElementComparator, primitive cells, volume scaling).
  T2  order-disorder match: E has partial occupancy and P is an ordered
      decoration of E (OrderDisorderElementComparator, supercells allowed),
      with equal reduced composition within comp_delta.
  T3  composition match: same element set and normalized composition within
      comp_delta (L1 distance), any structure.
The primary outcome is T2-or-better for structure-level predictions and T3 for
composition-level predictions (OQMD 2015 prototype list; see Methods).

Pinned and tested with pymatgen 2026.10.24 (Python 3.12). Run:
    python3 nc_realization_cohorts.py --selftest
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

DEFAULT_TOLS = {"ltol": 0.2, "stol": 0.3, "angle_tol": 5.0, "comp_delta": 0.02}


def _require_pymatgen():
    try:
        import pymatgen  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError("pymatgen is required: pip install pymatgen==2026.10.24") from exc


# ----------------------------------------------------------------------------
# Identity tiers
# ----------------------------------------------------------------------------
def composition_distance(comp_a, comp_b) -> float:
    """L1 distance between normalized fractional compositions; inf if the
    element sets differ."""
    _require_pymatgen()
    from pymatgen.core import Composition
    a = Composition(comp_a).fractional_composition
    b = Composition(comp_b).fractional_composition
    ea = {str(e) for e in a.elements}
    eb = {str(e) for e in b.elements}
    if ea != eb:
        return float("inf")
    return float(sum(abs(a[e] - b[e]) for e in ea))


def outcome_tier(pred, exp, tols: dict | None = None):
    """Return 'T1', 'T2', 'T3' or None for predicted structure `pred` versus
    experimental structure `exp` (pymatgen Structure objects). `exp` may be a
    composition string when only a composition-level report exists (then at most T3)."""
    _require_pymatgen()
    from pymatgen.analysis.structure_matcher import (ElementComparator,
                                                     OrderDisorderElementComparator,
                                                     StructureMatcher)
    from pymatgen.core import Structure
    t = dict(DEFAULT_TOLS, **(tols or {}))
    pred_comp = pred.composition if isinstance(pred, Structure) else pred
    exp_comp = exp.composition if isinstance(exp, Structure) else exp
    same_comp = composition_distance(pred_comp, exp_comp) <= t["comp_delta"]
    if isinstance(pred, Structure) and isinstance(exp, Structure) and same_comp:
        if exp.is_ordered:
            sm = StructureMatcher(ltol=t["ltol"], stol=t["stol"], angle_tol=t["angle_tol"],
                                  primitive_cell=True, scale=True, attempt_supercell=False,
                                  comparator=ElementComparator())
            if sm.fit(pred, exp):
                return "T1"
        else:
            sm = StructureMatcher(ltol=t["ltol"], stol=t["stol"], angle_tol=t["angle_tol"],
                                  primitive_cell=False, scale=True, attempt_supercell=True,
                                  comparator=OrderDisorderElementComparator())
            if sm.fit(pred, exp):
                return "T2"
    return "T3" if same_comp else None


def in_solid_solution_range(pred_comp, end_a, end_b, comp_delta: float = 0.02) -> bool:
    """T2 solid-solution rule: a report of a continuous solid solution between
    end-member compositions end_a and end_b (same structure type) matches a
    predicted composition that lies on the segment between them, within the
    L1 tolerance used elsewhere. Structure-type agreement is judged separately
    (T2 matcher or expert judgment)."""
    _require_pymatgen()
    import numpy as np
    from pymatgen.core import Composition
    comps = [Composition(x).fractional_composition for x in (pred_comp, end_a, end_b)]
    els = sorted({str(e) for c in comps for e in c.elements})
    p, a, b = (np.array([c[e] for e in els]) for c in comps)
    d = b - a
    lam = float(np.clip(np.dot(p - a, d) / np.dot(d, d), 0.0, 1.0)) if np.dot(d, d) > 0 else 0.0
    return float(np.abs(p - (a + lam * d)).sum()) <= comp_delta


TIER_RANK = {"T1": 3, "T2": 2, "T3": 1, None: 0}


def meets(tier, required: str) -> bool:
    """True if `tier` is at least as strict as `required` (T1 > T2 > T3)."""
    return TIER_RANK[tier] >= TIER_RANK[required]


# ----------------------------------------------------------------------------
# Running variable for the threshold design
# ----------------------------------------------------------------------------
def decomposition_energy(entries, target):
    """Decomposition energy (eV/atom) of `target` against the other entries of
    its phase diagram, computed from ONE release snapshot. Negative for entries on
    the hull (distance below the hull formed without them), positive and equal to
    the energy above hull otherwise."""
    _require_pymatgen()
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    pd = PhaseDiagram(entries)
    return float(pd.get_phase_separation_energy(target))


# ----------------------------------------------------------------------------
# Cohort logic from structured records
# ----------------------------------------------------------------------------
def build_cohort(pred_records, exp_records, required_tier: str, tier_fn=outcome_tier):
    """pred_records: dicts with keys pred_id, release_year, stable (bool),
        composition, structure (optional; None for composition-level predictions).
    exp_records: dicts with keys exp_id, first_report_year, composition,
        structure (optional).
    Returns one dict per prediction that is ever labeled stable:
        entry_year, preexisting (bool), realization_year (float or None)."""
    for r in pred_records:
        for k in ("pred_id", "release_year", "stable", "composition"):
            if k not in r:
                raise KeyError(f"prediction record missing {k}")
    for r in exp_records:
        for k in ("exp_id", "first_report_year", "composition"):
            if k not in r:
                raise KeyError(f"experimental record missing {k}")
    _require_pymatgen()
    from pymatgen.core import Composition
    by_sys = defaultdict(list)
    for r in exp_records:
        by_sys[Composition(r["composition"]).chemical_system].append(r)
    first_stable = {}
    for r in sorted(pred_records, key=lambda x: x["release_year"]):
        if r["stable"] and r["pred_id"] not in first_stable:
            first_stable[r["pred_id"]] = r
    out = []
    for pid, r in first_stable.items():
        sysk = Composition(r["composition"]).chemical_system
        pred_obj = r.get("structure") or r["composition"]
        matches = []
        for e in by_sys.get(sysk, []):
            exp_obj = e.get("structure") or e["composition"]
            tier = tier_fn(pred_obj, exp_obj)
            if meets(tier, required_tier):
                matches.append(float(e["first_report_year"]))
        entry = float(r["release_year"])
        pre = any(y <= entry for y in matches)
        later = [y for y in matches if y > entry]
        out.append({"pred_id": pid, "entry_year": entry, "preexisting": pre,
                    "realization_year": (min(later) if (later and not pre) else None)})
    return out


# ----------------------------------------------------------------------------
# Chemistry class (blocking and matching variable) and correction-driven flips
# ----------------------------------------------------------------------------
ANION_FAMILIES = {
    "oxide": {"O"},
    "halide": {"F", "Cl", "Br", "I"},
    "chalcogenide": {"S", "Se", "Te"},
    "pnictide": {"N", "P", "As", "Sb"},
    "hydride": {"H"},
    "carbide-boride": {"C", "B", "Si", "Ge"},
}


def chemistry_class(composition) -> str:
    """Deterministic class from the element set: the family of the most
    electronegative element (Pauling scale), so phosphates and sulfates count
    as oxides and thiophosphates as chalcogenides. 'intermetallic' if that
    element belongs to no listed family; 'element' for elemental entries. The
    rule is compositional, so it is computed for every entry before outcomes are
    known; it does not detect mixed-anion compounds, which fall in the class of
    their most electronegative anion."""
    _require_pymatgen()
    from pymatgen.core import Composition
    els = list(Composition(composition).elements)
    if len(els) == 1:
        return "element"
    top = max(els, key=lambda e: (e.X if e.X == e.X else -1.0))   # NaN-safe
    for fam, members in ANION_FAMILIES.items():
        if str(top) in members:
            return fam
    return "intermetallic"


def flip_is_correction_driven(energies_scheme_j, energies_scheme_k, target_id, tol=0.0):
    """Relabeling design, release pair (j, k). Inputs are the corrected total
    energies of ALL entries in the target's chemical system at release k,
    computed from the SAME raw energies under release j's and release k's
    energy-correction schemes: dicts {entry_id: (composition, total_energy)}.
    The flip of `target_id` to stable at k counts as correction-driven only if
    the target is stable (dHd <= tol) under scheme k and not stable under
    scheme j. Flips caused by new competing phases or recalculation fail this
    test and are excluded from group F."""
    _require_pymatgen()
    from pymatgen.analysis.phase_diagram import PDEntry
    from pymatgen.core import Composition

    def label(energies):
        ents = {i: PDEntry(Composition(c), e) for i, (c, e) in energies.items()}
        return decomposition_energy(list(ents.values()), ents[target_id]) <= tol

    if set(energies_scheme_j) != set(energies_scheme_k):
        raise ValueError("both schemes must cover the same release-k entries")
    return bool(label(energies_scheme_k) and not label(energies_scheme_j))


def load_snapshot(path: str):
    """Placeholder loader: refuses to run without the frozen snapshot file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"snapshot {path} not found; download and checksum it first (see runbook)")
    raise NotImplementedError("implement the database-specific parser for this snapshot")


# ----------------------------------------------------------------------------
# Self-test on toy structures (identity logic only; no real data)
# ----------------------------------------------------------------------------
def selftest() -> bool:
    _require_pymatgen()
    from pymatgen.analysis.phase_diagram import PDEntry
    from pymatgen.core import Composition, Lattice, Structure
    ok = True

    def check(name, cond, msg):
        nonlocal ok
        ok &= bool(cond)
        print(f"[{'PASS' if cond else 'FAIL'}] {name}: {msg}")

    a = 3.9
    # Ordered L1_0-type CuAu decoration of a conventional fcc cell
    l10 = Structure(Lattice.cubic(a), ["Cu", "Cu", "Au", "Au"],
                    [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    # Slightly strained copy of the same ordered structure (should be T1)
    l10b = Structure(Lattice.tetragonal(a * 1.01, a * 0.99), ["Cu", "Cu", "Au", "Au"],
                     [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    # Disordered fcc Cu0.5Au0.5 parent (should give T2 against the ordered decoration)
    fcc_dis = Structure(Lattice.cubic(a), [{"Cu": 0.5, "Au": 0.5}] * 4,
                        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    # Rock-salt and zinc-blende NaCl: same composition, different structure (T3)
    rs = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    zb = Structure.from_spacegroup("F-43m", Lattice.cubic(6.0), ["Na", "Cl"], [[0, 0, 0], [0.25, 0.25, 0.25]])

    check("T1_same_ordered", outcome_tier(l10, l10b) == "T1", f"{outcome_tier(l10, l10b)}")
    check("T2_order_disorder", outcome_tier(l10, fcc_dis) == "T2", f"{outcome_tier(l10, fcc_dis)}")
    check("T3_polymorph", outcome_tier(rs, zb) == "T3", f"{outcome_tier(rs, zb)}")
    check("None_other_composition", outcome_tier(rs, "NaCl2") is None, f"{outcome_tier(rs, 'NaCl2')}")
    check("T3_composition_only", outcome_tier("Cu2MnAl", "Cu2MnAl") == "T3", "composition-level")
    check("tier_order", meets("T1", "T2") and not meets("T3", "T2"), "T1 > T2 > T3")

    # Decomposition-energy sign convention on a toy binary (synthetic energies)
    ents = [PDEntry(Composition("A"), 0.0), PDEntry(Composition("B"), 0.0),
            PDEntry(Composition("AB"), -2.0),          # -1.0 eV/atom: on hull
            PDEntry(Composition("A3B"), -1.6)]         # -0.4 eV/atom: 0.1 above hull
    dAB = decomposition_energy(ents, ents[2])
    dA3B = decomposition_energy(ents, ents[3])
    # Without AB the hull at x=0.5 is the A3B-B tie line: -0.4 + 0.4/3 = -0.2667 eV/atom,
    # so dHd(AB) = -1.0 - (-0.2667) = -11/15 eV/atom. (An earlier hand expectation of -1.0
    # ignored A3B; the harness caught it. Logged in REDFLAGS.md.)
    check("dHd_stable_negative", abs(dAB - (-11.0 / 15.0)) < 1e-8, f"AB: {dAB:+.4f} eV/atom (expect -0.7333)")
    check("dHd_unstable_equals_ehull", abs(dA3B - 0.1) < 1e-8, f"A3B: {dA3B:+.4f} eV/atom (expect +0.1)")



    # Solid-solution rule (the TaCr2O6 case): rutile Ta(x)Cr(1-x)O2 reported for 0 <= x <= 0.65
    ok_ss = in_solid_solution_range("TaCr2O6", "CrO2", "Ta0.65Cr0.35O2")
    # x = 0.9 lies clearly outside. (x = 2/3, Ta2CrO6, is within the 0.02 L1
    # tolerance of the x = 0.65 end member, so composition alone cannot exclude it;
    # an earlier hand expectation missed this and the test caught it.)
    out_ss = in_solid_solution_range("Ta9CrO20", "CrO2", "Ta0.65Cr0.35O2")
    edge = in_solid_solution_range("Ta2CrO6", "CrO2", "Ta0.65Cr0.35O2")
    check("solid_solution_rule", ok_ss and not out_ss and edge,
          f"TaCr2O6 inside: {ok_ss}; Ta9CrO20 inside: {out_ss}; Ta2CrO6 within tolerance of the end member: {edge}")

    # Chemistry classes
    cls = {c: chemistry_class(c) for c in ("Cu2MnAl", "LiCoO2", "ZnSe", "Na3PS4", "CsPbBr3",
                                           "LiBH4", "LiFePO4", "TiC", "NiSb", "Fe")}
    expect = {"Cu2MnAl": "intermetallic", "LiCoO2": "oxide", "ZnSe": "chalcogenide",
              "Na3PS4": "chalcogenide", "CsPbBr3": "halide", "LiBH4": "hydride",
              "LiFePO4": "oxide", "TiC": "carbide-boride", "NiSb": "pnictide", "Fe": "element"}
    check("chemistry_class", cls == expect, f"{cls}")

    # Correction-driven flip: same raw energies, scheme k lowers AB by 0.3 eV/atom
    base = {"A": ("A", 0.0), "B": ("B", 0.0), "A3B": ("A3B", -1.6), "AB": ("AB", -0.4)}
    sch_k = dict(base); sch_k["AB"] = ("AB", -1.0)          # -0.5 eV/atom: on hull under scheme k
    check("flip_correction_driven", flip_is_correction_driven(base, sch_k, "AB"),
          "AB unstable under scheme j, stable under scheme k")
    check("flip_not_correction_driven", not flip_is_correction_driven(sch_k, sch_k, "AB"),
          "identical schemes: no correction-driven flip")

    # Cohort logic: entry, pre-existence, realization
    preds = [
        {"pred_id": "p1", "release_year": 2015, "stable": True, "composition": "NaCl", "structure": rs},
        {"pred_id": "p2", "release_year": 2016, "stable": False, "composition": "CuAu", "structure": l10},
        {"pred_id": "p2", "release_year": 2018, "stable": True, "composition": "CuAu", "structure": l10},
        {"pred_id": "p3", "release_year": 2019, "stable": True, "composition": "KBr", "structure": None},
    ]
    exps = [
        {"exp_id": "e1", "first_report_year": 1950, "composition": "NaCl", "structure": rs},
        {"exp_id": "e2", "first_report_year": 2021, "composition": "CuAu", "structure": fcc_dis},
    ]
    coh = {r["pred_id"]: r for r in build_cohort(preds, exps, required_tier="T2")}
    check("cohort_preexisting", coh["p1"]["preexisting"] and coh["p1"]["realization_year"] is None,
          "NaCl predicted in 2015 but reported in 1950 -> pre-existing")
    check("cohort_entry_first_stable", coh["p2"]["entry_year"] == 2018.0, "entry at first stable release")
    check("cohort_realized_T2", coh["p2"]["realization_year"] == 2021.0, "CuAu realized via T2 in 2021")
    check("cohort_unrealized", coh["p3"]["realization_year"] is None and not coh["p3"]["preexisting"],
          "KBr: no report -> right-censored")
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    sys.exit(0 if selftest() else 1)
