"""Outcome identity rules (Table 2) and the chemistry class used for blocking.

Tiers (a predicted entry P against an experimental record E):
  T1  structure match: E is ordered and P matches E with the pymatgen StructureMatcher
      (element comparator, primitive cells, volume scaling).
  T2  order-disorder match: E has partial occupancy and P is an ordered decoration of E
      (order-disorder comparator, supercells allowed), with matching composition.
  T3  composition match: same element set and normalized composition within comp_delta
      (L1 distance), any structure.
All cohorts in this study are composition-level, so T3 is the primary tier.
"""
from __future__ import annotations

DEFAULT_TOLS = {"ltol": 0.2, "stol": 0.3, "angle_tol": 5.0, "comp_delta": 0.02}


def _require_pymatgen():
    try:
        import pymatgen  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError("pymatgen is required; install the pinned requirements") from exc


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
    by the T2 structure matcher."""
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
