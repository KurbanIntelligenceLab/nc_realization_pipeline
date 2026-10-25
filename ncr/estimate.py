"""Certified estimates for every analysis, plus the checks that need no outcomes
(X10 window census, E6a premise check) and the validation summaries (X1, X3, X5)."""
from __future__ import annotations

import math
import os
import sys
from collections import defaultdict

import numpy as np

from . import config as C
from .sample import strata_for
from .snapshots import read_snapshot
from .util import decimal_year, log, parse_date, read_csv, read_json, sha256_file, write_csv, write_json

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import nc_realization_metrics as NCM                 # noqa: E402


def _y(s):
    d, _ = parse_date(s)
    return decimal_year(d) if d else math.nan


def _iv(r: dict) -> dict:
    return {"est": float(r["estimate"]), "lo": float(r["ci_low"]), "hi": float(r["ci_high"]), "se": float(r["se"])}


def _clip(d: dict) -> dict:
    """Shares live in [0, 1]; plotted intervals are truncated to that range and flagged."""
    out = {k: min(1.0, max(0.0, d[k])) for k in ("est", "lo", "hi")}
    out["truncated"] = any(abs(out[k] - d[k]) > 1e-12 for k in ("est", "lo", "hi"))
    return out


class Bundle:
    """A labeled population with its expert sample, ready for equation (1)."""

    def __init__(self, cfg, cohort_file, labels_file):
        self.cfg = cfg
        self.rows = {r["pred_id"]: r for r in read_csv(cohort_file)}
        frozen = labels_file + ".sha256"
        if not os.path.exists(frozen) or open(frozen, encoding="utf-8").read().strip() != sha256_file(labels_file):
            raise RuntimeError(f"{labels_file} is not frozen or changed after freezing")
        self.labs = read_csv(labels_file)
        self.name = os.path.splitext(os.path.basename(labels_file))[0]
        d = os.path.join(cfg["workdir"], "experts", self.name)
        key = read_csv(os.path.join(d, "sample_key.csv"))
        fin = {x["sample_id"]: x for x in read_csv(os.path.join(d, "final_labels.csv"))}
        if set(fin) != {k["sample_id"] for k in key}:
            raise RuntimeError("final_labels.csv does not cover the whole sample; finish adjudication")
        self.kappa = read_json(os.path.join(d, "kappa.json"))
        self.strata = strata_for(self.labs)
        self.idx = np.array([int(k["row_index"]) for k in key])
        self.N = len(self.labs)
        R = [self.rows[x["pred_id"]] for x in self.labs]
        self.R = R
        self.e = np.array([_y(r["entry_date"]) for r in R])
        self.m = np.array([_y(x["model_first_date"]) for x in self.labs])
        self.m_rec = np.array([_y(x["model_first_received_date"]) for x in self.labs])
        self.m_men = np.array([_y(x.get("model_first_mention_date", "")) for x in self.labs])
        self.h = np.array([_y(fin[k["sample_id"]]["final_first_report_date"]) for k in key])
        hr = [fin[k["sample_id"]].get("final_first_received_date", "") for k in key]
        self.h_rec = np.array([_y(s) if s else _y(fin[k["sample_id"]]["final_first_report_date"])
                               for s, k in zip(hr, key)])
        self.h_men = np.array([_y(fin[k["sample_id"]].get("final_first_mention_date", "")) for k in key])
        self.e_s = self.e[self.idx]

    @staticmethod
    def pre(first, e):
        return np.where(np.isfinite(first) & (first <= e), 1.0, 0.0)

    @staticmethod
    def within(first, lo, hi):
        return np.where(np.isfinite(first) & (first > lo) & (first <= hi), 1.0, 0.0)

    def cd(self, c, L, H, alpha=0.05):
        return NCM.census_difference(c, L, self.strata, self.idx, H, alpha)

    def errors_by_stratum(self, L, H):
        out = {}
        s = self.strata[self.idx]
        for k in np.unique(s):
            m = s == k
            x = int((L[self.idx][m] != H[m]).sum())
            n = int(m.sum())
            lo, hi = NCM.wilson_ci(x, n)
            out[str(k)] = {"N": int((self.strata == k).sum()), "n": n, "disagreements": x,
                           "rate": x / n, "ci_low": lo, "ci_high": hi}
        return out

    def hajek(self, mask, H):
        s = self.strata[self.idx]
        w = np.array([(self.strata == k).sum() / (s == k).sum() for k in s])
        mm = mask[self.idx]
        return float((w[mm] * H[mm]).sum() / w[mm].sum()) if mm.any() else math.nan


# ----------------------------------------------------------------------------
# E2 / E3 base rates (per source), naive baselines (X2)
# ----------------------------------------------------------------------------
def base_rates(cfg, b: Bundle, source: str, n_pop: int | None = None) -> dict:
    cutoff = decimal_year(C.cutoff(cfg))
    ml = source in cfg.get("ml_sources", [])
    horizons = cfg.get("ml_horizons_years", [2]) if ml else cfg.get("horizons_years", [3, 5, 10])
    alpha_h = 0.05 / len(horizons)
    two_phase = n_pop is not None and n_pop > b.N
    Lpre, Hpre = b.pre(b.m, b.e), b.pre(b.h, b.e_s)
    res = {"source": source, "N_labeled": b.N, "N_population": n_pop or b.N, "two_phase": two_phase,
           "kappa": b.kappa, "horizons": {}, "alpha_per_horizon": alpha_h}

    def share(c_mask, L, H, alpha):
        if two_phase:
            if not c_mask.all():
                raise NotImplementedError("two-phase domain estimates need complete follow-up for the whole "
                                          "phase-1 sample; shorten the horizon or extend the cutoff")
            return _iv(NCM.two_phase_stratified(L, b.strata, b.idx, H, n_pop, alpha))
        return _iv(b.cd(c_mask / c_mask.sum(), L, H, alpha))

    allmask = np.ones(b.N)
    res["pre_existence_all"] = share(allmask, Lpre, Hpre, 0.05)
    res["errors_by_stratum_pre"] = b.errors_by_stratum(Lpre, Hpre)
    for t in horizons:
        lag_sets = {}
        for lag in (1.0, 2.0):
            at = (b.e + t <= cutoff - lag).astype(float)
            lag_sets[lag] = at
        at = lag_sets[1.0]
        if at.sum() < 2:
            res["horizons"][str(t)] = {"note": "fewer than two entries with complete follow-up"}
            continue
        LY, HY = b.within(b.m, b.e, b.e + t), b.within(b.h, b.e_s, b.e_s + t)
        rec = {"n_at_risk": int(at.sum()),
               "realized_share": share(at, LY, HY, alpha_h),
               "pre_existence_within_R": share(at, Lpre, Hpre, alpha_h),
               "not_realized_share": share(at, 1 - Lpre - LY, 1 - Hpre - HY, alpha_h),
               "naive_model_only": float((at * LY).sum() / at.sum()),
               "naive_expert_only_hajek": b.hajek(at.astype(bool), HY),
               "errors_by_stratum_realized": b.errors_by_stratum(LY, HY)}
        try:
            rr = NCM.census_ratio(at, at, LY, 1 - Lpre, b.strata, b.idx, HY, 1 - Hpre, alpha_h)
            if two_phase:
                rr = two_phase_ratio_adjust(rr, b, LY, HY, 1 - Lpre, 1 - Hpre, n_pop, alpha_h)
            rec["realized_among_not_known"] = _iv(rr)
        except ValueError as exc:
            rec["realized_among_not_known"] = {"note": str(exc)}
        if not two_phase:
            LYr, HYr = b.within(b.m_rec, b.e, b.e + t), b.within(b.h_rec, b.e_s, b.e_s + t)
            rec["sensitivity_received_date"] = _iv(b.cd(at / at.sum(), LYr, HYr, alpha_h))
            at2 = lag_sets[2.0]
            if at2.sum() >= 2:
                rec["sensitivity_lag_2yr"] = _iv(b.cd(at2 / at2.sum(), LY, HY, alpha_h))
        db = cfg.get("_db_only", {}).get(b.name)
        if db is not None:
            Ldb = b.within(db, b.e, b.e + t)
            rec["naive_database_only"] = float((at * Ldb).sum() / at.sum())
        res["horizons"][str(t)] = rec
    return res


def two_phase_ratio_adjust(rr: dict, b, LY, HY, LD, HD, n_pop: int, alpha: float) -> dict:
    """Ratio of two shares for a two-phase cohort. census_ratio on the phase-1
    sample is conditionally (approximately) unbiased for the phase-1 ratio; the
    phase-1 sampling variance of a ratio of binary means, Y a subset of D, is
    (1 - n1/N) * pY (1 - R) / (n1 pD^2) by the delta method, with pY and pD
    estimated by equation (1). Checked by simulation in tests/test_pipeline.py."""
    n1 = b.N
    c = np.full(n1, 1.0 / n1)
    pY = min(max(b.cd(c, LY, HY)["estimate"], 0.0), 1.0)
    pD = min(max(b.cd(c, LD, HD)["estimate"], 1e-9), 1.0)
    R = rr["estimate"]
    v1 = (1.0 - n1 / n_pop) * pY * max(1.0 - R, 0.0) / (n1 * pD ** 2) * n_pop / max(n_pop - 1, 1)
    se = math.sqrt(rr["se"] ** 2 + v1)
    z = NCM.zcrit(alpha)
    return {"estimate": R, "se": se, "ci_low": R - z * se, "ci_high": R + z * se}


# ----------------------------------------------------------------------------
# E6 threshold design, X9 attention check, placebo, density and balance
# ----------------------------------------------------------------------------
def _smd(a, b):
    va, vb = np.var(a, ddof=1) if a.size > 1 else 0.0, np.var(b, ddof=1) if b.size > 1 else 0.0
    den = math.sqrt((va + vb) / 2)
    return 0.0 if den == 0 else float((a.mean() - b.mean()) / den)


def _binom_two_sided(k, n):
    """Exact two-sided binomial test of k successes in n trials at p = 1/2, computed in log
    space (math.comb(n, i) overflows a float for n above about 1,030). The p-value sums the
    probabilities of all outcomes no more likely than the observed one."""
    if n == 0:
        return 1.0
    c = math.lgamma(n + 1) + n * math.log(0.5)
    lp = [c - math.lgamma(i + 1) - math.lgamma(n - i + 1) for i in range(n + 1)]
    lk = lp[k]
    return float(min(1.0, sum(math.exp(v) for v in lp if v <= lk + 1e-9)))


def threshold(cfg, b: Bundle) -> dict:
    h = float(cfg.get("threshold", {}).get("horizon_years", 5))
    x = np.array([float(r["x_mev"]) for r in b.R])
    side = np.array([r["side"] for r in b.R])
    block = np.array([r["block"] for r in b.R])
    dis = np.array([r.get("side_disagrees", "0") == "1" for r in b.R])
    LY, HY = b.within(b.m, b.e, b.e + h), b.within(b.h, b.e_s, b.e_s + h)
    LM, HM = b.within(b.m_men, b.e, b.e + h), b.within(b.h_men, b.e_s, b.e_s + h)
    out = {"population": b.N, "kappa": b.kappa, "windows": {}}
    eff_of_int = cfg.get("threshold", {}).get("effect_of_interest")
    for w in cfg.get("threshold", {}).get("windows_mev", [10, 25, 50]):
        inw = np.abs(x) <= w
        A, B = inw & (side == "A"), inw & (side == "B")
        rec = {"N_A": int(A.sum()), "N_B": int(B.sum())}
        try:
            eff = NCM.threshold_effect(LY, b.strata, b.idx, HY, A, B, blocks=block)
            p = min(max(b.cd(NCM.coef_domain_mean(inw), LY, HY)["estimate"], 0.0), 1.0)
            rec.update({"effect": {k: eff[k] for k in ("estimate", "se", "se_certification_only", "ci_low", "ci_high")},
                        "blocks_used": eff["blocks_used"], "blocks_dropped": eff["blocks_dropped"],
                        "mde": NCM.threshold_mde(p, int(A.sum()), int(B.sum()), eff["se_certification_only"] ** 2),
                        "outcome_rate_in_window": p})
            if eff_of_int is not None:
                rec["mde_exceeds_effect_of_interest"] = rec["mde"] > float(eff_of_int)
            keep = ~dis
            e2 = NCM.threshold_effect(LY, b.strata, b.idx, HY, A & keep, B & keep, blocks=block)
            rec["sensitivity_excluding_label_disagreements"] = {k: e2[k] for k in ("estimate", "ci_low", "ci_high")}
            xr = np.array([float(r["x_recomputed_mev"]) if r.get("x_recomputed_mev", "") not in ("", None) else np.nan
                           for r in b.R])
            inr = np.isfinite(xr) & (np.abs(xr) <= w)
            if inr.any():
                e3 = NCM.threshold_effect(LY, b.strata, b.idx, HY, inr & (side == "A"), inr & (side == "B"), blocks=block)
                rec["sensitivity_recomputed_running_variable"] = {k: e3[k] for k in ("estimate", "ci_low", "ci_high")}
            at = NCM.threshold_effect(LM, b.strata, b.idx, HM, A, B, blocks=block)
            rec["attention_check"] = {k: at[k] for k in ("estimate", "se", "ci_low", "ci_high")}
        except ValueError as exc:
            rec["note"] = str(exc)
        out["windows"][str(w)] = rec
    # Balance in the primary window (check b), overall and for the blocked estimator
    inw = np.abs(x) <= 25
    out["balance_primary_window"] = _balance(b.R, inw, side, block)
    # Density (binomial) test within 5 meV/atom on each side
    nl, nr = int(((x > -5) & (x <= 0)).sum()), int(((x > 0) & (x <= 5)).sum())
    out["density_test"] = {"n_left": nl, "n_right": nr, "p_value": _binom_two_sided(nl, nl + nr)}
    # Placebo thresholds at +/-15 meV/atom with 10 meV/atom intervals
    pl = {}
    for c in (-15.0, 15.0):
        left, right = (x > c - 10) & (x <= c), (x > c) & (x <= c + 10)
        try:
            e = NCM.threshold_effect(LY, b.strata, b.idx, HY, left, right, blocks=block)
            pl[str(c)] = {k: e[k] for k in ("estimate", "ci_low", "ci_high")}
        except ValueError as exc:
            pl[str(c)] = {"note": str(exc)}
    out["placebo"] = pl
    # Fig. 3a bins (certified shares per 5 meV/atom bin)
    bins = []
    for lo in range(-50, 50, 5):
        dom = (x > lo) & (x <= lo + 5)
        if dom.sum() == 0:
            continue
        r = _iv(b.cd(NCM.coef_domain_mean(dom), LY, HY))
        bins.append({"center": lo + 2.5, **_clip(r), "n": int(dom.sum())})
    out["fig3a"] = {"window": 25, "bins": bins}
    return out


# ----------------------------------------------------------------------------
# E7 relabeling design
# ----------------------------------------------------------------------------
def relabel(cfg, b: Bundle) -> dict:
    h = float(cfg.get("threshold", {}).get("horizon_years", 5))
    grp = np.array([r["group"] for r in b.R])
    er = np.array([_y(r["release_r_date"]) for r in b.R])
    keys = np.array([r["cem_key"] for r in b.R])
    F, M = grp == "F", grp == "M"
    er_s = er[b.idx]
    LY, HY = b.within(b.m, b.e, b.e + h), b.within(b.h, b.e_s, b.e_s + h)
    LP, HP = b.within(b.m, er, b.e), b.within(b.h, er_s, b.e_s)
    out = {"counts": {str(g): int((grp == g).sum()) for g in np.unique(grp)},
           "minimum_F": int(cfg.get("relabel_min_F", 30))}
    if F.sum() < out["minimum_F"]:
        out["note"] = "group F below the pre-registered minimum; check (c) not assessable"
    for name, (L, H) in {"effect": (LY, HY), "pre_period_placebo": (LP, HP)}.items():
        try:
            e = NCM.threshold_effect(L, b.strata, b.idx, H, F, M, blocks=keys, weights="treated")
            out[name] = {k: e[k] for k in ("estimate", "se", "ci_low", "ci_high", "blocks_used", "blocks_dropped")}
        except ValueError as exc:
            out[name] = {"note": str(exc)}
    return out


# ----------------------------------------------------------------------------
# Checks without outcomes: X10 window census and E6a premise check
# ----------------------------------------------------------------------------
def _balance(rows, inw, side, block) -> dict:
    """Check (b) before any outcome: standardized differences of pre-treatment covariates
    between the sides, overall in the window and within each usable block."""
    covs = {"n_elements": np.array([float(r["n_elements"]) for r in rows]),
            "mean_electronegativity": np.array([float(r["mean_electronegativity"] or "nan") for r in rows]),
            "max_atomic_number": np.array([float(r["max_atomic_number"]) for r in rows])}
    for cls in sorted({r["chemistry_class"] for r in rows}):
        covs[f"class_{cls}"] = np.array([1.0 if r["chemistry_class"] == cls else 0.0 for r in rows])

    def smds(mask):
        A, B = mask & (side == "A"), mask & (side == "B")
        if A.sum() < 2 or B.sum() < 2:
            return None
        return {k: _smd(v[A & np.isfinite(v)], v[B & np.isfinite(v)]) for k, v in covs.items()}
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
    """Block labels for the threshold design. Schemes: 'class' (source x chemistry class, the
    pre-specified default); 'class_en' adds quintiles of mean electronegativity; 'class_en_nel'
    also splits by number of elements. Quintile edges come from the whole threshold population
    (no outcomes), so they are fixed before any outcome is seen."""
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


def step_x10(cfg: dict) -> dict:
    tc = cfg.get("threshold", {})
    p = float(C.need(cfg, "threshold.assumed_outcome_rate"))
    out = {}
    for src in tc.get("sources", []):
        from .cohorts import _slug
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
            mde = NCM.threshold_mde(p, nA, nB) if nA >= 2 and nB >= 2 else None
            rec[str(w)] = {"N_A": nA, "N_B": nB, "entries_in_usable_blocks": usable, "blocks": blocks,
                           "mde_lower_bound_without_certification_variance": mde,
                           "balance": _balance(rows, inw, side, block)}
            if mde is not None and tc.get("effect_of_interest") is not None:
                rec[str(w)]["mde_exceeds_effect_of_interest"] = mde > float(tc["effect_of_interest"])
        edges = en_quintile_edges(rows)
        by_window = {}
        for wv in tc.get("windows_mev", [10, 25, 50]):
            cw = {}
            inw_w = np.abs(x) <= wv
            for scheme in ("class", "class_en", "class_en_nel"):
                bk = block_keys(rows, scheme, edges)
                bal = _balance(rows, inw_w, side, bk)
                cw[scheme] = {"max_abs_blocked_smd": bal["max_abs_blocked_smd"], "passes": bal["check_b_passes_blocked"],
                              "blocked_smd": bal["blocked_smd"]}
            by_window[str(wv)] = cw
        cand = {}
        inw = np.abs(x) <= 25
        for scheme in ("class", "class_en", "class_en_nel"):
            bk = block_keys(rows, scheme, edges)
            bal = _balance(rows, inw, side, bk)
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
                            "mde_lower_bound": NCM.threshold_mde(p, uA, uB) if uA >= 2 and uB >= 2 else None}
        rec["block_scheme_candidates_25meV"] = {"electronegativity_quintile_edges": edges, "schemes": cand}
        rec["block_scheme_balance_by_window"] = by_window
        out[src] = rec
        log(f"x10: {src}: 25 meV window N_A={rec.get('25', {}).get('N_A')}, N_B={rec.get('25', {}).get('N_B')}, "
            f"MDE >= {rec.get('25', {}).get('mde_lower_bound_without_certification_variance')}")
        for sc, c in cand.items():
            log(f"x10: {src}: blocks '{sc}': max |blocked SMD| {c['max_abs_blocked_smd']}, passes {c['passes']}, "
                f"usable A/B {c['N_A_usable']}/{c['N_B_usable']}, MDE >= {c['mde_lower_bound']}")
    return out


def step_premise(cfg: dict) -> dict:
    """E6a: cross-database label disagreement by distance to the threshold (check a).
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
        rows = NCM.label_disagreement_by_bin(dist, sa, sb, edges)
        nmask, far = np.abs(dist) < 25, np.abs(dist) >= 100
        rate = lambda m: float((sa[m] != sb[m]).mean()) if m.any() else None
        lo25 = NCM.wilson_ci(int((sa[nmask] != sb[nmask]).sum()), int(nmask.sum()))[0] if nmask.any() else None
        out[key] = {"matched": len(common), "stable_in_a": int(sa.sum()), "stable_in_b": int(sb.sum()),
                    "bins": rows, "within_25_n": int(nmask.sum()), "within_25_rate": rate(nmask),
                    "within_25_wilson_low": lo25, "beyond_100_rate": rate(far),
                    "check_a_passes": bool(lo25 is not None and lo25 > 0.20 and rate(far) is not None
                                           and rate(far) < rate(nmask))}
        log(f"premise: {key}: matched {len(common)}; within 25 meV disagreement {rate(nmask)}; "
            f"beyond 100 {rate(far)}; check (a) passes = {out[key]['check_a_passes']}")
    return out


# ----------------------------------------------------------------------------
# X1 control panel, X3 residual-miss bound, X5 rerun agreement
# ----------------------------------------------------------------------------
def control_metrics(cfg, cohort_file, labels_file) -> dict:
    rows = {r["pred_id"]: r for r in read_csv(cohort_file)}
    labs = read_csv(labels_file)
    pos, dec, derr = [], [], []
    for x in labs:
        r = rows[x["pred_id"]]
        found = x["model_label"] == "report"
        if r["panel_kind"] in ("positive", "old"):
            pos.append(int(found))
            if found and r.get("known_first_report_date"):
                derr.append(_y(x["model_first_date"]) - _y(r["known_first_report_date"]))
        elif r["panel_kind"] == "decoy":
            dec.append(int(found))
    return NCM.control_panel_metrics(pos, dec, derr or None)


def residual_bound(cfg, b: Bundle, base: dict, horizon: str) -> dict | None:
    d = os.path.join(cfg["workdir"], "experts", b.name)
    p = os.path.join(d, "residual_audit_sheet.csv")
    if not os.path.exists(p):
        return None
    sheet = read_csv(p)
    if any(s["found_missed_report"] == "" for s in sheet):
        return {"note": "residual audit sheet not finished"}
    misses = sum(1 for s in sheet if s["found_missed_report"].strip().lower() in {"yes", "y", "1"})
    rec = base["horizons"].get(horizon, {})
    if "realized_share" not in rec:
        return {"note": f"no realized share at horizon {horizon}"}
    share_neg = 1.0 - rec["realized_share"]["est"] - rec["pre_existence_within_R"]["est"]
    return NCM.residual_miss_bound(rec["realized_share"]["hi"], max(share_neg, 0.0), misses, len(sheet))


def rerun_compare(labels_a: str, labels_b: str) -> dict:
    a = {x["pred_id"]: x for x in read_csv(labels_a)}
    b = {x["pred_id"]: x for x in read_csv(labels_b)}
    ids = sorted(set(a) & set(b))
    la = np.array([a[i]["model_label"] for i in ids])
    lb = np.array([b[i]["model_label"] for i in ids])
    same_date = np.array([a[i]["model_first_date"] == b[i]["model_first_date"] for i in ids])
    return {"n": len(ids), "label_agreement": float((la == lb).mean()), "kappa": NCM.cohen_kappa(la, lb),
            "label_and_date_agreement": float(((la == lb) & same_date).mean())}


# ----------------------------------------------------------------------------
# Design-stage coverage check on KNOWN compounds (no prediction outcomes are read)
# ----------------------------------------------------------------------------
def step_coverage_known(cfg: dict) -> dict:
    """Share of each release's ICSD-origin compositions that also appear in the experimental
    reports table (for example COD). It reads only compounds the source already listed as
    experimental, which are excluded from every cohort, so no prediction outcome is touched.
    Matching is by exact reduced formula, a conservative approximation of the T3 rule."""
    from pymatgen.core import Composition
    import warnings
    warnings.filterwarnings("ignore")
    from .cohorts import read_experiments
    exps = read_experiments(C.resolve(cfg, C.need(cfg, "experiments.path")))
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
        by_n = {}
        for n in sorted(set(comps.values())):
            ks = [k for k, v in comps.items() if v == n]
            by_n[str(n)] = {"n": len(ks), "found": sum(1 for k in ks if k in exp_red),
                            "share": sum(1 for k in ks if k in exp_red) / len(ks)}
        rel_date = r["date"]
        before = sum(1 for k in found if exp_first[k] <= rel_date)
        out[r["release_id"]] = {"icsd_origin_compositions": len(comps), "found_in_reports": len(found),
                                "share_found": len(found) / len(comps), "by_number_of_elements": by_n,
                                "found_with_report_dated_on_or_before_release": before}
        log(f"coverage-known: {r['release_id']}: {len(found)}/{len(comps)} ICSD-origin compositions "
            f"({len(found) / len(comps):.1%}) appear in the experimental reports table")
    return out
