#!/usr/bin/env python3
"""
nc_realization_metrics.py

Reference implementation (NumPy + Python standard library only) of every
statistical quantity used in the manuscript "Fifteen years of computational
materials predictions and their experimental realization".

Design principle
----------------
Every cohort entry i has a frozen language-model (LLM) label L_i in [0, 1] and an
unknown adjudicated expert label H_i in {0, 1}. Experts verify a stratified
simple random sample (without replacement) drawn AFTER the LLM labels are frozen.
Every reported quantity is a linear functional of the expert labels,

        theta_c = sum_i c_i * H_i ,        c_i known for every unit,

for example a cohort realization rate (c_i = 1/N), a rate in a domain such as
one side of the stability threshold (c_i = 1/N_A on the domain), a difference of
two domain rates (threshold design), a fixed-horizon cumulative realization
rate (c_i = 1/|R(t)| on the horizon risk set), or a vector of descriptive
regression coefficients (c_i = row i of (X^T X)^{-1} X^T).

The census difference estimator

        theta_hat = sum_i c_i L_i + sum_h (N_h / n_h) sum_{i in S_h} c_i (H_i - L_i)

is exactly design-unbiased for theta_c for ANY frozen L (Proposition 1 in the
Methods). Its variance is the stratified-total variance of z_i = c_i (H_i - L_i).

Integrity rules enforced in this file
-------------------------------------
* No function fabricates or imputes outcomes. Missing inputs raise.
* Anything produced by the synthetic smoke test is keyed with the prefix
  "SYNTHETIC_" and export_results() refuses to write such keys.
* The smoke test checks statistical PROPERTIES (bias, coverage, identities) of
  the estimators on simulated populations. It never produces a paper result.

Usage
-----
    python3 nc_realization_metrics.py --smoke      # run all self-checks
    python3 nc_realization_metrics.py --design     # print design half-widths
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from statistics import NormalDist

import numpy as np

SYNTH_PREFIX = "SYNTHETIC_"


# ----------------------------------------------------------------------------
# 0. Small utilities
# ----------------------------------------------------------------------------
def zcrit(alpha: float = 0.05) -> float:
    """Two-sided standard normal critical value."""
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    return NormalDist().inv_cdf(1.0 - alpha / 2.0)


def _as1d(x, name: str, dtype=float) -> np.ndarray:
    a = np.asarray(x, dtype=dtype)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if a.size == 0:
        raise ValueError(f"{name} is empty")
    if dtype is float and not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains NaN or inf; missing data must be resolved, not imputed")
    return a


# ----------------------------------------------------------------------------
# 1. Human-only simple random sample (the plan's +/-2.8 pp claim)
# ----------------------------------------------------------------------------
def srswor_halfwidth(p: float, n: int, N: int, alpha: float = 0.05) -> float:
    """Planned half-width of a binary-proportion CI from a simple random sample of
    n units drawn without replacement from N, with NO language-model labels.
    Var(mean) = (N - n)/(N - 1) * p(1 - p)/n."""
    if not (0 < n <= N):
        raise ValueError("need 0 < n <= N")
    fpc = (N - n) / (N - 1) if N > 1 else 0.0
    return zcrit(alpha) * math.sqrt(p * (1.0 - p) / n * fpc)


# ----------------------------------------------------------------------------
# 2. The census difference estimator (Proposition 1)
# ----------------------------------------------------------------------------
def census_difference(c_all, L_all, strata_all, sample_index, H_sample,
                      alpha: float = 0.05, variance: str = "floored"):
    """Design-unbiased estimate of theta = sum_i c_i H_i.

    Parameters
    ----------
    c_all : (N,) or (N, k) known coefficients for every unit (must not depend on H).
    L_all : (N,) frozen LLM labels or scores in [0, 1] for every unit.
    strata_all : (N,) stratum identifiers fixed before sampling (e.g. cohort x LLM label).
    sample_index : (n,) integer indices of the expert-verified units. Within each
        stratum these must be a simple random sample without replacement.
    H_sample : (n,) adjudicated expert labels for those units (0/1).
    variance : "plugin" uses the stratum sample variance of z = c (H - L).
        "floored" (default) additionally floors the per-stratum variance, in strata
        where L is constant and H is binary, at the Laplace-adjusted Bernoulli
        variance p~(1 - p~) * mean(c^2) with p~ = (x + 1)/(n_h + 2). The floor
        prevents the plug-in variance from collapsing when a stratum sample holds
        few or no disagreements. In simulation the plug-in interval undercovered
        badly when disagreements are rare (see scenario B3 in smoke_test) while the
        floored interval stayed at or above nominal; it is conservative when the
        LLM is very accurate (see smoke_test and REDFLAGS.md).

    Returns dict with estimate, se (or covariance for vector c), ci_low, ci_high.
    """
    L_all = _as1d(L_all, "L_all")
    N = L_all.size
    strata_all = np.asarray(strata_all)
    if strata_all.shape != (N,):
        raise ValueError("strata_all must have one entry per unit")
    c_all = np.asarray(c_all, dtype=float)
    vector = c_all.ndim == 2
    if c_all.shape[0] != N:
        raise ValueError("c_all must have one row per unit")
    if not np.all(np.isfinite(c_all)):
        raise ValueError("c_all contains NaN/inf")
    if np.any((L_all < 0) | (L_all > 1)):
        raise ValueError("L_all must lie in [0, 1]")
    sample_index = np.asarray(sample_index, dtype=int)
    H_sample = _as1d(H_sample, "H_sample")
    if sample_index.shape != H_sample.shape:
        raise ValueError("sample_index and H_sample must align")
    if np.unique(sample_index).size != sample_index.size:
        raise ValueError("sample_index contains duplicates (sampling must be without replacement)")
    if variance not in ("plugin", "floored"):
        raise ValueError("variance must be 'plugin' or 'floored'")

    z = zcrit(alpha)
    c2 = c_all if vector else c_all[:, None]
    est = (c2 * L_all[:, None]).sum(axis=0)
    k = c2.shape[1]
    cov = np.zeros((k, k))
    s_strata = strata_all[sample_index]
    for h in np.unique(strata_all):
        in_pop = strata_all == h
        Nh = int(in_pop.sum())
        m = s_strata == h
        nh = int(m.sum())
        if nh < 2:
            raise ValueError(f"stratum {h!r} has {nh} verified units; need at least 2")
        idx = sample_index[m]
        D = H_sample[m] - L_all[idx]
        zz = c2[idx] * D[:, None]                      # (nh, k)
        est = est + (Nh / nh) * zz.sum(axis=0)
        S = np.atleast_2d(np.cov(zz, rowvar=False, ddof=1)) if nh > 1 else np.zeros((k, k))
        if variance == "floored":
            Lh = L_all[in_pop]
            Hh = H_sample[m]
            if np.ptp(Lh) == 0.0 and np.all((Hh == 0) | (Hh == 1)):
                x = Hh.sum()
                pt = (x + 1.0) / (nh + 2.0)
                floor = pt * (1.0 - pt) * np.mean(c2[idx] ** 2, axis=0)
                diag = np.maximum(np.diag(S), floor)
                S = S.copy()
                np.fill_diagonal(S, diag)
        cov += Nh ** 2 * (1.0 - nh / Nh) * S / nh
    if vector:
        se = np.sqrt(np.diag(cov))
        return {"estimate": est, "cov": cov, "se": se,
                "ci_low": est - z * se, "ci_high": est + z * se}
    se = float(math.sqrt(cov[0, 0]))
    e = float(est[0])
    return {"estimate": e, "se": se, "ci_low": e - z * se, "ci_high": e + z * se}


def coef_mean(N: int) -> np.ndarray:
    return np.full(N, 1.0 / N)


def coef_domain_mean(in_domain) -> np.ndarray:
    d = np.asarray(in_domain, dtype=bool)
    nd = int(d.sum())
    if nd == 0:
        raise ValueError("empty domain")
    return d / nd


def coef_domain_difference(in_a, in_b) -> np.ndarray:
    """Coefficients for mean_H(A) - mean_H(B) (threshold design, local randomization)."""
    a = np.asarray(in_a, dtype=bool)
    b = np.asarray(in_b, dtype=bool)
    if np.any(a & b):
        raise ValueError("domains must be disjoint")
    return coef_domain_mean(a) - coef_domain_mean(b)


def coef_ols(X) -> np.ndarray:
    """Rows of (X^T X)^{-1} X^T: theta = sum_i c_i H_i is the finite-population OLS
    coefficient vector of H on X (descriptive linear-probability association)."""
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError("X must be 2-D")
    A = np.linalg.inv(X.T @ X)
    return (A @ X.T).T


def threshold_effect(L_all, strata_all, sample_index, H_sample, in_a, in_b,
                     alpha: float = 0.05, variance: str = "floored", blocks=None,
                     weights: str = "size"):
    """Local effect of the stable label in a threshold window (Assumption 1).
    census_difference() alone certifies the DESCRIPTIVE difference for the fixed
    window population; the causal estimand under as-if random assignment of side
    also carries assignment variance, so the interval adds the Neyman term
    p_A(1-p_A)/(N_A-1) + p_B(1-p_B)/(N_B-1). With `blocks` (e.g. source x
    chemistry class), the effect is the size-weighted average of within-block
    differences, w_b = N_b / sum N_b over blocks with at least 2 units per side;
    other blocks are dropped and reported. weights="treated" uses
    w_b = N_A,b / sum N_A,b instead, the average effect on side A (the relabeling
    design's effect on flipped entries, with CEM strata as blocks). The smoke test
    shows that the certification-only interval undercovers the causal effect."""
    if weights not in ("size", "treated"):
        raise ValueError("weights must be 'size' or 'treated'")
    in_a = np.asarray(in_a, dtype=bool)
    in_b = np.asarray(in_b, dtype=bool)
    N = in_a.size
    blocks = np.zeros(N, dtype=int) if blocks is None else np.asarray(blocks)
    used, dropped = [], []
    for b in np.unique(blocks[in_a | in_b]):
        nA, nB = int((in_a & (blocks == b)).sum()), int((in_b & (blocks == b)).sum())
        (used if (nA >= 2 and nB >= 2) else dropped).append(b)
    if not used:
        raise ValueError("no block has at least 2 units on each side")
    if weights == "size":
        sizes = {b: int(((in_a | in_b) & (blocks == b)).sum()) for b in used}
    else:
        sizes = {b: int((in_a & (blocks == b)).sum()) for b in used}
    tot = sum(sizes.values())
    c = np.zeros(N)
    v_rand = 0.0
    for b in used:
        w = sizes[b] / tot
        A, B = in_a & (blocks == b), in_b & (blocks == b)
        c += w * coef_domain_difference(A, B)
        pa = census_difference(coef_domain_mean(A), L_all, strata_all, sample_index,
                               H_sample, alpha, variance)["estimate"]
        pb = census_difference(coef_domain_mean(B), L_all, strata_all, sample_index,
                               H_sample, alpha, variance)["estimate"]
        pa, pb = min(max(pa, 0.0), 1.0), min(max(pb, 0.0), 1.0)
        v_rand += w * w * (pa * (1 - pa) / (A.sum() - 1) + pb * (1 - pb) / (B.sum() - 1))
    d = census_difference(c, L_all, strata_all, sample_index, H_sample, alpha, variance)
    se = math.sqrt(d["se"] ** 2 + v_rand)
    z = zcrit(alpha)
    e = d["estimate"]
    return {"estimate": e, "se": se, "se_certification_only": d["se"],
            "ci_low": e - z * se, "ci_high": e + z * se,
            "blocks_used": [str(x) for x in used], "blocks_dropped": [str(x) for x in dropped]}


def merge_small_strata(strata_all, min_size: int = 2):
    """Merge strata with fewer than min_size entries into the largest stratum,
    before sampling. Validity of equation (1) does not depend on how strata are
    formed, only on their being fixed before the draw; domains may cut across
    strata."""
    st = np.asarray(strata_all).astype(object)
    labels, counts = np.unique(st, return_counts=True)
    big = labels[np.argmax(counts)]
    small = set(labels[counts < min_size])
    return np.array([big if x in small else x for x in st], dtype=object)


def two_phase_stratified(L_s1, strata_s1, sample_index, H_sample, N_pop: int,
                         alpha: float = 0.05, variance: str = "floored"):
    """Cohort share for cohorts too large to label in full (for example GNoME).
    Phase 1: a simple random sample S1 of n1 entries from the N_pop-entry
    cohort is labeled by the model. Phase 2: a stratified expert sample within
    S1. Equation (1) applied to S1 (c_i = 1/n1) is unbiased for the S1 mean of
    H given S1 (Proposition 1), hence unbiased for the cohort share. Variance:
    (1 - n1/N) S_H^2 / n1 + Var_phase2, with S_H^2 estimated by p(1-p) for a
    binary label."""
    L_s1 = _as1d(L_s1, "L_s1")
    n1 = L_s1.size
    if not (2 <= n1 <= N_pop):
        raise ValueError("need 2 <= n1 <= N_pop")
    c = np.full(n1, 1.0 / n1)
    inner = census_difference(c, L_s1, strata_s1, sample_index, H_sample, alpha, variance)
    p = min(max(inner["estimate"], 0.0), 1.0)
    v1 = (1.0 - n1 / N_pop) * p * (1.0 - p) * N_pop / max(N_pop - 1, 1) / n1
    se = math.sqrt(inner["se"] ** 2 + v1)
    z = zcrit(alpha)
    e = inner["estimate"]
    return {"estimate": e, "se": se, "se_phase2": inner["se"], "ci_low": e - z * se, "ci_high": e + z * se}


def passage_in_document(passage: str, document_text: str, min_chars: int = 40) -> bool:
    """Grounding check for a model-positive label: the quoted passage must occur
    in the retrieved document after normalizing whitespace and case, and must be
    long enough to be specific. A label whose passage fails this check is
    recorded as an abstention, so a model cannot assert a report from memory."""
    norm = lambda t: " ".join(str(t).lower().split())
    p, d = norm(passage), norm(document_text)
    return len(p) >= min_chars and p in d


def census_ratio(c_num, c_den, L_num, L_den, strata_all, sample_index, H_num, H_den,
                 alpha: float = 0.05, variance: str = "floored"):
    """Ratio of two linear functionals, e.g. realization among predictions that
    were not already known: sum_i c_i Y_i / sum_i c_i (1 - P_i). Each total is
    estimated with equation (1); the variance uses the delta method with the
    joint covariance from the same stratified sample. L_num/L_den and
    H_num/H_den are the model and expert labels of the numerator and
    denominator indicators."""
    c_num = np.asarray(c_num, dtype=float)
    c_den = np.asarray(c_den, dtype=float)
    L_num = _as1d(L_num, "L_num")
    L_den = _as1d(L_den, "L_den")
    H_num = _as1d(H_num, "H_num")
    H_den = _as1d(H_den, "H_den")
    N = L_num.size
    if L_den.size != N or c_num.shape != (N,) or c_den.shape != (N,):
        raise ValueError("shape mismatch")
    idx = np.asarray(sample_index, dtype=int)
    est = np.array([(c_num * L_num).sum(), (c_den * L_den).sum()])
    cov = np.zeros((2, 2))
    strata_all = np.asarray(strata_all)
    s_strata = strata_all[idx]
    z = zcrit(alpha)
    for h in np.unique(strata_all):
        Nh = int((strata_all == h).sum())
        m = s_strata == h
        nh = int(m.sum())
        if nh < 2:
            raise ValueError(f"stratum {h!r} has {nh} verified units; need at least 2")
        ii = idx[m]
        zz = np.column_stack([c_num[ii] * (H_num[m] - L_num[ii]), c_den[ii] * (H_den[m] - L_den[ii])])
        est = est + (Nh / nh) * zz.sum(axis=0)
        S = np.cov(zz, rowvar=False, ddof=1)
        if variance == "floored":
            for col, (Lv, Hv, cv) in enumerate(((L_num, H_num[m], c_num), (L_den, H_den[m], c_den))):
                if np.ptp(Lv[strata_all == h]) == 0.0:
                    pt = (Hv.sum() + 1.0) / (nh + 2.0)
                    S[col, col] = max(S[col, col], pt * (1 - pt) * np.mean(cv[ii] ** 2))
        cov += Nh ** 2 * (1.0 - nh / Nh) * S / nh
    num, den = est
    if den <= 0:
        raise ValueError("estimated denominator is not positive")
    R = num / den
    var = (cov[0, 0] - 2 * R * cov[0, 1] + R * R * cov[1, 1]) / den ** 2
    se = math.sqrt(max(var, 0.0))
    return {"estimate": float(R), "se": se, "ci_low": R - z * se, "ci_high": R + z * se}


def cem_weights(treated, keys):
    """Coarsened exact matching for an average effect on the treated.
    treated: (N,) bool; keys: (N,) hashable stratum keys built only from
    pre-treatment database information (e.g. element count, chemistry class,
    dHd bin). Returns coefficients c with sum_i c_i Y_i = mean_T(Y) -
    sum_s (n_T,s / n_T) mean_C,s(Y), using only strata containing both groups,
    and the numbers of treated and control units dropped."""
    t = np.asarray(treated, dtype=bool)
    keys = np.asarray(keys)
    c = np.zeros(t.size)
    usable = [k for k in np.unique(keys) if t[keys == k].any() and (~t[keys == k]).any()]
    nT = int(sum(t[keys == k].sum() for k in usable))
    if nT == 0:
        raise ValueError("no stratum contains both treated and control units")
    for k in usable:
        m = keys == k
        nTs, nCs = int((t & m).sum()), int((~t & m).sum())
        c[t & m] = 1.0 / nT
        c[~t & m] = -(nTs / nT) / nCs
    dropped_T = int(t.sum()) - nT
    dropped_C = int((~t).sum()) - int(sum((~t[keys == k]).sum() for k in usable))
    return c, dropped_T, dropped_C


def weighted_auc(scores, labels, weights) -> float:
    """Design-weighted AUC for a stratified expert sample (weights N_k/n_k).
    An unweighted AUC on a sample that over-represents model positives can be
    biased for the population AUC when model labels correlate with the score."""
    s = _as1d(scores, "scores")
    y = np.asarray(labels).astype(int)
    w = _as1d(weights, "weights")
    pos, neg = y == 1, y == 0
    if not pos.any() or not neg.any():
        raise ValueError("AUC needs both classes")
    sp, wp, sn, wn = s[pos], w[pos], s[neg], w[neg]
    gt = (sp[:, None] > sn[None, :]).astype(float) + 0.5 * (sp[:, None] == sn[None, :])
    return float((wp[:, None] * wn[None, :] * gt).sum() / (wp.sum() * wn.sum()))


def wilson_ci(x: int, n: int, alpha: float = 0.05):
    if n <= 0:
        raise ValueError("n must be positive")
    z = zcrit(alpha)
    p = x / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    lo, hi = max(0.0, mid - half), min(1.0, mid + half)
    if x == 0:
        lo = 0.0          # exact bound; floating-point rounding otherwise leaves about 1e-18
    if x == n:
        hi = 1.0
    return lo, hi


def label_disagreement_by_bin(dist_to_hull, stable_a, stable_b, edges, alpha: float = 0.05):
    """Support check for the threshold design: among compounds computed by two
    databases, the share whose stable labels disagree, binned by |distance to
    hull| (meV/atom, from database A). Complete data, so no certification is
    needed; Wilson intervals per bin."""
    d = np.abs(_as1d(dist_to_hull, "dist_to_hull"))
    a = np.asarray(stable_a, dtype=bool)
    b = np.asarray(stable_b, dtype=bool)
    if not (a.shape == b.shape == d.shape):
        raise ValueError("inputs must align")
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (d >= lo) & (d < hi)
        n = int(m.sum())
        if n == 0:
            out.append({"bin": (lo, hi), "n": 0, "rate": None})
            continue
        x = int((a[m] != b[m]).sum())
        lo_ci, hi_ci = wilson_ci(x, n, alpha)
        out.append({"bin": (lo, hi), "n": n, "rate": x / n, "ci_low": lo_ci, "ci_high": hi_ci})
    return out


# ----------------------------------------------------------------------------
# Measurement validation helpers (experiments X1, X3, power for E6)
# ----------------------------------------------------------------------------
def control_panel_metrics(pred_pos_on_positives, pred_pos_on_decoys,
                          date_error_years=None, alpha: float = 0.05):
    """X1 control panel. Sensitivity on documented prediction-then-synthesis
    cases, specificity on decoy entries that should have no report, and the
    median absolute error of the extracted report date (years). These describe
    the pipeline for allocation and diagnosis only; they do not enter the
    certified estimates, which correct errors on the cohorts themselves."""
    p = np.asarray(pred_pos_on_positives).astype(int)
    d = np.asarray(pred_pos_on_decoys).astype(int)
    if p.size == 0 or d.size == 0:
        raise ValueError("need both positive and decoy panels")
    sens = p.mean(); spec = 1 - d.mean()
    out = {"sensitivity": float(sens), "sensitivity_ci": wilson_ci(int(p.sum()), p.size, alpha),
           "specificity": float(spec), "specificity_ci": wilson_ci(int((1 - d).sum()), d.size, alpha)}
    if date_error_years is not None:
        de = np.abs(_as1d(date_error_years, "date_error_years"))
        out["median_abs_date_error_years"] = float(np.median(de))
    return out


def residual_miss_bound(theta_hat_ci_high: float, share_expert_negative: float,
                        misses: int, n_audit: int, alpha: float = 0.05) -> dict:
    """X3 residual-miss audit. A third expert re-searches a random sample of
    n_audit entries that the adjudicated record labels as not realized and
    finds `misses` reports the panel missed. Because certified estimates are
    unbiased for the ADJUDICATED labels, residual misses bias realization
    downward. Upper bound on the true share: ci_high + share_neg * upper Wilson
    bound of the miss rate (a conservative, not exact, adjustment)."""
    if not (0 <= misses <= n_audit) or n_audit <= 0:
        raise ValueError("invalid audit counts")
    lo, hi = wilson_ci(misses, n_audit, alpha)
    return {"miss_rate": misses / n_audit, "miss_rate_ci": (lo, hi),
            "adjusted_upper_bound": float(theta_hat_ci_high + share_expert_negative * hi)}


def threshold_mde(p: float, n_a: int, n_b: int, cert_var: float = 0.0,
                  alpha: float = 0.05, power: float = 0.80) -> float:
    """Minimum detectable label effect (two-sided test) for the threshold
    design: (z_{1-alpha/2} + z_power) * sqrt(p(1-p)/(n_a-1) + p(1-p)/(n_b-1)
    + certification variance). With alpha = 0.05 and power 0.8 the multiplier
    is 2.80."""
    if n_a < 2 or n_b < 2:
        raise ValueError("need at least 2 entries per side")
    mult = NormalDist().inv_cdf(1 - alpha / 2) + NormalDist().inv_cdf(power)
    return mult * math.sqrt(p * (1 - p) / (n_a - 1) + p * (1 - p) / (n_b - 1) + cert_var)


# ----------------------------------------------------------------------------
# Table formatting conventions (Supplementary Information, "Table conventions")
# ----------------------------------------------------------------------------
def format_share(est: float, lo: float, hi: float, decimals: int = 1) -> str:
    """Format a share as a percentage with its 95% interval. The point estimate is
    rounded half-even; the interval is rounded OUTWARD (lower bound down, upper
    bound up), so rounding never narrows the interval and the rounded estimate
    always lies inside the rounded interval."""
    from decimal import Decimal, ROUND_HALF_EVEN, ROUND_FLOOR, ROUND_CEILING
    if not (lo <= est <= hi):
        raise ValueError("interval must contain the estimate")
    q = Decimal(1).scaleb(-decimals)
    to = lambda x, mode: Decimal(repr(float(100.0 * x))).quantize(q, rounding=mode)
    e, l, h = to(est, ROUND_HALF_EVEN), to(lo, ROUND_FLOOR), to(hi, ROUND_CEILING)
    return f"{e}% [{l}, {h}]"


def compare_rounded(a: float, b: float, decimals: int = 1) -> str:
    """Report the order of two shares at the published precision; equal rounded
    values are 'unresolved at this precision', never a ranking."""
    ra, rb = round(100 * a, decimals), round(100 * b, decimals)
    if ra == rb:
        return "unresolved at this precision"
    return "first larger" if ra > rb else "second larger"


# ----------------------------------------------------------------------------
# 3. Design planning: Neyman allocation and planned half-widths
# ----------------------------------------------------------------------------
def neyman_allocation(N_h, S_h, n_total: int, min_per_stratum: int = 30) -> np.ndarray:
    N_h = np.asarray(N_h, dtype=float)
    S_h = np.asarray(S_h, dtype=float)
    w = N_h * np.maximum(S_h, 1e-6)
    raw = n_total * w / w.sum()
    n_h = np.maximum(np.floor(raw).astype(int), min_per_stratum)
    n_h = np.minimum(n_h, N_h.astype(int))
    return n_h


def planned_halfwidth_llm_strata(N: int, p: float, sens: float, spec: float,
                                 n_total: int, alpha: float = 0.05,
                                 min_per_stratum: int = 30):
    """Planned half-width when strata are the LLM-positive and LLM-negative sets.
    Inside the LLM-positive stratum P(H=1) is the precision; inside the
    LLM-negative stratum it is the false-omission rate."""
    pos = p * sens + (1 - p) * (1 - spec)          # share labeled positive
    prec = p * sens / pos if pos > 0 else 0.0
    fom = p * (1 - sens) / (1 - pos) if pos < 1 else 0.0
    N_h = np.array([N * pos, N * (1 - pos)])
    p_h = np.array([prec, fom])
    S_h = np.sqrt(p_h * (1 - p_h))
    n_h = neyman_allocation(N_h, S_h, n_total, min_per_stratum)
    var = 0.0
    for Nh, nh, ph in zip(N_h, n_h, p_h):
        if nh <= 0:
            continue
        var += (Nh / N) ** 2 * (1 - nh / Nh) * ph * (1 - ph) / nh
    return {"halfwidth": zcrit(alpha) * math.sqrt(var), "n_h": n_h.tolist(),
            "precision": prec, "false_omission": fom, "share_llm_positive": pos}


# ----------------------------------------------------------------------------
# 4. Two-phase (double-sampling) difference estimator for credit shares
# ----------------------------------------------------------------------------
def two_phase_difference(N: int, L_phase1, H_phase2, L_phase2, alpha: float = 0.05):
    """Phase 1: SRSWOR of n1 newly reported compounds out of N in a year, all
    LLM-labeled (L_phase1). Phase 2: SRSWOR of n2 from phase 1, expert-labeled.
    theta_hat = mean(L | S1) + mean(H - L | S2) is unbiased for the population
    share; Var = (1 - n1/N) S_H^2 / n1 + (1 - n2/n1) S_D^2 / n2."""
    L1 = _as1d(L_phase1, "L_phase1")
    H2 = _as1d(H_phase2, "H_phase2")
    L2 = _as1d(L_phase2, "L_phase2")
    n1, n2 = L1.size, H2.size
    if L2.size != n2 or n2 < 2 or n2 > n1 or n1 > N:
        raise ValueError("inconsistent phase sizes")
    D2 = H2 - L2
    est = L1.mean() + D2.mean()
    var = (1 - n1 / N) * H2.var(ddof=1) / n1 + (1 - n2 / n1) * D2.var(ddof=1) / n2
    se = math.sqrt(var)
    z = zcrit(alpha)
    return {"estimate": float(est), "se": se, "ci_low": est - z * se, "ci_high": est + z * se}


def wls_trend(years, estimates, ses):
    """Inverse-variance weighted least-squares slope across independent yearly
    estimates (independent samples per year). Returns slope per year and SE."""
    x = _as1d(years, "years")
    y = _as1d(estimates, "estimates")
    s = _as1d(ses, "ses")
    if np.any(s <= 0):
        raise ValueError("standard errors must be positive")
    w = 1.0 / s ** 2
    xb = (w * x).sum() / w.sum()
    yb = (w * y).sum() / w.sum()
    sxx = (w * (x - xb) ** 2).sum()
    slope = (w * (x - xb) * (y - yb)).sum() / sxx
    return {"slope": float(slope), "se": float(math.sqrt(1.0 / sxx))}


# ----------------------------------------------------------------------------
# 5. Fixed-horizon cumulative realization and pre-existence indicators
# ----------------------------------------------------------------------------
def horizon_indicators(entry_year, first_report_year, horizon: float,
                       cutoff_year: float, lag_buffer: float = 1.0):
    """Risk set R(t): units whose window (entry, entry + t] closes no later
    than cutoff - lag_buffer. Outcome Y_i(t) = 1 if the FIRST matching report
    falls in (entry, entry + t]. Pre-existing units (first report on or before
    entry) stay in the denominator with Y_i(t) = 0: each prediction released as
    new is then exactly one of {already known, realized within t, not realized
    by t}, and every share is a linear functional (Proposition 1). The rate among
    non-pre-existing units is a ratio; use census_ratio(). first_report_year is
    np.nan for 'no report found by cutoff', an observed outcome, not missing data.
    Dates are decimal years."""
    e = _as1d(entry_year, "entry_year")
    f = np.asarray(first_report_year, dtype=float)
    if f.shape != e.shape:
        raise ValueError("shape mismatch")
    at_risk = e + horizon <= cutoff_year - lag_buffer
    y = np.where(np.isfinite(f) & (f > e) & (f <= e + horizon), 1.0, 0.0)
    return at_risk, y


def preexisting_indicator(entry_year, first_report_year):
    e = _as1d(entry_year, "entry_year")
    f = np.asarray(first_report_year, dtype=float)
    return np.where(np.isfinite(f) & (f <= e), 1.0, 0.0)


# ----------------------------------------------------------------------------
# 6. Agreement and discrimination
# ----------------------------------------------------------------------------
def cohen_kappa(a, b) -> float:
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape or a.size == 0:
        raise ValueError("label arrays must align and be non-empty")
    cats = np.unique(np.concatenate([a, b]))
    po = float(np.mean(a == b))
    pe = float(sum(np.mean(a == k) * np.mean(b == k) for k in cats))
    return 1.0 if pe == 1.0 else (po - pe) / (1.0 - pe)


def auc(scores, labels) -> float:
    s = _as1d(scores, "scores")
    y = np.asarray(labels).astype(int)
    pos, neg = s[y == 1], s[y == 0]
    if pos.size == 0 or neg.size == 0:
        raise ValueError("AUC needs both classes")
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    sv = allv[order]
    i = 0
    while i < sv.size:                      # average ranks for ties
        j = i
        while j + 1 < sv.size and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    rp = ranks[:pos.size].sum()
    return float((rp - pos.size * (pos.size + 1) / 2.0) / (pos.size * neg.size))


def bootstrap_auc_ci(scores, labels, B: int = 2000, alpha: float = 0.05, seed: int = 0):
    rng = np.random.default_rng(seed)
    s = _as1d(scores, "scores")
    y = np.asarray(labels).astype(int)
    vals = []
    for _ in range(B):
        idx = rng.integers(0, s.size, s.size)
        if y[idx].min() == y[idx].max():
            continue
        vals.append(auc(s[idx], y[idx]))
    lo, hi = np.quantile(vals, [alpha / 2, 1 - alpha / 2])
    return {"auc": auc(s, y), "ci_low": float(lo), "ci_high": float(hi)}


# ----------------------------------------------------------------------------
# 7. Sampling helper and export guard
# ----------------------------------------------------------------------------
def stratified_srswor(strata_all, n_h: dict, rng) -> np.ndarray:
    """Draw a stratified simple random sample without replacement.
    Record the seed and the hash of the frozen LLM labels before calling."""
    strata_all = np.asarray(strata_all)
    out = []
    for h, nh in n_h.items():
        pool = np.flatnonzero(strata_all == h)
        if nh > pool.size:
            raise ValueError(f"stratum {h!r}: requested {nh} > available {pool.size}")
        out.append(rng.choice(pool, size=nh, replace=False))
    return np.concatenate(out)


def export_results(results: dict, path: str) -> None:
    bad = [k for k in results if str(k).startswith(SYNTH_PREFIX)]
    if bad:
        raise PermissionError(f"refusing to export synthetic results: {bad[:3]}")
    with open(path, "w") as fh:
        json.dump(results, fh, indent=2, default=float)


# ----------------------------------------------------------------------------
# 8. Smoke test (properties only, on simulated populations)
# ----------------------------------------------------------------------------
def _sim_population(N, p, sens, spec, rng):
    H = (rng.random(N) < p).astype(float)
    L = np.where(H == 1, rng.random(N) < sens, rng.random(N) > spec).astype(float)
    return H, L


def _coverage(fn, R, truth):
    hits, ests = 0, []
    for r in range(R):
        out = fn(r)
        ests.append(out["estimate"])
        hits += out["ci_low"] <= truth <= out["ci_high"]
    ests = np.asarray(ests)
    return hits / R, ests.mean() - truth, ests.std(ddof=1) / math.sqrt(R)


def smoke_test(R: int = 2000, verbose: bool = True) -> dict:
    res = {}
    ok = True

    def check(name, cond, detail):
        nonlocal ok
        ok &= bool(cond)
        res[SYNTH_PREFIX + name] = detail
        if verbose:
            print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")

    # A. The plan's human-only figure
    hw = srswor_halfwidth(0.10, 400, 3231)
    check("plan_pm2.8pp_is_human_only", 0.027 <= hw <= 0.028,
          f"human-only half-width at p=0.10, n=400, N=3231 is {100*hw:.2f} pp")

    # B. Population mean, imperfect LLM, design-based coverage (population fixed)
    rng = np.random.default_rng(20261024)
    N = 3000
    H, L = _sim_population(N, p=0.08, sens=0.85, spec=0.97, rng=rng)
    strata = np.where(L == 1, "llm_pos", "llm_neg")
    truth = H.mean()
    n_h = {"llm_pos": 150, "llm_neg": 250}
    c = coef_mean(N)

    for vmode in ("plugin", "floored"):
        def run(r, vmode=vmode):
            g = np.random.default_rng(1000 + r)
            idx = stratified_srswor(strata, n_h, g)
            return census_difference(c, L, strata, idx, H[idx], variance=vmode)
        cov, bias, mcse = _coverage(run, R, truth)
        if vmode == "floored":
            check("mean_unbiased", abs(bias) < 3.5 * mcse, f"bias={bias:+.5f} (MC SE {mcse:.5f})")
            check("mean_coverage_floored", cov >= 0.94, f"coverage={cov:.3f}")
        else:
            res[SYNTH_PREFIX + "mean_coverage_plugin"] = f"coverage={cov:.3f} (reported, not gated)"
            if verbose:
                print(f"[INFO] mean_coverage_plugin: coverage={cov:.3f}")

    # B2. Rare-disagreement stress test: very accurate LLM, small strata
    H2, L2 = _sim_population(N, p=0.05, sens=0.97, spec=0.995, rng=np.random.default_rng(7))
    strata2 = np.where(L2 == 1, "llm_pos", "llm_neg")
    truth2 = H2.mean()
    n2 = {"llm_pos": 60, "llm_neg": 200}
    for vmode in ("plugin", "floored"):
        def run2(r, vmode=vmode):
            g = np.random.default_rng(5000 + r)
            idx = stratified_srswor(strata2, n2, g)
            return census_difference(coef_mean(N), L2, strata2, idx, H2[idx], variance=vmode)
        cov2, _, _ = _coverage(run2, R, truth2)
        if vmode == "floored":
            check("rare_disagreement_coverage_floored", cov2 >= 0.94,
                  f"coverage={cov2:.3f} (conservative when disagreements are rare)")
        else:
            res[SYNTH_PREFIX + "rare_disagreement_coverage_plugin"] = f"coverage={cov2:.3f}"
            if verbose:
                print(f"[INFO] rare_disagreement_coverage_plugin: coverage={cov2:.3f} "
                      "(plugin variance collapses when a stratum sample has no disagreements)")

    # B3. Rare-disagreement scenario with a small false-omission rate
    H4, L4 = _sim_population(N, p=0.02, sens=0.90, spec=0.99, rng=np.random.default_rng(4))
    strata4 = np.where(L4 == 1, "llm_pos", "llm_neg")
    truth4 = H4.mean()
    n4 = {"llm_pos": 60, "llm_neg": 300}
    for vmode in ("plugin", "floored"):
        def run4(r, vmode=vmode):
            g = np.random.default_rng(400 + r)
            idx = stratified_srswor(strata4, n4, g)
            return census_difference(coef_mean(N), L4, strata4, idx, H4[idx], variance=vmode)
        cov4, _, _ = _coverage(run4, R, truth4)
        if vmode == "floored":
            check("rare_omission_coverage_floored", cov4 >= 0.94, f"coverage={cov4:.3f}")
        else:
            res[SYNTH_PREFIX + "rare_omission_coverage_plugin"] = f"coverage={cov4:.3f}"
            if verbose:
                print(f"[INFO] rare_omission_coverage_plugin: coverage={cov4:.3f} "
                      "(plug-in interval far below nominal)")

    # C. Domain difference (threshold design): two sides of a window
    side = rng.random(N)
    in_a = side < 0.3
    in_b = (side >= 0.3) & (side < 0.6)
    Hd = H.copy()
    Hd[in_a] = (rng.random(in_a.sum()) < 0.12).astype(float)   # different rates by side
    Ld = np.where(Hd == 1, rng.random(N) < 0.85, rng.random(N) > 0.97).astype(float)
    strata_d = np.where(Ld == 1, "pos", "neg")
    strata_d = np.char.add(strata_d, np.where(in_a, "_A", np.where(in_b, "_B", "_O")))
    truth_d = Hd[in_a].mean() - Hd[in_b].mean()
    cd = coef_domain_difference(in_a, in_b)
    sizes = {h: int((strata_d == h).sum()) for h in np.unique(strata_d)}
    nd = {h: min(sizes[h], 40 if "_O" in h else 90) for h in sizes}

    def run_d(r):
        g = np.random.default_rng(9000 + r)
        idx = stratified_srswor(strata_d, nd, g)
        return census_difference(cd, Ld, strata_d, idx, Hd[idx])
    covd, biasd, mcsed = _coverage(run_d, R, truth_d)
    check("domain_difference_coverage", covd >= 0.94 and abs(biasd) < 3.5 * mcsed,
          f"coverage={covd:.3f}, bias={biasd:+.5f} (MC SE {mcsed:.5f})")

    # C2. Causal threshold effect: side assigned at random, population regenerated
    tau, pB, Nw = 0.04, 0.06, 1500
    hits_c, hits_m = 0, 0
    for r in range(R):
        g = np.random.default_rng(21000 + r)
        A = g.random(Nw) < 0.5
        Hw = (g.random(Nw) < np.where(A, pB + tau, pB)).astype(float)
        Lw = np.where(Hw == 1, g.random(Nw) < 0.85, g.random(Nw) > 0.97).astype(float)
        sw = np.char.add(np.where(Lw == 1, "pos", "neg"), np.where(A, "_A", "_B"))
        sizes = {h: int((sw == h).sum()) for h in np.unique(sw)}
        nw = {h: min(sizes[h], 400) for h in sizes}   # large expert sample: randomization variance dominates
        idx = stratified_srswor(sw, nw, g)
        o = threshold_effect(Lw, sw, idx, Hw[idx], A, ~A)
        z = zcrit(0.05)
        hits_c += o["ci_low"] <= tau <= o["ci_high"]
        hits_m += (o["estimate"] - z * o["se_certification_only"] <= tau <= o["estimate"] + z * o["se_certification_only"])
    check("threshold_causal_coverage", hits_c / R >= 0.94,
          f"coverage of the causal effect={hits_c / R:.3f}")
    res[SYNTH_PREFIX + "threshold_certification_only_coverage"] = f"{hits_m / R:.3f}"
    if verbose:
        print(f"[INFO] threshold_certification_only_coverage: {hits_m / R:.3f} "
              "(measurement-only interval ignores randomization variance)")

    # C3. Blocked threshold effect: side probability and base rate differ by block
    tau3, hits_b, est_b, est_u = 0.04, 0, [], []
    for r in range(R):
        g = np.random.default_rng(29000 + r)
        blk = g.integers(0, 2, 2000)
        pside = np.where(blk == 0, 0.3, 0.7)
        base = np.where(blk == 0, 0.02, 0.15)
        A = g.random(2000) < pside
        Hb = (g.random(2000) < base + tau3 * A).astype(float)
        Lb = np.where(Hb == 1, g.random(2000) < 0.85, g.random(2000) > 0.97).astype(float)
        sb = np.char.add(np.char.add(np.where(Lb == 1, "pos", "neg"), np.where(A, "_A", "_B")), blk.astype(str))
        sizes = {h: int((sb == h).sum()) for h in np.unique(sb)}
        idx = stratified_srswor(sb, {h: min(sizes[h], 120) for h in sizes}, g)
        o = threshold_effect(Lb, sb, idx, Hb[idx], A, ~A, blocks=blk)
        hits_b += o["ci_low"] <= tau3 <= o["ci_high"]
        est_b.append(o["estimate"])
        est_u.append(threshold_effect(Lb, sb, idx, Hb[idx], A, ~A)["estimate"])
    check("blocked_threshold_coverage", hits_b / R >= 0.94 and abs(np.mean(est_b) - tau3) < 0.004,
          f"coverage={hits_b / R:.3f}, mean estimate={np.mean(est_b):.4f} (true {tau3})")
    res[SYNTH_PREFIX + "unblocked_threshold_mean"] = f"{np.mean(est_u):.4f}"
    if verbose:
        print(f"[INFO] unblocked threshold mean estimate: {np.mean(est_u):.4f} "
              f"(true {tau3}; blocking removes confounding by block)")

    # C3b. Treated-weighted (relabeling) effect with block-varying effects
    hits_t, est_t = 0, []
    for r in range(min(R, 1500)):
        g = np.random.default_rng(31000 + r)
        blk = g.integers(0, 3, 2400)
        pA = np.array([0.1, 0.3, 0.6])[blk]
        A = g.random(2400) < pA
        eff = np.array([0.02, 0.05, 0.10])[blk]
        base = np.array([0.05, 0.10, 0.15])[blk]
        Ht = (g.random(2400) < base + eff * A).astype(float)
        Lt = np.where(Ht == 1, g.random(2400) < 0.85, g.random(2400) > 0.97).astype(float)
        st = merge_small_strata(np.char.add(np.char.add(np.where(Lt == 1, "p", "n"), np.where(A, "A", "B")), blk.astype(str)))
        sizes_t = {h: int((st == h).sum()) for h in np.unique(st)}
        idx = stratified_srswor(st, {h: min(sizes_t[h], 80) for h in sizes_t}, g)
        o = threshold_effect(Lt, st, idx, Ht[idx], A, ~A, blocks=blk, weights="treated")
        # superpopulation ATT target: block effects weighted by expected treated counts
        att_real = float(np.sum([eff[(blk == b) & A].sum() for b in range(3)]) / A.sum())
        hits_t += o["ci_low"] <= att_real <= o["ci_high"]
        est_t.append(o["estimate"] - att_real)
    check("treated_weighted_coverage", hits_t / min(R, 1500) >= 0.94 and abs(np.mean(est_t)) < 0.004,
          f"coverage={hits_t / min(R, 1500):.3f}, mean error={np.mean(est_t):+.4f}")

    # C4. CEM coefficients recover a known effect on the treated exactly
    keys = np.array(["a", "a", "a", "b", "b", "c"])
    trt = np.array([True, False, False, True, False, True])     # stratum c has no control
    base_k = {"a": 0.1, "b": 0.5, "c": 0.9}
    Yc = np.array([base_k[k] for k in keys]) + 0.2 * trt
    cc, dT, dC = cem_weights(trt, keys)
    check("cem_exact", abs(cc @ Yc - 0.2) < 1e-12 and dT == 1 and dC == 0,
          f"CEM difference={cc @ Yc:.3f} (true 0.2), dropped treated={dT}")

    # C5. Design-weighted AUC on a stratified expert sample
    ga = np.random.default_rng(31)
    Na = 4000
    Ha = (ga.random(Na) < 0.1).astype(float)
    score = Ha * 0.8 + ga.normal(0, 0.6, Na)
    hi = score > 0.8          # among negatives, the model errs more often for high-score entries
    La = (ga.random(Na) < np.where(Ha == 1, 0.85, 0.02 + 0.25 * hi)).astype(float)
    sa = np.where(La == 1, "pos", "neg")
    true_auc = auc(score, Ha)
    wa_list, ua_list = [], []
    Nh_a = {h: int((sa == h).sum()) for h in ("pos", "neg")}
    nh_a = {"pos": 200, "neg": 200}
    for r in range(min(R, 500)):
        g = np.random.default_rng(33000 + r)
        idx = stratified_srswor(sa, nh_a, g)
        w = np.array([Nh_a[h] / nh_a[h] for h in sa[idx]])
        wa_list.append(weighted_auc(score[idx], Ha[idx], w))
        ua_list.append(auc(score[idx], Ha[idx]))
    check("weighted_auc_unbiased", abs(np.mean(wa_list) - true_auc) < 0.01,
          f"weighted={np.mean(wa_list):.3f}, population={true_auc:.3f}")
    res[SYNTH_PREFIX + "unweighted_auc_mean"] = f"{np.mean(ua_list):.3f}"
    if verbose:
        print(f"[INFO] unweighted AUC on the stratified sample: {np.mean(ua_list):.3f} "
              f"(population {true_auc:.3f})")

    # D. Regression coefficients (vector c)
    X = np.column_stack([np.ones(N), rng.normal(size=N), (rng.random(N) < 0.4).astype(float)])
    lin = 0.06 + 0.03 * X[:, 1] + 0.05 * X[:, 2]
    Hr = (rng.random(N) < np.clip(lin, 0.001, 0.999)).astype(float)
    Lr = np.where(Hr == 1, rng.random(N) < 0.85, rng.random(N) > 0.96).astype(float)
    sr = np.where(Lr == 1, "pos", "neg")
    Cr = coef_ols(X)
    truth_r = Cr.T @ Hr
    nr = {"pos": 150, "neg": 250}
    hits = np.zeros(3)
    for r in range(R):
        g = np.random.default_rng(13000 + r)
        idx = stratified_srswor(sr, nr, g)
        o = census_difference(Cr, Lr, sr, idx, Hr[idx])
        hits += (o["ci_low"] <= truth_r) & (truth_r <= o["ci_high"])
    covr = hits / R
    check("ols_coefficient_coverage", bool(np.all(covr >= 0.94)),
          f"coverage per coefficient={np.round(covr, 3).tolist()}")

    # E. Two-phase credit-share estimator (population fixed, both phases resampled)
    Nc = 4000
    Hc = (rng.random(Nc) < 0.07).astype(float)
    Lc = np.where(Hc == 1, rng.random(Nc) < 0.8, rng.random(Nc) > 0.95).astype(float)
    truth_c = Hc.mean()

    def run_c(r):
        g = np.random.default_rng(17000 + r)
        s1 = g.choice(Nc, 600, replace=False)
        s2 = g.choice(s1, 120, replace=False)
        return two_phase_difference(Nc, Lc[s1], Hc[s2], Lc[s2])
    covc, biasc, mcsec = _coverage(run_c, R, truth_c)
    check("two_phase_coverage", 0.92 <= covc <= 0.975 and abs(biasc) < 3.5 * mcsec,
          f"coverage={covc:.3f}, bias={biasc:+.5f}")

    # F. Fixed-horizon logic (pre-existing units stay in the denominator with Y = 0)
    e = np.array([2015., 2015., 2020., 2024., 2015.])
    f = np.array([2018., np.nan, 2021., 2025., 2014.])
    at, y = horizon_indicators(e, f, horizon=5, cutoff_year=2026, lag_buffer=1)
    pre = preexisting_indicator(e, f)
    check("horizon_logic", at.tolist() == [True, True, True, False, True]
          and y.tolist() == [1, 0, 1, 1, 0] and pre.tolist() == [0, 0, 0, 0, 1],
          f"at_risk={at.tolist()}, y={y.tolist()}, preexisting={pre.tolist()}")

    # F1. Ratio estimator: realization among predictions not already known
    gr = np.random.default_rng(77)
    Nq = 3000
    P = (gr.random(Nq) < 0.15).astype(float)                       # already known
    Y = np.where(P == 0, (gr.random(Nq) < 0.08), 0).astype(float)   # realized after release
    anyrep = np.maximum(P, Y)
    Lany = np.where(anyrep == 1, gr.random(Nq) < 0.85, gr.random(Nq) > 0.97).astype(float)
    # model's timing call is wrong for 10% of its positives
    swap = (Lany == 1) & (gr.random(Nq) < 0.10)
    LP = np.where(Lany == 1, np.where(swap, 1 - P, P), 0.0)
    LY = np.where(Lany == 1, 1 - LP, 0.0)
    sq = np.where(Lany == 1, "pos", "neg")
    truth_q = Y.sum() / (1 - P).sum()
    cN = np.full(Nq, 1.0 / Nq)
    hits_q, est_q = 0, []
    for r in range(R):
        g = np.random.default_rng(23000 + r)
        idx = stratified_srswor(sq, {"pos": 150, "neg": 250}, g)
        o = census_ratio(cN, cN, LY, 1 - LP, sq, idx, Y[idx], 1 - P[idx])
        hits_q += o["ci_low"] <= truth_q <= o["ci_high"]
        est_q.append(o["estimate"])
    bias_q = np.mean(est_q) - truth_q
    check("ratio_coverage", hits_q / R >= 0.94 and abs(bias_q) < 0.002,
          f"coverage={hits_q / R:.3f}, bias={bias_q:+.5f} (ratio estimator, delta method)")

    # F2. Label-disagreement check    # F2. Label-disagreement check: two databases with independent calculation noise
    gg = np.random.default_rng(31)
    true_e = gg.uniform(-150, 150, 20000)                  # meV/atom
    ea = true_e + gg.normal(0, 30, true_e.size)
    eb = true_e + gg.normal(0, 30, true_e.size)
    bins = label_disagreement_by_bin(ea, ea <= 0, eb <= 0, [0, 25, 50, 100, 150])
    rates = [bb["rate"] for bb in bins]
    check("label_disagreement_decreases", rates[0] > rates[-1] and rates[0] > 0.2,
          f"disagreement by |dist| bin = {[round(r, 3) for r in rates]}")

    # G. Agreement and AUC sanity
    a = (rng.random(500) < 0.3).astype(int)
    check("kappa_identity", abs(cohen_kappa(a, a) - 1) < 1e-12, "kappa(a, a) = 1")
    b = (rng.random(500) < 0.3).astype(int)
    check("kappa_independent", abs(cohen_kappa(a, b)) < 0.15, f"kappa(independent)={cohen_kappa(a, b):+.3f}")
    check("auc_perfect", auc([0.1, 0.2, 0.8, 0.9], [0, 0, 1, 1]) == 1.0, "AUC perfect = 1")
    check("auc_ties", auc([0.5, 0.5], [0, 1]) == 0.5, "AUC all-tied = 0.5")

    # I. Exact first-principles checks (full enumeration)
    ex = exact_checks(verbose=False)
    check("exact_enumeration", all(ex.values()), f"{sum(ex.values())}/{len(ex)} exact identities hold: {[k for k,v in ex.items() if not v]}")

    # J. Measurement-validation helpers
    cp = control_panel_metrics([1] * 17 + [0] * 3, [0] * 98 + [1] * 2, date_error_years=[0, 0.1, 0.5, 1.0])
    check("control_panel", abs(cp["sensitivity"] - 0.85) < 1e-12 and abs(cp["specificity"] - 0.98) < 1e-12
          and cp["sensitivity_ci"][0] < 0.85 < cp["sensitivity_ci"][1],
          f"sens={cp['sensitivity']:.2f} {tuple(round(x, 3) for x in cp['sensitivity_ci'])}, spec={cp['specificity']:.2f}")
    rb = residual_miss_bound(0.08, 0.9, 2, 100)
    check("residual_miss_bound", rb["adjusted_upper_bound"] > 0.08 and rb["miss_rate"] == 0.02,
          f"miss rate 0.02, CI {tuple(round(x, 3) for x in rb['miss_rate_ci'])}, adjusted upper {rb['adjusted_upper_bound']:.3f}")
    m1 = threshold_mde(0.05, 500, 500); m2 = threshold_mde(0.05, 2000, 2000)
    mult = NormalDist().inv_cdf(0.975) + NormalDist().inv_cdf(0.8)
    check("threshold_mde", m2 < m1 and abs(mult - 2.8016) < 1e-3,
          f"MDE p=0.05: n=500/side {100*m1:.2f} pp, n=2000/side {100*m2:.2f} pp; multiplier {mult:.3f}")

    # K. Two-phase sampling for very large cohorts
    Nbig = 200000
    gpop = np.random.default_rng(55)
    Hbig = (gpop.random(Nbig) < 0.04).astype(float)
    Lbig = np.where(Hbig == 1, gpop.random(Nbig) < 0.85, gpop.random(Nbig) > 0.97).astype(float)
    truth_big = Hbig.mean()
    hits_tp, est_tp = 0, []
    Rtp = min(R, 1000)
    for r in range(Rtp):
        g = np.random.default_rng(41000 + r)
        s1 = g.choice(Nbig, 4000, replace=False)
        st1 = np.where(Lbig[s1] == 1, "pos", "neg")
        idx2 = stratified_srswor(st1, {"pos": 120, "neg": 280}, g)
        o = two_phase_stratified(Lbig[s1], st1, idx2, Hbig[s1][idx2], Nbig)
        hits_tp += o["ci_low"] <= truth_big <= o["ci_high"]
        est_tp.append(o["estimate"])
    check("two_phase_stratified_coverage", hits_tp / Rtp >= 0.94 and abs(np.mean(est_tp) - truth_big) < 0.002,
          f"coverage={hits_tp / Rtp:.3f}, bias={np.mean(est_tp) - truth_big:+.5f} (n1=4000 of N=200000)")

    # L. Grounding check for model-positive labels
    doc = "We report the solid-state synthesis of  Na3SbS4 and its   structure refined from X-ray diffraction data."
    ok_g = passage_in_document("the solid-state synthesis of Na3SbS4 and its structure refined", doc)
    bad_g = passage_in_document("the first synthesis of Na3SbS4 was reported in 2019 by our group", doc)
    check("grounding_check", ok_g and not bad_g, f"true passage found: {ok_g}; invented passage found: {bad_g}")

    # M. Table conventions: outward rounding keeps the estimate inside its interval
    gr = np.random.default_rng(99)
    ok_fmt = True
    for _ in range(2000):
        lo_, hi_ = sorted(gr.random(2))
        es_ = lo_ + (hi_ - lo_) * gr.random()
        txt = format_share(es_, lo_, hi_)
        e_, rest = txt.split("% [")
        l_, h_ = rest.rstrip("]").split(", ")
        ok_fmt &= float(l_) <= float(e_) <= float(h_) and float(l_) <= 100 * lo_ + 1e-9 and float(h_) >= 100 * hi_ - 1e-9
    check("format_share_outward", ok_fmt and format_share(0.12345, 0.10011, 0.14449) == "12.3% [10.0, 14.5]",
          f"example: {format_share(0.12345, 0.10011, 0.14449)}; 2000 random intervals contain their estimates and are never narrowed")
    check("compare_rounded_ties", compare_rounded(0.10049, 0.10001) == "unresolved at this precision"
          and compare_rounded(0.2, 0.1) == "first larger", "ties at one decimal are reported as unresolved")

    # H. Export guard
    blocked = False
    try:
        export_results({SYNTH_PREFIX + "x": 1}, "/tmp/should_not_exist.json")
    except PermissionError:
        blocked = True
    check("export_guard", blocked, "synthetic keys cannot be exported")

    res[SYNTH_PREFIX + "ALL_PASS"] = ok
    if verbose:
        print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return res


def compare_floor_variants(R: int = 1500, verbose: bool = True) -> dict:
    """Reproduces the choice of variance floor (Supplementary Note 4): coverage of
    nominal 95% intervals for a cohort share under four floors p~ = (x+a)/(n+b),
    in five simulated scenarios. 'plugin' means no floor."""
    variants = {"plugin": None, "jeffreys(0.5/1)": (0.5, 1.0),
                "laplace(1/2)": (1.0, 2.0), "agresti(1.92/3.84)": (1.92, 3.84)}
    scenarios = [("p.08 s.85 sp.97", 0.08, 0.85, 0.97, {"pos": 150, "neg": 250}),
                 ("p.05 s.97 sp.995", 0.05, 0.97, 0.995, {"pos": 60, "neg": 200}),
                 ("p.10 s.70 sp.90", 0.10, 0.70, 0.90, {"pos": 150, "neg": 250}),
                 ("p.02 s.90 sp.99", 0.02, 0.90, 0.99, {"pos": 60, "neg": 300})]
    out = {}
    z = zcrit(0.05)
    for si, (name, p, se_, sp, nh) in enumerate(scenarios):
        rng = np.random.default_rng(100 + si)
        N = 3000
        H, L = _sim_population(N, p, se_, sp, rng)
        st = np.where(L == 1, "pos", "neg")
        truth = H.mean()
        row = {}
        for vn, ab in variants.items():
            hits = 0
            for r in range(R):
                g = np.random.default_rng(10000 * (si + 1) + r)
                idx = stratified_srswor(st, nh, g)
                est = L.mean()
                var = 0.0
                for h in ("pos", "neg"):
                    m = st[idx] == h
                    Nh, n = int((st == h).sum()), int(m.sum())
                    D = H[idx][m] - L[idx][m]
                    est += (Nh / n) * D.sum() / N
                    s2 = D.var(ddof=1) / N ** 2
                    if ab is not None:
                        pt = (H[idx][m].sum() + ab[0]) / (n + ab[1])
                        s2 = max(s2, pt * (1 - pt) / N ** 2)
                    var += Nh ** 2 * (1 - n / Nh) * s2 / n
                hits += abs(est - truth) <= z * math.sqrt(var)
            row[vn] = round(hits / R, 3)
        out[name] = row
        if verbose:
            print(name, row)
    return out


def exact_checks(verbose: bool = True) -> dict:
    """First-principles verification by FULL ENUMERATION (no Monte Carlo):
    1. Proposition 1: E[theta_hat] = theta and E[v_hat] = Var(theta_hat) exactly,
       for scalar and vector weights, continuous model labels, unequal strata.
    2. Equation (2) identity: with c = 1/N and L constant within strata, equation
       (1) equals sum_k (N_k/N) * mean_{S_k}(H) for every sample.
    3. Two-phase estimator: unbiased, and the stated variance formula is exact.
    4. Blocked difference in means under complete randomization within blocks:
       unbiased for the block-size-weighted ATE; the Neyman variance estimator is
       conservative in expectation.
    5. Local-randomization bias with a linear outcome gradient and a uniform
       running variable: bias = -slope * w (side A below the threshold).
    Returns a dict of booleans."""
    from itertools import combinations
    rng = np.random.default_rng(2026)
    out = {}

    # ---- 1. Proposition 1 by enumeration
    sizes = (5, 6, 4)
    n_s = (2, 3, 2)
    N = sum(sizes)
    strata = np.concatenate([[k] * m for k, m in enumerate(sizes)])
    H = (rng.random(N) < 0.4).astype(float)
    L = rng.random(N)                                  # continuous model scores
    c = rng.normal(size=N)                             # arbitrary scalar weights
    C2 = rng.normal(size=(N, 2))                       # vector weights
    theta, theta2 = float(c @ H), C2.T @ H
    pools = [np.flatnonzero(strata == k) for k in range(len(sizes))]
    combos = [list(combinations(p, n)) for p, n in zip(pools, n_s)]
    ests, vars_, ests2, covs2 = [], [], [], []
    import itertools
    for pick in itertools.product(*combos):
        idx = np.array([i for grp in pick for i in grp])
        o = census_difference(c, L, strata, idx, H[idx], variance="plugin")
        ests.append(o["estimate"]); vars_.append(o["se"] ** 2)
        o2 = census_difference(C2, L, strata, idx, H[idx], variance="plugin")
        ests2.append(o2["estimate"]); covs2.append(o2["cov"])
    ests, vars_ = np.array(ests), np.array(vars_)
    out["prop1_unbiased"] = abs(ests.mean() - theta) < 1e-10
    out["prop1_variance_unbiased"] = abs(vars_.mean() - ests.var()) < 1e-10
    E2 = np.array(ests2)
    out["prop1_vector_unbiased"] = bool(np.allclose(E2.mean(axis=0), theta2, atol=1e-10))
    out["prop1_vector_cov_unbiased"] = bool(np.allclose(np.mean(covs2, axis=0),
                                                        np.cov(E2, rowvar=False, ddof=0), atol=1e-10))

    # ---- 2. Equation (2) identity
    Lb = np.array([float(k == 0) for k in strata])      # constant within strata
    ok = True
    for pick in itertools.islice(itertools.product(*combos), 50):
        idx = np.array([i for grp in pick for i in grp])
        e1 = census_difference(np.full(N, 1.0 / N), Lb, strata, idx, H[idx], variance="plugin")["estimate"]
        e2 = sum(len(p) / N * H[np.array(grp)].mean() for p, grp in zip(pools, pick))
        ok &= abs(e1 - e2) < 1e-12
    out["eq2_identity"] = bool(ok)

    # ---- 3. Two-phase estimator by enumeration
    Np, n1, n2 = 7, 4, 2
    Hp = (rng.random(Np) < 0.5).astype(float)
    Lp = rng.random(Np)
    tp_est, tp_var = [], []
    for s1 in combinations(range(Np), n1):
        for s2 in combinations(s1, n2):
            s1a, s2a = np.array(s1), np.array(s2)
            tp_est.append(Lp[s1a].mean() + (Hp[s2a] - Lp[s2a]).mean())
    tp_est = np.array(tp_est)
    S_H2 = Hp.var(ddof=1)
    S_D2 = (Hp - Lp).var(ddof=1)
    formula = (1 - n1 / Np) * S_H2 / n1 + (1 - n2 / n1) * S_D2 / n2
    out["two_phase_unbiased"] = abs(tp_est.mean() - Hp.mean()) < 1e-12
    out["two_phase_variance_formula_exact"] = abs(tp_est.var() - formula) < 1e-12

    # ---- 4. Blocked Neyman: two blocks, complete randomization within each
    blocks_def = [(6, 3), (5, 2)]                       # (N_b, N_A,b)
    Y1 = [rng.random(nb) < 0.5 for nb, _ in blocks_def]
    Y0 = [(rng.random(nb) < 0.3) for nb, _ in blocks_def]
    Ntot = sum(nb for nb, _ in blocks_def)
    ate = sum(nb / Ntot * (y1.astype(float) - y0.astype(float)).mean()
              for (nb, _), y1, y0 in zip(blocks_def, Y1, Y0))
    assigns = [list(combinations(range(nb), na)) for nb, na in blocks_def]
    dif, vhat = [], []
    for pick in itertools.product(*assigns):
        d, v = 0.0, 0.0
        for (nb, na), A, y1, y0 in zip(blocks_def, pick, Y1, Y0):
            a = np.zeros(nb, bool); a[list(A)] = True
            ya, yb = y1[a].astype(float), y0[~a].astype(float)
            w = nb / Ntot
            d += w * (ya.mean() - yb.mean())
            pa, pb = ya.mean(), yb.mean()
            v += w * w * (pa * (1 - pa) / (a.sum() - 1) + pb * (1 - pb) / ((~a).sum() - 1))
        dif.append(d); vhat.append(v)
    dif, vhat = np.array(dif), np.array(vhat)
    out["blocked_dim_unbiased"] = abs(dif.mean() - ate) < 1e-12
    out["neyman_conservative"] = vhat.mean() >= dif.var() - 1e-12

    # ---- 5. Bias of the local comparison under a linear gradient
    w, slope, tau = 25.0, 0.002, 0.03                  # meV/atom, per meV/atom, effect
    x = np.linspace(-w, w, 200001)                       # uniform running variable
    m = 0.10 + slope * x + tau * (x < 0)                 # side A (stable) is x < 0
    est = m[x < 0].mean() - m[x >= 0].mean()
    out["gradient_bias_equals_minus_slope_w"] = abs((est - tau) - (-slope * w)) < 1e-6

    if verbose:
        for k, v in out.items():
            print(f"[{'PASS' if v else 'FAIL'}] exact: {k}")
    return out


def design_table():
    """Design calculation (not a result): planned 95% half-widths for a pilot cohort
    of N units at an assumed realization rate p, for several LLM qualities."""
    N = 3231   # upper bound: Kirklin et al. 2015 list; the composition-novel subset is smaller
    print("Design calculation (illustrative; replace N and p with pilot values)")
    print(f"N={N}; human-only n=400, p=0.10: +/-{100*srswor_halfwidth(0.10, 400, N):.2f} pp")
    for p in (0.03, 0.10, 0.25):
        for sens, spec in ((0.70, 0.90), (0.85, 0.97), (0.95, 0.99)):
            out = planned_halfwidth_llm_strata(N, p, sens, spec, n_total=400)
            ho = srswor_halfwidth(p, 400, N)
            print(f"p={p:.2f} sens={sens:.2f} spec={spec:.2f}: stratified +/-{100*out['halfwidth']:.2f} pp "
                  f"(human-only +/-{100*ho:.2f} pp; n_h={out['n_h']})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--design", action="store_true")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--floors", action="store_true", help="compare variance-floor variants")
    ap.add_argument("--exact", action="store_true", help="exact enumeration checks of the theory")
    args = ap.parse_args()
    if args.design:
        design_table()
    if args.floors:
        compare_floor_variants()
    if args.exact:
        sys.exit(0 if all(exact_checks().values()) else 1)
    if args.smoke or not (args.smoke or args.design or args.floors or args.exact):
        out = smoke_test(R=args.reps)
        sys.exit(0 if out[SYNTH_PREFIX + "ALL_PASS"] else 1)
