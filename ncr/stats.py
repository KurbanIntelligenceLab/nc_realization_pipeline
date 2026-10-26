"""Statistics used in the analysis. Every outcome is observed for every unit (a complete
census of the cohort), so shares are exact for the cohort; intervals are Wilson intervals
(a superpopulation view), and the label effect carries only the assignment variance of a
randomized comparison within blocks."""
from __future__ import annotations

import math
from decimal import ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, Decimal
from statistics import NormalDist

import numpy as np


def zcrit(alpha: float = 0.05) -> float:
    """Two-sided standard normal critical value."""
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    return NormalDist().inv_cdf(1.0 - alpha / 2.0)


def wilson_ci(x: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson score interval for x successes in n trials."""
    if n <= 0:
        raise ValueError("n must be positive")
    z = zcrit(alpha)
    p = x / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    lo, hi = max(0.0, mid - half), min(1.0, mid + half)
    if x == 0:
        lo = 0.0          # exact; floating-point rounding otherwise leaves about 1e-18
    if x == n:
        hi = 1.0
    return lo, hi


def share(k: int, n: int, alpha: float = 0.05) -> dict:
    """A share with its Wilson interval and counts."""
    if n == 0:
        return {"est": None, "lo": None, "hi": None, "k": int(k), "n": 0}
    lo, hi = wilson_ci(int(k), int(n), alpha)
    return {"est": k / n, "lo": lo, "hi": hi, "k": int(k), "n": int(n)}


def format_share(est: float, lo: float, hi: float, decimals: int = 1) -> str:
    """Percentage with its interval. The estimate is rounded half-even and the interval
    outward (lower bound down, upper bound up), so rounding never narrows an interval."""
    if not (lo <= est <= hi):
        raise ValueError("interval must contain the estimate")
    q = Decimal(1).scaleb(-decimals)
    to = lambda x, mode: Decimal(repr(float(100.0 * x))).quantize(q, rounding=mode)  # noqa: E731
    return f"{to(est, ROUND_HALF_EVEN)}% [{to(lo, ROUND_FLOOR)}, {to(hi, ROUND_CEILING)}]"


def binomial_two_sided(k: int, n: int) -> float:
    """Exact two-sided binomial test of k successes in n trials at p = 1/2, in log space.
    The p-value sums the probabilities of all outcomes no more likely than the observed one."""
    if n == 0:
        return 1.0
    c = math.lgamma(n + 1) + n * math.log(0.5)
    lp = [c - math.lgamma(i + 1) - math.lgamma(n - i + 1) for i in range(n + 1)]
    lk = lp[k]
    return float(min(1.0, sum(math.exp(v) for v in lp if v <= lk + 1e-9)))


def standardized_difference(a: np.ndarray, b: np.ndarray) -> float:
    va = np.var(a, ddof=1) if a.size > 1 else 0.0
    vb = np.var(b, ddof=1) if b.size > 1 else 0.0
    den = math.sqrt((va + vb) / 2)
    return 0.0 if den == 0 else float((a.mean() - b.mean()) / den)


def blocked_difference(y, in_a, in_b, blocks, alpha: float = 0.05) -> dict:
    """Size-weighted average of within-block differences in means between sides A and B,
    over blocks with at least two units on each side; other blocks are dropped and reported.
    The variance is the plug-in conservative Neyman variance of a difference in means under
    complete randomization within blocks:
        sum_b w_b^2 [p_A,b (1 - p_A,b) / (N_A,b - 1) + p_B,b (1 - p_B,b) / (N_B,b - 1)].
    """
    y = np.asarray(y, dtype=float)
    in_a, in_b = np.asarray(in_a, dtype=bool), np.asarray(in_b, dtype=bool)
    blocks = np.asarray(blocks)
    used, dropped = [], []
    for b in np.unique(blocks[in_a | in_b]):
        na, nb = int((in_a & (blocks == b)).sum()), int((in_b & (blocks == b)).sum())
        (used if (na >= 2 and nb >= 2) else dropped).append(b)
    if not used:
        raise ValueError("no block has at least 2 units on each side")
    sizes = {b: int(((in_a | in_b) & (blocks == b)).sum()) for b in used}
    tot = sum(sizes.values())
    est, var = 0.0, 0.0
    for b in used:
        w = sizes[b] / tot
        A, B = in_a & (blocks == b), in_b & (blocks == b)
        pa, pb = float(y[A].mean()), float(y[B].mean())
        est += w * (pa - pb)
        var += w * w * (pa * (1 - pa) / (A.sum() - 1) + pb * (1 - pb) / (B.sum() - 1))
    se = math.sqrt(var)
    z = zcrit(alpha)
    return {"estimate": est, "se": se, "ci_low": est - z * se, "ci_high": est + z * se,
            "blocks_used": [str(x) for x in used], "blocks_dropped": [str(x) for x in dropped]}


def minimum_detectable_effect(p: float, n_a: int, n_b: int, alpha: float = 0.05, power: float = 0.80) -> float:
    """(z_{1-alpha/2} + z_power) * sqrt(p(1-p)/(n_a-1) + p(1-p)/(n_b-1)); multiplier 2.80
    at alpha = 0.05 and power 0.8."""
    if n_a < 2 or n_b < 2:
        raise ValueError("need at least 2 entries per side")
    mult = NormalDist().inv_cdf(1 - alpha / 2) + NormalDist().inv_cdf(power)
    return mult * math.sqrt(p * (1 - p) / (n_a - 1) + p * (1 - p) / (n_b - 1))


def chapman(n1: int, n2: int, m: int) -> dict:
    """Chapman capture-recapture estimate of the number of units captured by either source,
    with its estimated standard error. Assumes independent capture."""
    total = (n1 + 1) * (n2 + 1) / (m + 1) - 1
    var = (n1 + 1) * (n2 + 1) * (n1 - m) * (n2 - m) / ((m + 1) ** 2 * (m + 2))
    return {"total": total, "se": math.sqrt(var)}


def local_linear_jump(x, y, bandwidth: float, alpha: float = 0.05) -> dict:
    """Sharp local-linear estimate of the jump at x = 0 (side A is x <= 0), with a triangular
    kernel and heteroskedasticity-robust (HC1) standard errors."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    m = np.abs(x) <= bandwidth
    x, y = x[m], y[m]
    d = (x <= 0).astype(float)
    X = np.column_stack([np.ones_like(x), x, d, d * x])
    w = 1 - np.abs(x) / bandwidth
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    res = y - X @ beta
    inv = np.linalg.pinv((X * w[:, None]).T @ X)
    meat = (X * (w * res)[:, None]).T @ (X * (w * res)[:, None])
    n, k = X.shape
    V = inv @ meat @ inv * n / max(n - k, 1)
    se = float(math.sqrt(max(V[2, 2], 0.0)))
    z = zcrit(alpha)
    return {"bandwidth_mev": bandwidth, "n": int(n), "events": int(y.sum()), "estimate": float(beta[2]), "se": se,
            "ci_low": float(beta[2] - z * se), "ci_high": float(beta[2] + z * se)}


def label_disagreement_by_bin(dist_to_hull, stable_a, stable_b, edges, alpha: float = 0.05) -> list[dict]:
    """Share of compounds whose stable labels from two databases disagree, binned by the
    absolute distance to the threshold (meV/atom, from database A), with Wilson intervals."""
    d = np.abs(np.asarray(dist_to_hull, dtype=float))
    a, b = np.asarray(stable_a, dtype=bool), np.asarray(stable_b, dtype=bool)
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


def exact_checks() -> dict:
    """Identities verified by full enumeration on small populations (no Monte Carlo error):
    (1) the blocked difference in means is unbiased over all assignments within blocks;
    (2) its Neyman variance is conservative for the true randomization variance;
    (3) with a linear gradient and uniform spacing, the difference in means equals delta - beta*w."""
    from itertools import combinations
    out = {}
    rng = np.random.default_rng(3)
    blocks = [(4, 2), (5, 3)]                       # (units, treated) per block
    y1 = [rng.integers(0, 2, n).astype(float) for n, _ in blocks]
    y0 = [rng.integers(0, 2, n).astype(float) for n, _ in blocks]
    sizes = np.array([n for n, _ in blocks], dtype=float)
    w = sizes / sizes.sum()
    tau = float(sum(wb * (a - b).mean() for wb, a, b in zip(w, y1, y0)))
    ests, nvars = [], []
    assigns = [list(combinations(range(n), t)) for n, t in blocks]
    for a0 in assigns[0]:
        for a1 in assigns[1]:
            e, v = 0.0, 0.0
            for b, (n, t), A in zip(range(2), blocks, (a0, a1)):
                ta = np.zeros(n, dtype=bool)
                ta[list(A)] = True
                ya, yb = y1[b][ta], y0[b][~ta]
                e += w[b] * (ya.mean() - yb.mean())
                v += w[b] ** 2 * (ya.var(ddof=1) / ya.size + yb.var(ddof=1) / yb.size)
            ests.append(e)
            nvars.append(v)
    ests = np.array(ests)
    out["blocked_difference_unbiased"] = bool(abs(ests.mean() - tau) < 1e-12)
    out["neyman_variance_conservative"] = bool(np.mean(nvars) >= ests.var() - 1e-12)
    wd, beta, delta = 25.0, 0.002, 0.05
    xs = np.linspace(-wd, wd, 2001)
    ys = 0.1 + beta * xs + delta * (xs < 0)
    diff = ys[xs < 0].mean() - ys[xs > 0].mean()
    out["gradient_bias_equals_minus_beta_w"] = bool(abs(diff - (delta - beta * wd)) < 1e-3)
    return out


def simulate_blocked_design(reps: int = 2000, n: int = 2000, effect: float = 0.04, seed: int = 29000) -> dict:
    """Monte Carlo check of the blocked estimator when both the chance of the stable side and
    the base rate differ between two blocks (confounding by block), with a constant effect.
    Reports the mean blocked and unblocked estimates and the coverage of the 95% interval."""
    est_b, est_u, hits = [], [], 0
    for r in range(reps):
        g = np.random.default_rng(seed + r)
        blk = g.integers(0, 2, n)
        a = g.random(n) < np.where(blk == 0, 0.3, 0.7)
        y = (g.random(n) < np.where(blk == 0, 0.02, 0.15) + effect * a).astype(float)
        o = blocked_difference(y, a, ~a, blk)
        est_b.append(o["estimate"])
        est_u.append(float(y[a].mean() - y[~a].mean()))
        hits += o["ci_low"] <= effect <= o["ci_high"]
    return {"reps": reps, "true_effect": effect, "mean_blocked": float(np.mean(est_b)),
            "mean_unblocked": float(np.mean(est_u)), "coverage": hits / reps}
