#!/usr/bin/env python3
"""Figures for "Computational predictions of stable inorganic materials and their
later experimental synthesis".

Draws Fig. 2 (a, b), Fig. 3a and Supplementary Figs. S1 and S2 ONLY from a results
JSON written by nc_realization_metrics.export_results(). Nothing is estimated here.

Presentation rules enforced in code (each raises if violated):
  1. Every plotted share has a point estimate and a 95% interval (lo <= est <= hi).
  2. Share axes start at zero; the upper limit is min(1, 1.2 x largest upper limit).
  3. No smoothing, no interpolation between horizons, no truncated share axes.
  4. Colorblind-safe Okabe-Ito palette, with marker shapes as a second channel.
  5. Text is at least 7 pt at the final 180 mm width.
  6. Synthetic inputs are refused unless --selftest, which watermarks every panel.
After drawing, every plotted coordinate is read back and compared with the input.

Usage:
  python3 nc_realization_figures.py results.json outdir/
  python3 nc_realization_figures.py --selftest
"""
import argparse
import json
import os
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OKABE_ITO = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
             "blue": "#0072B2", "vermillion": "#D55E00", "grey": "#999999"}
MARKERS = ["o", "s", "^", "D", "v"]
MIN_PT = 7.0
WIDTH_MM = 180.0
SYNTH_PREFIX = "SYNTHETIC_"

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "pdf.fonttype": 42})


def _check_interval(d, where):
    for k in ("est", "lo", "hi"):
        if k not in d:
            raise ValueError(f"{where}: missing '{k}'; every plotted share needs an interval")
    if not (d["lo"] <= d["est"] <= d["hi"]):
        raise ValueError(f"{where}: interval does not contain the estimate")
    if not (0.0 <= d["lo"] and d["hi"] <= 1.0 + 1e-9):
        raise ValueError(f"{where}: share interval outside [0, 1]")


def _share_axis_upper(his):
    return min(1.0, 1.2 * max(his)) if his else 1.0


def _watermark(ax, synthetic):
    if synthetic:
        ax.text(0.5, 0.5, "SYNTHETIC TEST\nNOT DATA", transform=ax.transAxes, ha="center",
                va="center", fontsize=14, color="red", alpha=0.35, rotation=20)


def fig2(res, path, synthetic=False):
    a, b = res["fig2a"], res["fig2b"]
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(WIDTH_MM / 25.4, 2.8))
    sources = [row["source"] for row in a]
    y = np.arange(len(sources))
    left = np.zeros(len(sources))
    colors = [OKABE_ITO["grey"], OKABE_ITO["blue"], "#BFE3F7"]
    for part, col, lab in zip(("already_known", "realized", "not_realized"), colors,
                              ("already known", "structure recorded later", "no record")):
        vals = []
        for row in a:
            _check_interval(row[part], f"fig2a {row['source']} {part}")
            vals.append(row[part]["est"])
        vals = np.asarray(vals)
        axa.barh(y, vals, left=left, color=col, edgecolor="black", linewidth=0.4, label=lab)
        lo = np.array([row[part]["lo"] for row in a]); hi = np.array([row[part]["hi"] for row in a])
        axa.errorbar(left + vals, y, xerr=[vals - lo, hi - vals], fmt="none", ecolor="black",
                     elinewidth=0.6, capsize=1.5)
        left = left + vals
    if np.any(np.abs(left - 1.0) > 0.02):
        raise ValueError("fig2a: the three shares must sum to about 1 for each source")
    axa.set_yticks(y, sources); axa.set_xlim(0, 1); axa.set_xlabel("Share of predictions released as new")
    axa.set_title("a  status at the longest complete horizon", loc="left", fontsize=8)
    axa.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, frameon=False)
    his = []
    for k, src in enumerate(sorted({p["source"] for p in b})):
        pts = sorted([p for p in b if p["source"] == src], key=lambda p: p["t"])
        for p in pts:
            _check_interval(p, f"fig2b {src} t={p['t']}")
            his.append(p["hi"])
        t = np.array([p["t"] for p in pts]) + (k - 2) * 0.08
        est = np.array([p["est"] for p in pts])
        err = [est - np.array([p["lo"] for p in pts]), np.array([p["hi"] for p in pts]) - est]
        axb.errorbar(t, est, yerr=err, fmt=MARKERS[k % 5], linestyle="none", ms=4,
                     color=list(OKABE_ITO.values())[k % 6], capsize=2, label=src)
    axb.set_ylim(0, min(1.0, 1.7 * max(his)) if his else 1.0); axb.set_xlim(0, 11)   # headroom for the legend
    axb.set_xticks([0, 2, 3, 5, 10]); axb.set_xlabel("Years since release")
    axb.set_ylabel("Share with a recorded structure")
    axb.set_title("b  predictions not already known", loc="left", fontsize=8)
    axb.legend(loc="upper right", frameon=False, handletextpad=0.3)
    for ax in (axa, axb):
        _watermark(ax, synthetic)
    fig.tight_layout(); fig.savefig(path); return fig


def fig3a(res, path, synthetic=False):
    d = res["fig3a"]
    fig, ax = plt.subplots(figsize=(WIDTH_MM / 25.4 / 2, 2.6))
    w = d["window"]
    ax.axvspan(-w, w, color=OKABE_ITO["sky"], alpha=0.2, lw=0)
    ax.axvline(0, color="black", ls="--", lw=0.8)
    his = []
    for bn in d["bins"]:
        _check_interval(bn, f"fig3a bin {bn['center']}")
        his.append(bn["hi"])
        side_col = OKABE_ITO["vermillion"] if bn["center"] < 0 else OKABE_ITO["blue"]
        ax.errorbar(bn["center"], bn["est"], yerr=[[bn["est"] - bn["lo"]], [bn["hi"] - bn["est"]]],
                    fmt="o" if bn["center"] < 0 else "s", color=side_col, ms=3.5, capsize=1.5)
        ax.annotate(str(bn["n"]), (bn["center"], bn["hi"]), textcoords="offset points",
                    xytext=(0, 3), ha="center", va="bottom", rotation=90, fontsize=7)
    ax.set_ylim(0, 1.15 * _share_axis_upper(his)); ax.set_xlim(-2 * w - 3, 2 * w + 3)
    ax.set_xlabel(r"$\Delta H_\mathrm{d} - \tau_s$ displayed at release (meV per atom)")
    ax.set_ylabel("Share with a structure within 5 years")
    _watermark(ax, synthetic)
    fig.tight_layout(); fig.savefig(path); return fig


def figS1(premise, path):
    """Supplementary Fig. S1: cross-database label disagreement by distance to the threshold
    (check a), from the exported premise results only."""
    fig, ax = plt.subplots(figsize=(WIDTH_MM / 25.4 * 0.62, 2.6))
    pairs = list(premise.items())
    width = 0.38
    for k, (name, rec) in enumerate(pairs):
        bins = rec["bins"]
        xs = np.arange(len(bins)) + (k - 0.5) * width
        est = np.array([b["rate"] for b in bins]); lo = np.array([b["ci_low"] for b in bins]); hi = np.array([b["ci_high"] for b in bins])
        for b in bins:
            _check_interval({"est": b["rate"], "lo": b["ci_low"], "hi": b["ci_high"]}, f"figS1 {name} {b['bin']}")
        ax.bar(xs, est, width, color=list(OKABE_ITO.values())[k], edgecolor="black", linewidth=0.4,
               label=name.replace("oqmd_v1_0", "OQMD v1.0").replace("oqmd_v1_1", "OQMD v1.1").replace(" vs mp_2022_10_28", " vs Materials Project"))
        ax.errorbar(xs, est, yerr=[est - lo, hi - est], fmt="none", ecolor="black", elinewidth=0.6, capsize=1.5)
    labels = []
    for b in pairs[0][1]["bins"]:
        lo_, hi_ = b["bin"]
        labels.append(f"{lo_:g}-{hi_:g}" if hi_ < 1e6 else f">{lo_:g}")
    ax.axhline(0.20, color="black", ls="--", lw=0.8)
    ax.set_xticks(np.arange(len(labels)), labels)
    ax.set_xlabel("Distance to the OQMD threshold (meV per atom)")
    ax.set_ylabel("Share with discordant labels")
    ax.set_ylim(0, 1.2 * max(b["ci_high"] for _, r in pairs for b in r["bins"]))
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(path); return fig


def _min_font_pt(fig):
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()]
    return min(sizes) if sizes else MIN_PT


def _readback_fig2b(fig, res):
    axb = fig.axes[1]
    data = [np.asarray(c.lines[0].get_ydata(), dtype=float) for c in axb.containers]
    plotted = sorted(np.round(np.concatenate(data), 9).tolist())
    expected = sorted(np.round([p["est"] for p in res["fig2b"]], 9).tolist())
    return plotted == expected


def render_all(res, outdir, synthetic=False):
    if not synthetic and any(str(k).startswith(SYNTH_PREFIX) for k in res):
        raise PermissionError("refusing to plot synthetic results")
    os.makedirs(outdir, exist_ok=True)
    f2 = fig2(res, os.path.join(outdir, "fig2.pdf"), synthetic)
    f3 = fig3a(res, os.path.join(outdir, "fig3a.pdf"), synthetic)
    checks = {
        "fig2_axes_start_at_zero": f2.axes[0].get_xlim()[0] == 0 and f2.axes[1].get_ylim()[0] == 0,
        "fig3a_axis_starts_at_zero": f3.axes[0].get_ylim()[0] == 0,
        "min_font_at_least_7pt": min(_min_font_pt(f2), _min_font_pt(f3)) >= MIN_PT,
        "fig2b_points_match_input": _readback_fig2b(f2, res),
    }
    return checks


def _synthetic_results():
    rng = np.random.default_rng(0)
    srcs = ["OQMD 2015", "Materials Project", "AFLOW", "GNoME", "Alexandria"]
    a = []
    for s in srcs:
        k = rng.uniform(0.05, 0.3); r = rng.uniform(0.01, 0.15); n = 1 - k - r
        mk = lambda v: {"est": v, "lo": max(v - 0.02, 0), "hi": min(v + 0.02, 1)}
        a.append({"source": s, "already_known": mk(k), "realized": mk(r), "not_realized": mk(n)})
    b = [{"source": s, "t": t, "est": e, "lo": max(e - 0.01, 0), "hi": e + 0.01}
         for s in srcs[:3] for t, e in zip((3, 5, 10), sorted(rng.uniform(0.01, 0.12, 3)))]
    b += [{"source": s, "t": 2, "est": 0.02, "lo": 0.01, "hi": 0.03} for s in srcs[3:]]
    bins = [{"center": c, "est": 0.05, "lo": 0.03, "hi": 0.07, "n": 100} for c in range(-47, 50, 5)]
    return {"fig2a": a, "fig2b": b, "fig3a": {"window": 25, "bins": bins}}


def selftest():
    res = _synthetic_results()
    with tempfile.TemporaryDirectory() as d:
        checks = render_all(res, d, synthetic=True)
    bad = dict(res); bad["fig2b"] = [dict(res["fig2b"][0], lo=0.5)]
    try:
        fig2(bad, os.path.join(tempfile.gettempdir(), "x.pdf"), True); checks["rejects_bad_interval"] = False
    except ValueError:
        checks["rejects_bad_interval"] = True
    try:
        render_all({SYNTH_PREFIX + "x": 1, **res}, tempfile.gettempdir(), synthetic=False)
        checks["refuses_synthetic"] = False
    except PermissionError:
        checks["refuses_synthetic"] = True
    ok = all(checks.values())
    for k, v in checks.items():
        print(f"[{'PASS' if v else 'FAIL'}] {k}")
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?")
    ap.add_argument("outdir", nargs="?")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest or not args.results:
        sys.exit(0 if selftest() else 1)
    with open(args.results) as fh:
        res = json.load(fh)
    checks = render_all(res, args.outdir)
    print(checks)
    sys.exit(0 if all(checks.values()) else 1)
