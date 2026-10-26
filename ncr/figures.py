"""Figures 1 to 3 and Supplementary Figures 1 and 2, drawn only from exported results.

Presentation rules enforced in code (each raises if violated):
  1. Every estimated share has a point estimate and a 95% interval (lo <= est <= hi).
  2. Share axes start at zero; there is no smoothing and no interpolation between horizons.
  3. A muted, colorblind-safe palette in which every color keeps one meaning across figures
     (release, record source, or side of the threshold), with marker shapes as a second channel.
  4. Nature's specifications: 183 mm wide; text 5.5-7 pt (panel letters 8 pt bold); strokes
     0.25-1 pt; black text only; colors keyed inside the figure, not described in captions.
After drawing, every point of Fig. 1b is counted back against the exported estimates, and the end of
every curve in Fig. 2b is compared with the exported realized share.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

MIN_PT = 5.5
WIDTH_IN = 183.0 / 25.4                                     # Nature double column
INK = "#1A1A1A"
# Every color has one meaning across Figs 1-3.
COHORT_COLORS = ["#24548C", "#8DB3D8", "#D9822B"]          # OQMD v1.0, OQMD v1.1, Materials Project
STABLE, STABLE_L = "#B8463B", "#D99A92"                     # labeled stable (point, bar)
UNSTABLE, UNSTABLE_L = "#5F5F5F", "#D9D9D9"                 # labeled unstable (point, bar)
GREY, GREY_L, BAND = "#9B9B9B", "#D9D9D9", "#EEF2F7"
MARKERS = ["o", "s", "^", "D", "v"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "TeX Gyre Heros", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold",
    "mathtext.cal": "sans", "mathtext.sf": "sans",
    "font.size": 6.5, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2,
    "ytick.major.size": 2, "xtick.direction": "out", "ytick.direction": "out", "axes.spines.top": False,
    "axes.spines.right": False, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
    "ytick.color": INK, "text.color": INK, "legend.frameon": False, "legend.handlelength": 1.0,
    "legend.handletextpad": 0.4, "legend.borderaxespad": 0.2, "axes.labelpad": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600})


def _letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=8, fontweight="bold", ha="left", va="top")


def _errbar(ax, x, y, lo, hi, horizontal=False, **kw):
    err = [[y - lo], [hi - y]] if not horizontal else [[x - lo], [hi - x]]
    if horizontal:
        ax.errorbar(x, y, xerr=err, fmt="none", ecolor=kw.get("color", INK), elinewidth=0.6, capsize=0,
                    alpha=kw.get("alpha", 1.0))
    else:
        ax.errorbar(x, y, yerr=err, fmt="none", ecolor=kw.get("color", INK), elinewidth=0.6, capsize=0,
                    alpha=kw.get("alpha", 1.0))


def _check_interval(d, where):
    for k in ("est", "lo", "hi"):
        if k not in d:
            raise ValueError(f"{where}: missing '{k}'; every plotted share needs an interval")
    if not (d["lo"] <= d["est"] <= d["hi"]):
        raise ValueError(f"{where}: interval does not contain the estimate")
    if not (0.0 <= d["lo"] and d["hi"] <= 1.0 + 1e-9):
        raise ValueError(f"{where}: share interval outside [0, 1]")



# Outcome colors, used identically in Figs 1-3 and the Supplementary Figures.
KNOWN, LATER, NONE = "#8C8C8C", "#1A1A1A", "#E4E8EE"


def _ticks(axis, lo, hi, nbins=5):
    """Ticks at nice steps within [lo, hi], chosen from the data range, with the same number of
    decimals throughout the axis, zero written as 0 and negative values with a true minus sign."""
    vals = [v for v in matplotlib.ticker.MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
            if lo - 1e-12 <= v <= hi + 1e-12]
    dec = 0
    for v in vals:
        while dec < 4 and abs(round(v, dec) - v) > 1e-9:
            dec += 1
    axis.set_ticks(vals)
    axis.set_ticklabels(["0" if abs(v) < 1e-12 else f"{v:.{dec}f}".replace("-", "\u2212") for v in vals])


def _short(label):
    return label.split(" (")[0]


def _clean(ax, left=True):
    ax.tick_params(axis="y", length=0)
    if not left:
        ax.spines["left"].set_visible(False)


# --------------------------------------------------------------------------------- Figure 1
SRC_SHADES = {"COD only": "#3F3F3F", "both": "#8C8C8C", "ICSD-derived only": "#C6C6C6"}
LINESTYLES = [("-", INK), ("-", "#8C8C8C"), ((0, (3, 1.5)), INK)]


def _hull_panel(ax):
    hx = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    hy = np.array([0.0, -0.40, -0.56, -0.30, 0.0])
    hull = lambda x: np.interp(x, hx, hy)  # noqa: E731
    xs = np.linspace(0.04, 0.96, 300)
    ax.fill_between(xs, hull(xs), hull(xs) + 0.09, color=BAND, lw=0, zorder=0)
    ax.plot(hx, hy, color=INK, lw=0.8, zorder=1)
    ax.plot(hx[[0, -1]], hy[[0, -1]], "o", ms=2.2, color=INK, mew=0, zorder=3)
    ax.plot(hx[1:-1], hy[1:-1], "o", ms=3.8, color=STABLE, mew=0, zorder=3, label="Stable")
    nx = np.array([0.37, 0.62])
    ax.plot(nx, hull(nx) + 0.05, "s", ms=3.0, color=UNSTABLE, mew=0, zorder=3, label="Unstable, near hull")
    fx = np.array([0.2, 0.45, 0.82])
    fy = hull(fx) + np.array([0.26, 0.36, 0.2])
    ax.plot(fx, fy, "o", ms=2.8, mfc="white", mec=GREY, mew=0.6, zorder=3, label="Unstable")
    ax.plot([fx[1]] * 2, [hull(fx[1]), fy[1] - 0.035], color=INK, lw=0.5, ls=(0, (1.5, 1.2)), zorder=2)
    ax.text(fx[1] + 0.03, (hull(fx[1]) + fy[1]) / 2, r"$\Delta H_\mathrm{d}$", fontsize=6.5, va="center")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.7, 0.12)
    ax.set_xticks([0, 1], ["A", "B"])
    ax.set_yticks([0], ["0"])
    ax.set_xlabel("Composition")
    ax.set_ylabel("Formation energy")
    from matplotlib.patches import Patch
    hs, ls = ax.get_legend_handles_labels()
    ax.legend(hs + [Patch(color=BAND)], ls + ["Window"], loc="lower left", bbox_to_anchor=(-0.02, 1.0),
              ncol=4, columnspacing=0.6, handletextpad=0.15, handlelength=0.8)


def _swarm(values, half_height, rng, min_dx):
    """Vertical offsets that spread points sharing nearly the same x within a band."""
    v = np.asarray(values, dtype=float)
    order = np.argsort(v)
    offs = np.zeros(len(v))
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] - v[order[i]] < min_dx:
            j += 1
        k = j - i + 1
        if k == 1:
            offs[order[i]] = rng.uniform(-half_height, half_height)
        else:
            offs[order[i:j + 1]] = np.linspace(-half_height, half_height, k) + rng.uniform(-0.03, 0.03, k)
        i = j + 1
    return offs


def fig1(res, path):
    """Fig. 1. Every prediction of the three releases with a matching experimental structure
    record, placed at the date of its earliest record relative to the release: before it
    (already known; decades scale) or after it, within the longest horizon of complete
    follow-up (recorded later; years scale). Triangles mark the median age of the
    already-known records; counts give every prediction's outcome."""
    d = res["fig1"]
    led = d["ledger"]
    order = sorted(led, key=lambda k: led[k]["date"])
    fig = plt.figure(figsize=(WIDTH_IN, 2.45))
    bottom, height = 0.2, 0.62
    axl = fig.add_axes([0.145, bottom, 0.45, height])                  # before release, decades
    axr = fig.add_axes([0.607, bottom, 0.18, height], sharey=axl)     # after release, years
    rng = np.random.default_rng(11)
    ys = np.arange(len(order))[::-1].astype(float)
    band = 0.34
    fig._ledger_counts = {}
    age_max = 100.0
    for yv, k in zip(ys, order):
        v = led[k]
        e = v["date"]
        ages = e - np.array(v["known_years"])
        axl.scatter(-ages, yv + rng.uniform(-band, band, len(ages)), s=2.4, color=KNOWN,
                    alpha=0.55, lw=0, zorder=2)
        if v.get("age_median") is not None:
            axl.plot([-v["age_median"]], [yv + band + 0.1], "v", ms=4.2, color=INK, mew=0, zorder=4, clip_on=False)
        later = np.array(v["later_offsets"])
        axr.axvspan(0, v["t_max"], ymin=(yv + 0.5 - 0.44) / len(order), ymax=(yv + 0.5 + 0.44) / len(order),
                    color=BAND, lw=0, zorder=0)
        if len(later):
            axr.scatter(later, yv + _swarm(later, 0.28, rng, 0.25), s=13, color=LATER, lw=0, zorder=4)
        n_none = v["n"] - len(ages) - len(later)
        fig._ledger_counts[k] = (len(ages), len(later))
        frac = (yv + 0.5) / len(order)
        lines = [(f"{len(ages):,}", " already known"), (f"{len(later):,}", " recorded later"), (f"{n_none:,}", " no record")]
        for m, (num, lab) in enumerate(lines):
            yy = frac + (0.105 - m * 0.105) * 3 / len(order)
            axr.text(1.24, yy, num, transform=axr.transAxes, fontsize=6, ha="right", va="center",
                     fontweight="bold" if m == 1 else "normal")
            axr.text(1.24, yy, lab, transform=axr.transAxes, fontsize=6, ha="left", va="center")
    axr.axvline(0, color=INK, lw=0.9, zorder=3)                        # the release, once
    axl.set_xlim(-max(age_max, max(float(np.max(led[k]["date"] - np.array(led[k]["known_years"]))) for k in order
                                         if led[k]["known_years"]) + 3), 0)
    axl.set_xticks([-100, -75, -50, -25, 0], ["100", "75", "50", "25", "0"])
    axl.set_xlabel("Years before release")
    axr.set_xlim(0, 10.6)
    axr.set_xticks([0, 2, 4, 6, 8, 10])
    axr.set_xlabel("Years after release")
    axr.spines["left"].set_visible(False)
    axr.tick_params(axis="y", left=False, labelleft=False)
    axl.set_yticks(ys, [f"{_short(k)} ({k.split('(')[-1].rstrip(')')})\n{led[k]['n']:,} predictions" for k in order])
    _clean(axl, left=False)
    axl.set_ylim(-0.6, len(order) - 0.4)
    # axis break marks between the two time scales
    for ax, x0 in ((axl, 1.0), (axr, 0.0)):
        kw = {"transform": ax.transAxes, "color": INK, "lw": 0.6, "clip_on": False}
        ax.plot([x0 - 0.006, x0 + 0.006], [-0.025, 0.025], **kw)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([], [], ls="none", marker="o", ms=2.8, color=KNOWN, label="Earliest record before release (already known)"),
               Line2D([], [], ls="none", marker="v", ms=4.2, color=INK, label="Median age of these records"),
               Line2D([], [], ls="none", marker="o", ms=3.8, color=LATER, label="First record after release (recorded later)"),
               Patch(color=BAND, label="Complete follow-up")]
    axl.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.06), ncol=2, columnspacing=1.2,
               handlelength=1.0, labelspacing=0.3)
    fig.text(0.607, bottom + height + 0.012, "Release", fontsize=6, ha="center", va="bottom")
    fig.savefig(path)
    return fig


# --------------------------------------------------------------------------------- Figure 2
def fig2(res, path):
    """Fig. 2. a, how complete the experimental record is (capture-recapture coverage, records
    dated before 2017); b, predictions already known at release, by the sources that record
    them; c, predictions recorded later at every complete horizon."""
    a, c = res["fig2a"], res.get("fig2c", [])
    cov = res["fig2"]["coverage"]
    rz = res["fig2"]["realized_by_horizon"]
    order = [r["source"] for r in a][::-1]                          # release order, top to bottom
    fig = plt.figure(figsize=(WIDTH_IN, 2.35))
    bottom, height = 0.2, 0.64
    axa = fig.add_axes([0.14, bottom, 0.19, height])
    axb = fig.add_axes([0.385, bottom, 0.2, height])
    axc = fig.add_axes([0.77, bottom, 0.19, height])
    y = np.arange(len(order))[::-1].astype(float)
    # a: coverage of the record
    styles = {"COD": ("s", SRC_SHADES["COD only"], None), "ICSD-derived": ("s", "#A6A6A6", None),
              "either": ("D", INK, "white")}
    offs = {"COD": 0.2, "ICSD-derived": 0.0, "either": -0.2}
    for yy, src in zip(y, order):
        row = next((r for r in cov if r["source"] == src), None)
        if row is None:
            continue
        for part, (mk, col, face) in styles.items():
            d = row[part]
            _check_interval(d, f"fig2a {src} {part}")
            _errbar(axa, 100 * d["est"], yy + offs[part], 100 * d["lo"], 100 * d["hi"], horizontal=True, color=col)
            axa.plot([100 * d["est"]], [yy + offs[part]], mk, ms=3.6 if mk != "D" else 3.4, color=col, mew=0.8,
                     mfc=col if face is None else face, mec=col, zorder=4,
                     label={"COD": "COD", "ICSD-derived": "ICSD-derived", "either": "Either"}[part] if yy == y[0] else None)
    axa.set_yticks(y, order)
    _clean(axa, left=False)
    axa.set_ylim(-0.6, len(order) - 0.4)
    axa.set_xlim(0, 100)
    axa.set_xticks([0, 25, 50, 75, 100])
    axa.set_xlabel("Record coverage (%)")
    axa.legend(loc="lower left", bbox_to_anchor=(-0.04, 1.0), ncol=3, columnspacing=0.8, handletextpad=0.2)
    # b: already known at release, by source
    his = []
    for yy, src in zip(y, order):
        dd = next(r for r in a if r["source"] == src)["already_known"]
        _check_interval(dd, f"fig2b {src}")
        his.append(dd["hi"])
        parts = next(q for q in c if q["source"] == src)
        left = 0.0
        for part in ("COD only", "both", "ICSD-derived only"):
            w = 100 * parts[part] / dd["n"]
            axb.barh(yy, w, left=left, height=0.5, color=SRC_SHADES[part], edgecolor="white", linewidth=0.4,
                     label={"COD only": "COD only", "both": "Both", "ICSD-derived only": "ICSD-derived only"}[part]
                     if yy == y[0] else None)
            left += w
        _errbar(axb, 100 * dd["est"], yy, 100 * dd["lo"], 100 * dd["hi"], horizontal=True)
        axb.text(100 * dd["hi"] + 0.6, yy, f"{100 * dd['est']:.1f}%", fontsize=6, va="center", ha="left")
    axb.set_yticks(y, [""] * len(order))
    _clean(axb, left=False)
    axb.set_ylim(-0.6, len(order) - 0.4)
    axb.set_xlim(0, 100 * max(0.25, 1.2 * max(his)))
    axb.set_xlabel("Already known at release (%)")
    axb.legend(loc="lower left", bbox_to_anchor=(-0.04, 1.0), ncol=3, columnspacing=0.7, handlelength=0.8,
               handletextpad=0.3)
    # c: recorded later at every complete horizon
    rows = []
    for src in order:
        rows += sorted([r for r in rz if r["source"] == src], key=lambda r: r["t"])
    yc, yv, prev = [], 0.0, None
    for r in rows:
        if prev is not None and r["source"] != prev:
            yv -= 0.5
        yc.append(yv)
        yv -= 1.0
        prev = r["source"]
    fig._points = []
    his = []
    for r, yv in zip(rows, yc):
        _check_interval(r, f"fig2c {r['source']} {r['t']}")
        his.append(r["hi"])
        _errbar(axc, 100 * r["est"], yv, 100 * r["lo"], 100 * r["hi"], horizontal=True)
        zero = r["k"] == 0
        axc.plot([100 * r["est"]], [yv], "o", ms=3.6, color=LATER, mfc="white" if zero else LATER, mew=0.8,
                 zorder=4, clip_on=False)
        axc.text(100 * r["hi"] + 0.015, yv, f"{r['k']} of {r['n']:,}", fontsize=6, va="center", ha="left")
        fig._points.append((r["source"], r["t"], 100 * r["est"]))
    axc.set_yticks(yc, [f"{int(r['t'])} yr" for r in rows])
    for src in order:
        ys_ = [yv for r, yv in zip(rows, yc) if r["source"] == src]
        if ys_:
            name, yr = (src.split(" (") + [""])[:2]
            axc.text(-0.2, sum(ys_) / len(ys_), f"{name}\n({yr}" if yr else name, transform=axc.get_yaxis_transform(),
                     fontsize=6, ha="right", va="center")
    _clean(axc, left=False)
    axc.set_ylim(min(yc) - 0.6, 0.6)
    xmax = 100 * 1.4 * max(his)
    axc.set_xlim(0, xmax)
    _ticks(axc.xaxis, 0, xmax)
    axc.set_xlabel("Recorded later (%)")
    _letter(fig, 0.004, 0.975, "a")
    _letter(fig, 0.345, 0.975, "b")
    _letter(fig, 0.6, 0.975, "c")
    fig.savefig(path)
    return fig


# --------------------------------------------------------------------------------- Figure 3
def fig3(res, path):
    """Fig. 3. a, how a release labels a compound (schematic); b, label disagreement between
    databases by distance to the threshold; c, every entry with a structure recorded within five
    years (one marker each) above the number of entries per 5 meV/atom bin; d, the label effect in
    the primary window, the other windows and at placebo thresholds."""
    d = res["fig3"]
    w = d["window"]
    fig = plt.figure(figsize=(WIDTH_IN, 4.3))
    axs = fig.add_axes([0.07, 0.615, 0.33, 0.30])                     # a: what the label is
    axa = fig.add_axes([0.6, 0.615, 0.385, 0.33])                      # b: premise check
    axe = fig.add_axes([0.07, 0.33, 0.4, 0.13])                        # c: events
    axn = fig.add_axes([0.07, 0.085, 0.4, 0.22], sharex=axe)           # c: entries
    axc = fig.add_axes([0.745, 0.085, 0.24, 0.38])                     # d: label effect
    _hull_panel(axs)
    pairs = d.get("premise", [])
    his = []
    styles = [("o", INK), ("s", GREY)]
    for k, pr in enumerate(pairs):
        xs = np.arange(len(pr["bins"])) + (k - (len(pairs) - 1) / 2) * 0.22
        mk, col = styles[k % 2]
        for xv, bb in zip(xs, pr["bins"]):
            if not bb.get("n"):
                continue
            _check_interval({"est": bb["rate"], "lo": bb["ci_low"], "hi": bb["ci_high"]}, f"fig3a {pr['pair']}")
            his.append(bb["ci_high"])
            _errbar(axa, xv, 100 * bb["rate"], 100 * bb["ci_low"], 100 * bb["ci_high"], color=col)
        axa.plot(xs, [100 * bb["rate"] if bb.get("n") else np.nan for bb in pr["bins"]], mk, ms=3.4, color=col,
                 mew=0, ls="none", label=pr["pair"].split(" vs ")[0])
    if pairs:
        labs = [f"{b['bin'][0]:g}\u2013{b['bin'][1]:g}" if b["bin"][1] < 1e6 else f">{b['bin'][0]:g}"
                for b in pairs[0]["bins"]]
        axa.set_xticks(range(len(labs)), labs)
        axa.set_xlim(-0.5, len(labs) - 0.5)
        axa.tick_params(axis="x", labelsize=5.5)
    if pairs:
        inside = [i for i, b in enumerate(pairs[0]["bins"]) if b["bin"][1] <= w]
        if inside:
            axa.axvspan(min(inside) - 0.5, max(inside) + 0.5, color=BAND, lw=0, zorder=0)
    axa.axhline(20, color=INK, lw=0.5, ls=(0, (1.5, 1.5)))
    axa.text(len(pairs[0]["bins"]) - 0.55 if pairs else 5, 20.6, "Criterion", fontsize=6, ha="right", va="bottom")
    axa.set_ylim(0, 100 * 1.15 * max(his) if his else 50)
    axa.set_xlabel(r"Distance from the threshold (meV atom$^{-1}$)")
    axa.set_ylabel("Labels disagree (%)")
    axa.legend(loc="upper right")
    # b
    cent = np.array([bn["center"] for bn in d["bins"]])
    ks = np.array([bn["k"] for bn in d["bins"]])
    ns = np.array([bn["n"] for bn in d["bins"]])
    for ax in (axe, axn):
        ax.axvspan(-w, w, color=BAND, lw=0, zorder=0)
        ax.axvline(0, color=INK, ls=(0, (1.5, 1.5)), lw=0.5, zorder=1)
    for x, k in zip(cent, ks):
        for j in range(int(k)):
            axe.plot([x], [j + 1], "o" if x < 0 else "s", ms=3.4 if x < 0 else 3.0,
                     color=STABLE if x < 0 else UNSTABLE, mew=0, zorder=3)
    kmax = max(2, int(ks.max()))
    axe.set_ylim(0, kmax + 0.9)
    axe.set_yticks(range(0, kmax + 1))
    axe.set_ylabel("Events")
    axe.tick_params(axis="x", labelbottom=False, length=0)
    axe.spines["bottom"].set_visible(False)
    axn.bar(cent, ns, width=4.3, color=np.where(cent < 0, STABLE_L, UNSTABLE_L), edgecolor="none", zorder=2)
    axn.set_ylim(0, 1.12 * ns.max())
    _ticks(axn.yaxis, 0, 1.12 * ns.max(), nbins=4)
    axn.set_ylabel("Entries")
    axn.set_xlim(-2 * w, 2 * w)
    ticks = [-2 * w, -w, 0, w, 2 * w]
    axn.set_xticks(ticks, [f"{v:g}".replace("-", "\u2212") for v in ticks])
    axn.set_xlabel(r"Decomposition energy relative to the threshold (meV atom$^{-1}$)")
    ev = d.get("events") or {}
    if ev.get("N_A"):
        axe.plot([], [], "o", ms=3.4, color=STABLE, mew=0, label=f"Stable ({ev['events_A']} of {ev['N_A']:,} in window)")
        axe.plot([], [], "s", ms=3.0, color=UNSTABLE, mew=0,
                 label=f"Unstable ({ev['events_B']} of {ev['N_B']:,} in window)")
        axe.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, columnspacing=0.8, handletextpad=0.2)
    # c
    rows = [r for r in d.get("forest", []) if r["group"] in {"primary", "window", "placebo"}]
    yy, yv, prev = [], 0.0, None
    for r in rows:
        if prev is not None and r["group"] != prev:
            yv -= 0.45
        yy.append(yv)
        yv -= 0.8
        prev = r["group"]
    mde = d.get("mde")
    if mde:
        for m, xv in enumerate((-100 * mde, 100 * mde)):
            axc.axvline(xv, color=INK, lw=0.6, ls=(0, (1.5, 1.5)), zorder=1,
                        label=f"Smallest reliably detectable effect (\u00b1{100 * mde:.2f})" if m else None)
        axc.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=1.6)
    axc.axvline(0, color=INK, lw=0.5)
    for r, yv in zip(rows, yy):
        col = INK if r["group"] == "primary" else UNSTABLE
        if r.get("note") == "no events":
            axc.plot([0], [yv], "o", ms=3.0, mfc="white", mec=UNSTABLE, mew=0.6)
            continue
        _errbar(axc, 100 * r["est"], yv, 100 * r["lo"], 100 * r["hi"], horizontal=True, color=col)
        axc.plot([100 * r["est"]], [yv], "D" if r["group"] == "primary" else "o",
                 ms=3.6 if r["group"] == "primary" else 2.8, color=col, mew=0)
    U = r"meV atom$^{-1}$"
    lab = {"Primary (\u00b125 meV/atom)": f"\u00b125 {U} (primary)", "Window \u00b110 meV/atom": f"\u00b110 {U}",
           "Window \u00b150 meV/atom": f"\u00b150 {U}", "Placebo threshold at \u221215 meV/atom": f"Placebo at \u221215 {U}",
           "Placebo threshold at +15 meV/atom": f"Placebo at +15 {U}"}
    axc.set_yticks(yy, [lab.get(r["label"], r["label"]) for r in rows])
    _clean(axc, left=False)
    axc.set_ylim(min(yy) - 0.5, 0.55)
    _ticks(axc.xaxis, *axc.get_xlim(), nbins=6)
    axc.set_xlabel("Label effect (percentage points)")
    for x, yv, s_ in ((0.004, 0.99, "a"), (0.5, 0.99, "b"), (0.004, 0.525, "c"), (0.565, 0.525, "d")):
        _letter(fig, x, yv, s_)
    fig.savefig(path)
    return fig


# --------------------------------------------------------------------- Supplementary Figures
def figS1(res, path):
    """Supplementary Fig. 1. The dated releases with their longest complete follow-up, the
    publication years covered by each source of dated structure records, the reporting-lag
    buffer and the data cutoff."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    d = res["fig1"]
    fig = plt.figure(figsize=(WIDTH_IN * 0.7, 2.35))
    ax = fig.add_axes([0.33, 0.17, 0.63, 0.6])
    rel, src = d["releases"], d["sources"]
    ys_s = np.arange(len(src))[::-1].astype(float)
    ys_r = np.arange(len(rel))[::-1] + len(src) + 1.2
    x_lo, cut = 2010, d["cutoff"]
    ax.axvspan(cut - 1, cut, facecolor="#F2F2F2", lw=0, zorder=0, hatch="////", edgecolor="#BFBFBF")
    for r, yv in zip(rel, ys_r):
        ax.barh(yv, r["t_max"], left=r["date"], height=0.42, color="#4D4D4D", edgecolor="none", zorder=2)
        ax.plot([r["date"]], [yv], "D", ms=3.4, color=INK, mew=0, zorder=4)
    for s_, yv in zip(src, ys_s):
        x_end = min(s_["last_year"] + 1, cut)
        ax.barh(yv, x_end - x_lo, left=x_lo, height=0.42, color=GREY_L, edgecolor="none", zorder=2)
        ax.plot([x_lo - 0.18], [yv], "<", ms=4.6, color="#8C8C8C", mew=0, zorder=3, clip_on=False)
    ax.axvline(cut, color=INK, lw=0.6, ls=(0, (1.5, 1.5)), zorder=3)
    ax.set_yticks(list(ys_r) + list(ys_s), [f"{_short(r['label'])} ({r['n']:,} predictions)" for r in rel]
                  + [s_["name"] for s_ in src])
    _clean(ax, left=False)
    ax.tick_params(axis="y", pad=8)
    ax.set_xlim(x_lo, 2027.6)
    ax.set_ylim(-0.6, max(ys_r) + 0.6)
    ax.set_xticks([2010, 2014, 2018, 2022, 2026])
    ax.set_xlabel("Year")
    for yv, txt in ((max(ys_r) + 0.55, "Releases"), (max(ys_s) + 0.55, "Outcome sources")):
        xf, yf = fig.transFigure.inverted().transform(ax.transData.transform((x_lo, yv)))
        fig.text(0.012, yf, txt, fontsize=6, fontweight="bold", ha="left", va="bottom")
    handles = [Line2D([], [], ls="none", marker="D", ms=3.4, color=INK, label="Release"),
               Patch(color="#4D4D4D", label="Longest complete follow-up"),
               Patch(color=GREY_L, label="Dated records (continue before 2010)"),
               Patch(facecolor="#F2F2F2", edgecolor="#BFBFBF", hatch="////", label="Reporting-lag buffer (1 yr)"),
               Line2D([], [], color=INK, lw=0.6, ls=(0, (1.5, 1.5)), label="Data cutoff")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.5, 1.07), ncol=3, columnspacing=0.9,
              handlelength=1.3, labelspacing=0.35)
    fig.savefig(path)
    return fig


def figS2(res, path):
    """Supplementary Fig. 2. Standardized differences across the threshold for each blocking
    scheme, in the 10, 25 (primary) and 50 meV/atom windows, against the pre-specified limit."""
    from matplotlib.patches import Patch
    bal = res["fig3"].get("balance", {})
    covs = bal.get("covariates", [])
    wins = bal.get("windows") or {"25": bal.get("series", [])}
    names = {"n_elements": "Number of elements", "mean_electronegativity": "Mean electronegativity",
             "max_atomic_number": "Largest atomic number"}
    short = {"No blocks": "No blocks", "Chemistry class": "Chemistry class",
             "+ electronegativity": "+ Electronegativity", "+ number of elements (primary)": "+ Number of elements (used)"}
    styles = [("o", "#B3B3B3", "white"), ("s", "#9B9B9B", None), ("^", "#5F5F5F", None), ("o", INK, None)]
    keys = [w for w in ("10", "25", "50") if w in wins]
    fig = plt.figure(figsize=(WIDTH_IN, 2.15))
    left, right, gap = 0.17, 0.985, 0.025
    wpan = (right - left - gap * (len(keys) - 1)) / len(keys)
    vals = [abs(v) for w in keys for ser in wins[w] for v in ser["smd"].values() if v is not None]
    lim = max(0.45, 1.12 * max(vals)) if vals else 0.45
    ys = np.arange(len(covs))[::-1]
    fig._balance = {}
    for m, w in enumerate(keys):
        ax = fig.add_axes([left + m * (wpan + gap), 0.2, wpan, 0.56])
        ax.axvspan(-0.1, 0.1, color=BAND, lw=0, zorder=0)
        ax.axvline(0, color=INK, lw=0.5, zorder=1)
        nser = len(wins[w])
        for k, ser in enumerate(wins[w]):
            mk, col, face = styles[k % len(styles)]
            off = ((nser - 1) / 2 - k) * 0.17
            ax.plot([ser["smd"].get(cv) for cv in covs], ys + off, mk, ms=3.4, mec=col,
                    mfc=col if face is None else face, mew=0.7, ls="none",
                    label=short.get(ser["label"], ser["label"]) if m == 0 else None, zorder=3)
        used = wins[w][-1]["smd"]
        fig._balance[w] = max(abs(v) for v in used.values() if v is not None)
        ax.set_yticks(ys, [names.get(cv, cv) for cv in covs] if m == 0 else [""] * len(covs))
        _clean(ax, left=False)
        if m:
            ax.tick_params(axis="y", length=0)
        ax.set_xlim(-lim, lim)
        ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace("-", "\u2212")))
        ax.set_ylim(-0.6, len(covs) - 0.4)
        ax.set_title(f"\u00b1{w} meV atom$^{{-1}}$" + (" (primary)" if w == "25" else ""), fontsize=6.5, pad=3)
        if m == 1:
            ax.set_xlabel("Standardized difference (stable minus unstable side)")
        if m == 0:
            hs, ls_ = ax.get_legend_handles_labels()
            ax.legend(hs + [Patch(color=BAND)], ls_ + ["Pre-specified limit (|d| < 0.1)"], loc="lower left",
                      bbox_to_anchor=(-0.02, 1.13), ncol=5, columnspacing=0.9, handletextpad=0.15)
    fig.savefig(path)
    return fig


def _min_font_pt(fig):
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()]
    return min(sizes) if sizes else MIN_PT



def _max_font_pt(fig):
    return max(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip())


def _text_is_black(fig):
    from matplotlib.colors import to_hex
    return all(to_hex(t.get_color()) in (INK.lower(), "#000000", "#ffffff")
               for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip())


def _max_stroke_pt(fig):
    from matplotlib.lines import Line2D
    from matplotlib.collections import LineCollection
    widths = [ln.get_linewidth() for ln in fig.findobj(Line2D) if ln.get_linestyle() not in ("None", "none", "")]
    for lc in fig.findobj(LineCollection):
        widths += list(np.atleast_1d(lc.get_linewidth()))
    return max(widths) if widths else 0.0


def _ledger_matches(fig, res):
    """Every point of Fig. 1b is counted back and compared with the exported estimates."""
    want = {r["source"]: (r["already_known"]["k"], r["realized_simultaneous"]["k"]) for r in res["fig2a"]}
    got = getattr(fig, "_ledger_counts", {})
    return bool(got) and all(want.get(k) == v for k, v in got.items())


def _points_match(fig, res):
    """Every point of Fig. 2c equals an exported realized share, and none is missing."""
    got = sorted(getattr(fig, "_points", []))
    want = sorted((r["source"], r["t"], 100 * r["est"]) for r in res["fig2"]["realized_by_horizon"])
    return bool(got) and len(got) == len(want) and all(g[:2] == w[:2] and abs(g[2] - w[2]) < 1e-9 for g, w in zip(got, want))


def _layout_problems(fig):
    """Text that overlaps other text or runs past the figure edge, hyphens used as minus signs,
    and axes whose visible tick labels mix numbers of decimal places."""
    import itertools
    import re
    from matplotlib.text import Text
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    out, boxes = [], []
    for t in fig.findobj(Text):
        if not (t.get_visible() and t.get_text().strip()):
            continue
        bb = t.get_window_extent(r)
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5:
            out.append(f"clipped: {t.get_text()!r}")
        if re.search(r"(^|[\s(])-\d", t.get_text()):
            out.append(f"hyphen as minus: {t.get_text()!r}")
        boxes.append((t, bb))
    for (t1, b1), (t2, b2) in itertools.combinations(boxes, 2):
        if b1.overlaps(b2) and min(b1.x1, b2.x1) - max(b1.x0, b2.x0) > 1 and min(b1.y1, b2.y1) - max(b1.y0, b2.y0) > 1:
            out.append(f"overlap: {t1.get_text()!r} / {t2.get_text()!r}")
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            labs = [lb.get_text().replace("\u2212", "-") for tk, lb in zip(axis.get_majorticklocs(), axis.get_ticklabels())
                    if lb.get_text() and lo - 1e-9 <= tk <= hi + 1e-9]
            nums = [x for x in labs if re.fullmatch(r"-?\d+(\.\d+)?", x) and float(x) != 0]
            if len(nums) >= 2 and len({len(x.split(".")[1]) if "." in x else 0 for x in nums}) > 1:
                out.append(f"mixed decimals: {labs}")
    return out


def _balance_matches(fig, res):
    """Supplementary Fig. 2 draws every exported window, and for the scheme used its largest
    plotted difference equals the exported one."""
    bal = res["fig3"].get("balance", {})
    wins = bal.get("windows") or {"25": bal.get("series", [])}
    got = getattr(fig, "_balance", {})
    want = {w: max(abs(v) for v in ser[-1]["smd"].values() if v is not None) for w, ser in wins.items() if ser}
    return bool(want) and set(got) == set(want) and all(abs(got[w] - want[w]) < 1e-12 for w in want)


def render_all(res: dict, outdir: str) -> dict:
    """Draw every figure into outdir and return the presentation checks."""
    os.makedirs(outdir, exist_ok=True)
    figs = {"fig1": fig1(res, os.path.join(outdir, "fig1.pdf")),
            "fig2": fig2(res, os.path.join(outdir, "fig2.pdf")),
            "fig3": fig3(res, os.path.join(outdir, "fig3.pdf")),
            "figS1": figS1(res, os.path.join(outdir, "figS1.pdf")),
            "figS2": figS2(res, os.path.join(outdir, "figS2.pdf"))}
    f2, f3 = figs["fig2"], figs["fig3"]
    checks = {
        "fig2_axes_start_at_zero": bool(all(ax.get_xlim()[0] == 0 for ax in f2.axes[:3])),
        "fig3_axes_start_at_zero": bool(f3.axes[1].get_ylim()[0] == 0 and f3.axes[2].get_ylim()[0] == 0
                                        and f3.axes[3].get_ylim()[0] == 0),
        "min_font_at_least_5pt": bool(min(_min_font_pt(f) for f in figs.values()) >= MIN_PT),
        "max_font_at_most_8pt": bool(max(_max_font_pt(f) for f in figs.values()) <= 8.0),
        "only_black_or_white_text": bool(all(_text_is_black(f) for f in figs.values())),
        "strokes_at_most_1pt": bool(all(_max_stroke_pt(f) <= 1.0 for f in figs.values())),
        "fig1_points_match_estimates": _ledger_matches(figs["fig1"], res),
        "fig2c_points_match_estimates": _points_match(f2, res),
        "figS2_shows_every_window": _balance_matches(figs["figS2"], res),
        "no_overlap_clipping_hyphen_or_mixed_decimals": not any(_layout_problems(f) for f in figs.values()),
    }
    for f in figs.values():
        plt.close(f)
    return checks
