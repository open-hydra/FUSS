#!/usr/bin/env python3
"""plot_vv.py -- the figures of the pyrolysis, Darcy and ablating-wall V&V pages (docs/vv/) from a gate run.

    plot_vv.py [--runs DIR] [--out DIR] [--only PAGE[,PAGE...]] [--png DIR]

  --runs  the work root of the gate tests (their FUSS_GATE_WORK; default build/test/runs/gates), which holds one
          Gate_<Name>/ folder per gate
  --out   where the SVG files go (default docs/vv/images)
  --only  the pages to draw, by the name of their docs/vv/ file (default: every page below)
  --png   also write a PNG of every figure into this folder (a preview, drawn from the same figure)

Pages: pyrolysis-two-layer, pyrolysis-atc1, pyrolysis-atc2, pyrolysis-kinetics-0d, pyrolysis-energy-audit,
pyrolysis-literature, darcy-exact, ablating-wall.

The run outputs are read with the readers of each case's own verify.py and of test/pyrolysis/tools/ (imported, not
copied), and every reference is evaluated by the code the gate uses: the exact series and similarity solutions, the
0-D oracle, the committed PATO, FIAT, Amaryllis and digitised literature files. No number is typed into a figure.
Each figure prints the quantities its gate grades, recomputed from the same data, so that a figure and its page can
be compared. A page whose runs are missing is skipped with a SKIP line and the exit code is 1. The SVG files are
deterministic: a fixed hash salt and no date. Python 3.6 with matplotlib 3.3 or later; nothing runs the solver.
"""
import argparse
import contextlib
import glob
import importlib.util
import io
import json
import math
import os
import sys

sys.dont_write_bytecode = True          # importing the cases' scripts writes nothing into the source tree

import matplotlib                       # noqa: E402
matplotlib.use("svg")
import matplotlib.pyplot as plt         # noqa: E402
from matplotlib.lines import Line2D     # noqa: E402

TOOLS = os.path.dirname(os.path.abspath(__file__))
PYRO = os.path.dirname(TOOLS)
REPO = os.path.dirname(os.path.dirname(PYRO))
sys.path.insert(0, TOOLS)
import conv_rule as cr                  # noqa: E402

# ---------------------------------------------------------------------------------------------------------- style
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
BAND = "#b9b8b0"
plt.rcParams.update({
    "svg.hashsalt": "fuss-vv", "font.family": "DejaVu Sans", "font.size": 9, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8, "axes.titlesize": 9.5,
    "axes.titleweight": "normal", "axes.titlelocation": "left", "text.color": INK, "axes.labelcolor": INK,
    "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.axisbelow": True, "xtick.direction": "in", "ytick.direction": "in", "lines.linewidth": 1.5,
    "legend.frameon": True, "legend.framealpha": 0.92, "legend.edgecolor": GRID, "legend.fancybox": False,
    "figure.facecolor": "white", "savefig.facecolor": "white", "axes.facecolor": "white", "path.simplify": True,
})
FW = 7.5                                # figure width, inches


class Missing(Exception):
    pass


def need(*paths):
    for p in paths:
        if not os.path.exists(p):
            raise Missing(p)
    return paths[0] if len(paths) == 1 else paths


def load(name, rel):
    """the module at test/pyrolysis/<rel>, under a name of its own (several scripts are called verify.py)"""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, os.path.join(PYRO, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@contextlib.contextmanager
def quiet():
    """a verify.py helper that prints its own table: keep it out of this tool's output"""
    old = sys.stdout
    sys.stdout = io.StringIO()
    try:
        yield
    finally:
        sys.stdout = old


def gate(runs, g, *parts):
    return os.path.join(runs, g[0], g[1], *parts)


def save(fig, ctx, name):
    path = os.path.join(ctx["out"], name + ".svg")
    fig.savefig(path, format="svg", bbox_inches="tight", pad_inches=0.04, metadata={"Date": None, "Creator": None})
    if ctx["png"]:
        fig.savefig(os.path.join(ctx["png"], name + ".png"), format="png", dpi=150, bbox_inches="tight",
                    pad_inches=0.04)
    plt.close(fig)
    print("  wrote %s (%d bytes)" % (path, os.path.getsize(path)))


def source_legend(ax, entries, loc="best", **kw):
    """the legend of the sources: (label, line style, marker, filled) drawn in ink, independent of the colours"""
    hs = [Line2D([], [], color=INK2, ls=ls, marker=mk, ms=5, mfc=(INK2 if filled else "white"), mec=INK2, lw=1.4)
          for _lab, ls, mk, filled in entries]
    return ax.legend(hs, [e[0] for e in entries], loc=loc, **kw)


def end_labels(ax, items, gap=0.045, x=1.01):
    """direct labels right of the axes at the y of each curve's last value, spread to a minimum gap (axes fraction);
    items: [(y in data units, text, colour)]. Call after the limits are set."""
    to_ax = ax.transAxes.inverted()
    pts = []
    for y, text, col in items:
        ya = to_ax.transform(ax.transData.transform((ax.get_xlim()[1], y)))[1]
        pts.append([ya, text, col])
    pts.sort(key=lambda p: p[0])
    for i in range(1, len(pts)):
        pts[i][0] = max(pts[i][0], pts[i - 1][0] + gap)
    if pts and pts[-1][0] > 1.0:
        pts[-1][0] = 1.0
        for i in range(len(pts) - 2, -1, -1):
            pts[i][0] = min(pts[i][0], pts[i + 1][0] - gap)
    for ya, text, col in pts:
        ax.text(x, ya, text, transform=ax.transAxes, color=INK, fontsize=8, va="center", ha="left", clip_on=False)
        ax.plot([x - 0.008], [ya], transform=ax.transAxes, marker="s", ms=4, color=col, clip_on=False, ls="none")


def ladder(ax, x, e, p, color, marker, label, ref=None, offset=(0, 6), every=True, ms=6, mfc=None, tag=""):
    """one ladder of the convergence rule (tools/conv_rule.py) on log-log axes: the runs' errors e at x (coarsest
    first), the observed order of each pair written beside it, and the gate's window from the second-finest run, the
    errors an order in [BAND_LO p, BAND_HI p] would give at the finest run (shaded). ref: the refinement measure
    (cells, or 1/vnn) when x is not it; every: write the order of every pair (else of the finest pair only)."""
    ref = ref or x
    ax.plot(x, e, color=color, marker=marker, ms=ms, mfc=mfc or color, mec=color, mew=1.2, lw=1.3, label=label,
            zorder=3)
    lo, hi = cr.BAND_LO * p, cr.BAND_HI * p
    n = 24
    xs = [x[-2] * (x[-1] / float(x[-2])) ** (i / float(n - 1)) for i in range(n)]
    rs = [ref[-2] * (ref[-1] / float(ref[-2])) ** (i / float(n - 1)) for i in range(n)]
    ax.fill_between(xs, [e[-2] * (r / ref[-2]) ** (-hi) for r in rs], [e[-2] * (r / ref[-2]) ** (-lo) for r in rs],
                    color=color, alpha=0.16, lw=0, zorder=1)
    orders = []
    for i in range(len(x) - 1):
        o = cr.order(e[i], e[i + 1], ref[i], ref[i + 1])
        orders.append(o)
        if not every and i < len(x) - 2:
            continue
        xm, em = math.sqrt(x[i] * x[i + 1]), math.sqrt(e[i] * e[i + 1])
        ax.annotate(tag + "%.3f" % o, (xm, em), xytext=offset, textcoords="offset points", fontsize=7.5, color=INK2,
                    ha="center" if offset[0] == 0 else ("left" if offset[0] > 0 else "right"),
                    va="bottom" if offset[1] > 0 else ("top" if offset[1] < 0 else "center"))
    return orders


def log_time_axis(ax, t0, t1, label="t (s)"):
    """a log time axis with plain tick labels"""
    from matplotlib.ticker import FixedLocator, NullLocator
    ax.set_xscale("log")
    ax.set_xlim(t0, t1)
    ticks = [v for v in (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 60, 100, 300, 1000) if t0 <= v <= t1
             and not (v == 60 and t1 >= 100)]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels(["%g" % v for v in ticks])
    ax.set_xlabel(label)


def log_ticks(ax, axis="y"):
    """labels at 1, 2 and 5 x 10^k on a log axis that spans a decade or two (else only the powers of ten)"""
    from matplotlib.ticker import LogLocator, FuncFormatter, NullLocator

    def label(v, _pos):
        e = int(math.floor(math.log10(v) + 1e-9))
        m = int(round(v / 10.0 ** e))
        if -2 <= e <= 2:
            return "%g" % (m * 10.0 ** e)
        return r"$10^{%d}$" % e if m == 1 else r"$%d\times10^{%d}$" % (m, e)
    a = ax.yaxis if axis == "y" else ax.xaxis
    lo, hi = ax.get_ylim() if axis == "y" else ax.get_xlim()
    subs = (1.0, 2.0, 5.0) if hi / lo < 200.0 else (1.0,)
    a.set_major_locator(LogLocator(base=10.0, subs=subs))
    a.set_major_formatter(FuncFormatter(label))
    a.set_minor_locator(NullLocator())


def band_handle(label="order window", color=BAND, alpha=0.6):
    from matplotlib.patches import Patch
    return Patch(facecolor=color, alpha=alpha, edgecolor="none", label=label)


def ladder_legend(ax, **kw):
    """the legend of a ladder panel, with the shaded order window"""
    return ax.legend(handles=ax.get_legend_handles_labels()[0] + [band_handle()], **kw)


# ------------------------------------------------------------------------------------------------- two-layer slab
TWO = ("Gate_PyrolysisTwoLayer", "pyrolysis_two_layer")
TWO_MESHES = [  # run name of the 0.9 x limit pairs -> the page's own row names
    ("k20u", "ratio 20, equal cells"), ("k100u", "ratio 100, equal cells"),
    ("k20h5", "ratio 20, HIGH cells 5\u00d7 larger"), ("k20h02", "ratio 20, LOW cells 5\u00d7 larger"),
    ("k100h5", "ratio 100, HIGH cells 5\u00d7 larger"), ("k100h02", "ratio 100, LOW cells 5\u00d7 larger"),
    ("k100d2", r"ratio 100, 2-D, $\Delta y = 2\Delta x$")]
SCHEMES = [("RK3", "o"), ("RK2", "s"), ("euler", "^")]


def two_layer_at_sample(tlv, case):
    """(case.json, exact solution, time, x of the cell centres, T of the first row, error) at the field file nearest
    T_SAMPLE, as verify.py --ladder takes it"""
    cfg = json.load(open(os.path.join(case, "case.json")))
    ex = tlv.two_layer_exact.from_case(cfg)
    hit = sorted(((abs(t - tlv.T_SAMPLE), t, z) for t, _f, z in tlv.fields(case)), key=lambda h: h[:2])
    _d, t, z = hit[0]
    e, _x = tlv.error(case, ex, z, t)
    xc, rows = z[0]
    return cfg, ex, t, xc, rows[0], e


def page_two_layer(ctx):
    tlv = load("vv_two_layer", "two_layer/verify.py")
    base = gate(ctx["runs"], TWO)
    cases = need(*[os.path.join(base, n) for n in ("N020", "N040", "N080")])
    runs = [two_layer_at_sample(tlv, c) for c in cases]
    ns = [r[0]["ncells"] for r in runs]
    es = [r[5] for r in runs]
    print("  (a) max |T - T_exact| at the output nearest %g s: %s K; orders %s" % (
        tlv.T_SAMPLE, ", ".join("N %d %.4e" % (n, e) for n, e in zip(ns, es)),
        ", ".join("%.3f" % cr.order(a, b, m, n) for a, b, m, n in zip(es, es[1:], ns, ns[1:]))))
    cfg = runs[0][0]
    A, L = cfg["A"] * 1e3, cfg["L"] * 1e3
    names = [m["name"] for m in cfg["materials"]]

    # 1: the profile across the jump and its error, per mesh
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.1))
    ex, t = runs[-1][1], runs[-1][2]
    xf = [cfg["L"] * i / 400.0 for i in range(401)]
    a1.plot([v * 1e3 for v in xf], [ex.T(v, t) for v in xf], color=INK, lw=1.2, label="exact", zorder=2)
    marks = ["o", "s", "^"]
    for k, (c, ex_k, t_k, xc, T, e) in enumerate(runs):
        lab = "FUSS, %d cells" % c["ncells"]
        a1.plot([v * 1e3 for v in xc], T, ls="none", marker=marks[k], ms=4.5 - k, mfc="white" if k < 2 else SERIES[k],
                mec=SERIES[k], mew=1.0, color=SERIES[k], label=lab, zorder=3)
        a2.plot([v * 1e3 for v in xc], [Ti - ex_k.T(xi, t_k) for xi, Ti in zip(xc, T)], color=SERIES[k],
                marker=marks[k], ms=3.5, lw=1.0, label=lab)
    for ax in (a1, a2):
        ax.axvline(A, color=INK2, lw=0.8, ls=(0, (4, 3)))
        ax.set_xlim(0, L)
        ax.set_xlabel("x (mm)")
        ax.text(A / 2, 1.0, names[0], transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8,
                color=INK2)
        ax.text((A + L) / 2, 1.0, names[1], transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=8, color=INK2)
    a1.set_ylabel("T (K)")
    a1.legend(loc="upper left")
    a2.set_ylabel(r"$T - T_\mathrm{exact}$ (K)")
    a2.axhline(0.0, color=INK2, lw=0.6)
    a2.legend(loc="lower left")
    fig.tight_layout(w_pad=2.0)
    save(fig, ctx, "pyrolysis-two-layer-profiles")

    # 2: the ladder of (a) and the orders of (c)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.3), gridspec_kw={"width_ratios": [1, 1.25]})
    ladder(a1, ns, es, tlv.P_EXPECTED, SERIES[0], "o", r"max $|T - T_\mathrm{exact}|$ at 8 s", offset=(9, 0))
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xticks(ns)
    a1.set_xticklabels([str(n) for n in ns])
    a1.minorticks_off()
    a1.set_xlabel("cells N")
    log_ticks(a1)
    a1.set_ylabel("error (K)")
    ladder_legend(a1, loc="upper right")
    a1.set_title("(a) equal cells, ratio 20, RK3")
    rows = []
    for i, (key, lab) in enumerate(TWO_MESHES):
        for j, (sch, mk) in enumerate(SCHEMES):
            m1, m2 = (need(os.path.join(base, "sub_%s_%s_m%d" % (key, sch, k))) for k in (1, 2))
            r1, r2 = two_layer_at_sample(tlv, m1), two_layer_at_sample(tlv, m2)
            o = cr.order(r1[5], r2[5], r1[0]["ncells"], r2[0]["ncells"])
            rows.append((i, j, o, r1[5] > r2[5]))
            a2.plot([o], [i + (j - 1) * 0.22], ls="none", marker=mk, ms=5, color=SERIES[j], mfc=SERIES[j],
                    label=sch if i == 0 else None)
    lo, hi = cr.BAND_LO * tlv.P_EXPECTED, cr.BAND_HI * tlv.P_EXPECTED
    a2.axvspan(lo, hi, color=BAND, alpha=0.35, lw=0, zorder=0)
    a2.axvline(tlv.P_EXPECTED, color=INK2, lw=0.8, ls=(0, (4, 3)))
    a2.set_yticks(range(len(TWO_MESHES)))
    a2.set_yticklabels([lab for _k, lab in TWO_MESHES], fontsize=8)
    a2.invert_yaxis()
    a2.set_xlim(lo - 0.1, hi + 0.1)
    a2.set_xlabel("observed order, N to 2N")
    a2.grid(axis="y", lw=0)
    a2.legend(handles=a2.get_legend_handles_labels()[0] + [band_handle()], loc="center right",
              fontsize=7.5)
    a2.set_title("(c) 0.9 \u00d7 the nominal vnn limit")
    fig.tight_layout(w_pad=1.5)
    save(fig, ctx, "pyrolysis-two-layer-order")
    print("  (c) observed orders: min %.3f max %.3f, every error falls: %s" % (
        min(r[2] for r in rows), max(r[2] for r in rows), all(r[3] for r in rows)))


# ------------------------------------------------------------------------------------------- Ablation Test Case 1
PYR = ("Gate_Pyrolysis", "pyrolysis")


def mm(d):
    return "%g mm" % d


def page_atc1(ctx):
    v = load("vv_atc1", "ATC1/verify.py")
    case = need(gate(ctx["runs"], PYR, "ATC1"))
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    px = v.probe_x(pj)
    fiat = v.read_fiat(need(os.path.join(PYRO, "ATC1", "ref", "FIAT_T.dat")))
    series = {}                     # FUSS at the exact depths: the two bracketing probes, linear in x (verify.py)
    for d in v.DEPTHS_MM:
        lo, hi = pj["pairs"][str(d)]
        s_lo = v.read_series(need(os.path.join(case, "OUTPUT", lo + ".txt")))
        s_hi = v.read_series(need(os.path.join(case, "OUTPUT", hi + ".txt")))
        w = ((L - d * 1.0e-3) - px[lo]) / (px[hi] - px[lo])
        series[d] = [(t0, (1.0 - w) * v0 + w * v1) for (t0, v0), (t1, v1) in zip(s_lo, s_hi)]
    series[50] = v.read_series(need(os.path.join(case, "OUTPUT", "back.txt")))
    pato_dir = os.path.join(PYRO, "ATC1", "ref", "pato-298K")
    Ppl = v.read_cols(need(os.path.join(pato_dir, "Ta_plot")))
    Psu = v.read_cols(need(os.path.join(pato_dir, "Ta_surfacePatch")))
    pato = {d: [(r[0], r[c]) for r in rows] for d, (rows, c) in
            {1: (Ppl, 1), 2: (Ppl, 2), 4: (Ppl, 3), 8: (Ppl, 4), 16: (Ppl, 6), 50: (Psu, 2)}.items()}
    depths = v.DEPTHS_MM + [50]

    # the gate's two maxima on FIAT's 1 s grid (verify.py)
    spread = max((abs(v.interp_t(pato[d], row[0]) - row[v.FIAT_COL[d]]), d, row[0])
                 for d in pato for row in fiat if 0.0 <= row[0] <= 60.0)
    dev = {d: [] for d in depths}
    for d in depths:
        for row in fiat:
            t = row[0]
            if t < 1.0 or t > 60.0:
                continue
            f = v.interp_t(series[d], t)
            if f is None and series[d] and 0.0 <= t - series[d][-1][0] <= 0.2:
                f = series[d][-1][1]
            dev[d].append((t, abs(f - v.interp_t(pato[d], t))))
    worst = max((e, d, t) for d in depths for t, e in dev[d])
    env = [(row[0], max(abs(v.interp_t(pato[d], row[0]) - row[v.FIAT_COL[d]]) for d in depths))
           for row in fiat if 0.0 <= row[0] <= 60.0]
    print("  max |T_FUSS - T_PATO| %.2f K at %g mm, %g s; the PATO-FIAT spread %.2f K at %g mm, %g s"
          % (worst[0], worst[1], worst[2], spread[0], spread[1], spread[2]))

    # 1: the temperatures of the three codes at the six depths
    fig, ax = plt.subplots(figsize=(FW, 4.0))
    labels = []
    for i, d in enumerate(depths):
        col = SERIES[i]
        ax.plot([t for t, _ in series[d]], [T for _, T in series[d]], color=col, lw=1.5, zorder=3)
        pm = [(t, T) for t, T in pato[d] if t <= 60.0 + 1e-9 and abs(t / 2.0 - round(t / 2.0)) < 1e-6]
        ax.plot([t for t, _ in pm], [T for _, T in pm], ls="none", marker="o", ms=4, color=col, zorder=4)
        fm = [(r[0], r[v.FIAT_COL[d]]) for r in fiat if int(round(r[0])) % 2 == 1]
        ax.plot([t for t, _ in fm], [T for _, T in fm], ls="none", marker="s", ms=4, mfc="white", mec=col, mew=1.0,
                zorder=4)
        labels.append((series[d][-1][1], mm(d) if d != 50 else "50 mm (back face)", col))
    ax.set_xlim(0, 60)
    ax.set_ylim(250, 1700)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("T (K)")
    source_legend(ax, [("FUSS", "-", None, True), ("PATO 3.1", "none", "o", True), ("FIAT", "none", "s", False)],
                  loc="upper left", ncol=3)
    end_labels(ax, labels, gap=0.05)
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-atc1-temperatures")

    # 2: the deviation from PATO on FIAT's grid, and the tolerance it is held to
    fig, ax = plt.subplots(figsize=(FW, 3.4))
    ax.plot([t for t, _ in env], [e for _, e in env], color=INK2, lw=1.0, ls=(0, (1, 1.5)), zorder=2)
    for i, d in enumerate(depths):
        ax.plot([t for t, _ in dev[d]], [e for _, e in dev[d]], color=SERIES[i], lw=1.3, marker="o", ms=2.5, zorder=3,
                label=mm(d) if d != 50 else "50 mm")
    ax.axhline(spread[0], color=INK, lw=1.0, ls=(0, (5, 3)), zorder=2)
    ax.text(60, spread[0], "tolerance %.2f K  " % spread[0], ha="right", va="bottom", fontsize=8, color=INK)
    ax.plot([worst[2]], [worst[0]], ls="none", marker="o", ms=9, mfc="none", mec=INK, mew=1.0, zorder=5)
    ax.annotate("%.2f K (%g mm, %g s)" % worst, (worst[2], worst[0]), xytext=(-6, 8), textcoords="offset points",
                ha="right", fontsize=8, color=INK)
    ax.set_xlim(0, 61)
    ax.set_ylim(0, spread[0] * 1.15)
    ax.set_xlabel("t (s)")
    ax.set_ylabel(r"$|T_\mathrm{FUSS} - T_\mathrm{PATO}|$ (K)")
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([], [], color=INK2, lw=1.0, ls=(0, (1, 1.5))))
    l.append(r"$|T_\mathrm{PATO} - T_\mathrm{FIAT}|$, max over the depths")
    ax.legend(h, l, loc="center right", ncol=2)
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-atc1-deviation")


# ------------------------------------------------------------------------------- Ablation Test Case 2, chemistry off
DATC2 = ("Gate_DarcyAtc2ChemistryOff", "darcy_atc2_chemistry_off")


def grid_dev(series, pato, depths, interp, col0=1):
    """{depth: [(t, FUSS - PATO)]} on PATO's grid t = 0.1..60 s (the window of both ATC 2 verify.py)"""
    out = {}
    for c, d in enumerate(depths, col0):
        out[d] = []
        for row in pato:
            t = row[0]
            if t < 0.1 - 1e-9 or t > 60.0 + 1e-9:
                continue
            f = interp(series[d], t)
            if f is not None:
                out[d].append((t, f - row[c]))
    return out


def page_atc2(ctx):
    v1 = load("vv_atc2", "ATC2-chemistryOff/verify.py")
    v2 = load("vv_atc2_t2", "darcy/ATC2-chemistryOff-type2/verify.py")
    at = v2.abl_tools
    c1 = need(gate(ctx["runs"], PYR, "ATC2-chemistryOff"))
    c2 = need(gate(ctx["runs"], DATC2, "type2"))
    depths = v2.DEPTHS_MM
    pj1, pj2 = (json.load(open(os.path.join(c, "probes.json"))) for c in (c1, c2))
    s1 = v2.bracket_series(c1, pj1, v2.probe_x(pj1), pj1["L"], 1)
    s2 = v2.bracket_series(c2, pj2, v2.probe_x(pj2), pj2["L"], 1)
    p2 = v2.bracket_series(c2, pj2, v2.probe_x(pj2), pj2["L"], 2)
    P1 = v1.read_pato_plot(need(os.path.join(PYRO, "ATC2-chemistryOff", "ref", "Ta_plot")))
    ref2 = os.path.join(PYRO, "darcy", "ATC2-chemistryOff-type2", "ref")
    P2 = v2.read_pato_plot(need(os.path.join(ref2, "Ta_plot")))
    Pp = v2.read_pato_plot(need(os.path.join(ref2, "p_plot")))
    sp21, pk21 = at.spread_atc21()
    tol = sp21["T"][0]
    dp_tol = v2.DP_T_TERMS_PER_K * tol + sp21["mdot_g"][0] / pk21["mdot_g"]
    d1 = grid_dev(s1, P1, depths, v2.interp_t)
    d2 = grid_dev(s2, P2, depths, v2.interp_t)
    w1 = max((abs(e), d, t) for d in depths for t, e in d1[d])
    w2 = max((abs(e), d, t) for d in depths for t, e in d2[d])
    dp = {d: [(t, f - 0.0) for t, f in p2[d]] for d in depths}
    ddp = grid_dev(p2, Pp, depths, v2.interp_t)
    pmax = max(r[c] - v2.P_WALL for r in Pp if 0.1 - 1e-9 <= r[0] <= 60.0 + 1e-9 for c in range(1, len(depths) + 1))
    wp = max((abs(e), d, t) for d in depths for t, e in ddp[d])
    print("  (T) type 1: max |T_FUSS - T_PATO| %.2f K at %g mm, %.1f s; type 2: %.2f K at %g mm, %.1f s; "
          "tolerance %.2f K" % (w1 + w2 + (tol,)))
    print("  (p) max |dp_FUSS - dp_PATO| %.1f Pa at %g mm, %.1f s; PATO peak dp %.1f Pa: ratio %.4f, tolerance %.4f"
          % (wp[0], wp[1], wp[2], pmax, wp[0] / pmax, dp_tol))

    # 1: the temperatures of the type-2 run against PATO (type 1 differs by a few K: figure 2)
    fig, ax = plt.subplots(figsize=(FW, 4.0))
    labels = []
    for i, d in enumerate(depths):
        col = SERIES[i]
        ax.plot([t for t, _ in s2[d]], [T for _, T in s2[d]], color=col, lw=1.5, zorder=3)
        pm = [r for r in P2 if r[0] <= 60.0 + 1e-9 and abs(r[0] / 2.0 - round(r[0] / 2.0)) < 1e-6]
        ax.plot([r[0] for r in pm], [r[i + 1] for r in pm], ls="none", marker="o", ms=4, color=col, zorder=4)
        labels.append((s2[d][-1][1], mm(d), col))
    ax.set_xlim(0, 60)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("T (K)")
    ax.set_ylim(250, None)
    source_legend(ax, [("FUSS, type 2", "-", None, True), ("PATO 3.1", "none", "o", True)], loc="upper left",
                  ncol=2)
    end_labels(ax, labels, gap=0.05)
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-atc2-temperatures")

    # 2: the deviation from PATO of each closure against the borrowed tolerance (log time: the peaks are near 1 s)
    fig, axs = plt.subplots(1, 2, figsize=(FW, 3.3), sharey=True)
    for ax, dev, w, title in ((axs[0], d1, w1, "type 1"), (axs[1], d2, w2, "type 2 (Darcy)")):
        for i, d in enumerate(depths):
            ax.plot([t for t, _ in dev[d]], [abs(e) for _, e in dev[d]], color=SERIES[i], lw=1.2, label=mm(d))
        ax.axhline(tol, color=INK, lw=1.0, ls=(0, (5, 3)))
        ax.text(0.1, tol, " tolerance %.2f K" % tol, ha="left", va="bottom", fontsize=8, color=INK)
        ax.plot([w[2]], [w[0]], ls="none", marker="o", ms=9, mfc="none", mec=INK, mew=1.0, zorder=5)
        ax.annotate("%.2f K (%g mm, %g s)" % w, (w[2], w[0]), xytext=(8, 2), textcoords="offset points", fontsize=8,
                    color=INK, va="bottom")
        log_time_axis(ax, 0.1, 60)
        ax.set_title(title)
    axs[0].set_ylim(0, tol * 1.15)
    axs[0].set_ylabel(r"$|T_\mathrm{FUSS} - T_\mathrm{PATO}|$ (K)")
    axs[1].legend(loc="center right", ncol=1, fontsize=7.5)
    fig.tight_layout(w_pad=1.0)
    save(fig, ctx, "pyrolysis-atc2-deviation")

    # 3: the gas pressure rise of the type-2 run against PATO's
    fig, ax = plt.subplots(figsize=(FW, 3.6))
    labels = []
    for i, d in enumerate(depths):
        col = SERIES[i]
        ax.plot([t for t, _ in dp[d]], [f - v2.P_WALL for _, f in dp[d]], color=col, lw=1.5, zorder=3)
        pm = [r for r in Pp if r[0] <= 60.0 + 1e-9 and abs(r[0] / 2.0 - round(r[0] / 2.0)) < 1e-6]
        ax.plot([r[0] for r in pm], [r[i + 1] - v2.P_WALL for r in pm], ls="none", marker="o", ms=4, color=col,
                zorder=4)
        labels.append((dp[d][-1][1] - v2.P_WALL, mm(d), col))
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 1.12 * max(pmax, max(f - v2.P_WALL for d in depths for _t, f in dp[d])))
    ax.set_xlabel("t (s)")
    ax.set_ylabel(r"$p - p_w$ (Pa)")
    source_legend(ax, [("FUSS, type 2", "-", None, True), ("PATO 3.1", "none", "o", True)], loc="upper left")
    end_labels(ax, labels, gap=0.055)
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-atc2-pressure")


# ------------------------------------------------------------------------------------ kinetics and 0-D decomposition
D0D = ("Gate_Darcy0d", "darcy_0d")
LADDERS_0D = [  # label, gate, runs (coarsest first), the verify.py whose oracle step it takes
    ("adiabatic0d, type 1", PYR, ["adiabatic0d_T1000", "adiabatic0d_T1000_vnn0.125", "adiabatic0d_T1000_vnn0.0625"],
     "adiabatic0d/verify.py"),
    ("heated0d, type 1", PYR, ["heated0d", "heated0d_vnn0.0015625", "heated0d_vnn0.00078125"], "heated0d/verify.py"),
    ("adiabatic0d_type2", D0D, ["adiabatic0d_type2_vnn0.0625", "adiabatic0d_type2_vnn0.03125",
                                "adiabatic0d_type2_vnn0.015625"], "adiabatic0d/verify.py"),
    ("heated0d_type2", D0D, ["heated0d_type2", "heated0d_type2_vnn0.0015625", "heated0d_type2_vnn0.00078125"],
     "heated0d/verify.py")]


def ladder_0d(ctx, label, g, names, script):
    """the errors of one 0-D ladder against its oracle, as ode0d.check_ladder computes them: the oracle of the
    run's closure integrated once on the union of the runs' field-file times"""
    od = load("vv_ode0d", "tools/ode0d.py")
    dt = load("vv_" + script.split("/")[0], script).DT_ORACLE
    cases = [need(gate(ctx["runs"], g, n)) for n in names]
    variant = "type2" if od.input_value(cases[0], "gas-closure", "type1") == "type2" else "type1"
    objs = [od.Case(c, pato_tables=None) for c in cases]
    recs = [o.fuss_fields() for o in objs]
    times = sorted(set(r["t"] for rr in recs for r in rr))
    T, chi = objs[0].integrate(variant, times[-1], dt, times)
    orc = dict(zip(times, zip(T, chi)))
    vnn = [float(od.input_value(c, "vnn")) for c in cases]
    eT = [max(abs(r["T"] - orc[r["t"]][0]) for r in rr) for rr in recs]
    ec = [max(abs(r["chi"][k] - orc[r["t"]][1][k]) for r in rr for k in range(objs[0].nsub)) for rr in recs]
    print("  %s (%s oracle): vnn %s; max |T - T_or| %s K; max |chi - chi_or| %s" % (
        label, variant, " ".join("%g" % v for v in vnn), " ".join("%.4e" % e for e in eT),
        " ".join("%.3e" % e for e in ec)))
    return {"label": label, "vnn": vnn, "eT": eT, "ec": ec, "recs": recs, "orc": orc, "times": times,
            "nsub": objs[0].nsub}


def page_kinetics(ctx):
    kv = load("vv_kinetics", "kinetics_isothermal/verify.py")

    # 1: isothermal advancement against the closed-form propagator
    fig, ax = plt.subplots(figsize=(FW, 3.6))
    styles = ["-", (0, (5, 2.5)), (0, (1, 1.6))]
    marks = ["o", "s", "^"]
    worst = 0.0
    for i, name in enumerate(("T0700", "T0900", "T1200")):
        case = need(gate(ctx["runs"], PYR, name))
        _rv, reac = kv.read_pyrolysis_dat(need(os.path.join(case, "INPUT", "pyrolysis.dat")))
        T0 = kv.read_tec(os.path.join(case, "INPUT", "ic.tec"))[2]["T"][0]
        files = sorted(glob.glob(os.path.join(case, "OUTPUT", "field*.tec")))
        snaps = sorted((kv.read_tec(f)[1], kv.read_tec(f)[2]) for f in files)
        tt = [10 ** (-2 + 3.6 * k / 300.0) for k in range(301)]
        for r, rc in enumerate(reac):
            ax.plot(tt, [kv.chi_exact(rc, T0, t) for t in tt], color=SERIES[i], ls=styles[r], lw=1.3, zorder=2)
            ax.plot([t for t, _d in snaps], [d["chi_%d" % (r + 1)][0] for _t, d in snaps], ls="none", marker=marks[r],
                    ms=4.5, mfc="white" if r else SERIES[i], mec=SERIES[i], mew=1.0, zorder=3)
            worst = max([worst] + [abs(c - kv.chi_exact(rc, T0, t)) for t, d in snaps for c in d["chi_%d" % (r + 1)]])
        ax.text(1.01, kv.chi_exact(reac[0], T0, tt[-1]), "%g K" % T0, transform=ax.get_yaxis_transform(), fontsize=8,
                va="center", color=INK)
    print("  isothermal: max |chi - chi_exact| %.3e over every cell and field file" % worst)
    log_time_axis(ax, 0.01, 40)
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel(r"$\chi_r$")
    hs = [Line2D([], [], color=SERIES[i], lw=2, label="%s K" % T) for i, T in enumerate((700, 900, 1200))]
    hs += [Line2D([], [], color=INK2, ls=styles[r], marker=marks[r], ms=4.5, mfc=INK2 if r == 0 else "white",
                  mec=INK2, label="reaction %d" % (r + 1)) for r in range(3)]
    ax.legend(handles=hs, loc="upper left", ncol=2, title="closed form (line), FUSS (markers)", title_fontsize=8)
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-kinetics-0d-isothermal")

    lad = [ladder_0d(ctx, *spec) for spec in LADDERS_0D]

    # 2: the type-1 trajectories against the oracle, at each case's own vnn
    fig, axs = plt.subplots(2, 2, figsize=(FW, 5.0), sharex="col")
    for col, L in enumerate(lad[:2]):
        rr, orc = L["recs"][0], L["orc"]
        ts = L["times"]
        tlast = rr[-1]["t"]
        ot = [t for t in ts if t <= tlast + 1e-12]
        axs[0, col].plot(ot, [orc[t][0] for t in ot], color=INK, lw=1.2, zorder=2)
        axs[0, col].plot([r["t"] for r in rr], [r["T"] for r in rr], ls="none", marker="o", ms=4, color=SERIES[0],
                         zorder=3)
        for k in range(L["nsub"]):
            axs[1, col].plot(ot, [orc[t][1][k] for t in ot], color=SERIES[k], lw=1.2, ls=styles[k], zorder=2)
            axs[1, col].plot([r["t"] for r in rr], [r["chi"][k] for r in rr], ls="none", marker=marks[k], ms=4,
                             mfc="white" if k else SERIES[k], mec=SERIES[k], zorder=3, label="reaction %d" % (k + 1))
        axs[0, col].set_title("%s (type 1), vnn %g" % (L["label"].split(",")[0], L["vnn"][0]))
        axs[1, col].set_xlabel("t (s)")
        axs[1, col].set_ylim(-0.03, 1.03)
    axs[0, 0].set_ylabel("T (K)")
    axs[1, 0].set_ylabel(r"$\chi_r$")
    source_legend(axs[0, 1], [("oracle (tools/ode0d.py)", "-", None, True), ("FUSS", "none", "o", True)],
                  loc="lower right")
    axs[1, 1].legend(loc="lower right")
    fig.tight_layout(h_pad=0.8, w_pad=1.5)
    save(fig, ctx, "pyrolysis-kinetics-0d-trajectories")

    # 3: the ladders of the four cases
    fig, axs = plt.subplots(1, 2, figsize=(FW, 3.4))
    for i, L in enumerate(lad):
        ref = [1.0 / v for v in L["vnn"]]
        for ax, err in ((axs[0], L["eT"]), (axs[1], L["ec"])):
            ladder(ax, L["vnn"], err, 2, SERIES[i], marks[i % 3] if i < 3 else "D", L["label"], ref=ref,
                   offset=(9, 0), every=False)
    for ax, yl in ((axs[0], r"max $|T - T_\mathrm{or}|$ (K)"), (axs[1], r"max $|\chi - \chi_\mathrm{or}|$")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        log_ticks(ax)
        ax.set_xlabel("vnn")
        ax.set_ylabel(yl)
    ladder_legend(axs[1], loc="lower right", fontsize=7.5)
    fig.tight_layout(w_pad=1.5)
    save(fig, ctx, "pyrolysis-kinetics-0d-order")


# ------------------------------------------------------------------------------- energy audit and pure conduction
def page_energy_audit(ctx):
    """the order ladder of the pure-conduction case"""
    cv = load("vv_conduction", "conduction/verify.py")
    cases = [need(gate(ctx["runs"], PYR, n)) for n in ("conduction_N50", "conduction_graded100", "conduction_N200")]
    with quiet():
        res = [cv.errors(c) for c in cases]
    ns = [r[0] for r in res]
    fig, ax = plt.subplots(figsize=(FW * 0.62, 3.3))
    labs = ["(a) max over every depth and time", "(b) 1 mm, t = 1 s", r"(c) max at $\geq$ 2 mm"]
    for k in (1, 2, 3):
        e = [r[k] for r in res]
        o = ladder(ax, ns, e, cv.P_EXPECTED, SERIES[k - 1], ["o", "s", "^"][k - 1], labs[k - 1],
                   offset=[(9, 0), (9, -7), (9, 7)][k - 1], every=False, ms=[6, 7, 5][k - 1],
                   mfc=[None, "white", None][k - 1], tag=labs[k - 1][:4])
        print("  conduction %s: %s K, order %.3f" % (labs[k - 1].split(")")[0] + ")", " ".join("%.4f" % v for v in e),
                                                       o[-1]))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.set_xticklabels([str(n) for n in ns])
    ax.minorticks_off()
    ax.set_xlim(ns[0] / 1.25, ns[-1] * 1.6)
    ax.set_xlabel("cells N")
    log_ticks(ax)
    ax.set_ylabel(r"$|T - T_\mathrm{exact}|$ (K)")
    ladder_legend(ax, loc="lower left", fontsize=7.5)
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-energy-audit-conduction-order")


# ---------------------------------------------------------------------------------------------- literature benchmarks
BAUM = ("Gate_BenchBaumAtreya", "bench_baum_atreya")
ZHANG = ("Gate_BenchZhang2010", "bench_zhang2010")


def page_literature(ctx):
    bv = load("vv_baum", "lit/baum_atreya2007/verify.py")
    cases = [need(gate(ctx["runs"], BAUM, n)) for n in ("full_N100", "full_N200")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.4))
    marks = ["o", "s"]
    for k, c in enumerate(cases):
        m = json.load(open(os.path.join(c, "case.json")))
        exf, eL2, nsnap = bv.front_and_T(c, m)
        print("  Baum & Atreya %d cells: front error %.4f, L2 T error %.3f K (averages over %d snapshots, t >= %g s)"
              % (m["ncells"], exf, eL2, nsnap, bv.TMIN))
        snaps = []
        for f in glob.glob(os.path.join(c, "OUTPUT", "field*.tec")):
            dims, t, d = bv.dt.read_tec(f)
            ni = dims[0] - 1
            x = d["x"]
            xc = [0.5 * (x[i] + x[i + 1]) for i in range(ni)]
            h = xc[1] - xc[0]
            xf = sum(m["rv"] - r for r in d["rho_s"][:ni]) * h / (m["rv"] - m["rc"])    # the gate's mass-based front
            snaps.append((t, xf, xc, d["T"][:ni]))
        snaps.sort(key=lambda r: r[0])
        snaps = [r for r in snaps if r[0] > 0.0]
        a1.plot([r[0] for r in snaps], [r[1] * 1e3 for r in snaps], ls="none", marker=marks[k], ms=4.5,
                mfc="white" if k == 0 else SERIES[k], mec=SERIES[k], color=SERIES[k],
                label="FUSS, %d cells" % m["ncells"], zorder=3)
        t, _xf, xc, T = snaps[-1]
        a2.plot([v * 1e3 for v in xc], T, ls="none", marker=marks[k], ms=4, mfc="white" if k == 0 else SERIES[k],
                mec=SERIES[k], label="FUSS, %d cells" % m["ncells"], zorder=3)
    tend = snaps[-1][0]
    tt = [tend * i / 300.0 for i in range(301)]
    a1.plot(tt, [m["Cstar"] * math.sqrt(m["alpha_c"] * t) * 1e3 for t in tt], color=INK, lw=1.3, zorder=2,
            label=r"exact, $C^\star\sqrt{\alpha_c t}$")
    a1.axvspan(0, bv.TMIN, color=BAND, alpha=0.3, lw=0, zorder=0)
    a1.text(bv.TMIN / 2, 0.97, "not\naveraged", transform=a1.get_xaxis_transform(), ha="center", va="top",
            fontsize=7.5, color=INK2)
    a1.set_xlim(0, tend)
    a1.set_ylim(0, None)
    a1.set_xlabel("t (s)")
    a1.set_ylabel("front depth (mm)")
    a1.legend(loc="lower right")
    X = m["Cstar"] * math.sqrt(m["alpha_c"] * t)
    xs = [3.5 * X * i / 400.0 for i in range(401)]
    a2.plot([v * 1e3 for v in xs], [bv.dx.sharp_front_T(v, t, m["lam"], m["alpha_c"], m["alpha_v"], m["rv"], m["rc"],
                                                         m["cc"], m["cg"], m["Ts"], m["Tp"], m["T0"]) for v in xs],
            color=INK, lw=1.3, label="exact", zorder=2)
    a2.axvline(X * 1e3, color=INK2, lw=0.8, ls=(0, (4, 3)))
    a2.text(X * 1e3, 0.97, " front", transform=a2.get_xaxis_transform(), fontsize=7.5, color=INK2, va="top")
    a2.axvspan(3 * X * 1e3, 3.5 * X * 1e3, color=BAND, alpha=0.3, lw=0, zorder=0)
    a2.text(3.25 * X * 1e3, 0.97, r"$x > 3X$", transform=a2.get_xaxis_transform(), fontsize=7.5, color=INK2,
            ha="center", va="top")
    a2.set_xlim(0, 3.5 * X * 1e3)
    a2.set_xlabel("x from the heated surface (mm)")
    a2.set_ylabel("T (K)")
    a2.set_title("t = %.1f s" % t)
    a2.legend(loc="upper right", bbox_to_anchor=(0.84, 1.0))
    fig.tight_layout(w_pad=1.5)
    save(fig, ctx, "pyrolysis-literature-baum-atreya")

    zv = load("vv_zhang", "lit/zhang2010/verify.py")
    case = need(gate(ctx["runs"], ZHANG, "N050"))
    rT, rM = zv.read_ref("zhang_fig2_2_T.csv"), zv.read_ref("zhang_fig2_3_mloss.csv")
    sF = zv.fuss_series(case)
    planes = [("exposed", "exposed face"), ("middle", "mid-plane"), ("unexposed", "unexposed face")]
    codes = [("MATLAB", "o", True, 0), ("UEL", "s", False, 3), ("UMATHT", "^", False, 6)]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.5), gridspec_kw={"width_ratios": [1.15, 1]})
    for i, (key, lab) in enumerate(planes):
        col = SERIES[i]
        a1.plot([t for t, _ in sF[key]], [v for _, v in sF[key]], color=col, lw=1.5, zorder=3)
        for code, mk, filled, off in codes:
            ser = rT["%s_%s_K" % (code, key)]
            pts = [r for j, r in enumerate(ser) if j % 10 == off]
            a1.plot([t for t, _ in pts], [v for _, v in pts], ls="none", marker=mk, ms=4.5, color=col,
                    mfc=col if filled else "white", mec=col, zorder=4)
        mref = rT["MATLAB_%s_K" % key]
        dev = [(t, abs(v - zv.interp(mref, t))) for t, v in sF[key] if t > 20.0 and zv.interp(mref, t) is not None]
        other = [rT["%s_%s_K" % (c, key)] for c in ("UEL", "UMATHT") if rT.get("%s_%s_K" % (c, key))]
        spread = max(abs(v - zv.interp(o, t)) for o in other for t, v in mref
                     if t > 20.0 and zv.interp(o, t) is not None)
        print("  Zhang %-9s max |T - T_MATLAB| %.2f K (t > 20 s); the codes' spread %.2f K"
              % (key, max(v for _t, v in dev), spread))
        a2.plot([t for t, _ in dev], [v for _, v in dev], color=col, lw=1.2,
                label="%s, tolerance %.2f K" % (lab, spread))
        a2.axhline(spread, color=col, lw=1.0, ls=(0, (5, 3)))
    a1.set_xlim(0, 2000)
    a1.set_xlabel("t (s)")
    a1.set_ylabel("T (K)")
    hs = [Line2D([], [], color=SERIES[i], lw=2) for i in range(3)]
    l1 = a1.legend(hs, [lab for _k, lab in planes], loc="lower right")
    a1.add_artist(l1)
    source_legend(a1, [("FUSS", "-", None, True)] + [(c, "none", mk, f) for c, mk, f, _o in codes], loc="center right",
                  bbox_to_anchor=(1.0, 0.45))
    a2.set_xlim(0, 2000)
    a2.set_ylim(0, 7.0)
    a2.set_xlabel("t (s)")
    a2.set_ylabel(r"$|T_\mathrm{FUSS} - T_\mathrm{MATLAB}|$ (K)")
    a2.legend(loc="upper right", fontsize=7.5, bbox_to_anchor=(1.0, 0.7))
    fig.tight_layout(w_pad=1.5)
    save(fig, ctx, "pyrolysis-literature-zhang-temperature")

    exp = json.load(open(os.path.join(PYRO, "lit", "zhang2010", "expected.json")))
    pk = max(sF["mdot"], key=lambda r: r[1])
    refpk = [max(rM["%s_kg_s_m2" % c], key=lambda r: r[1]) for c, _m, _f, _o in codes]
    lo, hi = exp["mdot_lo"] * min(r[1] for r in refpk), exp["mdot_hi"] * max(r[1] for r in refpk)
    tmid = sum(r[0] for r in refpk) / len(refpk)
    print("  Zhang peak gas flux: FUSS %.4e at %.0f s; codes %s; window [%.4e, %.4e] x [%.0f, %.0f] s"
          % (pk[1], pk[0], ", ".join("%.4e at %.0f s" % r[::-1] for r in refpk), lo, hi, tmid - exp["dt_peak"],
             tmid + exp["dt_peak"]))
    fig, ax = plt.subplots(figsize=(FW, 3.3))
    ax.plot([t for t, _ in sF["mdot"]], [v * 1e3 for _, v in sF["mdot"]], color=SERIES[0], lw=1.6, zorder=3)
    for code, mk, filled, off in codes:
        ser = rM["%s_kg_s_m2" % code]
        pts = [r for j, r in enumerate(ser) if j % 4 == off % 4]
        ax.plot([t for t, _ in pts], [v * 1e3 for _, v in pts], ls="none", marker=mk, ms=4.5, color=INK2,
                mfc=INK2 if filled else "white", mec=INK2, zorder=4)
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((tmid - exp["dt_peak"], lo * 1e3), 2 * exp["dt_peak"], (hi - lo) * 1e3, facecolor=BAND,
                           alpha=0.45, edgecolor=INK2, lw=0.8, zorder=2))
    ax.plot([pk[0]], [pk[1] * 1e3], ls="none", marker="o", ms=9, mfc="none", mec=INK, mew=1.0, zorder=5)
    ax.annotate("FUSS peak %.3f g/(m$^2$s) at %.0f s" % (pk[1] * 1e3, pk[0]), (pk[0], pk[1] * 1e3), xytext=(70, 22),
                textcoords="offset points", fontsize=8, color=INK, va="bottom",
                arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.8})
    ax.set_xlim(0, 2000)
    ax.set_ylim(0, 1.25 * hi * 1e3)
    ax.set_xlabel("t (s)")
    ax.set_ylabel(r"$\dot m_{g}$ at the exposed face (g/(m$^2$s))")
    hs = [Line2D([], [], color=SERIES[0], lw=1.6)]
    hs += [Line2D([], [], ls="none", marker=mk, ms=4.5, color=INK2, mfc=INK2 if f else "white", mec=INK2)
           for _c, mk, f, _o in codes]
    from matplotlib.patches import Patch
    hs.append(Patch(facecolor=BAND, alpha=0.45, edgecolor=INK2, lw=0.8))
    ax.legend(hs, ["FUSS"] + [c for c, _m, _f, _o in codes] + ["peak window of the check"], loc="upper left")
    fig.tight_layout()
    save(fig, ctx, "pyrolysis-literature-zhang-massflux")


# ------------------------------------------------------------------------------------ exact benchmarks of type 2
EXACT = {"pressure_transient": (("Gate_BenchPressureTransient", "bench_pressure_transient"), ("N20", "N40", "N80")),
         "transpiration": (("Gate_BenchTranspiration", "bench_transpiration"), ("N20", "N40", "N80")),
         "sharp_front": (("Gate_BenchSharpFront", "bench_sharp_front"), ("N50", "N100", "N200"))}


def cells_axis(ax, ns, label="cells N"):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.set_xticklabels([str(n) for n in ns])
    ax.minorticks_off()
    ax.set_xlim(ns[0] / 1.25, ns[-1] * 1.6)
    ax.set_xlabel(label)
    log_ticks(ax)


def page_darcy_exact(ctx):
    vx = load("vv_verify_exact", "exact/verify_exact.py")
    ec = vx.ec
    dx = ec.dx
    marks = ["o", "s", "^"]

    def setup(kind):
        g, names = EXACT[kind]
        cases = [need(gate(ctx["runs"], g, n)) for n in names]
        res = vx.measure(kind, cases)
        order = json.load(open(os.path.join(PYRO, "exact", kind, "expected.json")))["order"]
        return cases, res, order

    def mesh_marker(ax, k, xs, ys, n, ms=4.0):
        ax.plot(xs, ys, ls="none", marker=marks[k], ms=ms, mfc="white" if k < 2 else SERIES[k], mec=SERIES[k],
                mew=1.0, label="FUSS, %d cells" % n, zorder=3)

    # 1: transient gas pressure, linear limit
    cases, res, order = setup("pressure_transient")
    print("  pressure_transient: max|p - p_ex|/dp %s" % ", ".join("N %d %.4e" % (r["n"], r["e"]) for r in res))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.2), gridspec_kw={"width_ratios": [1.3, 1]})
    for k, (c, r) in enumerate(zip(cases, res)):
        meta = json.load(open(os.path.join(c, "case.json")))
        snaps = ec.snapshots(c)
        t, _f, dims, d = min(snaps, key=lambda q: abs(q[0] / meta["tau_p"] - 0.5))
        xc, pp = ec.centres(dims, d), ec.row_i(dims, d, "p_pg")
        dp = meta["p0"] - meta["pw"]
        mesh_marker(a1, k, [x / meta["L"] for x in xc], [(v - meta["pw"]) / dp for v in pp], meta["ncells"],
                    ms=4.5 - k)
    xs = [meta["L"] * i / 300.0 for i in range(301)]
    a1.plot([x / meta["L"] for x in xs], [(dx.pressure_transient_linear(x, t, meta["L"], meta["D"], meta["p0"],
                                                                         meta["pw"]) - meta["pw"]) / dp for x in xs],
            color=INK, lw=1.2, label="exact series", zorder=2)
    a1.set_xlim(0, 1)
    a1.set_xlabel("x / L (sealed at 0, open at 1)")
    a1.set_ylabel(r"$(p - p_w)/(p_0 - p_w)$")
    a1.set_title(r"$t \approx 0.5\,L^2/D$")
    a1.legend(loc="lower left")
    ns = [r["n"] for r in res]
    ladder(a2, ns, [r["e"] for r in res], order["e"], SERIES[0], "o", r"max $|p - p_\mathrm{ex}|/\Delta p$",
           offset=(9, 0), every=False)
    cells_axis(a2, ns)
    a2.set_ylabel("error")
    ladder_legend(a2, loc="upper right", fontsize=7.5)
    fig.tight_layout(w_pad=1.5)
    save(fig, ctx, "darcy-exact-pressure-transient")

    # 2: steady transpiration
    cases, res, order = setup("transpiration")
    print("  transpiration: max|T - T_ex| %s K; mdot error %s" % (
        ", ".join("N %d %.4e" % (r["n"], r["e"]) for r in res), ", ".join("%.3e" % r["emdot"] for r in res)))
    fig, axs = plt.subplots(1, 3, figsize=(FW, 3.0), gridspec_kw={"width_ratios": [1.35, 1, 1]})
    for k, c in enumerate(cases):
        meta = json.load(open(os.path.join(c, "case.json")))
        t, _f, dims, d = ec.snapshots(c)[-1]
        mesh_marker(axs[0], k, [x * 1e3 for x in ec.centres(dims, d)], ec.row_i(dims, d, "T"), meta["ncells"],
                    ms=4.5 - k)
    md, A, B, Pe, lin = dx.transpiration_steady(meta["L"], meta["k"], meta["cpg"], meta["W"], meta["mu"], meta["K"],
                                                meta["pin"], meta["pout"], meta["h1"], meta["T1"], meta["h2"],
                                                meta["T2"])
    xs = [meta["L"] * i / 300.0 for i in range(301)]
    axs[0].plot([x * 1e3 for x in xs], [dx.transpiration_T(x, meta["L"], A, B, Pe, lin) for x in xs], color=INK,
                lw=1.2, label="exact", zorder=2)
    axs[0].set_xlim(0, meta["L"] * 1e3)
    axs[0].set_xlabel("x (mm), gas flowing to +x")
    axs[0].set_ylabel("T (K)")
    axs[0].set_title("steady state, Pe = %.3f" % Pe)
    axs[0].legend(loc="upper left", fontsize=7.5)
    ns = [r["n"] for r in res]
    ladder(axs[1], ns, [r["e"] for r in res], order["e"], SERIES[0], "o", r"max $|T - T_\mathrm{ex}|$",
           offset=(9, 0), every=False)
    ladder(axs[2], ns, [r["emdot"] for r in res], order["emdot"], SERIES[1], "s",
           r"$|\dot m - \dot m_\mathrm{ex}|/\dot m_\mathrm{ex}$", offset=(9, 0), every=False)
    for ax, yl in ((axs[1], "error (K)"), (axs[2], "relative error")):
        cells_axis(ax, ns)
        ax.set_ylabel(yl)
        ladder_legend(ax, loc="upper right", fontsize=7.5)
    fig.tight_layout(w_pad=1.0)
    save(fig, ctx, "darcy-exact-transpiration")

    # 3: decomposition front with transpiration through the char
    cases, res, order = setup("sharp_front")
    print("  sharp_front: L2 T error %s K; front error %s" % (
        ", ".join("N %d %.4f" % (r["n"], r["e"]) for r in res), ", ".join("%.3e" % r["exf"] for r in res)))
    fig, axs = plt.subplots(2, 2, figsize=(FW, 5.4))
    for k, c in enumerate(cases):
        meta, rows = ec.compare(c)
        snaps = [(t, e) for t, e in rows if isinstance(t, float)]
        mesh_marker(axs[0, 0], k, [t for t, _e in snaps], [e["xf"] * 1e3 for _t, e in snaps], meta["ncells"],
                    ms=4.5 - k)
        t, _f, dims, d = ec.snapshots(c)[-1]
        mesh_marker(axs[0, 1], k, [x * 1e3 for x in ec.centres(dims, d)], ec.row_i(dims, d, "T"), meta["ncells"],
                    ms=4.5 - k)
    tt = [t * i / 300.0 for i in range(301)]
    axs[0, 0].plot(tt, [dx.sharp_front_position(v, meta["lam"], meta["alpha_c"]) * 1e3 for v in tt], color=INK,
                   lw=1.2, label="exact", zorder=2)
    axs[0, 0].axvspan(0, 5.0, color=BAND, alpha=0.3, lw=0, zorder=0)
    axs[0, 0].text(2.5, 0.97, "not\naveraged", transform=axs[0, 0].get_xaxis_transform(), ha="center", va="top",
                   fontsize=7.5, color=INK2)
    axs[0, 0].set_xlim(0, t)
    axs[0, 0].set_ylim(0, None)
    axs[0, 0].set_xlabel("t (s)")
    axs[0, 0].set_ylabel("front depth (mm)")
    axs[0, 0].legend(loc="lower right", fontsize=7.5)
    xf = dx.sharp_front_position(t, meta["lam"], meta["alpha_c"])
    xs = [meta["L"] * 0.5 * i / 400.0 for i in range(401)]
    axs[0, 1].plot([x * 1e3 for x in xs], [dx.sharp_front_T(x, t, meta["lam"], meta["alpha_c"], meta["alpha_v"],
                                                            meta["rho_v"], meta["rho_c"], meta["cp"], meta["cpg"],
                                                            meta["Ts"], meta["Tp"], meta["T0"]) for x in xs],
                   color=INK, lw=1.2, label="exact", zorder=2)
    axs[0, 1].axvline(xf * 1e3, color=INK2, lw=0.8, ls=(0, (4, 3)))
    axs[0, 1].text(xf * 1e3, 0.97, " front", transform=axs[0, 1].get_xaxis_transform(), fontsize=7.5, color=INK2,
                   va="top")
    axs[0, 1].set_xlim(0, meta["L"] * 0.5 * 1e3)
    axs[0, 1].set_xlabel("x from the heated surface (mm)")
    axs[0, 1].set_ylabel("T (K)")
    axs[0, 1].set_title("t = %.0f s" % t)
    axs[0, 1].legend(loc="upper right", fontsize=7.5)
    ns = [r["n"] for r in res]
    ladder(axs[1, 0], ns, [r["e"] for r in res], order["e"], SERIES[0], "o", r"$L_2$ T error, time average",
           offset=(9, 0), every=False)
    ladder(axs[1, 1], ns, [r["exf"] for r in res], order["exf"], SERIES[1], "s", "front error, time average",
           offset=(9, 0), every=False)
    for ax, yl in ((axs[1, 0], "error (K)"), (axs[1, 1], "relative error")):
        cells_axis(ax, ns)
        ax.set_ylabel(yl)
        ladder_legend(ax, loc="upper right", fontsize=7.5)
    fig.tight_layout(h_pad=1.2, w_pad=1.5)
    save(fig, ctx, "darcy-exact-sharp-front")


# ------------------------------------------------------------------------------- ablating wall: ATC 2.1 against PATO
ATC21 = ("Gate_AblationAtc21", "ablation_atc21")


def page_ablating_wall(ctx):
    av = load("vv_atc21", "ablation/ATC21/verify.py")
    at = av.at
    case = need(gate(ctx["runs"], ATC21))
    ref = os.path.join(PYRO, "ablation", "ATC21", "ref")
    Pp, Ps, Pm = (at.pato_rows(need(os.path.join(ref, f))) for f in ("Ta_plot", "Ta_surfacePatch", "mass"))
    E = at.pato_rows(need(os.path.join(ref, "Amaryllis_Energy_TestCase_2.1.txt")))
    M = at.pato_rows(need(os.path.join(ref, "Amaryllis_Mass_TestCase_2.1.txt")))
    n = min(len(Pp), len(Ps), len(Pm), len(E), len(M))
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    px = {p["name"]: p["x"] for p in pj["probes"]}
    fT = {}                                   # verify.py: the two bracketing probes, linear in x
    for d in av.DEPTHS_MM:
        lo, hi = pj["pairs"][str(d)]
        a = at.pato_rows(need(os.path.join(case, "OUTPUT", lo + ".txt")))
        b = at.pato_rows(need(os.path.join(case, "OUTPUT", hi + ".txt")))
        w = ((L - d * 1e-3) - px[lo]) / (px[hi] - px[lo])
        fT[d] = [(r[0], (1 - w) * r[1] + w * q[1]) for r, q in zip(a, b)]
    seb = at.read_seb(case)
    faces = at.faces_306(case)
    key0 = sorted(faces)[0]
    kinks = [row[0] for row in at.read_env(need(os.path.join(case, "INPUT", faces[key0][3])))[1:-1]]
    s0 = [r for r in seb if (r["b"], r["i"], r["j"], r["k"], r["f"]) == key0]
    fTw = [(r["time"], r["T_w"]) for r in s0]
    fmg = [(r["time"], r["mdot_pg"]) for r in s0]
    fr = sorted(at.fronts(f, L, av.RHO_V, av.RHO_C) for f in glob.glob(os.path.join(case, "OUTPUT", "field*.tec")))
    # the gate's maxima and the PATO-Amaryllis spreads, on PATO's grid (verify.py)
    res, spread = {}, {}

    def upd(dct, k, dv, t, where=""):
        if dv is not None and abs(dv) > dct.get(k, (-1.0,))[0]:
            dct[k] = (abs(dv), t, where)
    fv, fc = [(t, v) for t, v, c in fr], [(t, c) for t, v, c in fr]
    for i in range(n):
        t = Pp[i][0]
        if t < 0.05:
            continue
        for k, d in enumerate(av.DEPTHS_MM):
            v = at.interp(fT[d], t)
            if v is not None:
                upd(res, "T", v - Pp[i][1 + k], t, "%d mm" % d)
            if E[i][2 + k] > 0:
                upd(spread, "T", Pp[i][1 + k] - E[i][2 + k], t, "%d mm" % d)
        v98, c02 = at.interp(fv, t), at.interp(fc, t)
        upd(res, "v98", None if v98 is None else v98 - Pm[i][3], t)
        upd(res, "c02", None if c02 is None else c02 - Pm[i][4], t)
        upd(spread, "v98", Pm[i][3] - M[i][3], t)
        upd(spread, "c02", Pm[i][4] - M[i][4], t)
        upd(spread, "Tw", Ps[i][1] - E[i][1], t)
        upd(spread, "mdot_g", Pm[i][1] - M[i][1], t)
        br = at.bracket(fTw, t)              # a grid time whose FUSS samples straddle an env breakpoint is skipped
        if br is not None and not any(br[0] < tk < br[1] or tk == t for tk in kinks):
            upd(res, "Tw", at.interp(fTw, t) - Ps[i][1], t)
            upd(res, "mdot_g", at.interp(fmg, t) - Pm[i][1], t)
    for q in ("T", "Tw", "mdot_g", "v98", "c02"):
        print("  (%s) max |FUSS - PATO| %.4g at %g s %s; the PATO-Amaryllis spread %.4g at %g s %s"
              % ((q,) + res[q] + spread[q]))
    rows = [i for i in range(n) if Pp[i][0] > 0.0]
    every = [i for i in rows if abs(Pp[i][0] / 4.0 - round(Pp[i][0] / 4.0)) < 1e-6]          # PATO markers
    every_a = [i for i in rows if abs((Pp[i][0] - 2.0) / 4.0 - round((Pp[i][0] - 2.0) / 4.0)) < 1e-6]   # Amaryllis
    srcs = [("FUSS", "-", None, True), ("PATO 3.1", "none", "o", True), ("Amaryllis", "none", "s", False)]

    # 1: in-depth temperatures
    fig, ax = plt.subplots(figsize=(FW, 4.0))
    labels = []
    for k, d in enumerate(av.DEPTHS_MM):
        col = SERIES[k]
        ax.plot([t for t, _ in fT[d]], [v for _, v in fT[d]], color=col, lw=1.5, zorder=3)
        ax.plot([Pp[i][0] for i in every], [Pp[i][1 + k] for i in every], ls="none", marker="o", ms=4, color=col,
                zorder=4)
        ax.plot([E[i][0] for i in every_a], [E[i][2 + k] for i in every_a], ls="none", marker="s", ms=4, mfc="white",
                mec=col, mew=1.0, zorder=4)
        labels.append((max(v for _, v in fT[d]), mm(d), col))
    tend = fT[av.DEPTHS_MM[0]][-1][0]
    ax.set_xlim(0, tend)
    ax.set_ylim(250, None)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("T (K)")
    source_legend(ax, srcs, loc="upper right")
    for y, text, col in labels:
        tmax = max(fT[int(text.split()[0])], key=lambda r: r[1])[0]
        ax.annotate(text, (tmax, y), xytext=(-4, 3), textcoords="offset points", ha="right", va="bottom",
                    fontsize=8, color=INK)
    fig.tight_layout()
    save(fig, ctx, "ablating-wall-atc21-temperatures")

    # 2: wall temperature and gas blowing rate (log time: the spreads of both peak in the first second)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.3))
    for ax, f, col_p, col_e, scale, yl in ((a1, fTw, (Ps, 1), (E, 1), 1.0, r"$T_w$ (K)"),
                                           (a2, fmg, (Pm, 1), (M, 1), 1e3, r"$\dot m_g$ (g/(m$^2$s))")):
        ax.plot([t for t, _ in f if t > 0], [v * scale for t, v in f if t > 0], color=SERIES[0], lw=1.5, zorder=3)
        logsel = []
        for g in [0.1 * (tend / 0.1) ** (j / 79.0) for j in range(80)]:      # marker times even in log t
            i = min(rows, key=lambda r: abs(math.log(Pp[r][0]) - math.log(g)))
            if i not in logsel:
                logsel.append(i)
        ip, ie = logsel[0::2], logsel[1::2]          # PATO and Amaryllis on alternate marker times
        ax.plot([col_p[0][i][0] for i in ip], [col_p[0][i][col_p[1]] * scale for i in ip], ls="none",
                marker="o", ms=4, color=INK2, zorder=4)
        ax.plot([col_e[0][i][0] for i in ie], [col_e[0][i][col_e[1]] * scale for i in ie], ls="none",
                marker="s", ms=4.5, mfc="white", mec=INK2, mew=1.0, zorder=4)
        log_time_axis(ax, 0.1, tend)
        ax.set_ylabel(yl)
        if ax is a2:
            ax.set_ylim(0, None)
    source_legend(a1, [("FUSS", "-", None, True), ("PATO 3.1", "none", "o", True), ("Amaryllis", "none", "s", False)],
                  loc="lower left")
    fig.tight_layout(w_pad=1.5)
    save(fig, ctx, "ablating-wall-atc21-wall")

    # 3: the virgin-98 % and char-2 % front depths
    fig, ax = plt.subplots(figsize=(FW, 3.3))
    for k, (col_f, col_r, lab) in enumerate(((1, 3, "virgin 98 %"), (2, 4, "char 2 %"))):
        col = SERIES[k]
        ax.plot([r[0] for r in fr], [r[col_f] * 1e3 for r in fr], color=col, lw=1.5, zorder=3, label=lab)
        ax.plot([Pm[i][0] for i in every], [Pm[i][col_r] * 1e3 for i in every], ls="none", marker="o", ms=4, color=col,
                zorder=4)
        ax.plot([M[i][0] for i in every_a], [M[i][col_r] * 1e3 for i in every_a], ls="none", marker="s", ms=4,
                mfc="white", mec=col, mew=1.0, zorder=4)
    ax.set_xlim(0, tend)
    ax.set_ylim(0, None)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("front depth below the surface (mm)")
    l1 = ax.legend(loc="upper left")
    ax.add_artist(l1)
    source_legend(ax, srcs, loc="lower right")
    fig.tight_layout()
    save(fig, ctx, "ablating-wall-atc21-fronts")


PAGES = [("pyrolysis-two-layer", page_two_layer), ("pyrolysis-atc1", page_atc1),
         ("pyrolysis-atc2", page_atc2),
         ("pyrolysis-kinetics-0d", page_kinetics),
         ("pyrolysis-energy-audit", page_energy_audit),
         ("pyrolysis-literature", page_literature),
         ("darcy-exact", page_darcy_exact),
         ("ablating-wall", page_ablating_wall)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", default=os.path.join(REPO, "build", "test", "runs", "gates"))
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "vv", "images"))
    ap.add_argument("--only", default="")
    ap.add_argument("--png", default="")
    a = ap.parse_args()
    want = [w for w in a.only.split(",") if w] or [p for p, _f in PAGES]
    unknown = [w for w in want if w not in dict(PAGES)]
    if unknown:
        sys.exit("plot_vv.py: unknown page %s; the pages are %s"
                 % (", ".join(unknown), ", ".join(p for p, _f in PAGES)))
    for d in (a.out, a.png):
        if d and not os.path.isdir(d):
            os.makedirs(d)
    ctx = {"runs": os.path.abspath(a.runs), "out": a.out, "png": a.png}
    status = 0
    for name, fn in PAGES:
        if name not in want:
            continue
        print("%s:" % name)
        try:
            fn(ctx)
        except Missing as m:
            print("  SKIP %s: no %s" % (name, m))
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
