#!/usr/bin/env python3
"""
verify.py -- gate of the Baum & Atreya (2007) charring-front benchmark (make_baum_atreya.py cases), 2026-10-13.

Usage: verify.py <full-case> [<full-case> ...] [--qs <qs-case>]
  full-cases: one mesh ladder (coarsest first) of --mode full runs; --qs: one --mode qs run (quasi-steady, the paper's
  beta = 0). The expected orders and the round-off floor of the full ledger in expected.json beside this script.
Exact solution (corrected front condition; darcy_exact.sharp_front_*): T(x, t), front X = C* sqrt(alpha_c t), and the
quasi-steady char pressure P^2(eta*) = 1 + 2 M int_0^eta* (T/Tp) d eta', M = (C*/2)(rv - rc) ac mu (R/W) Tp/(K_c ps^2)
(constant mu), constant in time at fixed eta* = x/sqrt(alpha_c t).
Checks, averages over 30 s <= t <= t_end, the first two by the convergence rule of tools/conv_rule.py (the error
falls on every rung and its order on the finest pair lies in [0.8 p, 1.2 p]):
  (front)  mass-based front error, p = 1 (a cell-by-cell front is first order by design)
  (T)      L2 T error (char + preheat, x <= 3X), p = 1 (the error follows the front)
  (p)      printed, not gated: <p^2 - ps^2> at the cell nearest 1 mm, over the probe times it lies in the char, against
           the exact (that cell moves with the mesh, 1.0 mm on N100, 1.1 mm on N200: not a fixed-point error)
  (mass)   every full run: |gas out + stored - stored(0) - solid lost| / solid lost <= ledger_full (round-off)
  (qs)     the quasi-steady run: |gas out - solid lost| / solid lost <= the one-step share of the run, (the largest
           solid loss between two budget rows) / (the steps per row) / (the total solid lost), the steps per row taken
           uniform (the run's steps, OUTPUT/residual-history.dat, over its budget intervals): the quasi-steady closure
           carries one step of gas (11 % before 2026-10-13, when it dropped the gas the stages released beyond the
           prediction)
FUSS_GATE_SELFTEST=1 empties the order windows and zeroes the ledger bounds. Python 3.6, no numpy.
"""
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import conv_rule as cr        # noqa: E402
import darcy_tools as dt      # noqa: E402
import darcy_exact as dx      # noqa: E402
import fuss_out as fo         # noqa: E402

R_UNIV = 8.314471469
AREA = 0.02                   # 2 x 0.1 m x 0.1 m: the cross-section of every 1-D pyrolysis case
TMIN = 30.0


def exact_pressure(m):
    """eta* nodes and P^2 on [0, C*] (Simpson)."""
    lam, ac, av = m["lam"], m["alpha_c"], m["alpha_v"]
    Cs = 2 * lam
    M = (Cs / 2) * (m["rv"] - m["rc"]) * ac * m["mu_p"] * (R_UNIV / m["W"]) * m["Tp"] / (m["Kc"] * m["ps"] ** 2)
    N = 2000
    etas = [Cs * i / N for i in range(N + 1)]
    G = lambda e: dx.sharp_front_T(e * math.sqrt(ac), 1.0, lam, ac, av, m["rv"], m["rc"], m["cc"], m["cg"], m["Ts"], m["Tp"], m["T0"]) / m["Tp"]
    I = [0.0]
    for i in range(N):
        a, b = etas[i], etas[i + 1]
        I.append(I[-1] + (b - a) / 6 * (G(a) + 4 * G(0.5 * (a + b)) + G(b)))
    return Cs, etas, [1 + 2 * M * v for v in I]


def p2_exact(x, t, m, Cs, etas, P2):
    e = x / math.sqrt(m["alpha_c"] * t)
    if e >= Cs:
        return P2[-1]
    k = min(int(e / Cs * (len(etas) - 1)), len(etas) - 2)
    w = (e - etas[k]) / (etas[k + 1] - etas[k])
    return (1 - w) * P2[k] + w * P2[k + 1]


def front_and_T(case, m):
    lam, ac, av = m["lam"], m["alpha_c"], m["alpha_v"]
    exf, eL2 = [], []
    for f in glob.glob(os.path.join(case, "OUTPUT", "field*.tec")):
        dims, t, d = dt.read_tec(f)
        if t < TMIN:
            continue
        ni = dims[0] - 1
        x = d["x"]
        xc = [0.5 * (x[i] + x[i + 1]) for i in range(ni)]
        h = xc[1] - xc[0]
        X = m["Cstar"] * math.sqrt(ac * t)
        xf = sum(m["rv"] - r for r in d["rho_s"][:ni]) * h / (m["rv"] - m["rc"])
        exf.append(abs(xf - X) / X)
        win = [i for i in range(ni) if xc[i] <= 3 * X]
        e = [d["T"][i] - dx.sharp_front_T(xc[i], t, lam, ac, av, m["rv"], m["rc"], m["cc"], m["cg"], m["Ts"], m["Tp"], m["T0"]) for i in win]
        eL2.append(math.sqrt(sum(v * v for v in e) / len(e)))
    return sum(exf) / len(exf), sum(eL2) / len(eL2), len(exf)


def char_pressure(case, m):
    Cs, etas, P2 = exact_pressure(m)
    ini = open(os.path.join(case, "input.ini")).read()
    if "[x10]" not in ini or not os.path.exists(os.path.join(case, "OUTPUT", "x10.txt")):
        return None, None          # a case without the 1 mm probe (not a make_baum_atreya.py case of 2026-10-13)
    i = int(ini.split("[x10]")[1].split("index-position = 1 ")[1].split()[0])
    dxc = m["L"] / m["ncells"]
    xc = (i - 0.5) * dxc
    rows = [r for r in fo.probe(case, "x10") if r[0] >= TMIN and xc + 1.5 * dxc <= Cs * math.sqrt(m["alpha_c"] * r[0])]
    ps = m["ps"]
    f = sum(r[2] ** 2 / ps ** 2 - 1.0 for r in rows) / len(rows)
    e = sum(p2_exact(xc, r[0], m, Cs, etas, P2) - 1.0 for r in rows) / len(rows)
    return f / e - 1.0, xc


def one_step_share(case):
    """(largest solid loss between two budget rows) / (steps per row) / (total solid lost), steps uniform per row"""
    b = fo.budget(case)
    with open(os.path.join(case, "OUTPUT", "residual-history.dat")) as fh:
        steps = sum(1 for line in fh if line.strip())
    loss = max(b[i - 1][2] - b[i][2] for i in range(1, len(b)))
    return loss / (steps / float(len(b) - 1)) / (b[0][2] - b[-1][2]), steps, len(b)


def ledger(case, with_storage):
    b = fo.budget(case)
    r0, r = b[0], b[-1]
    lost = r0[2] - r[2]
    stored = (r[7] - r0[7]) if with_storage else 0.0
    return (r[5] + stored - lost) / lost, r[0]


def main():
    argv = sys.argv[1:]
    qs = None
    if "--qs" in argv:
        qs = os.path.abspath(argv[argv.index("--qs") + 1])
    full = [os.path.abspath(a) for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or argv[i - 1] != "--qs")]
    if not full:
        sys.exit(__doc__)
    path = os.path.join(HERE, "expected.json")
    exp = json.load(open(path))
    sel = dt.selftest()
    lim_full = 0.0 if sel else exp["ledger_full"]
    ok = True
    meas = []
    print("Baum & Atreya 2007 (corrected): full type 2 %s%s" % (" ".join(os.path.basename(c) for c in full),
                                                            ", quasi-steady %s" % os.path.basename(qs) if qs else ""))
    for c in full:
        m = json.load(open(os.path.join(c, "case.json")))
        exf, eL2, ns = front_and_T(c, m)
        lg, tend = ledger(c, True)
        good = abs(lg) <= lim_full
        ok &= good
        meas.append({"n": m["ncells"], "exf": exf, "eL2": eL2, "ledger": lg})
        print("  %4d cells, c_g %g, C* %.6f: front error %.3e  L2 T error %.3f K  (%d snapshots)  (mass) ledger %+.2e at %.1f s  %s"
              % (m["ncells"], m["cg"], m["Cstar"], exf, eL2, ns, lg, tend, "PASS" if good else "FAIL"))
    b = meas[-1]
    m = json.load(open(os.path.join(full[-1], "case.json")))
    prel, xc = char_pressure(full[-1], m)
    b["p_rel"] = prel
    if prel is None:
        prel, xc = float("inf"), float("nan")
        print("  (p) the finest case has no x10 probe output")
    checks = []
    if len(meas) < 2:
        checks.append(("(front), (T): the convergence rule needs a ladder of 2 or more full runs", False))
    else:
        ns = [r["n"] for r in meas]
        for lab, key in (("(front) error", "exf"), ("(T) L2 error [K]", "eL2")):
            good, text = cr.check(lab, [r[key] for r in meas], ns, exp["order"][key], fmt="%.4e")
            checks.append((text, good))
    print("  (p) <p^2 - ps^2> at x = %.3f mm: %+.2f %% of exact (printed, not gated)" % (1e3 * xc, 100 * prel))
    if qs:
        lq, tq = ledger(qs, False)
        share, nst, nrow = one_step_share(qs)
        lim_qs = 0.0 if sel else share
        checks.append(("(qs) quasi-steady ledger %+.2e at %.1f s (|.| <= %.3e, the one-step share: %d steps over %d budget rows)"
                       % (lq, tq, lim_qs, nst, nrow), abs(lq) <= lim_qs))
    for text, good in checks:
        ok &= good
        print("  %s  %s" % (text, "PASS" if good else "FAIL"))
    print("baum_atreya2007: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
