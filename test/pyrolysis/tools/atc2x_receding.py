"""Plan 08 section 3.3: Ablation Test Cases 2.2 and 2.3 -- FUSS's ablating wall (BC 306) on the receding mesh (law
nodeshift) against a local PATO run (shared by test/pyrolysis/ablation/ATC22/verify.py and ATC23/verify.py).

PATO 3.1's AblationTestCase_2.x (movingMesh yes) at maxDeltaT 0.0005 is the oracle (ref/README.md of each case; at the
tutorial's 0.005 s PATO's own time error dominated ATC 2.1, plan 08 stage A). FUSS samples, on PATO's 0.1 s grid
t = 0.1 .. 120 s, from its own outputs:

  T at 1, 2, 4, 8, 12, 16, 24 mm below the INITIAL surface (DECISION X-D4): every cell moves under nodeshift, so a
      FUSS probe (a cell) cannot do it; each field file (one every 0.1 s) gives T at x = L0 - d, linear between the two
      CURRENT cell centres that bracket it; then linear in t. A depth is compared while both codes still have it: PATO
      writes 0 once its surface has passed the point, Amaryllis 0 or a stray value below 1 K (0.01191 at 8 mm, 40.7 s,
      of 2.2), so any value below 1 K counts as consumed; and FUSS stops comparing a depth once it lies in the last
      half cell, above the surface cell's centre (no two centres bracket it). There each code's in-cell rule decides
      the value, not the solution: PATO's cellPoint probe tracks the surface cell's centre value (2.2, 1 mm, 26 um under
      the surface: PATO 1445.5 K, FUSS's surface cell 1442.9 K, FUSS's centre-to-T_w line 1459.5 K; 2.3 at 22 um: 3073.7,
      3073.3, 3096.4 K), so the centre-to-wall line read 16-23 K of rule difference in the 0.3 s before each depth is
      consumed (measured 2026-10-13).
  T_w, mdot_g, mdot_c   seb.txt rows of the first 306 face, linear in t.
  s(t) = L0 - x_w, the virgin-98 % and char-2 % fronts   from each field file (the fronts by PATO's rule,
      abl_tools.fronts, depth from the INITIAL surface), linear in t.
  A grid time whose two FUSS samples straddle a breakpoint of the env table (0.1, 60, 60.1 s: a kink) is skipped for
  every quantity, as on ATC 2.1: across the switch-off at 60 s the samples are 0.1 s apart while the wall cools by
  140 K, so linear interpolation reads 56 K of sampling error into T at 12 mm (0.13 mm under the surface) on 2.2.

Checks (expected.json; DECISION 2026-10-13 P-D4 of plan 08):
  (T), (Tw)                      max |FUSS - PATO| <= measured + 2 K
  (mdot_g) (mdot_c) (s) (v98) (c02)   max |FUSS - PATO| <= measured + 10 % of PATO's peak of the quantity
  (seb)    plan 08 section 3.2 on every seb.txt row (abl_tools.seb_closure)
  (bal-m)  the type-2 mass M_s + M_g + M_gas_out + M_abl + G_abl of energy-budget.txt constant to 1e-12 of itself
  (bal-e)  E_tot + E_abl + Q_out - Q_src constant to 5e-10 of max(|E_tot(0)|, max |Q_out|): the remap's swept energy
           E_abl closes the budget of a receding wall. The drift is per-step round-off, not a time-step error: 5.6e-11
           (2.2) and 4.1e-11 (2.3) at 120 s, growing as t^0.5, and over the first 10 s of 2.2 it grows as the step
           shrinks (1.4e-11, 2.0e-11, 3.1e-11 at vnn 0.5, 0.25, 0.125); an O(dt) mismatch would shrink. A swept
           energy missing or valued inconsistently with the remap shows at 1e-3 and above.
  (R2)     printed, not gated: e_front = max_t |x_98,FUSS - x_98,PATO| against the FUSS cell width at the front then;
           plan 08 decision 4 adds a limited (MUSCL) reconstruction iff e_front > dx_local.
The 1.5x target of plan 06 section 4.4 (max |FUSS - PATO| <= 1.5 max |PATO - Amaryllis|) is printed per quantity and
recorded in the case README; it does not gate (P-D4).

FUSS_GATE_SELFTEST=1 sets every budget to 0. $FUSS_PATO_ATC2X_OUT (a PATO output/porousMat directory) replaces the
PATO files of ref/ (not the Amaryllis ones). Python 3.6, no numpy.
"""
import bisect
import datetime
import glob
import json
import os
import re
import sys

import abl_tools as at

DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]
L0 = 0.05                            # m, the slab (make_case.py ATC_L)
RHO_V, RHO_C = 280.0, 220.0          # TACOT virgin and char solid densities (INPUT/pyrolysis.dat)
T_SLACK = 2.0                        # K, P-D4
FRAC_SLACK = 0.10                    # of PATO's peak, P-D4
BAL_M, BAL_E = 1.0e-12, 5.0e-10
T_GONE = 1.0                         # K: a probe value below this is a consumed point (PATO 0, Amaryllis 0 or 0.0119)
QUANTS = ("T", "Tw", "mdot_g", "mdot_c", "s", "v98", "c02")
UNITS = {"T": "K", "Tw": "K", "mdot_g": "kg/m2/s", "mdot_c": "kg/m2/s", "s": "m", "v98": "m", "c02": "m"}


def read_field(path):
    """(t, x nodes of the j = k = 1 line, {var: cell values of block 1 along that line}) of a tecplot-ascii field file."""
    toks = open(path).read().split()
    k = toks.index("ZONE")
    names = " ".join(toks[:k]).split("=", 1)[1].split()
    s = k
    while not toks[s].startswith("SOLUTIONTIME="):
        s += 1
    zone = " ".join(toks[k:s + 1])
    I = int(re.search(r"\bI=(\d+)", zone).group(1))
    J = int(re.search(r"\bJ=(\d+)", zone).group(1))
    K = int(re.search(r"\bK=(\d+)", zone).group(1))
    t = float(re.search(r"SOLUTIONTIME=([-+0-9.eEdD]+)", zone).group(1).replace("D", "E"))
    s += 1
    nn = I * J * K
    nc = (I - 1) * max(J - 1, 1) * max(K - 1, 1)
    x = [float(v) for v in toks[s:s + I]]
    cells = {}
    for nm in ("T", "rho_s"):
        off = s + 3 * nn + (names.index(nm) - 3) * nc
        cells[nm] = [float(v) for v in toks[off:off + I - 1]]
    return t, x, cells


def t_at_depth(x, Tc, d):
    """T at x = L0 - d on the current mesh, linear between the two cell centres that bracket it, or None once the
    point lies above the surface cell's centre (or the surface has passed it)."""
    xd = L0 - d
    xc = [0.5 * (x[i] + x[i + 1]) for i in range(len(x) - 1)]
    if xd >= xc[-1]:
        return None
    i = bisect.bisect_right(xc, xd)
    return Tc[i - 1] + (Tc[i] - Tc[i - 1]) * (xd - xc[i - 1]) / (xc[i] - xc[i - 1])


def dx_at(x, xp):
    """The width of the current cell holding x = xp."""
    i = min(max(bisect.bisect_right(x, xp) - 1, 0), len(x) - 2)
    return x[i + 1] - x[i]


def run(here, label, amaryllis_tag, argv):
    record = "--record" in argv
    args = [a for a in argv if a != "--record"]
    case = os.path.abspath(args[0]) if args else here
    ref = os.path.join(here, "ref")
    pato = os.environ.get("FUSS_PATO_ATC2X_OUT") or ref
    exp_path = os.path.join(here, "expected.json")
    expected = json.load(open(exp_path)) if os.path.isfile(exp_path) else {}
    selftest = os.environ.get("FUSS_GATE_SELFTEST", "0") == "1"

    # ---- PATO and Amaryllis
    sub = (lambda f: os.path.join(pato, "scalar", f)) if os.path.isdir(os.path.join(pato, "scalar")) else \
          (lambda f: os.path.join(pato, f))
    Pp = at.pato_rows(sub("Ta_plot"))
    Ps = at.pato_rows(sub("Ta_surfacePatch"))
    Pm = at.pato_rows(os.path.join(pato, "mass"))
    E = at.pato_rows(os.path.join(ref, "Amaryllis_Energy_TestCase_%s.txt" % amaryllis_tag))
    M = at.pato_rows(os.path.join(ref, "Amaryllis_Mass_TestCase_%s.txt" % amaryllis_tag))
    n = min(len(Pp), len(Ps), len(Pm), len(E), len(M))

    # ---- FUSS: seb.txt (first 306 face) and the field files
    seb = at.read_seb(case)
    faces = at.faces_306(case)
    key0 = sorted(faces)[0]
    s0 = [r for r in seb if (r["b"], r["i"], r["j"], r["k"], r["f"]) == key0]
    fTw = [(r["time"], r["T_w"]) for r in s0]
    fmg = [(r["time"], r["mdot_pg"]) for r in s0]
    fmc = [(r["time"], r["mdot_ch"]) for r in s0]
    env = at.read_env(os.path.join(case, "INPUT", faces[key0][3]))
    kinks = [row[0] for row in env[1:-1]]
    fT = {d: [] for d in DEPTHS_MM}
    fs, fv, fc, fdx = [], [], [], []
    for path in glob.glob(os.path.join(case, "OUTPUT", "field*.tec")):
        t, x, cells = read_field(path)
        for d in DEPTHS_MM:
            fT[d].append((t, t_at_depth(x, cells["T"], d * 1e-3)))
        fs.append((t, L0 - x[-1]))
        tf, v98, c02 = at.fronts(path, L0, RHO_V, RHO_C)
        fv.append((t, v98))
        fc.append((t, c02))
        fdx.append((t, dx_at(x, L0 - v98) if v98 > 0 else None))
    for lst in [fs, fv, fc, fdx] + list(fT.values()):
        lst.sort(key=lambda p: p[0])

    def alive_interp(series, t):
        br = at.bracket(series, t)
        if br is None:
            return None
        lo = series[[p[0] for p in series].index(br[0])][1]
        hi = series[[p[0] for p in series].index(br[1])][1]
        if lo is None or hi is None:
            return None
        return at.interp([(br[0], lo), (br[1], hi)], t)

    # ---- maxima over PATO's grid
    res, spread, skipped = {}, {}, []

    def upd(dct, k, d, t, where=""):
        if d is not None and abs(d) > dct.get(k, (-1.0,))[0]:
            dct[k] = (abs(d), t, where)

    for i in range(n):
        t = Pp[i][0]
        if t < 0.05:
            continue
        for k, d in enumerate(DEPTHS_MM):
            if Pp[i][1 + k] > T_GONE and E[i][2 + k] > T_GONE:
                upd(spread, "T", Pp[i][1 + k] - E[i][2 + k], t, "%d mm" % d)
        upd(spread, "Tw", Ps[i][1] - E[i][1], t)
        for q, col in (("mdot_g", 1), ("mdot_c", 2), ("v98", 3), ("c02", 4), ("s", 5)):
            upd(spread, q, Pm[i][col] - M[i][col], t)
        brs = [at.bracket(fTw, t), at.bracket(fs, t)]
        if any(br is not None and any(br[0] < tk < br[1] or tk == t for tk in kinks) for br in brs):
            skipped.append(t)
            continue
        for k, d in enumerate(DEPTHS_MM):
            v = alive_interp(fT[d], t)
            if Pp[i][1 + k] > T_GONE and v is not None:
                upd(res, "T", v - Pp[i][1 + k], t, "%d mm" % d)
        if brs[0] is not None:
            upd(res, "Tw", at.interp(fTw, t) - Ps[i][1], t)
            upd(res, "mdot_g", at.interp(fmg, t) - Pm[i][1], t)
            upd(res, "mdot_c", at.interp(fmc, t) - Pm[i][2], t)
        for q, ser, col in (("s", fs, 5), ("v98", fv, 3), ("c02", fc, 4)):
            v = at.interp(ser, t)
            upd(res, q, None if v is None else v - Pm[i][col], t)
    peaks = {q: max(abs(r[col]) for r in Pm[:n]) for q, col in
             (("mdot_g", 1), ("mdot_c", 2), ("v98", 3), ("c02", 4), ("s", 5))}

    ok = True
    print("ATC %s gate: FUSS BC 306, receding mesh (nodeshift), vs PATO (%s), t = 0.1..120 s; %d grid times "
          "skipped (FUSS samples straddle an env breakpoint: %s)" % (label, pato, len(skipped),
                                                                    " ".join("%g" % t for t in skipped)))
    measured, budgets = {}, {}
    for q in QUANTS:
        m, t, where = res[q]
        measured[q] = {"max": m, "t": t, "where": where}
        if q in ("T", "Tw"):
            b = expected.get("budget", {}).get(q, m + T_SLACK)
        else:
            b = expected.get("budget", {}).get(q, m + FRAC_SLACK * peaks[q])
        if selftest:
            b = 0.0
        budgets[q] = b
        good = m <= b
        ok &= good
        tgt = 1.5 * spread[q][0]
        print("  (%-6s) max|FUSS - PATO| %.4g %s at %g s %s  budget %.4g  %s   [1.5x target %.4g: %s]" % (
            q, m, UNITS[q], t, where, b, "PASS" if good else "FAIL", tgt, "met" if m <= tgt else "NOT met"))

    # R2 (decision 4): the front error against the local cell size
    e_front, t_e, _ = res["v98"]
    dx_e = alive_interp(fdx, t_e)
    dxs = [p[1] for p in fdx if p[1] is not None]
    r2 = "MUSCL needed" if dx_e is not None and e_front > dx_e else "first order suffices"
    print("  (R2)     e_front = max|x98_FUSS - x98_PATO| = %.4g m at %g s; dx_local there %.4g m (front cells %.4g .. %.4g m)"
          " -> %s" % (e_front, t_e, dx_e if dx_e is not None else float("nan"), min(dxs), max(dxs), r2))

    # the budgets of energy-budget.txt close with the swept-out columns
    B = at.pato_rows(os.path.join(case, "OUTPUT", "energy-budget.txt"))
    tot = [r[2] + r[7] + r[5] + r[9] + r[10] for r in B]
    dm = max(abs(v - tot[0]) for v in tot) / tot[0]
    en = [r[1] + r[8] + r[3] - r[6] for r in B]
    scale = max(abs(B[0][1]), max(abs(r[3]) for r in B))
    de = max(abs(v - en[0]) for v in en) / scale
    lim_m, lim_e = (0.0, 0.0) if selftest else (BAL_M, BAL_E)
    good_m, good_e = dm <= lim_m, de <= lim_e
    ok &= good_m and good_e
    print("  (bal-m)  M_s + M_g + M_gas_out + M_abl + G_abl drift %.3e of the total (<= %.0e, %d rows)  %s"
          % (dm, lim_m, len(B), "PASS" if good_m else "FAIL"))
    print("  (bal-e)  E_tot + E_abl + Q_out - Q_src drift %.3e of max(|E0|, |Q_out|) (<= %.0e)  %s"
          % (de, lim_e, "PASS" if good_e else "FAIL"))

    tol = 0.0 if selftest else 1.0e-8
    good, lines, worst = at.seb_closure(case, tol if tol > 0 else 1e-30)
    ok &= good
    for l in lines:
        print("  " + l)

    if record:
        out = {
            "measured": measured,
            "budget": {q: (measured[q]["max"] + T_SLACK if q in ("T", "Tw") else measured[q]["max"] + FRAC_SLACK * peaks[q])
                       for q in measured},
            "pato_peaks": peaks,
            "spread_pato_amaryllis": {q: {"max": spread[q][0], "t": spread[q][1], "where": spread[q][2]} for q in spread},
            "target_1p5x_met": {q: measured[q]["max"] <= 1.5 * spread[q][0] for q in measured},
            "skipped_grid_times": skipped,
            "r2": {"e_front": e_front, "t": t_e, "dx_local": dx_e, "front_dx_min": min(dxs), "front_dx_max": max(dxs),
                   "decision": r2},
            "balances": {"mass": dm, "energy": de},
            "recorded": datetime.date.today().isoformat(),
            "note": "budget: T and Tw measured + 2 K, mdot_g, mdot_c, s and the fronts measured + 10 % of PATO's peak "
                    "(plan 08 P-D4); the 1.5x target is recorded, not gated; oracle ref/ = PATO 3.1 ATC 2.x (movingMesh "
                    "yes) at maxDeltaT 0.0005 (ref/README.md)",
        }
        with open(exp_path, "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("  recorded %s" % exp_path)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1
