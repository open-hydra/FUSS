#!/usr/bin/env python3
"""Plan 08 section 3.3, Ablation Test Case 2.1: the ablating wall BC 306 on a static mesh against a local PATO run.

The case is atc2x306 --case 2.1 of test/pyrolysis/tools/make_case.py: TACOT under the type-2 closure, face 2 BC 306 with
the equilibrium B' table and PATO's AblationTestCase_2.x boundary conditions (heating 0.1-60 s, chemistry off and pure
re-radiation to 300 K after 60.1 s), no-char = true. The oracle is PATO 3.1 on the same problem with movingMesh no
(ref/README.md): it zeroes B'c and mdot_c and keeps the table h_w, as FUSS's static 306 face does. PATO runs at
maxDeltaT 0.0005 s, a tenth of the tutorial's: at 0.005 s PATO's own time error reaches ~24 K at the wall after the
60 s switch-off (ref/pato-dt0.005/ keeps that run; ref/README.md has the step study).

FUSS samples, PATO's 0.1 s grid t = 0.1 .. 120 s:
  T at 1, 2, 4, 8, 12, 16, 24 mm   two bracketing index probes per depth (probes.json), linear in x then in t
  T_w, mdot_g                      seb.txt rows (the first 306 face), linear in t; a grid time whose two FUSS
                                   samples straddle a breakpoint of the env table (a kink of T_w: 0.1, 60, 60.1 s)
                                   is skipped -- linear interpolation across the switch-off at 60 s reads 72 K
                                   of sampling error into a wall that agrees to 0.1 K just before it
  virgin-98 %, char-2 % fronts     the field files (one every 0.1 s), PATO's own rule (abl_tools.fronts)
Checks (expected.json; DECISION 2026-10-13 P-D4 of plan 08):
  (T)      max over t and the 7 depths |T_FUSS - T_PATO| <= measured + 2 K
  (Tw)     max |T_w,FUSS - T_w,PATO| <= measured + 2 K
  (mdot_g) max |mdot_g,FUSS - mdot_g,PATO| <= measured + 10 % of PATO's peak mdot_g
  (v98), (c02)  the two fronts: max |depth_FUSS - depth_PATO| <= measured + 10 % of PATO's peak depth
  (static) mdot_ch, sdot and the recession are exactly 0 (no-char, static mesh)
  (seb)    plan 08 section 3.2 on every seb.txt row (abl_tools.seb_closure)
The 1.5x target of plan 06 section 4.4 (max |FUSS - PATO| <= 1.5 max |PATO - Amaryllis|, Amaryllis 2.1 in ref/) is
printed per quantity and recorded in the README; it does not gate (P-D4).

Usage: verify.py [case_dir] [--record]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets every budget to 0.
$FUSS_PATO_ATC21_OUT (a PATO output/porousMat directory) replaces the PATO files of ref/ (not the Amaryllis ones).
Python 3.6, no numpy.
"""
import datetime
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import abl_tools as at  # noqa: E402

DEPTHS_MM = [1, 2, 4, 8, 12, 16, 24]
RHO_V, RHO_C = 280.0, 220.0          # TACOT virgin and char solid densities (INPUT/pyrolysis.dat)
T_SLACK = 2.0                        # K, P-D4
FRAC_SLACK = 0.10                    # of PATO's peak, P-D4


def main():
    argv = sys.argv[1:]
    record = "--record" in argv
    args = [a for a in argv if a != "--record"]
    case = os.path.abspath(args[0]) if args else HERE
    ref = os.path.join(HERE, "ref")
    pato = os.environ.get("FUSS_PATO_ATC21_OUT") or ref     # a PATO output/porousMat directory overrides ref/
    exp_path = os.path.join(HERE, "expected.json")
    expected = json.load(open(exp_path)) if os.path.isfile(exp_path) else {}
    selftest = os.environ.get("FUSS_GATE_SELFTEST", "0") == "1"

    # ---- PATO and Amaryllis
    sub = (lambda f: os.path.join(pato, "scalar", f)) if os.path.isdir(os.path.join(pato, "scalar")) else \
          (lambda f: os.path.join(pato, f))
    Pp = at.pato_rows(sub("Ta_plot"))
    Ps = at.pato_rows(sub("Ta_surfacePatch"))
    Pm = at.pato_rows(os.path.join(pato, "mass"))
    E = at.pato_rows(os.path.join(ref, "Amaryllis_Energy_TestCase_2.1.txt"))
    M = at.pato_rows(os.path.join(ref, "Amaryllis_Mass_TestCase_2.1.txt"))
    n = min(len(Pp), len(Ps), len(Pm), len(E), len(M))

    # ---- FUSS
    pj = json.load(open(os.path.join(case, "probes.json")))
    L = pj["L"]
    px = {p["name"]: p["x"] for p in pj["probes"]}
    fT = {}
    for d in DEPTHS_MM:
        lo, hi = pj["pairs"][str(d)]
        a = at.pato_rows(os.path.join(case, "OUTPUT", lo + ".txt"))
        b = at.pato_rows(os.path.join(case, "OUTPUT", hi + ".txt"))
        w = ((L - d * 1e-3) - px[lo]) / (px[hi] - px[lo])
        if not -1e-9 <= w <= 1 + 1e-9:
            raise SystemExit("FAIL: probes %s/%s do not bracket %g mm" % (lo, hi, d))
        fT[d] = [(r[0], (1 - w) * r[1] + w * s[1]) for r, s in zip(a, b)]
    seb = at.read_seb(case)
    faces = at.faces_306(case)
    key0 = sorted(faces)[0]
    s0 = [r for r in seb if (r["b"], r["i"], r["j"], r["k"], r["f"]) == key0]
    fTw = [(r["time"], r["T_w"]) for r in s0]
    fmg = [(r["time"], r["mdot_pg"]) for r in s0]
    env = at.read_env(os.path.join(case, "INPUT", faces[key0][3]))
    kinks = [row[0] for row in env[1:-1]]
    fr = sorted(at.fronts(f, L, RHO_V, RHO_C) for f in glob.glob(os.path.join(case, "OUTPUT", "field*.tec")))
    fv = [(t, v) for t, v, c in fr]
    fc = [(t, c) for t, v, c in fr]

    # ---- maxima over PATO's grid
    res, spread = {}, {}
    skipped = []

    def upd(dct, k, d, t, where=""):
        if d is not None and abs(d) > dct.get(k, (-1.0,))[0]:
            dct[k] = (abs(d), t, where)

    for i in range(n):
        t = Pp[i][0]
        if t < 0.05:
            continue
        for k, d in enumerate(DEPTHS_MM):
            v = at.interp(fT[d], t)
            if v is not None:
                upd(res, "T", v - Pp[i][1 + k], t, "%d mm" % d)
            if E[i][2 + k] > 0:
                upd(spread, "T", Pp[i][1 + k] - E[i][2 + k], t, "%d mm" % d)
        br = at.bracket(fTw, t)
        if br is not None and any(br[0] < tk < br[1] or tk == t for tk in kinks):
            skipped.append(t)
        elif br is not None:
            upd(res, "Tw", at.interp(fTw, t) - Ps[i][1], t)
            upd(res, "mdot_g", at.interp(fmg, t) - Pm[i][1], t)
        upd(spread, "Tw", Ps[i][1] - E[i][1], t)
        upd(spread, "mdot_g", Pm[i][1] - M[i][1], t)
        v98, c02 = at.interp(fv, t), at.interp(fc, t)
        upd(res, "v98", None if v98 is None else v98 - Pm[i][3], t)
        upd(res, "c02", None if c02 is None else c02 - Pm[i][4], t)
        upd(spread, "v98", Pm[i][3] - M[i][3], t)
        upd(spread, "c02", Pm[i][4] - M[i][4], t)
    peaks = {"mdot_g": max(r[1] for r in Pm[:n]), "v98": max(r[3] for r in Pm[:n]), "c02": max(r[4] for r in Pm[:n])}
    units = {"T": "K", "Tw": "K", "mdot_g": "kg/m2/s", "v98": "m", "c02": "m"}

    ok = True
    print("ATC 2.1 gate: FUSS BC 306, static mesh, vs PATO (%s), t = 0.1..120 s; %d T_w/mdot_g grid times skipped "
          "(FUSS samples straddle an env breakpoint: %s)" % (pato, len(skipped), " ".join("%g" % t for t in skipped)))
    measured, budgets = {}, {}
    for q in ("T", "Tw", "mdot_g", "v98", "c02"):
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
            q, m, units[q], t, where, b, "PASS" if good else "FAIL", tgt, "met" if m <= tgt else "NOT met"))

    # static: nothing consumed, nothing receded
    worst_static = max(max(abs(r["mdot_ch"]), abs(r["Bch"])) for r in seb)
    st = [l.split() for l in open(os.path.join(case, "OUTPUT", "bc306-state.dat")) if not l.startswith("#")]
    rec = max(max(abs(float(r[5])), abs(float(r[7]))) for r in st)
    good = worst_static == 0.0 and rec == 0.0
    ok &= good
    print("  (static) max |mdot_ch|, |B'ch| over seb.txt %g; s_total, sdot in bc306-state.dat %g  %s" % (
        worst_static, rec, "PASS" if good else "FAIL"))

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
            "recorded": datetime.date.today().isoformat(),
            "note": "budget: T and Tw measured + 2 K, mdot_g and the fronts measured + 10 % of PATO's peak (plan 08 P-D4); "
                    "the 1.5x target is recorded, not gated; oracle ref/ = PATO 3.1 ATC 2.x with movingMesh no, DIC "
                    "preconditioners and maxDeltaT 0.0005 (ref/README.md), run 2026-10-13",
        }
        with open(exp_path, "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
            fh.write("\n")
        print("  recorded %s" % exp_path)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
