#!/usr/bin/env python3
"""
make_baum_atreya.py -- the Baum & Atreya (2007) sharp-front charring problem with Darcy gas flow, as a FUSS case.

H.R. Baum, A. Atreya, Proc. Combust. Inst. 31 (2007) 2633: a semi-infinite virgin solid at T0, its surface held at Ts from
t = 0, decomposing at Tp from rho_v to rho_c with the heat Q per kg of virgin; the gas flows through the char to the
surface at p_s (Darcy, no storage). The front sits at X = C* sqrt(alpha_c t); C* from the CORRECTED front condition (the
printed eq. 17 has lambda_c/lambda_v where rho_c c_c/(rho_v c_v) belongs; = Atreya & Baum 2002 eq. 27), computed by
darcy_exact.sharp_front_lambda (C* = 2 lam). Di Blasi property set; --cg 1475 = (rv cv - rc cc)/(rv - rc) makes Q
independent of T and puts the gas-enthalpy advection through the char into the exact solution (B = 5.36); --cg 0 is the
paper's own case (B = 0).

FUSS: type 2; --mode qs = quasi-steady (the paper's beta = 0), --mode full = with gas storage (eps_c 0.76, about 0.4 % of
the decomposed mass: the no-storage solution within that). The sharp front is zeroth-order kinetics from Tp with
A = sdot(100 s)/dx: one cell decomposes in one front-crossing time (instant kinetics release a whole cell's gas in one
step, 217 bar in the first cell), and A -> infinity as dx -> 0. Slab 40 mm, surface 302 open at p_s, back 301 sealed.
Probes (T, p_pg): the back cell (p there = p at the front: no flow beyond it), the surface cell, and the cells nearest
0.5, 1 and 2 mm (char pressure).

Usage: make_baum_atreya.py <case-dir> [--cg 1475] [--ncells 100] [--mode full|qs] [--ramp 0] [--tend 300] [--masterdir DIR]
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lit_case as lc        # noqa: E402
import darcy_exact as dx     # noqa: E402

RV, RC, CV, CC, KV, KC = 700.0, 140.0, 1400.0, 1100.0, 0.209, 0.071
Q, TP, T0, TS = 497302.6, 598.16, 298.16, 1298.16
W, MU, K_C, K_V, PS = 0.08, 2.3e-5, 1.0e-15, 1.0e-20, 101325.0
L, T_REF = 0.040, 100.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--cg", type=float, default=1475.0, help="gas cp [J/kg/K]: 1475 (Q flat, transpiration) or 0 (B&A)")
    ap.add_argument("--ncells", type=int, default=100)
    ap.add_argument("--mode", choices=["full", "qs"], default="full")
    ap.add_argument("--ramp", type=float, default=0.0, help="surface temperature ramp [s] (0 = the exact step)")
    ap.add_argument("--tend", type=float, default=300.0)
    ap.add_argument("--masterdir", default=None)
    a = ap.parse_args()
    case, nx, cg = os.path.abspath(a.case_dir), a.ncells, a.cg
    inp = os.path.join(case, "INPUT")
    os.makedirs(inp, exist_ok=True)
    name = "BA"
    xs = [L * i / nx for i in range(nx + 1)]
    av, ac = KV / (RV * CV), KC / (RC * CC)
    lam = dx.sharp_front_lambda(ac, av, KC, KV, RV, RC, CC, cg, TS, TP, T0, RV * Q)
    A_kin = lam * math.sqrt(ac / T_REF) / (L / nx)          # sdot(T_REF) / dx, sdot = (C*/2) sqrt(ac/t)
    lc.write_solid(os.path.join(inp, "virgin.dat"), name, 250, 1400, lambda T: CV, lambda T: CV * (T - TP) - Q, lambda T: KV)
    lc.write_solid(os.path.join(inp, "char.dat"), name, 250, 1400, lambda T: CC, lambda T: CC * (T - TP), lambda T: KC)
    lc.write_gas(os.path.join(inp, "gas.dat"), name, 250, 1400, W, lambda T: cg * (T - TP), lambda T: MU, [1.0e4, 1.0e8])
    lc.write_pyro(os.path.join(inp, "pyrolysis.dat"), name, [lc.phase(1000.0, 0.7, [lc.reaction(0.8, A_kin, 0.0, 0, 0, TP)])],
                  "0.3", "0.76", perm={"K_v": [K_V, 0, 0, 0, K_V, 0, 0, 0, K_V], "K_c": [K_C, 0, 0, 0, K_C, 0, 0, 0, K_C]},
                  provenance="test/pyrolysis/tools/make_baum_atreya.py (Baum & Atreya 2007, Di Blasi set)")
    qs = a.mode == "qs"
    title = ["Baum & Atreya 2007 sharp-front charring + Darcy (corrected front condition), Di Blasi set",
             "%d cells over %g mm, surface %.2f K (302, open at p_s%s), T0 %.2f K, zeroth-order decomposition from Tp %.2f K," % (
                 nx, L * 1e3, TS, ", %g s ramp" % a.ramp if a.ramp > 0 else "", T0, TP),
             "A %.6g /s (one cell per crossing time at %g s); c_g %g; K_c %g; %s; C* = %.6f" % (
                 A_kin, T_REF, cg, K_C, "quasi-steady" if qs else "gas storage (eps_c 0.76)", 2 * lam)]
    sec = lc.ini_common(title, a.tend, 0.5, 10.0, ["[FUSS-Pyrolysis]", "gas-closure = type2"],
                        ["[FUSS-Darcy]", "solver = tridiag", "tol = 1e-12", "max-iter = 100000", "p-init = %.1f" % PS,
                         "p-wall = %.1f" % PS, "open-bc-types = 302", "quasi-steady = %s" % ("true" if qs else "false")],
                        generator="test/pyrolysis/tools/make_baum_atreya.py")
    dxc = L / nx
    probes = [("back", nx, 0.05), ("surf", 1, 0.05)]
    probes += [("x%02d" % round(10 * d), int(round(d * 1e-3 / dxc + 0.5)), 0.05) for d in (0.5, 1.0, 2.0)]
    if a.ramp > 0:
        face1, tfiles = (302, lc.mc.TIMEFILE_LINE), {"Tw_time.dat": [(0.0, T0), (a.ramp, TS), (1.0e6, TS)]}
    else:
        face1, tfiles = (302, lc.mdc.VAL_FMT % TS), None
    lc.write_case(case, name, xs, face1, (301, lc.mdc.VAL_FMT % 0.0), [T0] * nx, sec, probes=probes, probe_vars="T p_pg",
                  masterdir=a.masterdir, time_files=tfiles,
                  meta={"kind": "baum_atreya", "rv": RV, "rc": RC, "cv": CV, "cc": CC, "kv": KV, "kc": KC, "Q": Q, "Tp": TP,
                        "T0": T0, "Ts": TS, "cg": cg, "W": W, "mu_p": MU, "nmu": 0, "Kc": K_C, "ps": PS, "L": L,
                        "ncells": nx, "t_end": a.tend, "A": A_kin, "t_ref": T_REF, "lam": lam, "Cstar": 2 * lam,
                        "alpha_v": av, "alpha_c": ac, "mode": a.mode, "ramp": a.ramp})
    print("make_baum_atreya.py: %s  (%d cells, c_g %g, %s, C* %.6f)" % (case, nx, cg, a.mode, 2 * lam))


if __name__ == "__main__":
    main()
