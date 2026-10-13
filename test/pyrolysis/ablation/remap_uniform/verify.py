#!/usr/bin/env python3
"""
Gate V2 (plan 08 section 2.9, stage C): a uniform type-2 slab receding at a forced rate stays uniform and loses exactly
what the receding face sweeps out.

    python3 verify.py RUN_DIR [--qs]      # --qs: the run had [FUSS-Darcy] quasi-steady = true

The case (generate.py) holds T0 = 300 K, below every TACOT reaction threshold, and nothing drives a flux; the face
recedes at the TEST ONLY sdot-override SDOT for 50 steps. Checks, from OUTPUT/field.tec (the last step),
energy-budget.txt (every row) and bc306-state.dat:

  (uniform)   T, rho_s, eta_prev, eta_st, T_prev spread <= 1e-13 of max |v|; T = T0 to 1e-13; p_pg spread <= 1e-12:
              the Darcy line solve on a moved mesh, whose exact answer is p_w, returns it only to its own round-off
              (1.5e-13 measured, both modes; a Laplacian line of N = 100 cells with one Dirichlet end has condition
              ~N^2, i.e. 1e4 eps ~ 1e-12). The quasi-steady run shows the floor is the solve's own: its solve never
              reads the remapped p (no storage term, the answer is p_w by construction), and its spread is 5.9e-14.
              The "h" column is the type-2 E, which holds G h_g with G = eta p: it inherits
              p's round-off times |G h_g| / E, about 4.5 at 300 K, so 5e-12 (6.8e-13 measured). chi <= 1e-14: they start
              at 0 and no reaction runs, but every remap returns F (1 - chi) to a few ulp, the clamp at 0 keeps the
              positive part, and that walk reaches 2e-15 in 50 steps.
  (closed)    E_tot, M_s, M_g = X0 (L0 - s(t)) / L0 with s = SDOT t, to 1e-12 at every row (the plan's bound).
  (closure)   E_tot + E_abl = E0 and M_s + M_abl = M_s0 to 1e-12; the type-2 mass M_s + M_g + M_gas_out + M_abl + G_abl
              = M_s0 + M_g0 to 1e-13 of that total (M_gas_out carries the solve's flux round-off, ~1e-15 kg against a
              stored gas of 7e-4 kg, so the gas alone is judged by (closed)). Quasi-steady holds no stored gas but the
              carry (the M_g column is then the equation-of-state gas, not a conserved quantity, as in plan 07's gate
              4.1 b1): M_s + M_gas_out + M_abl + G_abl = M_s0 to 1e-13, and the swept carry G_abl, zero here but for the
              ulp by which eta_st and eta_prev can differ, below 1e-13 of the swept equation-of-state gas.
  (recession) s_total of bc306-state.dat = SDOT t_end to 1e-13; the receding face's nodes at x = L0 - s_total to 1e-15 m.

The measured floors: T recovery in type 2 stops at |E(T) - E| <= 1e-13 |E| (NEWTON_REL_F), which puts E's closure at
2e-13 without quasi-steady; quasi-steady E (H_s alone) is piecewise linear in T and closes to 3e-15.

FUSS_GATE_SELFTEST=1 sets every tolerance to 0 (the proven-red hook).
"""
import os
import re
import sys

import numpy as np

T0 = 300.0
SDOT = 1.0e-5
L0 = 0.05
TOL_UNI, TOL_CHI, TOL_P, TOL_H = 1.0e-13, 1.0e-14, 1.0e-12, 5.0e-12
TOL_CLOSED, TOL_CLOSURE, TOL_MASS = 1.0e-12, 1.0e-12, 1.0e-13
TOL_S, TOL_X = 1.0e-13, 1.0e-15
if os.environ.get("FUSS_GATE_SELFTEST") == "1":
    TOL_UNI = TOL_CHI = TOL_P = TOL_H = TOL_CLOSED = TOL_CLOSURE = TOL_MASS = TOL_S = TOL_X = 0.0


def read_field(path):
    lines = open(path).read().split("\n")
    names = lines[0].split("=", 1)[1].split()
    I, J, K = (int(v) for v in re.search(r"I=(\d+), J=(\d+), K=(\d+)", lines[1]).groups())
    vals = np.array(" ".join(lines[2:]).split(), float)
    nn, nc = I * J * K, (I - 1) * (J - 1) * max(K - 1, 1)
    out, off = {}, 0
    for n, nm in enumerate(names):
        size = nn if n < 3 else nc
        out[nm] = vals[off:off + size]
        off += size
    return out


def main():
    run = sys.argv[1]
    qs = "--qs" in sys.argv[2:]
    bad = []

    def check(tag, ok, text):
        print("%-11s %s  %s" % ("(" + tag + ")", "ok  " if ok else "FAIL", text))
        if not ok:
            bad.append(tag)

    f = read_field(run + "/OUTPUT/field.tec")
    for nm in ("T", "rho_s", "eta_prev", "eta_st", "T_prev", "h"):
        v = f[nm]
        spread = (v.max() - v.min()) / abs(v).max()
        tol = TOL_H if nm == "h" else TOL_UNI
        check("uniform", spread <= tol, "%-8s spread %.2e (<= %.0e)" % (nm, spread, tol))
    dT = np.abs(f["T"] - T0).max() / T0
    check("uniform", dT <= TOL_UNI, "T - T0   max %.2e of T0 (<= %.0e)" % (dT, TOL_UNI))
    chis = [nm for nm in f if nm.startswith("chi_")]
    chimax = max(np.abs(f[nm]).max() for nm in chis) if chis else float("nan")
    check("uniform", len(chis) == 3 and chimax <= TOL_CHI, "chi      %d columns, max |chi| %.2e (<= %.0e)" % (len(chis), chimax, TOL_CHI))
    v = f["p_pg"]
    spread = (v.max() - v.min()) / abs(v).max()
    check("uniform", spread <= TOL_P, "p_pg     spread %.2e (<= %.0e)" % (spread, TOL_P))

    B = np.loadtxt(run + "/OUTPUT/energy-budget.txt", ndmin=2)
    if B.shape[1] != 11:
        check("closed", False, "energy-budget.txt has %d columns, the receding layout has 11" % B.shape[1])
        return 1
    t, E, Ms, Mout, Mg, Ea, Ma, Ga = B[:, 0], B[:, 1], B[:, 2], B[:, 5], B[:, 7], B[:, 8], B[:, 9], B[:, 10]
    frac = (L0 - SDOT * t) / L0
    for nm, x in (("E_tot", E), ("M_s", Ms), ("M_g", Mg)):
        err = np.abs(x - x[0] * frac).max() / abs(x[0])
        check("closed", err <= TOL_CLOSED, "%-5s vs X0 (L0 - s)/L0: max %.2e over %d rows (<= %.0e)" % (nm, err, len(t), TOL_CLOSED))
    for nm, x, a in (("E", E, Ea), ("M_s", Ms, Ma)):
        err = np.abs(x + a - x[0]).max() / abs(x[0])
        check("closure", err <= TOL_CLOSURE, "%-3s + %s_abl: max %.2e (<= %.0e)" % (nm, nm.split("_")[0], err, TOL_CLOSURE))
    if qs:
        err = np.abs(Ms + Mout + Ma + Ga - Ms[0]).max() / Ms[0]
        check("closure", err <= TOL_MASS, "M_s + M_gas_out + M_abl + G_abl: max %.2e of M_s0 (<= %.0e)" % (err, TOL_MASS))
        err = np.abs(Ga).max() / (Mg[0] * SDOT * t[-1] / L0)
        check("closure", err <= TOL_MASS, "swept carry G_abl: max %.2e of the swept eq.-of-state gas (<= %.0e)" % (err, TOL_MASS))
    else:
        tot0 = Ms[0] + Mg[0]
        err = np.abs(Ms + Mg + Mout + Ma + Ga - tot0).max() / tot0
        check("closure", err <= TOL_MASS, "M_s + M_g + M_gas_out + M_abl + G_abl: max %.2e of the total (<= %.0e)" % (err, TOL_MASS))
        print("(closure)   info  G_abl / swept stored gas = %.15f" % (Ga[-1] / (Mg[0] * SDOT * t[-1] / L0)))

    rows = [l.split() for l in open(run + "/OUTPUT/bc306-state.dat") if not l.startswith("#")]
    s_tot = np.array([float(r[5]) for r in rows])
    err = np.abs(s_tot - SDOT * t[-1]).max() / (SDOT * t[-1])
    check("recession", len(rows) > 0 and err <= TOL_S, "s_total %.16e vs SDOT t_end: %.2e (<= %.0e)" % (s_tot[0], err, TOL_S))
    x = f["x"]
    xs = x.max()
    err = abs(xs - (L0 - s_tot[0]))
    check("recession", err <= TOL_X, "surface x %.15e vs L0 - s_total: %.2e m (<= %.0e)" % (xs, err, TOL_X))

    print("V2 %s: %s" % ("quasi-steady" if qs else "stored gas", "PASS" if not bad else "FAIL (" + " ".join(sorted(set(bad))) + ")"))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
