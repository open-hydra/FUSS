#!/usr/bin/env python3
"""
Receding conjugate two-layer slab (plan 08 section 2.3 R5, stage B): the block interface stays conservative while one
block recedes.

Block 1 (k 50, 20 cells of 5 mm) recedes from its outer face at the TEST ONLY sdot-override of 1e-5 m/s (law nodeshift,
uniform compression toward its pinned face, the 101 interface); block 2 (k 5, 5 cells of 20 mm) does not move. Every
step after the move the ghost metrics are rebuilt (Setup_BC_Metrics): bc%Mg of block 2's interface faces must be the
CURRENT metric of block 1's compressed interface cells, or the two sides build different interface temperatures and
the interface leaks energy.

Checks:
  (pair)  at every solution output, for every 101 face, the flux it put into the residual and its partner's cancel:
          max |q_1 + q_2| / q_scale <= 1e-12, with q_scale = 2 k_1 A T_ref / dx_1 = 1000 W (k_1 50 W/m/K, A 5e-5 m2,
          T_ref 1000 K, dx_1 5 mm): the flux a temperature T_ref drives across half a cell, so the bound is a few
          hundred round-off units of the fluxes the face computes. (The plan's max|q| normalisation was dropped: early
          on the interface flux is ~3e-7 W, a difference of temperatures that agree to 7 digits, and even the correct
          scheme pairs only to 1e-7 of it; measured 2026-10-13.)
  (moved) the mesh moved: block 1's outer face is at x = s_total (bc306-state.dat) on every j (to the 15 digits of
          field.tec), its lines are compressed uniformly, the interface nodes (x = 0.1) and block 2 are unmoved.
The red control (run_physics_gates.sh) runs a copy with [FUSS-MeshMotion] skip-bc-metrics-refresh = true: (pair) must
fail (measured 8e-10 of q_scale at 12.5 s, against 4e-17 with the refresh).

Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets the bound to 0.
"""
import os
import re
import sys

Q_SCALE = 2.0 * 50.0 * 5.0e-5 * 1000.0 / 5.0e-3
L1 = 0.1


def read_zones(path):
    text = open(path).read().split("\n")
    names = text[0].split("=", 1)[1].split()
    zones, cur = [], None
    for line in text[1:]:
        if line.lstrip().startswith("ZONE"):
            m = re.search(r"I=(\d+), J=(\d+), K=(\d+)", line)
            cur = {"dim": tuple(int(x) for x in m.groups()), "v": []}
            zones.append(cur)
        elif cur is not None and line.strip():
            cur["v"] += [float(x) for x in line.split()]
    out = []
    for z in zones:
        I, J, K = z["dim"]
        nn = I * J * K
        x = z["v"][0:nn]
        out.append({"dim": z["dim"], "x": x})
    return out


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    tol = 0.0 if os.environ.get("FUSS_GATE_SELFTEST", "0") == "1" else 1.0e-12
    ok = True

    rows = [l.split() for l in open(os.path.join(case, "OUTPUT", "interface-flux.txt")) if not l.startswith("#")]
    by_t = {}
    for r in rows:
        by_t.setdefault(r[0], {})[tuple(int(x) for x in r[1:6])] = (tuple(int(x) for x in r[6:11]), float(r[11]))
    worst, worst_t, n_out = 0.0, None, 0
    for t, m in sorted(by_t.items(), key=lambda kv: float(kv[0])):
        n_out += 1
        for k, (partner, q) in m.items():
            if partner not in m:
                print("FAIL: face %s has no partner row %s at t = %s" % (k, partner, t))
                return 1
            mis = abs(q + m[partner][1]) / Q_SCALE
            if mis > worst:
                worst, worst_t = mis, float(t)
    good = n_out >= 4 and worst <= tol
    ok &= good
    print("(pair)  max |q_1 + q_2| / q_scale = %.3e at t = %s over %d outputs (<= %g, q_scale %g W)  %s"
          % (worst, worst_t, n_out, tol, Q_SCALE, "PASS" if good else "FAIL"))

    st = [l.split() for l in open(os.path.join(case, "OUTPUT", "bc306-state.dat")) if not l.startswith("#")]
    s_total = float(st[0][5])
    z = read_zones(os.path.join(case, "OUTPUT", "field.tec"))
    b1, b2 = z[0], z[1]
    I1, J1, K1 = b1["dim"]
    bad_surface = bad_interface = 0
    spread = 0.0
    for k in range(K1):
        for j in range(J1):
            line = [b1["x"][i + I1 * (j + J1 * k)] for i in range(I1)]
            if abs(line[0] - s_total) > 1.0e-17:     # field.tec prints 15 digits: 5e-19 m here
                bad_surface += 1
            if line[-1] != L1:
                bad_interface += 1
            dx = [line[i + 1] - line[i] for i in range(I1 - 1)]
            spread = max(spread, max(dx) - min(dx))
    I2, J2, K2 = b2["dim"]
    x2 = [b2["x"][i + I2 * (j + J2 * k)] for k in range(K2) for j in range(J2) for i in range(I2)]
    x2_expected = [L1 + 0.02 * i for k in range(K2) for j in range(J2) for i in range(I2)]
    b2_moved = max(abs(a - b) for a, b in zip(x2, x2_expected))
    good = s_total > 0 and bad_surface == 0 and bad_interface == 0 and spread < 1e-15 and b2_moved < 1e-15
    ok &= good
    print("(moved) s_total %.6e m; outer face at x = s_total on %d/%d node lines; interface at x = 0.1 on %d/%d; "
          "cell-size spread %.1e m; block 2 max |x - x0| %.1e  %s"
          % (s_total, J1 * K1 - bad_surface, J1 * K1, J1 * K1 - bad_interface, J1 * K1, spread, b2_moved,
             "PASS" if good else "FAIL"))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
