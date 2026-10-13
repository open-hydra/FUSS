#!/usr/bin/env python3
"""
make_zhang2010.py -- the decomposing-slab verification case of Zhang (2010) as a FUSS case.

Zhenyu Zhang, PhD thesis, Virginia Tech, 2010, sec. 2.4.1, Table 2.1, solved there by three codes (ABAQUS UMATHT, a 1-D
UEL element, Summers' MATLAB FEM): Henderson-form model (Summers 2010 eq. 12-14, h = Q + int_Tref^T C dT,
hg = int_Tref^T Cg dT, Tref = 273 K); 12.7 mm slab, rho 1700 -> 1255, Cp 1100, Cg 3000, k 0.30, A 5e28 1/s,
E 3.62e5 J/mol, n 4.6, Q = -2e5 J/kg; exposed face q_rad 25 kW/m2 - sigma Ts^4 (eps 1, no back-radiation term,
eq. 2.18) + 10 (300 - Ts); unexposed face adiabatic; T_ini 300 K; 2000 s.

FUSS: type 1 (gas leaves through the exposed face); virgin and char share h(T) = Q + Cp (T - 273), so h_g - h_bar is
Summers' (hg - h) term for term; the exposed face is ONE 305 line whose Tref reproduces the three terms of eq. 2.18
for every T (lit_case.tr_305, Tsur = 0). Probes: the exposed cell, the two cells around the mid-plane (even ncells),
the back cell; the surface temperatures come from wall*.tec.

Usage: make_zhang2010.py <case-dir> [--ncells 50] [--masterdir DIR]
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lit_case as lc   # noqa: E402

L = 0.0127
RHO_V, RHO_F, CP, CG, K = 1700.0, 1255.0, 1100.0, 3000.0, 0.30
Q, A, E, NORD, TREF = -2.0e5, 5.0e28, 3.62e5, 4.6, 273.0
QRAD, EPS, H, TINF, TINI = 25000.0, 1.0, 10.0, 300.0, 300.0
TEND = 2000.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--ncells", type=int, default=50)
    ap.add_argument("--masterdir", default=None)
    a = ap.parse_args()
    if a.ncells % 2:
        sys.exit("make_zhang2010.py: --ncells must be even (the mid-plane lies between two cell centres)")
    case, nx = os.path.abspath(a.case_dir), a.ncells
    inp = os.path.join(case, "INPUT")
    os.makedirs(inp, exist_ok=True)
    name = "ZHANG"
    hs = lambda T: Q + CP * (T - TREF)
    lc.write_solid(os.path.join(inp, "virgin.dat"), name, 250, 1200, lambda T: CP, hs, lambda T: K, 1.0, 1.0)
    lc.write_solid(os.path.join(inp, "char.dat"), name, 250, 1200, lambda T: CP, hs, lambda T: K, 1.0, 1.0)
    lc.write_gas(os.path.join(inp, "gas.dat"), name, 250, 1200, 0.03, lambda T: CG * (T - TREF), lambda T: 3.0e-5, [1.0e4, 1.0e6])
    F = (RHO_V - RHO_F) / RHO_V
    lc.write_pyro(os.path.join(inp, "pyrolysis.dat"), name, [lc.phase(2000.0, 0.85, [lc.reaction(F, A, E, NORD, 0, 0.0)])],
                  "0.1", "0.1", provenance="test/pyrolysis/tools/make_zhang2010.py (Zhang 2010 Table 2.1)")
    xs = [L * i / nx for i in range(nx + 1)]
    Tr = lc.tr_305(QRAD * EPS, H, EPS, TINF, Tsur=0.0)
    title = ["Zhang 2010 Table 2.1 decomposing slab (three-code verification), Henderson-form energy equation",
             "%d cells over 12.7 mm; exposed face 305 (h %g, eps %g, Tref %.6f = q_rad 25 kW/m2 + h (300 - T) - sigma T^4)," % (nx, H, EPS, Tr),
             "back face adiabatic; rho 1700 -> 1255, A 5e28, E 3.62e5, n 4.6; Q -2e5 J/kg at Tref 273 K; Cp 1100, Cg 3000, k 0.30"]
    sec = lc.ini_common(title, TEND, 0.5, 10.0, ["[FUSS-Pyrolysis]", "gas-closure = type1", "outflow-direction = -i",
                                                  "p-gas = 101325.0"], generator="test/pyrolysis/tools/make_zhang2010.py")
    mid = nx // 2
    probes = [("front", 1, 1.0), ("midlo", mid, 1.0), ("midhi", mid + 1, 1.0), ("back", nx, 1.0)]
    lc.write_case(case, name, xs, (305, lc.mdc.VAL3_FMT % (H, EPS, Tr)), (301, lc.mdc.VAL_FMT % 0.0), [TINI] * nx, sec,
                  probes=probes, masterdir=a.masterdir,
                  meta={"kind": "zhang2010", "L": L, "ncells": nx, "Tr": Tr, "rho_v": RHO_V, "rho_f": RHO_F, "t_end": TEND,
                        "area": 2 * lc.DY * lc.DZ})
    print("make_zhang2010.py: %s  (%d cells, Tref %.6f)" % (case, nx, Tr))


if __name__ == "__main__":
    main()
