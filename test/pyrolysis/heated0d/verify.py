#!/usr/bin/env python3
"""Gate 4.7 (verification record 2026-10-08): heated 0-D decomposition against the exact type-1 ODE.

The adiabatic0d slab with a uniform volumetric source q = 3.75e8 W/m3 (INPUT/st.tec) from 300 K, 3 s: the
state heats at ~650 K/s (the rate of ATC1's 1 mm probe in its first second) through the whole decomposition
and follows tools/ode0d.py (type1) exactly. Run at vnn 0.003125 (0.00625 before plan 14 S7; dt ~ 2.1 ms = ATC1's step): the splitting
error is 0.0016 K (second order: 15.8 K at vnn 0.5, 3.8 at 0.25, 0.59 at 0.1, 0.14 at 0.05, 0.034 at 0.025,
0.008 at 0.0125; record 2026-10-08; the vnn values before plan 14 S7, X/2 now).

Check (tools/ode0d.check_ladder, which states the rule in full): a ladder of runs of the case at vnn, vnn/2,
vnn/4 (the gates stage the vnn/2 and vnn/4 copies), coarsest first, against the event-located oracle (its
self-difference at dt 1e-4 s vs 5e-5 s is 1.7e-9 K, the RHS jumps at each of the ~1350 integer kelvins crossed):
  (a) max |T - T_or1| and (b) max |chi_r - chi_or1| fall on every rung and their order on the finest pair lies in
      [1.6, 2.4] (p = 2, tools/conv_rule.py)
  (c) the spatial non-uniformity of T and chi below the round-off bound of the run (N_steps x 3 x eps x max T)
  Reported, not gated: type2 - type1 (exact, -3.9 K at the end) and type1 with PATO's table hbar minus type1
  (the data difference, -0.15 K) when PATO's TACOT tables are on this machine.

Under gas-closure = type2 (plan 07 gate 4.1(d), test/pyrolysis/darcy/heated0d_type2: both 301 faces open at
p-wall, K = 1e-7 m2 so that the pressure excursion stays below 0.5 Pa and the fixed-p oracle applies) the
oracle is the type2 variant; (a) and (b) must fall on every rung, their observed order is printed (no expected
order is stated for the type-2 step), and the non-uniformity bound is the floor 1e-5 of tools/ode0d.py, set from
measured draws (the gas storage -eps p inside E follows the sub-pascal pressure excursion; the candidate physics
bound printed beside it is not adopted, see check_ladder).

Usage: verify.py case_dir [case_dir ...]   (coarsest first; one dir only: FAIL, the rule needs a ladder)
Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy (imports ../tools/ode0d.py, conv_rule.py).
FUSS_GATE_SELFTEST=1 empties the order window and the bounds (the proven-red hook).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import ode0d  # noqa: E402

DT_ORACLE = 1.0e-4


def main():
    cases = sys.argv[1:] or [HERE]
    for c in cases:
        q = ode0d.Case(c, pato_tables=None).q
        if q <= 0.0:
            print("heated0d: FAIL  %s/INPUT/st.tec missing or q <= 0 (q = %g)" % (c, q))
            return 1
    return 0 if ode0d.check_ladder("heated0d", cases, DT_ORACLE, pato_tables=ode0d.PATO_TACOT) else 1


if __name__ == "__main__":
    sys.exit(main())
