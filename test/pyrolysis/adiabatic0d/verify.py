#!/usr/bin/env python3
"""Gate 4.6 (verification record 2026-10-08): adiabatic 0-D decomposition against the exact type-1 ODE.

A 4-cell slab of real TACOT with every face adiabatic (BC 301, q = 0 on faces 1 and 2, 300 on 3/4) started at
T0 = 1000 K: the state stays uniform and follows the ODE of tools/ode0d.py (type1 variant) exactly; the only
FUSS error is the kinetics/energy splitting inside the RK3 stages (frozen stage temperature), second order in
dt here (0.051 K at vnn 0.5, 0.014 at 0.25, 0.0022 at 0.1, 0.0005 at 0.05; record 2026-10-08; the vnn values
before plan 14 S7, X/2 now).

Check (tools/ode0d.check_ladder, which states the rule in full): a ladder of runs of the case at vnn, vnn/2,
vnn/4 (the gate stages the vnn/2 and vnn/4 copies), coarsest first, against the event-located oracle:
  (a) max |T - T_or1| and (b) max |chi_r - chi_or1| fall on every rung and their order on the finest pair lies in
      [1.6, 2.4] (p = 2, tools/conv_rule.py)
  (c) the spatial non-uniformity of T and chi below the round-off bound of the run (N_steps x 3 x eps x max T)
  Reported, not gated: the exact type-2 (gas storage) minus type-1 difference (0.07 K here).

Under gas-closure = type2 (plan 07 gate 4.1(d), test/pyrolysis/darcy/adiabatic0d_type2: both 301 faces open
at p-wall through open-bc-types, K = 1e-7 m2, so that the pressure excursion stays at most 0.075 Pa and the
fixed-p oracle applies) the oracle is the type2 variant and (a), (b) take the same rule with p = 2, on the ladder
vnn/4, vnn/8, vnn/16: on the first two halvings from the committed vnn the T error does not yet fall (outside the
asymptotic range). At this excursion the case tests the kinetics/energy splitting under type 2, not the time order
of the pressure step. The non-uniformity bound is the floor 1e-5 of tools/ode0d.py, set from measured draws (the gas
storage -eps p inside E makes T follow the sub-pascal pressure excursion; the candidate physics bound printed beside
it is not adopted, see check_ladder).

Usage: verify.py case_dir [case_dir ...]   (coarsest first; one dir only: FAIL, the rule needs a ladder)
Exit 0 on PASS, 1 on FAIL. Python 3.6, no numpy (imports ../tools/ode0d.py, conv_rule.py).
FUSS_GATE_SELFTEST=1 empties the order window and the bounds (the proven-red hook).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import ode0d  # noqa: E402

DT_ORACLE = 2.0e-4


def main():
    cases = sys.argv[1:] or [os.path.join(HERE, "T1000")]
    return 0 if ode0d.check_ladder("adiabatic0d", cases, DT_ORACLE, pato_tables=None) else 1


if __name__ == "__main__":
    sys.exit(main())
