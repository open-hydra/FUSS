#!/usr/bin/env python3
"""Plan 08 gate 3.2: the surface energy balance of BC 306 closes, and its blowing correction sits at the Picard's
closed-form fixed point.

The case is the ATC 2.1 setup (atc2x306 of test/pyrolysis/tools/make_case.py, static mesh, no-char, lambda 0.5) with
the heating switched off at 1 s instead of 60 s (--heat-until 1.0), 2 s: seb.txt rows every 0.1 s cover the chemistry
on (0-1 s, the rhoeUeCH0 ramp included), the 1-1.1 s ramp (PATO's 0 < chemistryOn < 1: the chemistry-off balance) and
the chemistry off with re-radiation (1.1-2 s). On every row (abl_tools.seb_closure):
  (seb-a) |F| / max_k |term_k| <= 10 seb-tol (the Newton stops at seb-tol on this scaled residual)
  (seb-b) chemistry on, none of PATO's beta = 1 exceptions: |beta - zeta/(e^zeta - 1)| <= 10 seb-tol,
          zeta = 2 lambda (mdot_pg + mdot_ch)/rhoeUeCH0
  (seb-c) chemistry on: B'pg = mdot_pg/(rhoeUeCH0 beta), mdot_ch = B'ch rhoeUeCH0 beta to 1e-14 relative
plus (rows) the three regimes each have rows, and (off) in the chemistry-off rows the convective term is
hconv (T_edge - T_w) = 0 (hconv 0) and the gas-advection term is 0 (PATO's chemistryOn == 0 balance).
Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets the tolerances to 1e-30.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import abl_tools as at  # noqa: E402


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
    selftest = os.environ.get("FUSS_GATE_SELFTEST", "0") == "1"
    rows = at.read_seb(case)
    faces = at.faces_306(case)
    envs = {}
    count = {"on": 0, "ramp": 0, "off": 0}
    bad_off = 0
    for r in rows:
        ef = faces[(r["b"], r["i"], r["j"], r["k"], r["f"])][3]
        if ef not in envs:
            envs[ef] = at.read_env(os.path.join(case, "INPUT", ef))
        chem = at.env_at(envs[ef], r["time"])[4]
        reg = "on" if chem >= at.CHEM_ON else ("ramp" if chem > 0 else "off")
        count[reg] += 1
        if reg != "on" and (r["conv"] != 0.0 or r["adv_pg"] != 0.0 or r["adv_ch"] != 0.0):
            bad_off += 1
    ok = True
    good = min(count.values()) > 0
    ok &= good
    print("(rows) chemistry on %d, ramp %d, off %d rows  %s" % (count["on"], count["ramp"], count["off"], "PASS" if good else "FAIL"))
    good = bad_off == 0
    ok &= good
    print("(off)  ramp/off rows with a convective or advective term != 0: %d  %s" % (bad_off, "PASS" if good else "FAIL"))
    good, lines, _ = at.seb_closure(case, 1.0e-30 if selftest else 1.0e-8)
    ok &= good
    print("\n".join(lines))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
