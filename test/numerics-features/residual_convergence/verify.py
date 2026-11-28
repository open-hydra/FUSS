#!/usr/bin/env python3
"""
Residual convergence: does multigrid / IRS actually accelerate anything?

WHAT CHANGED AND WHY (plan 10, D7)
----------------------------------
These six cases used to measure ITERATIONS TO REACH res-threshold = 1e-8. That
is the natural statement, and it cost 1.21 million iterations across the family
-- `nominal` alone needed 369,089 -- which made the regression suite something
nobody would run casually.

Tightening the threshold does not help: convergence is geometric, so iterations
scale as log(1/tol) and going 1e-8 -> 1e-5 removes only about 37% of the work.

So the measurement is inverted. Every case now runs a FIXED budget of 30,000
fine-grid iterations with res-threshold set unreachably low, and what is
compared is the residual each one reaches. Same scientific claim, ~6x cheaper,
and deterministic by construction: no case can stop early, so no case can stop
at a different point than it did last time.

WHAT IS ASSERTED, AND WHAT IS NOT
---------------------------------
Asserted: the ORDERING. At equal fine-grid budget,

    nominal  >  mg1  >  mg2  >  mg3          more coarse pre-work, lower residual
    irs      <  nominal                      IRS accelerates
    mg_irs   <  everything                   both together

Not asserted: the residual VALUES. Those are already pinned bit-for-bit by the
regression manifest; repeating them here as thresholds would just be a second,
more brittle copy of the same baseline.

A caveat worth stating rather than hiding: the cases do not do equal TOTAL
work, because the multigrid ones additionally spend level2-iter iterations on
the coarse grid (2500 / 5000 / 7500 / 10000). That is deliberate -- a coarse
iteration is much cheaper than a fine one, and buying fine-grid progress with
coarse-grid work is precisely what multigrid is for. The comparison is at equal
FINE-grid budget.

PROVEN RED
----------
Setting mg3's level2-iter to 1 -- i.e. taking away its coarse-grid work --
fires two of the four checks: the budget check (it ran 30001 iterations against
an expected 37500) and the ordering check (mg2 5.9e-02 no longer exceeds mg3
1.8e-01). Checks 3 and 4 share their shape with check 2 and were not falsified
separately.

One attempt that did NOT work, recorded so nobody repeats it: setting
res-threshold to 1e-1 on `nominal` to make it stop early. `nominal` finishes the
budget at 1.8e-01, so the threshold was never reachable and nothing fired. If
you want to falsify the budget check directly, use a threshold above the
residual the case actually reaches.

WHY THIS FILE EXISTS AT ALL
---------------------------
Until 2026-09-10 the family's actual claim was asserted nowhere. It was read
off residual-history.dat by a human, if at all. That mattered: an uninitialised
`endsim` had been ending all four multigrid runs early, so they never reached
the fine grid and never produced a fine-grid solution -- and every artefact they
emitted was still perfectly reproducible, so the regression manifest was
entirely happy. A gate on the ordering would have caught it immediately.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

BUDGET = 30000          # fine-grid iterations, must match every input.ini
COARSE = {"nominal": 0, "irs": 0,
          "mg1": 2500, "mg2": 5000, "mg3": 7500, "mg_irs": 10000}
CASES = ["nominal", "irs", "mg1", "mg2", "mg3", "mg_irs"]


def final_residual(case):
    path = os.path.join(HERE, case, "OUTPUT", "residual-history.dat")
    with open(path) as fh:
        rows = [l.split() for l in fh if l.strip()]
    if not rows:
        raise RuntimeError("%s: empty residual history" % case)
    return len(rows), float(rows[-1][2])


def main():
    print()
    print("=" * 70)
    print(" RESIDUAL CONVERGENCE AT A FIXED %d-ITERATION FINE-GRID BUDGET" % BUDGET)
    print("=" * 70)
    print("  %-9s %-7s %-9s %s" % ("case", "coarse", "iters", "final residual"))

    res, iters = {}, {}
    for c in CASES:
        try:
            n, r = final_residual(c)
        except FileNotFoundError:
            print("  ERROR: %s has not been run" % c)
            return 1
        res[c], iters[c] = r, n
        print("  %-9s %-7d %-9d %.6e" % (c, COARSE[c], n, r))
    print()

    ok = True

    # 1. Every case must have run its whole budget. If one stopped early the
    #    comparison is between different amounts of work, and -- more to the
    #    point -- something is wrong with the stopping logic.
    for c in CASES:
        expect = BUDGET + COARSE[c]
        if abs(iters[c] - expect) > 1:
            print("  FAIL %s ran %d iterations, expected %d."
                  % (c, iters[c], expect))
            print("       A case that stops early has either hit res-threshold")
            print("       (it should be unreachable) or ended on an uninitialised")
            print("       flag -- see plan 10 D8.")
            ok = False
    if ok:
        print("  PASS every case ran its full budget")

    # 2. More coarse-grid pre-work leaves a lower residual.
    chain = ["nominal", "mg1", "mg2", "mg3"]
    bad = [(a, b) for a, b in zip(chain, chain[1:]) if not res[a] > res[b]]
    if bad:
        for a, b in bad:
            print("  FAIL %s (%.3e) should exceed %s (%.3e)"
                  % (a, res[a], b, res[b]))
        print("       Multigrid is not buying anything. Check that the coarse")
        print("       level is actually being run and prolongated.")
        ok = False
    else:
        print("  PASS nominal > mg1 > mg2 > mg3  (more coarse work, lower residual)")

    # 3. IRS accelerates.
    if not res["irs"] < res["nominal"]:
        print("  FAIL irs (%.3e) should be below nominal (%.3e)"
              % (res["irs"], res["nominal"]))
        ok = False
    else:
        print("  PASS irs < nominal  (%.1fx lower residual)"
              % (res["nominal"] / res["irs"]))

    # 4. Both together win.
    if not all(res["mg_irs"] < res[c] for c in CASES if c != "mg_irs"):
        print("  FAIL mg_irs (%.3e) is not the lowest residual" % res["mg_irs"])
        ok = False
    else:
        print("  PASS mg_irs lowest  (%.1e vs %.1e for nominal)"
              % (res["mg_irs"], res["nominal"]))

    print()
    print("  Result: %s" % ("PASS" if ok else "FAIL"))
    print()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
