"""conv_rule.py -- the convergence rule of plan 17 part B (P17-D1), for the verifies of the exact-solution cases.

A ladder of runs, coarsest first, each with its refinement measure r (cells, or 1/vnn) and its error e against the exact
solution: pass = every rung's error below the coarser one's, and the observed order on the finest pair,
log(e[-2] / e[-1]) / log(r[-1] / r[-2]), inside [BAND_LO p, BAND_HI p] for the scheme's expected order p (P17-D5).
No stored FUSS number enters. FUSS_GATE_SELFTEST=1 empties the band (the proven-red hook). Python 3.6, no numpy.
"""
import math
import os

BAND_LO, BAND_HI = 0.8, 1.2      # P17-D5, option U20


def selftest():
    return os.environ.get("FUSS_GATE_SELFTEST", "0") == "1"


def order(e1, e2, r1, r2):
    return math.log(e1 / e2) / math.log(float(r2) / r1) if e1 > 0 and e2 > 0 and r2 != r1 else float("nan")


def check(label, errors, refine, p, fmt="%.4e"):
    """(ok, text) for one quantity on one ladder: errors and refine coarsest first, p the expected order"""
    lo, hi = (1.0e9, -1.0e9) if selftest() else (BAND_LO * p, BAND_HI * p)
    falls = all(errors[i] > errors[i + 1] for i in range(len(errors) - 1))
    o = order(errors[-2], errors[-1], refine[-2], refine[-1])
    ok = falls and lo <= o <= hi
    text = "%s %s; falls on every rung %s; order on the finest pair %.3f in [%g, %g]" % (
        label, " -> ".join(fmt % e for e in errors), falls, o, lo, hi)
    return ok, text

