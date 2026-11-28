#!/usr/bin/env python3
"""Plan 08 gate 3.1 (b): with lambda = 0 the blowing correction is off, beta = 1 EXACTLY.

The case is the ATC 2.1 setup (atc2x306 of test/pyrolysis/tools/make_case.py, static mesh, no-char) with lambda = 0 on
the 306 line, 2 s, the chemistry on throughout (the rhoeUeCH0 ramp 0.003 -> 0.3 over the first 0.1 s included).
PATO's constantLambdaBlowingCorrectionModel returns 1 when lambda = 0, and so must FUSS on every face and output:
  (beta) beta == 1.0 bitwise in every seb.txt row
  (conv) the convective term equals rhoeUeCH0 (h_r - h_w), recomputed from the same row, to 1e-15 relative
  (seb)  plan 08 section 3.2 on every row (abl_tools.seb_closure; (b) has no row to test with lambda = 0)
Usage: verify.py [case_dir]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets the tolerances below zero.
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
    ok = True
    lams = sorted(set(v[0] for v in faces.values()))
    if lams != [0.0]:
        print("FAIL: the 306 faces carry lambda %s, the gate needs 0" % lams)
        return 1
    nb = sum(1 for r in rows if r["beta"] != 1.0)
    good = nb == 0 and not selftest
    ok &= good
    print("(beta) rows with beta != 1.0 (bitwise): %d of %d  %s" % (nb, len(rows), "PASS" if good else "FAIL"))
    worst = 0.0
    for r in rows:
        x = max(r["rhoeUeCH0"], 1.0e-6) * (r["h_r"] - r["h_w"])
        worst = max(worst, abs(r["conv"] - x) / max(abs(x), 1e-300))
    lim = -1.0 if selftest else 1.0e-15
    good = worst <= lim
    ok &= good
    print("(conv) max |conv - rhoeUeCH0 (h_r - h_w)| / |.| = %.3e  (<= %g)  %s" % (worst, lim, "PASS" if good else "FAIL"))
    good, lines, _ = at.seb_closure(case, 1.0e-30 if selftest else 1.0e-8)
    ok &= good
    print("\n".join(lines))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
