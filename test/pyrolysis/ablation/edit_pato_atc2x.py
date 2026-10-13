#!/usr/bin/env python3
"""edit_pato_atc2x.py <copy-of-AblationTestCase_2.x> --case {2.1,2.2,2.3} [--max-dt DT] [--ta-tol TOL] [--ncells N]
                   [--end-time T]

The PATO 3.1 oracle runs of plan 08 section 3.3, made from a COPY of $PATO_DIR/tutorials/1D/AblationTestCase_2.x
(an in-place ./Allrun writes 0/, constant/porousMat/polyMesh and output/ into PATO's tree):

  2.1  constant/porousMat/porousMatProperties: `movingMesh yes;` -> `movingMesh no;` (no recession), and in
       system/porousMat/fvSolution the preconditioner of Ta, p and Xsii DILU -> DIC: a static mesh makes those
       matrices symmetric and OpenFOAM aborts at step 1 on DILU ("Unknown symmetric matrix preconditioner DILU");
       same tolerances, so the converged answer is unchanged
  2.2  as shipped
  2.3  constant/porousMat/BoundaryConditions: h_r 1.5e6 -> 2.5e7 on the two heating rows (t = 0.1 and 60 s)
  --max-dt DT  system/controlDict maxDeltaT (the tutorial's 0.005 s). The ATC 2.1 oracle uses 0.0005: with the shipped
       step PATO's own time error is ~24 K at the wall 0.2 s after the 60 s switch-off and ~19 K at t = 0.2 s
       (ATC21/README.md, measured 2026-10-13)
  --ta-tol TOL  system/porousMat/fvSolution: the absolute tolerance of the Ta solve (shipped 1e-06, relTol 0). Below
       maxDeltaT ~0.0005 the per-step change of T falls under 1e-06 in the slow interior and the solve does zero
       iterations, i.e. T is not advanced: the 0.00025 ATC 2.1 run drifted 38 K between 40 and 60 s (30035 such solves).
       The reference step study uses 1e-12 there; count `No Iterations 0` in the log of any refined run.
  --ncells N  system/porousMat/blockMeshDict: N cells across the slab instead of 100, the grading (simpleGrading
       1 0.1 1, last/first cell 0.1) unchanged -- make_case.py --ncells N builds FUSS's mesh with the same rule
  --end-time T  system/controlDict endTime (120 s shipped): the mesh studies of the first seconds

Idempotent: an edit already made is left alone. Python 3.6, standard library only. Then:
  set +u; source /opt/anaconda3/etc/profile.d/conda.sh && conda activate pato
  cd <copy> && ./Allrun > log.Allrun 2>&1     # judge by the log (no FOAM FATAL) and output/porousMat/mass reaching 120
"""
import argparse
import os
import re
import sys


def edit(path, fn):
    with open(path) as fh:
        s = fh.read()
    t = fn(s)
    if t != s:
        with open(path, "w") as fh:
            fh.write(t)
        return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("copy")
    ap.add_argument("--case", required=True, choices=("2.1", "2.2", "2.3"))
    ap.add_argument("--max-dt", type=float, default=None)
    ap.add_argument("--ta-tol", type=float, default=None)
    ap.add_argument("--ncells", type=int, default=None)
    ap.add_argument("--end-time", type=float, default=None)
    a = ap.parse_args()
    c = a.copy
    for f in ("constant/porousMat/porousMatProperties", "system/porousMat/fvSolution",
              "constant/porousMat/BoundaryConditions", "system/controlDict", "system/porousMat/blockMeshDict"):
        if not os.path.isfile(os.path.join(c, f)):
            sys.exit("not an AblationTestCase_2.x copy: %s missing" % f)
    done = []
    if a.case == "2.1":
        if edit(os.path.join(c, "constant/porousMat/porousMatProperties"),
                lambda s: re.sub(r"^movingMesh(\s*)yes;", r"movingMesh\1no;", s, flags=re.M)):
            done.append("movingMesh no")

        def dic(s):
            n = s.count("preconditioner   DILU;")
            if n not in (0, 3):
                sys.exit("fvSolution: %d DILU preconditioners, expected the shipped 3 (Ta, p, Xsii)" % n)
            return s.replace("preconditioner   DILU;", "preconditioner   DIC;")
        if edit(os.path.join(c, "system/porousMat/fvSolution"), dic):
            done.append("DILU -> DIC for Ta, p, Xsii")
    if a.case == "2.3":
        def hr(s):
            t = re.sub(r"^(0\.1|60)(\s+101325\s+0\.3\s+)1\.5e6", lambda m: m.group(1) + m.group(2) + "2.5e7", s, flags=re.M)
            if t.count("2.5e7") != 2:
                sys.exit("BoundaryConditions: expected the two heating rows with h_r 1.5e6 (or already 2.5e7)")
            return t
        if edit(os.path.join(c, "constant/porousMat/BoundaryConditions"), hr):
            done.append("h_r 2.5e7 on the heating rows")
    if a.max_dt is not None:
        if edit(os.path.join(c, "system/controlDict"),
                lambda s: re.sub(r"^maxDeltaT(\s+)[^;]*;", lambda m: "maxDeltaT%s%s;" % (m.group(1), repr(a.max_dt)), s, flags=re.M)):
            done.append("maxDeltaT %s" % repr(a.max_dt))
    if a.ta_tol is not None:
        def tatol(s):
            i = s.index("Ta")
            j = s.index("}", i)
            blk = re.sub(r"tolerance(\s+)[^;]*;", lambda m: "tolerance%s%s;" % (m.group(1), repr(a.ta_tol)), s[i:j], count=1)
            return s[:i] + blk + s[j:]
        if edit(os.path.join(c, "system/porousMat/fvSolution"), tatol):
            done.append("Ta tolerance %s" % repr(a.ta_tol))
    if a.ncells is not None:
        def cells(s):
            t = re.sub(r"porousMat \(1 \d+ 1\) simpleGrading \(1 0\.1 1\)", "porousMat (1 %d 1) simpleGrading (1 0.1 1)" % a.ncells, s)
            if "porousMat (1 %d 1) simpleGrading (1 0.1 1)" % a.ncells not in t:
                sys.exit("blockMeshDict: the shipped hex line `porousMat (1 100 1) simpleGrading (1 0.1 1)` was not found")
            return t
        if edit(os.path.join(c, "system/porousMat/blockMeshDict"), cells):
            done.append("%d cells" % a.ncells)
    if a.end_time is not None:
        if edit(os.path.join(c, "system/controlDict"),
                lambda s: re.sub(r"^endTime(\s+)[^;]*;", lambda m: "endTime%s%s;" % (m.group(1), repr(a.end_time)), s, flags=re.M)):
            done.append("endTime %s" % repr(a.end_time))
    print("edit_pato_atc2x.py %s: %s" % (a.case, ", ".join(done) if done else "nothing to change"))


if __name__ == "__main__":
    main()
