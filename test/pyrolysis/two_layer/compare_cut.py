#!/usr/bin/env python3
"""The two-material slab as one block against the same slab cut at the material interface into two blocks (101).

    compare_cut.py <uncut case> <cut case> [--tol K] [--vars T,p_pg]

Pairs the field files of the two runs in time order and, per pair, compares the SOLUTIONTIME and every listed
cell variable cell by cell, the cells of both runs ordered by their centre x (the cut run's two zones joined).
Prints the largest |difference| of each variable and of the time; exit 0 when every one is <= --tol (default 0:
bytes of the printed values). Python 3.6, no numpy.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify  # noqa: E402  (read_field, fields)


def cells(path, names):
    """{name: [values in x order]}, the SOLUTIONTIME"""
    import re
    lines = open(path).read().split("\n")
    hdr = [t.strip('"') for t in re.split(r"[\s,]+", lines[0].split("=", 1)[-1]) if t.strip('"')]
    starts = [i for i, l in enumerate(lines) if l.strip().startswith("ZONE")]
    rows, time = [], None
    for zi, s in enumerate(starts):
        m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[s])
        ni, nj, nk = (int(v) for v in m.groups())
        time = float(re.search(r"SOLUTIONTIME=\s*([-+0-9.eEdD]+)", lines[s]).group(1).replace("D", "E"))
        end = starts[zi + 1] if zi + 1 < len(starts) else len(lines)
        vals = []
        for l in lines[s + 1:end]:
            for t in l.split():
                try:
                    vals.append(float(t))
                except ValueError:      # DATASETAUXDATA lines
                    pass
        nn, nc = ni * nj * nk, (ni - 1) * (nj - 1) * max(nk - 1, 1)
        d, off = {}, 0
        for nm in hdr:
            n = nn if nm in ("x", "y", "z") else nc
            d[nm] = vals[off:off + n]
            off += n
        x = d["x"]
        for c in range(nc):
            i = c % (ni - 1)
            r = c // (ni - 1)
            rows.append((r, 0.5 * (x[i] + x[i + 1])) + tuple(d[nm][c] for nm in names))
    rows.sort()
    return {nm: [row[2 + k] for row in rows] for k, nm in enumerate(names)}, time


def main():
    args = sys.argv[1:]
    tol, names, pos = 0.0, ["T"], []
    i = 0
    while i < len(args):
        if args[i] == "--tol":
            tol = float(args[i + 1]); i += 2
        elif args[i] == "--vars":
            names = args[i + 1].split(","); i += 2
        else:
            pos.append(args[i]); i += 1
    a, b = [sorted((t, f) for t, f, _z in verify.fields(c)) for c in pos]
    if len(a) != len(b) or not a:
        print("compare_cut: %d vs %d field files" % (len(a), len(b)))
        return 1
    worst = {nm: 0.0 for nm in names}
    wt = 0.0
    for (ta, fa), (tb, fb) in zip(a, b):
        va, _ = cells(fa, names)
        vb, _ = cells(fb, names)
        wt = max(wt, abs(ta - tb))
        for nm in names:
            if len(va[nm]) != len(vb[nm]):
                print("compare_cut: %s has %d cells, %s %d" % (fa, len(va[nm]), fb, len(vb[nm])))
                return 1
            worst[nm] = max(worst[nm], max(abs(x - y) for x, y in zip(va[nm], vb[nm])))
    ok = all(v <= tol for v in worst.values()) and wt <= tol
    print("compare_cut: %d field files; max |d time| %.3e s; %s (tol %g)  %s"
          % (len(a), wt, "; ".join("max |d %s| %.3e" % (nm, worst[nm]) for nm in names), tol, "ok" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
