#!/usr/bin/env python3
"""
Plan 12, FUSS standalone: the per-face blocks of the TEST ONLY coupled-env-file (`// face b i j k f`, hydra's env-dump
replayed). Five copies of run307 that differ only in INPUT/env307.dat (generate.py, faces/):

    python3 verify_faces.py DEST      # DEST/run307 (one table), DEST/same, DEST/ab, DEST/ab_rev, DEST/ba

  same    one block per face, both equal to run307's table   -> every OUTPUT file equals run307's, byte for byte
  ab_rev  ab's blocks in the other order                       -> every OUTPUT file equals ab's (a block goes to the face
                                                                  its header names, not to a place in the file)
  ab, ab_rev, ba                                               -> each face's seb.txt rows carry its own block's
                                                                  rhoeUeCH0 and h_r exactly, and alpha_qrad > 0 only where
                                                                  its block's q_rad is
  ba      ab's blocks swapped between the two faces            -> the j-mirror of ab: every field file and both faces'
                                                                  seb.txt rows, |ab(i,j) - ba(i,3-j)| <= MIRROR_TOL x the
                                                                  field's max (round-off: the mesh is j-symmetric and the
                                                                  solver is not exactly), while ba itself differs from
                                                                  ab by >= SWAP_MIN (the swap is not a no-op)

Exit 0 on PASS, 1 on FAIL. Python 3.6.
"""
import filecmp
import os
import sys

MIRROR_TOL = 1e-8      # measured 2026-10-22: 4.1e-11 (seb.txt mdot_pg); ab against its own mirror: 7.1
SWAP_MIN = 1e-3        # measured: 0.27 (T)
SKIP = ("x", "y", "z", "matID")


def outputs(d):
    return os.path.join(d, "OUTPUT")


def identical(a, b):
    fa, fb = sorted(os.listdir(outputs(a))), sorted(os.listdir(outputs(b)))
    same = [f for f in fa if f in fb and filecmp.cmp(os.path.join(outputs(a), f), os.path.join(outputs(b), f),
                                                     shallow=False)]
    return len(same), max(len(fa), len(fb)), fa == fb and len(same) == len(fa)


def blocks(d):
    """{(b, i, j, k, f): rows} of the run's INPUT/env307.dat."""
    out, key = {}, None
    for line in open(os.path.join(d, "INPUT", "env307.dat")):
        t = line.split()
        if t[:2] == ["//", "face"]:
            key = tuple(int(v) for v in t[2:7])
            out[key] = []
        elif t and not t[0].startswith("//"):
            out[key].append([float(v) for v in t])
    return out


def seb(d):
    lines = open(os.path.join(outputs(d), "seb.txt")).read().split("\n")
    names = [n for n in lines[0].lstrip("#").split() if not n.startswith("(")]   # "(SI)" closes the header
    rows = {}
    for line in lines[1:]:
        if line.strip() and not line.startswith("#"):
            v = line.split()
            rows.setdefault(tuple(int(x) for x in v[1:6]), []).append([float(x) for x in v[:len(names)]])
    return names, rows


def face_values(d):
    """Problems of run d: a face whose seb.txt rows do not carry its block's values."""
    bad = []
    names, rows = seb(d)
    c_rho, c_hr, c_q = names.index("rhoeUeCH0"), names.index("h_r"), names.index("alpha_qrad")
    for key, blk in sorted(blocks(d).items()):
        vals = set(tuple(r[1:]) for r in blk if r[1] > 0.0)
        if len(vals) != 1:
            bad.append("%s %s: the test needs one value set per block" % (os.path.basename(d), key))
            continue
        p_w, rho, h_r, q_rad = vals.pop()
        r = rows.get(key, [])
        if not r:
            bad.append("%s %s: no seb.txt rows" % (os.path.basename(d), key))
        elif any(x[c_rho] != rho or x[c_hr] != h_r for x in r):
            bad.append("%s %s: rhoeUeCH0/h_r not its block's (%r, %r)" % (os.path.basename(d), key, rho, h_r))
        elif any((x[c_q] > 0.0) != (q_rad > 0.0) for x in r):
            bad.append("%s %s: alpha_qrad > 0 does not follow its block's q_rad = %r" % (os.path.basename(d), key, q_rad))
    return bad


def field(path):
    lines = open(path).read().split("\n")
    names = lines[0].split("=", 1)[1].split()
    zone = lines[1]
    n = [int(zone.split(c + "=")[1].split(",")[0]) for c in ("I", "J", "K")]
    vals = [float(x) for x in lines[2:] if x.strip()]
    nn, nc = n[0] * n[1] * n[2], (n[0] - 1) * (n[1] - 1) * max(n[2] - 1, 1)
    out, o = {}, 0
    for name in names:
        m = nn if name in ("x", "y", "z") else nc
        out[name] = vals[o:o + m]
        o += m
    return out, n[0] - 1, n[1] - 1


def mirror(a, b):
    """(max mirrored difference, max direct difference of T, where), relative to each quantity's max over ab."""
    worst, direct, where = 0.0, 0.0, ""
    files = sorted(f for f in os.listdir(outputs(a)) if f.startswith("field") and f.endswith(".tec"))
    for f in files:
        fa, ni, nj = field(os.path.join(outputs(a), f))
        fb, _, _ = field(os.path.join(outputs(b), f))
        for name in fa:
            if name in SKIP:
                continue
            scale = max(abs(v) for v in fa[name])
            if scale == 0.0:
                continue
            sign = -1.0 if name == "u_pg_2" else 1.0
            for j in range(nj):
                for i in range(ni):
                    va = fa[name][j * ni + i]
                    e = abs(va - sign * fb[name][(nj - 1 - j) * ni + i]) / scale
                    if e > worst:
                        worst, where = e, "%s %s" % (f, name)
                    if name == "T":
                        direct = max(direct, abs(va - fb[name][j * ni + i]) / scale)
    names, ra = seb(a)
    _, rb = seb(b)
    keys = sorted(ra)
    cols = [c for c in range(names.index("T_w"), len(names))
            if names[c] not in ("iter", "F") and not names[c].startswith("nclamp")]
    for ka, kb in ((keys[0], keys[1]), (keys[1], keys[0])):
        for c in cols:
            scale = max(abs(x[c]) for x in ra[ka])
            if scale == 0.0:
                continue
            for xa, xb in zip(ra[ka], rb[kb]):
                e = abs(xa[c] - xb[c]) / scale
                if e > worst:
                    worst, where = e, "seb.txt %s" % names[c]
    return worst, direct, where


def main():
    dest = os.path.abspath(sys.argv[1])
    run = {r: os.path.join(dest, r) for r in ("run307", "same", "ab", "ab_rev", "ba")}
    bad = []
    n1, m1, ok1 = identical(run["same"], run["run307"])
    if not ok1:
        bad.append("same: %d of %d files equal run307's" % (n1, m1))
    n2, m2, ok2 = identical(run["ab_rev"], run["ab"])
    if not ok2:
        bad.append("ab_rev: %d of %d files equal ab's" % (n2, m2))
    for r in ("ab", "ab_rev", "ba"):
        bad += face_values(run[r])
    worst, direct, where = mirror(run["ab"], run["ba"])
    if not worst <= MIRROR_TOL:
        bad.append("ba is not the j-mirror of ab: %.3e (%s) > %.0e" % (worst, where, MIRROR_TOL))
    if not direct >= SWAP_MIN:
        bad.append("ba differs from ab by %.3e only (< %.0e): the swap changed nothing" % (direct, SWAP_MIN))
    print("coupled307 faces: same %d/%d = run307; ab_rev %d/%d = ab; face values %s; mirror ba/ab %.1e (%s; <= %.0e), "
          "swap %.2f (T)%s" % (n1, m1, n2, m2, "exact" if not any("rhoeUeCH0" in b or "alpha" in b for b in bad) else
                               "WRONG", worst, where, MIRROR_TOL, direct, "" if not bad else " -- " + "; ".join(bad)))
    print("PASS" if not bad else "FAIL")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
