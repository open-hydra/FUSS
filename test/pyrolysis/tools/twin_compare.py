#!/usr/bin/env python3
"""Compare the OUTPUT/ trees of two runs file by file, byte for byte.

    usage: twin_compare.py <case_a> <case_b> [--only-a NAME]... [--only-b NAME]... [--exclude PREFIX]...

Every file present in both OUTPUT/ trees must be byte-identical; a file present in one tree only must be named by
--only-a / --only-b (a path relative to OUTPUT/), and every name given must be absent from the other tree. Files whose
path starts with an --exclude PREFIX are left out of both trees (counted and named in the summary). For a file
that differs, the largest absolute difference between the numbers of the two files (token by token, same layout) is
printed, for the diagnosis. Exit 0 when the trees agree as stated, 1 otherwise. Python 3.6, no numpy.
"""
import os
import sys


def tree(case):
    root = os.path.join(case, "OUTPUT")
    out = set()
    for dp, _dn, fn in os.walk(root):
        for f in fn:
            out.add(os.path.relpath(os.path.join(dp, f), root))
    return root, out


def numdiff(pa, pb):
    ta, tb = open(pa).read().split(), open(pb).read().split()
    if len(ta) != len(tb):
        return "token counts %d vs %d" % (len(ta), len(tb))
    worst, where = 0.0, None
    for i, (x, y) in enumerate(zip(ta, tb)):
        if x == y:
            continue
        try:
            d = abs(float(x.rstrip(",")) - float(y.rstrip(",")))
        except ValueError:
            return "token %d: '%s' vs '%s'" % (i, x, y)
        if d > worst:
            worst, where = d, i
    return "max |difference| %.3e (token %s)" % (worst, where)


def main():
    args = sys.argv[1:]
    only = {"a": set(), "b": set()}
    excl = []
    pos = []
    i = 0
    while i < len(args):
        if args[i] in ("--only-a", "--only-b"):
            only[args[i][-1]].add(args[i + 1])
            i += 2
        elif args[i] == "--exclude":
            excl.append(args[i + 1])
            i += 2
        else:
            pos.append(args[i])
            i += 1
    if len(pos) != 2:
        sys.exit(__doc__)
    ra, fa = tree(pos[0])
    rb, fb = tree(pos[1])
    left = {f for f in fa | fb if any(f.startswith(e) for e in excl)}
    fa, fb = fa - left, fb - left
    bad = 0
    same = 0
    for f in sorted(fa & fb):
        a, b = open(os.path.join(ra, f), "rb").read(), open(os.path.join(rb, f), "rb").read()
        if a == b:
            same += 1
        else:
            bad += 1
            print("DIFFERS   %s: %s" % (f, numdiff(os.path.join(ra, f), os.path.join(rb, f))))
    for side, files, other in (("a", fa - fb, pos[0]), ("b", fb - fa, pos[1])):
        for f in sorted(files):
            if f in only[side]:
                print("only-%s    %s (expected)" % (side, f))
            else:
                bad += 1
                print("only-%s    %s (NOT expected)" % (side, f))
    for side, files in (("a", fa - fb), ("b", fb - fa)):
        for f in sorted(only[side] - files):
            bad += 1
            print("missing   %s was expected in %s only" % (f, side))
    tail = ", %d excluded (%s*)" % (len(left), "*, ".join(excl)) if excl else ""
    print("twin_compare: %d identical file(s), %d problem(s)%s" % (same, bad, tail))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
