#!/usr/bin/env python3
"""Rewrite every record of one face of one block of a FUSS bc.txt to another BC type.

    usage: set_face_bc.py <bc.txt> <block> <face> <type> [<line> ...]

The records of the face keep their cell indices; their own lines (one for 101, 201, 103, 301-307; the "n1 n2" line
and n1 + n2 donor lines for 102, 104; none for 0, 200, 300) are replaced by the <line> arguments, written after every
record of the face, in the order given. The other records are copied as they are. The file is rewritten in place.

Used by the gates that build twins of a committed case at run time (run_physics_gates.sh: pyrolysis_103):
    set_face_bc.py INPUT/bc.txt 1 2 103 "1 1 1 1 1 1 1 1 1"     a coupled face, standalone (no partner)
    set_face_bc.py INPUT/bc.txt 1 1 300                         a symmetry face
    set_face_bc.py INPUT/bc.txt 1 2 104 "1 1" "1 1 1 1 1.0" "1 1 1 1 1.0"   a chimera coupled face
Python 3.6, no numpy.
"""
import sys

ONE_LINE = {101, 201, 103} | set(range(301, 308))
NO_LINE = {0, 200, 300}
CHIMERA = {102, 104}


def records(lines):
    """[(header, [own lines])] in file order; refuses a type it does not know"""
    out, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        head = lines[i]
        t = int(head.split()[5])
        i += 1
        own = []
        if t in ONE_LINE:
            own = [lines[i]]
            i += 1
        elif t in CHIMERA:
            n1, n2 = (int(v) for v in lines[i].split()[:2])
            own = lines[i:i + 1 + n1 + n2]
            i += 1 + n1 + n2
        elif t not in NO_LINE:
            sys.exit("set_face_bc.py: unknown BC type %d in %s" % (t, head.strip()))
        out.append((head, own))
    return out


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    path, blk, face, code = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    own_new = [l + "\n" for l in sys.argv[5:]]
    if code in ONE_LINE and len(own_new) != 1:
        sys.exit("set_face_bc.py: type %d takes one line" % code)
    if code in NO_LINE and own_new:
        sys.exit("set_face_bc.py: type %d takes no line" % code)
    lines = open(path).read().split("\n")
    lines = [l + "\n" for l in lines]
    n = 0
    text = []
    for head, own in records(lines):
        f = head.split()
        if int(f[0]) == blk and int(f[4]) == face:
            text.append("%8s%8s%8s%8s%8s%8d\n" % (f[0], f[1], f[2], f[3], f[4], code))
            text.extend(own_new)
            n += 1
        else:
            text.append(head)
            text.extend(own)
    if n == 0:
        sys.exit("set_face_bc.py: no record on face %d of block %d" % (face, blk))
    open(path, "w").write("".join(text))
    print("set_face_bc.py: %d record(s) of block %d face %d -> %d" % (n, blk, face, code))


if __name__ == "__main__":
    main()
