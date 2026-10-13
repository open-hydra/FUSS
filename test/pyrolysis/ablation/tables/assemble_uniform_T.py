#!/usr/bin/env python3
"""Assemble the B' table FUSS reads (plan 08 sections 2.2 and 3.4) from a createBprimeTable run.

createBprimeTable (PATO 3.1) writes, besides its own table, the uniform-T scans of every pressure
and B'g of its inputs: <run>/p_<p>/Bc_Bg<v> (rows `T B'c`) and <run>/p_<p>/hw_Bg<v> (rows `T h_w`),
nTemp rows each. This script writes them in PATO's 5-column layout, `p B'g B'c T h_w`, pressure
blocks in the order of `pressures`, B'g sub-blocks in the order of `BgList`, T ascending in each,
under the two header lines FUSS_Lib_Bprime reads (`// kind:` and `// axes:`).

Why not the generator's own table: it places its rows at log-uniform B'c, so across the oxidation
plateau one sub-block jumps 1093 -> 2187 K, and a PATO ATC 2.2 run driven by it misses the
on-the-fly run by 15.7 K at the wall; this assembly of the same run misses it by 0.16 K (README.md).

Usage: assemble_uniform_T.py <createBprimeTable run dir> <out file>
       (the run dir holds createBprimeTable.inputs and the p_<p>/ scan directories)
Python 3.6, no third-party modules.
"""
import re
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    run, out = sys.argv[1], sys.argv[2]
    inp = open(run + "/createBprimeTable.inputs").read()
    pressures = [float(x) for x in re.search(r"pressures\s*\(([^)]*)\)", inp).group(1).split()]
    bgs = [float(x) for x in re.search(r"BgList\s*\(([^)]*)\)", inp).group(1).split()]
    ntemp = int(re.search(r"nTemp\s+(\d+)", inp).group(1))

    def rd(f):
        return [[float(x) for x in l.split()] for l in open(f) if l.strip() and not l.startswith("//")]

    n = 0
    with open(out, "w") as o:
        o.write("// kind: equilibrium\n// axes: p Bprime-g T\n")
        for p in pressures:
            src = "%s/p_%g" % (run, p)
            for bg in bgs:
                bc = rd(src + "/Bc_Bg%g" % bg)
                hw = rd(src + "/hw_Bg%g" % bg)
                if not (len(bc) == len(hw) == ntemp):
                    sys.exit("%s: Bg %g has %d/%d rows, expected nTemp %d" % (src, bg, len(bc), len(hw), ntemp))
                for a, b in zip(bc, hw):
                    if a[0] != b[0]:
                        sys.exit("%s: Bg %g: T %g in Bc_Bg, %g in hw_Bg" % (src, bg, a[0], b[0]))
                    o.write("%g %g %g %g %g\n" % (p, bg, a[1], a[0], b[1]))
                    n += 1
    print("%s: %d rows (%d pressures x %d B'g x %d T)" % (out, n, len(pressures), len(bgs), ntemp))


if __name__ == "__main__":
    main()
