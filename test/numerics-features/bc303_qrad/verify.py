#!/usr/bin/env python3
"""Gate: the radiative flux read for a type-303 face must reach the flux.

Checks, on every cell of every 303 face (B1F2 and B1F4 here), the identity the BC
defines (Lib_BC_Fluxes_Wall_HeatTransfer.f90):

    qw = hconv*(Tref - Tw) + qrad      =>     qw - hconv*(Tref - Tw) - qrad = 0

with hconv, qrad, Tref taken from INPUT/bc.txt.  Exact to round-off, so the
tolerance is 1e-9 relative to qrad.  Before the dispatch fix the residual is -qrad
(the dispatch passed the never-assigned bc%qw, which happened to be zero).
Exit 0 on pass, 1 on fail.  Usage: verify.py [case_dir]  (default: this dir)
"""
import os, re, sys

CASE = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
WALL = os.path.join(CASE, "OUTPUT", "wall.tec")
BC   = os.path.join(CASE, "INPUT", "bc.txt")
REL_TOL = 1e-9

def parse_wall_tec(path):
    variables, zones, hdr, raw = [], {}, None, []
    for line in open(path):
        s = line.strip()
        if re.match(r"(?i)variables", s):
            variables = re.findall(r'"([^"]+)"', s); continue
        if re.match(r"(?i)zone", s):
            if hdr: zones[hdr["name"]] = (hdr, raw)
            raw = []
            g = lambda pat: re.search(pat, s, re.I)
            hdr = {"name": g(r'T\s*=\s*(\w+)').group(1),
                   "ni": int(g(r'\bI\s*=\s*(\d+)').group(1)), "nj": int(g(r'\bJ\s*=\s*(\d+)').group(1)),
                   "nk": int(g(r'\bK\s*=\s*(\d+)').group(1)),
                   "cc": tuple(map(int, g(r'\[(\d+)-(\d+)\]\s*=\s*CELLCENTERED').groups()))}
            continue
        try: raw.append(float(s))
        except ValueError: pass
    if hdr: zones[hdr["name"]] = (hdr, raw)
    return variables, zones

def zone_vars(variables, hdr, raw):
    ni, nj, nk = hdr["ni"], hdr["nj"], hdr["nk"]
    n_nodal = ni*nj*nk; n_cell = max(ni-1,1)*max(nj-1,1)*max(nk-1,1)
    out, idx = {}, 0
    for v, name in enumerate(variables, 1):
        n = n_cell if hdr["cc"][0] <= v <= hdr["cc"][1] else n_nodal
        out[name] = raw[idx:idx+n]; idx += n
    return out

def bc303_params(path):
    """Return {face: (hconv, qrad, Tref)} for every 303 record, asserting they agree per face."""
    lines = open(path).read().splitlines(); params = {}
    for n, l in enumerate(lines):
        f = l.split()
        if len(f) == 6 and f[5] == "303":
            b, face = int(f[0]), int(f[4])
            vals = tuple(float(x) for x in lines[n+1].replace(",", " ").split())
            key = f"B{b}F{face}"
            if key in params and params[key] != vals:
                sys.exit(f"FAIL: {key} has non-uniform 303 parameters; this gate assumes one (hconv,qrad,Tref) per face")
            params[key] = vals
    return params

def main():
    if not os.path.isfile(WALL): sys.exit(f"FAIL: {WALL} not found")
    params = bc303_params(BC)
    if not params: sys.exit("FAIL: no type-303 face in bc.txt")
    variables, zones = parse_wall_tec(WALL)
    worst, npts, ok = 0.0, 0, True
    print(f"BC 303 gate: qw - hconv*(Tref - Tw) - qrad on every 303 face cell  (tol {REL_TOL:g} rel.)")
    for zone, (hconv, qrad, Tref) in sorted(params.items()):
        if zone not in zones: sys.exit(f"FAIL: zone {zone} not in wall.tec")
        zv = zone_vars(variables, *zones[zone])
        res = [qw - hconv*(Tref - tw) - qrad for tw, qw in zip(zv["Tw"], zv["qw"])]
        m = max(abs(r) for r in res); worst = max(worst, m); npts += len(res)
        status = "ok" if m <= REL_TOL*abs(qrad) else "FAIL"
        if status == "FAIL": ok = False
        print(f"  {zone}: {len(res):4d} cells  hconv={hconv:g} qrad={qrad:g} Tref={Tref:g}  max|residual| = {m:.3e} W/m2  [{status}]")
    print(f"  {npts} cells, worst residual {worst:.3e} W/m2 vs budget {REL_TOL*abs(qrad):.3e}")
    if not ok:
        print("FAIL: the qrad read from bc.txt did not reach the wall flux"); sys.exit(1)
    print("PASS")

if __name__ == "__main__":
    main()
