#!/usr/bin/env python3
"""
make_tacot_tables.py -- write the FUSS pyrolysis material tables from PATO 3.1's TACOT data.

Implements plan 06 (phase 2, pyrolysis on a static mesh) section 3.2, the DECIDED block
"who generates the files": until ATLAS grows a [GPB-Pyrolysis] builder, the four files a
pyrolysis case needs in INPUT/ are written by this script from the PATO 3.1 material
directory data/Materials/Composites/TACOT:

    pyrolysis.dat   phase constants and kinetics, PATO key names, one `material <name>` block
    virgin.dat      T cp h ki kj kk emissivity absorptivity   (nvar = 7), 1 K rows
    char.dat        same layout, same T range as virgin.dat
    gas.dat         T p M h mu                                (nvar = 4), one zone per p-level

Tecplot layout (what ORION's tec_read_points_multivars accepts, plan 06 section 3.2):
line 1 `TITLE = "FUSS pyrolysis table"`, line 2 `VARIABLES = "T" ...`, then per zone
`ZONE T="<material>"` (gas: `ZONE T="<material> p=<p>"`), `I=<rows>, F=POINT`, and the rows,
blank-separated reals; nothing after the last zone's last row.

Resampling (plan 06 section 3.2, DECIDED 2026-10-06): every column except the virgin/char h is
resampled linearly onto integer T; the virgin/char h is rebuilt from PATO's linear cp on each
PATO interval [T_a, T_b] as
    h(T) = h_a + int_{T_a}^{T} cp_lin dT' + r_ab (T - T_a)/(T_b - T_a),
    r_ab = h_b - h_a - (cp_a + cp_b)(T_b - T_a)/2,
so h hits PATO's nodes exactly and its slope is PATO's cp. The gas h stays linear.

Usage:
    make_tacot_tables.py --pato-dir <PATO 3.1 root> --out <dir>
        [--isothermal | --constant-gas M MU H | --zero-emissivity] [--material NAME]

    --isothermal       plan 06 section 4.2 test tables: solid cp = 1000, h = 1000 (T - 298.15),
                       k = 0.4; gas h = 1000 (T - 298.15), M = 0.022, mu = 3e-5; pyrolysis.dat =
                       TACOT constants with eps_g_v = eps_g_c = 0.8. The solid h and the gas h
                       come from ONE function and ONE format string (h_iso, FMT), so the two files
                       carry byte-identical h columns at equal T.
    --constant-gas     gas.dat whose every row of every p-level carries the constants M, mu, h
                       (plan 07's synthetic gas).
    --zero-emissivity  a full copy of PATO's TACOT directory into <dir> with column 8 (emissivity)
                       of `virgin` and `char` set to 0.0, every other byte unchanged (PATO cross-run
                       of plan 06 section 4.5). Writes no FUSS table.

Self-checks (TACOT datum checks of the plan) are printed as PASS/FAIL lines; exit 1 on any FAIL.
Python 3.6, no third-party modules.
"""

import argparse
import datetime
import math
import os
import re
import shutil
import sys

FMT = "%.9E"        # the ONE format string for every tabulated value (10 significant digits)
TFMT = "%.1f"       # T as an integer-valued real, e.g. 256.0
TOL_INT = 1e-9      # |T - round(T)| tolerance for "integer-valued"

TACOT_REL = os.path.join("data", "Materials", "Composites", "TACOT")

SOLID_VARS = 'VARIABLES = "T" "cp" "h" "ki" "kj" "kk" "emissivity" "absorptivity"'
GAS_VARS = 'VARIABLES = "T" "p" "M" "h" "mu"'
TITLE = 'TITLE = "FUSS pyrolysis table"'

# --isothermal constants (plan 06 section 4.2)
ISO_CP = 1000.0
ISO_K = 0.4
ISO_M = 0.022
ISO_MU = 3.0e-5
ISO_EPS_G = "0.8"


def h_iso(T):
    """The ONE isothermal enthalpy function, shared by the solid and the gas writers."""
    return 1000.0 * (T - 298.15)


# ----------------------------------------------------------------------------- helpers
def die(msg):
    sys.stderr.write("[ERROR] make_tacot_tables: %s\n" % msg)
    sys.exit(1)


def is_int_valued(x):
    return abs(x - round(x)) <= TOL_INT


def first_integer_at_or_above(x):
    """First integer inside [x, ...): ceil(x), with a tolerance for x already integer."""
    if is_int_valued(x):
        return int(round(x))
    return int(math.ceil(x))


def last_integer_at_or_below(x):
    if is_int_valued(x):
        return int(round(x))
    return int(math.floor(x))


def read_pato_table(path, ncol):
    """Read a PATO `//`-commented numeric table. Returns the rows (lists of floats) in file order."""
    rows = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("//"):
                continue
            toks = s.split()
            try:
                vals = [float(t) for t in toks]
            except ValueError:
                die("%s line %d: non-numeric data line: %r" % (path, ln, line.rstrip("\n")))
            if len(vals) != ncol:
                die("%s line %d: expected %d columns, found %d" % (path, ln, ncol, len(vals)))
            rows.append(vals)
    if not rows:
        die("%s: no data rows" % path)
    return rows


def group_by_pressure(rows):
    """Split rows into consecutive runs of equal p (column 0). Returns [(p, [row, ...]), ...]."""
    groups = []
    for r in rows:
        if groups and groups[-1][0] == r[0]:
            groups[-1][1].append(r)
        else:
            groups.append((r[0], [r]))
    return groups


def check_ascending(Ts, what):
    for a, b in zip(Ts, Ts[1:]):
        if not b > a:
            die("%s: T not strictly ascending (%g then %g)" % (what, a, b))


def lerp(Tn, vn, T, bisect=True):
    """Linear interpolation of the nodes (Tn ascending, vn) at T, T inside [Tn[0], Tn[-1]].

    bisect=True uses a binary search for the interval, bisect=False a linear scan: the same
    arithmetic through two independent lookups (the self-check compares the two).
    """
    n = len(Tn)
    if T < Tn[0] or T > Tn[-1]:
        die("interpolation outside the table: T = %g not in [%g, %g]" % (T, Tn[0], Tn[-1]))
    if bisect:
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if Tn[mid] <= T:
                lo = mid
            else:
                hi = mid
        a = lo
    else:
        a = 0
        while a < n - 2 and Tn[a + 1] <= T:
            a += 1
    b = a + 1
    if T == Tn[a]:
        return vn[a]
    if T == Tn[b]:
        return vn[b]
    return vn[a] + (vn[b] - vn[a]) * (T - Tn[a]) / (Tn[b] - Tn[a])


def h_from_cp(Tn, cpn, hn, T):
    """Plan 06 section 3.2 interval rule: h on [T_a, T_b] from PATO's linear cp, pinned at the nodes."""
    n = len(Tn)
    lo, hi = 0, n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if Tn[mid] <= T:
            lo = mid
        else:
            hi = mid
    a, b = lo, lo + 1
    if T == Tn[a]:
        return hn[a]
    if T == Tn[b]:
        return hn[b]
    Ta, Tb = Tn[a], Tn[b]
    cpa, cpb = cpn[a], cpn[b]
    dT = Tb - Ta
    r_ab = hn[b] - hn[a] - 0.5 * (cpa + cpb) * dT
    x = T - Ta
    integral = cpa * x + 0.5 * (cpb - cpa) * x * x / dT
    return hn[a] + integral + r_ab * x / dT


def residual_slopes(Tn, cpn, hn):
    """|r_ab|/(T_b - T_a) per PATO interval, J/kg/K (plan: <= 6.9 for TACOT)."""
    out = []
    for a in range(len(Tn) - 1):
        dT = Tn[a + 1] - Tn[a]
        r = hn[a + 1] - hn[a] - 0.5 * (cpn[a] + cpn[a + 1]) * dT
        out.append(abs(r) / dT)
    return out


def fmt_row(T, values):
    return TFMT % T + " " + " ".join(FMT % v for v in values)


def write_table(path, varline, zones):
    """zones: list of (title, rows); rows: list of (T, [values...])."""
    lines = [TITLE, varline]
    for title, rows in zones:
        lines.append('ZONE T="%s"' % title)
        lines.append("I=%d, F=POINT" % len(rows))
        for T, vals in rows:
            lines.append(fmt_row(T, vals))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def read_fuss_table(path, nvar):
    """Re-read a written table the way the ORION reader does (zone count by substring, rows list-directed)."""
    zones = []
    cur = None
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if "ZONE" in s or "Zone" in s:
                cur = {"title": s, "rows": []}
                zones.append(cur)
                continue
            toks = s.split()
            try:
                vals = [float(t) for t in toks]
            except ValueError:
                continue
            if len(vals) != nvar + 1:
                die("%s: row with %d values, expected %d" % (path, len(vals), nvar + 1))
            cur["rows"].append(vals)
    return zones


# ----------------------------------------------------------------------------- constantProperties
STMT_RE = re.compile(r"^\s*(\S+)\s+(?:\S+\s+\[[^\]]*\]\s+)?([-+0-9.eE]+)\s*;")
KEY_SCALAR = ("nSolidPhases", "eps_g_v", "eps_g_c")
KEY_PHASE_RE = re.compile(r"^(rhoI|epsI|nPyroReac)\[(\d+)\]$")
KEY_REAC_RE = re.compile(r"^(F|A|E|m|n|T|h)\[(\d+)\]\[(\d+)\]$")


def parse_constant_properties(path):
    """Return an ordered dict-like list of (key, value-text) for the pyrolysis.dat grammar keys only."""
    with open(path) as fh:
        text = fh.read()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//[^\n]*", "", text)
    found = {}
    order = []
    for ln, line in enumerate(text.split("\n"), 1):
        m = STMT_RE.match(line)
        if not m:
            continue
        key, val = m.group(1), m.group(2)
        if key in KEY_SCALAR or KEY_PHASE_RE.match(key) or KEY_REAC_RE.match(key):
            if key in found:
                die("%s: key %s appears twice" % (path, key))
            found[key] = val
            order.append(key)
    if "nSolidPhases" not in found:
        die("%s: nSolidPhases not found" % path)
    return found


def build_constants(found, path):
    """Validate the grammar (every key exactly once) and return the structured constants."""
    nph = int(found["nSolidPhases"])
    if nph < 1:
        die("%s: nSolidPhases must be >= 1" % path)
    for k in ("eps_g_v", "eps_g_c"):
        if k not in found:
            die("%s: %s not found" % (path, k))
    phases = []
    for i in range(1, nph + 1):
        ph = {}
        for k in ("rhoI", "epsI", "nPyroReac"):
            key = "%s[%d]" % (k, i)
            if key not in found:
                die("%s: %s not found" % (path, key))
            ph[k] = found[key]
        nre = int(ph["nPyroReac"])
        reacs = []
        for j in range(1, nre + 1):
            rc = {}
            for k in ("F", "A", "E", "m", "n", "T", "h"):
                key = "%s[%d][%d]" % (k, i, j)
                if key not in found:
                    die("%s: %s not found" % (path, key))
                rc[k] = found[key]
            reacs.append(rc)
        ph["reactions"] = reacs
        phases.append(ph)
    # keys present in the file but outside the declared counts are an error (exactly once each)
    used = set(["nSolidPhases", "eps_g_v", "eps_g_c"])
    for i, ph in enumerate(phases, 1):
        for k in ("rhoI", "epsI", "nPyroReac"):
            used.add("%s[%d]" % (k, i))
        for j in range(len(ph["reactions"])):
            for k in ("F", "A", "E", "m", "n", "T", "h"):
                used.add("%s[%d][%d]" % (k, i, j + 1))
    extra = sorted(set(found) - used)
    if extra:
        die("%s: keys outside the declared phase/reaction counts: %s" % (path, ", ".join(extra)))
    return {"nSolidPhases": nph, "eps_g_v": found["eps_g_v"], "eps_g_c": found["eps_g_c"],
            "phases": phases}


def write_pyrolysis_dat(path, name, const, provenance, eps_override=None):
    lines = ["// FUSS pyrolysis constants (plan 06 section 3.2 grammar): one `key value` per line,",
             "// `//` comments allowed. Written by make_tacot_tables.py from",
             "// %s" % provenance,
             "// on %s%s" % (datetime.date.today().isoformat(),
                             "" if eps_override is None else
                             " with eps_g_v = eps_g_c = %s (--isothermal)" % eps_override),
             "material %s" % name,
             "nSolidPhases %d" % const["nSolidPhases"]]
    for i, ph in enumerate(const["phases"], 1):
        lines.append("rhoI[%d] %s" % (i, ph["rhoI"]))
        lines.append("epsI[%d] %s" % (i, ph["epsI"]))
        lines.append("nPyroReac[%d] %s" % (i, ph["nPyroReac"]))
    eg_v = const["eps_g_v"] if eps_override is None else eps_override
    eg_c = const["eps_g_c"] if eps_override is None else eps_override
    lines.append("eps_g_v %s" % eg_v)
    lines.append("eps_g_c %s" % eg_c)
    for i, ph in enumerate(const["phases"], 1):
        for j, rc in enumerate(ph["reactions"], 1):
            for k in ("F", "A", "E", "m", "n", "T", "h"):
                lines.append("%s[%d][%d] %s" % (k, i, j, rc[k]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def derived_densities(const):
    rho_v = 0.0
    rho_c = 0.0
    for ph in const["phases"]:
        eps_rho = float(ph["epsI"]) * float(ph["rhoI"])
        rho_v += eps_rho
        sumF = sum(float(rc["F"]) for rc in ph["reactions"])
        rho_c += eps_rho * (1.0 - sumF)
    return rho_v, rho_c


# ----------------------------------------------------------------------------- solid tables
def load_solid(path):
    """PATO virgin/char: p T cp h ki kj kk emissivity absorptivity, several p sub-tables that must agree."""
    rows = read_pato_table(path, 9)
    groups = group_by_pressure(rows)
    ref = groups[0][1]
    for p, g in groups[1:]:
        if len(g) != len(ref):
            die("%s: p-level %g has %d rows, p-level %g has %d" % (path, p, len(g), groups[0][0], len(ref)))
        for ra, rb in zip(ref, g):
            if ra[1:] != rb[1:]:
                die("%s: the p-levels %g and %g differ at T = %g; refusing (plan 06 section 3.2)"
                    % (path, groups[0][0], p, ra[1]))
    Tn = [r[1] for r in ref]
    check_ascending(Tn, path)
    cols = {"cp": [r[2] for r in ref], "h": [r[3] for r in ref],
            "ki": [r[4] for r in ref], "kj": [r[5] for r in ref], "kk": [r[6] for r in ref],
            "emissivity": [r[7] for r in ref], "absorptivity": [r[8] for r in ref]}
    return Tn, cols, [p for p, _ in groups]


def resample_solid(Tn, cols, isothermal):
    Tlo = first_integer_at_or_above(Tn[0])
    Thi = last_integer_at_or_below(Tn[-1])
    rows = []
    for k in range(Tlo, Thi + 1):
        T = float(k)
        if isothermal:
            vals = [ISO_CP, h_iso(T), ISO_K, ISO_K, ISO_K,
                    lerp(Tn, cols["emissivity"], T), lerp(Tn, cols["absorptivity"], T)]
        else:
            vals = [lerp(Tn, cols["cp"], T),
                    h_from_cp(Tn, cols["cp"], cols["h"], T),
                    lerp(Tn, cols["ki"], T), lerp(Tn, cols["kj"], T), lerp(Tn, cols["kk"], T),
                    lerp(Tn, cols["emissivity"], T), lerp(Tn, cols["absorptivity"], T)]
        rows.append((T, vals))
    return rows


# ----------------------------------------------------------------------------- gas table
def load_gas(path):
    """PATO gasProperties: p T M h_g nu, p-levels listed descending; returned sorted by ascending p."""
    rows = read_pato_table(path, 5)
    groups = group_by_pressure(rows)
    ps = [p for p, _ in groups]
    if len(set(ps)) != len(ps):
        die("%s: a p-level appears in two separate blocks" % path)
    groups.sort(key=lambda g: g[0])
    out = []
    for p, g in groups:
        Tn = [r[1] for r in g]
        check_ascending(Tn, "%s p=%g" % (path, p))
        out.append((p, Tn, {"M": [r[2] for r in g], "h": [r[3] for r in g], "mu": [r[4] for r in g]}))
    return out


def resample_gas(levels, mode, const_gas=None):
    """mode: 'tacot' | 'isothermal' | 'constant'. Returns [(title-suffix p, rows)], common T range."""
    Tlo = max(first_integer_at_or_above(Tn[0]) for _, Tn, _ in levels)
    Thi = min(last_integer_at_or_below(Tn[-1]) for _, Tn, _ in levels)
    zones = []
    for p, Tn, cols in levels:
        rows = []
        for k in range(Tlo, Thi + 1):
            T = float(k)
            if mode == "isothermal":
                vals = [p, ISO_M, h_iso(T), ISO_MU]
            elif mode == "constant":
                M, MU, H = const_gas
                vals = [p, M, H, MU]
            else:
                vals = [p, lerp(Tn, cols["M"], T), lerp(Tn, cols["h"], T), lerp(Tn, cols["mu"], T)]
            rows.append((T, vals))
        zones.append((p, rows))
    return zones


# ----------------------------------------------------------------------------- zero emissivity
def zero_emissivity_copy(src_dir, out_dir):
    """Copy the whole TACOT directory; in `virgin` and `char` set the 8th token of every data line to 0.0."""
    os.makedirs(out_dir, exist_ok=True)
    edited = {}
    for name in sorted(os.listdir(src_dir)):
        src = os.path.join(src_dir, name)
        dst = os.path.join(out_dir, name)
        if os.path.isdir(src):
            if os.path.exists(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            continue
        if name in ("virgin", "char"):
            with open(src, newline="") as fh:
                text = fh.read()
            out_lines = []
            n_edit = 0
            for line in text.splitlines(True):
                parts = re.split(r"(\s+)", line)      # tokens and separators alternate, all bytes kept
                toks = [i for i, t in enumerate(parts) if t and not t.isspace()]
                first = parts[toks[0]] if toks else ""
                is_data = False
                try:
                    float(first)
                    is_data = True
                except ValueError:
                    is_data = False
                if is_data:
                    if len(toks) != 9:
                        die("%s: data line with %d tokens, expected 9" % (src, len(toks)))
                    parts[toks[7]] = "0.0"
                    n_edit += 1
                out_lines.append("".join(parts))
            with open(dst, "w", newline="") as fh:
                fh.write("".join(out_lines))
            shutil.copystat(src, dst)
            edited[name] = n_edit
        else:
            shutil.copy2(src, dst)
    return edited


# ----------------------------------------------------------------------------- checks
class Checker(object):
    def __init__(self):
        self.fails = 0

    def check(self, ok, label, detail):
        print("%s  %s: %s" % ("PASS" if ok else "FAIL", label, detail))
        if not ok:
            self.fails += 1


def validate_material(name):
    if not name or len(name) > 16:
        die("material name must be 1..16 characters (FUSS clen = 16): %r" % name)
    if re.search(r"\s", name):
        die("material name must not contain whitespace: %r" % name)
    for bad in ("I=", "ZONE", "Zone", "DATA"):
        if bad in name:
            die("material name must not contain %r (it breaks the ORION table reader): %r" % (bad, name))


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage="%(prog)s --pato-dir DIR --out DIR [--material NAME] "
                                       "[--isothermal | --constant-gas M MU H | --zero-emissivity]")
    ap.add_argument("--pato-dir", required=True, help="PATO 3.1 root (holds data/Materials/Composites/TACOT)")
    ap.add_argument("--out", required=True, help="output directory (a case's INPUT/, or the TACOT copy)")
    ap.add_argument("--material", default="TACOT", help="material name (ZONE label, `material` block)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--isothermal", action="store_true", help="plan 06 section 4.2 test tables")
    g.add_argument("--constant-gas-set", action="store_true", help=argparse.SUPPRESS)
    g.add_argument("--zero-emissivity", action="store_true",
                   help="copy PATO's TACOT directory with emissivity = 0 in virgin and char")
    # `--constant-gas M MU H`: pulled out of argv by hand so that a negative H such as -7.09e6
    # (which argparse would take for an option flag) is accepted.
    argv = list(sys.argv[1:])
    const_gas = None
    if "--constant-gas" in argv:
        i = argv.index("--constant-gas")
        vals = argv[i + 1:i + 4]
        if len(vals) != 3:
            die("--constant-gas needs three values: M MU H")
        try:
            const_gas = [float(v) for v in vals]
        except ValueError:
            die("--constant-gas values must be reals: %r" % vals)
        argv = argv[:i] + ["--constant-gas-set"] + argv[i + 4:]
    args = ap.parse_args(argv)
    args.constant_gas = const_gas

    tacot = os.path.join(args.pato_dir, TACOT_REL)
    for f in ("constantProperties", "virgin", "char", "gasProperties"):
        if not os.path.isfile(os.path.join(tacot, f)):
            die("%s not found under %s" % (f, tacot))
    validate_material(args.material)
    os.makedirs(args.out, exist_ok=True)
    chk = Checker()

    if args.zero_emissivity:
        edited = zero_emissivity_copy(tacot, args.out)
        print("zero-emissivity copy of %s -> %s" % (tacot, args.out))
        for name in ("virgin", "char"):
            orig = read_pato_table(os.path.join(tacot, name), 9)
            copy = read_pato_table(os.path.join(args.out, name), 9)
            same_rows = len(orig) == len(copy)
            others_equal = same_rows and all(a[:7] == b[:7] and a[8] == b[8] for a, b in zip(orig, copy))
            eps_zero = same_rows and all(b[7] == 0.0 for b in copy)
            chk.check(same_rows and others_equal and eps_zero, "%s emissivity column zeroed" % name,
                      "%d data lines edited, every other column equal, %d rows" % (edited[name], len(copy)))
            with open(os.path.join(tacot, name), "rb") as fa, open(os.path.join(args.out, name), "rb") as fb:
                ba, bb = fa.read(), fb.read()
            chk.check(len(ba) == len(bb), "%s byte count preserved" % name, "%d -> %d bytes" % (len(ba), len(bb)))
        copied = sorted(os.listdir(args.out))
        chk.check(set(copied) >= set(os.listdir(tacot)), "directory copied whole", ", ".join(copied))
        sys.exit(1 if chk.fails else 0)

    mode = "isothermal" if args.isothermal else ("constant" if args.constant_gas else "tacot")
    prov = os.path.join(tacot, "constantProperties")
    found = parse_constant_properties(prov)
    const = build_constants(found, prov)
    rho_v, rho_c = derived_densities(const)

    Tv, cv, pv = load_solid(os.path.join(tacot, "virgin"))
    Tc, cc, pc = load_solid(os.path.join(tacot, "char"))
    gas_levels = load_gas(os.path.join(tacot, "gasProperties"))

    virgin_rows = resample_solid(Tv, cv, args.isothermal)
    char_rows = resample_solid(Tc, cc, args.isothermal)
    gas_zones = resample_gas(gas_levels, mode, args.constant_gas)

    name = args.material
    out = args.out
    write_pyrolysis_dat(os.path.join(out, "pyrolysis.dat"), name, const, prov,
                        eps_override=ISO_EPS_G if args.isothermal else None)
    write_table(os.path.join(out, "virgin.dat"), SOLID_VARS, [(name, virgin_rows)])
    write_table(os.path.join(out, "char.dat"), SOLID_VARS, [(name, char_rows)])
    write_table(os.path.join(out, "gas.dat"), GAS_VARS,
                [("%s p=%s" % (name, FMT % p), rows) for p, rows in gas_zones])

    print("make_tacot_tables.py: mode = %s, material = %s, PATO = %s" % (mode, name, tacot))
    print("  virgin p-levels %s, char p-levels %s, gas p-levels (ascending) %s"
          % (pv, pc, [p for p, _ in gas_zones]))
    print("  wrote %s" % ", ".join(os.path.join(out, f) for f in ("pyrolysis.dat", "virgin.dat", "char.dat", "gas.dat")))
    rs_v = residual_slopes(Tv, cv["cp"], cv["h"])
    rs_c = residual_slopes(Tc, cc["cp"], cc["h"])
    print("  info: max |r_ab|/(T_b - T_a) = %.3f (virgin) %.3f (char) J/kg/K (plan: <= 6.9)" % (max(rs_v), max(rs_c)))

    # ------------------------------------------------------------------ self-checks
    print("self-checks (%s):" % ("real TACOT datum checks" if mode == "tacot" else "structural checks, mode " + mode))
    vz = read_fuss_table(os.path.join(out, "virgin.dat"), 7)
    cz = read_fuss_table(os.path.join(out, "char.dat"), 7)
    gz = read_fuss_table(os.path.join(out, "gas.dat"), 4)

    def h_at(rows, T):
        for r in rows:
            if r[0] == T:
                return r[2]
        return None

    nv, nc = len(vz[0]["rows"]), len(cz[0]["rows"])
    chk.check(len(vz) == 1 and len(cz) == 1, "virgin/char one zone each", "%d / %d zones" % (len(vz), len(cz)))
    Tlo_v, Thi_v = vz[0]["rows"][0][0], vz[0]["rows"][-1][0]
    Tlo_c, Thi_c = cz[0]["rows"][0][0], cz[0]["rows"][-1][0]
    chk.check(nv == nc and Tlo_v == Tlo_c and Thi_v == Thi_c, "virgin/char rows and T range agree",
              "%d rows %g..%g K / %d rows %g..%g K" % (nv, Tlo_v, Thi_v, nc, Tlo_c, Thi_c))
    ngz = len(gz)
    ngr = [len(z["rows"]) for z in gz]
    chk.check(ngz == 5 and all(n == ngr[0] for n in ngr), "gas zones x rows",
              "%d zones x %s rows (%g..%g K)" % (ngz, ngr, gz[0]["rows"][0][0], gz[0]["rows"][-1][0]))
    gp = [z["rows"][0][1] for z in gz]
    p_const = all(all(r[1] == z["rows"][0][1] for r in z["rows"]) for z in gz)
    chk.check(p_const and all(b > a for a, b in zip(gp, gp[1:])), "gas p constant per zone, ascending across zones",
              "%s" % gp)
    step_ok = True
    for z in vz + cz + gz:
        Ts = [r[0] for r in z["rows"]]
        if not is_int_valued(Ts[0]) or any(b - a != 1.0 for a, b in zip(Ts, Ts[1:])):
            step_ok = False
    chk.check(step_ok, "every zone's T starts at an integer and steps by exactly 1 K", "checked %d zones" % (len(vz) + len(cz) + len(gz)))

    if mode == "tacot":
        chk.check(nv == 3745 and nc == 3745, "virgin/char 3745 rows (256..4000 K)", "%d / %d" % (nv, nc))
        chk.check(ngz == 5 and ngr[0] == 3776, "5 gas zones x 3776 rows (200..3975 K)", "%d x %s" % (ngz, ngr))
        hc298 = h_at(cz[0]["rows"], 298.0)
        hv298 = h_at(vz[0]["rows"], 298.0)
        chk.check(hc298 is not None and abs(hc298) <= 1.0, "abs(h_c(298)) <= 1 J/kg", "%r" % hc298)
        chk.check(hv298 is not None and abs(hv298 - (-8.57e5)) <= 1.0, "h_v(298) = -8.57e5 +- 1 J/kg", "%r" % hv298)
        chk.check(abs(rho_v - 280.0) <= 1e-9, "rho_v = sum eps_I rho_I = 280", "%.12g" % rho_v)
        chk.check(abs(rho_c - 220.0) <= 1e-9, "rho_c = rho_v - sum eps_I rho_I F = 220", "%.12g" % rho_c)
        hg = None
        for z in gz:
            if abs(z["rows"][0][1] - 1.01325e5) <= 1e-6 * 1.01325e5:
                for r in z["rows"]:
                    if r[0] == 300.0:
                        hg = r[3]
        chk.check(hg is not None and abs(hg - (-7.09e6)) <= 1.0, "gas h(300 K, 1.01325e5 Pa) = -7.09e6 J/kg", "%r" % hg)
        # resampled cp vs PATO's own linear interpolant at the integers (independent interval lookup)
        for label, Tn, cols, rows in (("virgin", Tv, cv, virgin_rows), ("char", Tc, cc, char_rows)):
            dmax = 0.0
            for T, vals in rows:
                dmax = max(dmax, abs(vals[0] - lerp(Tn, cols["cp"], T, bisect=False)))
            chk.check(dmax == 0.0, "%s: max |resampled cp - PATO linear interpolant at the integers| = 0" % label, "%r" % dmax)
            # node pinning of h
            dn = max(abs(h_from_cp(Tn, cols["cp"], cols["h"], Tk) - hk)
                     for Tk, hk in zip(Tn, cols["h"]) if is_int_valued(Tk))
            chk.check(dn == 0.0, "%s: h equals PATO's h at every integer PATO node" % label, "%r" % dn)
        for label, Tn, cols, rows, zz in (("virgin", Tv, cv, virgin_rows, vz), ("char", Tc, cc, char_rows, cz)):
            rt = 0.0
            for (T, vals), r in zip(rows, zz[0]["rows"]):
                for a, b in zip(vals, r[1:]):
                    rt = max(rt, abs(a - b) / max(abs(a), 1e-300))
            print("info  %s: file round-trip max relative error (%s) = %.3e" % (label, FMT, rt))
        for z in gz:
            p = z["rows"][0][1]
            lev = [l for l in gas_levels if l[0] == p][0]
            dmax = max(abs(r[3] - lerp(lev[1], lev[2]["h"], r[0], bisect=False)) for r in z["rows"]
                       if is_int_valued(r[0]) and int(r[0]) % 25 == 0)
            chk.check(dmax == 0.0, "gas p=%s: h equals PATO's h at the 25 K nodes" % (FMT % p), "%r" % dmax)
    else:
        # the identity gate 4.2 rests on: virgin/char h strings == gas h strings at equal T
        def h_strings(path, hcol):
            d = {}
            with open(path) as fh:
                for line in fh:
                    toks = line.split()
                    if not toks:
                        continue
                    try:
                        T = float(toks[0])
                    except ValueError:
                        continue
                    d.setdefault(toks[0], toks[hcol])
            return d
        if mode == "isothermal":
            hv = h_strings(os.path.join(out, "virgin.dat"), 2)
            hc = h_strings(os.path.join(out, "char.dat"), 2)
            hg_s = h_strings(os.path.join(out, "gas.dat"), 3)
            common = sorted(set(hv) & set(hg_s), key=float)
            ident = all(hv[k] == hg_s[k] and hc[k] == hg_s[k] for k in common)
            chk.check(ident and len(common) > 0, "isothermal: solid h and gas h byte-identical at equal T",
                      "%d common T rows (%s..%s)" % (len(common), common[0] if common else "-", common[-1] if common else "-"))
            cp_const = all(vals[0] == ISO_CP and vals[2] == ISO_K for _, vals in virgin_rows + char_rows)
            chk.check(cp_const, "isothermal: cp = 1000, k = 0.4 on every solid row", "")
            with open(os.path.join(out, "pyrolysis.dat")) as fh:
                txt = fh.read()
            chk.check("eps_g_v 0.8\n" in txt and "eps_g_c 0.8\n" in txt, "isothermal: eps_g_v = eps_g_c = 0.8", "")
        else:
            M, MU, H = args.constant_gas
            ok = all(all(r[2] == M and r[3] == H and r[4] == MU for r in z["rows"]) for z in gz)
            chk.check(ok, "constant-gas: every row carries M, h, mu = %r %r %r" % (M, H, MU), "")

    sys.exit(1 if chk.fails else 0)


if __name__ == "__main__":
    main()
