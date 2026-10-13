#!/bin/bash
#===============================================================================
#         FILE: test/run_restart_gate.sh
#        USAGE: ./test/run_restart_gate.sh [-w WORKDIR]
#  DESCRIPTION: A restarted run must reproduce a continuous one across the
#               restart point (plan 05 section 6).
#
#               For each case it runs 2N iterations straight through, then N
#               iterations followed by a restart for another N, and compares the
#               two solutions field by field.
#
#               Five cases, because they fail for unrelated reasons:
#
#                 qvol_Plate         static mesh WITH a volumetric source. The
#                                    source file name used to be assigned only
#                                    on the newrun branch of
#                                    Setup_Input_Solution, so a restart opened
#                                    an undefined filename, the read failed and
#                                    qvol was silently zeroed for the rest of
#                                    the run. Nothing reported it. This case is
#                                    the gate on that.
#
#                 receding/N040      MOVING mesh. The solution file has to carry
#                                    the current node positions, or the restart
#                                    silently resets the geometry to the initial
#                                    condition and continues on the wrong mesh.
#                                    Uses law = translation: it displaces
#                                    incrementally, so the mesh IS the whole
#                                    state. The oscillatory `prescribed` law
#                                    anchors to a cached reference mesh and is
#                                    refused on restart by input validation --
#                                    see Check_Mesh_Motion_Compatibility.
#
#                 pyrolysis/energy_audit  static mesh with the reaction state
#                                    chi_* written to and read back from the
#                                    solution file (plan 06 section 3.6); the
#                                    gate compares every chi_* column too.
#                 pyrolysis/darcy/mass_balance_1d  the type-2 gas closure (plan 07,
#                                    gate 4.1(b2), time-accurate): the gas
#                                    pressure p_pg, the stored-gas state
#                                    eta_st and the two history fields of the
#                                    extrapolated coefficients (eta_prev,
#                                    T_prev) cross the restart; the gate
#                                    compares them and the gas velocity. T
#                                    is uniform and constant there and pi is a
#                                    constant source, so the two history
#                                    fields cannot matter in this case:
#                 pyrolysis/darcy/heated0d_type2  heated, decomposing, type 2:
#                                    the production of each step is evaluated
#                                    at T* = T^n + (T^n - T^(n-1))/2 and the
#                                    storage coefficient extrapolated from
#                                    eta^(n-1), so a restart that loses or
#                                    shifts T_prev or eta_prev moves the
#                                    answer. A T_prev column holding the
#                                    in-memory value (two steps old at the end
#                                    of a step) instead of TO failed it
#                                    (2026-10-09: rel|dT| 2.7e-8). Fixed, the
#                                    case sits at ~1e-12 (T) to ~6e-12 (chi),
#                                    above the ~1e-15 of the others: the type-2
#                                    inversion stops at |f| <= 1e-13 |E|, and
#                                    the restart's rounded guess can stop one
#                                    Newton iterate apart.
#
#               WHY THE TOLERANCE IS NOT ZERO
#               ----------------------------
#               Bit-identity was the first criterion tried and it is wrong: the
#               restart goes through the ASCII solution file, which carries 15
#               significant decimal digits, so the state is perturbed by ~5e-16
#               relative on the way in. Both cases then differ by a few times
#               1e-15 relative -- the round-trip floor, amplified slightly by
#               the remaining steps -- and not by anything structural.
#
#               The budget below is 1e-11 relative: four orders above that
#               measured floor, and many orders below what losing state costs.
#               For scale, the two defects this gate exists for produce
#               relative differences of 1.6e-1 and 4.6e-3 respectively, both
#               verified by reintroducing them.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK="${TMPDIR:-/tmp}/fuss-restart-gate.$$"

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR]"; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$WORK"
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS restart gate"
echo "  work dir : $WORK"
echo

# $1 = case path relative to the repo root, $2 = iterations per half
run_case () {
  local case_rel=$1 half=$2
  local src="$ROOT/$case_rel"
  local name; name=$(echo "$case_rel" | tr '/' '_')
  local cont="$WORK/$name.cont" rst="$WORK/$name.rst"

  rm -rf "$cont" "$rst"
  for d in "$cont" "$rst"; do
    mkdir -p "$d"
    cp -r "$src/INPUT" "$src/input.ini" "$src/FUSS.sh" "$d/" 2>/dev/null
    [ -d "$src/MESH" ] && cp -r "$src/MESH" "$d/"
    # the BC 302 time file of a pyrolysis case lives at the case root
    cp "$src"/*.dat "$d/" 2>/dev/null || true
    # The copies sit one level deeper than any real case, so point FUSS.sh at
    # the master binary by absolute path instead of counting ../ levels.
    sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
    mkdir -p "$d/OUTPUT" "$d/bin"
  done

  set_iters () {  # $1 = dir, $2 = iterations, $3 = newrun
    python3 - "$1/input.ini" "$2" "$3" <<'PY'
import re, sys
p, it, newrun = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(p).read()
def setkey(s, key, val):
    pat = re.compile(r"^(%s\s*=).*$" % re.escape(key), re.M)
    return pat.sub(r"\1 " + val, s) if pat.search(s) else s
s = setkey(s, "newrun", newrun)
s = setkey(s, "iter-threshold", it)
# Convergence and wall-clock limits must not end a leg early, or the two runs
# would stop at different points and the comparison would be meaningless.
s = setkey(s, "res-threshold", "1e-30")
s = setkey(s, "time-threshold", "1e30")
s = setkey(s, "sol-overwrite", "true")
if "iter-threshold" not in s:
    s = s.replace("[FUSS-Parameters]", "[FUSS-Parameters]\niter-threshold = " + it, 1)
open(p, "w").write(s)
PY
  }

  set_iters "$cont" "$((2*half))" "true"
  ( cd "$cont" && ./FUSS.sh -p 1 solve ) > "$cont/run.log" 2>&1 \
    || { bad "$case_rel: continuous run failed (see $cont/run.log)"; return; }

  set_iters "$rst" "$half" "true"
  ( cd "$rst" && ./FUSS.sh -p 1 solve ) > "$rst/run1.log" 2>&1 \
    || { bad "$case_rel: first leg failed (see $rst/run1.log)"; return; }

  set_iters "$rst" "$half" "false"
  ( cd "$rst" && ./FUSS.sh -p 1 solve ) > "$rst/run2.log" 2>&1 \
    || { bad "$case_rel: restarted leg failed (see $rst/run2.log)"; return; }

  local out rc
  out=$(python3 - "$cont/OUTPUT/field.tec" "$rst/OUTPUT/field.tec" <<'PY'
import re, sys

def read_tec(path):
    lines = open(path).read().split("\n")
    names = [t.strip('"') for t in lines[0].split("=", 1)[-1].split() if t.strip('"')]  # VARIABLES names, quoted (ORION <= v1.6) or bare (ORION >= v1.7.0 writes `VARIABLES = x y z T ...`)
    m = re.search(r"I=(\d+),\s*J=(\d+),\s*K=(\d+)", lines[1])
    ni, nj, nk = (int(x) for x in m.groups())
    vals = []
    for line in lines[2:]:
        s = line.strip()
        if not s:
            continue
        try:
            vals.append(float(s))
        except ValueError:
            pass
    n_nodal = ni * nj * nk
    n_cell = (ni - 1) * (nj - 1) * max(nk - 1, 1)
    out, off = {}, 0
    for name in names:
        n = n_nodal if name in ("x", "y", "z") else n_cell
        out[name] = vals[off:off + n]
        off += n
    return out

REL_TOL = 1.0e-11   # see the header: the ASCII round-trip floor is ~1e-15

a, b = read_tec(sys.argv[1]), read_tec(sys.argv[2])
worst, ok = [], True
compared = set()
# Keys from the header itself, not a literal list: a misspelt literal would pass
# silently. The pyrolysis state (chi_*) and the gas flux travel through the same
# file (plan 06 section 3.6) and are compared whenever they are present; so
# is the type-2 state (p_pg, eta_st, eta_prev, T_prev; plan 07). Under type 2
# the gas flux mdot_pg and the velocity u_pg_* are DERIVED from differences of
# p (Darcy), and their round-trip floor is the 15-digit quantum of p over the
# cell-to-cell difference: ~5e-12 on mass_balance_1d (dp_cell/p ~ 2e-4), ~1e-8
# on heated0d_type2 (dp_cell/p ~ 1e-7, measured 1.4e-10). They are printed as
# info there, not gated; type 1's mdot_pg (an integral of the production) stays
# gated.
type2 = "p_pg" in a
derived = [k for k in a if type2 and (k == "mdot_pg" or k.startswith("u_pg_"))]
keys = [k for k in a if k not in derived and (k in ("T", "qvol", "x", "y", "z") or k.startswith("chi_")
        or k == "mdot_pg" or k in ("p_pg", "eta_st", "eta_prev", "T_prev"))]
info = []
for key in keys + derived:
    if key not in a or key not in b:
        continue
    compared.add(key)
    if len(a[key]) != len(b[key]):
        print("SIZE-MISMATCH %s" % key)
        sys.exit(2)
    d = max((abs(u - v) for u, v in zip(a[key], b[key])), default=0.0)
    if key.startswith("u_pg_"):
        # a velocity component is scaled by the whole vector: the y and z
        # components of a 1-D case are round-off, and scaled by themselves
        # their round-trip difference reads as O(1)
        scale = max((abs(v) for k2 in a if k2.startswith("u_pg_") for v in a[k2]), default=0.0)
    else:
        scale = max((abs(v) for v in a[key]), default=0.0)
    rel = d / scale if scale > 0.0 else d
    if key in derived:
        info.append((key, rel))
        continue
    worst.append((key, rel))
    if rel > REL_TOL:
        ok = False

# Reported so a zero source is visible as a zero source rather than as a
# coincidental match: two runs that BOTH lost qvol would agree perfectly.
qmax = max((abs(v) for v in a.get("qvol", [0.0])), default=0.0)

if "T" not in compared or "x" not in compared:
    print("MISSING-NAMES: compared %s -- the VARIABLES header was not parsed" % sorted(compared))
    sys.exit(2)
print("  ".join("rel|d%s|=%.2e" % kv for kv in worst) + "   qvol_peak=%.3e" % qmax
      + ("   info (derived from p, not gated): " + "  ".join("rel|d%s|=%.2e" % kv for kv in info) if info else ""))
sys.exit(0 if ok else 1)
PY
)
  rc=$?
  printf '        %s\n' "$out"
  if [ $rc -eq 0 ]; then
    pass "$case_rel: restart reproduces the continuous run (within the ASCII floor)"
  else
    bad "$case_rel: restarted run differs from the continuous run"
  fi
}

run_case "test/steady-state/qvol_Plate"      200
run_case "test/moving-mesh/receding/N040"    200
# pyrolysis (plan 06): the reaction state chi_* must travel through the restart
# file with T; the budget accumulators restart from zero by design and are not
# compared. 100 iterations per half: the legs run to ~34 s of simulated time
# (the case itself stops at 10 s; the time file holds 1644 K to 60 s).
run_case "test/pyrolysis/energy_audit"       100
# Darcy (plan 07 section 3.8): 4.1(b2) is time-accurate (quasi-steady = false),
# so the gas pressure carries history across the restart; 200 steps per half
# run to ~17 s, far past the case's own 2 s.
run_case "test/pyrolysis/darcy/mass_balance_1d" 200
# the history fields of the type-2 step (T_prev, eta_prev) under heating and
# kinetics: 200 steps per half (~0.4 s each, 300 K to ~560 K, reacting)
run_case "test/pyrolysis/darcy/heated0d_type2"  200
# quasi-steady type 2 with threshold kinetics (2026-10-13): the gas the stages
# released beyond the prediction waits in (eta_st - eta_prev) p_pg, all three
# restart columns; 800 steps per half (~136 s each, inside the 300 s case) cross many cell fronts
run_case "test/pyrolysis/lit/baum_atreya2007/qs_N100" 800

echo
if [ "$fail" -eq 0 ]; then
  echo "restart gate passed"
else
  echo "$fail restart check(s) FAILED"
  exit 1
fi
