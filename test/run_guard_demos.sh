#!/bin/bash
#===============================================================================
#         FILE: test/run_guard_demos.sh
#        USAGE: ./test/run_guard_demos.sh [-w WORKDIR]
#  DESCRIPTION: Every moving-mesh guard is run against the configuration it is
#               supposed to reject, and must actually reject it.
#
#               A guard that has never been seen to fire is a comment. Two of
#               these were written and did NOT fire the first time -- the
#               IRS/multigrid checks read obj_mesh_motion%enabled, which
#               Check_Input runs too early to see -- and that was only caught by
#               trying them. Hence this file.
#
#               Both the error TEXT and the EXIT CODE are checked. The codes
#               are 1 for input validation and 2 for geometry/mesh, set by
#               plan 10 D2; before that every abort used Fortran's bare `stop`,
#               which exits 0, so a guard that fired and a run that succeeded
#               were indistinguishable by status and only the text could be
#               trusted. Checking both means a guard that stops for the WRONG
#               reason no longer passes just because the message matches.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BASE="$ROOT/test/moving-mesh/receding/N040"
BASE_PYRO="$ROOT/test/pyrolysis/energy_audit"
WORK="${TMPDIR:-/tmp}/fuss-guard-demos.$$"

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

echo "FUSS guard demonstrations"
echo "  base case : test/moving-mesh/receding/N040 (mesh motion), test/pyrolysis/energy_audit (pyrolysis)"
echo "  work dir  : $WORK"
echo

# $1 = short name, $2 = expected exit code, $3 = expected substring of the
# error, $4.. = ini edits (each "key=value", or "RAW:<text>" appended verbatim)
demo () {
  local name=$1 expect_rc=$2 expect=$3; shift 3
  local d="$WORK/$name"

  rm -rf "$d"; mkdir -p "$d"
  cp -r "$BASE/INPUT" "$BASE/MESH" "$BASE/input.ini" "$BASE/FUSS.sh" "$d/"
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"

  # Keep every demo to a couple of steps: the guards all fire at setup or on
  # the first mesh update, so a long run would only waste time.
  python3 - "$d/input.ini" "iter-threshold" "3" <<'PY'
import re, sys
p, k, v = sys.argv[1:4]
s = open(p).read()
pat = re.compile(r"^(%s\s*=).*$" % re.escape(k), re.M)
s = pat.sub(r"\1 " + v, s) if pat.search(s) else s
open(p, "w").write(s)
PY

  local kv sec
  for kv in "$@"; do
    if [ "${kv:0:4}" = "RAW:" ]; then
      printf '\n%s\n' "${kv:4}" >> "$d/input.ini"
      continue
    fi
    # SEC:<section>:key=value puts a NEW key under a named section; a bare
    # key=value edits an existing key anywhere, or adds it to [FUSS-MeshMotion].
    sec="FUSS-MeshMotion"
    if [ "${kv:0:4}" = "SEC:" ]; then
      kv=${kv:4}; sec=${kv%%:*}; kv=${kv#*:}
    fi
    python3 - "$d/input.ini" "${kv%%=*}" "${kv#*=}" "$sec" <<'PY'
import re, sys
p, k, v, sec = sys.argv[1], sys.argv[2].strip(), sys.argv[3].strip(), sys.argv[4].strip()
s = open(p).read()
pat = re.compile(r"^(%s\s*=).*$" % re.escape(k), re.M)
if pat.search(s):
    s = pat.sub(lambda m: m.group(1) + " " + v, s)
else:
    # New keys must land in the section that owns them. Dropping them into an
    # unrelated section is not an error the registry reports for a FUSS- section
    # it knows -- wait, it IS fatal ("unknown option") since d4f153c -- so the
    # demo would exit 1 with the WRONG message and the text check would catch
    # it; before d4f153c the key was simply never read and the guard appeared
    # not to fire. Either way: name the section. That happened while writing this.
    hdr = "[%s]" % sec
    if hdr not in s:
        sys.exit("demo input has no section %s to receive %s" % (hdr, k))
    s = s.replace(hdr, "%s\n%s = %s" % (hdr, k, v), 1)
open(p, "w").write(s)
PY
  done

  ( cd "$d" && ./FUSS.sh -p 1 solve ) > "$d/run.log" 2>&1
  local rc=$?

  if ! grep -qF "$expect" "$d/run.log"; then
    bad "$name: expected \"$expect\", not found in $d/run.log"
    tail -4 "$d/run.log" | sed 's/^/          /'
    return
  fi
  if [ "$rc" -ne "$expect_rc" ]; then
    bad "$name: refused correctly but exited $rc, expected $expect_rc"
    return
  fi
  pass "$name -> exit $rc, $(grep -m1 -F "$expect" "$d/run.log" | sed 's/^ *//')"
}

# ---------------------------------------------------------------------------
# Input-validation guards (fire before the first step).
# ---------------------------------------------------------------------------
demo irs        1 "not compatible with implicit residual smoothing" "irs=true"
# `irs-beta > 0` alone also enables IRS (Assign_Setup promotes
# it), but that promotion used to run AFTER Check_Input, so this exact input
# passed the guard above and smoothed the ALE-remapped residual anyway.
demo irs_beta   1 "not compatible with implicit residual smoothing" "SEC:FUSS-Numerics:irs-beta=0.5"
demo multigrid  1 "not compatible with multigrid" \
    "RAW:[FUSS-Multigrid]
levels = 2
level1-iter = 3
level2-iter = 3"
demo primitive  1 "requires integration-variables = cons"           "integration-variables=prim"
demo steady     1 "requires a time-accurate run"                    "time-accurate=false"
demo restart    1 "cannot be restarted"                             "law=prescribed" "newrun=false"

# ---------------------------------------------------------------------------
# Runtime geometry guards (fire on the first mesh update).
#
# too-fast: rigid translation never tangles a mesh -- volumes are preserved --
#           so this isolates the remap's stability limit on its own. dx/dt for
#           this case is about 0.029 m/s, so 1 m/s oversteps it by ~35x.
# tangled:  an oscillation whose amplitude dwarfs the cell size folds cells
#           inside out, which is a different failure and a different message.
# ---------------------------------------------------------------------------
demo too_fast 2 "mesh moved too far in one time step" "vel=1.0 0.0 0.0"
demo tangled  2 "mesh update failed" \
              "law=prescribed" "amp=0.05 0.05 0.0" "kx=500.0 500.0 0.0" \
              "ky=500.0 500.0 0.0" "omega=1.0"

# ---------------------------------------------------------------------------
# Input validation: options nothing reads (plan 10 D1).
# ---------------------------------------------------------------------------
# Injected into an EXISTING FUSS section, not appended as a new one. Appending
# "[FUSS-Numerics]\nvnnn = ..." does not work: a repeated section header is
# ignored wholesale, so the key never reaches the option loop and this demo
# silently tested nothing the first time it was written. The repeated-header
# case is its own check below.
demo unknown_key 1 "unknown option in input.ini" "vnnn=99.0"

demo dup_section 1 "repeated section in input.ini" \
    "RAW:[FUSS-Numerics]
vnn = 0.123"

# The OTHER direction, and the one that actually carries risk. An unknown-key
# check that is too eager rejects working input, and a FUSS input.ini is full of
# sections belonging to ATLAS -- [GRIB-*], [GPB-*], [ICB-*], [BCB-*] and
# arbitrarily named per-BC blocks. Those must pass through untouched.
demo_ok () {
  local name=$1; shift
  local d="$WORK/$name"
  rm -rf "$d"; mkdir -p "$d"
  cp -r "$BASE/INPUT" "$BASE/MESH" "$BASE/input.ini" "$BASE/FUSS.sh" "$d/"
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  sed -i 's/^iter-threshold = .*/iter-threshold = 3/' "$d/input.ini"
  mkdir -p "$d/OUTPUT" "$d/bin"
  printf '\n%s\n' "$1" >> "$d/input.ini"

  ( cd "$d" && ./FUSS.sh -p 1 solve ) > "$d/run.log" 2>&1
  local rc=$?
  if [ "$rc" -eq 0 ]; then
    pass "$name -> accepted, as it must be (exit 0)"
  else
    bad "$name: a section FUSS does not own was rejected (exit $rc)"
    grep -m2 -i "unknown option" "$d/run.log" | sed 's/^/          /'
  fi
}

demo_ok foreign_section "[GRIB-SomethingElse]
whatever = 1
[Twall-made-up]
type = wall
q = 0.0"

# ---------------------------------------------------------------------------
# Pyrolysis guards (plan 06 section 3.7), on test/pyrolysis/energy_audit.
# Same contract as demo(): both the text and the exit code (1 = input) are
# checked. Edit forms, in addition to demo()'s key=value / SEC: / RAW:
#   RM:<path>                       remove a file of the copy (relative path)
#   PYRO:<key> <value>              edit a key of INPUT/pyrolysis.dat
#   BCFACE:<face>:<type>:<line2>    retype every bc.txt record of that face and
#                                   replace its value line (`-` = no value line)
# A bare key=value that is not already in the ini goes to [FUSS-Pyrolysis].
# Guard 8 (MPI, `mpi_size_ > 1`) has no demo: it needs -m 2 (documented).
# ---------------------------------------------------------------------------
demo_pyro () {
  local name=$1 expect_rc=$2 expect=$3; shift 3
  local d="$WORK/$name"

  rm -rf "$d"; mkdir -p "$d"
  cp -r "$BASE_PYRO/INPUT" "$BASE_PYRO/MESH" "$BASE_PYRO/input.ini" "$BASE_PYRO/FUSS.sh" "$d/"
  cp "$BASE_PYRO"/*.dat "$d/" 2>/dev/null || true
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"
  sed -i 's/^iter-threshold = .*/iter-threshold = 3/' "$d/input.ini"

  local kv sec
  for kv in "$@"; do
    case "${kv%%:*}" in
      RAW)  printf '\n%s\n' "${kv#RAW:}" >> "$d/input.ini"; continue ;;
      RM)   rm -f "$d/${kv#RM:}"; continue ;;
      PYRO) python3 - "$d/INPUT/pyrolysis.dat" "${kv#PYRO:}" <<'PYEDIT'
import re, sys
p, edit = sys.argv[1], sys.argv[2]
key, val = edit.split(None, 1)
s = open(p).read()
pat = re.compile(r"^(%s)\s+.*$" % re.escape(key), re.M)
if not pat.search(s):
    sys.exit("pyrolysis.dat has no key %s" % key)
open(p, "w").write(pat.sub(lambda m: m.group(1) + " " + val, s))
PYEDIT
            continue ;;
      BCFACE) python3 - "$d/INPUT/bc.txt" "${kv#BCFACE:}" <<'PYBC'
import sys
p, spec = sys.argv[1], sys.argv[2]
face, typ, line2 = spec.split(":", 2)
face, typ = int(face), int(typ)
lines = open(p).read().split("\n")
out, i = [], 0
WITH_LINE2 = {101, 103, 301, 302, 303, 304, 305}
while i < len(lines):
    f = lines[i].split()
    if len(f) == 6 and all(t.lstrip("-").isdigit() for t in f):
        rec = [int(t) for t in f]
        has2 = rec[5] in WITH_LINE2
        if rec[4] == face:
            out.append("%8d%8d%8d%8d%8d%8d" % (rec[0], rec[1], rec[2], rec[3], rec[4], typ))
            if line2 != "-":
                out.append(" " + line2)
            i += 2 if has2 else 1
            continue
    out.append(lines[i]); i += 1
open(p, "w").write("\n".join(out))
PYBC
            continue ;;
    esac
    sec="FUSS-Pyrolysis"
    if [ "${kv:0:4}" = "SEC:" ]; then
      kv=${kv:4}; sec=${kv%%:*}; kv=${kv#*:}
    fi
    python3 - "$d/input.ini" "${kv%%=*}" "${kv#*=}" "$sec" <<'PYINI'
import re, sys
p, k, v, sec = sys.argv[1], sys.argv[2].strip(), sys.argv[3].strip(), sys.argv[4].strip()
s = open(p).read()
pat = re.compile(r"^(%s\s*=).*$" % re.escape(k), re.M)
if pat.search(s):
    s = pat.sub(lambda m: m.group(1) + " " + v, s)
else:
    hdr = "[%s]" % sec
    if hdr not in s:
        sys.exit("demo input has no section %s to receive %s" % (hdr, k))
    s = s.replace(hdr, "%s\n%s = %s" % (hdr, k, v), 1)
open(p, "w").write(s)
PYINI
  done

  ( cd "$d" && ./FUSS.sh -p 1 solve ) > "$d/run.log" 2>&1
  local rc=$?

  if ! grep -qF "$expect" "$d/run.log"; then
    bad "$name: expected \"$expect\", not found in $d/run.log"
    tail -4 "$d/run.log" | sed 's/^/          /'
    return
  fi
  if [ "$rc" -ne "$expect_rc" ]; then
    bad "$name: refused correctly but exited $rc, expected $expect_rc"
    return
  fi
  pass "$name -> exit $rc, $(grep -m1 -F "$expect" "$d/run.log" | sed 's/^ *//')"
}

demo_pyro pyro_irs          1 "pyrolysis is not compatible with implicit residual smoothing" "irs=true"
demo_pyro pyro_irs_beta     1 "pyrolysis is not compatible with implicit residual smoothing" "SEC:FUSS-Numerics:irs-beta=0.5"
demo_pyro pyro_multigrid    1 "pyrolysis is not compatible with multigrid" \
    "RAW:[FUSS-Multigrid]
levels = 2
level1-iter = 3
level2-iter = 3"
demo_pyro pyro_prim         1 "pyrolysis requires integration-variables = cons"  "integration-variables=prim"
demo_pyro pyro_type2        1 "type-2 gas closure is plan 07"                    "gas-closure=type2"
demo_pyro pyro_steady       1 "pyrolysis requires a time-accurate run"           "time-accurate=false"
demo_pyro pyro_moving       1 "pyrolysis on a moving mesh is plan 08" \
    "RAW:[FUSS-MeshMotion]
law = translation
vel = 1.0e-4 0.0 0.0"
demo_pyro pyro_no_outflow   1 "outflow-direction must be set"                    "outflow-direction=none"
# face 5 (type 0, no value line) becomes a 103 self-connection with its 9-integer line
demo_pyro pyro_connection   1 "pyrolysis phase 2 is single-block: no 101/102/103 faces" "BCFACE:5:103:1 1 1 1 6 1 0 0 1"
# the outflow face 2 (302 with a time file) becomes a 300 symmetry plane, no value line
demo_pyro pyro_outflow_face 1 "outflow face must be a wall BC (301-305)"         "BCFACE:2:300:-"
demo_pyro pyro_missing_virgin 1 "[ERROR] pyrolysis tables"                       "RM:INPUT/virgin.dat"
# F[2][1] = 0.9 makes the phase's sum of F 1.15 > 1
demo_pyro pyro_bad_F        1 "[ERROR] pyrolysis tables"                         "PYRO:F[2][1] 0.9"

echo
if [ "$fail" -eq 0 ]; then
  echo "all guards fired as intended"
else
  echo "$fail guard(s) did NOT fire -- they are decoration, not protection"
  exit 1
fi
