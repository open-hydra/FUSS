#!/bin/bash
#===============================================================================
#         FILE: test/run_guard_demos.sh
#        USAGE: ./test/run_guard_demos.sh [-w WORKDIR] [-m NRANKS]
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
#
#               -m N (plan 11) runs every demo with the mpi preset's bin/FUSS-mpi on
#               N ranks through $I_MPI_ROOT/bin/mpirun, under timeout 900: each
#               demo must keep its text AND its code there too (a stop on only
#               some of the ranks used to make the launcher report 255). A base
#               that MPI refuses as a whole -- the ablating wall -- is refused by
#               name with code 1 before the demo's own check; that outcome is
#               counted apart ("MPI") and is the first demonstration of that
#               refusal. Chimera (102) and multigrid run on N ranks since the merge
#               of upstream's chimera-mpi (2026-10-17), pyrolysis since plan 11
#               stages 3a (type 1) and 3b (type 2; 2026-10-24).
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BASE="$ROOT/test/moving-mesh/receding/N040"
BASE_PYRO="$ROOT/test/pyrolysis/energy_audit"
BASE_ABL="$ROOT/test/pyrolysis/ablation/demo_base"
BASE_ABL_PYRO="$ROOT/test/pyrolysis/ablation/demo_base_pyro"
BASE_ABL_2BLK="$ROOT/test/pyrolysis/ablation/demo_base_2blk"
BASE_R5="$ROOT/test/moving-mesh/conjugate_receding"
WORK="${TMPDIR:-/tmp}/fuss-guard-demos.$$"
NMPI=0

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -m | --mpi)  NMPI=$2; shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR] [-m NRANKS]"; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

# How every demo runs the solver: the default executable, or (-m) the mpi preset's on N ranks
# through Intel MPI's own launcher, never whatever mpirun comes first on PATH.
SOLVE="./FUSS.sh -p 1 solve"
if [ "$NMPI" -gt 0 ]; then
  [ -n "${I_MPI_ROOT:-}" ] || { echo "-m needs I_MPI_ROOT (the Intel MPI of the mpi preset)"; exit 1; }
  export FUSS_EXE=FUSS-mpi FUSS_MPIRUN="$I_MPI_ROOT/bin/mpirun"
  SOLVE="timeout 900 ./FUSS.sh -p 1 -m $NMPI solve"
fi

mkdir -p "$WORK"
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

# A helper called before it is defined, or a missing FUSS.sh, exits 127 and used to be
# skipped without a FAIL (the irs_beta demo, found 2026-10-16): count it.
set -E
trap '[ $? -eq 127 ] && bad "line $LINENO: exit 127 -- a demo helper or the solver script was not found"' ERR

# Under -m N > 1: a base that MPI refuses as a whole is stopped by name with code 1 before the demo's
# own check. That is the expected outcome there; it is counted apart, never as a PASS.
nmpi_ref=0
MPI_REFUSALS=("the ablating wall (306) is OpenMP-only until plan 11")
mpi_refused () {   # $1 = log, $2 = exit code: true when an MPI refusal stopped the run with code 1
  [ "$NMPI" -gt 1 ] && [ "$2" -eq 1 ] || return 1
  local t
  for t in "${MPI_REFUSALS[@]}"; do
    if grep -qF "$t" "$1"; then MPI_REF_TEXT=$t; return 0; fi
  done
  return 1
}
refused () { printf '  \033[0;33mMPI\033[0m   %s\n' "$1"; nmpi_ref=$((nmpi_ref+1)); }

echo "FUSS guard demonstrations"
echo "  base case : test/moving-mesh/receding/N040 (mesh motion), test/pyrolysis/energy_audit (pyrolysis),"
echo "              test/pyrolysis/ablation/demo_base{,_pyro} (the ablating wall 306)"
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

  ( cd "$d" && $SOLVE ) > "$d/run.log" 2>&1
  local rc=$?

  if ! grep -qF "$expect" "$d/run.log" || [ "$rc" -ne "$expect_rc" ]; then
    if mpi_refused "$d/run.log" "$rc"; then refused "$name -> exit 1, $MPI_REF_TEXT"; return; fi
  fi
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

  ( cd "$d" && $SOLVE ) > "$d/run.log" 2>&1
  local rc=$?
  if [ "$rc" -ne 0 ] && mpi_refused "$d/run.log" "$rc"; then refused "$name -> exit 1, $MPI_REF_TEXT"; return; fi
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

# hydra's own sections pass too: a hydra case's input.ini carries [HYDRA-Parameters] and
# [HYDRA-IO], which only hydra reads (until 2026-10-14 FUSS counted every HYDRA- section as
# its own and refused these keys, so hydra-MF stopped at setup). [HYDRA-Multigrid] stays
# FUSS's: FUSS reads its multigrid levels there under hydra, so a misspelt key is refused.
demo_ok hydra_sections "[HYDRA-Parameters]
N-coupling = 1
FUSS-iter = 1
MOSE-iter = 1
time-accurate = false
[HYDRA-IO]
coupling-diter = 10
coupling-wall = true"

demo hydra_multigrid_typo 1 "unknown option in input.ini" \
    "RAW:[HYDRA-Multigrid]
levelz = 2"

# ---------------------------------------------------------------------------
# Pyrolysis guards (plan 06 section 3.7), on test/pyrolysis/energy_audit.
# Same contract as demo(): both the text and the exit code (1 = input) are
# checked. Edit forms, in addition to demo()'s key=value / SEC: / RAW:
#   RM:<path>                       remove a file of the copy (relative path)
#   PYRO:<key> <value>              edit a key of INPUT/pyrolysis.dat
#   BCFACE:<face>:<type>:<line2>    retype every bc.txt record of that face and
#                                   replace its value line (`-` = no value line)
#   DELKEY:<key>                    delete every `key = ...` line of input.ini
#   PYROADD:<line>                  append a line to INPUT/pyrolysis.dat
#   FILE:<path>:<text>              replace a file of the copy with <text> (plan 08: env rows, B' rows)
#   BCBLK:<block>:<face>:<type>:<line2>  as BCFACE, for one block's records only
#   PY:<python code>                run in the copy's directory (an edit no other form can make)
# A bare key=value that is not already in the ini goes to [FUSS-Pyrolysis].
# Guard 8 (MPI, `mpi_size_ > 1`) has no demo: it needs -m 2 (documented).
# demo_on BASE NAME ... is the same on another base (the ablating-wall bases).
# ---------------------------------------------------------------------------
demo_pyro () { demo_on "$BASE_PYRO" "$@"; }
demo_on () {
  local base=$1 name=$2 expect_rc=$3 expect=$4; shift 4
  local d="$WORK/$name"

  rm -rf "$d"; mkdir -p "$d"
  cp -r "$base/INPUT" "$base/MESH" "$base/input.ini" "$base/FUSS.sh" "$d/"
  cp "$base"/*.dat "$d/" 2>/dev/null || true
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"
  sed -i 's/^iter-threshold = .*/iter-threshold = 3/' "$d/input.ini"

  local kv sec rest
  for kv in "$@"; do
    case "${kv%%:*}" in
      RAW)  printf '\n%s\n' "${kv#RAW:}" >> "$d/input.ini"; continue ;;
      RM)   rm -f "$d/${kv#RM:}"; continue ;;
      FILE) rest=${kv#FILE:}; printf '%s\n' "${rest#*:}" > "$d/${rest%%:*}"; continue ;;
      PY)   ( cd "$d" && python3 -c "${kv#PY:}" ) || echo "PY edit failed" >> "$d/edit.err"; continue ;;
      BCBLK) python3 - "$d/INPUT/bc.txt" "${kv#BCBLK:}" <<'PYBB'
import sys
p, spec = sys.argv[1], sys.argv[2]
blk, face, typ, line2 = spec.split(":", 3)
blk, face, typ = int(blk), int(face), int(typ)
lines = open(p).read().split("\n")
out, i = [], 0
WITH_LINE2 = {101, 103, 301, 302, 303, 304, 305, 306, 307}
while i < len(lines):
    f = lines[i].split()
    if len(f) == 6 and all(t.lstrip("-").isdigit() for t in f):
        rec = [int(t) for t in f]
        has2 = rec[5] in WITH_LINE2
        if rec[0] == blk and rec[4] == face:
            out.append("%8d%8d%8d%8d%8d%8d" % (rec[0], rec[1], rec[2], rec[3], rec[4], typ))
            if line2 != "-":
                out.append(" " + line2)
            i += 2 if has2 else 1
            continue
    out.append(lines[i]); i += 1
open(p, "w").write("\n".join(out))
PYBB
            continue ;;
      DELKEY) sed -i "/^${kv#DELKEY:}[[:space:]]*=/d" "$d/input.ini"; continue ;;
      PYROADD) printf '%s\n' "${kv#PYROADD:}" >> "$d/INPUT/pyrolysis.dat"; continue ;;
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
WITH_LINE2 = {101, 103, 301, 302, 303, 304, 305, 306, 307}
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

  ( cd "$d" && $SOLVE ) > "$d/run.log" 2>&1
  local rc=$?

  if ! grep -qF "$expect" "$d/run.log" || [ "$rc" -ne "$expect_rc" ]; then
    if mpi_refused "$d/run.log" "$rc"; then refused "$name -> exit 1, $MPI_REF_TEXT"; return; fi
  fi
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

# `irs-beta > 0` alone also enables IRS (Assign_Setup promotes
# it), but that promotion used to run AFTER Check_Input, so this exact input
# passed the `irs` guard (demo irs, above) and smoothed the ALE-remapped residual anyway.
# The base says `irs = false`, which beside irs-beta > 0 is refused as a contradiction
# (irs_contradiction below), so this demo deletes it: what it tests is the promotion.
demo_on "$BASE" irs_beta 1 "not compatible with implicit residual smoothing" \
    "DELKEY:irs" "SEC:FUSS-Numerics:irs-beta=0.5"
# (Until 2026-10-16 this call stood above demo_on's definition: bash ran nothing, exit 127, and
# no FAIL was counted.)
demo_pyro pyro_irs          1 "pyrolysis is not compatible with implicit residual smoothing" "irs=true"
demo_pyro pyro_irs_beta     1 "pyrolysis is not compatible with implicit residual smoothing" "DELKEY:irs" "SEC:FUSS-Numerics:irs-beta=0.5"
demo_pyro pyro_multigrid    1 "pyrolysis is not compatible with multigrid" \
    "RAW:[FUSS-Multigrid]
levels = 2
level1-iter = 3
level2-iter = 3"
demo_pyro pyro_prim         1 "pyrolysis requires integration-variables = cons"  "integration-variables=prim"
demo_pyro pyro_steady       1 "pyrolysis requires a time-accurate run"           "time-accurate=false"
demo_pyro pyro_moving       1 "pyrolysis on a moving mesh needs law = nodeshift" \
    "RAW:[FUSS-MeshMotion]
law = translation
vel = 1.0e-4 0.0 0.0"
demo_pyro pyro_no_outflow   1 "outflow-direction must be set"                    "outflow-direction=none"
# face 5 (type 0, no value line) becomes a 103 self-connection with its 9-integer line
demo_pyro pyro_connection   1 "pyrolysis: 102, 103 and 104 faces are refused (a 101 face needs gas-closure = type2)" "BCFACE:5:103:1 1 1 1 6 1 0 0 1"
# a block connection under type 1: its outflow line integral cannot cross it (plan 08 section 2.9 (4))
demo_pyro pyro_type1_multiblock 1 "gas-closure = type1 is single-block" "BCFACE:5:101:1 1 1 1 6 1 0 0 1"
demo_pyro pyro_chimera_104  1 "pyrolysis: 102, 103 and 104 faces are refused (a 101 face needs gas-closure = type2)" "BCFACE:5:104:1 1
 1 1 1 1 1.0
 1 1 1 1 1.0"
# the outflow face 2 (302 with a time file) becomes a 300 symmetry plane, no value line
demo_pyro pyro_outflow_face 1 "outflow face must be a wall BC (301-307)"         "BCFACE:2:300:-"
demo_pyro pyro_missing_virgin 1 "[ERROR] pyrolysis tables"                       "RM:INPUT/virgin.dat"
# F[2][1] = 0.9 makes the phase's sum of F 1.15 > 1
demo_pyro pyro_bad_F        1 "[ERROR] pyrolysis tables"                         "PYRO:F[2][1] 0.9"

# ---------------------------------------------------------------------------
# Darcy closure guards (plan 07 section 3.8), on the same pyrolysis base with
# gas-closure = type2. The base sets p-gas, which type 2 refuses (guard 13), so
# every other type-2 demo deletes it first; the bc-level demos also add the
# permeability tensors the table check requires, so that each demo stops for
# its own reason. darcy_open_no_gas is deferred to the first two-material case
# (plan 07 section 3.8: one material cannot have a cell without gas zones).
# ---------------------------------------------------------------------------
DARCY_T2="gas-closure=type2"
DARCY_K1="PYROADD:K_v 1e-11 0 0 0 1e-11 0 0 0 1e-11"
DARCY_K2="PYROADD:K_c 2e-11 0 0 0 2e-11 0 0 0 2e-11"
# guard 12: a [FUSS-Darcy] section that nothing would read -- without pyrolysis, and under type 1
demo darcy_no_pyro 1 "[FUSS-Darcy] is read only with [FUSS-Pyrolysis] gas-closure = type2" \
    "RAW:[FUSS-Darcy]
solver = tridiag"
demo_pyro darcy_type1        1 "[FUSS-Darcy] is read only with [FUSS-Pyrolysis] gas-closure = type2" \
    "RAW:[FUSS-Darcy]
solver = tridiag"
# guard 13: the type-1 pressure key under type 2
demo_pyro darcy_p_gas_key    1 "p-gas is refused with gas-closure = type2" "$DARCY_T2"
# guards 14-16: the hand-checked [FUSS-Darcy] ranges and lists
demo_pyro darcy_bad_vnn      1 "pseudo-vnn must be in (0, 1]" "$DARCY_T2" "DELKEY:p-gas" \
    "RAW:[FUSS-Darcy]
pseudo-vnn = 1.5"
demo_pyro darcy_bad_open_type 1 "open-bc-types: bad token \"300\"" "$DARCY_T2" "DELKEY:p-gas" \
    "RAW:[FUSS-Darcy]
open-bc-types = 300"
demo_pyro darcy_axes_follow_grid 1 "axes-follow-grid = true is not built" "$DARCY_T2" "DELKEY:p-gas" \
    "RAW:[FUSS-Darcy]
axes-follow-grid = true"
# the table check: type 2 needs the permeability tensors
demo_pyro darcy_missing_K    1 "key K_v, K_c: required by gas-closure = type2" "$DARCY_T2" "DELKEY:p-gas"
# bc-level (Check_Pyrolysis_BC_Topology): the base's faces are 301/302/300/0, none open by default
demo_pyro darcy_sealed_with_source 1 "no face is open to the gas while the pyrolysis gas is produced" \
    "$DARCY_T2" "DELKEY:p-gas" "$DARCY_K1" "$DARCY_K2"
demo_pyro darcy_sealed_quasi_steady 1 "a sealed domain has no unique quasi-steady pressure" \
    "$DARCY_T2" "DELKEY:p-gas" "$DARCY_K1" "$DARCY_K2" \
    "RAW:[FUSS-Darcy]
quasi-steady = true"
# face 5 (type 0) becomes a 102 chimera face with one donor per ghost row, volume fraction 1, so
# that IO_BC's no-donor abort (104: upstream 6ad97e0, 102 too since chimera-mpi's 4634cd0) does not
# stop the run first. The Darcy text precedes plan 06 guard 9's (which refuses 101/102/103/104 on
# every pyrolysis run until plan 08)
demo_pyro darcy_chimera      1 "a chimera (102, 104) face has no pressure condition" \
    "$DARCY_T2" "DELKEY:p-gas" "$DARCY_K1" "$DARCY_K2" "BCFACE:5:102:1 1
 1 1 1 1 1.0
 1 1 1 1 1.0"
# the same with a 104 face (chimera and coupled wall, upstream 6ad97e0)
demo_pyro darcy_chimera_104  1 "a chimera (102, 104) face has no pressure condition" \
    "$DARCY_T2" "DELKEY:p-gas" "$DARCY_K1" "$DARCY_K2" "BCFACE:5:104:1 1
 1 1 1 1 1.0
 1 1 1 1 1.0"
# the per-face p_w of a 303 line (plan 07 section 4.3) must be positive
demo_pyro darcy_bad_pw       1 "the fourth value of a 303 line (p_w) must be > 0" \
    "$DARCY_T2" "DELKEY:p-gas" "$DARCY_K1" "$DARCY_K2" "BCFACE:1:303:0.0, 0.0, 300.0, -1.0"
# probe variables: p_pg exists only under type 2, and an unknown name is refused (it used to be ignored)
demo_pyro probe_p_pg_type1   1 "probe variable p_pg needs [FUSS-Pyrolysis] gas-closure = type2" "variables=T p_pg"
demo_pyro probe_unknown_var  1 "unknown probe variable \"q\"" "variables=T q"

# ---------------------------------------------------------------------------
# Ablating-wall guards (plan 08 stage A: the bc.txt line, the env and B' tables, the setup rules, the restart
# state, the solve), on test/pyrolysis/ablation/demo_base (kinematic 306, legacy material) and demo_base_pyro
# (type 2, the toy B' table). abl_base_ok and abl_base_pyro_ok first prove that both bases run as committed, so
# that every demo below stops for its own reason. The MPI refusal ("OpenMP-only until plan 11") needs -m 2 and
# has no demo, as guard 8.
# ---------------------------------------------------------------------------
ABL_ENV_HDR="// t p_w rhoeUeCH0 h_r q_rad chem hconv T_edge"
demo_ok_on () {   # $1 = base, $2 = name: the base must run to completion (exit 0, no [ERROR])
  local base=$1 name=$2 d="$WORK/$2"
  rm -rf "$d"; mkdir -p "$d"
  cp -r "$base/INPUT" "$base/MESH" "$base/input.ini" "$base/FUSS.sh" "$d/"
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"
  ( cd "$d" && $SOLVE ) > "$d/run.log" 2>&1
  local rc=$?
  if [ "$rc" -ne 0 ] && mpi_refused "$d/run.log" "$rc"; then refused "$name -> exit 1, $MPI_REF_TEXT"; return; fi
  if [ "$rc" -ne 0 ] || grep -q '\[ERROR\]' "$d/run.log"; then
    bad "$name: the base itself does not run (exit $rc, see $d/run.log)"
    return
  fi
  pass "$name -> exit 0, the base runs as committed"
}
demo_ok_on "$BASE_ABL"      abl_base_ok
demo_ok_on "$BASE_ABL_PYRO" abl_base_pyro_ok

# bc.txt line 2 of a 306 record
demo_on "$BASE_ABL" abl_bad_line    1 "a 306 line is \"lambda  T_inf  <bprime-file>  <env-file>\"" "BCFACE:2:306:0.5, 300.0"
demo_on "$BASE_ABL" abl_bad_lambda  1 "a 306 line needs lambda >= 0 and T_inf > 0" "BCFACE:2:306:-1.0, 300.0, none, env_demo.dat"
demo_on "$BASE_ABL" abl_bad_tinf    1 "a 306 line needs lambda >= 0 and T_inf > 0" "BCFACE:2:306:0.0, 0.0, none, env_demo.dat"
# the env table (time_table_t%load)
demo_on "$BASE_ABL" abl_env_missing   1 "env table INPUT/env_demo.dat: cannot open" "RM:INPUT/env_demo.dat"
demo_on "$BASE_ABL" abl_env_one_row   1 "env table INPUT/env_demo.dat: fewer than two rows" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.0 0.0 0.0 0.0 10.0 400.0"
demo_on "$BASE_ABL" abl_env_order     1 "t must increase strictly from row to row" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.0 0.0 0.0 0.0 10.0 400.0
0.0 101325.0 0.0 0.0 0.0 0.0 10.0 400.0"
demo_on "$BASE_ABL" abl_env_short_row 1 "line 2: expected 8 numbers (t and 7 columns)" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.0 0.0 0.0 0.0 10.0
1.0 101325.0 0.0 0.0 0.0 0.0 10.0"
demo_on "$BASE_ABL" abl_env_not_finite 1 "value is not finite" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.0 0.0 0.0 0.0 NaN 400.0
1.0 101325.0 0.0 0.0 0.0 0.0 10.0 400.0"
# the env rows (Check_Ablation), one demo per condition of the one rule
for spec in "pw:0.0 0.0 0.0 0.0 0.0 0.0 10.0 400.0" "rhoe:0.0 101325.0 -1.0 0.0 0.0 0.0 10.0 400.0" \
            "hconv:0.0 101325.0 0.0 0.0 0.0 0.0 -10.0 400.0" "chem_neg:0.0 101325.0 0.0 0.0 0.0 -0.5 10.0 400.0" \
            "chem_gt1:0.0 101325.0 0.0 0.0 0.0 1.5 10.0 400.0" "tedge:0.0 101325.0 0.0 0.0 0.0 0.0 10.0 0.0"; do
  demo_on "$BASE_ABL" "abl_env_row_${spec%%:*}" 1 "needs p_w > 0, rhoeUeCH0 >= 0, hconv >= 0, 0 <= chem <= 1 and T_edge > 0" \
      "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
${spec#*:}
1.0e6 101325.0 0.0 0.0 0.0 0.0 10.0 400.0"
done
demo_on "$BASE_ABL" abl_chem_without_pyro 1 "without [FUSS-Pyrolysis] (kinematic mode) has no pyrolysis gas and no B' table" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.1 1.0e6 0.0 1.0 0.0 300.0
1.0e6 101325.0 0.1 1.0e6 0.0 1.0 0.0 300.0"
# setup rules
# multigrid needs an even cell count across, which the 1-D demo base lacks (Check_Multigrid would stop first): the
# 40 x 4 receding base, its mesh motion removed and its face 1 (302) made a kinematic 306 face
demo_on "$BASE" abl_multigrid   1 "the ablating wall (306) is not compatible with multigrid" \
    "DELKEY:law" "DELKEY:vel" "DELKEY:gcl-tolerance" "BCFACE:1:306:0.0, 300.0, none, env_demo.dat" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.0 0.0 0.0 0.0 10.0 400.0
1.0e6 101325.0 0.0 0.0 0.0 0.0 10.0 400.0" \
    "RAW:[FUSS-Multigrid]
levels = 2
level1-iter = 3
level2-iter = 3"
demo_on "$BASE_ABL" abl_306_translation 1 "the ablating wall (306) cannot ride a translation mesh" \
    "law=translation" "SEC:FUSS-MeshMotion:vel=1.0e-5 0.0 0.0"
# [FUSS-Ablation] registry rules
demo_on "$BASE_ABL" abl_bad_seb_tol      1 "seb-tol must be > 0"      "SEC:FUSS-Ablation:seb-tol=0.0"
demo_on "$BASE_ABL" abl_bad_seb_max_iter 1 "seb-max-iter must be > 0" "SEC:FUSS-Ablation:seb-max-iter=0"
demo_on "$BASE_ABL" abl_bad_oor          1 "bprime-out-of-range must be one of" "SEC:FUSS-Ablation:bprime-out-of-range=maybe"
# the type-1 line integral ends at the outflow face: a 306 face elsewhere is refused (energy_audit: outflow +i, face 2)
demo_pyro abl_type1_not_outflow 1 "under gas-closure = type1 a 306 face must be the outflow face" \
    "BCFACE:1:306:0.0, 300.0, none, env_demo.dat" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 0.0 0.0 0.0 0.0 10.0 400.0
1.0e6 101325.0 0.0 0.0 0.0 0.0 10.0 400.0"
# the B' table (FUSS_Lib_Bprime) on the pyrolysis base
ABL_B_HDR="// kind: equilibrium
// axes: p Bprime-g T"
demo_on "$BASE_ABL_PYRO" abl_bprime_missing 1 "B' table INPUT/bprime_demo.dat: cannot open" "RM:INPUT/bprime_demo.dat"
demo_on "$BASE_ABL_PYRO" abl_bprime_order   1 "p must decrease from block to block" \
    "FILE:INPUT/bprime_demo.dat:$ABL_B_HDR
50000.0 10.0 0.1 300.0 -700000.0
50000.0 10.0 0.1 4000.0 3000000.0
50000.0 0.001 0.1 300.0 -700000.0
50000.0 0.001 0.1 4000.0 3000000.0
200000.0 10.0 0.1 300.0 -700000.0
200000.0 10.0 0.1 4000.0 3000000.0
200000.0 0.001 0.1 300.0 -700000.0
200000.0 0.001 0.1 4000.0 3000000.0"
demo_on "$BASE_ABL_PYRO" abl_bprime_bg_order 1 "B'g must decrease from sub-block to sub-block" \
    "FILE:INPUT/bprime_demo.dat:$ABL_B_HDR
200000.0 0.001 0.1 300.0 -700000.0
200000.0 0.001 0.1 4000.0 3000000.0
200000.0 10.0 0.1 300.0 -700000.0
200000.0 10.0 0.1 4000.0 3000000.0"
demo_on "$BASE_ABL_PYRO" abl_bprime_T_order 1 "T must not decrease inside an increasing (p, B'g) sub-table" \
    "FILE:INPUT/bprime_demo.dat:$ABL_B_HDR
200000.0 10.0 0.1 300.0 -700000.0
200000.0 10.0 0.1 4000.0 3000000.0
200000.0 10.0 0.1 3000.0 2000000.0
200000.0 0.001 0.1 300.0 -700000.0
200000.0 0.001 0.1 4000.0 3000000.0"
demo_on "$BASE_ABL_PYRO" abl_bprime_kind 1 "unknown kind \"magic\"" \
    "FILE:INPUT/bprime_demo.dat:// kind: magic
200000.0 10.0 0.1 300.0 -700000.0
200000.0 10.0 0.1 4000.0 3000000.0"
demo_on "$BASE_ABL_PYRO" abl_bprime_one_row 1 "fewer than two data rows" \
    "FILE:INPUT/bprime_demo.dat:$ABL_B_HDR
200000.0 10.0 0.1 300.0 -700000.0"
demo_on "$BASE_ABL_PYRO" abl_bprime_none_chem 1 "has chem > 0 and the face has no B' table (none)" \
    "BCFACE:2:306:0.5, 300.0, none, env_demo.dat"
demo_on "$BASE_ABL_PYRO" abl_pw_outside_gas 1 "lies outside the gas table of the face's material" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 2.0e6 0.1 2.0e6 0.0 1.0 0.0 300.0
1.0e6 2.0e6 0.1 2.0e6 0.0 1.0 0.0 300.0"
# exit 3: the solve. A p_w inside the gas table but outside the B' table's p range clamps the lookup;
# bprime-out-of-range = abort stops there. rhoeUeCH0 1e9 with two Newton iterations cannot converge.
demo_on "$BASE_ABL_PYRO" abl_bprime_abort 3 "the B' table query left the table (bprime-out-of-range = abort)" \
    "SEC:FUSS-Ablation:bprime-out-of-range=abort" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 3.0e5 0.1 2.0e6 0.0 1.0 0.0 300.0
1.0e6 3.0e5 0.1 2.0e6 0.0 1.0 0.0 300.0"
demo_on "$BASE_ABL_PYRO" seb_fail 3 "SEB did not converge" "SEC:FUSS-Ablation:seb-max-iter=2" \
    "FILE:INPUT/env_demo.dat:$ABL_ENV_HDR
0.0 101325.0 1.0e9 2.0e6 0.0 1.0 0.0 300.0
1.0e6 101325.0 1.0e9 2.0e6 0.0 1.0 0.0 300.0"

# ---------------------------------------------------------------------------
# Receding mesh (plan 08 stage B, law = nodeshift): the bc-level rules (Check_Ablation_BC_Topology), the node-map
# rules (Setup_Recession: (a) a closed 306 patch, (b) equal receding lines across a 101 interface), and the stops of
# a moving surface (exit 2). demo_base recedes at sdot-override 1e-5 m/s; demo_base_2blk is two blocks, the second
# taller; conjugate_receding (R5) is two blocks that share their 101 face.
# ---------------------------------------------------------------------------
ABL_306_LINE="0.0, 300.0, none, env_demo.dat"
demo_ok_on "$BASE_ABL_2BLK" abl_base_2blk_ok
demo_on "$BASE_ABL" abl_two_faces      1 "carries 306 faces on two face ids" "BCFACE:3:306:$ABL_306_LINE"
demo_on "$BASE_ABL" abl_opposite_306   1 "are opposite faces" "BCFACE:1:306:$ABL_306_LINE"
demo_on "$BASE_ABL" abl_nodeshift_no306 1 "law = nodeshift moves the ablating faces (306), and bc.txt has none" \
    "BCFACE:2:302:    3.000000E+02,"
demo_on "$BASE_ABL" abl_receding_coupled 1 "recedes (law = nodeshift) and owns a coupled face 103" \
    "BCFACE:1:103:       1       1       1       1       1       1       0       0       1"
# a second material in the receding block: phase.txt, a copy of the properties zone, the last cell's matID
demo_on "$BASE_ABL" abl_multimat       1 "holds more than one material" \
    "PY:p = 'INPUT/phase.txt'; open(p, 'a').write('mat2 1\n'); s = open('INPUT/properties.dat').read().split('\n'); open('INPUT/properties.dat', 'a').write('\n'.join(['ZONE T=\"mat2\"'] + s[3:]) + '\n'); L = open('INPUT/ic.tec').read().rstrip('\n').split('\n'); L[-1] = ' 2.000000000000000E+00'; open('INPUT/ic.tec', 'w').write('\n'.join(L) + '\n')"
# a volumetric source in the receding block: INPUT/st.tec from ic.tec's nodes, qvol 1e3 W/m3 in every cell
demo_on "$BASE_ABL" abl_qvol           1 "has a volumetric source qvol" \
    "PY:L = open('INPUT/ic.tec').read().split('\n'); nn = 21 * 2 * 2; out = [' VARIABLES =\"x\" \"y\" \"z\" qvol', L[1].replace('[4-5]=CELLCENTERED', '[4]=CELLCENTERED')] + L[2:2 + 3 * nn] + [' 1.000000000000000E+03'] * 20; open('INPUT/st.tec', 'w').write('\n'.join(out) + '\n')"
# the node map on two blocks: the bottom of block 1 alone ablating leaves the shared corner in a block without 306;
# both bottoms ablating make the corner end lines of 0.02 and 0.04 m
demo_on "$BASE_ABL_2BLK" abl_open_patch     1 "has no 306 face containing it" \
    "BCBLK:1:1:301:    0.000000E+00," "BCBLK:1:3:306:$ABL_306_LINE"
demo_on "$BASE_ABL_2BLK" abl_lline_mismatch 1 "on lines of different length" \
    "BCBLK:1:1:301:    0.000000E+00," "BCFACE:3:306:$ABL_306_LINE"
# exit 2, the moving surface. 5.3e-3 m/s moves the 2.5 mm surface cell of demo_base by twice its size in the first
# 0.9375 s step: the remap's sweep ratio stops it while 45 mm of the 50 mm line remain (min-thickness 25 mm).
demo_on "$BASE_ABL" recession_too_fast 2 "the mesh moved too far in one time step" "sdot-override=5.3e-3"
# min-thickness = the line itself: consumed at the first step, whatever the rate
demo_on "$BASE_ABL" block_consumed     2 "ablating block consumed" "SEC:FUSS-Ablation:min-thickness=0.05"
# R5 with the TEST ONLY perturbation of one shared node at step 2: the interface assertion stops it
demo_on "$BASE_R5" abl_node_mismatch   2 "the mesh tore at a block interface" "SEC:FUSS-Ablation:perturb-interface-node=true"
# Pyrolysis on the receding mesh (plan 08 stage C): the remap's fixed-size work vectors hold 64 sub-phases (reactions);
# 65 are refused at the first remap (64 run: checked 2026-10-13). The conservation audit (error stop 2) cannot fire on a
# correct remap; its demonstration is a scratch build (plan-bucket/records/2026-10-13-plan08-stageC/).
REAC65="s = open('INPUT/pyrolysis.dat').read(); s = s.replace('nPyroReac[2] 0', 'nPyroReac[2] 65\\n' + ''.join('F[2][%d] 0.005\\nA[2][%d] 12000\\nE[2][%d] 71130.89\\nm[2][%d] 3\\nn[2][%d] 0\\nT[2][%d] 333.3\\nh[2][%d] -4e6\\n' % ((r,) * 7) for r in range(1, 66))); open('INPUT/pyrolysis.dat', 'w').write(s)"
demo_on "$BASE_ABL_PYRO" abl_remap_subphases 1 "pyrolysis remap: more than 64 sub-phases" \
    "RAW:[FUSS-MeshMotion]
law = nodeshift" "PY:$REAC65"
# The coupled ablating wall 307 (plan 08 stage D, section 2.5 ii): line 2 = the 103 connectivity record, T_inf and the
# B' table, no lambda; its exchange values come from hydra or, standalone, from the TEST ONLY coupled-env-file.
L307="1 1 1 1 1 1 0 0 1 300.0 bprime_demo.dat"
CENV="0.0 101325.0 0.1 2000000.0 0.0
1000000.0 101325.0 0.1 2000000.0 0.0"
demo_on "$BASE_ABL_PYRO" abl307_lambda       1 "(11 values, 12 read): a coupled face takes no lambda" \
    "BCFACE:2:307:1 1 1 1 1 1 0 0 1 0.5 300.0 bprime_demo.dat"
demo_on "$BASE_ABL_PYRO" abl307_no_table     1 "a 307 line needs T_inf > 0 and a B' table" "BCFACE:2:307:1 1 1 1 1 1 0 0 1 300.0 none"
demo_on "$BASE_ABL_PYRO" abl307_no_exchange  1 "no exchange values (p_w = 0) at this step" "BCFACE:2:307:$L307"
demo_on "$BASE_ABL_PYRO" abl307_env_no307    1 "coupled-env-file is set, but bc.txt has no coupled ablating wall (307)" \
    "FILE:INPUT/cenv.dat:$CENV" "SEC:FUSS-Ablation:coupled-env-file=cenv.dat"
demo_on "$BASE_ABL_PYRO" abl307_env_bad_row  1 "needs p_w > 0, rhoeUeCH >= 0 and q_rad >= 0" "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:0.0 0.0 0.1 2000000.0 0.0
1000000.0 101325.0 0.1 2000000.0 0.0" "SEC:FUSS-Ablation:coupled-env-file=cenv.dat"
demo_on "$BASE_ABL" abl307_without_pyro      1 "a coupled ablating wall (307) runs with the chemistry on, which needs [FUSS-Pyrolysis]" \
    "PY:import shutil; shutil.copy('$BASE_ABL_PYRO/INPUT/bprime_demo.dat', 'INPUT/bprime_demo.dat')" \
    "BCFACE:2:307:$L307" "SEC:FUSS-MeshMotion:law=static"
demo_on "$BASE_ABL_PYRO" abl307_receding     1 "recedes (law = nodeshift) and owns a coupled face 307" \
    "BCFACE:1:307:$L307" "FILE:INPUT/cenv.dat:$CENV" "SEC:FUSS-Ablation:coupled-env-file=cenv.dat" \
    "RAW:[FUSS-MeshMotion]
law = nodeshift"
# Plan 12: coupled-env-file split into one block per 307 face by `// face b i j k f` lines (hydra's env-dump replayed).
# The base's 307 face is (1 20 1 1 2); FBLK is its block, opened by the setup row hydra writes (no values yet).
FBLK="// face 1 20 1 1 2
0.0 0.0 0.0 0.0 0.0
1e-12 101325.0 0.1 2000000.0 0.0
1000000.0 101325.0 0.1 2000000.0 0.0"
CENV_KEY="SEC:FUSS-Ablation:coupled-env-file=cenv.dat"
demo_on "$BASE_ABL_PYRO" abl307_faces_header  1 'a block header is "// face b i j k f": five integers' "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:// face 1 20 1 1
$CENV" "$CENV_KEY"
demo_on "$BASE_ABL_PYRO" abl307_faces_orphan  1 'a row before the first "// face" line' "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:$CENV
$FBLK" "$CENV_KEY"
demo_on "$BASE_ABL_PYRO" abl307_faces_no307   1 'names no coupled ablating wall (307) of bc.txt' "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:$FBLK
// face 1 1 1 1 1
$CENV" "$CENV_KEY"
demo_on "$BASE_ABL_PYRO" abl307_faces_twice   1 'a second block for the 307 face "// face 1 20 1 1 2"' "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:$FBLK
$FBLK" "$CENV_KEY"
demo_on "$BASE_ABL_PYRO" abl307_faces_missing 1 'no "// face" block for the 307 face of bc.txt record 1 ("// face 1 1 1 1 1")' \
    "BCFACE:1:307:$L307" "BCFACE:2:307:$L307" "FILE:INPUT/cenv.dat:$FBLK" "$CENV_KEY"
demo_on "$BASE_ABL_PYRO" abl307_faces_short   1 'has fewer than two rows' "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:// face 1 20 1 1 2
0.0 101325.0 0.1 2000000.0 0.0" "$CENV_KEY"
demo_on "$BASE_ABL_PYRO" abl307_faces_bad_row 1 'p_w = 0 only on the first row of a block' "BCFACE:2:307:$L307" \
    "FILE:INPUT/cenv.dat:// face 1 20 1 1 2
0.0 101325.0 0.1 2000000.0 0.0
1000000.0 0.0 0.1 2000000.0 0.0" "$CENV_KEY"
# hydra's setup row, then no values until t = 1e6 s: the first step has none
demo_on "$BASE_ABL_PYRO" abl307_faces_early   1 "precedes the first exchange row of the face's coupled-env-file block" \
    "BCFACE:2:307:$L307" "FILE:INPUT/cenv.dat:// face 1 20 1 1 2
0.0 0.0 0.0 0.0 0.0
1000000.0 101325.0 0.1 2000000.0 0.0" "$CENV_KEY"
# the wall smoothing's stencil is centred on its node (plan 08 section 2.10)
demo_on "$BASE_ABL" abl_bad_width            1 "smoothing-width = 4 must be odd" "SEC:FUSS-Ablation:smoothing-width=4"

# the restart state (Setup_Ablation): a 5-step run of the kinematic base, OUTPUT/bc306-state.dat edited, newrun = false
demo_abl_restart () {   # $1 = name, $2 = expected text, $3 = python edit of the state file (or RM)
  local name=$1 expect=$2 edit=$3 d="$WORK/$1"
  rm -rf "$d"; mkdir -p "$d"
  cp -r "$BASE_ABL/INPUT" "$BASE_ABL/MESH" "$BASE_ABL/input.ini" "$BASE_ABL/FUSS.sh" "$d/"
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"
  ( cd "$d" && $SOLVE ) > "$d/run1.log" 2>&1
  local rc1=$?
  if [ "$rc1" -ne 0 ]; then
    if mpi_refused "$d/run1.log" "$rc1"; then refused "$name -> exit 1, $MPI_REF_TEXT"; return; fi
    bad "$name: the first leg failed (see $d/run1.log)"; return
  fi
  if [ "$edit" = "RM" ]; then
    rm -f "$d/OUTPUT/bc306-state.dat"
  else
    python3 -c "$edit" "$d/OUTPUT/bc306-state.dat" || { bad "$name: the state-file edit failed"; return; }
  fi
  sed -i 's/^newrun = .*/newrun = false/' "$d/input.ini"
  ( cd "$d" && $SOLVE ) > "$d/run.log" 2>&1
  local rc=$?
  if ! grep -qF "$expect" "$d/run.log" || [ "$rc" -ne 1 ]; then
    if mpi_refused "$d/run.log" "$rc"; then refused "$name -> exit 1, $MPI_REF_TEXT"; return; fi
  fi
  if ! grep -qF "$expect" "$d/run.log"; then
    bad "$name: expected \"$expect\", not found in $d/run.log"; tail -4 "$d/run.log" | sed 's/^/          /'; return
  fi
  if [ "$rc" -ne 1 ]; then bad "$name: refused correctly but exited $rc, expected 1"; return; fi
  pass "$name -> exit $rc, $(grep -m1 -F "$expect" "$d/run.log" | sed 's/^ *//')"
}
ABL_ROWS='import sys; p = sys.argv[1]; L = open(p).read().splitlines(); h = [l for l in L if l.startswith("#")]; r = [l for l in L if not l.startswith("#")]'
demo_abl_restart abl_restart_no_state "OUTPUT/bc306-state.dat is missing" RM
demo_abl_restart abl_restart_bad_row  "expected \"b i j k f s_total Tw_prev sdot beta Bpg_lag\"" \
    "$ABL_ROWS; open(p, 'w').write('\n'.join(h + [r[0].split()[0] + ' x'] + r[1:]) + '\n')"
demo_abl_restart abl_restart_foreign_row "no 306 or 307 face of bc.txt matches this row" \
    "$ABL_ROWS; t = r[0].split(); t[1] = '7'; open(p, 'w').write('\n'.join(h + [' '.join(t)] + r[1:]) + '\n')"
demo_abl_restart abl_restart_dup_row  "a second row for the same face" \
    "$ABL_ROWS; open(p, 'w').write('\n'.join(h + r + r) + '\n')"
demo_abl_restart abl_restart_no_row   "no row for the 306 face of bc.txt record" \
    "$ABL_ROWS; open(p, 'w').write('\n'.join(h) + '\n')"

# ---------------------------------------------------------------------------
# Plan 10 (2026-10-15): latent defects closed. Same contract (text and exit code).
# ---------------------------------------------------------------------------
BASE_RAD="$ROOT/test/steady-state/radiative_Rod"
BASE_PROBE="$ROOT/test/transient/oscillating_Rod"

# D12: a wall radiating far above Tref (Tref = 1 K) behind a weak solid link
# (k = 0.005: G = 2k/dx = 10 W/m2/K). The wall-temperature fixed point has slope
# -1.39 at the first call, the iterates oscillate with growing amplitude and the
# cap of 100 is reached. Until 2026-10-15 both loops (304, 305) warned and the run
# continued on the unconverged flux, exit 0.
RAD_K='p = "INPUT/properties.dat"; L = open(p).read().split("\n")
for n, l in enumerate(L):
    f = l.split()
    if len(f) == 5 and f[0][0].isdigit(): f[3] = "0.005"; L[n] = " ".join(f)
open(p, "w").write("\n".join(L))'
demo_on "$BASE_RAD" radiation_noconv     3 "BC_Wall_Radiation: no convergence in 100 iterations" \
    "BCFACE:2:304:   0.980000E+00,    0.100000E+01," "PY:$RAD_K" "SEC:FUSS-Parameters:iter-threshold=3"
demo_on "$BASE_RAD" radiation_noconv_305 3 "BC_Wall_RadiationConvection: no convergence in 100 iterations" \
    "BCFACE:2:305:   0.000000E+00,    0.980000E+00,    0.100000E+01," "PY:$RAD_K" "SEC:FUSS-Parameters:iter-threshold=3"

# D13 follow-up: an explicit irs = false beside irs-beta > 0 (the base says irs = false).
demo irs_contradiction 1 "irs = false contradicts irs-beta > 0" "SEC:FUSS-Numerics:irs-beta=0.5"

# D5: without TecIO (USE_TECIO = false, this build) ORION stops at the first binary
# Tecplot write or read with a bare stop, exit 0, after the setup.
demo binary_without_tecio    1 "sol-format = tecplot binary needs TecIO" "sol-format=tecplot binary"
demo binary_ic_without_tecio 1 "ic-format = tecplot binary needs TecIO"  "SEC:FUSS-IO:ic-format=tecplot binary"

# D18: a string longer than its variable was cut without a message. A 17-character
# probe name became 16 characters, its own section was looked up under the cut name,
# and the run stopped later on "unknown probe variable". Now refused by name.
PROBE_NAME='s = open("input.ini").read(); s = s.replace("probe1 = p1", "probe1 = NAME").replace("[p1]", "[NAME]")
s = s.replace("time-threshold = 160.0", "time-threshold = 1.0"); open("input.ini", "w").write(s)'
demo_on "$BASE_PROBE" probe_name_long 1 "[FUSS-Probes] probe1 has 17 characters, the limit is 16" \
    "PY:${PROBE_NAME//NAME/abcdefghijklmnopq}"
# ... and up to 16 characters are written in full (OUTPUT/<name>.txt; before, the
# path itself was cut at 16: `d01_lo` wrote OUTPUT/d01_lo.tx), with any number of
# variables (a ninth used to write past an 8-element array): 9 x T -> 10 columns.
PROBE_NINE='import re; s = open("input.ini").read()
s = re.sub(r"(\[abcdefghijklmnop\][^\[]*?variables = )T", r"\1T T T T T T T T T", s, flags=re.S); open("input.ini", "w").write(s)'
demo_on "$BASE_PROBE" probe_name_16 0 "Elapsed" "PY:${PROBE_NAME//NAME/abcdefghijklmnop}" "PY:$PROBE_NINE"
f="$WORK/probe_name_16/OUTPUT/abcdefghijklmnop.txt"
if [ -s "$f" ] && [ "$(awk 'NR == 1 {print NF}' "$f")" = "10" ]; then
  pass "probe_name_16 -> OUTPUT/abcdefghijklmnop.txt written in full, 9 variables + time"
else
  bad "probe_name_16: OUTPUT/abcdefghijklmnop.txt missing or not 10 columns ($(ls "$WORK/probe_name_16/OUTPUT" | tr '\n' ' '))"
fi

# D19's sibling: a 301/302 time file in a steady run, whose time stays -1 (before the
# file's first row: the series returned an undefined value there). Refused unless a
# clock runs (FUSS or hydra time-accurate).
demo_on "$BASE_PROBE" tfile_steady 1 "needs a time-accurate run" "time-accurate=false"

# D23: a bc.txt record that names a face or block that does not exist indexed n_bf and
# the blocks out of bounds; an unreadable record printed a line and went on.
BC_FACE7='p = "INPUT/bc.txt"; L = open(p).read().split("\n"); f = L[0].split(); f[4] = "7"
L[0] = "".join("%8s" % t for t in f); open(p, "w").write("\n".join(L))'
demo_on "$BASE" bc_face_out_of_range 1 "bc.txt record 1: block 1, face 7" "PY:$BC_FACE7"
BC_CONN9='p = "INPUT/bc.txt"; L = open(p).read().split("\n")
k = next(n for n, l in enumerate(L) if l.split()[-1:] == ["101"]); f = L[k + 1].split(); f[0] = "9"
L[k + 1] = "".join("%8s" % t for t in f); open(p, "w").write("\n".join(L))'
demo_on "$ROOT/test/numerics-features/conjugate_two_layer" bc_connection_out_of_range 1 \
    "a 101 line is \"bs is js ks fs d11 d12 d21 d22\" with a block and face of this mesh" "PY:$BC_CONN9"

echo
[ "$NMPI" -gt 1 ] && echo "$nmpi_ref demo(s) refused as a whole by an MPI guard (the ablating wall), code 1"
if [ "$fail" -eq 0 ]; then
  echo "all guards fired as intended"
else
  echo "$fail guard(s) did NOT fire -- they are decoration, not protection"
  exit 1
fi
