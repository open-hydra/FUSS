#!/bin/bash
#===============================================================================
#         FILE: test/run_mpi_gate.sh
#        USAGE: ./test/run_mpi_gate.sh [-w WORKDIR] [-c CASE:NRANKS]...
#  DESCRIPTION: A case run on N MPI ranks must give the outputs of the same case
#               run on one rank (plan 11 section 2.5). Both runs use the mpi
#               preset's bin/FUSS-mpi at OMP_NUM_THREADS=1, the N-rank one through
#               $I_MPI_ROOT/bin/mpirun under timeout 900 (a deadlock is a FAIL,
#               not a hang).
#
#               Rules: every OUTPUT file byte-identical, except
#                 - wall*.tec, compared with -0.0 equal to +0.0: the wall values
#                   are summed over the ranks, each entry from the rank that owns
#                   it, and a -0.0 of the owner becomes +0.0 once another rank's
#                   zero is added (the decision of 2026-10-06 keeps the sum);
#                 - residual-history.dat, not compared: the residual is an MPI sum
#                   whose rounding depends on the rank count.
#               So a case must stop on a fixed iteration count or time, never on
#               res-threshold (conjugate_two_layer_fixed exists for that).
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK="${TMPDIR:-/tmp}/fuss-mpi-gate.$$"
CASES=()

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -c | --case) CASES+=("$2"); shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR] [-c CASE:NRANKS]..."; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
# the static legacy cases of plan 11 stage 1: three blocks on four ranks (one idle), and two blocks with
# IRS (the R exchange) on two
[ "${#CASES[@]}" -gt 0 ] || CASES=("test/transient/multimat_Plate:4" "test/numerics-features/conjugate_two_layer_fixed:2")

[ -n "${I_MPI_ROOT:-}" ] || { echo "needs I_MPI_ROOT (the Intel MPI of the mpi preset)"; exit 1; }
export FUSS_EXE=FUSS-mpi FUSS_MPIRUN="$I_MPI_ROOT/bin/mpirun"
"$ROOT/test/check_build_config.sh" --preset mpi || { echo "refusing to grade a build that does not match the preset"; exit 1; }

mkdir -p "$WORK"
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS MPI gate: N ranks against one"
echo "  binary   : $ROOT/bin/FUSS-mpi ($(md5sum "$ROOT/bin/FUSS-mpi" | cut -c1-8))"
echo "  work dir : $WORK"
echo

for spec in "${CASES[@]}"; do
  case_rel=${spec%%:*}; n=${spec##*:}
  src="$ROOT/$case_rel"; name=$(echo "$case_rel" | tr '/' '_')
  for r in 1 "$n"; do
    d="$WORK/$name.m$r"
    rm -rf "$d"; mkdir -p "$d"
    cp -r "$src/INPUT" "$src/input.ini" "$src/FUSS.sh" "$d/"
    [ -d "$src/MESH" ] && cp -r "$src/MESH" "$d/"
    cp "$src"/*.dat "$d/" 2>/dev/null || true
    sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
    mkdir -p "$d/OUTPUT" "$d/bin"
    ( cd "$d" && timeout 900 ./FUSS.sh -p 1 -m "$r" solve ) > "$d/run.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ] || grep -q "\[ERROR\]\|NaN or T<0" "$d/run.log"; then
      bad "$case_rel: the $r-rank run failed (exit $rc$([ $rc -eq 124 ] && echo ', timed out'); see $d/run.log)"
      continue 2
    fi
  done
  out=$(python3 - "$WORK/$name.m1/OUTPUT" "$WORK/$name.m$n/OUTPUT" <<'PY'
import os, sys
a, b = sys.argv[1], sys.argv[2]
fa = sorted(f for f in os.listdir(a) if f != 'residual-history.dat')
fb = sorted(f for f in os.listdir(b) if f != 'residual-history.dat')
if fa != fb:
    print('different file sets: %s' % sorted(set(fa) ^ set(fb))); sys.exit(1)
same = signed_zero = 0; bad = []
for f in fa:
    x = open(os.path.join(a, f), 'rb').read(); y = open(os.path.join(b, f), 'rb').read()
    if x == y:
        same += 1; continue
    if f.startswith('wall'):
        nz = lambda s: s.replace(b'-0.000000000000000E+000', b' 0.000000000000000E+000')   # the sign's column
        if nz(x) == nz(y):
            signed_zero += 1; continue
    bad.append(f)
print('%d files: %d byte-identical, %d wall files equal up to the sign of zero%s'
      % (len(fa), same, signed_zero, '' if not bad else '; DIFFERENT: ' + ', '.join(bad[:6])))
sys.exit(1 if bad else 0)
PY
)
  if [ $? -eq 0 ]; then
    pass "$case_rel: $n ranks = 1 rank ($out)"
  else
    bad "$case_rel: $n ranks differ from 1 rank ($out)"
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "MPI gate passed"
else
  echo "$fail case(s) FAILED"
  exit 1
fi
