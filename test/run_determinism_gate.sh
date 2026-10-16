#!/bin/bash
#===============================================================================
#         FILE: test/run_determinism_gate.sh
#        USAGE: ./test/run_determinism_gate.sh [-p NTHREADS] [-m NRANKS] [-a] [-x] [-X] [-c CASE]...
#  DESCRIPTION: Run the whole suite TWICE with identical settings and require
#               the two manifests to be identical.
#
#               This is not a formality. It is the gate that found the
#               uninitialised `endmg` in Mod_Explicit: on every single-grid run
#               that flag was read without ever being assigned, so whatever sat
#               on the stack decided whether a step counted as an output step.
#               The solution fields were bit-identical either way -- only the
#               NUMBERING of the snapshot files moved, 1..101 one run and
#               10..1010 the next -- which is precisely why it survived so long.
#               It looked like a quirk of the test harness.
#
#               The same flag pair, read on the other branch, can end a
#               multigrid run early and report it as normal completion.
#
#               Every comparison gate in this repository rests on the same
#               executable producing the same bytes twice. That assumption is
#               cheap to check and expensive to be wrong about, so check it.
#
#               A difference here does NOT mean the physics changed. Read the
#               diff before concluding anything: a manifest entry can move
#               because a file was renamed, not because its contents changed.
#
#               -x: a SEPARATE section, run after the normal
#               two-pass comparison above and only when asked for, that checks
#               determinism across THREAD COUNTS rather than across repeated
#               runs at a fixed count. It solves moving-mesh/freestream at
#               -p 1 and -p 8 in a scratch copy (never the in-tree OUTPUT/,
#               same reasoning as run_physics_gates.sh's `stage`) and requires
#               OUTPUT/field.tec to be byte-identical between the two. This is
#               the moving-mesh claim ("identical at 1 and 8 threads"),
#               which until 2026-09-15 had only ever been checked by hand
#               (field.tec md5 identical, max |T-T0| 9.663e-13 K at both
#               thread counts). Only field.tec is asserted, not the whole
#               regression manifest: wall.tec and residual-history.dat carry
#               OpenMP-reduction residual norms, whose summation order can
#               legitimately move across thread counts in the last bit even
#               when the FIELD -- the thing the claim is actually about -- is
#               unaffected, so routing this through run_regression.sh's
#               manifest diff would risk a false FAIL that this gate did not
#               intend to make. Absent -x, behaviour is exactly as before.
#
#               Since plan 07 the -x section runs a LIST of cases, each at
#               -p 1 and -p 8 with field.tec compared byte for byte:
#               moving-mesh/freestream (above), and the Darcy gates 4.1(b)
#               pyrolysis/darcy/mass_balance_1d (the line-implicit tridiag
#               solve, lines distributed over threads) and 4.3
#               pyrolysis/darcy/tensor_linear_2d (the damped Jacobi sweep,
#               cells distributed over threads, and a run that stops on its
#               convergence residual). -c CASE (path under test/, repeatable)
#               replaces the list and implies -x; -X runs the cross-thread
#               section alone, without the two full passes (development use:
#               the two passes are the gate).
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NTHREADS=1
EXTRA=""
XTHREAD=0
XONLY=0
XCASES=()
NMPI=0

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -m | --mpi)      NMPI=$2; shift 2 ;;
    -a | --all)      EXTRA="-a"; shift ;;
    -x | --cross-thread) XTHREAD=1; shift ;;
    -X | --cross-thread-only) XTHREAD=1; XONLY=1; shift ;;
    -c | --case) XTHREAD=1; XCASES+=("$2"); shift 2 ;;
    -h | --help) echo "usage: $0 [-p NTHREADS] [-m NRANKS] [-a] [-x] [-X] [-c CASE]..."; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
[ "${#XCASES[@]}" -gt 0 ] || XCASES=(moving-mesh/freestream pyrolysis/darcy/mass_balance_1d pyrolysis/darcy/tensor_linear_2d)

WORK="${TMPDIR:-/tmp}/fuss-determinism.$$"
mkdir -p "$WORK"

# -m N (plan 11): both passes and the cross-thread runs use the mpi preset's bin/FUSS-mpi on N ranks
# through Intel MPI's own launcher (run_regression.sh -m does the same for the passes).
EXE=FUSS ; MPIOPT="" ; TMO=""
if [ "$NMPI" -gt 0 ]; then
  [ -n "${I_MPI_ROOT:-}" ] || { echo "-m needs I_MPI_ROOT (the Intel MPI of the mpi preset)"; exit 1; }
  EXE=FUSS-mpi ; MPIOPT="-m $NMPI"
  [ "$NMPI" -gt 1 ] && TMO="timeout 900"
  export FUSS_EXE=FUSS-mpi FUSS_MPIRUN="$I_MPI_ROOT/bin/mpirun"
fi

echo "FUSS determinism gate"
echo "  threads : $NTHREADS"
[ "$NMPI" -gt 0 ] && echo "  ranks   : $NMPI (bin/FUSS-mpi)"
echo "  work    : $WORK"
echo

for pass in A B; do
  [ "$XONLY" -eq 1 ] && break
  printf '  pass %s ... ' "$pass"
  if "$ROOT/test/run_regression.sh" -p "$NTHREADS" $MPIOPT $EXTRA -o "$WORK/$pass" \
       > "$WORK/$pass.log" 2>&1; then
    echo "done ($(wc -l < "$WORK/$pass/manifest.txt") artefacts)"
  elif [ "$NMPI" -gt 1 ]; then
    # On N ranks the cases MPI refuses as a whole (pyrolysis, the ablating wall, chimera, multigrid)
    # exit 1 by name; any other case that did not run is a failure
    crashed=0
    for lg in $(grep -l . "$WORK/$pass"/*.log); do
      grep -q "FAILED (see $lg)" "$WORK/$pass.log" || continue
      if ! grep -qE "pyrolysis is OpenMP-only|the ablating wall \(306\) is OpenMP-only|chimera faces \(102\) are OpenMP-only|multigrid is OpenMP-only" "$lg"; then
        echo; echo "  FAIL  $(basename "$lg" .log) did not run and was not refused by an MPI guard (see $lg)"
        crashed=$((crashed+1))
      fi
    done
    [ "$crashed" -eq 0 ] || exit 1
    echo "done ($(wc -l < "$WORK/$pass/manifest.txt") artefacts; $(grep -c 'FAILED (see' "$WORK/$pass.log") cases refused by an MPI guard)"
  else
    echo "FAILED -- see $WORK/$pass.log"
    exit 1
  fi
done

echo
if [ "$XONLY" -eq 1 ]; then
  echo "  (-X: the two full passes skipped)"
elif diff -q "$WORK/A/manifest.txt" "$WORK/B/manifest.txt" > /dev/null; then
  echo "  PASS  two identical runs produced identical manifests"
else
  echo "  FAIL  the same binary produced different artefacts on two runs"
  echo
  diff "$WORK/A/manifest.txt" "$WORK/B/manifest.txt" | head -40
  echo
  echo "  Until this passes, every bit-identity gate in test/ is meaningless."
  exit 1
fi

# ---------------------------------------------------------------------------
# -x: cross-thread check, see the header note above.
# ---------------------------------------------------------------------------
if [ "$XTHREAD" -eq 1 ]; then
  echo
  XWORK="$WORK/xthread"
  mkdir -p "$XWORK"

  # Snapshot bin/FUSS ONCE, up front, so the -p 1 and -p 8 runs are guaranteed
  # to grade the same executable regardless of what else touches
  # $ROOT/bin/FUSS meanwhile: bin/ can be rebuilt by a concurrent build
  # and was observed CHURNING mid-run (same reasoning as
  # run_physics_gates.sh's BINSNAP; see that file's header for the failure it
  # was written to avoid).
  cp "$ROOT/bin/$EXE" "$XWORK/FUSS.snapshot"

  xfail=0
  for xc in "${XCASES[@]}"; do
    echo "  cross-thread check: $xc, -p 1 vs -p 8"
    SRC="$ROOT/test/$xc"
    xname=$(echo "$xc" | tr '/' '_')
    if [ ! -f "$SRC/FUSS.sh" ]; then
      echo "  FAIL  no case at test/$xc"
      xfail=$((xfail+1)); continue
    fi
    for nt in 1 8; do
      d="$XWORK/$xname/p$nt"
      rm -rf "$d"; mkdir -p "$d"
      cp -r "$SRC/INPUT" "$SRC/MESH" "$SRC/input.ini" "$SRC/FUSS.sh" "$d/"
      cp "$SRC"/*.dat "$d/" 2>/dev/null || true
      sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
      mkdir -p "$d/OUTPUT" "$d/bin"
      # Work copy's mtime is pushed into the future so FUSS.sh's own runtime
      # `-nt` copy check (which still points at the LIVE $ROOT/bin/FUSS via the
      # sedded MASTERDIR) is a guaranteed no-op and cannot overwrite the
      # snapshot with a possibly-different, possibly-mid-write binary.
      cp "$XWORK/FUSS.snapshot" "$d/bin/$EXE"
      touch -d '+1 hour' "$d/bin/$EXE"

      printf '    -p %d ... ' "$nt"
      if ( cd "$d" && $TMO ./FUSS.sh -p "$nt" $MPIOPT solve ) > "$d/run.log" 2>&1; then
        echo "done"
      else
        echo "FAILED -- see $d/run.log"
        exit 1
      fi
    done

    m1=$(md5sum "$XWORK/$xname/p1/OUTPUT/field.tec" | cut -d' ' -f1)
    m8=$(md5sum "$XWORK/$xname/p8/OUTPUT/field.tec" | cut -d' ' -f1)
    echo "    field.tec md5  -p 1: $m1"
    echo "    field.tec md5  -p 8: $m8"
    if [ "$m1" = "$m8" ]; then
      echo "  PASS  $xc field.tec is bit-identical at 1 and 8 threads"
    else
      echo "  FAIL  $xc field.tec differs between 1 and 8 threads -- see $XWORK/$xname"
      xfail=$((xfail+1))
    fi
  done
  [ "$xfail" -eq 0 ] || exit 1
fi
