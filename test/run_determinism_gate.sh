#!/bin/bash
#===============================================================================
#         FILE: test/run_determinism_gate.sh
#        USAGE: ./test/run_determinism_gate.sh [-p NTHREADS] [-a] [-x]
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
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NTHREADS=1
EXTRA=""
XTHREAD=0

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -a | --all)      EXTRA="-a"; shift ;;
    -x | --cross-thread) XTHREAD=1; shift ;;
    -h | --help) echo "usage: $0 [-p NTHREADS] [-a] [-x]"; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

WORK="${TMPDIR:-/tmp}/fuss-determinism.$$"
mkdir -p "$WORK"

echo "FUSS determinism gate"
echo "  threads : $NTHREADS"
echo "  work    : $WORK"
echo

for pass in A B; do
  printf '  pass %s ... ' "$pass"
  if "$ROOT/test/run_regression.sh" -p "$NTHREADS" $EXTRA -o "$WORK/$pass" \
       > "$WORK/$pass.log" 2>&1; then
    echo "done ($(wc -l < "$WORK/$pass/manifest.txt") artefacts)"
  else
    echo "FAILED -- see $WORK/$pass.log"
    exit 1
  fi
done

echo
if diff -q "$WORK/A/manifest.txt" "$WORK/B/manifest.txt" > /dev/null; then
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
  echo "  cross-thread check: moving-mesh/freestream, -p 1 vs -p 8"
  SRC="$ROOT/test/moving-mesh/freestream"
  XWORK="$WORK/xthread"
  mkdir -p "$XWORK"

  # Snapshot bin/FUSS ONCE, up front, so the -p 1 and -p 8 runs are guaranteed
  # to grade the same executable regardless of what else touches
  # $ROOT/bin/FUSS meanwhile: bin/ can be rebuilt by a concurrent build
  # and was observed CHURNING mid-run (same reasoning as
  # run_physics_gates.sh's BINSNAP; see that file's header for the failure it
  # was written to avoid).
  cp "$ROOT/bin/FUSS" "$XWORK/FUSS.snapshot"

  for nt in 1 8; do
    d="$XWORK/p$nt"
    rm -rf "$d"; mkdir -p "$d"
    cp -r "$SRC/INPUT" "$SRC/MESH" "$SRC/input.ini" "$SRC/FUSS.sh" "$d/"
    sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
    mkdir -p "$d/OUTPUT" "$d/bin"
    # Work copy's mtime is pushed into the future so FUSS.sh's own runtime
    # `-nt` copy check (which still points at the LIVE $ROOT/bin/FUSS via the
    # sedded MASTERDIR) is a guaranteed no-op and cannot overwrite the
    # snapshot with a possibly-different, possibly-mid-write binary.
    cp "$XWORK/FUSS.snapshot" "$d/bin/FUSS"
    touch -d '+1 hour' "$d/bin/FUSS"

    printf '    -p %d ... ' "$nt"
    if ( cd "$d" && ./FUSS.sh -p "$nt" solve ) > "$d/run.log" 2>&1; then
      echo "done"
    else
      echo "FAILED -- see $d/run.log"
      exit 1
    fi
  done

  m1=$(md5sum "$XWORK/p1/OUTPUT/field.tec" | cut -d' ' -f1)
  m8=$(md5sum "$XWORK/p8/OUTPUT/field.tec" | cut -d' ' -f1)
  echo "    field.tec md5  -p 1: $m1"
  echo "    field.tec md5  -p 8: $m8"
  if [ "$m1" = "$m8" ]; then
    echo "  PASS  freestream field.tec is bit-identical at 1 and 8 threads"
  else
    echo "  FAIL  freestream field.tec differs between 1 and 8 threads -- see $XWORK"
    exit 1
  fi
fi
