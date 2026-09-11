#!/bin/bash
#===============================================================================
#         FILE: test/run_init_sweep.sh
#        USAGE: ./test/run_init_sweep.sh [-p NTHREADS]
#  DESCRIPTION: Detector for the "read a variable that was never assigned"
#               defect family (plan 10, D8).
#
#               Builds the solver TWICE with opposite garbage-fill policies and
#               runs the whole suite under each. If any artefact differs, some
#               code path read uninitialised memory, and the differing case
#               says where to look.
#
#                 A   -finit-local-zero
#                       every local set to zero / .false.
#                 B   -finit-logical=true -finit-integer=-999999
#                     -finit-real=nan -finit-derived
#                       every local set to something conspicuous
#
#               THREE instances of this defect have been found in FUSS, each by
#               accident rather than by a test:
#
#                 blk%vol/M/dl ghost entries   results depended on the memory
#                                              layout of a derived type; 254 of
#                                              262 artefacts flipped
#                 coarse SOLUTIONTIME          1.33e+180 in file headers
#                 endsim / endmg               output file NAMES varied run to
#                                              run with identical solutions;
#                                              multigrid runs ended early
#
#               Three is a pattern. This is the check that would have caught
#               all three in one pass.
#
#               KNOWN LIMIT, measured not assumed: these flags apply to LOCAL
#               variables. gfortran does not reliably apply -finit-real to
#               allocatable components, which is how the first instance above
#               evaded an earlier -finit-real=snan attempt. A clean sweep is
#               therefore evidence about locals, not a proof about the whole
#               program. Allocated arrays are covered instead by explicit
#               initialisation in Mod_Allocate_Data.
#
#               PROVEN RED. Reverting the endsim/endmg initialisation in
#               Mod_Explicit and re-running just multimat_Plate under both
#               builds reproduces the original defect exactly:
#
#                 pass A (-finit-local-zero, endmg = .false.)  files 1, 2, 3
#                 pass B (-finit-logical=true, endmg = .true.) files 10, 20, 30
#
#               which is also the cleanest available explanation of why the
#               original bug looked like run-to-run randomness: the two values
#               the stack happened to hold are exactly the two the flags force.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NTHREADS=1

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -h | --help) echo "usage: $0 [-p NTHREADS]"; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

WORK="${TMPDIR:-/tmp}/fuss-init-sweep.$$"
mkdir -p "$WORK"

FLAGS_A="-finit-local-zero"
FLAGS_B="-finit-logical=true -finit-integer=-999999 -finit-real=nan -finit-derived"

# Base flags copied from the normal RELEASE build. gcc-ar / gcc-ranlib are not
# optional: with -flto and plain ar the FUSSL archive gets an empty symbol index
# and the link fails with undefined references that look like missing sources.
BASE="-fno-underscoring -march=native -cpp -ffree-line-length-none"

echo "FUSS uninitialised-read sweep"
echo "  work    : $WORK"
echo "  threads : $NTHREADS"
echo

# CLEANING UP AFTER THIS SCRIPT IS NOT OPTIONAL.
#
# CMakeLists.txt hands add_subdirectory an explicit binary dir for ORION and
# FiNeR: "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/". That path is fixed, so
# building from ANY build directory writes those subprojects into $ROOT/build
# -- while their .mod files go to the sweep's own module dir. The result is that
# $ROOT/build is left referring to module files that are not there, and the next
# ordinary `cmake --build build` fails with errors in ORION, PENF, FACE and
# FiNeR that look nothing like the cause. That is not hypothetical; it happened
# the first time this script was run.
#
# So the trap below rebuilds $ROOT/build from scratch. It costs a couple of
# minutes and it is the difference between a diagnostic and a trap for the next
# person.
SAVED="$WORK/FUSS.original"
[ -f "$ROOT/bin/FUSS" ] && cp "$ROOT/bin/FUSS" "$SAVED"

restore () {
  echo
  echo "  restoring $ROOT/build (the sweep clobbers ORION/FiNeR in it) ..."
  rm -rf "$ROOT/build"
  cmake -S "$ROOT" -B "$ROOT/build" \
        -DCMAKE_BUILD_TYPE=RELEASE \
        -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
        -DCMAKE_C_COMPILER=/usr/bin/gcc \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DORION_PATH="$ROOT/lib/ORION/" \
        -DFINER_PATH="$ROOT/lib/third_party/FiNeR/" \
        > "$WORK/restore.log" 2>&1 \
    && cmake --build "$ROOT/build" -j 8 >> "$WORK/restore.log" 2>&1 \
    && echo "  restored." \
    || echo "  RESTORE FAILED -- see $WORK/restore.log; rebuild build/ by hand."
}
trap restore EXIT

for pass in A B; do
  case $pass in
    A) EXTRA="$FLAGS_A" ;;
    B) EXTRA="$FLAGS_B" ;;
  esac

  bdir="$WORK/build-$pass"
  printf '  pass %s  configure+build ... ' "$pass"
  # The compiler MUST be named explicitly. A fresh configure on this machine
  # finds Intel's ifx first and fails in SetFortranFlags.cmake with "No compile
  # flags were found" -- and the -finit-* flags are gfortran-specific anyway.
  if ! cmake -S "$ROOT" -B "$bdir" \
        -DCMAKE_BUILD_TYPE=RELEASE \
        -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
        -DCMAKE_C_COMPILER=/usr/bin/gcc \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DCMAKE_AR=/usr/bin/gcc-ar-12 \
        -DCMAKE_RANLIB=/usr/bin/gcc-ranlib-12 \
        -DORION_PATH="$ROOT/lib/ORION/" \
        -DFINER_PATH="$ROOT/lib/third_party/FiNeR/" \
        -DCMAKE_Fortran_FLAGS="$BASE $EXTRA" \
        > "$WORK/$pass.cmake.log" 2>&1; then
    echo "FAILED to configure -- see $WORK/$pass.cmake.log"; exit 1
  fi
  if ! cmake --build "$bdir" -j 8 > "$WORK/$pass.build.log" 2>&1; then
    echo "FAILED to build -- see $WORK/$pass.build.log"; exit 1
  fi
  echo "done"

  # No binary shuffling needed: src/app/CMakeLists.txt sets
  # CMAKE_RUNTIME_OUTPUT_DIRECTORY to ${CMAKE_SOURCE_DIR}/bin, so an
  # out-of-source build still writes $ROOT/bin/FUSS -- which is exactly what
  # every case's FUSS.sh copies. Build then run, in that order, per pass.
  # `touch` so the -nt check in FUSS.sh always refreshes the case-local copy.
  touch "$ROOT/bin/FUSS"

  printf '  pass %s  running suite ... ' "$pass"
  if "$ROOT/test/run_regression.sh" -p "$NTHREADS" -o "$WORK/$pass" \
       > "$WORK/$pass.run.log" 2>&1; then
    echo "$(wc -l < "$WORK/$pass/manifest.txt") artefacts"
  else
    echo "FAILED -- see $WORK/$pass.run.log"; exit 1
  fi
done

echo
if diff -q "$WORK/A/manifest.txt" "$WORK/B/manifest.txt" > /dev/null; then
  echo "  PASS  opposite initialisation policies gave identical artefacts"
  echo "        (no local variable is read before being assigned)"
  exit 0
else
  echo "  FAIL  the two builds disagree -- something reads uninitialised memory"
  echo
  diff "$WORK/A/manifest.txt" "$WORK/B/manifest.txt" \
    | grep '^[<>]' | cut -f1 | sed 's/^[<>] //' | sort -u | sed 's/^/          /'
  echo
  echo "  Those cases are where to look. Note the symptom may be subtler than a"
  echo "  wrong number: the endsim/endmg instance changed only output file NAMES"
  echo "  while every solution stayed bit-identical."
  exit 1
fi
