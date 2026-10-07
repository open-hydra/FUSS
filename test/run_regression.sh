#!/bin/bash
#===============================================================================
#         FILE: test/run_regression.sh
#        USAGE: ./test/run_regression.sh [-p N] [-o OUTDIR] [-c CASE]...
#  DESCRIPTION: Run the FUSS test cases and record a fingerprint manifest of
#               their solution files, so two builds can be compared exactly.
#
#               Used for two distinct jobs:
#                 1. measuring the run-to-run / thread-to-thread determinism
#                    floor (run it twice with identical settings and diff)
#                 2. the phase regression gates -- e.g. plan 05 section 5.2
#                    "zero mesh velocity must be bit-identical to Phase 0"
#
#               SOLUTIONTIME in the .tec headers is NOT stripped: for a steady
#               run it carries the converged iteration count, so a change in it
#               is a real result change, not formatting noise.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NTHREADS=1
export LC_ALL=C   # `find | sort` below orders the manifest rows; a shell with another collation moves rows (multimat_Plate field1/field10/field2) and -k reports a false difference
OUTDIR=""
ALL=0
FRESH=0
CHECK=0
BASELINE="${FUSS_BASELINE:-$ROOT/test/baseline/manifest.txt}"   # env override exists so the -k path can be TESTED
CASES=()

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -o | --out)      OUTDIR=$2;   shift 2 ;;
    -c | --case)     CASES+=("$2"); shift 2 ;;
    -a | --all)      ALL=1; shift ;;
    -B | --fresh-build)    FRESH=1; shift ;;
    -k | --check-baseline) CHECK=1; shift ;;
    -h | --help)
      echo "usage: $0 [-p NTHREADS] [-o OUTDIR] [-a] [-B] [-k] [-c CASE]..."
      echo "  with no -c, every test/**/FUSS.sh case is run"
      echo "  -a also runs cases marked .slow (grid-convergence studies);"
      echo "     without it they are skipped, so keep -a consistent between"
      echo "     the two manifests you intend to diff"
      echo "  -B wipes build/ and rebuilds from the preset before running."
      echo "     REQUIRED when recording a baseline -- see the note in the script."
      echo "  -k after the run, diff the manifest against test/baseline/manifest.txt"
      echo "     and exit 1 on any difference (only meaningful for a full run)"
      exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

# ---------------------------------------------------------------------------
# WHY -B EXISTS. On 2026-09-14 an INCREMENTAL relink of the tree -- one changed
# module, `cmake --build build` on an existing build dir -- produced a binary
# whose results differed from a fresh preset build of the SAME SOURCE on 90 of
# 290 artefacts, including static single-grid cases. The fresh build matched a
# fresh build of the previous commit exactly, so the source change was inert;
# the incremental binary was the odd one out. The build dir was wiped before
# it could be examined, so the mechanism is unknown (LTO whole-program codegen
# is one candidate). What is known: two fresh preset builds agree with each
# other, and an incremental one may not. So a baseline is recorded only from a
# fresh build, and this script can do that itself so nobody has to remember.
# ---------------------------------------------------------------------------
if [ "$FRESH" -eq 1 ]; then
  echo "fresh build requested: wiping build/ and configuring from the preset"
  ( cd "$ROOT" && rm -rf build && cmake --preset default > /dev/null 2>&1 \
      && cmake --build build -j 8 > /dev/null 2>&1 ) \
    || { echo "fresh build FAILED"; exit 1; }
fi

if [ ${#CASES[@]} -eq 0 ]; then
  # deterministic ordering matters: the manifest is diffed line by line
  while IFS= read -r c; do
    # A .slow marker means the case exists to be MEASURED (a refinement study
    # with its own verify.py), not to be fingerprinted for bit-identity, and it
    # costs minutes rather than seconds. Excluding it by default keeps this
    # script something people actually run.
    if [ "$ALL" -eq 0 ] && [ -f "$ROOT/$c/.slow" ]; then continue; fi
    CASES+=("$c")
  done < <(
    cd "$ROOT" && find test -name FUSS.sh -printf '%h\n' | sort
  )
fi

# Refuse to fingerprint anything from a build that was not configured the way
# CMakePresets.json says. A mis-configured build is silent -- it compiles, runs,
# and agrees with the baseline on all but one artefact -- so this has to be
# checked rather than remembered. See test/check_build_config.sh.
if ! "$ROOT/test/check_build_config.sh"; then
  echo
  echo "refusing to record a manifest from a build that does not match the preset"
  exit 1
fi

[ -n "$OUTDIR" ] || OUTDIR="$ROOT/test/.regression/omp${NTHREADS}"
mkdir -p "$OUTDIR"
MANIFEST="$OUTDIR/manifest.txt"
: > "$MANIFEST"

echo "FUSS regression run"
echo "  root     : $ROOT"
echo "  threads  : $NTHREADS"
echo "  cases    : ${#CASES[@]}"
echo "  manifest : $MANIFEST"
echo

fail=0
for case_dir in "${CASES[@]}"; do
  abs="$ROOT/$case_dir"
  if [ ! -x "$abs/FUSS.sh" ]; then
    echo "SKIP  $case_dir (no FUSS.sh)"
    continue
  fi

  # Start from a clean slate so a stale result cannot masquerade as a fresh one.
  #
  # SOLUTION/ is deliberately NOT touched. It holds committed REFERENCE data,
  # and the solver never writes there -- every output path in Write_vtk_tec is
  # rooted at 'OUTPUT/'. An earlier version of this script wiped it too, which
  # quietly deleted ten tracked files and 341k lines of reference solutions the
  # first time the suite was run.
  rm -rf "$abs/OUTPUT" "$abs/bin"

  printf '%-52s' "RUN   $case_dir"
  log="$OUTDIR/$(echo "$case_dir" | tr '/' '_').log"
  if ( cd "$abs" && ./FUSS.sh -p "$NTHREADS" solve ) > "$log" 2>&1; then
    printf 'ok\n'
  else
    printf 'FAILED (see %s)\n' "$log"
    fail=$((fail+1))
    echo "$case_dir	<RUN-FAILED>" >> "$MANIFEST"
    continue
  fi

  # Fingerprint every solution artefact, in sorted order.
  #
  # OUTPUT/ only. SOLUTION/ is static reference data checked into the
  # repository, so fingerprinting it would add the same constant lines to every
  # manifest -- and, because the earlier version of this script emptied it
  # before this find ran, no existing baseline contains those lines. Keeping to
  # OUTPUT keeps new manifests comparable with the old ones.
  while IFS= read -r f; do
    rel=${f#"$abs"/}
    sum=$(md5sum "$f" | cut -d' ' -f1)
    echo "$case_dir	$rel	$sum" >> "$MANIFEST"
  done < <(find "$abs/OUTPUT" -type f 2>/dev/null | sort)
done

echo
echo "wrote $(wc -l < "$MANIFEST") manifest lines to $MANIFEST"
[ "$fail" -eq 0 ] || { echo "$fail case(s) failed to run"; exit 1; }

# Compare against the committed baseline. This is the regression GATE; the
# manifest above is only its input. A baseline that lives in a temp directory
# is not a baseline -- every one recorded before 2026-09-14 was lost that way.
if [ "$CHECK" -eq 1 ]; then
  if [ ! -f "$BASELINE" ]; then
    echo "no committed baseline at $BASELINE"; exit 1
  fi
  if diff "$BASELINE" "$MANIFEST" > "$OUTDIR/baseline.diff"; then
    echo "baseline: IDENTICAL to $BASELINE"
  else
    echo "baseline: DIFFERS from $BASELINE in these cases:"
    grep '^[<>]' "$OUTDIR/baseline.diff" | cut -f1 | sed 's/^[<>] //' | sort -u | sed 's/^/    /'
    echo "  full diff: $OUTDIR/baseline.diff"
    echo "  If the change is intended, re-record with: $0 -B -o <dir> && cp <dir>/manifest.txt $BASELINE"
    exit 1
  fi
fi
