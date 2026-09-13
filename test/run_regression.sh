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
OUTDIR=""
ALL=0
CASES=()

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -o | --out)      OUTDIR=$2;   shift 2 ;;
    -c | --case)     CASES+=("$2"); shift 2 ;;
    -a | --all)      ALL=1; shift ;;
    -h | --help)
      echo "usage: $0 [-p NTHREADS] [-o OUTDIR] [-a] [-c CASE]..."
      echo "  with no -c, every test/**/FUSS.sh case is run"
      echo "  -a also runs cases marked .slow (grid-convergence studies);"
      echo "     without it they are skipped, so keep -a consistent between"
      echo "     the two manifests you intend to diff"
      exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

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
