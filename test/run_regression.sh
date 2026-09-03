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
CASES=()

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -o | --out)      OUTDIR=$2;   shift 2 ;;
    -c | --case)     CASES+=("$2"); shift 2 ;;
    -h | --help)
      echo "usage: $0 [-p NTHREADS] [-o OUTDIR] [-c CASE]..."
      echo "  with no -c, every test/**/FUSS.sh case is run"
      exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

if [ ${#CASES[@]} -eq 0 ]; then
  # deterministic ordering matters: the manifest is diffed line by line
  while IFS= read -r c; do CASES+=("$c"); done < <(
    cd "$ROOT" && find test -name FUSS.sh -printf '%h\n' | sort
  )
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

  # Start from a clean slate so a stale SOLUTION cannot masquerade as a result.
  rm -rf "$abs/SOLUTION" "$abs/OUTPUT" "$abs/bin"
  mkdir -p "$abs/SOLUTION"

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
  while IFS= read -r f; do
    rel=${f#"$abs"/}
    sum=$(md5sum "$f" | cut -d' ' -f1)
    echo "$case_dir	$rel	$sum" >> "$MANIFEST"
  done < <(find "$abs/SOLUTION" "$abs/OUTPUT" -type f 2>/dev/null | sort)
done

echo
echo "wrote $(wc -l < "$MANIFEST") manifest lines to $MANIFEST"
[ "$fail" -eq 0 ] || { echo "$fail case(s) failed to run"; exit 1; }
