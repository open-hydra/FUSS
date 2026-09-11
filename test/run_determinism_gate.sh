#!/bin/bash
#===============================================================================
#         FILE: test/run_determinism_gate.sh
#        USAGE: ./test/run_determinism_gate.sh [-p NTHREADS] [-a]
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
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NTHREADS=1
EXTRA=""

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -a | --all)      EXTRA="-a"; shift ;;
    -h | --help) echo "usage: $0 [-p NTHREADS] [-a]"; exit 1 ;;
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
