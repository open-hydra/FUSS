#!/bin/bash
#===============================================================================
#         FILE: test/run_morph_guards.sh
#        USAGE: ./test/run_morph_guards.sh
#  DESCRIPTION: Architecture guards for the MORPH mesh component (plan 05 5.0,
#               plan 09 section 2), which FUSS builds from the submodule
#               lib/MORPH (or the tree MORPH_PATH names: the environment, else
#               build/CMakeCache.txt).
#
#               The entire value of MORPH being a separate component rests on
#               one rule -- it must never depend on FUSS -- and on it staying
#               true as later phases add code. A design rule nobody checks
#               decays within two phases, so it is checked here, mechanically,
#               in milliseconds.
#
#               FUSS's side: exactly one FUSS file names MORPH types (the
#               adapter). MORPH's side: MORPH's own guards
#               (lib/MORPH/test/run_morph_guards.sh) on the tree FUSS builds:
#               no FUSS module, no stop, no I/O, no module but its own, a
#               standalone -Wall-clean compile, and the unit programs at 1 and
#               8 threads with the same output.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MORPH=${MORPH_PATH:-$(sed -n 's/^MORPH_PATH:[A-Z]*=//p' "$ROOT/build/CMakeCache.txt" 2>/dev/null)}
MORPH=${MORPH:-$ROOT/lib/MORPH}
MORPH=${MORPH%/}

fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "MORPH architecture guards"
echo "  component: $MORPH"
echo

if [ ! -x "$MORPH/test/run_morph_guards.sh" ]; then
  bad "no MORPH at $MORPH: run git submodule update --init lib/MORPH (or ./install.sh update)"
  echo; echo "1 MORPH guard(s) FAILED"; exit 1
fi

# ---------------------------------------------------------------------------
# (b) Exactly one FUSS file may name MORPH types (the adapter, plan 09
#     section 7). Two or more and the seam is no longer a single file.
# ---------------------------------------------------------------------------
adapters=$(grep -rIl --include='*.f90' -E '^[[:space:]]*use[[:space:]]+Morph_' "$ROOT/src" 2>/dev/null || true)
n_adapters=$(printf '%s' "$adapters" | grep -c . || true)
if [ "$n_adapters" -eq 1 ]; then
  pass "exactly one FUSS file names MORPH: $(basename "$adapters")"
else
  bad "$n_adapters FUSS files name MORPH types; the seam must be a single adapter:"
  echo "$adapters" | sed 's/^/          /'
fi

# ---------------------------------------------------------------------------
# MORPH's guards, on the tree FUSS builds (MORPH's own wall-clean-app needs a
# build of MORPH's driver, which FUSS does not configure: MORPH's ctest has it)
# ---------------------------------------------------------------------------
echo
if ! bash "$MORPH/test/run_morph_guards.sh" no-fuss no-stop no-io lib-isolated wall-clean units; then
  fail=$((fail+1))
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "all MORPH guards passed"
else
  echo "$fail MORPH guard set(s) FAILED"
  exit 1
fi
