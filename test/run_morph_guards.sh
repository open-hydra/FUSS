#!/bin/bash
#===============================================================================
#         FILE: test/run_morph_guards.sh
#        USAGE: ./test/run_morph_guards.sh
#  DESCRIPTION: Architecture guards for the MORPH mesh component (plan 05 5.0,
#               plan 09 section 2).
#
#               The entire value of MORPH being a separate component rests on
#               one rule -- it must never depend on FUSS -- and on it staying
#               true as later phases add code. A design rule nobody checks
#               decays within two phases, so it is checked here, mechanically,
#               in milliseconds.
#
#               Also builds MORPH standalone and runs its unit tests.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MORPH="$ROOT/src/lib/morph"
FC=${FC:-gfortran}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "MORPH architecture guards"
echo "  component: $MORPH"
echo

# ---------------------------------------------------------------------------
# (a) MORPH must not use any FUSS module. This is THE rule.
# ---------------------------------------------------------------------------
hits=$(grep -rIl --include='*.f90' -E '^[[:space:]]*use[[:space:]]+FUSS_' "$MORPH" 2>/dev/null || true)
if [ -z "$hits" ]; then
  pass "MORPH uses no FUSS module"
else
  bad "MORPH depends on FUSS -- extraction is no longer a directory move:"
  echo "$hits" | sed 's/^/          /'
fi

# ---------------------------------------------------------------------------
# (b) Exactly one FUSS file may name MORPH types (the adapter).
#     Absent until the adapter lands; zero is acceptable, two or more is not.
# ---------------------------------------------------------------------------
adapters=$(grep -rIl --include='*.f90' -E '^[[:space:]]*use[[:space:]]+Morph_' "$ROOT/src/lib" 2>/dev/null \
           | grep -v "^$MORPH/" || true)
n_adapters=$(printf '%s' "$adapters" | grep -c . || true)
if [ "$n_adapters" -le 1 ]; then
  if [ "$n_adapters" -eq 0 ]; then
    pass "no FUSS file names MORPH yet (adapter not landed)"
  else
    pass "exactly one FUSS file names MORPH: $(basename "$adapters")"
  fi
else
  bad "$n_adapters FUSS files name MORPH types; the seam must be a single adapter:"
  echo "$adapters" | sed 's/^/          /'
fi

# ---------------------------------------------------------------------------
# (c) MORPH must never abort. A library that stops cannot be embedded.
#     `error stop` in the test programs is fine -- those are executables.
# ---------------------------------------------------------------------------
stops=$(grep -rn --include='*.f90' -E '^[[:space:]]*stop\b' "$MORPH" 2>/dev/null \
        | grep -v '/test/' || true)
if [ -z "$stops" ]; then
  pass "MORPH contains no 'stop' outside its test programs"
else
  bad "MORPH aborts the program instead of returning a status:"
  echo "$stops" | sed 's/^/          /'
fi

# ---------------------------------------------------------------------------
# (d) MORPH compiles standalone, with no FUSS sources on the path.
#     If this fails, something crossed the seam that (a) did not catch.
# ---------------------------------------------------------------------------
# Dependency order matters: Morph_Metrics uses Morph_GCL for the cell volume,
# so GCL must be compiled first or its .mod does not exist yet.
SRC="
$MORPH/base/Morph_Types_m.f90
$MORPH/metrics/Morph_GCL.f90
$MORPH/metrics/Morph_Metrics.f90
$MORPH/motion/Morph_Motion_m.f90
$MORPH/motion/Morph_Motion_Static.f90
$MORPH/motion/Morph_Motion_Prescribed.f90
$MORPH/quality/Morph_Quality.f90
$MORPH/Morph_API.f90
"
if ( cd "$WORK" && $FC -c -O2 -fopenmp -fimplicit-none -Wall \
       -Wno-unused-dummy-argument $SRC ) > "$WORK/build.log" 2>&1; then
  pass "MORPH compiles standalone (-Wall clean)"
else
  bad "MORPH does not compile standalone:"
  sed 's/^/          /' "$WORK/build.log" | head -20
fi

# ---------------------------------------------------------------------------
# (e) Unit tests: the GCL identity, single cell and whole block.
# ---------------------------------------------------------------------------
for prog in test_morph_gcl test_morph_block; do
  if ( cd "$WORK" && $FC -O2 -fopenmp -fimplicit-none -o "$prog" \
         $SRC "$MORPH/test/$prog.f90" ) >> "$WORK/build.log" 2>&1; then
    for nt in 1 8; do
      if ( cd "$WORK" && OMP_NUM_THREADS=$nt "./$prog" ) > "$WORK/$prog.$nt.out" 2>&1; then
        pass "$prog (OMP_NUM_THREADS=$nt)"
      else
        bad "$prog failed with OMP_NUM_THREADS=$nt:"
        sed 's/^/          /' "$WORK/$prog.$nt.out" | grep -vE '^\s+#|Error termination' | head -15
      fi
    done
    if ! diff -q "$WORK/$prog.1.out" "$WORK/$prog.8.out" >/dev/null 2>&1; then
      bad "$prog output differs between 1 and 8 threads (race in the metric refresh?)"
    fi
  else
    bad "could not build $prog"
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "all MORPH guards passed"
else
  echo "$fail MORPH guard(s) FAILED"
  exit 1
fi
