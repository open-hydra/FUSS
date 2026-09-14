#!/bin/bash
#===============================================================================
#         FILE: test/run_runtime_ini_gate.sh
#        USAGE: ./test/run_runtime_ini_gate.sh [-w WORKDIR]
#  DESCRIPTION: input.ini is re-read every `ini-diter` iterations WHILE THE RUN
#               IS IN PROGRESS. This gates what that re-read is allowed to do
#               (plan 10, D3).
#
#               It used to re-assign every registered parameter and then compute
#               Validate_Registry() into a variable it never looked at. So an
#               out-of-range value was accepted in silence, `newrun` and `law`
#               were re-assignable without the state derived from them being
#               re-derived, and -- worst -- editing `irs = true` into a running
#               moving-mesh case reached exactly the configuration that is
#               refused at startup, because the guards only run at setup.
#
#               The policy now: only runtime-mutable parameters are applied.
#               Anything else is IGNORED, the run continues on the value it
#               started with, and the attempt is appended to
#               OUTPUT/<prefix>runtime-ini-ignored.txt -- a file created ONLY
#               when something was ignored, so its existence is the signal.
#
#               HOW THESE TESTS WORK, since the defect only exists mid-flight:
#               each case starts a run configured to continue effectively
#               forever, waits until it is demonstrably running, edits
#               input.ini, and lets the run finish. `time-threshold` is
#               runtime-mutable, so lowering it is BOTH the positive test and
#               the mechanism that stops the run.
#
#               NOT `iter-threshold`, which is what this script used at first:
#               it is copied into domain%itermax at setup and the stopping test
#               reads itermax, so lowering it mid-run does nothing and the run
#               hung. That is why it is classed immutable -- see RUNTIME_MUTABLE
#               in Backend_INI.f90.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BASE="$ROOT/test/moving-mesh/receding/N040"
WORK="${TMPDIR:-/tmp}/fuss-runtime-ini.$$"

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR]"; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$WORK"
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS runtime input.ini gate"
echo "  work dir : $WORK"
echo

REPORT="OUTPUT/runtime-ini-ignored.txt"

# $1 = case name, $2.. = sed expressions applied to input.ini MID-RUN
start_case () {
  local name=$1; shift
  local d="$WORK/$name"

  rm -rf "$d"; mkdir -p "$d"
  cp -r "$BASE/INPUT" "$BASE/MESH" "$BASE/input.ini" "$BASE/FUSS.sh" "$d/"
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"

  # Run effectively forever, and re-read input.ini often.
  sed -i 's/^iter-threshold = .*/iter-threshold = 100000000/' "$d/input.ini"
  sed -i 's/^time-threshold = .*/time-threshold = 1e30/'      "$d/input.ini"
  sed -i 's/^ini-diter     = .*/ini-diter     = 20/'          "$d/input.ini"

  ( cd "$d" && ./FUSS.sh -p 1 solve > run.log 2>&1 ) &
  local runner=$!

  # Wait until it is demonstrably past the first re-read, so the edit lands
  # mid-flight rather than racing setup.
  local waited=0
  while [ ! -s "$d/OUTPUT/residual-history.dat" ] && [ $waited -lt 200 ]; do
    sleep 0.1; waited=$((waited+1))
  done
  sleep 0.5

  local e
  for e in "$@"; do sed -i "$e" "$d/input.ini"; done
  # Always stop it: time-threshold is read live every step, so this is also the
  # positive control -- if the mutable path were broken the run would hang here
  # and the test would time out rather than quietly pass.
  sed -i 's/^time-threshold = .*/time-threshold = 1e-9/' "$d/input.ini"

  wait $runner
  echo "$d"
}

# ---------------------------------------------------------------------------
d=$(start_case mutable_applies)
if grep -q "Time of operation" "$d/run.log"; then
  pass "a runtime-mutable change takes effect (time-threshold lowered, run stopped)"
else
  bad "mutable_applies: run did not terminate after time-threshold was lowered"
fi
if [ ! -f "$d/$REPORT" ]; then
  pass "no report file when only mutable parameters changed"
else
  bad "mutable_applies: a report file was written when nothing should have been ignored"
  sed 's/^/          /' "$d/$REPORT"
fi

# ---------------------------------------------------------------------------
d=$(start_case immutable_ignored 's/^vnn = .*/vnn = 0.123/')
if [ -f "$d/$REPORT" ] && grep -q "vnn" "$d/$REPORT"; then
  pass "an immutable change is ignored and named: $(grep -m1 vnn "$d/$REPORT" | sed 's/^ *//')"
else
  bad "immutable_ignored: vnn change was not reported"
fi
if grep -qE "Delta t = 0\.868029E-01" "$d/run.log"; then
  pass "...and the running value really was kept (dt unchanged)"
else
  bad "immutable_ignored: dt changed, so the edit was applied after all"
  grep -m1 "Delta t" "$d/run.log" | sed 's/^/          /'
fi

# ---------------------------------------------------------------------------
# The one that matters most: this combination is REFUSED at startup, so
# reaching it mid-run was a way around the guard entirely.
d=$(start_case guard_bypass 's/^irs = .*/irs = true/')
if [ -f "$d/$REPORT" ] && grep -q "irs" "$d/$REPORT"; then
  pass "enabling IRS mid-run on a moving mesh is ignored (the setup guard cannot be bypassed)"
else
  bad "guard_bypass: irs change was not reported"
fi

# ---------------------------------------------------------------------------
# A mutable parameter that IS read live must actually change behaviour, and
# must not be reported as ignored.
d=$(start_case mutable_not_reported 's/^res-threshold = .*/res-threshold = 1e-12/')
if [ ! -f "$d/$REPORT" ]; then
  pass "a mutable change is applied silently, with no report"
else
  bad "mutable_not_reported: res-threshold was treated as immutable"
  sed 's/^/          /' "$d/$REPORT"
fi

# ---------------------------------------------------------------------------
# A mutable parameter edited to an OUT-OF-RANGE value must be rejected, the old
# value kept, and the rejection reported. Validation is per parameter, through
# Validate_Param -- whole-registry validation cannot be re-run after setup.
d=$(start_case invalid_value 's/^res-threshold = .*/res-threshold = -1.0/')
if [ -f "$d/$REPORT" ] && grep -q "res-threshold must be" "$d/$REPORT" \
   && grep -q 'rejected "-1.0"' "$d/$REPORT"; then
  pass "an out-of-range mutable value is rejected: $(grep -m1 'res-threshold must' "$d/$REPORT" | sed 's/^ *//')"
else
  bad "invalid_value: the validation result was not acted on"
  [ -f "$d/$REPORT" ] && sed 's/^/          /' "$d/$REPORT"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "runtime ini gate passed"
else
  echo "$fail check(s) FAILED"
  exit 1
fi
