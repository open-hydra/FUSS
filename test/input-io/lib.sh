#!/bin/bash
# Helpers of the input-io tests (input checks, stops, restart, decomposition map). Source it from a test script:
#   BIN=<FUSS binary>   WORK=<work folder>   FUSS_CASES=<folder with steady-state/, transient/, numerics-features/>
# The test cases are copied (cp -rL) into WORK; the source tree is never written. One thread, one core.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES="${FUSS_CASES:-$(cd "$HERE/.." && pwd)}"
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg LC_ALL=C.UTF-8 LANG=C.UTF-8 PYTHONIOENCODING=utf-8 I_MPI_PIN=off
ulimit -s unlimited 2>/dev/null || true
NFAIL=0

# prep <case> <dir>: a copy of the case (path relative to test/) with an empty OUTPUT
prep() {
  rm -rf "$2"; mkdir -p "$2"; cp -rL "$CASES/$1/." "$2/"; rm -rf "$2/OUTPUT" "$2/bin"; mkdir -p "$2/OUTPUT"
}

# setkey <dir> <key> <value> [section]: replace a key of input.ini, or add it under the section (default FUSS-Parameters)
setkey() {
  local w=$1 k=$2 v=$3 s=${4:-FUSS-Parameters}
  if grep -qE "^$k[[:space:]]*=" "$w/input.ini"; then sed -i -E "s|^$k[[:space:]]*=.*|$k = $v|" "$w/input.ini"
  else sed -i -E "s|^\[$s\]|[$s]\n$k = $v|" "$w/input.ini"; fi
}

# delkey <dir> <key>: remove a key of input.ini
delkey() { sed -i -E "/^$2[[:space:]]*=/d" "$1/input.ini"; }

# run <dir> [ranks]: run the binary in the folder (mpirun for more than one rank); status in <dir>/exit.code, log in run.log
run() {
  local w=$1 np=${2:-1}
  if [[ $np -gt 1 ]]; then (cd "$w" && mpirun -np "$np" "$BIN" > run.log 2>&1); echo $? > "$w/exit.code"
  else (cd "$w" && "$BIN" > run.log 2>&1); echo $? > "$w/exit.code"; fi
  return 0
}
code() { cat "$1/exit.code"; }

# Report: CASE <name> PASS|FAIL <detail>; the test fails if one case does
case_result() { # case_result <name> <0 pass | 1 fail> <detail>
  if [[ $2 -eq 0 ]]; then echo "CASE $1 PASS $3"; else echo "CASE $1 FAIL $3"; NFAIL=$((NFAIL+1)); fi
}

# expect_stop <name> <dir> <regex>: the run in <dir> ended with exit code 1, the log matches the regex, and no solution was written
expect_stop() {
  local n=$1 w=$2 re=$3 c; c=$(code "$w")
  if [[ $c != 1 ]]; then case_result "$n" 1 "exit code $c, expected 1 ($w/run.log)"; return; fi
  if ! grep -qE -- "$re" "$w/run.log"; then case_result "$n" 1 "message '$re' not in $w/run.log"; return; fi
  if ls "$w"/OUTPUT/field* > /dev/null 2>&1; then case_result "$n" 1 "a solution was written"; return; fi
  case_result "$n" 0 "exit 1, '$re'"
}

# expect_ok <name> <dir>: exit code 0
expect_ok() { local c; c=$(code "$2"); if [[ $c == 0 ]]; then case_result "$1" 0 "exit 0"; else case_result "$1" 1 "exit code $c ($2/run.log)"; fi; }

# ndiff <dirA> <dirB>: number of files of OUTPUT that differ between the two runs, or are in one only (prints the names)
ndiff() {
  local a=$1/OUTPUT b=$2/OUTPUT n=0 f
  while read -r f; do
    if [[ ! -f $b/$f ]] || ! cmp -s "$a/$f" "$b/$f"; then n=$((n+1)); echo "  differs: $f" >&2; fi
  done < <(cd "$a" && find . -type f | sort)
  while read -r f; do [[ -f $a/$f ]] || { n=$((n+1)); echo "  only in $b: $f" >&2; }; done < <(cd "$b" && find . -type f | sort)
  echo $n
}

finish() { if [[ $NFAIL -gt 0 ]]; then echo "FAIL: $NFAIL case(s)"; exit 1; fi; echo "PASS"; exit 0; }

mkdir -p "$WORK"
