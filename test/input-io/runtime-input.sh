#!/bin/bash
# RuntimeInput: input.ini read again during the run (every ini-diter iterations) keeps to the rules of the registry.
# An error of the reading does what ini-error says: with warn (the default) the reading is not taken, one warning gives
# the iteration and the error, and the run goes on with every value it had (residual history byte for byte that of the run
# without the change); with stop the run stops with the message and exit code 1, as at set-up.
# Base: residual_convergence/nominal (about 1 ms per iteration), 2000 iterations, numbered solutions every 500, input.ini
# read again at iteration 1500; the change is written once the solution of iteration 500 exists (a host too slow to poll
# makes the test say so and fail).
#   usage: runtime-input.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

base() {  # base <name> [ini-diter]
  local w="$WORK/$1"
  prep numerics-features/residual_convergence/nominal "$w"
  setkey "$w" iter-threshold 2000; setkey "$w" res-threshold 1e-30
  setkey "$w" sol-diter 500 FUSS-IO; setkey "$w" sol-overwrite false FUSS-IO; setkey "$w" ini-diter "${2:-1500}" FUSS-IO
}

# change <name> <shell command>: run the case; once OUTPUT/field1.tec exists, apply the command (it edits input.ini)
change() {
  local w="$WORK/$1" cmd=$2 pid
  ( cd "$w" && "$BIN" > run.log 2>&1; echo $? > exit.code ) &
  pid=$!
  until [[ -f $w/OUTPUT/field1.tec ]] || ! kill -0 $pid 2>/dev/null; do sleep 0.02; done
  if [[ -f $w/OUTPUT/field1.tec ]]; then ( cd "$w" && eval "$cmd" ); fi
  wait $pid
}
# atomic rewrite of input.ini through a copy
edit() { ( cd "$1" && cp input.ini input.new && eval "$2" && mv input.new input.ini ); }

# control: the case without a change; and a run that starts with a rule broken stops at set-up
base control; run "$WORK/control"; expect_ok control "$WORK/control"
base start-bad; setkey "$WORK/start-bad" vnn -1 FUSS-Numerics; run "$WORK/start-bad"; expect_stop start-bad "$WORK/start-bad" 'vnn must be > 0'

# ini-diter reads that change nothing: the run reads input.ini every 100 iterations and writes the files of the control
base unchanged 100; run "$WORK/unchanged"
n=$(ndiff "$WORK/control" "$WORK/unchanged"); [[ $(code "$WORK/unchanged") == 0 && $n == 0 ]] && case_result unchanged 0 "same files" || case_result unchanged 1 "exit $(code "$WORK/unchanged"), $n files differ"

# ini-error = warn: the reading is not taken, the run goes on with its values
warn_case() {  # warn_case <name> <edit command> <ini-error at the start>
  base "$1"; [[ -n ${3:-} ]] && setkey "$WORK/$1" ini-error "$3" FUSS-IO
  change "$1" "$2"
  local w="$WORK/$1" n; n=$(ndiff "$WORK/control" "$w")
  local nn; nn=$(grep -c 'NOT TAKEN' "$w/run.log")
  if [[ $(code "$w") != 0 ]]; then case_result "$1" 1 "exit code $(code "$w")"
  elif [[ $nn != 1 ]]; then case_result "$1" 1 "NOT TAKEN written $nn times (a host too slow to poll?)"
  elif grep -q '\[ERROR\]' "$w/run.log"; then case_result "$1" 1 "an [ERROR] line"
  elif [[ $(cmp "$WORK/control/OUTPUT/residual-history.dat" "$w/OUTPUT/residual-history.dat" 2>&1 | wc -l) != 0 ]]; then case_result "$1" 1 "residual history differs from the run without the change"
  else case_result "$1" 0 "not taken once, residual history unchanged"; fi
}
warn_case warn-vnn       'sed -i -E "s|^vnn.*|vnn = -1|" input.ini'
warn_case warn-empty     ': > input.ini'
warn_case warn-comments  'printf "; nothing\n" > input.ini'
warn_case warn-bad-value 'sed -i -E "s|^vnn.*|vnn = 0.4|; s|^sol-diter.*|sol-diter = seven|" input.ini'
# the value of ini-error in force when the reading starts decides, also when the file writes another one
warn_case warn-overruled 'sed -i -E "s|^vnn.*|vnn = -1|; s|^ini-error.*|ini-error = stop|" input.ini'

# ini-error = stop: the run stops at the reading of iteration 1500, with the message and exit code 1
stop_case() {  # stop_case <name> <edit command> <message regex>
  base "$1"; setkey "$WORK/$1" ini-error stop FUSS-IO; change "$1" "$2"
  local w="$WORK/$1"
  if [[ $(code "$w") != 1 ]]; then case_result "$1" 1 "exit code $(code "$w"), expected 1 (a host too slow to poll?)"
  elif ! grep -qE -- "$3" "$w/run.log"; then case_result "$1" 1 "message '$3' not in run.log"
  elif grep -q 'NOT TAKEN' "$w/run.log"; then case_result "$1" 1 "NOT TAKEN with ini-error = stop"
  elif [[ -f $w/OUTPUT/field4.tec ]]; then case_result "$1" 1 "the run went on past the reading"
  else case_result "$1" 0 "stopped, '$3'"; fi
}
stop_case stop-vnn   'sed -i -E "s|^vnn.*|vnn = -1|" input.ini'  'vnn must be > 0'
stop_case stop-empty ': > input.ini'                             'input.ini.*is empty'

finish
