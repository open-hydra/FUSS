#!/bin/bash
# MPIRuntimeInput: rank 0 reads input.ini and passes its text to every rank, at set-up and at every reading during the run.
# Two ranks run from two folders (MPMD): the folder of rank 1 holds the case without input.ini, or an input.ini with other
# values, or one that is rewritten during the run. The run must write the files of the run of one folder on two ranks,
# byte for byte, and end (a rank that read a different sol-diter, ini-diter or res-diter would call its collectives in
# another order and hang).
#   usage: mpi-runtime-input.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

base() {  # base <name>: nominal (two blocks), 600 iterations, numbered solutions every 100, input.ini read at 400
  local w="$WORK/$1"
  prep numerics-features/residual_convergence/nominal "$w"
  setkey "$w" iter-threshold 600; setkey "$w" res-threshold 1e-30
  setkey "$w" sol-diter 100 FUSS-IO; setkey "$w" sol-overwrite false FUSS-IO; setkey "$w" ini-diter 400 FUSS-IO
}
# mpmd <dirA> <dirB>: rank 0 in A, rank 1 in B; the status in A/exit.code
mpmd() {
  ( timeout 120 mpirun -n 1 -wdir "$1" "$BIN" : -n 1 -wdir "$2" "$BIN" > "$1/run.log" 2>&1; echo $? > "$1/exit.code" )
}

# reference: one folder, two ranks
base ref; run "$WORK/ref" 2; expect_ok ref "$WORK/ref"

# 1. rank 1 has no input.ini: the run is the reference
base missing-a; cp -a "$WORK/missing-a" "$WORK/missing-b"; rm "$WORK/missing-b/input.ini"
mpmd "$WORK/missing-a" "$WORK/missing-b"
n=$(ndiff "$WORK/ref" "$WORK/missing-a"); [[ $(code "$WORK/missing-a") == 0 && $n == 0 ]] && case_result rank1-missing 0 "same files as the one-folder run" || case_result rank1-missing 1 "exit $(code "$WORK/missing-a"), $n files differ"

# 2. rank 1 has an input.ini with other values (sol-diter, ini-diter, iter-threshold): rank 0's values are the run's
base other-a; cp -a "$WORK/other-a" "$WORK/other-b"
setkey "$WORK/other-b" sol-diter 7 FUSS-IO; setkey "$WORK/other-b" ini-diter 13 FUSS-IO; setkey "$WORK/other-b" iter-threshold 40
mpmd "$WORK/other-a" "$WORK/other-b"
n=$(ndiff "$WORK/ref" "$WORK/other-a"); [[ $(code "$WORK/other-a") == 0 && $n == 0 ]] && case_result rank1-other 0 "same files as the one-folder run" || case_result rank1-other 1 "exit $(code "$WORK/other-a"), $n files differ"

# 3. the input.ini of rank 1 is rewritten during the run (sol-diter 7): no effect, no hang
base during-a; cp -a "$WORK/during-a" "$WORK/during-b"
( cd "$WORK/during-a" && timeout 120 mpirun -n 1 -wdir "$WORK/during-a" "$BIN" : -n 1 -wdir "$WORK/during-b" "$BIN" > run.log 2>&1; echo $? > exit.code ) &
pid=$!
until [[ -f $WORK/during-a/OUTPUT/field1.tec ]] || ! kill -0 $pid 2>/dev/null; do sleep 0.01; done
( cd "$WORK/during-b" && sed 's|^sol-diter.*|sol-diter = 7|' input.ini > input.new && mv input.new input.ini )
wait $pid
n=$(ndiff "$WORK/ref" "$WORK/during-a"); [[ $(code "$WORK/during-a") == 0 && $n == 0 ]] && case_result rank1-rewritten 0 "same files, no hang" || case_result rank1-rewritten 1 "exit $(code "$WORK/during-a"), $n files differ"

finish
