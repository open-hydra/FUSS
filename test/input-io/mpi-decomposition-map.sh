#!/bin/bash
# MPIDecompositionMap: with MPI the map of INPUT/ (key decomposition-map) decides which rank owns which block: the two blocks of
# nominal reversed on 2 ranks, and on 4 ranks (two ranks without a block idle, with a warning), give the files of the serial run
# byte for byte; a map written for 3 ranks, run on 2, stops every rank.
#   usage: mpi-decomposition-map.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
R=numerics-features/residual_convergence/nominal
mk() {  # mk <name> <rank of block 1> <rank of block 2> [no key]
  prep $R "$WORK/$1"; setkey "$WORK/$1" iter-threshold 50
  python3 "$(dirname "${BASH_SOURCE[0]}")/write_map.py" "$WORK/$1/INPUT/ic.tec" "$WORK/$1/INPUT/decomposition.map" "$2" "$3"
  [[ -z ${4:-} ]] && setkey "$WORK/$1" decomposition-map decomposition.map FUSS-IO
}
mk serial 0 0 nokey; run "$WORK/serial"; expect_ok serial "$WORK/serial"
mk rev2 1 0; run "$WORK/rev2" 2
n=$(ndiff "$WORK/serial" "$WORK/rev2"); if [[ $(code "$WORK/rev2") == 0 && $n == 0 ]] && grep -q 'partition read from' "$WORK/rev2/run.log"; then case_result reversed-2-ranks 0 "same files as the serial run"; else case_result reversed-2-ranks 1 "exit $(code "$WORK/rev2"), $n files differ"; fi
mk four 3 1; run "$WORK/four" 4
n=$(ndiff "$WORK/serial" "$WORK/four"); if [[ $(code "$WORK/four") == 0 && $n == 0 ]] && grep -q 'gives no block to 2 of 4 ranks' "$WORK/four/run.log"; then case_result idle-ranks-4 0 "same files, warning"; else case_result idle-ranks-4 1 "exit $(code "$WORK/four"), $n files differ, or no warning"; fi
mk three 0 2; run "$WORK/three" 2
if [[ $(code "$WORK/three") == 1 ]] && grep -qE '\[RANK +0\] STOP' "$WORK/three/run.log" && grep -qE '\[RANK +1\] STOP' "$WORK/three/run.log" && grep -q 'outside the ranks of the run' "$WORK/three/run.log"; then case_result map-3-ranks-on-2 0 "every rank stops"; else case_result map-3-ranks-on-2 1 "exit $(code "$WORK/three")"; fi
finish
