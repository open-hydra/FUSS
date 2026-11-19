#!/bin/bash
# MPIRestart: a run on two ranks split in two (N iterations, then a restart that runs N more) writes the files of the run in
# one go, byte for byte (vtk raw): the counters are written by the root and read by every rank. nominal (two blocks) and
# mg_irs (two blocks, two grid levels: the restart goes on at level 1).
#   usage: mpi-restart.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
R=numerics-features/residual_convergence

mpi_split() {  # mpi_split <name> <case> <N> [level2-iter]
  local n=$1 c=$2 N=$3 a="$WORK/$1-A" b="$WORK/$1-B" d key
  prep "$c" "$a"; prep "$c" "$b"
  for d in "$a" "$b"; do
    setkey "$d" res-threshold 1e-30; setkey "$d" sol-format "vtk raw" FUSS-IO
    [[ -n ${4:-} ]] && setkey "$d" level2-iter "$4" FUSS-Multigrid
  done
  if [[ -n ${4:-} ]]; then key="level1-iter FUSS-Multigrid"; else key="iter-threshold FUSS-Parameters"; fi
  setkey "$a" ${key% *} $((2*N)) ${key#* }; run "$a" 2
  setkey "$b" ${key% *} "$N" ${key#* }; run "$b" 2; setkey "$b" newrun false; run "$b" 2
  if [[ $(code "$a") != 0 || $(code "$b") != 0 ]]; then case_result "$n" 1 "exit A $(code "$a") B $(code "$b")"; return; fi
  local nd; nd=$(ndiff "$a" "$b")
  if [[ $nd == 0 ]]; then case_result "$n" 0 "0 differing files"; else case_result "$n" 1 "$nd differing files"; fi
}
mpi_split nominal $R/nominal 100
mpi_split mg      $R/mg_irs  100 300
finish
