#!/bin/bash
# run_case.sh <case dir> <work dir> <FUSS binary> [--iters N] [--verify] [--mpi N] [--none-inert]
#             [--dt-face-k-inert]
# Runs a test case in a fresh copy (the source tree is not written), 1 thread.
#   --iters N     stop after N iterations (iter-threshold and every level<m>-iter)
#   --verify      then run the case's verify.py
#   --mpi N       run serial and on N ranks, and require identical output files
#   --none-inert  run again with every FAS key spelled out under cycle = none, and require
#                 identical output files: the FAS keys must be inert outside a cycle
#   --dt-face-k-inert  run again with dt-face-k = true and require identical output files: on one
#                 material at constant k the face conductivities equal the cell's, so must the dt
set -u
case_dir=$(cd "$1" && pwd); work=$2; bin=$3; shift 3
iters=""; verify=0; np=0; none_inert=0; dtk_inert=0
while [[ $# -gt 0 ]]; do
  case $1 in
    --iters) iters=$2; shift 2 ;;
    --verify) verify=1; shift ;;
    --mpi) np=$2; shift 2 ;;
    --none-inert) none_inert=1; shift ;;
    --dt-face-k-inert) dtk_inert=1; shift ;;
    *) echo "run_case.sh: unknown option $1"; exit 2 ;;
  esac
done
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg
ulimit -s unlimited 2>/dev/null

prepare() {  # prepare <dir>: copy of the case with a clean OUTPUT and the iteration cap
  rm -rf "$1"; mkdir -p "$1"
  cp -r "$case_dir/." "$1/"; rm -rf "$1/OUTPUT" "$1/bin"; mkdir -p "$1/OUTPUT"
  if [[ -n $iters ]]; then
    sed -i -E '/^(iter-threshold|level[0-9]+-iter)[[:space:]]*=/d' "$1/input.ini"
    sed -i -E "s/^\[FUSS-Parameters\]/[FUSS-Parameters]\niter-threshold = $iters/" "$1/input.ini"
    sed -i -E "s/^levels[[:space:]]*=[[:space:]]*([0-9]+)/levels = \1\nLEVELS_ITER/" "$1/input.ini"
    local nl=$(sed -n -E 's/^levels[[:space:]]*=[[:space:]]*([0-9]+).*/\1/p' "$1/input.ini")
    local lines=""; for ((m=1; m<=${nl:-0}; m++)); do lines+="level$m-iter = $iters\n"; done
    sed -i "s/^LEVELS_ITER$/${lines%\\n}/" "$1/input.ini"
  fi
}

run() {  # run <dir> <ranks>: the solver must end with status 0 and write a residual history
  (cd "$1" && if [[ $2 -gt 1 ]]; then mpirun -np $2 "$bin"; else "$bin"; fi) > "$1/run.log" 2>&1
  local st=$?
  if [[ $st -ne 0 ]]; then echo "FAIL: FUSS exit status $st ($1/run.log)"; tail -n 5 "$1/run.log"; exit 1; fi
  [[ -s "$1/OUTPUT/residual-history.dat" ]] || { echo "FAIL: no residual history in $1/OUTPUT"; exit 1; }
}

prepare "$work"
run "$work" 1
if [[ $np -gt 1 ]]; then
  prepare "$work-mpi"
  run "$work-mpi" "$np"
  for f in $(cd "$work/OUTPUT" && ls); do
    cmp -s "$work/OUTPUT/$f" "$work-mpi/OUTPUT/$f" || { echo "FAIL: $f differs between 1 process and $np ranks"; exit 1; }
  done
  echo "1 process and $np ranks: identical output files"
fi
if [[ $none_inert -eq 1 ]]; then
  prepare "$work-none"
  sed -i -E '/^(cycle|pre-iter|post-iter|coarse-iter)[[:space:]]*=/d' "$work-none/input.ini"
  sed -i -E 's/^\[FUSS-Multigrid\]/[FUSS-Multigrid]\ncycle = none\npre-iter = 7\npost-iter = 5\ncoarse-iter = 3/' "$work-none/input.ini"
  run "$work-none" 1
  for f in $(cd "$work/OUTPUT" && ls); do
    cmp -s "$work/OUTPUT/$f" "$work-none/OUTPUT/$f" || { echo "FAIL: $f differs when the FAS keys are set under cycle = none"; exit 1; }
  done
  echo "FAS keys inert under cycle = none: identical output files"
fi
if [[ $dtk_inert -eq 1 ]]; then
  prepare "$work-dtk"
  sed -i -E '/^dt-face-k[[:space:]]*=/d' "$work-dtk/input.ini"
  sed -i -E 's/^\[FUSS-Numerics\]/[FUSS-Numerics]\ndt-face-k = true/' "$work-dtk/input.ini"
  run "$work-dtk" 1
  for f in $(cd "$work/OUTPUT" && ls); do
    cmp -s "$work/OUTPUT/$f" "$work-dtk/OUTPUT/$f" || { echo "FAIL: $f differs with dt-face-k = true"; exit 1; }
  done
  echo "dt-face-k inert on one material at constant k: identical output files"
fi
if [[ $verify -eq 1 ]]; then
  (cd "$work" && python3 -B verify.py) || exit 1
fi
echo "PASS"
