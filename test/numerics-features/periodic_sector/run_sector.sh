#!/bin/bash
# run_sector.sh <periodic_sector dir> <work dir> <FUSS binary> [--mpi N]
# PeriodicSector: one sector of a helical annular wall with its side faces periodic (201), as two blocks (sector) and as
# one block periodic with itself (sector1), against the full ring of 6 sectors joined by block connections (101), cut in
# the same blocks (ring, ring1), for diffusive-metric = stencil and mean: the fields agree to 1e-11 (relative) after 100
# RK3 iterations, and symmetry (300) on the side faces does not. The periodic faces paired in another order (flip) stop
# the run at set-up. With --mpi N the sector runs again on N ranks (the periodic pair across two ranks) and must give the
# output files of the serial run.
set -u
src=$(cd "$1" && pwd); work=$2; bin=$3; shift 3
np=0
[[ ${1:-} == --mpi ]] && np=$2
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg
ulimit -s unlimited 2>/dev/null
rm -rf "$work"; mkdir -p "$work"

run() {  # run <dir> <ranks>
  (cd "$1" && if [[ $2 -gt 1 ]]; then mpirun -np $2 "$bin"; else "$bin"; fi) > "$1/run.log" 2>&1 \
    || { echo "FAIL: $1: solver exit status $?"; tail -n 5 "$1/run.log"; exit 1; }
}

for metric in stencil mean; do
  for k in ring ring1 sector sector1 sym; do
    python3 "$src/make_sector.py" case $k "$work/$metric-$k" --metric $metric || exit 1
    run "$work/$metric-$k" 1
  done
  grep -A7 'Periodic faces (201)' "$work/$metric-sector/run.log" | head -4
  python3 "$src/compare_sector.py" "$work/$metric-ring" "$work/$metric-sector"  1e-11 || { echo "FAIL: $metric: sector differs from the ring"; exit 1; }
  python3 "$src/compare_sector.py" "$work/$metric-ring1" "$work/$metric-sector1" 1e-11 || { echo "FAIL: $metric: sector1 differs from the ring"; exit 1; }
  python3 "$src/compare_sector.py" "$work/$metric-ring" "$work/$metric-sym" 1e-6 --differ || { echo "FAIL: $metric: symmetric sides give the ring"; exit 1; }
done

python3 "$src/make_sector.py" case flip "$work/flip" || exit 1
(cd "$work/flip" && "$bin" > run.log 2>&1) && { echo "FAIL: flip: the run did not stop"; exit 1; }
grep -q 'not images by a rigid motion' "$work/flip/run.log" || { echo "FAIL: flip: no misfit message"; tail -n 3 "$work/flip/run.log"; exit 1; }
echo "flip: stopped at set-up: $(grep -o 'not images by a rigid motion (misfit [^)]*)' "$work/flip/run.log" | head -1)"

if [[ $np -gt 1 ]]; then
  python3 "$src/make_sector.py" case sector "$work/stencil-sector-mpi" || exit 1
  run "$work/stencil-sector-mpi" "$np"
  for f in $(cd "$work/stencil-sector/OUTPUT" && ls); do
    cmp -s "$work/stencil-sector/OUTPUT/$f" "$work/stencil-sector-mpi/OUTPUT/$f" \
      || { echo "FAIL: $f differs between 1 process and $np ranks"; exit 1; }
  done
  echo "1 process and $np ranks: identical output files"
fi
echo PASS
