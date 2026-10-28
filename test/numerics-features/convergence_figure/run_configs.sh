#!/bin/bash
# run_configs.sh <runs dir> [-c <case,...>] [-k <config,...>] [-j <parallel runs>]
#
# Residual convergence of none / IRS / V / W / V + IRS / W + IRS on four cases of test/,
# the runs behind the figures of docs/vv/fas-cycles.md (then: python3 plot_convergence.py <runs dir>).
#   plate        residual_convergence/fas_v   two-material plate, 2 blocks (101), 2 levels
#   liner        multimat_liner               four-layer liner wall, 2 blocks (101), 4 levels
#   chimera_mg   chimera_mg                   two overlapping blocks (102), 3 levels
#   anisotropic  anisotropic_Plate            cells of aspect ratio 100, 3 levels
# Every run: RK3, cons, res-threshold 1e-10, residual every iteration, one thread. Without IRS
# vnn 0.5; with IRS vnn 1.0 and irs-beta 0.5. none and IRS run on one grid (iter-threshold
# 600000); the cycles on the case's levels with a 1-iteration descent on each coarse level
# (level1-iter 100000). The keys are written into the copy with set_keys.py.
# Each run lands in <runs dir>/<case>_<config> with run.log, run.err and run.time
# ("exit <status> <seconds> s"). The long single-grid runs start first (the liner ones take
# about 40 and 30 minutes); -j runs that many at a time; -c and -k select cases and configurations
# (none irs v w v_irs w_irs), e.g. -c liner -k v,w,v_irs,w_irs to redo the liner cycles.
# The solver is bin/FUSS of this repository, or $FUSS: a path, or a command such as
# FUSS="taskset -c 0-3 /path/to/FUSS".
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)
FUSS=${FUSS:-$root/bin/FUSS}
cases=plate,liner,chimera_mg,anisotropic; cfgs=none,irs,v,w,v_irs,w_irs; jobs=1
[[ $# -ge 1 && $1 != -* ]] || { sed -n '2,19p' "$0"; exit 2; }
out=$1; shift
while [[ $# -gt 0 ]]; do
  case $1 in
    -c) cases=$2; shift 2 ;;
    -k) cfgs=$2; shift 2 ;;
    -j) jobs=$2; shift 2 ;;
    *) echo "run_configs.sh: unknown option $1"; exit 2 ;;
  esac
done
command -v ${FUSS%% *} > /dev/null || { echo "run_configs.sh: no solver at $FUSS (set FUSS=)"; exit 2; }
mkdir -p "$out"; out=$(cd "$out" && pwd)
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M FUSS here
ulimit -s unlimited 2>/dev/null

T=$root/test/numerics-features
C='FUSS-Parameters:res-threshold=1e-10;FUSS-Numerics:time-scheme=RK3;FUSS-Numerics:integration-variables=cons'
C="$C;FUSS-IO:res-diter=1;FUSS-IO:shell-diter=10000;FUSS-IO:sol-diter=100000000;FUSS-IO:ini-diter=100000000"
NOIRS='FUSS-Numerics:vnn=0.5;FUSS-Numerics:irs-beta=0'
IRS='FUSS-Numerics:vnn=1.0;FUSS-Numerics:irs-beta=0.5'
ONE='FUSS-Multigrid:levels=1;FUSS-Multigrid:cycle=none;FUSS-Parameters:iter-threshold=600000'
for m in 1 2 3 4; do ONE="$ONE;FUSS-Multigrid:level$m-iter="; done   # the case's level budgets go

run_one() {  # run_one <case dir> <run dir> <keys>: a copy of the case with the keys set, run serial
  local src=$1 d=$2 keys=$3 s st
  rm -rf "$d"; mkdir -p "$d"; cp -r "$src/." "$d/"; rm -rf "$d/OUTPUT" "$d/bin"; mkdir -p "$d/OUTPUT"
  python3 "$here/set_keys.py" "$d/input.ini" "$keys" || { echo "exit 2 0 s" > "$d/run.time"; return; }
  s=$(date +%s)
  (cd "$d" && $FUSS > run.log 2> run.err); st=$?
  echo "exit $st $(( $(date +%s) - s )) s" > "$d/run.time"
  echo "$(basename "$d"): exit $st, $(( $(date +%s) - s )) s, last residual line: $(tail -1 "$d/OUTPUT/residual-history.dat" 2>/dev/null)"
}
export -f run_one

: > "$out/jobs.txt"
for cfg in none irs v w v_irs w_irs; do       # single-grid runs first: they are the long ones
  [[ ",$cfgs," == *",$cfg,"* ]] || continue
  for spec in "liner|$T/multimat_liner|4" "plate|$T/residual_convergence/fas_v|2" \
              "anisotropic|$T/anisotropic_Plate|3" "chimera_mg|$T/chimera_mg|3"; do
    IFS='|' read -r name dir L <<< "$spec"
    [[ ",$cases," == *",$name,"* ]] || continue
    LV="FUSS-Multigrid:levels=$L;FUSS-Multigrid:level1-iter=100000;FUSS-Parameters:iter-threshold="
    for m in $(seq 2 "$L"); do LV="$LV;FUSS-Multigrid:level$m-iter=1"; done
    case $cfg in
      none)  keys="$C;$ONE;$NOIRS" ;;
      irs)   keys="$C;$ONE;$IRS" ;;
      v)     keys="$C;$LV;FUSS-Multigrid:cycle=v;$NOIRS" ;;
      w)     keys="$C;$LV;FUSS-Multigrid:cycle=w;$NOIRS" ;;
      v_irs) keys="$C;$LV;FUSS-Multigrid:cycle=v;$IRS" ;;
      w_irs) keys="$C;$LV;FUSS-Multigrid:cycle=w;$IRS" ;;
    esac
    echo "$dir|$out/${name}_$cfg|$keys" >> "$out/jobs.txt"
  done
done
echo "run_configs.sh: $(wc -l < "$out/jobs.txt") runs, $jobs at a time, solver $FUSS, start $(date)"
tr '|' '\0' < "$out/jobs.txt" | tr '\n' '\0' | xargs -0 -n 3 -P "$jobs" bash -c 'run_one "$@"' _
echo "run_configs.sh: done $(date)"
