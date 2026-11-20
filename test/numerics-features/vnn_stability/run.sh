#!/usr/bin/env bash
#===============================================================================
#  VnnStability: the stability limit of the time step, vnn = 0.5 for forward Euler and 0.628 for RK3
#
#  run.sh <FUSS binary> <work dir>: on 16 x 16 x 16 cells of constant properties, walls 302 at 1000 K, T = 1000 K plus a random
#  field of +- 50 K at the start (the highest modes of the cells are in it), the local time step of Mod_dt (the diffusion limit
#  vnn, with the factor (2 + m + 1/m)/4 on skewed cells), serial, one thread, no IRS, res-threshold 1e-30 and a budget of N = 500
#  iterations (FUSS_VNN_ITER). The error of a mode grows by |1 - 4 vnn| per Euler step on the highest mode of the lattice, and
#  by the third-order polynomial of RK3 (amplitude 1 at 0.628), so just under the limit the highest modes do not decay and
#  the run cannot converge in N; just over it they grow by 2 - 8 % per iteration and the state check stops the run early:
#    cart    Cartesian cells
#    skew    rhombohedral cells, edges 115 degrees apart in the three pairs (as numerics-features/skewed_3D)
#  and for each of the two meshes
#    euler 0.50  and  RK3 0.60   must reach the budget (status `cap`): exit 0, N iterations, a finite field within
#                                twice the amplitude of the start (the noise must not grow)
#    euler 0.51  and  RK3 0.64   must stop (status `diverged`): exit 1, "Integration failed" in the log, before N, no solution
#  A run that ends before N with exit 0 has converged (status `converged`) and is a failure of the test, as is any other exit.
#  Over the four limits the stops of the campaign were at iterations 80 - 250 (cart RK3 0.64: 80, skew euler 0.51: 252).
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN=$(readlink -f "${1:?usage: run.sh <FUSS binary> <work dir>}")
WORK=${2:?usage: run.sh <FUSS binary> <work dir>}
CAP=${FUSS_VNN_ITER:-500}
MC="python3 -I -B $HERE/../../verification/convergence-order/make_case.py"
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg LC_ALL=C.UTF-8 LANG=C.UTF-8 PYTHONIOENCODING=utf-8
ulimit -s unlimited 2>/dev/null || true

rm -rf "$WORK"; mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
status=0
fail() { echo "  FAIL: $*"; status=1; }

# name mesh scheme vnn expected
RUNS=("cart-euler-0.50 cart euler 0.50 cap"  "cart-euler-0.51 cart euler 0.51 diverged"
      "cart-rk3-0.60   cart RK3   0.60 cap"  "cart-rk3-0.64   cart RK3   0.64 diverged"
      "skew-euler-0.50 skew115 euler 0.50 cap" "skew-euler-0.51 skew115 euler 0.51 diverged"
      "skew-rk3-0.60   skew115 RK3 0.60 cap"  "skew-rk3-0.64   skew115 RK3 0.64 diverged")

for r in "${RUNS[@]}"; do
  read -r name mesh scheme vnn want <<< "$r"
  (
    d="$WORK/$name"
    $MC case "$d" --sol const --dim 3 --mesh "$mesh" --n 16 --scheme "$scheme" --vnn "$vnn" --ic noise --noise 50 \
        --cap "$CAP" --res 1e-30 --resditer 1 --shellditer 100 > "$d.gen.log" 2>&1 || { echo 99 > "$d.exit"; exit; }
    ( cd "$d" && "$BIN" > run.log 2>&1 ); echo $? > "$d.exit"
  ) &
done
wait

for r in "${RUNS[@]}"; do
  read -r name mesh scheme vnn want <<< "$r"
  d="$WORK/$name"
  ex=$(cat "$d.exit" 2> /dev/null || echo 98)
  line=$($MC classify "$d" "$ex" "$CAP")
  got=$(sed -E 's/.* status=([^ ]+).*/\1/' <<< "$line")
  it=$(sed -E 's/.* iterations=([^ ]+).*/\1/' <<< "$line")
  dev=$(sed -E 's/.* max-dev=([^ ]+).*/\1/' <<< "$line")
  fin=$(sed -E 's/.* finite=([^ ]+).*/\1/' <<< "$line")
  printf '    %-16s %-6s vnn %s  exit %s  %s\n' "$name" "$scheme" "$vnn" "$ex" "${line#RUN }"
  [[ $got == "$want" ]] || { fail "$name: the run ended as '$got' (exit $ex), expected '$want' ($d/run.log)"; continue; }
  if [[ $want == cap ]]; then
    [[ $ex -eq 0 && $it -eq $CAP ]] || fail "$name: exit $ex after $it iterations, expected exit 0 after $CAP"
    [[ $fin == True ]] || fail "$name: the field is not finite"
    awk -v d="$dev" 'BEGIN { exit !(d + 0 <= 50.0) }' || fail "$name: the noise grew to $dev K (start: 50 K)"
  else
    [[ $ex -eq 1 ]] || fail "$name: exit status $ex, expected 1"
    grep -q 'Integration failed' "$d/run.log" || fail "$name: no 'Integration failed' message in $d/run.log"
    [[ $it -lt $CAP ]] || fail "$name: stopped only at the budget"
    [[ ! -e "$d/OUTPUT/field.tec" ]] || fail "$name: a solution was written by a run that stopped"
  fi
done

if [[ $status -ne 0 ]]; then
  echo "VnnStability: FAIL"
  exit 1
fi
echo "VnnStability: ok: euler 0.50 and RK3 0.60 reach $CAP iterations with a finite field, euler 0.51 and RK3 0.64 stop (exit 1), Cartesian and skewed cells"
exit 0
