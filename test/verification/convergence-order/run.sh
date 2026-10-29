#!/usr/bin/env bash
#===============================================================================
#  ConvergenceOrder and MultigridSolution (verification tier)
#
#  run.sh <FUSS binary> <work dir> [order|fas]: the part given (both with none) writes its decks with make_case.py (no ATLAS)
#  into <work dir>, solves each to res-threshold 1e-10 (serial, one thread, default RK3, vnn 0.4, diffusive-metric = stencil)
#  and compares the cell values with the exact solution (volume-weighted RMS and maximum over the cells).
#
#  order   the observed order of the discretisation on meshes that are not smooth: six cases, each at N = 16, 32 and 64
#          cells per side, in parallel; a case passes when the finest pair (32 -> 64) shows p(L2) >= 1.8 and p(Linf) >= 1.6
#          (ORDER_L2, ORDER_LINF). The order of the wall flux is not asserted (p(L2) of qw 0.8 - 1.0 on rand, 1.4 - 1.6 on jump
#          and shear, measured 2026-10-29 with `make_case.py error` on the runs of this test).
#            pois2-rand   Poisson T = 1000 + 100 sin(pi x) sin(pi y), source in INPUT/st.tec; interior nodes moved by 0.25 of
#                         the cell size, hash-seeded, scaled with h; 2 x 2 blocks joined by connections (they change nothing in
#                         the error: the one-block errors of the campaign are the same to 7 digits)
#            pois2-jump   the same, cells 4 times wider beyond x = 0.5, sheared x += 0.3 y: a jump of the cell size inside one
#                         block (diffusive-metric = mean is first order here)
#            pois2-shear  the same on a plate sheared by x += 0.5 y, 2 x 2 blocks
#            harm2-*      harmonic T = 1000 + 100 sin(pi x) sinh(pi y) / sinh(pi), the same three meshes (no source)
#  fas     the multigrid solution is the single-grid one: the harm2-rand and pois2-jump plates at N = 32 on 3 levels, FAS W and
#          FAS W with IRS (irs-beta 0.5, vnn 1.0), 2 x 2 blocks, converge to res-threshold 1e-10 in a few hundred iterations,
#          with no fallback of the IRS or of the cycle, and their cell values differ from the single-grid run of the same mesh
#          by at most 1e-6 K (FUSS_FAS_DT; the campaign measured <= 7e-9 K on 40 plates).
#  With no part both run. The exit status is 0 when every check passes.
#  FUSS_ORDER_METRIC=mean runs the order part with diffusive-metric = mean, whose order is 1 or less on the random and jump
#  meshes: the test must then fail (the control that it can).
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN=$(readlink -f "${1:?usage: run.sh <FUSS binary> <work dir> [order|fas]}")
WORK=${2:?usage: run.sh <FUSS binary> <work dir> [order|fas]}
PARTS=${3:-"order fas"}
ORDER_L2=${ORDER_L2:-1.8}
ORDER_LINF=${ORDER_LINF:-1.6}
FAS_DT=${FUSS_FAS_DT:-1e-6}
METRIC=${FUSS_ORDER_METRIC:-stencil}     # `mean` is first order on random and jump meshes: the control that the test can fail
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg LC_ALL=C.UTF-8 LANG=C.UTF-8 PYTHONIOENCODING=utf-8
ulimit -s unlimited 2>/dev/null || true
MC="python3 -I -B $HERE/make_case.py"

rm -rf "$WORK"; mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
status=0
fail() { echo "  FAIL: $*"; status=1; }

# solve <dir> <make_case.py case options...>: write the deck, run the solver, require exit 0 and a residual at the threshold.
# Prints nothing; on a failure writes the reason to <dir>/FAILED (read by the caller after `wait`).
solve() {
  local d=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  $MC case "$d" "$@" > "$d/gen.log" 2>&1 || { echo "the deck was not written ($d/gen.log)" > "$d/FAILED"; return; }
  ( cd "$d" && "$BIN" > run.log 2>&1 ); local st=$?
  [[ $st -eq 0 ]] || { echo "the solver exited with $st ($d/run.log)" > "$d/FAILED"; return; }
  local last it
  last=$(tail -n 1 "$d/OUTPUT/residual-history.dat" | awk '{print $3}')
  it=$(tail -n 1 "$d/OUTPUT/residual-history.dat" | awk '{print $1}')
  awk -v r="$last" 'BEGIN { exit !(r + 0 <= 1.0000001e-10) }' \
    || { echo "not converged: last residual $last after $it iterations ($d/run.log)" > "$d/FAILED"; return; }
  echo "$it" > "$d/iterations"
  rm -rf "$d/INPUT" "$d/MESH"
}

# errors <dir>: "l2 linf" of the run, or nothing
errors() {
  local err
  err=$($MC error "$1" --no-walls 2> /dev/null) || return 1
  sed -E 's/.* l2=([^ ]+).* linf=([^ ]+).*/\1 \2/' <<< "$err"
}

#-------------------------------------------------------------------------------------------------- order
spec() {   # <case> -> "sol mesh blocks"
  case $1 in
    pois2-rand)  echo "pois2 rand 2,2,1" ;;   pois2-jump)  echo "pois2 jump 1,1,1" ;;   pois2-shear) echo "pois2 shear 2,2,1" ;;
    harm2-rand)  echo "harm2 rand 2,2,1" ;;   harm2-jump)  echo "harm2 jump 1,1,1" ;;   harm2-shear) echo "harm2 shear 2,2,1" ;;
  esac
}

if [[ " $PARTS " == *" order "* ]]; then
  CASES=(pois2-rand pois2-jump pois2-shear harm2-rand harm2-jump harm2-shear)
  NS="16 32 64"
  echo "ConvergenceOrder: ${#CASES[@]} cases x N = $NS, in parallel"
  for c in "${CASES[@]}"; do
    read -r sol mesh blocks <<< "$(spec "$c")"
    (
      for n in $NS; do
        # the budget is about twice the iterations the campaign needed (1.6e3, 6.3e3, 2.5e4 at N = 16, 32, 64)
        solve "$WORK/order/$c/N$n" --sol "$sol" --mesh "$mesh" --n "$n" --blocks "$blocks" --metric "$METRIC" \
              --cap $((15 * n * n)) --resditer 10 --shellditer 1000
        [[ -e "$WORK/order/$c/N$n/FAILED" ]] && break
      done
    ) &
  done
  wait
  for c in "${CASES[@]}"; do
    triples=(); ok=1
    for n in $NS; do
      d="$WORK/order/$c/N$n"
      if [[ -e "$d/FAILED" ]]; then fail "$c N=$n: $(cat "$d/FAILED")"; ok=0; break; fi
      e=$(errors "$d") || { fail "$c N=$n: no error line (the field was not written)"; ok=0; break; }
      read -r l2 linf <<< "$e"
      echo "    $c N=$n  L2 $l2  Linf $linf  ($(cat "$d/iterations") iterations)"
      triples+=("$n:$l2:$linf")
    done
    [[ $ok == 1 ]] || continue
    line=$($MC orders "$c" "$ORDER_L2" "$ORDER_LINF" "${triples[@]}") \
      || fail "$c: observed order of the finest pair below p(L2) $ORDER_L2 or p(Linf) $ORDER_LINF"
    echo "  $line"
  done
fi

#-------------------------------------------------------------------------------------------------- fas
if [[ " $PARTS " == *" fas "* ]]; then
  echo "MultigridSolution: single grid, FAS W and FAS W + IRS on 3 levels, N = 32, 2 x 2 blocks"
  FM=(harm2-rand pois2-jump)
  for c in "${FM[@]}"; do
    read -r sol mesh _ <<< "$(spec "$c")"
    base=(--sol "$sol" --mesh "$mesh" --n 32 --blocks 2,2,1 --resditer 1 --shellditer 100)
    ( solve "$WORK/fas/$c/single" "${base[@]}" --cap 15000 ) &
    ( solve "$WORK/fas/$c/w" "${base[@]}" --levels 3 --cycle w --cap 3000 ) &
    ( solve "$WORK/fas/$c/wirs" "${base[@]}" --levels 3 --cycle w --irs 0.5 --vnn 1.0 --cap 3000 ) &
  done
  wait
  for c in "${FM[@]}"; do
    d="$WORK/fas/$c"
    if [[ -e "$d/single/FAILED" ]]; then fail "$c single grid: $(cat "$d/single/FAILED")"; continue; fi
    for v in w wirs; do
      if [[ -e "$d/$v/FAILED" ]]; then fail "$c $v: $(cat "$d/$v/FAILED")"; continue; fi
      # no fallback: the IRS and the cycle report their fallbacks in the log ("IRS: n ...", "FAS: n prolongated / n coarse")
      if grep -qE '^(IRS|FAS): [1-9][0-9]* ' "$d/$v/run.log"; then
        fail "$c $v: a fallback was reported: $(grep -E '^(IRS|FAS): [1-9][0-9]* ' "$d/$v/run.log" | head -n 1)"; continue
      fi
      dt=$($MC diff "$d/$v" "$d/single") || { fail "$c $v: no field to compare"; continue; }
      echo "    $c $v: $(cat "$d/$v/iterations") iterations (single grid $(cat "$d/single/iterations")), max |T - T_single| = $dt K"
      awk -v x="$dt" -v m="$FAS_DT" 'BEGIN { exit !(x + 0 <= m + 0) }' \
        || fail "$c $v: the multigrid solution differs from the single-grid one by $dt K, more than $FAS_DT K"
    done
  done
fi

if [[ $status -ne 0 ]]; then
  echo "verification: FAIL"
  exit 1
fi
echo "verification: ok: ${PARTS}"
exit 0
