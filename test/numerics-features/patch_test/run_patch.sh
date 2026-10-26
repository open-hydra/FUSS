#!/bin/bash
# run_patch.sh <patch_test dir> <work dir> <FUSS binary>
# Linear-field patch test of the diffusive fluxes (see patch_case.py): one forward-Euler step from the exact field
# T = A + g.x on non-smooth meshes, the flux sum of every cell and the wall output must vanish to round-off.
#   part 1, default metric (stencil): 8 families (para rand jump fold2 sym wallmix box3d step), every group <= 1e-11
#   part 2, metric-mean (diffusive-metric = mean): para exact (<= 1e-11), and `step` (one jump of the cell size 1.5)
#           shows the bias (1+r)^2/(4r) - 1 of the mean metric: 6.4201e-3 of k|g|A_cell in the interior group, +- 1e-7
set -u
src=$(cd "$1" && pwd); work=$2; bin=$3
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg
ulimit -s unlimited 2>/dev/null
rm -rf "$work"; mkdir -p "$work"
fail=0
for fam in para rand jump fold2 sym wallmix box3d step; do
  w=$work/$fam
  python3 "$src/patch_case.py" case $fam "$w" > /dev/null || exit 1
  (cd "$w" && "$bin" > run.log 2>&1) || { echo "FAIL: $fam: solver exit status $?"; tail -n 5 "$w/run.log"; exit 1; }
  python3 "$src/patch_case.py" check "$w" 1e-11 | tee "$w/check.txt" | grep -v ' ok$' && fail=1
done
echo "part 1 (stencil): $( [[ $fail -eq 0 ]] && echo 'all groups <= 1e-11' || echo FAIL )"
mfail=0
for fam in para step; do
  w=$work/mean-$fam
  python3 "$src/patch_case.py" case $fam "$w" --metric mean > /dev/null || exit 1
  (cd "$w" && "$bin" > run.log 2>&1) || { echo "FAIL: mean $fam: solver exit status $?"; tail -n 5 "$w/run.log"; exit 1; }
  python3 "$src/patch_case.py" check "$w" 1e-11 > "$w/check.txt"
  if [[ $fam == para ]]; then
    grep -q FAIL "$w/check.txt" && { grep FAIL "$w/check.txt"; mfail=1; }
  else
    v=$(awk '$1=="RESID" && $3=="interior" {print $4}' "$w/check.txt")
    awk -v v="$v" 'BEGIN { d = v - 6.4201e-3; if (d < 0) d = -d; exit !(d <= 1e-7) }' \
      || { echo "FAIL: mean step interior $v, expected 6.4201e-3 +- 1e-7"; mfail=1; }
    echo "mean step interior: $v (expected 6.4201e-3)"
  fi
done
echo "part 2 (metric-mean): $( [[ $mfail -eq 0 ]] && echo PASS || echo FAIL )"
[[ $fail -eq 0 && $mfail -eq 0 ]] && echo PASS || exit 1
