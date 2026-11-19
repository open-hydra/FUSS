#!/bin/bash
# run_split.sh <connection_split dir> <work dir> <FUSS binary>
# ConnectionSplit: the rand mesh as one block and as two blocks joined by a block connection (101), 50 RK3
# iterations: the fields agree to 1e-12 (relative), and the two orders of the connection records give the same bits.
set -u
src=$(cd "$1" && pwd); work=$2; bin=$3
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg
ulimit -s unlimited 2>/dev/null
rm -rf "$work"; mkdir -p "$work"
for k in one two two-reversed; do
  python3 "$src/make_split.py" case $k "$work/$k" || exit 1
  (cd "$work/$k" && "$bin" > run.log 2>&1) || { echo "FAIL: $k: solver exit status $?"; tail -n 5 "$work/$k/run.log"; exit 1; }
done
python3 "$src/compare_split.py" "$work/one" "$work/two" 1e-12 || { echo "FAIL: two blocks differ from one block"; exit 1; }
for f in $(cd "$work/two/OUTPUT" && ls); do
  cmp -s "$work/two/OUTPUT/$f" "$work/two-reversed/OUTPUT/$f" || { echo "FAIL: $f differs between the two record orders"; exit 1; }
done
echo "record orders: identical output files"
echo PASS
