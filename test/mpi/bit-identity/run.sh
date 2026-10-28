#!/usr/bin/env bash
#===============================================================================
#  The same output bits on any layout of MPI ranks and OpenMP threads
#  (needs-mpi tier)
#
#  usage: run.sh <FUSS binary> <work dir>
#
#  A run must write the same bits with one thread, with N threads, with N ranks
#  and with ranks and threads together. Two cases (make_case.py bit), each
#  solved on 1x1 (reference), 1x4, 2x1, 4x1 and 2x2 (MPI ranks x OpenMP threads);
#  every file the run writes in OUTPUT/ (field, diagnostic and wall solutions in
#  vtk raw, residual-history.dat) must be identical byte for byte to that of 1x1,
#  and the list of files the same, and the log must show the ranks and threads
#  asked:
#    steady     6 blocks: a 2 x 2 box of 10 x 10 x 3 cells around a vertex the
#               four blocks share (101 connections, walls 302 / 301 / 303 / 305,
#               symmetry in z) and a chimera pair of 16 x 16 x 3 and 15 x 15 x 3
#               cells overlapping by 3 x 3 cells, two chimera faces meeting at an
#               edge in both blocks; the mesh sheared (x += 0.3 y); steady, local
#               time step, cons, IRS (irs-beta = 0.5), 30 iterations;
#    transient  the same, time-accurate RK3, prim, 20 steps.
#  On 2 and 4 ranks the blocks of the chimera pair sit on different ranks, and so
#  does the vertex of the box (LPT split, Mod_MPI.f90).
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$(readlink -f "${1:?usage: run.sh <FUSS binary> <work dir>}")"
WORK="${2:?usage: run.sh <FUSS binary> <work dir>}"
MPIEXEC=${FUSS_MPIEXEC:-mpirun}
LAYOUTS="1x4 2x1 4x1 2x2"

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M MPLBACKEND=Agg

rm -rf "$WORK"; mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
cd "$WORK"

fail() { echo "MPIBitIdentity: FAIL: $*"; exit 1; }

python3 -B "$HERE/make_case.py" bit --outdir "$WORK/case" > /dev/null || fail "case generation"

# write_ini <case> <file>
write_ini() {
  local accurate=false iters=30 vars=cons irs=0.5
  if [[ $1 == transient ]]; then accurate=true; iters=20; vars=prim; irs=0.0; fi
  printf '%s\n' '[FUSS-Parameters]' 'newrun = true' 'res-threshold = 1e-300' "iter-threshold = $iters" '' \
    '[FUSS-Numerics]' 'time-scheme = RK3' 'vnn = 0.5' "time-accurate = $accurate" \
    "integration-variables = $vars" "irs-beta = $irs" '' \
    '[FUSS-IO]' 'sol-diter = 1000000' 'sol-overwrite = true' 'sol-format = vtk raw' 'ini-diter = 1000000' > "$2"
}

# run <case> <ranks> <threads>: solve and keep the folder as runs/<case>/<ranks>x<threads>
run() {
  local c="$1" np="$2" nt="$3" d="runs/$1/$2x$3"
  rm -rf "$d"; mkdir -p "$d/OUTPUT"
  cp -r case/INPUT case/MESH "$d/"
  write_ini "$c" "$d/input.ini"
  if ! (cd "$d" && OMP_NUM_THREADS="$nt" "$MPIEXEC" -np "$np" "$BIN" > solve.log 2>&1); then
    tail -20 "$d/solve.log"; fail "solver exited non-zero ($c ${np}x$nt)"
  fi
  [[ -s "$d/OUTPUT/residual-history.dat" ]] || fail "no residual history ($c ${np}x$nt)"
  ls "$d"/OUTPUT/vtk/fieldBlock*.vts > /dev/null 2>&1 || fail "no field output ($c ${np}x$nt)"
  grep -q "Number of threads --> *$nt\$" "$d/solve.log" || fail "the run did not use $nt threads ($c ${np}x$nt)"
  grep -q "Number of ranks *--> *$np\$" "$d/solve.log" || fail "the run did not use $np ranks ($c ${np}x$nt)"
}

# same <case> <layout>: every file of runs/<case>/1x1/OUTPUT identical to that of <layout>, and the same list of files
same() {
  local c="$1" b="$2" n=0 f
  diff <(cd "runs/$c/1x1/OUTPUT" && find . -type f | sort) <(cd "runs/$c/$b/OUTPUT" && find . -type f | sort) > /dev/null \
    || { echo "  different file lists: $c 1x1 vs $b"; return 1; }
  while IFS= read -r f; do
    n=$((n + 1))
    cmp -s "runs/$c/1x1/OUTPUT/$f" "runs/$c/$b/OUTPUT/$f" || { echo "  differs: $c ${f#./} (1x1 vs $b)"; return 1; }
  done < <(cd "runs/$c/1x1/OUTPUT" && find . -type f | sort)
  [[ $n -gt 0 ]] || { echo "  nothing compared: $c 1x1 vs $b"; return 1; }
  echo "  identical: $c 1x1 vs $b ($n files)"
}

status=0
for c in steady transient; do
  run "$c" 1 1
  for cfg in $LAYOUTS; do
    run "$c" "${cfg%x*}" "${cfg#*x}"
    same "$c" "$cfg" || status=1
  done
done

if [[ $status -ne 0 ]]; then
  echo "MPIBitIdentity: FAIL"
  exit 1
fi
echo "MPIBitIdentity: ok: steady (cons, IRS) and time-accurate (prim) on a box of 4 blocks and a chimera pair: 1x1, 1x4, 2x1, 4x1 and 2x2 bit-identical"
exit 0
