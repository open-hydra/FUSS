#!/usr/bin/env bash
#===============================================================================
#  MPI calls only from the main thread (needs-mpi tier)
#
#  usage: run.sh <FUSS binary> <work dir>      (FUSS_MPI_THREAD_PROBE = the probe library, set by CTest)
#
#  FUSS initialises MPI with MPI_THREAD_FUNNELED (Mod_MPI): the MPI standard then allows MPI calls only from the thread that
#  called MPI_INIT_THREAD. The halo exchange of the temperature, the IRS residual exchange and the chimera exchange run
#  inside the OpenMP parallel regions of the Runge-Kutta stages, so every MPI call there must be made by that thread. A
#  probe library preloaded into the solver (mpi_thread_probe.c, built with the tests) counts the MPI calls made by any other
#  thread; this test asserts that the count is zero on every rank for
#      2x2, 2x4 and 2x2 with implicit residual smoothing
#  (MPI ranks x OpenMP threads), that the probe did intercept the halo exchange on every rank (MPI_STARTALL, MPI_ISEND and
#  MPI_IRECV of the chimera exchange: otherwise the zero would mean nothing), and that the largest tag of the point-to-point
#  calls is at most 32767, the least MPI_TAG_UB the MPI standard allows (MPI 4.1, Sec. 3.2.3). The case is that of
#  MPIBitIdentity (make_case.py bit: a box of 4 blocks joined by connections and a chimera pair, 6 blocks), steady, cons, 10
#  iterations; the run with residual smoothing has irs-beta = 0.5.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$(readlink -f "${1:?usage: run.sh <FUSS binary> <work dir>}")"
WORK="${2:?usage: run.sh <FUSS binary> <work dir>}"
MPIEXEC=${FUSS_MPIEXEC:-mpirun}
PROBE=${FUSS_MPI_THREAD_PROBE:-}

fail() { echo "MPIThreadFunnel: FAIL: $*"; exit 1; }
[[ -n "$PROBE" && -f "$PROBE" ]] || fail "probe library not found (FUSS_MPI_THREAD_PROBE='$PROBE')"
PROBE="$(readlink -f "$PROBE")"

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M MPLBACKEND=Agg

rm -rf "$WORK"; mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
cd "$WORK"
python3 -B "$HERE/../bit-identity/make_case.py" bit --outdir "$WORK/case" > /dev/null || fail "case generation"

# write_ini <irs-beta> <file>
write_ini() {
  printf '%s\n' '[FUSS-Parameters]' 'newrun = true' 'res-threshold = 1e-300' 'iter-threshold = 10' '' \
    '[FUSS-Numerics]' 'time-scheme = RK3' 'vnn = 0.5' 'time-accurate = false' 'integration-variables = cons' \
    "irs-beta = $1" '' \
    '[FUSS-IO]' 'sol-diter = 1000000' 'sol-overwrite = true' 'sol-format = vtk raw' 'ini-diter = 1000000' > "$2"
}

status=0
for run in 2x2 2x4 2x2-irs; do
  cfg=${run%%-*}
  np=${cfg%x*}; nt=${cfg#*x}
  d="runs/$run"
  rm -rf "$d"; mkdir -p "$d/OUTPUT"
  cp -r case/INPUT case/MESH "$d/"
  if [[ $run == *-irs ]]; then write_ini 0.5 "$d/input.ini"; else write_ini 0.0 "$d/input.ini"; fi
  if ! (cd "$d" && OMP_NUM_THREADS="$nt" "$MPIEXEC" -np "$np" env LD_PRELOAD="$PROBE" "$BIN" > solve.log 2> probe.log); then
    tail -n 20 "$d/solve.log" "$d/probe.log"; fail "solver exited non-zero ($run)"
  fi
  # with a single thread every call would come from the main thread: the check means something only if the launcher
  # passed OMP_NUM_THREADS to the ranks
  grep -qE "Number of threads *--> *$nt\$" "$d/solve.log" || fail "solver did not run with $nt threads ($run)"
  grep -qE "Number of ranks *--> *$np\$" "$d/solve.log" || fail "solver did not run with $np ranks ($run)"
  if [[ $run == *-irs ]]; then
    grep -qE 'Implicit residual smoothing +[^ ]' "$d/solve.log" && ! grep -qE 'Implicit residual smoothing +none' "$d/solve.log" \
      || fail "implicit residual smoothing not active ($run)"
  else
    grep -qE 'Implicit residual smoothing +none' "$d/solve.log" || fail "implicit residual smoothing active in a run without it ($run)"
  fi
  lines=$(grep -c '^MPI-THREAD-PROBE ' "$d/probe.log")
  [[ "$lines" -eq "$np" ]] || fail "expected $np probe reports, found $lines ($run)"
  # the exchanges must have been seen on every rank: MPI_STARTALL (connections) and MPI_ISEND (chimera) at every stage
  seen=$(grep '^MPI-THREAD-PROBE ' "$d/probe.log" | grep ' mpi_startall_:[1-9]' | grep -c ' mpi_isend_:[1-9]')
  [[ "$seen" -eq "$np" ]] || fail "probe did not intercept the halo exchange on every rank ($run)"
  bad=$(grep '^MPI-THREAD-PROBE ' "$d/probe.log" | grep -vc ' other=0 ')
  if [[ "$bad" -ne 0 ]]; then
    echo "  $run: MPI calls from threads other than the main thread:"
    grep '^MPI-THREAD-PROBE ' "$d/probe.log" | grep -v ' other=0 ' | sed 's/^/    /'
    status=1
  else
    echo "  $run: all MPI calls from the main thread on every rank ($(grep '^MPI-THREAD-PROBE ' "$d/probe.log" | grep -o ' calls=[0-9]*' | tr -d '\n'))"
  fi
  maxtag=$(grep '^MPI-THREAD-PROBE ' "$d/probe.log" | grep -o ' maxtag=-\?[0-9]*' | cut -d= -f2 | sort -n | tail -1)
  [[ -n "$maxtag" ]] || fail "probe did not report the largest tag ($run)"
  if [[ "$maxtag" -gt 32767 ]]; then
    echo "  $run: MPI tag $maxtag above 32767, the least MPI_TAG_UB the MPI standard allows"
    status=1
  else
    echo "  $run: largest MPI tag $maxtag"
  fi
done

if [[ $status -ne 0 ]]; then
  echo "MPIThreadFunnel: FAIL"
  exit 1
fi
echo "MPIThreadFunnel: ok: MPI called only by the main thread, tags at most 32767 (2x2, 2x4, 2x2 with IRS)"
exit 0
