#!/usr/bin/env bash
# ChimeraEdge: a linear field on a sheared chimera pair, with walls consistent with it, is preserved to round-off, also at the
# cells at the ends of the chimera faces (make_case.py edge, chimera_edge_check.py).   usage: chimera_edge.sh <FUSS binary> <work dir>
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$(readlink -f "${1:?usage: chimera_edge.sh <FUSS binary> <work dir>}")"
WORK="${2:?usage: chimera_edge.sh <FUSS binary> <work dir>}"
ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg

rm -rf "$WORK"; mkdir -p "$WORK/OUTPUT"; cd "$WORK"
python3 -B "$HERE/make_case.py" edge --outdir "$WORK" > /dev/null || { echo "ChimeraEdge: FAIL: case generation"; exit 1; }
printf '%s\n' '[FUSS-Parameters]' 'newrun = true' 'res-threshold = 1e-300' 'iter-threshold = 5' '' \
  '[FUSS-Numerics]' 'time-scheme = RK3' 'vnn = 0.5' 'time-accurate = false' 'integration-variables = prim' '' \
  '[FUSS-IO]' 'sol-diter = 1000000' 'sol-overwrite = true' 'sol-format = tecplot ascii' 'ini-diter = 1000000' > input.ini
"$BIN" > run.log 2>&1 || { tail -5 run.log; echo "ChimeraEdge: FAIL: solver exited non-zero"; exit 1; }
python3 -B "$HERE/chimera_edge_check.py" "$WORK" 1e-9 || { echo "ChimeraEdge: FAIL"; exit 1; }
echo "ChimeraEdge: ok"
