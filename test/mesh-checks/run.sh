#!/usr/bin/env bash
#===============================================================================
#  MetricChecks: the checks of the mesh at set-up, and the type of the mesh
#
#  run.sh <FUSS binary> <work dir> [part ...]: every part writes its deck with make_case.py into <work dir>/<part>, runs the
#  solver once (one thread) and checks the exit status, the messages and, where there is one, the field:
##    thin3d    two Cartesian 3-D blocks side by side, much thinner in z than wide, block 1 one cell thick in k, block 2 with
#              four cells, adiabatic lateral faces, walls 302 at the k faces with T linear in z: the steady field is linear in
#              z to 1e-9 K in both (when the type of the mesh came from block 1 it was a 2-D mesh: the step ignored the k
#              direction of block 2 and the run left the properties table)
#    mixed2d   a plane block and an axisymmetric wedge in one 2-D mesh: set-up stops
#    fold2d, fold3d  one node moved across the next cell: set-up stops, with block, cell and corner
#    chimturn  a chimera face whose edges along it grow 4 times across the first cell: the first ghost layer behind it turns
#              over, set-up stops (exit status not 0), with the ghost cell, the face and the ratio of the first two cells
#    chimwall  the same mesh with a wall (302) where the chimera face was: it runs
#  With no part, every part runs. The exit status is 0 when every part passes.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN=$(readlink -f "${1:?FUSS binary}")
WORK=${2:?work dir}
shift 2
PARTS=(thin3d mixed2d fold2d fold3d chimturn chimwall)
[[ $# -gt 0 ]] && PARTS=("$@")
export OMP_NUM_THREADS=1 KMP_STACKSIZE=100M MPLBACKEND=Agg
ulimit -s unlimited 2>/dev/null

status=0
nfail=0
fail() { echo "  FAIL: $*"; status=1; nfail=$((nfail + 1)); }

# expect_stop <part> <message regexp>...: the run stops (exit status not 0) and its log holds each message
expect_stop() {
  local part=$1 st=$2; shift 2
  [[ $st -ne 0 ]] || fail "$part: the run did not stop"
  for re in "$@"; do
    grep -qE "$re" "$WORK/$part/run.log" || fail "$part: no message matching '$re' in $WORK/$part/run.log"
  done
}

# field_error <part>: max |T - (A + GZ z)| over the cells of OUTPUT/field.tec (the thin3d field), K
field_error() {
  python3 -I - "$WORK/$1/OUTPUT/field.tec" "$HERE" <<'PY'
import sys
sys.path.insert(0, sys.argv[2])
import make_case as mc
zones, names, vals, i = {}, [], [], 0
text = open(sys.argv[1]).read().split('\n')
cur = None
for line in text:
    s = line.strip()
    if s.upper().startswith('VARIABLES'):
        import re
        names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', s.split('=', 1)[1])]
    elif s.upper().startswith('ZONE'):
        import re
        cur = re.search(r'T\s*=\s*"?([^",\s]+)', s).group(1)
        dims = tuple(int(re.search(r'\b%s\s*=\s*(\d+)' % c, s).group(1)) for c in 'IJK')
        zones[cur] = [dims, []]
    elif s and not s[0].isalpha() and cur:
        zones[cur][1] += [float(x) for x in s.split()]
worst = 0.0
for title, (dims, v) in zones.items():
    ni, nj, nk = dims
    nn = ni * nj * nk
    nc = (ni - 1) * (nj - 1) * (nk - 1)
    iT = names.index('T')
    off = 3 * nn + (iT - 3) * nc
    zs = v[2 * nn:3 * nn]
    for k in range(nk - 1):
        for j in range(nj - 1):
            for i in range(ni - 1):
                zc = sum(zs[(k + c) * ni * nj + (j + b) * ni + i + a] for a in (0, 1) for b in (0, 1) for c in (0, 1)) / 8
                t = v[off + (k * (nj - 1) + j) * (ni - 1) + i]
                worst = max(worst, abs(t - (mc.A + mc.GZ * zc)))
print('%.3e' % worst)
PY
}

for part in "${PARTS[@]}"; do
  echo "=== $part"
  before=$nfail
  rm -rf "$WORK/$part"; mkdir -p "$WORK/$part"
  python3 -I "$HERE/make_case.py" "$part" "$WORK/$part" > /dev/null || { fail "$part: no deck"; continue; }
  (cd "$WORK/$part" && "$BIN" > run.log 2>&1); st=$?
  case $part in
    thin3d)
      [[ $st -eq 0 ]] || fail "thin3d: the run stopped (exit $st), see $WORK/thin3d/run.log"
      err=$(field_error thin3d); echo "  max |T - exact| = $err K (a block one cell thick in k next to a 3-D block, T linear in z)"
      awk -v e="$err" 'BEGIN { exit !(e + 0 <= 1e-9) }' || fail "thin3d: the field is not the linear one (block 1 decoupled?)" ;;
    mixed2d)
      expect_stop mixed2d $st 'Mesh type: 1 of the 2 blocks of the 2-D mesh are axisymmetric wedges' ;;
    fold2d|fold3d)
      expect_stop $part $st 'Mesh.* block 1, cell .*, corner node .*: the cell is folded or turned against its block' ;;
    chimturn)
      expect_stop chimturn $st 'Metric tensor: block 1, cell .* ghost cell behind face 2' 'r = dl_2/dl_1 of the first two cells across the face' ;;
    chimwall)
      [[ $st -eq 0 ]] || fail "chimwall: the run stopped (exit $st), see $WORK/chimwall/run.log"
      [[ -s $WORK/chimwall/OUTPUT/residual-history.dat ]] || fail "chimwall: no residual history" ;;
    *) fail "unknown part $part" ;;
  esac
  [[ $nfail -eq $before ]] && echo "=== $part: PASS" || echo "=== $part: FAIL"
done
[[ $status -eq 0 ]] && echo "MetricChecks: all parts passed" || echo "MetricChecks: FAILED"
exit $status
