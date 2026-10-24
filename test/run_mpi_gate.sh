#!/bin/bash
#===============================================================================
#         FILE: test/run_mpi_gate.sh
#        USAGE: ./test/run_mpi_gate.sh [-w WORKDIR] [-c CASE:NRANKS]...
#  DESCRIPTION: A case run on N MPI ranks must give the outputs of the same case
#               run on one rank (plan 11 section 2.5). Both runs use the mpi
#               preset's bin/FUSS-mpi at OMP_NUM_THREADS=1, the N-rank one through
#               $I_MPI_ROOT/bin/mpirun under timeout 900 (a deadlock is a FAIL,
#               not a hang).
#
#               Rules: every OUTPUT file byte-identical, except
#                 - wall*.tec, compared with -0.0 equal to +0.0: the wall values
#                   are summed over the ranks, each entry from the rank that owns
#                   it, and a -0.0 of the owner becomes +0.0 once another rank's
#                   zero is added (the decision of 2026-10-06 keeps the sum);
#                 - residual-history.dat, not compared: the residual is an MPI sum
#                   whose rounding depends on the rank count.
#               So a case must stop on a fixed iteration count or time, never on
#               res-threshold (conjugate_two_layer_fixed exists for that).
#
#               A moving mesh (plan 11 stage 2): in BOTH runs every wall.tec patch
#               must lie on the face nodes of the field.tec written at the same
#               output (string equality: one writer). Before stage 2 the patches
#               kept the setup surface, on one rank too.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK="${TMPDIR:-/tmp}/fuss-mpi-gate.$$"
CASES=()

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -c | --case) CASES+=("$2"); shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR] [-c CASE:NRANKS]..."; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
# plan 11 stage 1, static: three blocks on four ranks (one idle), two blocks with IRS (the R exchange) on
# two; stage 2, a moving mesh: the 2x2 and 4x1 splits on four ranks (every block connection between two
# ranks), one moving block on two ranks (one idle); upstream's chimera-mpi (merged 2026-10-17): two blocks
# joined by chimera faces, each donor on the other rank, and two multigrid levels, on two ranks; stage 3a:
# type-1 pyrolysis (single-block, guard 9) on four and on two ranks, every other rank idle through the
# collectives -- the energy budget, the probes and the output files must not move; stage 3b: type 2
# split across the ranks -- twenty blocks in a row with the kinetics (the chi, p_pg, Gamma and h_g Gamma
# halos, the correction group of blocks in a row), a 2 x 2 rotated-tensor slab under the Jacobi sweep
# (the pressure stencil pg(3:6), the sweep control over thousands of sweeps), two blocks side by side
# (a correction group whose blocks share their line positions)
[ "${#CASES[@]}" -gt 0 ] || CASES=("test/transient/multimat_Plate:4" "test/numerics-features/conjugate_two_layer_fixed:2"
  "test/moving-mesh/freestream_4blk:4" "test/moving-mesh/conservation_4blk:4" "test/moving-mesh/receding_4blk:4"
  "test/moving-mesh/receding/N040:2" "test/numerics-features/chimera_fixed:2"
  "test/numerics-features/residual_convergence/mg1:2"
  "test/pyrolysis/energy_audit:4" "test/pyrolysis/heated0d:2"
  "test/pyrolysis/lit/baum_atreya2007/full_N100_20blk:4" "test/pyrolysis/darcy/tensor_linear_2d_2x2:4"
  "test/pyrolysis/darcy/mass_balance_1d_2x:2")

[ -n "${I_MPI_ROOT:-}" ] || { echo "needs I_MPI_ROOT (the Intel MPI of the mpi preset)"; exit 1; }
export FUSS_EXE=FUSS-mpi FUSS_MPIRUN="$I_MPI_ROOT/bin/mpirun"
"$ROOT/test/check_build_config.sh" --preset mpi || { echo "refusing to grade a build that does not match the preset"; exit 1; }

mkdir -p "$WORK"
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS MPI gate: N ranks against one"
echo "  binary   : $ROOT/bin/FUSS-mpi ($(md5sum "$ROOT/bin/FUSS-mpi" | cut -c1-8))"
echo "  work dir : $WORK"
echo

for spec in "${CASES[@]}"; do
  case_rel=${spec%%:*}; n=${spec##*:}
  src="$ROOT/$case_rel"; name=$(echo "$case_rel" | tr '/' '_')
  for r in 1 "$n"; do
    d="$WORK/$name.m$r"
    rm -rf "$d"; mkdir -p "$d"
    cp -r "$src/INPUT" "$src/input.ini" "$src/FUSS.sh" "$d/"
    [ -d "$src/MESH" ] && cp -r "$src/MESH" "$d/"
    cp "$src"/*.dat "$d/" 2>/dev/null || true
    sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
    mkdir -p "$d/OUTPUT" "$d/bin"
    ( cd "$d" && timeout 900 ./FUSS.sh -p 1 -m "$r" solve ) > "$d/run.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ] || grep -q "\[ERROR\]\|NaN or T<0" "$d/run.log"; then
      bad "$case_rel: the $r-rank run failed (exit $rc$([ $rc -eq 124 ] && echo ', timed out'); see $d/run.log)"
      continue 2
    fi
  done
  out=$(python3 - "$WORK/$name.m1/OUTPUT" "$WORK/$name.m$n/OUTPUT" "$src/input.ini" <<'PY'
import os, re, sys
a, b, ini = sys.argv[1], sys.argv[2], sys.argv[3]

def zones(path):
    # name -> (I, J, K, the x, y, z strings); every zone, nodal coordinates first
    txt = open(path).read(); out = {}
    for part in re.split(r'\n\s*ZONE', txt)[1:]:
        head, body = part.split('\n', 1)
        name = re.search(r'T\s*=\s*"?([^",]+)', head).group(1).strip()
        I, J, K = (int(x) for x in re.search(r'I=(\d+),\s*J=(\d+),\s*K=(\d+)', head).groups())
        v = body.split(); n = I * J * K
        out[name] = (I, J, K, v[0:n], v[n:2*n], v[2*n:3*n])
    return out

def wall_on_mesh(d):
    # every wall patch B<b>F<f> against face f of block b of the field written with it
    bad = 0; checked = 0
    for fw in sorted(os.listdir(d)):
        m = re.match(r'wall(\d*)\.tec$', fw)
        if not m: continue
        ff = os.path.join(d, 'field%s.tec' % m.group(1))
        if not os.path.exists(ff): continue
        F = list(zones(ff).values())
        for name, (I, J, K, x, y, z) in zones(os.path.join(d, fw)).items():
            bb, f = (int(t) for t in re.match(r'B(\d+)F(\d+)', name).groups())
            FI, FJ, FK, fx, fy, fz = F[bb - 1]
            ax = (f - 1) // 2; hi = f % 2 == 0
            for k in range(K):
                for j in range(J):
                    for i in range(I):
                        ijk = [i, j, k]; ijk[ax] = ([FI, FJ, FK][ax] - 1) if hi else 0
                        p = ijk[0] + FI * (ijk[1] + FJ * ijk[2]); q = i + I * (j + J * k)
                        checked += 1
                        if (x[q], y[q], z[q]) != (fx[p], fy[p], fz[p]): bad += 1
    return bad, checked

moving = re.search(r'^\s*law\s*=\s*(\w+)', open(ini).read(), re.M)
moving = moving is not None and moving.group(1) != 'static'
if moving:
    for d in (a, b):
        nb, nc = wall_on_mesh(d)
        if nb:
            print('wall.tec patches off the mesh: %d of %d nodes (%s)' % (nb, nc, os.path.basename(os.path.dirname(d)))); sys.exit(1)
fa = sorted(f for f in os.listdir(a) if f != 'residual-history.dat')
fb = sorted(f for f in os.listdir(b) if f != 'residual-history.dat')
if fa != fb:
    print('different file sets: %s' % sorted(set(fa) ^ set(fb))); sys.exit(1)
same = signed_zero = 0; bad = []
for f in fa:
    x = open(os.path.join(a, f), 'rb').read(); y = open(os.path.join(b, f), 'rb').read()
    if x == y:
        same += 1; continue
    if f.startswith('wall'):
        nz = lambda s: s.replace(b'-0.000000000000000E+000', b' 0.000000000000000E+000')   # the sign's column
        if nz(x) == nz(y):
            signed_zero += 1; continue
    bad.append(f)
print('%d files: %d byte-identical, %d wall files equal up to the sign of zero%s%s'
      % (len(fa), same, signed_zero, '; wall patches on the mesh' if moving else '',
         '' if not bad else '; DIFFERENT: ' + ', '.join(bad[:6])))
sys.exit(1 if bad else 0)
PY
)
  if [ $? -eq 0 ]; then
    pass "$case_rel: $n ranks = 1 rank ($out)"
  else
    bad "$case_rel: $n ranks differ from 1 rank ($out)"
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "MPI gate passed"
else
  echo "$fail case(s) FAILED"
  exit 1
fi
