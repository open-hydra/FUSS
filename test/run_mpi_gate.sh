#!/bin/bash
#===============================================================================
#         FILE: test/run_mpi_gate.sh
#        USAGE: ./test/run_mpi_gate.sh [-w WORKDIR] [-r] [-c CASE:NRANKS[:verify]]... [-b FUSS_BINARY]
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
#                 - residual-history.dat is compared byte for byte too (plan 14 S5: the
#                   residual is summed per block, ours' pyrolysis/Darcy terms included,
#                   and over the blocks in block order, e4b5cf1); a difference names the
#                   first differing iteration, the largest relative difference and the
#                   row counts. -r is the default; FUSS_MPIGATE_RESIDUAL=0 skips it.
#               Before the cases: test/check_mpi_tags.py (every point-to-point tag a
#               named constant, distinct, 0..7000; plan 14 D4).
#               So a case must stop on a fixed iteration count or time, never on
#               res-threshold (conjugate_two_layer_fixed exists for that).
#
#               A moving mesh (plan 11 stage 2): in BOTH runs every wall.tec patch
#               must lie on the face nodes of the field.tec written at the same
#               output (string equality: one writer). Before stage 2 the patches
#               kept the setup surface, on one rank too.
#
#               CASE:N:verify (plan 11 stage 4) also runs the case's own verify.py
#               on the N-rank output, which must pass (the receding conjugate
#               slabs: the interface fluxes of the two sides pair to 1e-12).
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK="${TMPDIR:-/tmp}/fuss-mpi-gate.$$"
CASES=()
BIN="${FUSS_GATE_BIN:-}"

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -c | --case) CASES+=("$2"); shift 2 ;;
    -r | --residual) export FUSS_MPIGATE_RESIDUAL=1; shift ;;
    -b | --binary) BIN=$2; shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR] [-r] [-c CASE:NRANKS[:verify]]... [-b FUSS_BINARY]"; exit 1 ;;
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
# (a correction group whose blocks share their line positions); stage 4: the ablating wall -- the receding
# conjugate slabs on two ranks with their interface pairing (verify.py), ATC 2.2 cut in four with the
# receding block on the last rank (the rates' allgather, the faces' state to root), the smoothed slab whose
# two blocks share the receding surface (the nodal rates across ranks, the interface assertion), the coupled
# ablating wall 307 cut in four with its faces on the last rank (hydra's contract: the exchange values exist
# on the rank of the face's block only, as the TEST ONLY coupled-env-file now writes them); plan 14 S9d Part B:
# the periodic sectors (test/moving-mesh/periodic) on two ranks, blocks A and B on different ranks: the periodic pair
# and the block connection between them across the ranks, the per-step check of the pair, the joined seam of a
# receding wall (its line directions from the rank of the representative), the Darcy rows across the pair
[ "${#CASES[@]}" -gt 0 ] || CASES=("test/transient/multimat_Plate:4" "test/numerics-features/conjugate_two_layer_fixed:2"
  "test/moving-mesh/freestream_4blk:4" "test/moving-mesh/conservation_4blk:4" "test/moving-mesh/receding_4blk:4"
  "test/moving-mesh/receding/N040:2" "test/numerics-features/chimera_fixed:2"
  "test/numerics-features/residual_convergence/mg1:2"
  "test/pyrolysis/energy_audit:4" "test/pyrolysis/heated0d:2"
  "test/pyrolysis/lit/baum_atreya2007/full_N100_20blk:4" "test/pyrolysis/darcy/tensor_linear_2d_2x2:4"
  "test/pyrolysis/darcy/mass_balance_1d_2x:2"
  "test/moving-mesh/conjugate_receding:2:verify" "test/moving-mesh/conjugate_receding_t2:2:verify"
  "test/pyrolysis/ablation/atc22_4blk:4" "test/pyrolysis/ablation/smoothing/w5_2blk:2"
  "test/pyrolysis/ablation/coupled307/run307_4blk:4"
  "test/moving-mesh/periodic/translate:2" "test/moving-mesh/periodic/prescribed:2"
  "test/moving-mesh/periodic/receding:2" "test/moving-mesh/periodic/receding1:2"
  "test/pyrolysis/darcy/periodic_static:2" "test/pyrolysis/darcy/periodic_static1:2"
  "test/pyrolysis/ablation/periodic_receding:2")

export FUSS_MPIGATE_RESIDUAL=${FUSS_MPIGATE_RESIDUAL:-1}   # residual-history.dat compared by default (plan 14 S5)
[ -n "${I_MPI_ROOT:-}" ] || { echo "needs I_MPI_ROOT (the Intel MPI of the mpi preset)"; exit 1; }
export FUSS_EXE=FUSS-mpi FUSS_MPIRUN="$I_MPI_ROOT/bin/mpirun"
"$ROOT/test/check_build_config.sh" --preset mpi || { echo "refusing to grade a build that does not match the preset"; exit 1; }

mkdir -p "$WORK"
# -b PATH (or FUSS_GATE_BIN): every run takes that executable, copied once into $WORK/.master/bin/FUSS-mpi, which is
# then the MASTERDIR of the cases' FUSS.sh; without it, bin/FUSS-mpi of this tree as before
MASTERDIR=$ROOT
if [ -n "$BIN" ]; then
  if [ ! -f "$BIN" ] || [ ! -x "$BIN" ]; then echo "no FUSS executable at $BIN"; exit 1; fi
  MASTERDIR="$WORK/.master"
  mkdir -p "$MASTERDIR/bin" && cp "$BIN" "$MASTERDIR/bin/$FUSS_EXE" || exit 1
fi
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS MPI gate: N ranks against one"
echo "  binary   : $MASTERDIR/bin/FUSS-mpi ($(md5sum "$MASTERDIR/bin/FUSS-mpi" | cut -c1-8))"
echo "  work dir : $WORK"
echo
if python3 -I "$ROOT/test/check_mpi_tags.py" --dir "$ROOT/src" > "$WORK/mpi-tags.log" 2>&1; then
  pass "MPI tags: $(tail -1 "$WORK/mpi-tags.log")"
else
  bad "MPI tags: a point-to-point tag is not a distinct named constant (see $WORK/mpi-tags.log)"
fi

for spec in "${CASES[@]}"; do
  IFS=: read -r case_rel n check <<< "$spec"
  src="$ROOT/$case_rel"; name=$(echo "$case_rel" | tr '/' '_')
  for r in 1 "$n"; do
    d="$WORK/$name.m$r"
    rm -rf "$d"; mkdir -p "$d"
    cp -r "$src/INPUT" "$src/input.ini" "$src/FUSS.sh" "$d/"
    [ -d "$src/MESH" ] && cp -r "$src/MESH" "$d/"
    cp "$src"/*.dat "$d/" 2>/dev/null || true
    sed -i "s|^MASTERDIR=.*|MASTERDIR=$MASTERDIR|" "$d/FUSS.sh"
    mkdir -p "$d/OUTPUT" "$d/bin"
    ( cd "$d" && timeout 900 ./FUSS.sh -p 1 -m "$r" solve ) > "$d/run.log" 2>&1
    rc=$?
    # every stop exits 1; a NaN or out-of-table temperature (Check_And_Fix_State) also prints two [ERROR] lines
    # and "[RANK n] STOP: NaN or T outside the properties table"
    if [ $rc -ne 0 ] || grep -q "\[ERROR\]\|NaN or T outside the properties table" "$d/run.log"; then
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
        nz = lambda s: re.sub(rb'-(0\.0+E\+000)', rb' \1', s)   # the sign's column (15 digits, or 17 from ORION a978930 on)
        if nz(x) == nz(y):
            signed_zero += 1; continue
    bad.append(f)
resid = ''
if os.environ.get('FUSS_MPIGATE_RESIDUAL') == '1':
    ha, hb = (open(os.path.join(d, 'residual-history.dat'), 'rb').read() for d in (a, b))
    if ha == hb:
        resid = '; residual-history.dat byte-identical'
    else:
        ra = [l.split() for l in ha.decode().split('\n') if l.strip() and not l.lstrip().startswith('#')]
        rb = [l.split() for l in hb.decode().split('\n') if l.strip() and not l.lstrip().startswith('#')]
        first = next((x[0] for x, y in zip(ra, rb) if x != y), None)
        rel = 0.0
        for x, y in zip(ra, rb):
            for p, q in zip(x[1:], y[1:]):
                p, q = float(p), float(q)
                if p != q: rel = max(rel, abs(p - q) / max(abs(p), abs(q)))
        bad.append('residual-history.dat (first difference at iteration %s, max rel %.1e, %d vs %d rows)'
                   % (first, rel, len(ra), len(rb)))
print('%d files: %d byte-identical, %d wall files equal up to the sign of zero%s%s%s'
      % (len(fa), same, signed_zero, '; wall patches on the mesh' if moving else '', resid,
         '' if not bad else '; DIFFERENT: ' + ', '.join(bad[:6])))
sys.exit(1 if bad else 0)
PY
)
  if [ $? -eq 0 ]; then
    pass "$case_rel: $n ranks = 1 rank ($out)"
  else
    bad "$case_rel: $n ranks differ from 1 rank ($out)"
  fi
  if [ "${check:-}" = "verify" ]; then
    if python3 "$src/verify.py" "$WORK/$name.m$n" > "$WORK/$name.m$n/verify.log" 2>&1; then
      pass "$case_rel: verify.py on the $n-rank output ($(grep -m1 -E '\(pair\)' "$WORK/$name.m$n/verify.log" | tr -s ' ' | cut -c1-110))"
    else
      bad "$case_rel: verify.py FAILED on the $n-rank output (see $WORK/$name.m$n/verify.log)"
    fi
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "MPI gate passed"
else
  echo "$fail case(s) FAILED"
  exit 1
fi
