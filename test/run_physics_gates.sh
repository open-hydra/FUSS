#!/bin/bash
#===============================================================================
#         FILE: test/run_physics_gates.sh
#        USAGE: ./test/run_physics_gates.sh [-p NTHREADS] [-w WORKDIR] [-c CASE] [-a]
#  DESCRIPTION: Run the solver on a fresh copy of each physics-claim case and
#               grade it with that case's own verify.py.
#
#               WHY THIS EXISTS
#               ---------------
#               Until 2026-09-15 no script invoked any of these verify.py
#               scripts. The only automated correctness check in this tree was
#               md5 identity to test/baseline/manifest.txt -- and a manifest
#               is re-recordable on request (run_regression.sh -B). A change
#               that breaks conservation, or halves bc%Mg, or removes mg3's
#               coarse work, and is then baselined, passes every other gate in
#               test/ without complaint. This script is the one that actually
#               checks the physics, not just that the bytes did not move.
#
#               Every case is copied into $WORK first (see `stage`, below) so
#               that the in-tree OUTPUT/ is never disturbed and the committed
#               SOLUTION/ reference data is never at risk -- run_regression.sh
#               already lost ten tracked files once by wiping SOLUTION/, and
#               that mistake is not repeated here.
#
#               CASES, AND WHAT EACH ONE PROVES
#               --------------------------------
#                 freestream            -- a uniform field with
#                                        arbitrary prescribed mesh motion must
#                                        stay uniform to ~1e-9 K (measured:
#                                        ~1e-12 K); anything else is a GCL
#                                        violation, geometry vs. flux disagree.
#                 conservation          -- total volumetric-
#                                        enthalpy energy on a moving mesh with
#                                        adiabatic walls is conserved to 1e-10
#                                        relative; the reference cell volumes
#                                        are recomputed independently from the
#                                        node coordinates, not read back from
#                                        the solver.
#                 conjugate_two_layer   -- a two-material slab with
#                                        UNEQUAL cell spacing across the block
#                                        interface (dx2/dx1=4) must reproduce
#                                        the exact (piecewise-linear) interface
#                                        temperature to 1e-6 K. This is the
#                                        only case in the suite where the
#                                        answer depends on bc%Mg being the
#                                        NEIGHBOUR's metric and not the local
#                                        one -- halving bc%Mg moves the answer
#                                        by 28 K here and changes nothing on
#                                        any equal-spacing case.
#                 receding              -- a receding surface
#                                        modelled as rigid translation must
#                                        show first-order decay-length error
#                                        against the closed form
#                                        d_num/d = Pe/ln(1+Pe), the exact root
#                                        of the steady 3-term recurrence for
#                                        donor-cell upwinding on the swept
#                                        volume. Default grids N040/N080 (see
#                                        NOTE ON -a below); -a adds the .slow
#                                        N160/N320.
#                 residual_convergence  -- at a FIXED 30000-
#                                        iteration fine-grid budget, more
#                                        coarse-grid pre-work must leave a
#                                        lower residual (nominal>mg1>mg2>mg3),
#                                        IRS must accelerate, and both
#                                        together must win. Six sub-cases,
#                                        ~205k iterations total -- this is the
#                                        expensive one (see RUNTIME below).
#
#               Reserved for future groups -- append here as they land, do NOT
#               start a second gate script: one physics gate for every group,
#               not one per phase:
#                 bc303_qrad            -- the radiative flux read for
#                                        a type-303 face must reach both the flux
#                                        and the wall output: verify.py checks the
#                                        BC's own identity qw - hconv*(Tref-Tw) =
#                                        qrad on all 160 face cells to 1e-9 rel.
#                                        Red on 6cad426 (residual = qrad), green
#                                        after the two-site fix. 200 iterations.
#                 pyrolysis   -- plan 06 phase 2, six sub-checks: kinetics_isothermal
#                                T0700/T0900/T1200 (closed-form chi, 1e-8; T exact),
#                                energy_audit (energy and solid-mass balances 1e-10,
#                                independent state sum, vnn 0.5/0.25/0.125 first-order
#                                refinement, the no-kinetics variant exact zeros, 1 vs 8
#                                threads byte-identical), ATC1 vs FIAT, and
#                                ATC2-chemistryOff vs the bundled PATO zero-emissivity
#                                run (or $FUSS_PATO_ATC2_OUT). ~3 min at -p 1.
#                 darcy_*     -- plan 07 phase 2b (gas-closure = type2), cases under
#                                test/pyrolysis/darcy/: sealed_nosource (p stays p-init
#                                exactly, zero flux), mass_balance_1d (steady face fluxes
#                                pi x: b2 time-accurate 1e-6, b1 quasi-steady sedded copy
#                                1e-9), sealed_profile (closed-form p^2 profile, observed
#                                order >= 1.9 on N020..N160), klinkenberg (the same grids with
#                                K + (mu/p) B: exact (p+b)^2 profile, b = mu B/K = p_w/2, order
#                                >= 1.9 and a shift >= 100 x the error), tensor_linear_2d (rotated K:
#                                |u_y|/|u_x| = K_yx/K_xx to 1e-10; the diagonal-K red control
#                                written into a staged copy), 0d (adiabatic0d/heated0d under
#                                type 2 vs the exact type-2 ODE), atc2_chemistry_off (type 2
#                                vs the bundled PATO run, T and p, plus the type-1 case on
#                                the same mesh for the sign check). ~2 min at -p 1.
#                 bench_*     -- 2026-10-13, benchmarks of the type-2 closure, ladders
#                                generated into $WORK by test/pyrolysis/tools/make_*.py:
#                                exact (test/pyrolysis/exact/, verify_exact.py):
#                                pressure_transient (series, order >= 1.8 on N020..N080),
#                                transpiration (exponential T + Darcy mdot, order >= 1.8),
#                                sharp_front (similarity, front order in [0.7, 1.4]);
#                                literature (test/pyrolysis/lit/): baum_atreya (Baum &
#                                Atreya 2007 corrected: front, T, char pressure, gas ledgers
#                                of both type-2 modes), zhang2010 (three-code slab, BC 305,
#                                within the codes' spread). ~4 min at -p 1.
#                 ablation    -- add gate_ablation()   + list it
#
#               NOTE ON -a (receding). verify.py's order/decay/pointwise-error
#               checks all pass on N040/N080 alone (order 0.881 on that pair,
#               inside the [0.85,1.15] band -- measured 2026-09-15; N040 is
#               2s, N080 15s at -p 1), so that pair is the default. N160/N320
#               are the ones marked .slow in the source tree (N320 is ~20
#               minutes per its own header) and are added by -a. GRIDS is
#               overridden by sedding the WORK COPY of generate.py, not by
#               editing verify.py or by any importlib trick -- same discipline
#               run_guard_demos.sh already uses on FUSS.sh's MASTERDIR, and it
#               keeps the invocation the same `python3 verify.py` a human uses.
#
#               HOW A CASE'S FUSS.sh FINDS THE BINARY. FUSS.sh resolves
#               MASTERDIR=../../.. (or ../../../.. one level deeper) relative
#               to its own location, so a copy under $WORK cannot find
#               $ROOT/bin/FUSS by relative path any more. `stage` (below) seds
#               that one line to an absolute path, exactly as
#               run_guard_demos.sh's demo() already does; nothing about
#               FUSS.sh's own logic changes. `stage` never copies a case's
#               bin/ directory: FUSS.sh only refreshes bin/FUSS "if $MASTER is
#               newer than $LOCAL" (bash `-nt`), which is also true whenever
#               $LOCAL is simply absent -- but a STALE local FUSS.sh binary
#               with a newer mtime than $ROOT/bin/FUSS would silently win and
#               the gate would grade the wrong executable. Always start bin/
#               empty.
#
#               THE FUSS_GATE_SELFTEST HOOK (kept, not removed -- documented
#               here per the acceptance note that a gate never seen to fail is
#               not a gate). FUSS_GATE_SELFTEST=1 does two things, neither of
#               which touches a checked-in verify.py or FUSS.sh:
#                 - freestream: a COPY of verify.py is written into the work
#                   dir with TOL_ABS sedded to 0.0, and that copy is run
#                   instead of the original -- proves the harness turns a
#                   non-zero verify.py exit into a FAIL line and a non-zero
#                   script exit (the "plumbing" proof, cheap: ~2s).
#                 - residual_convergence: the WORK COPY of mg3/input.ini has
#                   level2-iter sedded from 7500 to 1, i.e. mg3's coarse-grid
#                   work is taken away -- a real physics-assertion red, not
#                   just plumbing, input-only, no solver code touched. The
#                   original record had this firing the budget check AND the ordering
#                   check together (30001 vs expected 37500, mg2 no longer
#                   exceeding mg3). Re-measured here on 2026-09-15 against
#                   THIS binary (see the provenance note below): the budget
#                   check still fires, but on 60000 iterations against the
#                   same expected 37500 -- not 30001 -- and mg1>mg2>mg3 still
#                   held (mg3's residual fell to 4.96e-04, well BELOW mg2's,
#                   rather than rising above it), so the ordering check did
#                   NOT fire this time. Disagreeing with that earlier
#                   record is reported here rather than silently assumed to
#                   match. 60000 is exactly 2x mg3's level1-iter (30000) and
#                   the shared iter-threshold (30000) -- consistent with the
#                   near-zero coarse budget letting the multigrid cycle
#                   restart the fine grid for a second full pass rather than
#                   terminating after one, which would also explain the
#                   unexpectedly LOW residual. A concrete, testable
#                   hypothesis about the cycle-restart logic, not chased
#                   further here (out of this script's scope) -- either way the
#                   mechanism this hook exists to prove -- an input-only
#                   change reaching a genuine verify.py assertion, no solver
#                   code touched -- is demonstrated regardless of which of
#                   the four checks fires.
#               Measured green/red runs and the runtime are recorded at the
#               bottom of this header once both were executed (2026-09-15).
#
#               RUNTIME (measured 2026-09-15, -p 1, default set, this script's
#               own `time`, wall clock; solver-self-reported per case in
#               parens; binary md5 deea5f54... stable start-to-end):
#                 freestream 4s (0.070 min), conservation 3s (0.052 min),
#                 conjugate_two_layer 3s (0.054 min), receding N040+N080 18s
#                 (0.039+0.257 min), residual_convergence 6m16s total across
#                 its 6 sub-cases (nominal 0.849, irs 1.196, mg1 0.937,
#                 mg2 0.934, mg3 1.005, mg_irs 1.352 min) -- ~205k iterations
#                 dominate, as expected. TOTAL: real 6m49s.
#               This machine was under HEAVY external load during the
#               measurement (load average ~97 on 96 cores, other users'
#               jobs, unrelated to this script); the numbers
#               above are still representative because they are the solver's
#               own self-timed durations, not wall-clock deltas taken from
#               outside a contended scheduler.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NTHREADS=1
ALL=0
WORK=""
ONLY=""
SELFTEST="${FUSS_GATE_SELFTEST:-0}"

while test $# -gt 0; do
  case $1 in
    -p | --parallel) NTHREADS=$2; shift 2 ;;
    -w | --work)     WORK=$2;     shift 2 ;;
    -c | --case)     ONLY=$2;     shift 2 ;;
    -a | --all)      ALL=1; shift ;;
    -h | --help)
      echo "usage: $0 [-p NTHREADS] [-w WORKDIR] [-c CASE] [-a]"
      echo "  -c CASE   run only this case (freestream, conservation, bc303_qrad,"
      echo "            conjugate_two_layer, receding, residual_convergence, pyrolysis,"
      echo "            darcy_sealed_nosource, darcy_mass_balance_1d, darcy_sealed_profile,"
      echo "            darcy_tensor_linear_2d, darcy_0d, darcy_atc2_chemistry_off)"
      echo "  -a        also run receding's .slow N160/N320 grids"
      echo "  -w WORKDIR  scratch dir for the copied cases (default: a fresh"
      echo "              dir under \$TMPDIR/fuss-physics-gates.\$\$)"
      echo "  FUSS_GATE_SELFTEST=1  proven-red hook, see the file header"
      exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

[ -n "$WORK" ] || WORK="${TMPDIR:-/tmp}/fuss-physics-gates.$$"
mkdir -p "$WORK"

# Refuse to grade physics against a build that was not configured the way
# CMakePresets.json says -- same reasoning as run_regression.sh: a
# mis-configured build (serial, or TecIO-linked) is not loud on its own.
if ! "$ROOT/test/check_build_config.sh"; then
  echo
  echo "refusing to run physics gates against a build that does not match the preset"
  exit 1
fi

# Snapshot bin/FUSS EXACTLY ONCE, here, rather than letting `stage` re-read
# $ROOT/bin/FUSS for every case. bin/ can be rebuilt by a concurrent build while
# this script runs, and was observed CHURNING mid-run
# (caught at 0 bytes, mid-write, between two cases of the same invocation).
# Without this, two cases in the SAME run could silently grade two different
# executables. One snapshot makes "one invocation grades one binary" true by
# construction, not by luck.
BINSNAP="$WORK/.bin-snapshot"
mkdir -p "$BINSNAP"
cp "$ROOT/bin/FUSS" "$BINSNAP/FUSS"
BINMD5=$(md5sum "$BINSNAP/FUSS" | cut -d' ' -f1)

fail=0
KEYNUM=""
pass_line() { printf '  \033[0;32mPASS\033[0m  %-24s %s\n' "$1" "$2"; }
fail_line() { printf '  \033[0;31mFAIL\033[0m  %-24s %s\n' "$1" "$2"; fail=$((fail+1)); }

echo "FUSS physics gates"
echo "  root      : $ROOT"
echo "  threads   : $NTHREADS"
echo "  work      : $WORK"
echo "  binary    : $BINMD5  (bin/FUSS, snapshotted once at start)"
[ "$ALL" -eq 1 ] && echo "  receding  : all four grids (-a)"
[ "$SELFTEST" = "1" ] && echo "  SELFTEST  : FUSS_GATE_SELFTEST=1 -- tolerances deliberately broken"
echo

# ---------------------------------------------------------------------------
# stage SRC DEST -- minimal case copy: INPUT, MESH, input.ini, FUSS.sh, and a
# fresh OUTPUT/bin. See the header note on why SOLUTION/ (reference data) and
# any leftover OUTPUT/ are never copied.
#
# bin/FUSS is copied from $BINSNAP (the ONE snapshot taken at script start,
# not a fresh read of $ROOT/bin/FUSS per case -- see the snapshot comment
# above for why: bin/ can be rebuilt by a concurrent build, and it
# was observed CHURNING mid-run, 0 bytes mid-write
# at one sample. Reading it once and reusing that copy for every case is what
# makes "one invocation grades one binary" true regardless of what else
# touches $ROOT/bin/FUSS while this script runs. The work copy's mtime is
# pushed an hour into the future so FUSS.sh's own runtime
# `[[ "$MASTER" -nt "$LOCAL" ]] && cp` check (which still points at the LIVE,
# possibly-churning $ROOT/bin/FUSS via the sedded MASTERDIR) is a guaranteed
# no-op and can never overwrite the snapshot. MASTERDIR is still sedded to an
# absolute path (same discipline run_guard_demos.sh already uses) so the
# script's own logic is otherwise untouched.
# ---------------------------------------------------------------------------
stage () {
  local src=$1 dest=$2
  rm -rf "$dest"
  mkdir -p "$dest"
  cp -r "$src/INPUT" "$src/MESH" "$src/input.ini" "$src/FUSS.sh" "$dest/"
  # Case-root files some cases need at run time: the BC 302 time file
  # (Tw_time.dat, resolved relative to the run directory) and probes.json.
  cp "$src"/*.dat "$src"/*.json "$dest/" 2>/dev/null || true
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$dest/FUSS.sh"
  mkdir -p "$dest/OUTPUT" "$dest/bin"
  cp "$BINSNAP/FUSS" "$dest/bin/FUSS"
  touch -d '+1 hour' "$dest/bin/FUSS"
}

# solve DEST NTHREADS -- run FUSS.sh in DEST, log to DEST/run.log, return rc
solve () {
  local dest=$1 nt=$2 rc
  ( cd "$dest" && ./FUSS.sh -p "$nt" solve ) > "$dest/run.log" 2>&1
  rc=$?
  return $rc
}

# ---------------------------------------------------------------------------
# freestream
# ---------------------------------------------------------------------------
gate_freestream () {
  local src="$ROOT/test/moving-mesh/freestream" dest="$WORK/freestream" rc
  stage "$src" "$dest"
  solve "$dest" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then
    KEYNUM="solver exit $rc  (see $dest/run.log)"
    return 1
  fi

  local vpy="$src/verify.py"
  if [ "$SELFTEST" = "1" ]; then
    cp "$src/verify.py" "$dest/verify.selftest.py"
    sed -i 's/^TOL_ABS = .*/TOL_ABS = 0.0     # FUSS_GATE_SELFTEST: tightened to prove the gate can fail/' \
      "$dest/verify.selftest.py"
    vpy="$dest/verify.selftest.py"
  fi

  ( cd "$dest" && python3 "$vpy" ) > "$dest/verify.log" 2>&1
  local rc=$?
  KEYNUM=$(grep -F -m1 'max |T - T0|' "$dest/verify.log" | sed 's/^ *//')
  [ -n "$KEYNUM" ] || KEYNUM="(see $dest/verify.log)"
  return $rc
}

# ---------------------------------------------------------------------------
# conservation
# ---------------------------------------------------------------------------
gate_conservation () {
  local src="$ROOT/test/moving-mesh/conservation" dest="$WORK/conservation" rc
  stage "$src" "$dest"
  solve "$dest" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then
    KEYNUM="solver exit $rc  (see $dest/run.log)"
    return 1
  fi

  ( cd "$dest" && python3 "$src/verify.py" ) > "$dest/verify.log" 2>&1
  local rc=$?
  KEYNUM=$(grep -F -m1 'relative imbalance' "$dest/verify.log" | sed 's/^ *//')
  [ -n "$KEYNUM" ] || KEYNUM="(see $dest/verify.log)"
  return $rc
}

# ---------------------------------------------------------------------------
# bc303_qrad
# ---------------------------------------------------------------------------
gate_bc303_qrad () {
  local src="$ROOT/test/numerics-features/bc303_qrad" dest="$WORK/bc303_qrad" rc
  stage "$src" "$dest"
  solve "$dest" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then
    KEYNUM="solver exit $rc  (see $dest/run.log)"
    return 1
  fi

  # verify.py takes the case directory as its argument (it reads INPUT/bc.txt
  # for hconv/qrad/Tref and OUTPUT/wall.tec), so it runs unmodified on the copy.
  python3 "$src/verify.py" "$dest" > "$dest/verify.log" 2>&1
  local rc=$?
  KEYNUM=$(grep -F -m1 'worst residual' "$dest/verify.log" | sed 's/^ *//')
  [ -n "$KEYNUM" ] || KEYNUM="(see $dest/verify.log)"
  return $rc
}

# ---------------------------------------------------------------------------
# pyrolysis (plan 06 phase 2). One PASS/FAIL line for the group; the per-case
# verify logs sit under $WORK/pyrolysis/. The energy_audit variants are made
# from the STAGED copy (vnn sedded, the reaction lines dropped from
# pyrolysis.dat), so the gate never needs PATO or the table generator.
# ---------------------------------------------------------------------------
gate_pyrolysis () {
  local base="$ROOT/test/pyrolysis" dest="$WORK/pyrolysis" rc sub v
  rm -rf "$dest"; mkdir -p "$dest"
  KEYNUM=""

  # 4.2 isothermal kinetics against the closed-form propagator
  for sub in T0700 T0900 T1200; do
    stage "$base/kinetics_isothermal/$sub" "$dest/$sub"
    solve "$dest/$sub" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$sub: solver exit $rc  (see $dest/$sub/run.log)"; return 1; fi
    if ! python3 "$base/kinetics_isothermal/verify.py" "$dest/$sub" > "$dest/$sub/verify.log" 2>&1; then
      KEYNUM="$sub: $(grep -F -m1 'max |chi' "$dest/$sub/verify.log" | sed 's/^ *//')"; return 1
    fi
  done
  local chi_err; chi_err=$(grep -F -m1 'max |chi' "$dest/T1200/verify.log" | awk '{print $6}')

  # 4.3 energy audit: the case, two vnn refinements, the no-kinetics variant, 8 threads
  stage "$base/energy_audit" "$dest/energy_audit"
  solve "$dest/energy_audit" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="energy_audit: solver exit $rc  (see $dest/energy_audit/run.log)"; return 1; fi
  for v in 0.25 0.125; do
    stage "$base/energy_audit" "$dest/energy_audit_vnn$v"
    sed -i "s/^vnn = .*/vnn = $v/" "$dest/energy_audit_vnn$v/input.ini"
    solve "$dest/energy_audit_vnn$v" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="energy_audit vnn=$v: solver exit $rc"; return 1; fi
  done
  stage "$base/energy_audit" "$dest/energy_audit_nokin"
  python3 - "$dest/energy_audit_nokin/INPUT/pyrolysis.dat" <<'PYNOKIN'
import re, sys
p = sys.argv[1]
out = []
for line in open(p):
    key = line.split("//", 1)[0].split()
    k = key[0] if key else ""
    if re.match(r"^nPyroReac\[\d+\]$", k):
        out.append(k + " 0\n")
    elif re.match(r"^[FAEmnTh]\[\d+\]\[\d+\]$", k):
        continue
    else:
        out.append(line)
open(p, "w").write("".join(out))
PYNOKIN
  solve "$dest/energy_audit_nokin" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="energy_audit no-kinetics: solver exit $rc"; return 1; fi
  stage "$base/energy_audit" "$dest/energy_audit_p8"
  solve "$dest/energy_audit_p8" 8; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="energy_audit -p 8: solver exit $rc"; return 1; fi
  if ! python3 "$base/energy_audit/verify.py" "$dest/energy_audit" \
         --vnn025 "$dest/energy_audit_vnn0.25" --vnn0125 "$dest/energy_audit_vnn0.125" \
         --nokin "$dest/energy_audit_nokin" --p8 "$dest/energy_audit_p8" > "$dest/energy_audit/verify.log" 2>&1; then
    KEYNUM="energy_audit: $(grep -E -m1 'FAIL' "$dest/energy_audit/verify.log" | sed 's/^ *//')"; return 1
  fi
  local ebal; ebal=$(grep -F -m1 '(i)' "$dest/energy_audit/verify.log" | awk '{print $5}')

  # 4.6-4.8 exact-solution gates of the verification record 2026-10-08: adiabatic 0-D decomposition and
  # heated 0-D decomposition against the exact type-1 ODE (tools/ode0d.py), pure conduction on the ATC1 mesh
  # against the exact slab series (tools/slab_series.py)
  local sub src
  for sub in adiabatic0d/T1000 heated0d conduction/graded100; do
    src="$base/$sub"; stage "$src" "$dest/$(echo "$sub" | tr / _)"
    solve "$dest/$(echo "$sub" | tr / _)" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$sub: solver exit $rc  (see $dest/$(echo "$sub" | tr / _)/run.log)"; return 1; fi
    if ! python3 "$base/${sub%%/*}/verify.py" "$dest/$(echo "$sub" | tr / _)" > "$dest/$(echo "$sub" | tr / _)/verify.log" 2>&1; then
      KEYNUM="$sub: $(grep -E -m1 'FAIL' "$dest/$(echo "$sub" | tr / _)/verify.log" | sed 's/^ *//')"; return 1
    fi
  done
  local ad0 ht0 cnd
  ad0=$(grep -F -m1 'max |T - T_or1|' "$dest/adiabatic0d_T1000/verify.log" | awk '{print $6}')
  ht0=$(grep -F -m1 'max |T - T_or1|' "$dest/heated0d/verify.log" | awk '{print $6}')
  cnd=$(grep -F -m1 'max |T - exact|, all' "$dest/conduction_graded100/verify.log" | awk -F'= ' '{print $2}' | awk '{print $1}')

  # 4.4 ATC1 against FIAT (bundled)
  stage "$base/ATC1" "$dest/ATC1"
  solve "$dest/ATC1" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="ATC1: solver exit $rc  (see $dest/ATC1/run.log)"; return 1; fi
  if ! python3 "$base/ATC1/verify.py" "$dest/ATC1" > "$dest/ATC1/verify.log" 2>&1; then
    KEYNUM="ATC1: $(grep -E -m1 'FAIL|max over the gated' "$dest/ATC1/verify.log" | sed 's/^ *//')"; return 1
  fi
  # the gate is the move of the 100-cell discretisation error against the Richardson anchor (2026-10-08); FIAT is informational
  local atc1 atc1f
  atc1=$(grep -F -m1 'max over the gated points' "$dest/ATC1/verify.log" | awk -F'= ' '{print $2}' | awk '{print $1}')
  atc1f=$(grep -F -m1 'max over the six probes vs FIAT' "$dest/ATC1/verify.log" | awk -F': ' '{print $2}' | awk '{print $1}')

  # 4.5 ATC2-chemistryOff against PATO ($FUSS_PATO_ATC2_OUT, else the bundled ref/)
  stage "$base/ATC2-chemistryOff" "$dest/ATC2-chemistryOff"
  solve "$dest/ATC2-chemistryOff" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="ATC2-chemistryOff: solver exit $rc  (see $dest/ATC2-chemistryOff/run.log)"; return 1; fi
  if ! python3 "$base/ATC2-chemistryOff/verify.py" "$dest/ATC2-chemistryOff" > "$dest/ATC2-chemistryOff/verify.log" 2>&1; then
    KEYNUM="ATC2-chemistryOff: $(grep -F -m1 'max over' "$dest/ATC2-chemistryOff/verify.log" | sed 's/^ *//')"; return 1
  fi
  local atc2; atc2=$(grep -F -m1 'max over' "$dest/ATC2-chemistryOff/verify.log" | awk '{print $6}')
  grep -q '^SKIP' "$dest/ATC2-chemistryOff/verify.log" && atc2="SKIP"

  KEYNUM="chi err ${chi_err:-?}  energy balance ${ebal:-?}  0-D adiabatic ${ad0:-?} K  0-D heated ${ht0:-?} K  conduction ${cnd:-?} K  ATC1 move vs f* ${atc1:-?} K (vs FIAT ${atc1f:-?} K, info)  ATC2off max|dT| ${atc2:-?} K"
  return 0
}

# ---------------------------------------------------------------------------
# darcy (plan 07 phase 2b, gas-closure = type2): one gate per verification
# group of plan 07 section 4, cases under test/pyrolysis/darcy/. Every
# verify.py there honours FUSS_GATE_SELFTEST=1 (tolerances to 0 or below).
# ---------------------------------------------------------------------------
DARCY="$ROOT/test/pyrolysis/darcy"

# darcy_verify NAME SRCVERIFY DEST [ARGS...] -- run a verify.py on DEST, log to
# DEST/verify.log (or DEST/verify-NAME.log), return its exit code
darcy_verify () {
  local name=$1 vpy=$2 dest=$3; shift 3
  python3 "$vpy" "$dest" "$@" > "$dest/verify$name.log" 2>&1
}

# 4.1 (a) sealed, no source: p stays p-init exactly, zero flux and velocity
gate_darcy_sealed_nosource () {
  local src="$DARCY/sealed_nosource" dest="$WORK/darcy_sealed_nosource" rc
  stage "$src" "$dest"
  solve "$dest" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="solver exit $rc  (see $dest/run.log)"; return 1; fi
  darcy_verify "" "$src/verify.py" "$dest"; rc=$?
  KEYNUM="max|p-p0|/p0 $(grep -F -m1 '(a) max |p - p0|' "$dest/verify.log" | awk '{print $(NF-3)}')  max|mdot| $(grep -F -m1 '(c) max |mdot_pg|' "$dest/verify.log" | awk '{print $(NF-3)}')"
  return $rc
}

# 4.1 (b) steady gas mass balance on a 1-D line: (b2) the case as committed,
# time-accurate; (b1) its staged copy with quasi-steady = true
gate_darcy_mass_balance_1d () {
  local src="$DARCY/mass_balance_1d" dest="$WORK/darcy_mass_balance_1d" rc v b1 b2
  for v in b2 b1; do
    stage "$src" "$dest/$v"
    [ "$v" = "b1" ] && sed -i 's/^quasi-steady = .*/quasi-steady = true/' "$dest/$v/input.ini"
    solve "$dest/$v" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$v: solver exit $rc  (see $dest/$v/run.log)"; return 1; fi
    if ! darcy_verify "" "$src/verify.py" "$dest/$v"; then
      KEYNUM="$v: $(grep -E -m1 'FAIL' "$dest/$v/verify.log" | sed 's/^ *//')"; return 1
    fi
  done
  b2=$(grep -F -m1 '(i)   max_f' "$dest/b2/verify.log" | awk '{print $(NF-3)}')
  b1=$(grep -F -m1 '(i)   max_f' "$dest/b1/verify.log" | awk '{print $(NF-3)}')
  KEYNUM="max_f|mdot-pi x|/(pi L): b2 ${b2:-?} (tol 1e-6)  b1 ${b1:-?} (tol 1e-9)"
  return 0
}

# 4.2 analytic steady pressure profile on four grids, observed order >= 1.9
gate_darcy_sealed_profile () {
  local src="$DARCY/sealed_profile" dest="$WORK/darcy_sealed_profile" rc g
  rm -rf "$dest"; mkdir -p "$dest"
  for g in N020 N040 N080 N160; do
    stage "$src/$g" "$dest/$g"
    solve "$dest/$g" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$g: solver exit $rc  (see $dest/$g/run.log)"; return 1; fi
  done
  darcy_verify "" "$src/verify.py" "$dest"; rc=$?
  KEYNUM=$(grep -F -m1 'order on the finest pair' "$dest/verify.log" | sed 's/^ *//')
  [ -n "$KEYNUM" ] || KEYNUM="(see $dest/verify.log)"
  return $rc
}

# Klinkenberg oracle (plan 07, closed 2026-10-13): the sealed_profile grids with
# [FUSS-Darcy] klinkenberg = true and B_v = B_c = B I, B = b K/mu, b = p_w/2;
# p + b then obeys the plain equation, (p+b)^2 = (p_w+b)^2 + (pi R T mu/(W K))(L^2 - x^2)
gate_darcy_klinkenberg () {
  local src="$DARCY/sealed_profile" dest="$WORK/darcy_klinkenberg" rc g
  rm -rf "$dest"; mkdir -p "$dest"
  for g in N020 N040 N080 N160; do
    stage "$src/$g" "$dest/$g"
    python3 - "$dest/$g" "$dest" "$ROOT/test/pyrolysis/tools" > "$dest/$g.edit.log" 2>&1 <<'PY' || { KEYNUM="$g: Klinkenberg edit failed (see $dest/$g.edit.log)"; return 1; }
import json, os, sys
case, parent, tools = sys.argv[1:4]
sys.path.insert(0, tools)
import darcy_tools as dt
meta = dt.case_meta(case)
W, mu, h = dt.gas_constants(case)
b = 0.5 * meta["p0"]
B = b * meta["K"] / mu
with open(os.path.join(case, "INPUT", "pyrolysis.dat"), "a") as fh:
    fh.write("B_v %.15g 0 0 0 %.15g 0 0 0 %.15g\nB_c %.15g 0 0 0 %.15g 0 0 0 %.15g\n" % ((B,) * 6))
ini = os.path.join(case, "input.ini")
s = open(ini).read()
assert s.count("[FUSS-Darcy]\n") == 1
open(ini, "w").write(s.replace("[FUSS-Darcy]\n", "[FUSS-Darcy]\nklinkenberg = true\n"))
json.dump({"b": b, "B": B}, open(os.path.join(parent, "klinkenberg.json"), "w"))
print("b %.6g Pa, B %.6g m2/s" % (b, B))
PY
    solve "$dest/$g" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$g: solver exit $rc  (see $dest/$g/run.log)"; return 1; fi
  done
  darcy_verify "" "$src/verify_klinkenberg.py" "$dest"; rc=$?
  KEYNUM="$(grep -F -m1 '(a) observed orders' "$dest/verify.log" | sed -E 's/^ *\(a\) //; s/  (PASS|FAIL)$//')  $(grep -F -m1 '(b) Klinkenberg' "$dest/verify.log" | sed -E 's/^ *\(b\) Klinkenberg moves the profile by ([^ ]+) Pa.*/shift \1 Pa/')"
  return $rc
}

# 4.3 orthotropic permeability: |u_y|/|u_x| = K_yx/K_xx in every cell; the red
# control (K diagonal, written into the staged copy) must give u_y = 0 to the
# solver tolerance
gate_darcy_tensor_linear_2d () {
  local src="$DARCY/tensor_linear_2d" dest="$WORK/darcy_tensor_linear_2d" rc
  stage "$src" "$dest/rotated"
  solve "$dest/rotated" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="rotated: solver exit $rc  (see $dest/rotated/run.log)"; return 1; fi
  if ! darcy_verify "" "$src/verify.py" "$dest/rotated"; then
    KEYNUM="rotated: $(grep -E -m1 'FAIL' "$dest/rotated/verify.log" | sed 's/^ *//')"; return 1
  fi
  stage "$src" "$dest/control"
  python3 - "$dest/control" <<'PYCTRL'
import json, os, sys
case = sys.argv[1]
meta = json.load(open(os.path.join(case, "case.json")))
k1, k2 = meta["K1"], meta["K2"]
K = [k1, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, k2]
p = os.path.join(case, "INPUT", "pyrolysis.dat")
out, seen = [], set()
for s in open(p):
    key = s.split("//", 1)[0].split()
    k = key[0] if key else ""
    if k in ("K_v", "K_c"):
        out.append("%s %s\n" % (k, " ".join("%.15g" % v for v in K)))
        seen.add(k)
    else:
        out.append(s)
if seen != {"K_v", "K_c"}:
    sys.exit("control: K_v/K_c lines not found in %s" % p)
open(p, "w").write("".join(out))
meta["identity"], meta["K"] = True, K
json.dump(meta, open(os.path.join(case, "case.json"), "w"), indent=1, sort_keys=True)
PYCTRL
  solve "$dest/control" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="control: solver exit $rc  (see $dest/control/run.log)"; return 1; fi
  if ! darcy_verify "" "$src/verify.py" "$dest/control" --control; then
    KEYNUM="control: $(grep -E -m1 'FAIL' "$dest/control/verify.log" | sed 's/^ *//')"; return 1
  fi
  KEYNUM="ratio err $(grep -F -m1 '(a) max | |u_y|/|u_x|' "$dest/rotated/verify.log" | awk '{print $(NF-3)}') (tol 1e-10)  control max|u_y|/max|u_x| $(grep -F -m1 '(a) max |u_y| / max |u_x|' "$dest/control/verify.log" | awk '{print $(NF-3)}')"
  return 0
}

# 4.1 (d) the two 0-D decompositions under type 2 against the exact type-2 ODE
gate_darcy_0d () {
  local dest="$WORK/darcy_0d" rc sub ad ht
  rm -rf "$dest"; mkdir -p "$dest"
  for sub in adiabatic0d_type2 heated0d_type2; do
    stage "$DARCY/$sub" "$dest/$sub"
    solve "$dest/$sub" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$sub: solver exit $rc  (see $dest/$sub/run.log)"; return 1; fi
    if ! darcy_verify "" "$DARCY/$sub/verify.py" "$dest/$sub"; then
      KEYNUM="$sub: $(grep -E -m1 'FAIL' "$dest/$sub/verify.log" | sed 's/^ *//')"; return 1
    fi
  done
  ad=$(grep -F -m1 'max |T - T_or2|' "$dest/adiabatic0d_type2/verify.log" | awk '{print $6}')
  ht=$(grep -F -m1 'max |T - T_or2|' "$dest/heated0d_type2/verify.log" | awk '{print $6}')
  KEYNUM="0-D adiabatic ${ad:-?} K (tol 0.1)  0-D heated ${ht:-?} K (tol 0.02) vs the exact type-2 ODE"
  return 0
}

# 4.4 ATC2-chemistryOff under type 2 against the bundled PATO run (T and p
# probed), and against the type-1 case on the same mesh (the sign check)
gate_darcy_atc2_chemistry_off () {
  local dest="$WORK/darcy_atc2_chemistry_off" rc
  rm -rf "$dest"; mkdir -p "$dest"
  stage "$ROOT/test/pyrolysis/ATC2-chemistryOff" "$dest/type1"
  solve "$dest/type1" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="type1: solver exit $rc  (see $dest/type1/run.log)"; return 1; fi
  stage "$DARCY/ATC2-chemistryOff-type2" "$dest/type2"
  solve "$dest/type2" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="type2: solver exit $rc  (see $dest/type2/run.log)"; return 1; fi
  darcy_verify "" "$DARCY/ATC2-chemistryOff-type2/verify.py" "$dest/type2" --type1 "$dest/type1"; rc=$?
  if grep -q '^SKIP' "$dest/type2/verify.log"; then KEYNUM="SKIP (no PATO probe files)"; return 0; fi
  KEYNUM="max|dT| $(grep -F -m1 '(T)   max over' "$dest/type2/verify.log" | awk '{print $7}') K (budget $(grep -F -m1 '(T)   max over' "$dest/type2/verify.log" | sed -E 's/.*budget ([0-9.]+) K.*/\1/'), B_T 22.75)  dp ratio $(grep -F -m1 '        ratio' "$dest/type2/verify.log" | awk '{print $2}') (tol 0.15)  type2-type1 hot max $(grep -F -m1 '(sign) max' "$dest/type2/verify.log" | awk '{print $10}') K"
  return $rc
}

# ---------------------------------------------------------------------------
# conjugate_two_layer (needs generate.py alongside verify.py:
# verify.py imports the analytic profile and geometry constants from it)
# ---------------------------------------------------------------------------
gate_conjugate_two_layer () {
  local src="$ROOT/test/numerics-features/conjugate_two_layer" dest="$WORK/conjugate_two_layer" rc
  stage "$src" "$dest"
  cp "$src/generate.py" "$src/verify.py" "$dest/"
  solve "$dest" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then
    KEYNUM="solver exit $rc  (see $dest/run.log)"
    return 1
  fi

  python3 "$dest/verify.py" > "$dest/verify.log" 2>&1
  local rc=$?
  KEYNUM=$(grep -F -m1 'interface error' "$dest/verify.log" | sed 's/^ *//')
  [ -n "$KEYNUM" ] || KEYNUM="(see $dest/verify.log)"
  return $rc
}

# ---------------------------------------------------------------------------
# receding. verify.py locates itself via __file__ (HERE), so it
# is invoked in place in $dest; GRIDS is overridden by sedding the work copy
# of generate.py (see the header note on -a).
# ---------------------------------------------------------------------------
gate_receding () {
  local src="$ROOT/test/moving-mesh/receding" dest="$WORK/receding"
  rm -rf "$dest"
  mkdir -p "$dest"
  cp "$src/generate.py" "$src/verify.py" "$dest/"

  local grids
  if [ "$ALL" -eq 1 ]; then
    grids="40 80 160 320"
    sed -i 's/^GRIDS = .*/GRIDS = [40, 80, 160, 320]/' "$dest/generate.py"
  else
    grids="40 80"
    sed -i 's/^GRIDS = .*/GRIDS = [40, 80]/' "$dest/generate.py"
  fi

  local n nd rc
  for n in $grids; do
    nd=$(printf 'N%03d' "$n")
    stage "$src/$nd" "$dest/$nd"
    solve "$dest/$nd" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then
      KEYNUM="$nd: solver exit $rc  (see $dest/$nd/run.log)"
      return 1
    fi
  done

  python3 "$dest/verify.py" > "$dest/verify.log" 2>&1
  local rc=$? decay_row order_row reldiff order_val
  # Pulled from the two per-grid TABLES, not the wrapped PASS/FAIL prose
  # (verify.py splits "convergence order N.NNN on the finest pair (first
  # order," across two print() calls, which truncates under a single grep).
  # Decay table rows: N Pe d_measured d_predicted rel.diff pts  (6 fields).
  # Pointwise table rows: N max|dT| L2 order budget meas/pred at-xi/d (7).
  decay_row=$(awk '/^  [0-9]+ / && NF==6 {line=$0} END{print line}' "$dest/verify.log")
  order_row=$(awk '/^  [0-9]+ / && NF==7 {line=$0} END{print line}' "$dest/verify.log")
  reldiff=$(awk '{print $5}' <<< "$decay_row")
  order_val=$(awk '{print $4}' <<< "$order_row")
  KEYNUM="order=${order_val:-?} (band 0.85-1.15)  finest-grid decay rel.diff=${reldiff:-?} (tol 0.02)"
  return $rc
}

# ---------------------------------------------------------------------------
# residual_convergence. Six sub-cases sharing one verify.py,
# which locates itself via __file__ and reads $dest/<subcase>/OUTPUT.
# ---------------------------------------------------------------------------
gate_residual_convergence () {
  local src="$ROOT/test/numerics-features/residual_convergence" dest="$WORK/residual_convergence"
  rm -rf "$dest"
  mkdir -p "$dest"
  cp "$src/verify.py" "$dest/"

  local sub rc
  for sub in nominal irs mg1 mg2 mg3 mg_irs; do
    stage "$src/$sub" "$dest/$sub"
    if [ "$SELFTEST" = "1" ] && [ "$sub" = "mg3" ]; then
      # The proven-red case: take away mg3's coarse-grid work.
      # Fires the budget check (measured here: 60000 iterations vs an
      # expected 37500 -- see the file header for how this differs from the
      # original recorded numbers). Input-only, no solver code touched.
      sed -i 's/^level2-iter = .*/level2-iter = 1   ; FUSS_GATE_SELFTEST: mg3 loses its coarse work/' \
        "$dest/$sub/input.ini"
    fi
    solve "$dest/$sub" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then
      KEYNUM="$sub: solver exit $rc  (see $dest/$sub/run.log)"
      return 1
    fi
  done

  python3 "$dest/verify.py" > "$dest/verify.log" 2>&1
  local rc=$?
  KEYNUM=$(grep -E '^  (nominal|irs|mg1|mg2|mg3|mg_irs) ' "$dest/verify.log" \
             | awk '{printf "%s=%s ", $1, $4}')
  [ -n "$KEYNUM" ] || KEYNUM="(see $dest/verify.log)"
  return $rc
}

# ---------------------------------------------------------------------------
# bench_* (2026-10-13): the benchmarks of the type-2 closure, against exact
# solutions (test/pyrolysis/exact/) and the literature (test/pyrolysis/lit/).
# The mesh ladders are generated into $WORK by the committed generators
# (test/pyrolysis/tools/make_*.py); the in-tree case of each benchmark is its
# coarsest member, for the regression manifest. Every verify honours
# FUSS_GATE_SELFTEST=1 (budgets to 0, order windows empty).
# ---------------------------------------------------------------------------
TOOLS="$ROOT/test/pyrolysis/tools"

# gen_case GENERATOR DEST ARGS... -- write a case with a committed generator
# (absolute MASTERDIR), then give it the snapshot binary as `stage` does
gen_case () {
  local gen=$1 dest=$2; shift 2
  rm -rf "$dest"
  python3 "$TOOLS/$gen" "$dest" "$@" --masterdir "$ROOT" > "$dest.gen.log" 2>&1 || return 1
  mkdir -p "$dest/OUTPUT" "$dest/bin"
  cp "$BINSNAP/FUSS" "$dest/bin/FUSS"
  touch -d '+1 hour' "$dest/bin/FUSS"
}

# bench_ladder KIND NAME N... -- generate and run an exact-solution ladder, then verify_exact.py
bench_ladder () {
  local kind=$1 name=$2; shift 2
  local dest="$WORK/$name" n rc cases=()
  rm -rf "$dest"; mkdir -p "$dest"
  for n in "$@"; do
    gen_case make_exact_case.py "$dest/N$n" --kind "$kind" --tables "$ROOT/test/pyrolysis/exact/tables" --ncells "$n" \
      || { KEYNUM="N$n: generator failed (see $dest/N$n.gen.log)"; return 1; }
    solve "$dest/N$n" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="N$n: solver exit $rc  (see $dest/N$n/run.log)"; return 1; fi
    cases+=("$dest/N$n")
  done
  python3 "$ROOT/test/pyrolysis/exact/verify_exact.py" "$kind" "${cases[@]}" > "$dest/verify.log" 2>&1; rc=$?
  KEYNUM=$(grep -E '(finest|order)' "$dest/verify.log" | sed -E 's/  (PASS|FAIL)$//; s/ on the finest (grid|pair)//' | tr -s ' ' | paste -sd ';' -)
  return $rc
}

# exact transient Darcy pressure (frozen T): 2nd order in space at vnn 0.5
gate_bench_pressure_transient () { bench_ladder pressure_transient bench_pressure_transient 20 40 80; }

# exact steady transpiration (exponential T, Darcy mdot): 2nd order since the
# open-face gas enthalpy is taken at the wall temperature (was 1st order)
gate_bench_transpiration () { bench_ladder transpiration bench_transpiration 20 40 80; }

# exact sharp front with transpiration through the char: 1st order by design
gate_bench_sharp_front () { bench_ladder sharp_front bench_sharp_front 50 100 200; }

# Baum & Atreya 2007 (corrected): full type 2 on N100/N200 (front, T, char
# pressure, gas-mass ledger to round-off; the surface steps to 1298 K at t = 0,
# which stopped the run before the storage-coefficient limiter) and the
# quasi-steady N100 gas ledger (one carried step; 11 % lost before the carry)
gate_bench_baum_atreya () {
  local dest="$WORK/bench_baum_atreya" c rc
  rm -rf "$dest"; mkdir -p "$dest"
  gen_case make_baum_atreya.py "$dest/full_N100" --mode full --cg 1475 --ncells 100 || { KEYNUM="generator failed"; return 1; }
  gen_case make_baum_atreya.py "$dest/full_N200" --mode full --cg 1475 --ncells 200 || { KEYNUM="generator failed"; return 1; }
  gen_case make_baum_atreya.py "$dest/qs_N100" --mode qs --cg 0 --ncells 100 || { KEYNUM="generator failed"; return 1; }
  for c in full_N100 full_N200 qs_N100; do
    solve "$dest/$c" "$NTHREADS"; rc=$?
    if [ "$rc" -ne 0 ]; then KEYNUM="$c: solver exit $rc  (see $dest/$c/run.log)"; return 1; fi
  done
  python3 "$ROOT/test/pyrolysis/lit/baum_atreya2007/verify.py" "$dest/full_N100" "$dest/full_N200" --qs "$dest/qs_N100" \
    > "$dest/verify.log" 2>&1; rc=$?
  KEYNUM=$(grep -E '^  \((front|T|p|qs)\)' "$dest/verify.log" | sed -E 's/  (PASS|FAIL)$//; s/^  //' | tr -s ' ' | paste -sd ';' -)
  return $rc
}

# Zhang 2010 three-code decomposing slab (BC 305 front face, type 1): within
# the codes' own spread of the MATLAB code (7.4 K off before the BC 305 fix)
gate_bench_zhang2010 () {
  local dest="$WORK/bench_zhang2010" rc
  rm -rf "$dest"; mkdir -p "$dest"
  gen_case make_zhang2010.py "$dest/N050" --ncells 50 || { KEYNUM="generator failed"; return 1; }
  solve "$dest/N050" "$NTHREADS"; rc=$?
  if [ "$rc" -ne 0 ]; then KEYNUM="solver exit $rc  (see $dest/N050/run.log)"; return 1; fi
  python3 "$ROOT/test/pyrolysis/lit/zhang2010/verify.py" "$dest/N050" > "$dest/verify.log" 2>&1; rc=$?
  KEYNUM=$(grep -E '^  \((T|mdot)\)' "$dest/verify.log" | sed -E 's/\(T\) +([a-z]+) +max\|T - T_MATLAB\| ([0-9.]+) K.*/\1 \2 K/; s/^  //; s/  (PASS|FAIL)$//' \
           | cut -c1-60 | paste -sd ';' -)
  return $rc
}

# ---------------------------------------------------------------------------
CASES=(freestream conservation bc303_qrad conjugate_two_layer receding residual_convergence pyrolysis
       darcy_sealed_nosource darcy_mass_balance_1d darcy_sealed_profile darcy_klinkenberg darcy_tensor_linear_2d darcy_0d
       darcy_atc2_chemistry_off
       bench_pressure_transient bench_transpiration bench_sharp_front bench_baum_atreya bench_zhang2010)

if [ -n "$ONLY" ]; then
  case " ${CASES[*]} " in
    *" $ONLY "*) CASES=("$ONLY") ;;
    *) echo "unknown case: $ONLY (known: ${CASES[*]})" >&2; exit 1 ;;
  esac
fi

for c in "${CASES[@]}"; do
  KEYNUM=""
  if "gate_$c"; then
    pass_line "$c" "$KEYNUM"
  else
    fail_line "$c" "$KEYNUM"
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "all physics gates passed"
else
  echo "$fail physics gate(s) FAILED -- see the per-case verify.log / run.log under $WORK"
  exit 1
fi
