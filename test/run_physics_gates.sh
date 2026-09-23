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
#                 pyrolysis   -- add gate_pyrolysis()  + list it
#                 darcy       -- add gate_darcy()      + list it
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
      echo "            conjugate_two_layer, receding, residual_convergence)"
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
      sed -i 's/^level2-iter = .*/level2-iter = 1   # FUSS_GATE_SELFTEST: mg3 loses its coarse work/' \
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
CASES=(freestream conservation bc303_qrad conjugate_two_layer receding residual_convergence)

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
