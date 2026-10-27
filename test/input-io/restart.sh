#!/bin/bash
# Restart: a run split in two (N iterations, then a restart that runs N more) writes the files of the run in one go (2N
# iterations), byte for byte: the solution stores the counters that decide the next iteration (CYCLE and FUSS_ITERATION,
# FUSS_OUTPUT_NUMBER, FUSS_GRID_LEVEL, FUSS_LAST_OUTPUT_TIME; vtk field data, Tecplot data set auxiliary data) and the
# restart reads them back. Also: the stops of a solution without a time or with a counter that is not finite, a solution
# without counters (legacy path), the rows written after the solution, a time-accurate run started from a steady
# solution, the source INPUT/st.* read with ic-format.
#   usage: restart.sh <FUSS binary> <work folder> [core|szplt|all]     (one CASE line per check, the battery reads them)
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); GROUP=${3:-core}; source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

R=numerics-features/residual_convergence
iters() {  # iters <dir> <n>: the iteration budget of the run (level 1 with multigrid levels)
  if grep -qE '^levels[[:space:]]*=[[:space:]]*[2-9]' "$1/input.ini"; then setkey "$1" level1-iter "$2" FUSS-Multigrid
  else setkey "$1" iter-threshold "$2"; fi
}
fmt() { setkey "$1" sol-format "$2" FUSS-IO; }

# split <name> <case> <sol-format> <N> [extra edits as 'key value section' triples applied to every run]
split() {
  local n=$1 c=$2 f=$3 N=$4 a="$WORK/$1-A" b="$WORK/$1-B"; shift 4
  prep "$c" "$a"; prep "$c" "$b"
  for d in "$a" "$b"; do setkey "$d" res-threshold 1e-30; fmt "$d" "$f"; done
  while [[ $# -gt 0 ]]; do setkey "$a" "$1" "$2" "$3"; setkey "$b" "$1" "$2" "$3"; shift 3; done
  iters "$a" $((2*N)); run "$a"
  iters "$b" "$N"; run "$b"; setkey "$b" newrun false; run "$b"
  if [[ $(code "$a") != 0 || $(code "$b") != 0 ]]; then case_result "$n" 1 "exit A $(code "$a") B $(code "$b")"; return; fi
  local nd; nd=$(ndiff "$a" "$b")
  if [[ $nd == 0 ]]; then case_result "$n" 0 "0 differing files"; else case_result "$n" 1 "$nd differing files"; fi
}

# split_time <name> <case> <sol-format> <T>: A runs to time T in one go, B to T/2 and then (restart) to T
split_time() {
  local n=$1 c=$2 f=$3 T=$4 a="$WORK/$1-A" b="$WORK/$1-B"; shift 4
  prep "$c" "$a"; prep "$c" "$b"
  for d in "$a" "$b"; do fmt "$d" "$f"; setkey "$d" time-threshold "$T"; done
  while [[ $# -gt 0 ]]; do setkey "$a" "$1" "$2" "$3"; setkey "$b" "$1" "$2" "$3"; shift 3; done
  run "$a"
  setkey "$b" time-threshold "$(awk -v t="$T" 'BEGIN{print t/2}')"; run "$b"; setkey "$b" newrun false; setkey "$b" time-threshold "$T"; run "$b"
  if [[ $(code "$a") != 0 || $(code "$b") != 0 ]]; then case_result "$n" 1 "exit A $(code "$a") B $(code "$b")"; return; fi
  local nd; nd=$(ndiff "$a" "$b")
  if [[ $nd == 0 ]]; then case_result "$n" 0 "0 differing files"; else case_result "$n" 1 "$nd differing files"; fi
}

if [[ $GROUP == core || $GROUP == all ]]; then
  split nominal-vtk-raw     $R/nominal "vtk raw" 200
  split nominal-tecplot     $R/nominal "tecplot ascii" 200
  # multigrid: the restart goes on at level 1 (FUSS_GRID_LEVEL), not from the coarsest level
  split mg-vtk-raw          $R/mg_irs  "vtk raw" 200 level2-iter 300 FUSS-Multigrid
  # FAS cycles: no coarse state between cycles. A run cut by its iteration budget inside a cycle restarts from the cut
  # state, not from the cycle: the split is on a cycle boundary (a cycle takes 3 iterations here)
  split fas-vtk-raw         $R/fas_w   "vtk raw" 99
  # time-based probes and outputs (sol-dtime), time-accurate: the split at half the time threshold
  split_time probes-vtk-raw transient/oscillating_Rod "vtk raw" 160
  # numbered solutions: the number goes on, no doubled rows (the split on 22, between the outputs at 20 and 25)
  split numbered-tecplot    transient/variable_cp_Plate "tecplot ascii" 22

  # numbered solutions, the split on an output iteration: the output due at the last iteration is written (commit 6)
  split numbered-due-tecplot transient/variable_cp_Plate "tecplot ascii" 20

  # the first part ends on a due output: by sol-diter with sol-overwrite = true, with two grid levels, by sol-dtime only
  split due-overwrite-vtk-raw $R/nominal "vtk raw" 200 sol-diter 100 FUSS-IO
  split due-mg-vtk-raw        $R/mg_irs  "vtk raw" 200 level2-iter 300 FUSS-Multigrid sol-diter 100 FUSS-IO
  split_time due-dtime-vtk-raw transient/oscillating_Rod "vtk raw" 160 sol-dtime 0.01 FUSS-IO

  # a solution without counters (written before them): the restart goes on from its time, counts from there
  prep $R/nominal "$WORK/legacy"; setkey "$WORK/legacy" res-threshold 1e-30; iters "$WORK/legacy" 100; run "$WORK/legacy"
  sed -i '/^DATASETAUXDATA/d' "$WORK/legacy/OUTPUT/field.tec"
  setkey "$WORK/legacy" newrun false; iters "$WORK/legacy" 100; run "$WORK/legacy"
  last=$(tail -1 "$WORK/legacy/OUTPUT/residual-history.dat" | awk '{print $1}')
  if [[ $(code "$WORK/legacy") == 0 && $last == 200 ]]; then case_result legacy-no-counters 0 "restart from 100, last iteration 200"
  else case_result legacy-no-counters 1 "exit $(code "$WORK/legacy"), last iteration $last"; fi

  # rows after the solution: the run went on past the solution the restart reads (field2 at iteration 20): its rows are dropped
  a="$WORK/rows-A"; b="$WORK/rows-B"
  for d in "$a" "$b"; do prep $R/nominal "$d"; setkey "$d" res-threshold 1e-30; setkey "$d" sol-diter 10 FUSS-IO; setkey "$d" sol-overwrite false FUSS-IO; done
  iters "$a" 40; run "$a"
  iters "$b" 25; run "$b"; rm "$b/OUTPUT/field.tec"; setkey "$b" newrun false; iters "$b" 20; run "$b"
  nd=$(ndiff "$a" "$b"); if [[ $(code "$b") == 0 && $nd == 0 ]]; then case_result rows-after-solution 0 "0 differing files"
  else case_result rows-after-solution 1 "exit $(code "$b"), $nd differing files"; fi

  # a time-accurate run started from a steady solution starts at t = 0 with new counters: it is the run from the same
  # field as the initial condition
  prep steady-state/temperature_Plate "$WORK/tfs-B"; setkey "$WORK/tfs-B" iter-threshold 50; run "$WORK/tfs-B"
  prep steady-state/temperature_Plate "$WORK/tfs-C"; cp "$WORK/tfs-B/OUTPUT/field.tec" "$WORK/tfs-C/INPUT/ic.tec"
  for d in "$WORK/tfs-B" "$WORK/tfs-C"; do setkey "$d" time-accurate true FUSS-Numerics; setkey "$d" vnn 0.25 FUSS-Numerics; setkey "$d" iter-threshold 5; done
  setkey "$WORK/tfs-B" newrun false; run "$WORK/tfs-B"; run "$WORK/tfs-C"
  first=$(awk '$1==1' "$WORK/tfs-B/OUTPUT/residual-history.dat" | tail -1 | awk '{print $2}')
  cmp -s "$WORK/tfs-B/OUTPUT/field.tec" "$WORK/tfs-C/OUTPUT/field.tec" && same=1 || same=0
  if [[ $(code "$WORK/tfs-B") == 0 && $same == 1 ]] && awk -v t="$first" 'BEGIN{exit !(t+0 > 0 && t+0 < 10)}'; then case_result transient-from-steady 0 "starts at t = 0 (first row time $first), field = run from ic.tec"
  else case_result transient-from-steady 1 "exit $(code "$WORK/tfs-B"), first row time $first, same field $same"; fi

  # the source INPUT/st.* is read with ic-format: qvol_Plate restarted in vtk raw (ic-format tecplot ascii, st.tec)
  prep steady-state/qvol_Plate "$WORK/src-A"; prep steady-state/qvol_Plate "$WORK/src-B"
  for d in "$WORK/src-A" "$WORK/src-B"; do setkey "$d" res-threshold 1e-30; fmt "$d" "vtk raw"; done
  iters "$WORK/src-A" 200; run "$WORK/src-A"
  iters "$WORK/src-B" 100; run "$WORK/src-B"; setkey "$WORK/src-B" newrun false; run "$WORK/src-B"
  nd=$(ndiff "$WORK/src-A" "$WORK/src-B")
  if [[ $(code "$WORK/src-A") == 0 && $(code "$WORK/src-B") == 0 && $nd == 0 ]]; then case_result source-ic-format 0 "0 differing files"
  else case_result source-ic-format 1 "exit $(code "$WORK/src-B"), $nd differing files"; fi

  # the .vtm of the field, wall and diagnostic files list their block files relative to their own folder (vtk/<name>.vts)
  prep $R/nominal "$WORK/vtm"; iters "$WORK/vtm" 10; fmt "$WORK/vtm" "vtk raw"; run "$WORK/vtm"
  bad=$(grep -h 'file="' "$WORK"/vtm/OUTPUT/*.vtm | grep -vc 'file="vtk/[A-Za-z0-9_-]*\.vts"')
  nv=$(grep -h 'file="' "$WORK"/vtm/OUTPUT/*.vtm | wc -l)
  if [[ $(code "$WORK/vtm") == 0 && $bad == 0 && $nv -gt 0 ]]; then case_result vtm-entries 0 "$nv entries relative to the .vtm"; else case_result vtm-entries 1 "exit $(code "$WORK/vtm"), $bad entries of $nv not vtk/<name>.vts"; fi

  # a vtk solution without TIME, and one with a counter that is not finite, stop at set-up; the vts in ascii are editable
  prep $R/nominal "$WORK/bad"; setkey "$WORK/bad" res-threshold 1e-30; iters "$WORK/bad" 10; fmt "$WORK/bad" "vtk ascii"; run "$WORK/bad"
  cp -a "$WORK/bad" "$WORK/bad-notime"; cp -a "$WORK/bad" "$WORK/bad-nan"
  for f in "$WORK"/bad-notime/OUTPUT/vtk/field*.vts; do python3 - "$f" notime <<'PY'
import re, sys
f, mode = sys.argv[1:3]; t = open(f).read()
t = re.sub(r'<DataArray[^>]*Name="TIME"[^>]*>.*?</DataArray>', '', t, flags=re.S)
open(f, 'w').write(t)
PY
  done
  for f in "$WORK"/bad-nan/OUTPUT/vtk/field*.vts; do python3 - "$f" nan <<'PY'
import re, sys
f = sys.argv[1]; t = open(f).read()
t = re.sub(r'(<DataArray[^>]*Name="FUSS_GRID_LEVEL"[^>]*>)\s*[-0-9.eE+]+', r'\1 NaN', t, flags=re.S)
open(f, 'w').write(t)
PY
  done
  for n in bad-notime bad-nan; do setkey "$WORK/$n" newrun false; run "$WORK/$n"; done
  grep -q "OUTPUT" /dev/null
  c=$(code "$WORK/bad-notime"); if [[ $c == 1 ]] && grep -q 'the solution holds no time' "$WORK/bad-notime/run.log"; then case_result vtk-no-time 0 "stop at set-up"; else case_result vtk-no-time 1 "exit $c"; fi
  c=$(code "$WORK/bad-nan"); if [[ $c == 1 ]] && grep -q 'FUSS_GRID_LEVEL is not finite' "$WORK/bad-nan/run.log"; then case_result vtk-counter-nan 0 "stop at set-up"; else case_result vtk-counter-nan 1 "exit $c"; fi
fi

if [[ $GROUP == szplt || $GROUP == all ]]; then
  # .szplt: the zone time (signed) and the counters in the binary aux data
  split nominal-szplt $R/nominal "tecplot binary" 200
fi

finish
