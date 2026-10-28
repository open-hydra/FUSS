#!/bin/bash
# MPIExitCodes: with two ranks every error stop ends the run with exit code 1, the STOP line of both ranks where the
# decision is collective, and no hang. The errors: an input error found by rank 0 (broadcast), a bc.txt error (every
# rank reads it), an initial-condition read error on rank 1 only (MPMD: its folder has no ic.tec), a value of rank 1's
# initial condition that is not finite, the reading of input.ini during the run with ini-error = stop, a thin block next to a connection.
#   usage: mpi-exit-codes.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

base() { prep numerics-features/residual_convergence/nominal "$WORK/$1"; setkey "$WORK/$1" iter-threshold 20; }
both() {  # both <name> <dir> <regex> [solutions-ok]: exit 1 from mpirun, the regex in the log, a STOP line of rank 0 and of rank 1, no solution
  local n=$1 w="$WORK/$2" re=$3 c; c=$(code "$w")
  if [[ $c != 1 ]]; then case_result "$n" 1 "exit code $c, expected 1 ($w/run.log)"
  elif ! grep -qE -- "$re" "$w/run.log"; then case_result "$n" 1 "message '$re' not in run.log"
  elif ! grep -qE '\[RANK +0\] STOP' "$w/run.log" || ! grep -qE '\[RANK +1\] STOP' "$w/run.log"; then case_result "$n" 1 "a rank wrote no STOP line"
  elif [[ ${4:-} != solutions-ok ]] && ls "$w"/OUTPUT/field* > /dev/null 2>&1; then case_result "$n" 1 "a solution was written"
  else case_result "$n" 0 "exit 1, '$re', STOP of both ranks"; fi
}
mpmd() { ( timeout 120 mpirun -n 1 -wdir "$1" "$BIN" : -n 1 -wdir "$2" "$BIN" > "$1/run.log" 2>&1; echo $? > "$1/exit.code" ); }

base control; run "$WORK/control" 2; expect_ok control "$WORK/control"

# 1. an input error: rank 0 checks the registry and broadcasts its verdict
base input; setkey "$WORK/input" vnn -1 FUSS-Numerics; run "$WORK/input" 2; both input input 'vnn must be > 0'

# 2. a bc.txt line that cannot be read: every rank reads the file
base bc; sed -i '1s/302/303/; 2s/.*/ 750.0, 0.0,/' "$WORK/bc/INPUT/bc.txt"; run "$WORK/bc" 2; both bc bc 'bc.txt line 2 \(record 1, type 303\).*short'

# 3. the initial condition of rank 1 cannot be read (MPMD: no ic.tec in its folder): the ranks decide together
base ic-a; cp -a "$WORK/ic-a" "$WORK/ic-b"; rm "$WORK/ic-b/INPUT/ic.tec"; mpmd "$WORK/ic-a" "$WORK/ic-b"
both ic-read ic-a 'reading input file.*\(rank 1\)'

# 4. a cell of the initial condition of rank 1 is not finite
base nan-a; cp -a "$WORK/nan-a" "$WORK/nan-b"
python3 - "$WORK/nan-b/INPUT/ic.tec" <<'PY'
import sys
f = sys.argv[1]; L = open(f).read().split('\n')
# find the first value after the nodal coordinates: the header zone lines give I, J, K of the first zone
import re
m = re.search(r'I=(\d+), J=(\d+), K=(\d+)', L[1]); ni, nj, nk = map(int, m.groups())
L[2 + 3 * ni * nj * nk + 5] = ' NaN'
open(f, 'w').write('\n'.join(L))
PY
mpmd "$WORK/nan-a" "$WORK/nan-b"; both ic-nan nan-a 'holds a value that is not finite.*\(rank 1\)'

# 5. input.ini read again with ini-error = stop and a rule broken: both ranks stop at the reading
base runtime; setkey "$WORK/runtime" iter-threshold 600; setkey "$WORK/runtime" ini-error stop FUSS-IO
setkey "$WORK/runtime" ini-diter 400 FUSS-IO; setkey "$WORK/runtime" sol-diter 100 FUSS-IO; setkey "$WORK/runtime" sol-overwrite false FUSS-IO
setkey "$WORK/runtime" res-threshold 1e-30
( cd "$WORK/runtime" && timeout 120 mpirun -np 2 "$BIN" > run.log 2>&1; echo $? > exit.code ) &
pid=$!
until [[ -f $WORK/runtime/OUTPUT/field1.tec ]] || ! kill -0 $pid 2>/dev/null; do sleep 0.01; done
( cd "$WORK/runtime" && sed 's|^vnn.*|vnn = -1|' input.ini > input.new && mv input.new input.ini )
wait $pid
both runtime-stop runtime 'vnn must be > 0' solutions-ok

# 6. a divergence found inside the OpenMP loop of the time step: the rank that finds it ends with exit code 1 (stop_from_thread,
# no MPI call there) and the launcher ends the other rank, no hang
base diverge; setkey "$WORK/diverge" vnn 0.7 FUSS-Numerics; setkey "$WORK/diverge" iter-threshold 400
( cd "$WORK/diverge" && timeout 120 mpirun -np 2 "$BIN" > run.log 2>&1; echo $? > exit.code )
c=$(code "$WORK/diverge")
if [[ $c == 1 ]] && grep -q 'STOP: NaN or T outside the properties table' "$WORK/diverge/run.log"; then case_result diverge 0 "exit 1, STOP line"; else case_result diverge 1 "exit code $c"; fi

# 7. a block one cell across a connection (mesh-checks/make_case.py thin1s, one block per rank): the rank of block 1 finds
# the connection, the rank of block 2 the sheared thin block; both ranks stop together at the set-up
rm -rf "$WORK/thin"; mkdir -p "$WORK/thin"; python3 -I "$(dirname "${BASH_SOURCE[0]}")/../mesh-checks/make_case.py" thin1s "$WORK/thin" > /dev/null
run "$WORK/thin" 2; both thin thin 'block 2 across it has 1 cell\(s\) across its face 1, fewer than the 2'
grep -q 'block 2 has one cell in direction 1.*not orthogonal' "$WORK/thin/run.log" \
  && case_result thin-lines 0 "message of the rank of block 2" || case_result thin-lines 1 "no 'not orthogonal' message in $WORK/thin/run.log"

# 8. a block number or a cell index of bc.txt outside its range: every rank reads the file, so every rank stops at the same record,
# at the reading and before the table of the owner ranks is indexed with the number (partner block 9 of 2: the MPI library
# ended the run with exit 6; donor block 9: exit 6; a donor cell outside its block: no message, exit 0); chimera deck with its control
chim() { prep numerics-features/chimera "$WORK/$1"; setkey "$WORK/$1" iter-threshold 20; }
chim range-ctl; run "$WORK/range-ctl" 2; expect_ok range-ctl "$WORK/range-ctl"
base range-partner-block; sed -i '22s/.*/       9       1       1       1       1       1       0       0       1/' "$WORK/range-partner-block/INPUT/bc.txt"
run "$WORK/range-partner-block" 2
both range-partner-block range-partner-block 'bc.txt line 22 \(record 11, type 101\): partner block bs = 9 is outside the admissible range 1 to 2'
base range-partner-index; sed -i '22s/.*/       2       1      99       1       1       1       0       0       1/' "$WORK/range-partner-index/INPUT/bc.txt"
run "$WORK/range-partner-index" 2
both range-partner-index range-partner-index 'bc.txt line 22 \(record 11, type 101\): partner cell index js of block 2 = 99 is outside the admissible range 1 to 10'
chim range-donor-block; sed -i '23s/^ *2 /       9 /' "$WORK/range-donor-block/INPUT/bc.txt"
run "$WORK/range-donor-block" 2
both range-donor-block range-donor-block 'bc.txt line 23 \(record 11, type 102\): donor 1, block bD = 9 is outside the admissible range 1 to 2'
chim range-donor-index; sed -i '23s/^ *2 *1 *1 /       2       1      99 /' "$WORK/range-donor-index/INPUT/bc.txt"
run "$WORK/range-donor-index" 2
both range-donor-index range-donor-index 'bc.txt line 23 \(record 11, type 102\): donor 1, cell index jD of block 2 = 99 is outside the admissible range 1 to 7'

finish
