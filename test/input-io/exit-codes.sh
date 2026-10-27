#!/bin/bash
# ExitCodes: every error that stops the run ends it with exit code 1, a message that names the file, the line or the key,
# and no solution written. Base cases: temperature_Plate (60 x 40, steady, wall temperatures), oscillating_Rod (probes),
# residual_convergence/nominal (two blocks).
#   usage: exit-codes.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

plate() { prep steady-state/temperature_Plate "$WORK/$1"; setkey "$WORK/$1" iter-threshold 5; }
rod()   { prep transient/oscillating_Rod "$WORK/$1"; setkey "$WORK/$1" iter-threshold 5; }
nominal() { prep numerics-features/residual_convergence/nominal "$WORK/$1"; setkey "$WORK/$1" iter-threshold 5; }
go() { run "$WORK/$1"; expect_stop "$1" "$WORK/$1" "$2"; }

# control
plate control; run "$WORK/control"; expect_ok control "$WORK/control"

# bc.txt: a record line that is not six integers; a property line that is short; a property line with a null value; a file with fewer records
plate bc-record;  sed -i '3s/.*/ 1 x 2 1 1 302/' "$WORK/bc-record/INPUT/bc.txt"; go bc-record 'bc.txt.*(line 3|record|Boundary conditions number)'
plate bc-short;   sed -i '1s/302/303/; 2s/.*/ 750.0, 0.0,/' "$WORK/bc-short/INPUT/bc.txt"; go bc-short 'bc.txt line 2 \(record 1, type 303\).*short'
plate bc-null;    sed -i '1s/302/303/; 2s/.*/ 750.0,, 273.15,/' "$WORK/bc-null/INPUT/bc.txt"; go bc-null 'bc.txt line 2 \(record 1, type 303\).*qrad is empty'
plate bc-text;    sed -i '1s/302/303/; 2s/.*/ 750.0, 0.0, abc,/' "$WORK/bc-text/INPUT/bc.txt"; go bc-text 'bc.txt line 2 \(record 1, type 303\)'
plate bc-count;  head -n -2 "$WORK/bc-count/INPUT/bc.txt" > "$WORK/bc-count/INPUT/bc.new"; mv "$WORK/bc-count/INPUT/bc.new" "$WORK/bc-count/INPUT/bc.txt"; go bc-count 'Boundary conditions number'
plate bc-missing; rm "$WORK/bc-missing/INPUT/bc.txt"; go bc-missing 'Boundary condition file not found'

# time tables of a wall temperature (302): decreasing times, a missing file, an empty file
tt() {  # tt <name> <file content or -> <regex>
  plate "$1"; sed -i '2s/.*/ INPUT\/Twall.txt,/' "$WORK/$1/INPUT/bc.txt"
  [[ $2 != - ]] && printf "$2" > "$WORK/$1/INPUT/Twall.txt"
  go "$1" "$3"
}
tt table-decreasing '0.0 300\n1.0 310\n0.5 320\n' 'time table .*Twall.txt.*times must increase'
tt table-missing    -                                  'time table .*Twall.txt.*not found'
tt table-empty      ''                                 'time table .*Twall.txt.*empty'
tt table-bad-row    '0.0 300\n1.0 abc\n'               'time table .*Twall.txt.*cannot read row 2'

# the probes: a block that does not exist; an index outside the block; an unknown variable
rod probe-block; sed -i -E 's|^index-position = 1 90 1 1|index-position = 5 90 1 1|' "$WORK/probe-block/input.ini"; go probe-block 'probe block 5 does not exist'
rod probe-index; sed -i -E 's|^index-position = 1 90 1 1|index-position = 1 900 1 1|' "$WORK/probe-index/input.ini"; go probe-index 'probe index-position 900 outside block 1'
rod probe-var;   sed -i -E '0,/^variables = T/s//variables = T Z/' "$WORK/probe-var/input.ini"; go probe-var 'unknown probe variable "Z"'

# FAS cycle with one level; multigrid levels that the blocks cannot divide are checked at commit 9 (Lib_Multigrid)
nominal fas-levels; printf '\n[FUSS-Multigrid]\ncycle = w\n' >> "$WORK/fas-levels/input.ini"; go fas-levels 'cycle needs levels >= 2'

# the initial condition: a file that is missing; a restart without the solution
plate ic-missing; rm "$WORK/ic-missing/INPUT/ic.tec"; go ic-missing 'reading input file'
plate restart-missing; setkey "$WORK/restart-missing" newrun false; go restart-missing 'reading input file'

# values that are not finite: a wall temperature, a time table, a cell of the initial condition
plate bc-nan;  sed -i '2s/.*/ NaN,/' "$WORK/bc-nan/INPUT/bc.txt"; go bc-nan 'bc.txt line 2 \(record 1, type 302\).*T_wall is not finite'
tt table-nan '0.0 300\n1.0 NaN\n' 'time table "INPUT/Twall.txt", row 2.*not finite'
plate ic-nan
python3 - "$WORK/ic-nan/INPUT/ic.tec" <<'PY'
import sys
f = sys.argv[1]; L = open(f).read().split('\n')
# header (2 lines), then the coordinates of the 61 x 41 x 2 nodes, then the 60 x 40 temperatures of the cells
k = 2 + 3 * 61 * 41 * 2 + 5
L[k] = ' NaN'
open(f, 'w').write('\n'.join(L))
PY
go ic-nan 'ic.tec holds a value that is not finite.*T of block 1, cell i, j, k = 6, 1, 1'

# the stop of a rank that is not rank 0 is tested with MPI (mpi-exit-codes.sh)
finish
