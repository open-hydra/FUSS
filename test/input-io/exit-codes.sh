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
# a code that FUSS does not know is read without a data line: its data line is read as the next record, and the stop names it
plate bc-code;    sed -i '1s/302$/399/' "$WORK/bc-code/INPUT/bc.txt"; go bc-code 'bc.txt line 2 \(record 2\).*record 1 has code 399'

# bc.txt: a block number, a face or a cell index outside the block, in the record itself (nominal, 2 blocks of 100 x 10 x 1 cells:
# line 21 is the header of record 11, a connection of block 1 face 2, line 22 its partner), in the partner of a connection (101),
# in the donor of a chimera record (102, chimera deck: line 21 the header of record 11, line 23 its first donor). The stop is at the
# reading, with the line, the record, the value and the range: before, a block that does not exist ended the run in the MPI
# library (exit 6) or with a segmentation fault, and a donor cell outside its block was read past the array without a message.
# The donors of a coupled wall (104) are cells of the other solver: not tested.
chim() { prep numerics-features/chimera "$WORK/$1"; setkey "$WORK/$1" iter-threshold 5; }
nominal bc-range-ctl; run "$WORK/bc-range-ctl"; expect_ok bc-range-ctl "$WORK/bc-range-ctl"
chim bc-range-chim-ctl; run "$WORK/bc-range-chim-ctl"; expect_ok bc-range-chim-ctl "$WORK/bc-range-chim-ctl"
nominal bc-own-block;  sed -i '21s/.*/       7     100       1       1       2     101/' "$WORK/bc-own-block/INPUT/bc.txt"
go bc-own-block 'bc.txt line 21 \(record 11, type 101\): block b = 7 is outside the admissible range 1 to 2'
nominal bc-own-index;  sed -i '21s/.*/       1     100      11       1       2     101/' "$WORK/bc-own-index/INPUT/bc.txt"
go bc-own-index 'bc.txt line 21 \(record 11, type 101\): cell index j of block 1 = 11 is outside the admissible range 1 to 10'
nominal bc-own-face;   sed -i '21s/.*/       1     100       1       1       7     101/' "$WORK/bc-own-face/INPUT/bc.txt"
go bc-own-face 'bc.txt line 21 \(record 11, type 101\): face f of block 1 = 7 is outside the admissible range 1 to 6'
nominal bc-partner-block; sed -i '22s/.*/       9       1       1       1       1       1       0       0       1/' "$WORK/bc-partner-block/INPUT/bc.txt"
go bc-partner-block 'bc.txt line 22 \(record 11, type 101\): partner block bs = 9 is outside the admissible range 1 to 2'
nominal bc-partner-index; sed -i '22s/.*/       2       1      99       1       1       1       0       0       1/' "$WORK/bc-partner-index/INPUT/bc.txt"
go bc-partner-index 'bc.txt line 22 \(record 11, type 101\): partner cell index js of block 2 = 99 is outside the admissible range 1 to 10'
nominal bc-partner-face; sed -i '22s/.*/       2       1       1       1       0       1       0       0       1/' "$WORK/bc-partner-face/INPUT/bc.txt"
go bc-partner-face 'bc.txt line 22 \(record 11, type 101\): partner face fs of block 2 = 0 is outside the admissible range 1 to 6'
chim bc-donor-block;  sed -i '23s/^ *2 /       9 /' "$WORK/bc-donor-block/INPUT/bc.txt"
go bc-donor-block 'bc.txt line 23 \(record 11, type 102\): donor 1, block bD = 9 is outside the admissible range 1 to 2'
chim bc-donor-index;  sed -i '23s/^ *2 *1 *1 /       2       1      99 /' "$WORK/bc-donor-index/INPUT/bc.txt"
go bc-donor-index 'bc.txt line 23 \(record 11, type 102\): donor 1, cell index jD of block 2 = 99 is outside the admissible range 1 to 7'
chim bc-donor-late;   sed -i '24s/^ *2 *2 *1 /       2     101       1 /' "$WORK/bc-donor-late/INPUT/bc.txt"
go bc-donor-late 'bc.txt line 24 \(record 11, type 102\): donor 2, cell index iD of block 2 = 101 is outside the admissible range 1 to 100'

# time tables of a wall temperature (302): decreasing times, a missing file, an empty file, a row that is not a number
tt() {  # tt <name> <file content or -> <regex>
  plate "$1"; sed -i '2s/.*/ INPUT\/Twall.txt,/' "$WORK/$1/INPUT/bc.txt"
  [[ $2 != - ]] && printf "$2" > "$WORK/$1/INPUT/Twall.txt"
  go "$1" "$3"
}
tt table-decreasing '0.0 300\n1.0 310\n0.5 320\n' 'time table .*Twall.txt.*times must increase'
tt table-missing    -                                  'time table .*Twall.txt.*not found'
tt table-empty      ''                                 'time table .*Twall.txt.*empty'
tt table-bad-row    '0.0 300\n1.0 abc\n'               'time table .*Twall.txt.*cannot read row 2'
tt table-bad-line   '0.0 300\n\n1.0 abc\n'             'time table .*Twall.txt.*cannot read row 2 \(line 3\)'
# blank lines, also at the end, are not rows (as MOSE)
plate table-blank; sed -i '2s/.*/ INPUT\/Twall.txt,/' "$WORK/table-blank/INPUT/bc.txt"
printf '0.0 300\n\n1.0 310\n  \n\n' > "$WORK/table-blank/INPUT/Twall.txt"; run "$WORK/table-blank"; expect_ok table-blank "$WORK/table-blank"

# the probes: a block that does not exist; an index outside the block; an unknown variable
rod probe-block; sed -i -E 's|^index-position = 1 90 1 1|index-position = 5 90 1 1|' "$WORK/probe-block/input.ini"; go probe-block 'probe block 5 does not exist'
rod probe-index; sed -i -E 's|^index-position = 1 90 1 1|index-position = 1 900 1 1|' "$WORK/probe-index/input.ini"; go probe-index 'probe index-position 900 outside block 1'
rod probe-var;   sed -i -E '0,/^variables = T/s//variables = T Z/' "$WORK/probe-var/input.ini"; go probe-var 'unknown probe variable "Z"'

# FAS cycle with one level; multigrid levels that the blocks cannot divide or leave thinner than 2 cells; a divergence
nominal fas-levels; printf '\n[FUSS-Multigrid]\ncycle = w\n' >> "$WORK/fas-levels/input.ini"; go fas-levels 'cycle needs levels >= 2'

nominal mg-divisible; printf '\n[FUSS-Multigrid]\nlevels = 4\nlevel1-iter = 5\nlevel2-iter = 5\nlevel3-iter = 5\nlevel4-iter = 5\n' >> "$WORK/mg-divisible/input.ini"; go mg-divisible 'multigrid with 4 levels'
nominal diverge; setkey "$WORK/diverge" vnn 0.7 FUSS-Numerics; setkey "$WORK/diverge" iter-threshold 400; go diverge 'NaN or T outside the properties table'

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

# an initial condition with T only (no matID): matID was read out of the bounds of the variables
plate ic-1var
python3 - "$WORK/ic-1var/INPUT/ic.tec" <<'PY'
import sys
f = sys.argv[1]; L = open(f).read().split('\n')
L[0] = L[0].replace(' matID', ''); L[1] = L[1].replace('[4-5]=CELLCENTERED', '[4]=CELLCENTERED')
n = 2 + 3 * 61 * 41 * 2 + 60 * 40          # header, nodes, temperatures: the matID values follow
open(f, 'w').write('\n'.join(L[:n]) + '\n')
PY
go ic-1var 'ic.tec: block 1 holds 1 variable\(s\); FUSS needs T and matID'

# properties.dat: a value that is not finite and a density of zero, in the row of 1 K that no cell of the plate reaches
plate prop-nan;  sed -i -E '0,/^1\.0 [0-9.]+ /s//1.0 NaN /' "$WORK/prop-nan/INPUT/properties.dat"
go prop-nan 'properties.dat, material mat, row T = 1 K: Cp is not finite'
plate prop-zero; sed -i -E '0,/^1\.0 ([0-9.]+) [0-9.]+ /s//1.0 \1 0.0 /' "$WORK/prop-zero/INPUT/properties.dat"
go prop-zero 'properties.dat, material mat, row T = 1 K: Density = +0.000E\+00 must be > 0'

# the stop of a rank that is not rank 0 is tested with MPI (mpi-exit-codes.sh)
finish
