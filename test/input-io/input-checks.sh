#!/bin/bash
# InputChecks: the values of input.ini are checked, and every error stops the run at set-up with exit code 1, the message
# naming the key and its value, and no solution written. Base case: oscillating_Rod (two probes), 5 iterations.
#   usage: input-checks.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

base() {  # base <name>: a fresh copy of the base case in $WORK/<name>
  prep transient/oscillating_Rod "$WORK/$1"; setkey "$WORK/$1" iter-threshold 5
}

# control: the unchanged case runs
base control; run "$WORK/control"; expect_ok control "$WORK/control"

# values out of their rule or that do not read as their type
t() {  # t <name> <section> <key> <value> <message regex>
  base "$1"; setkey "$WORK/$1" "$3" "$4" "$2"; run "$WORK/$1"; expect_stop "$1" "$WORK/$1" "$5"
}
t vnn-negative      FUSS-Numerics vnn -1          'vnn.*> 0'
t vnn-nan           FUSS-Numerics vnn NaN         'vnn'
t vnn-word          FUSS-Numerics vnn abc         'vnn.*abc'
t sol-diter-real    FUSS-IO       sol-diter 2.5   'sol-diter.*2.5.*not an integer'
t sol-diter-exp     FUSS-IO       sol-diter 1e4   'sol-diter.*1e4'
t sol-overwrite-yes FUSS-IO       sol-overwrite yes 'sol-overwrite.*yes'
t ini-error-word    FUSS-IO       ini-error maybe 'ini-error'
t time-scheme-bad   FUSS-Numerics time-scheme RK9 'time-scheme'
t res-threshold-0   FUSS-Parameters res-threshold 0 'res-threshold.*> 0'

# levels is read outside the registry: its value is checked too
base levels-real; printf '\n[FUSS-Multigrid]\nlevels = 2.0\n' >> "$WORK/levels-real/input.ini"; run "$WORK/levels-real"
expect_stop levels-real "$WORK/levels-real" 'levels = 2.0.*not an integer'

# the probes
base probe-index; sed -i -E 's|^index-position = 1 90 1 1|index-position = 1 -1 1 1|' "$WORK/probe-index/input.ini"
run "$WORK/probe-index"; expect_stop probe-index "$WORK/probe-index" 'index-position'
base probe-gap; sed -i -E 's|^probe2 = p2|probe3 = p2|' "$WORK/probe-gap/input.ini"
run "$WORK/probe-gap"; expect_stop probe-gap "$WORK/probe-gap" 'probe3.*probe2 is missing.*without gaps'

# input.ini that cannot be used: missing, empty, comment lines only
base missing; rm "$WORK/missing/input.ini"; run "$WORK/missing"; expect_stop missing "$WORK/missing" 'input.ini.*cannot be opened'
base empty; : > "$WORK/empty/input.ini"; run "$WORK/empty"; expect_stop empty "$WORK/empty" 'input.ini.*is empty'
base comments; printf '; nothing here\n; a = b\n' > "$WORK/comments/input.ini"; run "$WORK/comments"
expect_stop comments "$WORK/comments" 'input.ini'

# a comment line that holds "=" is a comment: the run is the one without it (byte for byte), and nothing is reported as not read
base eqline-a; run "$WORK/eqline-a"
base eqline-b; sed -i -E 's|^\[FUSS-Numerics\]|[FUSS-Numerics]\n; a comment with a = sign|' "$WORK/eqline-b/input.ini"; run "$WORK/eqline-b"
n=$(ndiff "$WORK/eqline-a" "$WORK/eqline-b"); if [[ $(code "$WORK/eqline-b") == 0 && $n == 0 ]] && ! grep -q 'is not read' "$WORK/eqline-b/run.log"; then
  case_result eqline 0 "same files, no warning"; else case_result eqline 1 "exit $(code "$WORK/eqline-b"), $n files differ, or a warning"; fi

# a line "key =" that the reader does not give is reported with its line number, and the run is unchanged
base unread; sed -i -E 's|^\[FUSS-Numerics\]|[FUSS-Numerics]\nfoo =|' "$WORK/unread/input.ini"; run "$WORK/unread"
n=$(ndiff "$WORK/eqline-a" "$WORK/unread")
if [[ $(code "$WORK/unread") == 0 && $n == 0 ]] && grep -q 'input.ini, line [0-9]*: foo of \[FUSS-Numerics\] is not read' "$WORK/unread/run.log"; then
  case_result unread 0 "warning with the line, same files"; else case_result unread 1 "exit $(code "$WORK/unread"), $n files differ, or no warning"; fi

finish
