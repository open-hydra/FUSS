#!/bin/bash
# ThinBlocks: a block one cell across stops the set-up with exit code 1 and a message that names the block (before the ghost
# metrics read node planes outside a block), and the thin blocks that are exact still run. Decks from make_case.py, plane
# mesh, 8 cells along the long side, 1 cell across (0.0625 thick), one Euler step.
#   thin1, thin1s, thin1j  a 1 x 8 (8 x 1) block joined by a connection (101) to an 8 x 8 block: the block across the
#                          connection has 1 cell across its face, fewer than the gc = 2 of its ghost cells (thin1s: sheared,
#                          so the thin block also fails the second rule)
#   thin1ws, thin1jws      a sheared 1 x 8 and a sheared 8 x 1 block, walls (302) only: the lines of the thin direction are
#                          not orthogonal to the others, the face gradient takes no difference along it (qw error 0.16-0.30)
#   thin1w, thin1jw        the same blocks unsheared: they run (exact), as do the k direction of every plane case and
#                          thin3d (MetricChecks: one cell in k of a 3-D mesh, Cartesian)
#   usage: thin-blocks.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/../input-io/lib.sh"
GEN="$(dirname "${BASH_SOURCE[0]}")/make_case.py"

deck() { rm -rf "$WORK/$1"; mkdir -p "$WORK/$1"; python3 -I "$GEN" "$1" "$WORK/$1" > /dev/null || { case_result "$1" 1 "no deck"; return 1; }; run "$WORK/$1"; }

# stops: exit 1, the message of the rule, a STOP line, no solution
for p in thin1 thin1j; do
  deck $p && expect_stop $p "$WORK/$p" 'block 2 across it has 1 cell\(s\) across its face [13], fewer than the 2 that the ghost cells'
done
deck thin1s  && expect_stop thin1s  "$WORK/thin1s"  'block 2 across it has 1 cell\(s\) across its face 1, fewer than the 2 that the ghost cells'
grep -q 'block 2 has one cell in direction 1 and its lines in direction 2 are not orthogonal' "$WORK/thin1s/run.log" \
  && case_result thin1s-lines 0 "second rule too" || case_result thin1s-lines 1 "no 'not orthogonal' message in $WORK/thin1s/run.log"
deck thin1ws  && expect_stop thin1ws  "$WORK/thin1ws"  'block 1 has one cell in direction 1 and its lines in direction 2 are not orthogonal'
deck thin1jws && expect_stop thin1jws "$WORK/thin1jws" 'block 1 has one cell in direction 2 and its lines in direction 1 are not orthogonal'

# controls: the unsheared thin blocks run
for p in thin1w thin1jw; do deck $p && expect_ok $p "$WORK/$p"; done

finish
