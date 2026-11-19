#!/bin/bash
# DecompositionMap: the key decomposition-map of [FUSS-IO] names a file of INPUT/ (in the format of the map of the ATLAS
# mesh decomposer: one line per block with 13 integers: new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank) that gives the MPI
# rank of every block. One rank here: a map with every block on rank 0 is the run without a map, byte for byte; every error
# of the map stops the run (exit code 1, the message names the file and the line); a map in INPUT/ is read only when the key
# names it.
#   usage: decomposition-map.sh <FUSS binary> <work folder>
BIN=$(readlink -f "$1"); WORK=$(readlink -m "$2"); source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
R=numerics-features/residual_convergence/nominal

base() { prep $R "$WORK/$1"; setkey "$WORK/$1" iter-threshold 20; }
# mapfile <dir> <rank of block 1> <rank of block 2>: INPUT/decomposition.map of the two blocks of nominal
mapfile() { python3 "$(dirname "${BASH_SOURCE[0]}")/write_map.py" "$1/INPUT/ic.tec" "$1/INPUT/decomposition.map" "$2" "$3"; }
withmap() { setkey "$WORK/$1" decomposition-map "${2:-decomposition.map}" FUSS-IO; }

base plain; run "$WORK/plain"; expect_ok plain "$WORK/plain"
base allzero; mapfile "$WORK/allzero" 0 0; withmap allzero; run "$WORK/allzero"
n=$(ndiff "$WORK/plain" "$WORK/allzero"); if [[ $(code "$WORK/allzero") == 0 && $n == 0 ]] && grep -q 'partition read from' "$WORK/allzero/run.log"; then case_result map-all-rank-0 0 "same files as the run without a map"; else case_result map-all-rank-0 1 "exit $(code "$WORK/allzero"), $n files differ"; fi

e() {  # e <name> <regex> <python edit of the map lines L>
  base "$1"; mapfile "$WORK/$1" 0 0; withmap "$1"
  python3 - "$WORK/$1/INPUT/decomposition.map" "$3" <<'PY'
import sys
f, code = sys.argv[1:3]; L = [l for l in open(f).read().split('\n') if l.strip() != '' and l[0] != '#']
ns = {'L': L}; exec(code, ns); open(f, 'w').write('# edited\n' + '\n'.join(ns['L']) + '\n# end of block lines\n')
PY
  run "$WORK/$1"; expect_stop "$1" "$WORK/$1" "$2"
}
e map-one-line  'decomposition.map: 1 block lines for the 2 blocks' 'L = L[:1]'
e map-order     'lists block 2 where block 1 is expected'           'L = [L[1], L[0]]'
e map-size      'block 1 is 7 x'                                    'w = L[0].split(); w[8] = "7"; L[0] = " ".join(w)'
e map-cells     'cells in the map'                                  'w = L[0].split(); w[11] = "5"; L[0] = " ".join(w)'
e map-rank      'is on rank 3, outside the ranks of the run'         'w = L[1].split(); w[12] = "3"; L[1] = " ".join(w)'
e map-not-13    'is not 13 integers'                                 'L[1] = "2 0 1 2"'
e map-words     'is not 13 integers'                                 'L[1] = L[1].replace(L[1].split()[1], "x", 1)'

base nofile; withmap nofile nothere.map; run "$WORK/nofile"; expect_stop map-missing "$WORK/nofile" 'nothere.map: cannot be opened'
base isdir; mkdir "$WORK/isdir/INPUT/adir"; withmap isdir adir; run "$WORK/isdir"; expect_stop map-directory "$WORK/isdir" 'is a directory, not a file'
base folder; mapfile "$WORK/folder" 0 0; withmap folder "INPUT/decomposition.map"; run "$WORK/folder"; expect_stop map-folder "$WORK/folder" 'without a folder'

# a map in INPUT/ is read only when the key names it: a wrong map without the key changes nothing
base unnamed; mapfile "$WORK/unnamed" 0 7; run "$WORK/unnamed"
n=$(ndiff "$WORK/plain" "$WORK/unnamed"); [[ $(code "$WORK/unnamed") == 0 && $n == 0 ]] && case_result map-unnamed 0 "not read" || case_result map-unnamed 1 "exit $(code "$WORK/unnamed"), $n files differ"
finish
