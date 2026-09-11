#!/bin/bash
#===============================================================================
#         FILE: test/run_guard_demos.sh
#        USAGE: ./test/run_guard_demos.sh [-w WORKDIR]
#  DESCRIPTION: Every moving-mesh guard is run against the configuration it is
#               supposed to reject, and must actually reject it.
#
#               A guard that has never been seen to fire is a comment. Two of
#               these were written and did NOT fire the first time -- the
#               IRS/multigrid checks read obj_mesh_motion%enabled, which
#               Check_Input runs too early to see -- and that was only caught by
#               trying them. Hence this file.
#
#               Matching is on the ERROR TEXT, not the exit status: Fortran's
#               bare `stop` exits 0, so a guard that fires and a run that
#               succeeds are indistinguishable by exit code alone. That is worth
#               knowing before trusting any "it failed, good" check here.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BASE="$ROOT/test/moving-mesh/receding/N040"
WORK="${TMPDIR:-/tmp}/fuss-guard-demos.$$"

while test $# -gt 0; do
  case $1 in
    -w | --work) WORK=$2; shift 2 ;;
    -h | --help) echo "usage: $0 [-w WORKDIR]"; exit 1 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$WORK"
fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "FUSS moving-mesh guard demonstrations"
echo "  base case : test/moving-mesh/receding/N040"
echo "  work dir  : $WORK"
echo

# $1 = short name, $2 = expected substring of the error, $3.. = ini edits
# (each edit is "section-less key = value" applied by regex, or "+SECTION:line")
demo () {
  local name=$1 expect=$2; shift 2
  local d="$WORK/$name"

  rm -rf "$d"; mkdir -p "$d"
  cp -r "$BASE/INPUT" "$BASE/MESH" "$BASE/input.ini" "$BASE/FUSS.sh" "$d/"
  sed -i "s|^MASTERDIR=.*|MASTERDIR=$ROOT|" "$d/FUSS.sh"
  mkdir -p "$d/OUTPUT" "$d/bin"

  # Keep every demo to a couple of steps: the guards all fire at setup or on
  # the first mesh update, so a long run would only waste time.
  python3 - "$d/input.ini" "iter-threshold" "3" <<'PY'
import re, sys
p, k, v = sys.argv[1:4]
s = open(p).read()
pat = re.compile(r"^(%s\s*=).*$" % re.escape(k), re.M)
s = pat.sub(r"\1 " + v, s) if pat.search(s) else s
open(p, "w").write(s)
PY

  local kv
  for kv in "$@"; do
    if [ "${kv:0:4}" = "RAW:" ]; then
      printf '\n%s\n' "${kv:4}" >> "$d/input.ini"
      continue
    fi
    python3 - "$d/input.ini" "${kv%%=*}" "${kv#*=}" <<'PY'
import re, sys
p, k, v = sys.argv[1], sys.argv[2].strip(), sys.argv[3].strip()
s = open(p).read()
pat = re.compile(r"^(%s\s*=).*$" % re.escape(k), re.M)
if pat.search(s):
    s = pat.sub(lambda m: m.group(1) + " " + v, s)
else:
    # New keys belong to the motion section. Dropping them into an unrelated
    # section is not an error the registry reports -- it simply never reads
    # them -- so the guard would appear not to fire when in fact it was never
    # given the input that should trigger it. That happened while writing this.
    s = s.replace("[FUSS-MeshMotion]", "[FUSS-MeshMotion]\n%s = %s" % (k, v), 1)
open(p, "w").write(s)
PY
  done

  ( cd "$d" && ./FUSS.sh -p 1 solve ) > "$d/run.log" 2>&1
  if grep -qF "$expect" "$d/run.log"; then
    pass "$name -> refused with: $(grep -m1 -F "$expect" "$d/run.log" | sed 's/^ *//')"
  else
    bad "$name: expected \"$expect\", not found in $d/run.log"
    tail -4 "$d/run.log" | sed 's/^/          /'
  fi
}

# ---------------------------------------------------------------------------
# Input-validation guards (fire before the first step).
# ---------------------------------------------------------------------------
demo irs        "not compatible with implicit residual smoothing" "irs=true"
demo multigrid  "not compatible with multigrid" \
    "RAW:[FUSS-Multigrid]
levels = 2
level1-iter = 3
level2-iter = 3"
demo primitive  "requires integration-variables = cons"           "integration-variables=prim"
demo steady     "requires a time-accurate run"                    "time-accurate=false"
demo restart    "cannot be restarted"                             "law=prescribed" "newrun=false"

# ---------------------------------------------------------------------------
# Runtime geometry guards (fire on the first mesh update).
#
# too-fast: rigid translation never tangles a mesh -- volumes are preserved --
#           so this isolates the remap's stability limit on its own. dx/dt for
#           this case is about 0.029 m/s, so 1 m/s oversteps it by ~35x.
# tangled:  an oscillation whose amplitude dwarfs the cell size folds cells
#           inside out, which is a different failure and a different message.
# ---------------------------------------------------------------------------
demo too_fast "mesh moved too far in one time step" "vel=1.0 0.0 0.0"
demo tangled  "mesh update failed" \
              "law=prescribed" "amp=0.05 0.05 0.0" "kx=500.0 500.0 0.0" \
              "ky=500.0 500.0 0.0" "omega=1.0"

echo
if [ "$fail" -eq 0 ]; then
  echo "all guards fired as intended"
else
  echo "$fail guard(s) did NOT fire -- they are decoration, not protection"
  exit 1
fi
