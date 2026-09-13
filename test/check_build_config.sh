#!/bin/bash
#===============================================================================
#         FILE: test/check_build_config.sh
#        USAGE: ./test/check_build_config.sh [--quiet]
#  DESCRIPTION: Assert that build/ was configured the way CMakePresets.json says.
#
#               WHY THIS EXISTS
#               ---------------
#               build/ was once rebuilt with a hand-written cmake line instead of
#               the preset. CMake's own defaults are the OPPOSITE of this tree's
#               on two options that change results:
#
#                   USE_OPENMP   preset true    CMake default false
#                   USE_TECIO    preset false   CMake default true
#
#               The resulting binary configured, compiled and ran without a
#               single complaint, and differed from the recorded baseline in
#               exactly ONE artefact -- convective_Plate's wall.tec, with
#               field.tec identical. It cost two full regression runs to find
#               and explain.
#
#               Nothing about that failure was loud. A serial binary is not an
#               error, TecIO is not an error, and 289 of 290 artefacts agreeing
#               looks like success. So the check is mechanical and runs before
#               the suite, rather than depending on anyone remembering.
#
#               Quick tell if you are ever debugging this by hand: the correct
#               binary is ~1.09 MB; a TecIO-linked one is ~3.0 MB.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
QUIET=0
[ "${1:-}" = "--quiet" ] && QUIET=1

CACHE="$ROOT/build/CMakeCache.txt"
PRESET="$ROOT/CMakePresets.json"

if [ ! -f "$PRESET" ]; then
  echo "  build-config: no CMakePresets.json, nothing to check against" >&2
  exit 0
fi
if [ ! -f "$CACHE" ]; then
  echo "  build-config: build/ is not configured. Run: cmake --preset default" >&2
  exit 1
fi

python3 - "$PRESET" "$CACHE" "$QUIET" <<'PY'
import json, re, sys

preset_path, cache_path, quiet = sys.argv[1], sys.argv[2], sys.argv[3] == "1"

with open(preset_path) as fh:
    doc = json.load(fh)

want = {}
for p in doc.get("configurePresets", []):
    if p.get("name") == "default":
        want = dict(p.get("cacheVariables", {}))
        break
if not want:
    print("  build-config: no 'default' preset, nothing to check")
    sys.exit(0)

# Read the cache as NAME:TYPE=VALUE
have = {}
with open(cache_path) as fh:
    for line in fh:
        m = re.match(r"^([A-Za-z_0-9]+):[A-Z]+=(.*)$", line.strip())
        if m:
            have[m.group(1)] = m.group(2)

def norm(v):
    v = str(v).strip()
    low = v.lower()
    if low in ("on", "true", "1", "yes"):   return "true"
    if low in ("off", "false", "0", "no"):  return "false"
    return v.rstrip("/")          # paths in the preset may carry a trailing /

bad = []
for k, v in want.items():
    if isinstance(v, dict):       # {"type": ..., "value": ...} form
        v = v.get("value", "")
    if k not in have:
        bad.append((k, norm(v), "<absent from cache>"))
    elif norm(have[k]) != norm(v):
        bad.append((k, norm(v), norm(have[k])))

if bad:
    print("  build-config: build/ DOES NOT match CMakePresets.json")
    for k, w, h in bad:
        print("      %-28s preset=%-12s cache=%s" % (k, w, h))
    print("  Results from this build are not comparable with any recorded baseline.")
    print("  Fix with:  rm -rf build && cmake --preset default && cmake --build build -j 8")
    sys.exit(1)

if not quiet:
    print("  build-config: matches CMakePresets.json (%d options checked)" % len(want))
sys.exit(0)
PY
