#!/bin/bash
#===============================================================================
#         FILE: test/check_build_config.sh
#        USAGE: ./test/check_build_config.sh [--quiet] [--preset NAME]
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
NAME=default
while test $# -gt 0; do
  case $1 in
    --quiet)  QUIET=1; shift ;;
    --preset) NAME=$2; shift 2 ;;
    *) echo "usage: $0 [--quiet] [--preset NAME]" >&2; exit 1 ;;
  esac
done

PRESET="$ROOT/CMakePresets.json"

if [ ! -f "$PRESET" ]; then
  echo "  build-config: no CMakePresets.json, nothing to check against" >&2
  exit 0
fi

# Each preset its own tree (plan 11): the cache checked is <binaryDir>/CMakeCache.txt of
# that preset, its cache variables merged with those of the presets it inherits.
python3 - "$PRESET" "$ROOT" "$NAME" "$QUIET" <<'PY'
import json, os, re, sys

preset_path, root, name, quiet = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == "1"

with open(preset_path) as fh:
    doc = json.load(fh)
presets = {p.get("name"): p for p in doc.get("configurePresets", [])}

def resolved(n, seen=()):
    p = presets.get(n)
    if p is None or n in seen:
        return {}, None
    cv, bdir = {}, None
    inh = p.get("inherits", [])
    for parent in ([inh] if isinstance(inh, str) else inh):
        pcv, pdir = resolved(parent, seen + (n,))
        cv.update(pcv)
        bdir = pdir or bdir
    cv.update(p.get("cacheVariables", {}))
    return cv, p.get("binaryDir", bdir)

if name not in presets:
    print("  build-config: no '%s' preset, nothing to check" % name)
    sys.exit(0)
want, bdir = resolved(name)
bdir = (bdir or "${sourceDir}/build").replace("${sourceDir}", root)
cache_path = os.path.join(bdir, "CMakeCache.txt")
if not os.path.isfile(cache_path):
    print("  build-config: %s is not configured. Run: cmake --preset %s" % (os.path.relpath(bdir, root), name))
    sys.exit(1)

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

# FindMPI moves the directories in which it finds mpi.mod / mpif.h out of
# MPI_<lang>_ADDITIONAL_INCLUDE_DIRS into MPI_<lang>_MODULE_DIR / _F77_HEADER_DIR and
# rewrites the cache entry without them: a preset directory counts as used when it is
# in any of the three (plan 11 stage 0)
def mpi_dirs(k):
    lang = k[len("MPI_"):-len("_ADDITIONAL_INCLUDE_DIRS")]
    got = [have.get(k, "")] + [have.get("MPI_%s_%s" % (lang, x), "") for x in ("MODULE_DIR", "F77_HEADER_DIR")]
    return set(norm(d) for g in got for d in g.split(";") if d and not d.endswith("-NOTFOUND"))

bad = []
for k, v in want.items():
    if isinstance(v, dict):       # {"type": ..., "value": ...} form
        v = v.get("value", "")
    if k.startswith("MPI_") and k.endswith("_ADDITIONAL_INCLUDE_DIRS"):
        missing = set(norm(d) for d in str(v).split(";") if d) - mpi_dirs(k)
        if missing:
            bad.append((k, ";".join(sorted(missing)), "not in the cache's include/module/header dirs"))
    elif k not in have:
        bad.append((k, norm(v), "<absent from cache>"))
    elif norm(have[k]) != norm(v):
        bad.append((k, norm(v), norm(have[k])))

tree = os.path.relpath(bdir, root)
if bad:
    print("  build-config: %s/ DOES NOT match CMakePresets.json (preset %s)" % (tree, name))
    for k, w, h in bad:
        print("      %-28s preset=%-12s cache=%s" % (k, w, h))
    print("  Results from this build are not comparable with any recorded baseline.")
    print("  Fix with:  rm -rf %s && cmake --preset %s && cmake --build %s -j 8" % (tree, name, tree))
    sys.exit(1)

if not quiet:
    print("  build-config: %s/ matches preset %s of CMakePresets.json (%d options checked)" % (tree, name, len(want)))
sys.exit(0)
PY
