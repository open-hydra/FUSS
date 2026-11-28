#!/bin/bash
#===============================================================================
#         FILE: test/moving-mesh/receding/run.sh
#  DESCRIPTION: Run the whole grid sequence for gate 5.3 and report the order.
#               The two finest grids are the slow part; N320 is ~20 minutes.
#===============================================================================
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

for d in "$HERE"/N*/; do
  n=$(basename "$d")
  printf '%-8s ' "$n"
  ( cd "$d" && rm -rf OUTPUT && mkdir -p OUTPUT bin && ./FUSS.sh -p 1 solve ) > "$d/run.log" 2>&1 \
    && echo ok || { echo "FAILED (see $d/run.log)"; exit 1; }
done

python3 "$HERE/verify.py"
