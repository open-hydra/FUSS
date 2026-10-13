#!/usr/bin/env python3
"""Plan 08 section 3.3, Ablation Test Case 2.3: the ablating wall BC 306 on the receding mesh (law nodeshift) against a
local PATO run at maxDeltaT 0.0005 (ref/). The checks, the sampling and the tolerances are those of
test/pyrolysis/tools/atc2x_receding.py.

Usage: verify.py [case_dir] [--record]   Exit 0 on PASS, 1 on FAIL. FUSS_GATE_SELFTEST=1 sets every budget to 0.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import atc2x_receding  # noqa: E402

if __name__ == "__main__":
    sys.exit(atc2x_receding.run(HERE, "2.3", "2.3", sys.argv[1:]))
