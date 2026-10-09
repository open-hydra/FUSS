#!/usr/bin/env python3
"""Plan 07 gate 4.1(d): the adiabatic0d case under gas-closure = type2 (see ../../adiabatic0d/verify.py, which
detects the closure from input.ini and compares against the type2 variant of tools/ode0d.py)."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
case = sys.argv[1] if len(sys.argv) > 1 else HERE
sys.exit(subprocess.call([sys.executable, os.path.join(HERE, "..", "..", "adiabatic0d", "verify.py"), case]))
