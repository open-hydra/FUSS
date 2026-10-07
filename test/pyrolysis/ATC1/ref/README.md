# Reference data for the ATC1 gate (plan 06 section 4.4)

- `FIAT_T.dat`, `FIAT_pyrolysisFront.dat`: copies of PATO 3.1's `tutorials/1D/AblationTestCase_1.0/data/ref/FIAT/{T,pyrolysisFront}`
  (the FIAT solution of Ablation Test Case 1; `T` columns `t, T(0 mm), T(1), T(2), T(4), T(8), T(16), T(50 mm)`, 0..60 s).
- `pato-298K/`: PATO 3.1 run of the unmodified tutorial with the initial and the first ramp value at 298 K
  (`origin.0/porousMat/Ta`: `internalField uniform 298`, table entry `(0 298)`), 2026-10-07; `Ta_plot` columns
  `t` then 1, 2, 4, 8, 12, 16, 24 mm below the heated face, `Ta_surfacePatch` columns `t`, surface, back face;
  `spread_298.txt` = `tools/measure_spread.py` on it (13.18 K max over the six gated probes).
