# Adiabatic 0-D Decomposition

The physics, the reference solution and the results are on the V&V page:
[Pyrolysis Kinetics and 0-D Decomposition](../../../docs/vv/pyrolysis-kinetics-0d.md).

## Setup (`T1000`)

| Parameter | Value |
|---|---|
| Mesh | one block, 4 cells over 4 mm |
| Both i-faces | adiabatic (`301`, `q = 0`); faces 3/4 also `301` |
| Material | real TACOT tables, unchanged |
| Initial condition | 1000 K |
| Pyrolysis | `[FUSS-Pyrolysis] gas-closure = type1`, `outflow-direction = +i`, `p-gas = 101325.0` |
| Time scheme | RK3, `vnn = 0.25`, time-accurate |
| Integration variables | `cons` |
| Stop | `time-threshold = 20.0`, a field file every 0.5 s |

Nothing varies along `i`, so every cell stays uniform and follows the zero-dimensional decomposition
ODE of `test/pyrolysis/tools/ode0d.py` (type-1 variant) exactly, except for the kinetics/energy
splitting inside the RK3 stages (frozen stage temperature), which is second order in `dt`. Generated
by `test/pyrolysis/tools/make_case.py --T0 1000` (see the `input.ini` header).

The type-2 twin is `test/pyrolysis/darcy/adiabatic0d_type2/` (both `301` faces open at `p-wall`,
`K = 1e-7` m², so the pressure excursion stays below 0.1 Pa and the fixed-`p` oracle applies); its
`verify.py` calls this directory's `verify.py`, which detects the closure.

## Check

Over every `OUTPUT/field<i>.tec`:

- `max |T - T_or1| <= TOL_T` (0.10 K; measured 0.051 K at `vnn = 0.5`, now run at `0.25`);
- `max |chi_r - chi_or1| <= TOL_CHI` (2e-3);
- spatial non-uniformity of `T` and `chi` `<= TOL_UNIF` (1e-6).

Reported, not gated: the exact type-2 (gas storage) minus type-1 difference, and the FUSS error at
the last file.

## Run

```bash
cd test/pyrolysis/adiabatic0d/T1000 && ./FUSS.sh solve
cd .. && python3 verify.py T1000
```

Physics gate: `test/run_physics_gates.sh -c pyrolysis` (together with `kinetics_isothermal`,
`energy_audit`, `heated0d`, `conduction`, `ATC1` and `ATC2-chemistryOff`; not in CTest).
