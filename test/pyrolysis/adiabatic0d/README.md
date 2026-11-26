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

`verify.py` takes a ladder of runs of the case, coarsest first: `T1000` at `vnn = 0.25` and two copies
at `vnn = 0.125` and `0.0625` (the gate stages them). Over every `OUTPUT/field<i>.tec` of each run
(the oracle integrated once on the union of the runs' field times, `tools/ode0d.check_ladder`):

- `max |T - T_or1|` and `max |chi_r - chi_or1|` fall on every rung and their order on the finest
  pair lies in [1.6, 2.4] (`tools/conv_rule.py`, p = 2): 0.0510, 0.0137, 0.00347 K, order 1.981;
  chi 1.22e-3, 3.43e-4, 8.97e-5, order 1.934;
- the spatial non-uniformity of each run below its round-off bound, `N_steps x 3 x eps x max T` for
  `T` (6.3e-11 K at 96 steps) and `N_steps x 3 x eps` for `chi`.

Reported, not gated: the exact type-2 (gas storage) minus type-1 difference (0.0696 K).

## Run

```bash
cd test/pyrolysis/adiabatic0d/T1000 && ./FUSS.sh solve
# two copies with vnn = 0.125 and 0.0625 in input.ini, solved likewise, then
cd .. && python3 verify.py T1000 <copy at 0.125> <copy at 0.0625>
```

Physics gate: `test/run_physics_gates.sh -c pyrolysis` (CTest `Gate_Pyrolysis`; together with
`kinetics_isothermal`, `energy_audit`, `heated0d`, `conduction`, `ATC1` and `ATC2-chemistryOff`).
