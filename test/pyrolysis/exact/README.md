# Exact-Solution Benchmarks of the Type-2 Closure

The physics, the reference solutions and the results are on the V&V page:
[Type-2 (Darcy) Gas Closure: Exact and Analytic Benchmarks](../../../docs/vv/darcy-exact.md#exact-solution-benchmarks-of-the-type-2-closure).

Three problems with closed-form or numerically-exact solutions, independent of the steady-state
cases of `../darcy/`: transient pressure diffusion, steady transpiration with Robin heat exchange, and
a decomposition front with transpiration through the char. Each is a mesh ladder (coarsest member
in-tree); `tables/` holds the `tools/make_tacot_tables.py --isothermal` set they start from.
`tools/make_exact_case.py` writes the cases, `tools/darcy_exact.py` the exact solutions.

## Cases

| Case | What it proves |
|---|---|
| `pressure_transient/N020` (+ N040, N080) | a slab sealed at one end, a small pressure step at the other: the linearised Darcy equation (`∂p/∂t = D ∂²p/∂x²`) against its separation-of-variables series; second order in space, the backward-Euler step first order in time |
| `transpiration/N020` (+ N040, N080) | both faces open, Robin heat exchange at both: the steady temperature profile and mass flux against the fixed-point solution of the coupled energy/Darcy system; second order |
| `sharp_front/N050` (+ N100, N200) | a decomposing semi-infinite slab with instantaneous decomposition and gas transpiration through the char: the Neumann-type similarity solution of `tools/darcy_exact.py`, first order by design (a cell-by-cell front) |

## Check

`verify_exact.py <kind> <case> [<case> ...]` (`kind`: `pressure_transient`, `transpiration` or
`sharp_front`) runs the mesh ladder (coarsest first), compares each grid with the exact solution of
`darcy_exact.py`, and passes by the convergence rule of `tools/conv_rule.py`: each error falls on
every grid and its observed order on the finest pair lies in [0.8 p, 1.2 p], p from
`<kind>/expected.json` (the expected orders only: pressure 2; transpiration T and mass flux 2; sharp
front, front and L2 T, 1). Today: 2.079; 2.003 and 1.978; 1.037 and 0.966. `FUSS_GATE_SELFTEST=1`
empties the order windows.

Only the coarsest grid of each kind ships in-tree (`pressure_transient/N020`, `transpiration/N020`,
`sharp_front/N050`); the finer grids of a ladder are written by the gate script at run time with
`tools/make_exact_case.py`.

## Run

```bash
cd test/pyrolysis/exact/pressure_transient/N020
./FUSS.sh solve
cd ../.. && python3 verify_exact.py pressure_transient pressure_transient/N020
```

Physics gates (not in CTest): `test/run_physics_gates.sh -c bench_pressure_transient`,
`-c bench_transpiration`, `-c bench_sharp_front` (each generates and runs its own N020/N040/N080, or
N050/N100/N200, ladder into a work directory).
