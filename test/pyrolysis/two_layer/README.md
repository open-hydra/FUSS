# two_layer: a face between two materials inside a block under pyrolysis

A slab of two materials in one block, no kinetics, constant properties: the composite-slab heat equation.

| | material 1 `HIGH` (x < A) | material 2 `LOW` (x > A) |
|---|---|---|
| `k_s` [W/(m K)] | 8.0 | 0.4 |
| `rho_s` [kg/m3] | 2800 | 280 |
| `cp_s` [J/(kg K)] | 4000 | 1000 |

`L = 0.01 m`, `A = L/2` on a cell face; face 1 (x = 0) adiabatic (301), face 2 (x = L) at 1000 K from t = 0 (302,
constant), `T(x, 0) = 300 K`; gas closure type 1 (no kinetics: no gas). The conductivity ratio is 20, the
diffusivity of `HIGH` half that of `LOW`.

- `N020/`: the committed deck (20 cells, RK3, `vnn = 0.05`, 8 s), written by
  `test/pyrolysis/tools/make_two_layer.py N020 --ncells 20`. Edit the generator, not the deck: the gate checks that a
  fresh generation equals it.
- `verify.py <case>`: the error of every field file against the exact solution
  (`test/pyrolysis/tools/two_layer_exact.py`, the eigenfunction series of the composite slab) and whether the run
  reached its end with `T` inside `[T0 - 1, TW + 1]`.
- `verify.py --ladder <N020> <N040> <N080>`: the errors at 8 s and the observed orders; pass when every order lies in
  `[1.8, 2.2]` around the scheme's second order (the reasoning is in the docstring).
- `verify.py --decay <case>`: the check of a run at the nominal `vnn` limit: it reaches its end and its error falls
  from the first output to the last (at the limit the shortest mode is not damped, for one material as for two).
- `compare_cut.py <uncut> <cut> --tol K`: the slab as one block against the slab cut at the interface into two
  blocks (`test/tools/split_blocks.py <case> 2 1`), cell by cell.

Gate: `test/run_physics_gates.sh -c pyrolysis_two_layer`: the ladder above; seven meshes (`make_two_layer.py`
`--k-ratio 20|100`, `--cells-high/--cells-low` 8/40 and 40/8, and `--ny 4 --dy 5e-4` at ratio 100) at the nominal
limits of RK3, RK2 and euler (`--decay`) and at 0.9 x those limits (second order on N and 2N); the cut of the type-2
variant `make_two_layer.py --type2`.
