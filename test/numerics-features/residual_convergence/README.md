# Residual Convergence

Six convergence strategies on one two-material plate, and three FAS cycles. The physics and the
results are on the V&V pages: [Residual Convergence](../../../docs/vv/residual-convergence.md) and
[Multigrid Cycles (FAS)](../../../docs/vv/fas-cycles.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | two blocks of 100 × 10 cells, 1.0 m × 0.1 m × 0.01 m, joined at x = 0.5 m (101) |
| Faces | 302 at 1000 K (x = 0) and 2000 K (x = 1 m); `null` elsewhere |
| Materials | k = 10 and 100 W/mK; ρ = 8000 kg/m³, c<sub>p</sub> = 450 J/kgK |
| Initial condition | 1500 K |
| Time scheme | RK3, steady (local time step) |
| Integration variables | `prim` (the `fas_*` cases: `cons`) |
| Stop | `res-threshold = 1e-8` |

| Case | `vnn` | `irs-beta` | `levels` | `level2-iter` | `cycle` |
|---|---|---|---|---|---|
| `nominal` | 0.25 | 0 | 1 | — | — |
| `irs` | 1.0 | 0.5 | 1 | — | — |
| `mg1` | 0.25 | 0 | 2 | 25 000 | — |
| `mg2` | 0.25 | 0 | 2 | 50 000 | — |
| `mg3` | 0.25 | 0 | 2 | 75 000 | — |
| `mg_irs` | 1.0 | 0.5 | 2 | 100 000 | — |
| `fas_v` | 0.25 | 0 | 2 | 1 | `v` |
| `fas_w` | 0.25 | 0 | 2 | 1 | `w` |
| `fas_w_irs` | 1.0 | 0.5 | 2 | 1 | `w` |

`level1-iter = 1000000` is only a cap: the fine level stops on `res-threshold`.

## Run

```bash
cd test/numerics-features/residual_convergence/<case> && ./FUSS.sh solve   # one case
cd test/numerics-features/residual_convergence && python3 plot.py         # the figure, from SOLUTION/
```

`SOLUTION/<case>.dat` is the `OUTPUT/residual-history.dat` of that case run to the threshold,
thinned to every 100th row. The first 100 rows, the 300 rows after each change of grid level, and
the first and last rows are kept, so the last row is the converged iteration. Re-run the case for
the full history.

CTest runs each case for 300 iterations (`ResidualConvergence_<case>`, label `smoke`);
`MultigridNoneInert` checks that the FAS keys change nothing under `cycle = none`.
