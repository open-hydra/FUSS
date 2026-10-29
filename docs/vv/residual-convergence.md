# Residual Convergence – Numerical Acceleration Strategies

Steady-state heat conduction across a two-material plate in the X direction. The physical problem is simple enough to admit an exact analytical solution, which makes it an ideal vehicle for assessing the convergence behaviour of different numerical acceleration strategies available in FUSS: baseline explicit RK3, implicit residual smoothing (IRS), and multigrid. Six configurations are compared — varying the IRS flag, the multigrid level-2 iteration budget, and their combination — all driven to the same residual threshold of 1 × 10⁻⁸ on the same 100 × 10 cells-per-block mesh.

---

## Problem setup

A rectangular plate (1.0 m × 0.1 m × 0.01 m) is partitioned into two material zones in the X direction:

| Block | X extent (m) | Material | k (W/mK) | ρ (kg/m³) | c_p (J/kgK) |
|---|---|---|---|---|---|
| Block 1 | 0.0 – 0.5 | mat1 | 10.0 | 8000 | 450 |
| Block 2 | 0.5 – 1.0 | mat2 | 100.0 | 8000 | 450 |

The two blocks are joined at x = 0.5 m by a block connection (101). No heat crosses the other faces (`null`).

**Boundary conditions**

| Face | Location | Condition | Value |
|---|---|---|---|
| Face 1 – Block 1 | x = 0.0 m | Prescribed temperature (302) – cold wall | T = 1000 K |
| Face 2 – Block 2 | x = 1.0 m | Prescribed temperature (302) – hot wall  | T = 2000 K |
| All other faces   | y, z directions | `null` (no flux) | — |

**Initial condition**

| Domain | T_init |
|---|---|
| All blocks | 1500 K (uniform) |

## Analytical solution

At steady state, the heat flux q is uniform and determined by the total thermal resistance of the two layers in series:

$$q = \frac{T_\text{hot} - T_\text{cold}}{L_1/k_1 + L_2/k_2} = \frac{2000 - 1000}{0.5/10 + 0.5/100} = \frac{1000}{0.055} \approx 18182\ \text{W/m}^2$$

The interface temperature is:

$$T_\text{int} = T_\text{cold} + q\,\frac{L_1}{k_1} = 1000 + 18182 \times 0.05 \approx 1909\ \text{K}$$

The temperature profile is piecewise linear in each block:

$$T(x) = \begin{cases} 1000 + 1818\,x & 0 \le x \le 0.5\ \text{m} \\ 2000 - 182\,(1-x) & 0.5 \le x \le 1.0\ \text{m} \end{cases}$$

## Cases and numerical setup

All cases use RK3 time-marching in pseudo-time (time-accurate = false) with primitive integration variables on the same 100 × 10 cells-per-block mesh (2 000 cells total) and are run until the residual has fallen by res-threshold = 1 × 10⁻⁸ from its first value on the finest level (the test is relative since October 2026, see [Time Integration](../theory/time-integration.md)). The cases differ in IRS activation, multigrid activation, and the level-2 iteration budget:

| Case | VNN | IRS | IRS β | Multigrid levels | level2-iter |
|---|---|---|---|---|---|
| nominal | 0.25 | No  | —   | — | — |
| irs     | 1.0 | Yes | 0.5 | — | — |
| mg1     | 0.25 | No  | —   | 2 |  25 000 |
| mg2     | 0.25 | No  | —   | 2 |  50 000 |
| mg3     | 0.25 | No  | —   | 2 |  75 000 |
| mg_irs  | 1.0 | Yes | 0.5 | 2 | 100 000 |

IRS (Implicit Residual Smoothing) smooths the residual before the update, which allows a larger VNN. Multigrid is grid sequencing: level 2 (50 × 5 cells per block) runs `level2-iter` iterations, or until its residual has fallen by `res-threshold` from its first value on level 2; its solution is prolongated once to level 1, which runs until its residual has fallen by `res-threshold` from the first level-1 residual (`level1-iter = 1 000 000` is only a cap). The coarse level is not visited again.

## Results

Residual convergence histories for all six cases are shown below.

![Residual history comparison](images/residual_convergence.svg)

Summary of the runs. Each ends when its residual has fallen by 10⁻⁸ from the first residual of its finest level: for `nominal` and `irs` the residual of the start, for the multigrid runs the first residual of level 1 after the prolongation, which is already about 60 (`mg_irs`: 74) times smaller than the first one of level 2. The multigrid runs therefore end about two orders deeper, and the iterations of the stop compare different depths; the last column compares the strategies at the same depth, the first level-1 iteration whose residual is at most 3 × 10⁻⁵:

| Case | Initial residual | First level-1 residual | Final residual | Iterations | Iterations to 3 × 10⁻⁵ |
|---|---|---|---|---|---|
| nominal | 6.88 × 10² | — | 6.88 × 10⁻⁶ | 236 927 | 207 060 |
| irs     | 2.00 × 10³ | — | 2.00 × 10⁻⁵ |  60 855 |  58 802 |
| mg1     | 4.87 × 10² | 8.20 | 8.20 × 10⁻⁸ | 251 958 | 132 206 |
| mg2     | 4.87 × 10² | 8.23 | 8.23 × 10⁻⁸ | 180 608 |  61 094 |
| mg3     | 4.87 × 10² | 8.23 | 8.23 × 10⁻⁸ | 163 642 |  67 584 |
| mg_irs  | 1.42 × 10³ | 19.2 | 1.92 × 10⁻⁷ |  55 523 |  29 903 |

**Key observations** (at the same depth, last column)

- **IRS** alone reduces the iteration count by ≈ 3.5× relative to the nominal RK3 run (207 060 → 58 802), because the larger allowable VNN = 1.0 effectively damps low-frequency error components more aggressively per iteration.
- **Multigrid** alone gains 1.6× (mg1, 25 000 coarse iterations) to 3.4× (mg2, 50 000). Level 2 now stops on its own test as well: in mg3 after 64 579 of its 75 000 iterations, when its residual has fallen by 10⁻⁸, so mg3 hands over a coarse solution as converged as mg2's and needs a few more level-1 iterations (67 584 against 61 094). The best multigrid-only case is about as fast as IRS alone.
- **Combining multigrid and IRS** (mg_irs, level 2 converged after 16 556 iterations) yields the best result of all configurations: 29 903 iterations, half those of IRS alone and 6.9× fewer than the nominal baseline.

All six runs approach the same solution; how close depends on where each stops: the final fields differ from the piecewise-linear profile above by at most 4.4 × 10⁻³ K (`nominal`, the shallowest stop relative to the field), 3.2 × 10⁻³ K (`irs`), 5.3 × 10⁻⁵ K (mg1–mg3) and 3.1 × 10⁻⁵ K (`mg_irs`). The iteration counts, the residuals, the final fields and the figure (with its data in `SOLUTION/*.dat`) are those of one re-measurement of the six cases as shipped, none of them capped, with the relative `res-threshold` (branch `post-merge/sync-followups`, ifx 2024.1 with `USE_IPO` on, one thread, 2026-10-29). With the absolute test of the release before, the same decks ran 369 498 (nominal), 99 415 (irs), 294 644, 223 371, 212 414 (mg1–mg3) and 79 816 (mg_irs) iterations, to a residual of 10⁻⁸ in every case. The histories in `SOLUTION/*.dat` are thinned to every 100th row (the first 100 rows, the 300 rows after each multigrid level change, the first and the last row are kept, so the last row of each file is the converged iteration of the table); the full history is reproduced by re-running the case.

```bash
cd test/numerics-features/residual_convergence/<case> && ./FUSS.sh solve   # one case
cd test/numerics-features/residual_convergence && python3 plot.py         # the figure, from SOLUTION/
# SOLUTION/<case>.dat: the OUTPUT/residual-history.dat of that case run to the threshold, thinned to every 100th row
```
- All six cases end on the same relative test, and none of the acceleration strategies changes the converged solution: the fields approach the exact profile as the residual falls.