# Residual Convergence – Numerical Acceleration Strategies

Steady-state heat conduction across a two-material plate in the X direction. The physical problem is simple enough to admit an exact analytical solution, which makes it an ideal vehicle for assessing the convergence behaviour of different numerical acceleration strategies available in FUSS: baseline explicit RK3, implicit residual smoothing (IRS), and multigrid. Six configurations are compared — varying the IRS flag, the multigrid level-2 iteration budget, and their combination — all driven to the same residual threshold of 1 × 10⁻⁸ on the same 100 × 10 cells-per-block mesh.

---

## Problem setup

A rectangular plate (1.0 m × 0.1 m × 0.01 m) is partitioned into two material zones in the X direction:

| Block | X extent (m) | Material | k (W/mK) | ρ (kg/m³) | c_v (J/kgK) |
|---|---|---|---|---|---|
| Block 1 | 0.0 – 0.5 | mat1 | 10.0 | 8000 | 450 |
| Block 2 | 0.5 – 1.0 | mat2 | 100.0 | 8000 | 450 |

The two blocks are connected at x = 0.5 m via a `connection` (continuity) interface. No heat flux crosses the remaining faces (2-D plane geometry).

**Boundary conditions**

| Face | Location | Condition | Value |
|---|---|---|---|
| Face 1 – Block 1 | x = 0.0 m | Prescribed temperature – cold wall | T = 1000 K |
| Face 2 – Block 2 | x = 1.0 m | Prescribed temperature – hot wall  | T = 2000 K |
| All other faces   | y, z directions | Adiabatic (2-D plane) | q = 0 |

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

All cases use RK3 time-marching in pseudo-time (time-accurate = false) with primitive integration variables on the same 100 × 10 cells-per-block mesh (2 000 cells total) and are run until the energy residual drops below res-threshold = 1 × 10⁻⁸. The cases differ in IRS activation, multigrid activation, and the level-2 iteration budget:

| Case | VNN | IRS | IRS β | Multigrid levels | level2-iter |
|---|---|---|---|---|---|
| nominal | 0.5 | No  | —   | — | — |
| irs     | 2.0 | Yes | 0.5 | — | — |
| mg1     | 0.5 | No  | —   | 2 |  25 000 |
| mg2     | 0.5 | No  | —   | 2 |  50 000 |
| mg3     | 0.5 | No  | —   | 2 |  75 000 |
| mg_irs  | 2.0 | Yes | 0.5 | 2 | 100 000 |

IRS (Implicit Residual Smoothing) allows a larger effective CFL number by filtering the explicit residual before the update step. The multigrid strategy uses 2 levels; `level2-iter` is the maximum number of iterations spent on the coarse-grid correction per multigrid cycle.

## Results

Residual convergence histories for all six cases are shown below.

![Residual history comparison](images/residual_convergence.svg)

Summary of iterations required to reach res-threshold = 1 × 10⁻⁸:

| Case | Initial residual | Final residual | Iterations |
|---|---|---|---|
| nominal | 6.88 × 10² | ~1.0 × 10⁻⁸ | 369 089 |
| irs     | 2.00 × 10³ | ~1.0 × 10⁻⁸ |  99 305 |
| mg1     | 4.87 × 10² | ~1.0 × 10⁻⁸ | 294 094 |
| mg2     | 4.87 × 10² | ~1.0 × 10⁻⁸ | 219 100 |
| mg3     | 4.87 × 10² | ~1.0 × 10⁻⁸ | 144 106 |
| mg_irs  | 1.42 × 10³ | ~1.0 × 10⁻⁸ |  80 603 |

**Key observations**

- **IRS** alone reduces the iteration count by ≈ 3.7× relative to the nominal RK3 run (369 089 → 99 305), because the larger allowable VNN = 2.0 effectively damps low-frequency error components more aggressively per iteration.
- **Multigrid** alone improves with a larger coarse-grid iteration budget: increasing `level2-iter` from 25 000 (mg1) to 75 000 (mg3) cuts the total iteration count from 294 094 to 144 106. Even the best multigrid-only case (mg3) remains less efficient than IRS alone (99 305).
- **Combining multigrid and IRS** (mg_irs) yields the best result of all configurations, converging in 80 603 iterations — about 19 % faster than IRS alone and 4.6× faster than the nominal baseline.
- All six cases reach the same residual threshold of 1 × 10⁻⁸, confirming that none of the acceleration strategies compromises the converged solution.