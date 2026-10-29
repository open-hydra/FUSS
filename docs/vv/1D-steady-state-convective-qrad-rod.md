# 1-D Steady State Thermal Analysis – Convective Rod with Imposed Flux

Steady conduction in a rod with a prescribed temperature at one end and, at the other, a convective wall that also receives a prescribed heat flux (BC 303: `hconv`, `Tref`, `qrad`). With a constant conductivity the solution is linear and known in closed form, so the test checks that both parts of the BC 303 flux reach the solution and the wall output.

---

## Problem setup

A rod of length L = 0.1 m (0.01 m thick) at T₀ = 300 K at x = 0. At x = L the heat flux into the rod is

$$q = h\,(T_\text{ref} - T_L) + q_\text{rad}$$

and in steady state it crosses the rod to the cold end: $k\,(T_L - T_0)/L = q$. Hence

$$T_L = \frac{k\,T_0/L + h\,T_\text{ref} + q_\text{rad}}{k/L + h}, \qquad T(x) = T_0 + (T_L - T_0)\,\frac{x}{L}$$

**Boundary conditions**

| Face | Condition | Value |
|---|---|---|
| Face 1 (x = 0) | Prescribed temperature (302) | T = 300 K |
| Face 2 (x = L) | Convection + imposed flux (303) | h = 500 W/m²K, T_ref = 400 K, q_rad = 2 × 10⁴ W/m² |
| Faces 3–6 | `null` (adiabatic) | — |

**Material properties**: k = 55.563 W/mK, ρ = 7850 kg/m³, c_p = 460 J/kgK (constant).

## Numerical setup

| Parameter | Value |
|---|---|
| Mesh | 100 × 2 cells, uniform |
| Time scheme | RK3, VNN 1.0, steady (local time step) |
| Integration variables | prim |
| Implicit residual smoothing | β = 0.5 |
| Initial condition | 300 K |
| Residual threshold | 1 × 10⁻⁸ of the first residual (reached at iteration 43 477) |

## Results and verification

`verify.py` compares the field and the wall output with the exact solution: T_L = 366.311 K, flux through the rod 36 844.44 W/m² (without q_rad, T_L would be 347.365 K).

| Quantity | Error | Tolerance |
|---|---|---|
| Field, max \|T − T_exact\| | 2.9 × 10⁻⁶ K | 10⁻³ K |
| Face 2, max \|T_w − T_L\| | 2.7 × 10⁻⁶ K | 10⁻³ K |
| Face 2, max \|q_w − q\| / q | 3.6 × 10⁻⁸ | 10⁻⁴ |
| Face 1, max \|q_w + q\| / q (the flux leaves the rod) | 8.8 × 10⁻⁸ | 10⁻⁴ |

The finite-volume scheme is exact for a linear profile, so the errors are those of the residual convergence. Before October 2026 BC 303 dropped `qrad` (the flux and the wall output read a field that only BC 301 sets): the run ended at T_L = 347.37 K and this test failed.

```bash
cd test/steady-state/convective_qrad_Rod
./FUSS.sh solve && python3 verify.py
```
