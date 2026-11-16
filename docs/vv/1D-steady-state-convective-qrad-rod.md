# 1-D Convective Rod with Imposed Flux

A rod is held at 300 K at one end. The other
end is a convective wall that also receives a prescribed heat flux. The test checks that both parts
of the BC 303 flux reach the solution and the wall output.

Setup and run: `test/steady-state/convective_qrad_Rod/README.md`.

---

## Problem

A rod of length L = 0.1 m, k = 55.563 W/(m·K), T<sub>0</sub> = 300 K at x = 0. At x = L the heat
flux into the rod is

$$q = h\,(T_\text{ref} - T_L) + q_\text{rad}$$

with h = 500 W/(m²·K), T<sub>ref</sub> = 400 K and q<sub>rad</sub> = 2 × 10⁴ W/m². In steady state
this flux crosses the rod, $k\,(T_L - T_0)/L = q$, so

$$T_L = \frac{k\,T_0/L + h\,T_\text{ref} + q_\text{rad}}{k/L + h}, \qquad T(x) = T_0 + (T_L - T_0)\,\frac{x}{L}$$

That gives T<sub>L</sub> = 366.311 K and q = 36 844.44 W/m². Without q<sub>rad</sub>, T<sub>L</sub>
would be 347.365 K.

## Result

| Quantity | Error | Tolerance |
|---|---|---|
| Field, max \|T − T<sub>exact</sub>\| | 2.9 × 10⁻⁶ K | 10⁻³ K |
| Face 2, max \|T<sub>w</sub> − T<sub>L</sub>\| | 2.7 × 10⁻⁶ K | 10⁻³ K |
| Face 2, max \|q<sub>w</sub> − q\| / q | 3.6 × 10⁻⁸ | 10⁻⁴ |
| Face 1, max \|q<sub>w</sub> + q\| / q | 8.8 × 10⁻⁸ | 10⁻⁴ |

On face 1 the heat leaves the rod, so q<sub>w</sub> = −q. The scheme is exact for a linear
profile: the errors are those of the iteration at its stopping threshold.
