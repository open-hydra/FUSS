# 1-D Radiative Rod

A rod is held at 1000 K at one end and radiates to 300 K at the other. The test checks the radiation wall (BC 304).

**Reference**: NAFEMS P16, *Benchmark Tests for Thermal Analysis*, Test 9 (iii), 1986.
Setup and run: `test/steady-state/radiative_Rod/README.md`.

---

## Problem

A rod of length L = 0.1 m, k = 55.563 W/(m·K).

| Face | Condition |
|---|---|
| x = 0 | T = 1000 K (302) |
| x = L | radiation, ε = 0.98, T<sub>ref</sub> = 300 K (304) |
| lateral | adiabatic |

Without a source the steady profile is linear. The free-end temperature $T_L$ balances conduction
and radiation:

$$\frac{k\,(T_0 - T_L)}{L} = \varepsilon\,\sigma\left(T_L^4 - T_\text{ref}^4\right)$$

Its root is 926.9707 K; NAFEMS rounds it to 927 K.

## Result

| Location | NAFEMS | FUSS | Error |
|---|---|---|---|
| x = 0.1 m | 927 K | 926.971 K | 0.0032 % |

The tolerance is 1 %. The difference from NAFEMS is the rounding of the reference: FUSS agrees with
the exact root to the last printed digit.
