# 2-D Convective Plate

A plate is heated on one edge, cooled by convection on two, and insulated on the last. The test checks the convective wall (BC 303) and an
adiabatic wall.

**Reference**: NAFEMS P16, *Benchmark Tests for Thermal Analysis*, Test T4, 1986.
Setup and run: `test/steady-state/convective_Plate/README.md`.

---

## Problem

A plate of 0.6 m × 1.0 m, k = 52 W/(m·K), no source.

| Edge | Condition |
|---|---|
| left (x = 0) | adiabatic |
| right (x = 0.6 m) | convection, h = 750 W/(m²·K), T<sub>∞</sub> = 0 °C |
| bottom (y = 0) | T = 100 °C |
| top (y = 1.0 m) | convection, h = 750 W/(m²·K), T<sub>∞</sub> = 0 °C |

On the convective edges $-k\,\partial T/\partial n = h\,(T - T_\infty)$. The reference is the wall
temperature at point A (x = 0.6 m, y = 0.2 m), point E of the benchmark: 18.25 °C.

## Result

Point A lies on the right wall, between two wall cells; the value is their mean.

| Location | NAFEMS | FUSS | Error |
|---|---|---|---|
| A (0.6 m, 0.2 m) | 18.25 °C | 18.26 °C | 0.004 % |

The difference is 0.012 K, against a tolerance of 1 % (in kelvin).
