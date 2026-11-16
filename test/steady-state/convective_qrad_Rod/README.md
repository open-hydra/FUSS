# Convective Rod with Imposed Flux

The physics, the reference solution and the results are on the V&V page:
[Convective Rod with Imposed Flux](../../../docs/vv/1D-steady-state-convective-qrad-rod.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | 100 × 2 cells, 0.1 m × 0.01 m, uniform |
| Face 1 (x = 0) | `T = 300` (302) |
| Face 2 (x = L) | `hconv = 500`, `Tref = 400`, `qrad = 2e4` (303) |
| Faces 3–6 | `null` |
| Material | k = 55.563 W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK |
| Initial condition | 300 K |
| Time scheme | RK3, `vnn = 1.0`, steady |
| Integration variables | `prim` |
| IRS | `irs-beta = 0.5` |
| Stop | `res-threshold = 1e-8` (iteration 43 477) |

## Check

`verify.py` compares the field and the wall output of both faces with the exact linear solution. Tolerances: 10⁻³ K on temperatures, 10⁻⁴ relative on fluxes.

## Run

```bash
cd test/steady-state/convective_qrad_Rod
./FUSS.sh solve && python3 verify.py
```

CTest: `ConvectiveQradRod` (labels fast, validation).
