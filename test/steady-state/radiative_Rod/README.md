# Radiative Rod

The physics, the reference solution and the results are on the V&V page:
[Radiative Rod](../../../docs/vv/1D-steady-state-radiative-rod.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | 100 × 2 cells, 0.1 m × 0.01 m, uniform |
| Face 1 (x = 0) | `type = wall`, `T = 1000` (302) |
| Face 2 (x = L) | `type = wall`, `eps = 0.98`, `Tref = 300` (304) |
| Faces 3–6 | `null` |
| Material | k = 55.563 W/mK, ρ = 7850 kg/m³, c<sub>p</sub> = 460 J/kgK |
| Initial condition | 1000 K |
| Time scheme | RK3, `vnn = 1.0`, steady (local time step) |
| Integration variables | `prim` |
| IRS | `irs-beta = 0.5` |
| Stop | `res-threshold = 1e-8` (iteration 56 642) |

## Check

`verify.py` reads the wall temperature of face 2 and passes below 1 % from the NAFEMS value, 927 K.

## Run

```bash
cd test/steady-state/radiative_Rod
./FUSS.sh solve && python3 verify.py
```

CTest: `RadiativeRod` (labels fast, validation).
