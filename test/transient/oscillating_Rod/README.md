# Oscillating Rod

The physics, the reference solution and the results are on the V&V page:
[Oscillating Rod](../../../docs/vv/1D-transient-oscillating-rod.md).

## Setup

| Parameter | Value |
|---|---|
| Mesh | 100 × 2 cells, 0.1 m × 0.01 m (Δx = 1 mm) |
| Face 1 (x = 0) | `type = wall`, `T = 2000` |
| Face 2 (x = 0.1 m) | `type = wall`, `T-time-file = Tw_time.dat` (one period, `periodic`) |
| Faces 3–6 | `null` |
| Material | k = 1 W/mK, ρ = 2300 kg/m³, c<sub>p</sub> = 985 J/kgK |
| Initial condition | 2000 K |
| Time scheme | RK3, `vnn = 0.25`, `time-accurate = true` (Δt = 0.545 s) |
| Integration variables | `prim` |
| End time | 160 s (two periods); field every 10 s |
| Probes | `p1` at (1, 90, 1, 1), x ≈ 0.0895 m; `p2` at (1, 91, 1, 1), x ≈ 0.0905 m; `dtime = 0.1` |

`bc_time.py` writes `Tw_time.dat`: T(t) = 2000 + 1000 sin(2π t / 80 s), one row per second.

## Check

`verify.py` takes the mean of the two probes, compares its peak rise with 94.69 K and passes below 1 %. It writes the two figures of the V&V page.

## Run

```bash
cd test/transient/oscillating_Rod
./FUSS.sh solve && python3 verify.py
```

CTest: `OscillatingRod` (labels fast, validation).
