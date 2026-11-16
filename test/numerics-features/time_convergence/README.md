# Parallel Timing

One 1000-iteration run in eight parallel layouts. The results are on the V&V page:
[Parallel Scaling](../../../docs/vv/parallel-scaling.md).

## Setup (`simdir/`)

| Parameter | Value |
|---|---|
| Mesh | four blocks of 1000 × 20 cells, 4.0 m × 0.2 m × 0.01 m, joined by `connection` (101) |
| Faces | 302 at 1000 K (x = 0) and 2000 K (x = 4 m); `null` elsewhere |
| Materials | blocks 1–2: k = 10 W/mK; blocks 3–4: k = 100 W/mK; ρ = 8000 kg/m³, c<sub>p</sub> = 450 J/kgK |
| Initial condition | 1500 K |
| Time scheme | RK3, `vnn = 0.75`, steady (local time step) |
| Integration variables | `prim` |
| IRS | `irs-beta = 0.5` |
| Stop | `iter-threshold = 1000` (the residual falls from 3.9 × 10³ to 5.6) |

The timings need an idle node with 10 to 40 cores.

## Run

```bash
cd test/numerics-features/time_convergence/simdir
./FUSS.sh -p 10 solve > ../SOLUTION/openMP-10      # 1 x 10
./FUSS.sh -m 2 -p 5 solve > ../SOLUTION/MPI-2x5     # 2 x 5
...
cd .. && python3 plot.py                            # time_convergence.svg
```

`plot.py` reads the total time of each log in `SOLUTION/` (the `Elapsed` line of the timing
summary, or `Time of operation`) and draws the figure.

CTest runs 100 iterations (`TimeConvergence`, label `smoke`).
