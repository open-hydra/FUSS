# Convergence figure: none / IRS / V / W / V + IRS / W + IRS on four cases

The runs and the two figures of [Multigrid Cycles (FAS)](../../../docs/vv/fas-cycles.md)
(`docs/vv/images/fas-cycles-convergence-{iterations,time}.png`). Nothing here is a CTest: the
single-grid liner runs take about 40 and 30 minutes.

| Case | Directory | Levels |
|---|---|---|
| `plate` | `residual_convergence/fas_v`: two-material plate, 2 blocks (101) | 2 |
| `liner` | `multimat_liner`: four-layer liner wall, 2 blocks (101), material jumps inside the blocks (`bc<m>.txt` for 4 levels) | 4 |
| `chimera_mg` | `chimera_mg`: two overlapping blocks (102) | 3 |
| `anisotropic` | `anisotropic_Plate`: cells of aspect ratio 100 | 3 |

```bash
./run_configs.sh /path/to/runs -j 2             # bin/FUSS of the repository, 2 runs at a time
FUSS=/path/to/FUSS ./run_configs.sh /path/to/runs -c liner,plate   # another solver, two cases
./run_configs.sh /path/to/runs -c liner -k v,w,v_irs,w_irs           # only the liner cycles
python3 plot_convergence.py /path/to/runs       # runs/convergence_{iterations,time}.png + summary
```

`run_configs.sh` copies each case, sets the keys with `set_keys.py` (RK3, `cons`,
`res-threshold = 1e-10`, `res-diter = 1`; `vnn = 0.5` without IRS, `vnn = 1.0` and
`irs-beta = 0.5` with it; `none` and `IRS` on one grid, the cycles on the case's levels with
`level<m>-iter = 1` on the coarse levels) and runs it with one thread in
`<runs>/<case>_<config>`, writing `run.log`, `run.err` and `run.time` (`exit <status> <seconds> s`).
`FUSS` may carry a prefix, e.g. `FUSS="taskset -c 0-3 /path/to/FUSS"`. The wall time of the
figure is the `Elapsed` line of each run's timing summary: one core per run, so the times are
only as reproducible as the machine is idle.
