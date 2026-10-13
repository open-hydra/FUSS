# Performance Notes

This page lists what the code already does for speed and the known costs that have not been addressed. Gains are quoted only where they were measured; none of the open items has been tried.

---

## Measuring

- The end-of-run summary prints the iterations, the solver time, the time per iteration and the elapsed time (`Mod_Timers`).
- `timer-diter = N` (`[FUSS-IO]`) prints, every N iterations, the iteration time on the slowest, fastest and mean rank, the load imbalance, and the fractions of time spent waiting for the ghost exchange and for collectives.
- Scaling measurements: [Parallel Scaling Analysis](../vv/parallel-scaling.md).

---

## In the code

| Item | Where |
|------|-------|
| `RELEASE` flags: CMake's default optimisation plus the first accepted of each list: loop unrolling (`-funroll-loops` / `-unroll`), inlining (`-finline-functions` / `-inline`), interprocedural optimisation (`-ipo` / `-flto`); host architecture in every build type (`-march=native` for GNU, `-xHost` for Intel). No `-ffast-math`. | `cmake/SetFortranFlags.cmake` |
| `diffusive-metric = inverse-mean` stores the inverse metric of each cell (`Minv`) at set-up, so a face needs one 3×3 inversion instead of three. `mean` needs none. | `Mod_Metrics`, `Lib_Diffusive` |
| One OpenMP parallel region per iteration around all Runge–Kutta stages, with orphaned `!$omp do` loops inside. | `Mod_Explicit` |
| MPI ghost and chimera-donor exchange with persistent, non-blocking requests; the local boundary conditions are processed while the messages are in flight. | `Lib_Ghost` (`Fill_Ghost_Cell`), `Mod_GhostExchange` |
| After set-up each rank frees the large arrays of the blocks it does not own (keeping the chimera donors it needs). | `Mod_Allocate_Data` (`deallocate_remote_computation_data`) |
| The IRS / no-IRS choice is made once per stage, outside the cell loops. | `Lib_Newstate` |
| `co_T` finds the row of $T(h)$ by bisection when the material's Energy column does not decrease (checked at load), by a scan from $T_\text{min}$ otherwise: the same row, so the results are bit-identical. Serial solver time with `cons`, against the scan it replaced: ×3–4 at 1000–2000 K on a 1–3000 K table, ×2.0–2.2 under FAS cycles with IRS, ×1.24 at 300–373 K on a 1–1000 K table; `prim` unchanged (it calls `co_T` only in chimera ghosts and multigrid restriction). | `Lib_Solid`, `Load_Table` |

---

## Open items

### Property look-ups

- `Diffusive_Flux` calls `co_Kappa` once per interior face, at the mean of the two cell temperatures, so each cell's table is interpolated on up to $2\,n_\text{dir}$ faces per stage. Caching $k$ per cell and averaging it at the face would change the discretisation (mean of $k$ instead of $k$ at the mean temperature), and so the results.
- `Compute_dt` interpolates $k$, $\rho$ and $c_p$ in every cell at every iteration.

### Flux loops

- `Diffusive_Flux` receives the stencil as ten scalar temperatures and two 3×3 metrics per face. Whether the compiler vectorises the face loops has not been checked (`-fopt-info-vec-missed` for GNU, `-qopt-report` for Intel).
- `Fluxes_blk` sweeps each block once per direction: three loop nests, each reading `T` and updating `R`.

### Data layout

The metric is stored per cell as a derived type holding `c(3,3)` (`FUSS_tensor_3D_type`), and the cell lengths as `c(3)` (`FUSS_vector_3D_type`): an array of structures, with strided access when a loop reads one component. Moving to arrays indexed `(d1, d2, i, j, k)` would touch the metric set-up, the ghost cells, the boundary conditions and the fluxes.
