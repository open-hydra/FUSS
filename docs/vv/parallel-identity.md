# Parallel Identity: the Same Bits on Any Rank/Thread Layout

A distributed-memory run must write exactly the bits a single-rank run writes: block partitioning,
halo exchange, chimera donors and the pyrolysis/Darcy reductions must never change the converged
state, only how it is computed. These checks run the general (non-pyrolysis) bit-identity and
MPI-threading harness once, and then `test/run_mpi_gate.sh`'s one-rule gate — N ranks reproduce 1 rank
— across every case touched by this branch's features: single- and multi-block pyrolysis kinetics,
the type-2 (Darcy) closure split across ranks, and the coupled ablating wall.

Setup and run: `test/mpi/README.md`; the general gate is `test/run_mpi_gate.sh` (no case directory of
its own).

---

## Bit identity across MPI ranks and OpenMP threads

`test/mpi/bit-identity`: two cases (`make_case.py bit`, steady and transient) on a sheared mesh
($x \mathrel{+}= 0.3y$): a 2×2 box of 10×10×3-cell blocks around a shared vertex (BC 302/301/303/305,
symmetry in $z$) plus a chimera pair of 16×16×3 and 15×15×3 cells overlapping by 3×3 cells, two
chimera faces meeting at an edge in both blocks (six blocks total). Steady: local time step, `cons`,
IRS (`irs-beta = 0.5`), 30 iterations. Transient: the same mesh, time-accurate RK3, `prim`, 20 steps.
Each is solved on five layouts — 1×1 (the reference), 1×4, 2×1, 4×1, 2×2 (MPI ranks × OpenMP threads)
— splitting the chimera pair and the shared vertex across ranks on 2 and 4 ranks (the least-loaded
block partition, `Mod_MPI.f90`).

**Pass criteria.** Every file `OUTPUT/` writes (field, diagnostic and wall solutions, Tecplot or raw,
`residual-history.dat`) is byte-identical to the 1×1 run on every other layout, and the set of files
written is the same; the run log must report the ranks and threads asked for.

## MPI calls only from the main thread

`test/mpi/thread-funnel`: FUSS initialises MPI with `MPI_THREAD_FUNNELED`, so every MPI call inside
the OpenMP-parallel Runge-Kutta stages (the temperature halo exchange, the IRS residual exchange, the
chimera exchange) must be issued by the thread that called `MPI_INIT_THREAD`. A probe library
(`mpi_thread_probe.c`, a `PMPI` interposer preloaded into the solver) counts MPI calls made by any
other thread, on the `bit-identity` steady case (6 blocks) at three layouts: 2×2, 2×4, and 2×2 with
IRS.

**Pass criteria.** The probe's count of non-main-thread MPI calls is exactly 0 on every rank at every
layout; the probe did intercept the halo exchange (`MPI_STARTALL`, `MPI_ISEND`, `MPI_IRECV` of the
chimera exchange) on every rank, so that a zero count is not vacuous; and the largest point-to-point
tag used is $\le32767$, the least `MPI_TAG_UB` the MPI standard guarantees (MPI 4.1, Sec. 3.2.3).

---

## `test/run_mpi_gate.sh`: N ranks reproduce 1 rank

The gate runs each listed case on $N$ MPI ranks (Intel MPI, `OMP_NUM_THREADS=1`) and on 1, and
requires every `OUTPUT/` file to be byte-identical between them, with two documented exceptions: wall
files, where a `-0.0` owned by one rank can become `+0.0` once another rank's exact zero is added in
the cross-rank sum (the sum itself is kept); and `residual-history.dat`, compared byte for byte too
(the residual is summed per block, this branch's pyrolysis/Darcy terms included, then over the blocks
in block order). A case must stop on a fixed iteration count or elapsed time, never on
`res-threshold` (a residual summed over ranks cannot be the stopping test for an N-vs-1 comparison).

Of the gate's full case list, the ones that exercise this branch's pyrolysis and ablation features on
a fixed mesh (no mesh motion):

| Case | Ranks | What it exercises |
|---|---:|---|
| `test/pyrolysis/energy_audit` | 4 | single-block type-1 kinetics, with every other rank idle through the collectives; the energy budget, probes and field files must not move |
| `test/pyrolysis/heated0d` | 2 | single-block type-1 kinetics with a volumetric source |
| `test/pyrolysis/lit/baum_atreya2007/full_N100_20blk` | 4 | type-2 kinetics split across twenty blocks in a row (one correction group): the $\chi$ halo, the $p_{pg}$ halo, the neighbour's $\Gamma$ and $h_g\Gamma$ all change across every cut as the charring front crosses it |
| `test/pyrolysis/darcy/tensor_linear_2d_2x2` | 4 | the rotated-tensor Darcy slab cut into 2×2 blocks: the connection rows' tangential pressure stencil, and the damped-Jacobi sweep control over thousands of sweeps per step, across ranks |
| `test/pyrolysis/darcy/mass_balance_1d_2x` | 2 | two blocks side by side sharing one line-implicit correction group: the group's summed line adds the two blocks' partial sums in a different order from a single block's running sum |
| `test/pyrolysis/ablation/coupled307/run307_4blk` | 4 | the coupled ablating wall 307 cut in four, both faces on the last rank: the per-rank exchange values, and the per-face outputs gathered to root |
| `test/pyrolysis/darcy/periodic_static`, `periodic_static1` | 2 | the periodic faces (BC 201) of a type-2 Darcy slab split so that the periodic pair crosses ranks: the pair's per-step check and the Darcy rows across the pair |

(The gate's other pyrolysis/ablation entries are out of scope for this page.)

**Pass criteria.** Every `OUTPUT/` file of each case above is byte-identical between the listed rank
count and 1 rank, under the wall-file and residual-history rules stated above. The per-case
`README.md` files quote the specific byte-for-byte file counts and, where measured independently of
an exact match, the relative difference against a single-block/single-rank run of the same binary
(for example `darcy/mass_balance_1d_2x`: $T$ relative $10^{-14}$, $p_{pg}$ relative $7.8\times10^{-14}$
against its own four-line, single-block source — a cross-check of the case construction, not of the
rank count).
