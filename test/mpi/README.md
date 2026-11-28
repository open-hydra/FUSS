# MPI Bit Identity and Thread Funnelling

The physics, the reference solution and the results are on the V&V page:
[Parallel Identity: the Same Bits on Any Rank/Thread Layout](../../docs/vv/parallel-identity.md).

## `bit-identity/`

A box of six blocks (a 2×2 box of 10×10×3 cells around a shared vertex, joined by `101`
connections, walls `302`/`301`/`303`/`305`, symmetry in `z`, sheared `x += 0.3 y`) and a chimera pair
(16×16×3 and 15×15×3 cells overlapping by 3×3, two chimera faces meeting at an edge in both blocks):
two cases (`make_case.py bit`, `steady` and `transient`), each solved at 1×1 (reference), 1×4, 2×1,
4×1 and 2×2 (MPI ranks × OpenMP threads). Every `OUTPUT/` file (field, diagnostic and wall solutions
in `vtk raw`, `residual-history.dat`) must be byte-identical to the 1×1 run, file for file.

## `thread-funnel/`

The same case (`steady`, 10 iterations; a second run with `irs-beta = 0.5`), at 2×2, 2×4 and 2×2 with
IRS. FUSS initialises MPI with `MPI_THREAD_FUNNELED`: a probe library preloaded into the solver
(`mpi_thread_probe.c`, built with the tests) counts the MPI calls made by any thread other than the
one that called `MPI_INIT_THREAD`. The test asserts that count is zero on every rank, that the probe
did intercept the halo exchange (`MPI_STARTALL`, `MPI_ISEND`, `MPI_IRECV` of the chimera exchange,
otherwise the zero would mean nothing), and that the largest point-to-point tag is at most 32767.

## Check

Both scripts are self-contained: they generate the case(s), run every layout, and diff the output
files themselves (`run.sh` prints `FAIL: ...` and exits 1 on the first difference).

## Run

```bash
cmake --build build -j 8   # or build-mpi/ for the mpi preset
test/mpi/bit-identity/run.sh    build/bin/FUSS   /tmp/fuss-bit-identity
test/mpi/thread-funnel/run.sh   build/bin/FUSS   /tmp/fuss-thread-funnel   # needs FUSS_MPI_THREAD_PROBE
```

CTest (`needs-mpi` tier, registered only in a tree built with `USE_MPI`): `MPIBitIdentity`,
`MPIThreadFunnel`. The pyrolysis, Darcy and ablating-wall MPI checks are run by
`test/run_mpi_gate.sh`, the CTest test `Gate_MPI` of the same tree (`ctest -L mpi-gate`; see the V&V page for the case list).
