# Parallel Execution

FUSS runs on OpenMP threads, on MPI ranks, or on both (one MPI rank per socket or node, threads inside it). The build options are `--use-openmp` and `--use-mpi` of `install.sh build` (see [Installation](../getting-started/installation.md)).

## OpenMP

The cell loops of every block are shared among the threads. The number of threads is `OMP_NUM_THREADS`:

```bash
./FUSS.sh -p 8 solve          # 8 threads, 1 process
```

`FUSS.sh` also sets `ulimit -s unlimited` and `KMP_STACKSIZE=100M`, which the Intel runtime needs for the per-thread temporaries.

The metrics of the ghost cells of the symmetry faces (300) are filled one layer at a time, every face of a layer before the next: on a one-cell-thick direction the second layer of one face is a first-layer ghost cell of the opposite face, which another thread would otherwise be writing in the same loop. The set-up of the metrics is the same for any number of threads.

## MPI

```bash
./FUSS.sh -m 4 solve          # 4 ranks (mpirun -np 4)
./FUSS.sh -m 4 -p 6 solve     # 4 ranks x 6 threads
```

**Partition.** The unit of work is a whole block: at setup the blocks are dealt to the ranks, the largest first, each to the rank with the least cells so far. The log reports the result:

```
  FUSS MPI partition: 3 blocks over 2 ranks, balance 100.0% of ideal
```

**Decomposition map.** `decomposition-map = <file>` in `[FUSS-IO]` takes the rank of every block from a file of `INPUT/` (the name only, with the phase prefix of `bc.txt`; default `none`: the partition above). The file has the format of the map of the ATLAS mesh decomposer: lines starting with `#` are comments; then one line per block, blocks 1 to *n* in order, with 13 integers `new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank`, of which FUSS uses `new` (the block), `ni nj nk` and `cells` (checked against the mesh) and `rank` (from 0); the first `#` line after the block lines ends them. Rank 0 reads and checks the file and passes the result to the other ranks. A file that does not exist or is a directory, a line that is not 13 integers, blocks out of order or missing, sizes that differ from the mesh, a rank outside the run, or (with at least as many blocks as ranks) a highest rank below the last, which says the map is for another number of ranks, stops every rank with exit code 1. With fewer blocks than ranks the ranks without a block idle, with a warning. A map in `INPUT/` is read only when the key names it.

A block is never split. With more ranks than blocks the extra ranks have no work, and the largest block caps the balance: to use more ranks, split the mesh into more blocks.

**Communication.** At every stage each rank exchanges the ghost cells of the block connections (BC 101) that cross ranks and the chimera donor cells (BC 102) of its receivers. With multigrid the partition is made once, on the finest level, and every level has its own exchange schedule. Output goes through rank 0: it gathers the temperature, material and source fields and the wall values, and writes the files. MPI is initialised with `MPI_THREAD_FUNNELED`: every MPI call is made by the primary thread of each rank, the other threads of the rank work on the packing, the unpacking and the boundary conditions between the calls.

**Stops.** An error found by rank 0 (the input file, the material tables, the boundary conditions) is passed to the other ranks, which stop with it; an error that every rank finds on the same data (a line of `bc.txt`, a time table, a probe, the FAS keys) stops every rank after a barrier and `MPI_FINALIZE`; an error that only some ranks can find (the initial condition each rank reads, a probe file) is reduced over the ranks first. Each rank writes `[RANK r] STOP: <message>`, the run ends with exit code 1, and the output of every rank is flushed before it ends. See [Running FUSS](using.md#exit-codes).

**Reproducibility.** The results do not depend on the number of OpenMP threads, of MPI ranks or on how the blocks are split over the ranks: any layout of ranks × threads writes identical files, bit for bit, including `residual-history.dat`, whose sum is taken in a fixed order (over the cells of each block, line by line, then over the blocks in turn). The CTest label `needs-mpi` checks it on the chimera and multigrid cases with 2 ranks and, on a case with a box of four blocks and a chimera pair, on 1×1, 1×4, 2×1, 4×1 and 2×2 ranks × threads (see [Testing](../development/testing.md)).

**Limits.** BC 104 (coupling to the fluid solver across non-matching faces, in hydra) reads the partner cells in place: hydra stops when they sit on another rank, so in practice it runs on one rank.
