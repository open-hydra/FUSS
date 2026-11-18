# Parallel Execution

FUSS runs on OpenMP threads, on MPI ranks, or on both (one MPI rank per socket or node, threads inside it). The build options are `--use-openmp` and `--use-mpi` of `install.sh build` (see [Installation](../getting-started/installation.md)).

## OpenMP

The cell loops of every block are shared among the threads. The number of threads is `OMP_NUM_THREADS`:

```bash
./FUSS.sh -p 8 solve          # 8 threads, 1 process
```

`FUSS.sh` also sets `ulimit -s unlimited` and `KMP_STACKSIZE=100M`, which the Intel runtime needs for the per-thread temporaries.

## MPI

```bash
./FUSS.sh -m 4 solve          # 4 ranks (mpirun -np 4)
./FUSS.sh -m 4 -p 6 solve     # 4 ranks x 6 threads
```

**Partition.** The unit of work is a whole block: at setup the blocks are dealt to the ranks, the largest first, each to the rank with the least cells so far. The log reports the result:

```
  FUSS MPI partition: 3 blocks over 2 ranks, balance 100.0% of ideal
```

A block is never split. With more ranks than blocks the extra ranks have no work, and the largest block caps the balance: to use more ranks, split the mesh into more blocks.

**Communication.** At every stage each rank exchanges the ghost cells of the block connections (BC 101) that cross ranks and the chimera donor cells (BC 102) of its receivers. With multigrid the partition is made once, on the finest level, and every level has its own exchange schedule. Output goes through rank 0: it gathers the temperature, material and source fields and the wall values, and writes the files. MPI is initialised with `MPI_THREAD_FUNNELED`: every MPI call is made by the primary thread of each rank, the other threads of the rank work on the packing, the unpacking and the boundary conditions between the calls.

**Reproducibility.** The results do not depend on the number of OpenMP threads, of MPI ranks or on how the blocks are split over the ranks: any layout of ranks × threads writes identical files, bit for bit, including `residual-history.dat`, whose sum is taken in a fixed order (over the cells of each block, line by line, then over the blocks in turn). The CTest label `needs-mpi` checks it on the chimera and multigrid cases with 2 ranks and, on a case with a box of four blocks and a chimera pair, on 1×1, 1×4, 2×1, 4×1 and 2×2 ranks × threads (see [Testing](../development/testing.md)).

**Limits.** BC 104 (coupling to the fluid solver across non-matching faces, in hydra) reads the partner cells in place: hydra stops when they sit on another rank, so in practice it runs on one rank.
