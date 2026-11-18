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

**Communication.** At every stage each rank exchanges the ghost cells of the block connections (BC 101) that cross ranks and the chimera donor cells (BC 102) of its receivers. With multigrid the partition is made once, on the finest level, and every level has its own exchange schedule. Output goes through rank 0: it gathers the temperature, material and source fields and the wall values, and writes the files.

**Reproducibility.** For a given block split, the results do not depend on the number of ranks: 1 process and N ranks write identical files, bit for bit. The CTest label `needs-mpi` checks it on the chimera and multigrid cases (see [Testing](../development/testing.md)).

**Limits.** BC 104 (coupling to the fluid solver across non-matching faces, in hydra) reads the partner cells in place: hydra stops when they sit on another rank, so in practice it runs on one rank.
