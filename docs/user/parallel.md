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

**Decomposition map.** `decomposition-map = <file>` in `[FUSS-IO]` takes the rank of every block from a file in `INPUT/` (the name only, with the phase prefix of `bc.txt`). The default, `none`, is the partition above.

- The format is that of the ATLAS mesh decomposer: `#` comment lines, then one line per block, blocks 1 to *n* in order, with 13 integers `new parent i0 i1 j0 j1 k0 k1 ni nj nk cells rank`.
- FUSS uses `new` (the block), `ni nj nk` and `cells` (checked against the mesh) and `rank` (from 0). The first `#` line after the block lines ends them.
- Rank 0 reads and checks the file. A missing file, a malformed line, blocks out of order or missing, sizes that differ from the mesh, a rank outside the run, or a map written for another number of ranks stops every rank with exit code 1.
- With fewer blocks than ranks the ranks without a block idle, with a warning.
- A map in `INPUT/` is read only when the key names it.

A block is never split. With more ranks than blocks the extra ranks have no work, and the largest block caps the balance: to use more ranks, split the mesh into more blocks.

**Communication.**

- At every stage each rank exchanges the ghost cells of the block connections (101) that cross ranks and the chimera donors (102) of its receivers.
- With multigrid the partition is made once, on the finest level; every level has its own exchange schedule.
- Output goes through rank 0, which gathers the fields and the wall values and writes the files.
- MPI runs in `MPI_THREAD_FUNNELED` mode: only the primary thread of a rank calls MPI. The other threads pack, unpack and apply boundary conditions between the calls.

**Stops.** Every rank stops together, with `[RANK r] STOP: <message>` and exit code 1, after flushing its output (see [Running FUSS](using.md#exit-codes)).

**Reproducibility.** **The results do not depend on the number of threads, of ranks, or on the split of the blocks over the ranks**: the layouts tested write identical files, bit for bit. The residual is summed in a fixed order. The CTest label `needs-mpi` checks this on the chimera and multigrid cases with 2 ranks, and on a box of four blocks with a chimera pair at 1×1, 1×4, 2×1, 4×1 and 2×2 ranks × threads (see [Testing](../development/testing.md)).

**Limits.** BC 104 (coupling to the fluid solver across non-matching faces, in hydra) reads the partner cells in place: hydra stops when they sit on another rank, so in practice it runs on one rank.
