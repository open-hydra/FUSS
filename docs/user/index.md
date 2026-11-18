# User Guide

How to set up, run and restart a FUSS case.

<div class="grid cards" markdown>

-   :material-play-circle-outline:{ .lg .middle } __Running FUSS__

    ---

    Case directory, mesh, launching the solver, log, output files, probes, restart.

    [:octicons-arrow-right-24: Running FUSS](using.md)

-   :material-file-document-outline:{ .lg .middle } __Input File__

    ---

    `input.ini`: syntax, the sections FUSS reads, defaults and validation.

    [:octicons-arrow-right-24: Input file](input.md) · [Parameter registry](registry.md)

-   :material-database-outline:{ .lg .middle } __Initial Conditions__

    ---

    Initial field `INPUT/ic.*`, heat source `INPUT/st.*`, materials (`phase.txt`, `properties.dat`), restart.

    [:octicons-arrow-right-24: Initial conditions](initial-conditions.md)

-   :material-arrow-left-right:{ .lg .middle } __Boundary Conditions__

    ---

    `INPUT/bc.txt` format, BC codes, wall model, time tables, ATLAS keywords.

    [:octicons-arrow-right-24: Boundary conditions](boundary-conditions.md)

-   :material-server-network:{ .lg .middle } __Parallel Execution__

    ---

    OpenMP threads, MPI ranks, block partition.

    [:octicons-arrow-right-24: Parallel execution](parallel.md)

</div>

If you are new to FUSS, start from [Getting Started](../getting-started/index.md).
