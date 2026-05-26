# Getting Started

Welcome to the FUSS getting-started guide. This section will help you install FUSS and run your first examples.

!!! info "FUSS within Hydra"
    FUSS is the solid heat-conduction solver of **Hydra**. The pre-processor **ATLAS** (distributed separately from FUSS) is typically used to generate the input files that FUSS reads: initial conditions, boundary conditions, and material property tables. See the [overview page](../overview.md#hydra-cfd-suite) for an overview of the ecosystem. FUSS also exposes a built-in mesh generator (GRIB) for simple Cartesian geometries, so simple cases can be set up without ATLAS.

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Build the API and the executable from source

    [:octicons-arrow-right-24: Install FUSS](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

    Get up and running in minutes

    [:octicons-arrow-right-24: Quick start tutorial](quick-start.md)

</div>

## Scope of This Section

This getting started guide covers:

1. **Installation**
   - Step-by-step build instructions
   - Build options and customisation

2. **Quick Start**
   - Run your first case