# Input File

FUSS is configured through a single **INI-format** file called `input.ini`, located in the case root directory.

The file is organised into parameter blocks — each block controls a specific aspect of the simulation (run control, numerics, mesh generation, materials, initial conditions, boundary conditions, output). Parameters not specified take their default values.

## File Structure

```ini
[SECTION-NAME]
parameter = value
```

Example:

```ini
[FUSS-Parameters]
res-threshold = 1e-8
iter-threshold = 500000

[FUSS-Numerics]
time-scheme = RK3
vnn = 2.0
time-accurate = false
integration-variables = cons
irs = true
irs-beta = 0.5
```

!!! warning "Default values"
    If a parameter is not specified, the **default value** is used silently.

!!! warning "Parameter names are case-sensitive"
    Incorrect names or values are ignored and the default is used instead. Check the [Parameter Registry](registry.md) for the exact spelling.

---

## Sections

| Section | Description | Reference |
|---------|-------------|-----------|
| `[FUSS-Parameters]` | Simulation control: convergence thresholds, restart | [→](registry.md#fuss-parameters) |
| `[FUSS-Numerics]` | Time scheme, VNN, IRS, integration variables | [→](registry.md#fuss-numerics) |
| `[FUSS-Multigrid]` | Multigrid acceleration settings | [→](registry.md#fuss-multigrid) |
| `[FUSS-IO]` | Output formats, frequencies, and overwrite policy | [→](registry.md#fuss-io) |
| `[FUSS-Probes]` | Point probes for time-history recording | [→](registry.md#fuss-probes) |
| `[GRIB-meshgen]` | Built-in mesh generator global settings | [User Guide](using.md#mesh) |
| `[GRIB-Block*]` | Per-block mesh geometry and resolution | [User Guide](using.md#mesh) |
| `[GPB-Phase1]` | Material property definitions | [Initial Conditions](initial-conditions.md#materials) |
| `[ICB-Block*]` | Per-block initial conditions | [Initial Conditions](initial-conditions.md) |
| `[BCB-Block*]` | Per-block face boundary condition assignments | [Boundary Conditions](boundary-conditions.md) |
| Named BC sections | BC type and parameters (temperature, flux, convection, radiation) | [Boundary Conditions](boundary-conditions.md) |

---

## Parameter Registry

The complete list of all parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code at every release — do not edit it manually.
