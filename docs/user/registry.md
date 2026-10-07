# FUSS Input Parameters


## FUSS-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false |  no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 |  no | Residual convergence threshold |
| time-threshold | 1e30 | > 0 |  no | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 |  no | Maximum number of iterations |

## FUSS-IO

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| ic-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw |  no | Initial condition (INPUT/ic.*) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk raw |  no | Solution (OUTPUT/field.*) format |
| sol-diter | 1000000000 | > 0 |  no | Solution output iter frequency |
| sol-dtime | 1e30 | > 0 |  no | Solution output time frequency |
| sol-overwrite | true | true, false |  no | Overwrite solution files |
| res-diter | 1 | > 0 |  no | Residual history iter frequency |
| shell-diter | 1 | > 0 |  no | Shell update iter frequency |
| ini-diter | 10000 | > 0 |  no | input.ini update iter frequency |
| timer-diter | 0 | >= 0 |  no | Wall-clock timing report iter frequency (0 = off) |

## FUSS-Probes

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| probe1 | probe-placeholder |  |  no | Probe file name |

## probe-placeholder

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| variables | none |  |  no | Probe variables to write |
| dtime | 1e30 | > 0 |  no | Probe output time frequency |
| diter | 1000000000 | > 0 |  no | Probe output iter frequency |
| index-position | 0 0 0 0 | >= 0 |  no | Probe location by index |
| position | 0.0 0.0 0.0 |  |  no | Probe location by coordinates |

## FUSS-Numerics

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| time-scheme | euler | euler, RK2, RK3 |  no | Time integration solver |
| vnn | 0.3 | > 0 |  no | VNN parameter |
| vnn-rise-threshold | 0 | >= 0 |  no | VNN rise threshold |
| time-accurate | .false. | logical | yes | Time accurate switch |
| integration-variables | cons | cons ,  prim |  no | Integration variables (cons/prim) |
| irs | .false. | logical |  no | Implicit Residual Smoothing |
| irs-beta | 0.0 | >= 0 |  no | IRS beta parameter |

## FUSS-Multigrid

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| level1-iter | 0 | >= 0 |  no | Iterations for multigrid level 1 |
| level2-iter | 0 | >= 0 |  no | Iterations for multigrid level 2 |

## FUSS-MeshMotion

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| law | static | static ,  prescribed ,  translation |  no | Mesh motion law |
| vel | 0.0 | real |  no | Rigid-translation velocity [m/s] |
| amp | 0.0 | >= 0 |  no | Prescribed-motion displacement amplitude per coordinate [m] |
| kx | 0.0 | real |  no | Prescribed-motion wavenumbers multiplying x [1/m] |
| ky | 0.0 | real |  no | Prescribed-motion wavenumbers multiplying y [1/m] |
| kz | 0.0 | real |  no | Prescribed-motion wavenumbers multiplying z [1/m] |
| omega | 0.0 | real |  no | Prescribed-motion angular frequency [1/s] |
| gcl-tolerance | 1.0e-10 | > 0 |  no | Max relative discrete-GCL residual before the run is stopped |
| taper-to-boundary | .false. | logical |  no | Taper prescribed motion to zero on the block boundary |
