# FUSS Input Parameters


## FUSS-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false |  no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 |  no | Residual convergence threshold |
| time-threshold | 1e30 | > 0 |  no | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 |  no | Maximum number of iterations (re-read at runtime) |

## FUSS-IO

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| ic-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk binary, vtk raw |  no | Initial condition (INPUT/ic.*) format |
| sol-format | tecplot ascii | tecplot ascii, tecplot binary, vtk ascii, vtk binary, vtk raw |  no | Solution (OUTPUT/field.*) format |
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
| dt-face-k | .false. | logical |  no | Cell time step from the larger face conductivity of each direction, evaluated as the fluxes evaluate it, instead of the cell k: needed where a block holds a material jump, at which the first cells of the less conductive material exchange heat through up to their neighbour's k. Changes no step on one material at constant k; a second pass over the cells with up to 9 property lookups each (re-read at runtime) |
| time-accurate | .false. | logical | yes | Time accurate switch |
| integration-variables | cons | cons ,  prim |  no | Integration variables (cons/prim) |
| irs | .false. | logical |  no | Implicit Residual Smoothing (switched on only at setup) |
| irs-beta | 0.0 | >= 0 |  no | IRS beta parameter |

## FUSS-Multigrid

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| level1-iter | 0 | >= 0 |  no | Iterations for multigrid level 1 (re-read at runtime) |
| level2-iter | 0 | >= 0 |  no | Iterations for multigrid level 2 (re-read at runtime) |
| prolongated-diter | 0 | >= 0 |  no | Coarse levels: least level iterations between two writes of field/wall-prolongated, which are always written at the end of a level (0 = with every field write; re-read at runtime) |
| cycle | none | none, v, w |  no | FAS cycle on the levels below the running one: none = grid sequencing only, v, w (steady runs, integration-variables cons). w is the robust choice: v without IRS can limit-cycle on 4 or more levels |
| pre-iter | 2 | >= 0 |  no | FAS: smoothing iterations on a level before its coarse-grid visit |
| post-iter | 1 | >= 0 |  no | FAS: smoothing iterations on a level after the coarse-grid correction |
| coarse-iter | 10 | > 0 |  no | FAS: smoothing iterations on the coarsest level of a cycle |
