# FUSS Input Parameters


## FUSS-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false |  no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 |  no | Residual convergence threshold (not used under law = nodeshift: a receding surface is never steady, the run ends on its time or iteration limit) |
| time-threshold | 1e30 | > 0 |  no | Maximum simulation time |
| iter-threshold | 1000000000 | > 0 |  no | Maximum number of iterations (re-read at runtime); with more than one multigrid level each level stops at its level<m>-iter instead |

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
| decomposition-map | none |  |  no | File of INPUT/ (name only, with the phase prefix) that gives the MPI rank of every block, in the format of the ATLAS mesh decomposer; none: blocks balanced by cell count |
| ini-error | warn | warn, stop |  no | input.ini read again with an error: warn (the reading is not taken, the run goes on) or stop |
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
| vnn | 0.15 | > 0 |  no | Diffusion number: dt = vnn rho cp / sum_d k_d |grad xi_d|^2 (limits 0.5 euler/RK2, 0.628 RK3; MOSE's convention; until October 2026 twice this value) |
| vnn-rise-threshold | 0 | >= 0 |  no | VNN rise threshold |
| time-accurate | .false. | logical | yes | Time accurate switch |
| integration-variables | cons | cons ,  prim |  no | Integration variables (cons/prim) |
| irs-beta | 0.0 | >= 0 |  no | Implicit residual smoothing coefficient: 0 = off, > 0 switches IRS on (only at set-up; the value is re-read at runtime) |
| diffusive-metric | stencil | stencil, mean |  no | Face metric of the diffusive fluxes: stencil takes the face Jacobian from the same differences of the cell centres as the face gradient, exact for a linear field on any mesh except at the chimera faces (102), which take the two-sided balance with the cell metrics with every value; mean takes the mean of the two cell metrics at the interior faces (the block connections the two-sided balance, the walls the cell metric) and overstates the gradient by (1+r)^2/(4r) at a cell-size ratio r; read at set-up only |

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

## FUSS-MeshMotion

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| law | static | static ,  prescribed ,  translation ,  nodeshift |  no | Mesh motion law (nodeshift: the ablating faces 306 recede, plan 08) |
| vel | 0.0 |  |  no | Rigid-translation velocity [m/s] |
| amp | 0.0 | >= 0 |  no | Prescribed-motion displacement amplitude per coordinate [m] |
| kx | 0.0 |  |  no | Prescribed-motion wavenumbers multiplying x [1/m] |
| ky | 0.0 |  |  no | Prescribed-motion wavenumbers multiplying y [1/m] |
| kz | 0.0 |  |  no | Prescribed-motion wavenumbers multiplying z [1/m] |
| omega | 0.0 |  |  no | Prescribed-motion angular frequency [1/s] |
| gcl-tolerance | 1.0e-10 | > 0 |  no | Max relative discrete-GCL residual before the run is stopped |
| taper-to-boundary | .false. | logical |  no | Taper prescribed motion to zero on the block boundary |
| skip-bc-metrics-refresh | .false. | logical |  no | TEST ONLY: do not rebuild the ghost-cell and face metrics after the mesh moves |

## FUSS-Pyrolysis

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| gas-closure | type1 | type1 ,  type2 |  no | Pyrolysis-gas closure: type1 = quasi-steady 1-D outflow (type2 = Darcy, plan 07) |
| outflow-direction | none | none ,  +i ,  -i ,  +j ,  -j ,  +k ,  -k |  no | Index direction along which the pyrolysis gas leaves the (single) block |
| p-gas | 101325.0 | > 0 |  no | Prescribed pyrolysis-gas pressure [Pa] at which the gas table is evaluated (type1) |
| hbar-guard-report | .false. | logical |  no | Print at setup the materials whose hbar guard fired (rho_v = rho_c, no net decomposition) |

## FUSS-Darcy

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| solver | tridiag | tridiag ,  pseudo |  no | Pressure solver: tridiag = line-implicit sweeps along the open-face direction of each block (pseudo on a block whose open faces lie in several directions, or none); pseudo = damped Jacobi everywhere |
| tol | 1.0e-8 | > 0 |  no | Stop test of the pressure sweeps: max cell mass residual <= tol * max(|face flux|, |source|, |storage|) |
| max-iter | 100000 | > 0 |  no | Sweep cap per material step; reached -> the run stops (exit code 1) |
| p-ref | 101325.0 | > 0 |  no | Pressure scale [Pa] of the reported per-sweep change (not the stop test) |
| p-init | 101325.0 | > 0 |  no | Uniform gas pressure [Pa] on newrun; ignored on restart |
| p-wall | 101325.0 | > 0 |  no | Dirichlet gas pressure [Pa] on an open face whose bc.txt line carries no p_w of its own |
| open-bc-types | 303 305 306 |  |  no | Blank-separated BC codes whose faces are open to the gas (p = p_w); from {301,302,303,304,305,306,103} |
| quasi-steady | .false. | logical |  no | Drop the gas storage: d(eta p)/dt from the pressure equation and eps (rho_g h_g - p) from the energy E; the gas released but not yet expelled is carried to the next step as (eta_st - eta_prev) p |
| klinkenberg | .false. | logical |  no | K_eff = K + (mu/p) B with the tensors B_v, B_c of pyrolysis.dat (thesis eq. 3.15) |
| pseudo-vnn | 0.5 | > 0 |  no | Damping of the Jacobi (pseudo) sweep: p <- p + pseudo-vnn * residual / diagonal; (0, 1] |
| axes-follow-grid | .false. | logical |  no | Rotate the permeability tensor with the grid lines (not built: false only) |
| pi-source | 0.0 | >= 0 |  no | TEST ONLY: constant gas production [kg/m3/s] added to the kinetics in every cell |

## FUSS-Ablation

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| seb-tol | 1.0e-8 | > 0 |  no | Stop test of the surface energy balance of BC 306: |dT_w| <= seb-tol*T_w, |F| <= seb-tol*max|term| and |d(beta)| <= seb-tol |
| seb-max-iter | 50 | > 0 |  no | Iteration cap of the surface energy balance; reached -> the run stops (exit code 1) |
| no-char | .false. | logical |  no | TEST ONLY: B'_ch = mdot_ch = 0 on every 306 face; h_w is still read from the B' table |
| bprime-out-of-range | clamp | clamp ,  abort |  no | A B' table query outside the table: clamp (PATO's behaviour, counted in seb.txt) or abort (the run stops (exit code 1) at the first p or B'g outside) |
| min-thickness | 0.0 | >= 0 |  no | law = nodeshift: a receding line shorter than this stops the run (exit code 1) [m]; 0 = 10 x the smallest initial cell on the receding lines, at most half the shortest line, resolved at a new run's setup and kept by a restart |
| sdot-override | -1.0 |  |  no | TEST ONLY: the recession rate of every 306 face [m/s] (kinematic recession); < 0 = off |
| dump-interface-flux | .false. | logical |  no | TEST ONLY: write the flux of every 101 face to OUTPUT/<prefix>interface-flux.txt at each output |
| perturb-interface-node | .false. | logical |  no | TEST ONLY: at step 2 move the first shared 101 node of the lowest receding block by 1e-12 m in x (the interface-node assertion must stop the run) |
| smoothing-width | 5 | >= 1 |  no | law = nodeshift: nodes of the stencil that averages the surface recession rate along each surface grid line (thesis Eq. 5.1, truncated at a line end; odd; 1 = no smoothing) |
| coupled-env-file | none |  |  no | TEST ONLY: every coupled ablating wall (307) takes its exchange values from this INPUT/ file (rows t p_w rhoeUeCH h_r q_rad, linear in t) instead of hydra; lines // face b i j k f split it into one block per 307 face (a replay of hydra's env-dump); none = from hydra |
