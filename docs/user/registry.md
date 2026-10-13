# FUSS Input Parameters


## FUSS-Parameters

| Parameter | Default | Allowed | Required | Description |
|-----------|---------|---------|----------|-------------|
| newrun | true | true , false |  no | Start a new simulation (false = restart) |
| res-threshold | 1e-10 | > 0 |  no | Residual convergence threshold (not used under law = nodeshift: a receding surface is never steady, the run ends on its time or iteration limit) |
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
| law | static | static ,  prescribed ,  translation ,  nodeshift |  no | Mesh motion law (nodeshift: the ablating faces 306 recede, plan 08) |
| vel | 0.0 | real |  no | Rigid-translation velocity [m/s] |
| amp | 0.0 | >= 0 |  no | Prescribed-motion displacement amplitude per coordinate [m] |
| kx | 0.0 | real |  no | Prescribed-motion wavenumbers multiplying x [1/m] |
| ky | 0.0 | real |  no | Prescribed-motion wavenumbers multiplying y [1/m] |
| kz | 0.0 | real |  no | Prescribed-motion wavenumbers multiplying z [1/m] |
| omega | 0.0 | real |  no | Prescribed-motion angular frequency [1/s] |
| gcl-tolerance | 1.0e-10 | > 0 |  no | Max relative discrete-GCL residual before the run is stopped |
| taper-to-boundary | .false. | logical |  no | Taper prescribed motion to zero on the block boundary |
| skip-bc-metrics-refresh | .false. | logical |  no | TEST ONLY: do not rebuild the ghost-cell metrics after the mesh moves |

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
| max-iter | 100000 | > 0 |  no | Sweep cap per material step; reached -> error stop 3 |
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
| seb-max-iter | 50 | > 0 |  no | Iteration cap of the surface energy balance; reached -> error stop 3 |
| no-char | .false. | logical |  no | TEST ONLY: B'_ch = mdot_ch = 0 on every 306 face; h_w is still read from the B' table |
| bprime-out-of-range | clamp | clamp ,  abort |  no | A B' table query outside the table: clamp (PATO's behaviour, counted in seb.txt) or abort (error stop 3 at the first p or B'g outside) |
| min-thickness | 0.0 | >= 0 |  no | law = nodeshift: a receding line shorter than this stops the run (error stop 2) [m]; 0 = 10 x the smallest initial cell on the receding lines, at most half the shortest line, resolved at a new run's setup and kept by a restart |
| sdot-override | -1.0 | real |  no | TEST ONLY: the recession rate of every 306 face [m/s] (kinematic recession); < 0 = off |
| dump-interface-flux | .false. | logical |  no | TEST ONLY: write the flux of every 101 face to OUTPUT/<prefix>interface-flux.txt at each output |
| perturb-interface-node | .false. | logical |  no | TEST ONLY: at step 2 move the first shared 101 node of the lowest receding block by 1e-12 m in x (the interface-node assertion must stop the run) |
| smoothing-width | 5 | >= 1 |  no | law = nodeshift: nodes of the stencil that averages the surface recession rate along each surface grid line (thesis Eq. 5.1, truncated at a line end; odd; 1 = no smoothing) |
| coupled-env-file | none |  |  no | TEST ONLY: every coupled ablating wall (307) takes its exchange values from this INPUT/ file (rows t p_w rhoeUeCH h_r q_rad, linear in t) instead of hydra; none = from hydra |
