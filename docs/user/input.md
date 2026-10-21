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
    In a FUSS section an unknown (misspelt) key, a repeated section, an invalid value or a string longer than the
    variable that holds it is refused at setup with the key named. Check the [Parameter Registry](registry.md) for the
    exact spelling. `irs-beta > 0` turns IRS on by itself; `irs = false` beside it is refused as a contradiction.

---

## Sections

| Section | Description | Reference |
|---------|-------------|-----------|
| `[FUSS-Parameters]` | Simulation control: convergence thresholds, restart | [→](registry.md#fuss-parameters) |
| `[FUSS-Numerics]` | Time scheme, VNN, IRS, integration variables | [→](registry.md#fuss-numerics) |
| `[FUSS-Multigrid]` | Multigrid acceleration settings | [→](registry.md#fuss-multigrid) |
| `[FUSS-IO]` | Output formats, frequencies, and overwrite policy | [→](registry.md#fuss-io) |
| `[FUSS-Probes]` | Point probes for time-history recording | [→](registry.md#fuss-probes) |
| `[FUSS-Pyrolysis]` | Pyrolysis of a charring solid (its presence switches the model on) | [→](#pyrolysis) |
| `[FUSS-Darcy]` | Gas pressure solver of the type-2 gas closure (read only with `gas-closure = type2`) | [→](#type-2-gas-closure-darcy) |
| `[GRIB-meshgen]` | Built-in mesh generator global settings | [User Guide](using.md#mesh) |
| `[GRIB-Block*]` | Per-block mesh geometry and resolution | [User Guide](using.md#mesh) |
| `[GPB-Phase1]` | Material property definitions | [Initial Conditions](initial-conditions.md#materials) |
| `[ICB-Block*]` | Per-block initial conditions | [Initial Conditions](initial-conditions.md) |
| `[BCB-Block*]` | Per-block face boundary condition assignments | [Boundary Conditions](boundary-conditions.md) |
| Named BC sections | BC type and parameters (temperature, flux, convection, radiation) | [Boundary Conditions](boundary-conditions.md) |

---

## Parameter Registry

The complete list of all parameters, defaults, and allowed values is in the **[Parameter Registry](registry.md)**. It is generated automatically from the source code at every release — do not edit it manually.

---

## Pyrolysis

A section `[FUSS-Pyrolysis]` switches on the pyrolysis model of a charring solid (phase 2 of the porous-pyrolysis programme, plan 06): the solid is a set of phases, each decomposing through Arrhenius reactions, the stored energy `E(T, χ)` is the integrated variable and the pyrolysis gas leaves through one face of the (single) block. Without the section the solver runs the legacy conduction model unchanged, bit for bit.

```ini
[FUSS-Pyrolysis]
gas-closure       = type1        ; type1: quasi-steady 1-D outflow; type2: Darcy pressure equation (below)
outflow-direction = +i           ; type1 only: +i, -i, +j, -j, +k, -k, the gas leaves through that end of the block
p-gas             = 101325.0     ; type1 only: [Pa] pressure at which the gas table is evaluated
hbar-guard-report = false        ; report the materials with rho_v = rho_c at setup
```

Requirements enforced at setup (`error stop 1` with a message): `integration-variables = cons`, `time-accurate = true`, no IRS, one grid level, `law = static` or `law = nodeshift` (the receding mesh of an ablating wall BC 306: every step the pyrolysis state is remapped conservatively onto the moved cells, at most 64 reactions), no chimera (102, 104) or coupled (103) face, OpenMP only; under `type1` also one block without block connections (101: its outflow line integral does not cross block faces), an `outflow-direction` and a wall BC (301–307) on the outflow end of the block. Under `type2` blocks may be connected (101): the reaction state, the gas pressure and both conductivities cross the connection, and its gas mass and enthalpy flows are conservative (the two blocks evaluate one face flux).

`type1` is exact only for one-dimensional, quasi-steady gas flow with an impermeable back face. The gas crosses an interior face with the enthalpy at the mean temperature of the two cells and leaves the block at the wall temperature of the outflow face (second order in space; until 2026-10-15 every face took the upstream cell's temperature, first order: ATC1's 100-cell mesh was 7 K too hot at 1 mm, 1 s, now 1.9 K). `type2` solves the gas momentum (Darcy) and the gas storage instead; see [Type-2 gas closure](#type-2-gas-closure-darcy).

### Material files

Besides `INPUT/<prefix>properties.dat` (still loaded, not used by the pyrolysis path), every material of `phase.txt` needs four more files in `INPUT/`. Until ATLAS can write them, `test/pyrolysis/tools/make_tacot_tables.py --pato-dir <PATO 3.1 root> --out INPUT` writes them from PATO's TACOT data, and `make_case.py` builds complete cases.

| file | content |
|------|---------|
| `<prefix>pyrolysis.dat` | plain text, `//` comments, one `material <name>` block per material in `phase.txt` order, `key value` lines with PATO's `constantProperties` names: `nSolidPhases`, `rhoI[i]`, `epsI[i]`, `nPyroReac[i]`, `eps_g_v`, `eps_g_c`, and per reaction `F[i][j] A[i][j] E[i][j] m[i][j] n[i][j] T[i][j] h[i][j]` (`h` is parsed and ignored: the decomposition heat is in the enthalpy tables). Each key exactly once per block. Type 2 also needs the permeability tensors `K_v` and `K_c` [m²] of the virgin and the charred material, nine reals each, row-major, in the mesh frame (symmetric, positive diagonal; an isotropic `K` is `K 0 0 0 K 0 0 0 K`), blended with the charring degree like the conductivity; `B_v`, `B_c` (same layout) are read when `[FUSS-Darcy] klinkenberg = true`. Type 1 ignores the four keys. |
| `<prefix>virgin.dat` | Tecplot POINT table, one zone per material: `T cp h ki kj kk emissivity absorptivity`, integer `T` rows 1 K apart (TACOT: 256–4000 K); `h` absolute, `h(298.15 K)` = the enthalpy of formation; `ki = kj = kk` (orthotropic conductivity is not supported) |
| `<prefix>char.dat` | same layout and the same `T` range as `virgin.dat` |
| `<prefix>gas.dat` | one zone per (material, pressure level), material-major, `p` ascending: `T p M h mu` (`M` kg/mol, `h` J/kg, `mu` Pa·s), integer `T` rows 1 K apart, the same range in every zone of a material (TACOT: 200–3975 K, five levels) |

Tecplot layout of the three tables: `TITLE = "FUSS pyrolysis table"`, `VARIABLES = "T" ...`, then per zone `ZONE T="<material>"` (gas: `ZONE T="<material> p=<p>"`), `I=<rows>, F=POINT` and the rows; nothing after the last row.

### Output

The solution file gains the columns `chi_1 … chi_nsub` (reaction advancements, read back on restart; `chi` when the material has one reaction), `rho_s` and `mdot_pg` (gas mass flux through the cell's downstream face, kg m⁻² s⁻¹); `k`, `rho`, `cs` carry the blended solid conductivity, the total density and the blended specific heat, and `h` the stored energy `E` [J/m³]. `OUTPUT/<prefix>energy-budget.txt` records, at every solution write, `t E_tot M_s Q_out Q_gas_out M_gas_out Q_src` (the discrete energy and solid-mass balances close to round-off; `test/pyrolysis/energy_audit/verify.py` checks them). On a receding mesh (`law = nodeshift`) it gains `E_abl M_abl`, the energy and solid mass the receding faces swept out (counted from the start of the run, like `Q_*`), so that `E_tot + E_abl + Q_out − Q_src` and `M_s + M_abl + M_gas_out` stay constant.

### Type-2 gas closure (Darcy)

`gas-closure = type2` (plan 07) solves the gas pressure `p_pg` of the pyrolysis gas from the gas mass balance with Darcy's law, `∂(ε_g ρ_g)/∂t + ∇·(ρ_g u_g) = π`, `u_g = −(K/(μ ε_g))·∇p` (`ρ_g = p W/(ℜ T)`, `W` and `μ` from `gas.dat`, `K` the blended permeability tensor), once per material step: the conductance `ρ_g K/μ` at the start of the step, the storage coefficient `ε_g W/(ℜT)` extrapolated to its end (limited: never below half its start value) and the gas production to its middle. The pressure update itself is backward Euler, so a transient pressure field is first order in time (exact pressure-diffusion benchmark, `test/pyrolysis/exact/pressure_transient`); `T` and `χ` keep the second order of the 0-D gates. The stored gas mass `ε_g ρ_g` is carried as a state and corrected at the end of the step so that the gas mass balance closes to round-off; when the correction would take a cell's stored gas below half its start value (threshold kinetics), the reactions that ran supply the difference. The stored energy becomes `E = H_s(T, χ) + ε_g (ρ_g h_g − p)`, and the gas enthalpy is advected with the Darcy mass flux. The gas leaves through the **open** faces, where `p` is prescribed, with the enthalpy `h_g(p_w, T_w)` at the wall temperature `T_w` that the face's BC computes (the cell temperature on a face without one, 103); every other face is impermeable.

```ini
[FUSS-Pyrolysis]
gas-closure = type2

[FUSS-Darcy]
solver        = tridiag       ; tridiag (line-implicit along the open-face direction) | pseudo (damped Jacobi)
tol           = 1.0e-8        ; sweep stop: max cell mass residual <= tol * the largest flux/source/storage term
p-init        = 101325.0      ; [Pa] uniform initial gas pressure (new run)
p-wall        = 101325.0      ; [Pa] pressure on an open face whose bc.txt line gives none
open-bc-types = 303 305 306   ; BC codes whose faces are open to the gas
quasi-steady  = false         ; true: no gas storage (neither in the mass balance nor in E)
```

All keys and defaults: [Parameter Registry](registry.md#fuss-darcy).

- **Open faces.** By default the convective (303) and radiative-convective (305) faces are open, `p = p_w`; the value line of a 303 or 305 face takes an optional fourth number `p_w` [Pa] (`hconv, qrad, Tref, p_w`; without it `p-wall` applies). Any wall code of `{301, 302, 303, 304, 305, 306, 103}` can be listed in `open-bc-types`. An ablating wall (306, 307) is open whatever the list holds, at the `p_w` of its env table or exchange.
- **Refused at setup:** `[FUSS-Darcy]` without `gas-closure = type2`; `p-gas` with `type2`; `pseudo-vnn` outside (0, 1]; an unknown `open-bc-types` token; `axes-follow-grid = true` (not built); a 102 or 104 (chimera) face; no open face while gas is produced (kinetics on, or `pi-source > 0`); no open face with `quasi-steady = true`; a material without `K_v`/`K_c`; a `p_w ≤ 0` in `bc.txt`.
- **Solver.** `tridiag` solves each grid line along the direction of a block's open faces exactly. A block without an
  open face whose 101 face ends the lines of a line-implicit block continues them.
  - The lines of a block, and the blocks across a connection, are coupled by line Jacobi with the conductance across the
    lines on the diagonal.
  - After each sweep, an additive correction removes the error that is uniform across the lines. Line-implicit blocks
    joined through 101 faces, in a row along their lines or side by side, share one such correction, so a domain cut
    into blocks converges as the uncut one (to round-off). Lines that carry the same data (a one-dimensional case cut
    into identical lines, in one block or several) are therefore solved exactly in one sweep.
  - A two-dimensional block takes more sweeps (hundreds on the 40 × 20 slab of `test/pyrolysis/ablation/smoothing`).
  - Before 2026-10-13 the flow across the lines sat whole on the right-hand side, and the sweeps diverged on a
    two-dimensional block whenever the storage `V/dt` fell below the cross conductance.
  - A block whose open faces lie in several directions, or none, falls back to `pseudo`. `max-iter` sweeps without convergence stop the run (`error stop 3`). `OUTPUT/<prefix>darcy.txt` records per step `iter time block sweeps resid dp_last`.
- **Quasi-steady.** With no stored gas, the difference between the gas the pressure solve expels (production predicted at the mid-step temperature) and the gas the RK stages release waits as a carry, `(eta_st − eta_prev)·p_pg` (restart columns), and joins the next step's source (never below zero): the gas mass is exact with a one-step lag. Before 2026-10-13 the difference was dropped, and threshold kinetics lost 5–11 % of the gas (Baum & Atreya benchmark, `test/pyrolysis/lit/baum_atreya2007`).
- **Output.** The solution file gains `p_pg`, the cell gas velocity `u_pg_1 u_pg_2 u_pg_3` [m/s] and the state of the scheme `eta_prev`, `eta_st`, `T_prev` (read back on restart with `p_pg`; a type-2 run restarted from a file without them, such as a type-1 solution, starts the gas at `p-init` and the equation of state with a NOTICE, and its first step has no extrapolation history); `mdot_pg` is the magnitude of the cell's Darcy mass flux `|ρ_g (K/μ)·∇p|` [kg m⁻² s⁻¹] and `u_pg` the interstitial velocity. `energy-budget.txt` gains an eighth column `M_g`, the stored gas mass; on a receding mesh, `E_abl M_abl G_abl` follow it (`G_abl` the stored gas swept out; under `quasi-steady` the swept carry), so that `M_s + M_g + M_gas_out + M_abl + G_abl` stays constant (quasi-steady: without the equation-of-state gas `M_g`). A probe may record `variables = T p_pg`.

## Ablating wall (306, 307) and recession

A face of type 306 solves the surface energy balance of the material response for its wall temperature. The balance
covers conduction, the convection of the boundary layer blocked by the blowing gas, the advection of the pyrolysis gas and
of the char, absorbed and emitted radiation, and the B′ table of the surface chemistry (plan 08; thesis eqs. 3.46–3.50,
PATO's `BprimeBoundaryConditions`). Type 307 is the same balance on a face coupled to MOSE, with its boundary-layer
values from hydra. For the `bc.txt` lines, the env table and the B′ table, see
[Boundary Conditions](boundary-conditions.md#ablating-wall-306-and-coupled-ablating-wall-307).

```ini
[FUSS-MeshMotion]
law = nodeshift                ; the 306 faces recede; static = they blow and consume no char

[FUSS-Ablation]
seb-tol             = 1.0e-8   ; stop test of the surface energy balance
seb-max-iter        = 50       ; reached -> error stop 3
smoothing-width     = 5        ; nodes of the recession smoothing along the surface lines (odd; 1 = off)
min-thickness       = 0.0      ; [m] a receding line shorter than this stops the run; 0 = 10 x the smallest cell, at most half the line
bprime-out-of-range = clamp    ; clamp | abort
```

All keys: [Parameter Registry](registry.md#fuss-ablation).

- **Static mesh** (`law = static`). The face blows its pyrolysis gas and consumes no char (`B′_c = ṁ_c = ṡ = 0`, the
  wall enthalpy still from the table), as PATO does on a static mesh.
- **Recession** (`law = nodeshift`). Each face recedes at `ṡ = ṁ_c/ρ_s`.
  - Every surface node takes the area-weighted rate of its faces, smoothed along the surface grid lines over
    `smoothing-width` nodes (thesis Eq. 5.1). The stencil is truncated at a line end and crosses block connections.
  - Every grid line from the receding face to the opposite face is compressed uniformly: a graded mesh keeps its
    grading, and the opposite face stays put.
  - A block recedes from one face id. Its 306 patch must be closed across block connections, and blocks sharing a
    receding node must give their lines the same length.
  - The run stops (`error stop 2`) when a line would fall below `min-thickness` ("ablating block consumed"), when a cell
    would sweep out more than its volume in one step, or when the two sides of a block connection disagree on a node.
  - A receding run never stops on `res-threshold`: it ends on `time-threshold` or `iter-threshold`.
- **Pyrolysis on the receding mesh.** Every step the pyrolysis state is remapped conservatively onto the moved cells:
  - what is remapped: the energy, the reactive solid masses and, under type 2, the stored gas, the gas pressure and the
    histories of the scheme;
  - the receding face sweeps out the wall state: the solid at `T_w` (the face's last balance; on the first step of a
    new run, before any balance, the surface cell's `T`) and, under type 2, `p_w` with the stored gas `eta_st·p_w`,
    `eta_st` being the surface cell's storage coefficient (equal to `η(T_w, p_w)·p_w` to first order). A receding face
    without `[FUSS-Pyrolysis]` sweeps the cell's own state instead (`test/pyrolysis/ablation/receding_front`: about
    9 % of the wall's rise on its Robin twin);
  - an in-code conservation audit stops a broken remap (`error stop 2`, `[ERROR] remap conservation`);
  - at most 64 reactions.
- **Coupled ablating wall (307).** It runs on a static mesh with `lambda = 0`, its values written by hydra before every
  step. A standalone run reads them from the TEST ONLY `coupled-env-file`.
- **Output.**
  - `wall.tec` carries seven columns in a run with a 306/307 face (two otherwise): `T_w`, the conducted flux `q_w`, the
    recession rate `ṡ`, `ṁ_pg`, `ṁ_ch`, `B′` and the recession `s_total`.
  - `OUTPUT/<prefix>seb.txt` holds every term of each face's last balance, at every solution output.
  - `OUTPUT/<prefix>bc306-state.dat` holds each face's restart state (recession, wall temperature, recession rate,
    blowing correction, lagged B′_pg) and the resolved `min-thickness`. It is rewritten at every output, and a restart
    needs it.
  - `energy-budget.txt` gains the swept-out columns (see [Pyrolysis](#output)).
- **Refused at setup:**
  - multigrid, or MPI (OpenMP only until plan 11);
  - `law = translation` or `prescribed` with a 306 face, or `law = nodeshift` without one;
  - a receding block that owns a coupled face (103, 104, 307), several materials or a volumetric source;
  - an even `smoothing-width`;
  - a chemistry-on env row without `[FUSS-Pyrolysis]` or without a B′ table;
  - a 307 face without pyrolysis or without a B′ table, or with a `lambda`;
  - `coupled-env-file` without a 307 face.
- **TEST ONLY keys:**
  - `sdot-override`: imposes the recession rate;
  - `no-char`: makes `B′_c = ṁ_c = 0`;
  - `coupled-env-file`: 307 standalone;
  - `dump-interface-flux`: writes `OUTPUT/<prefix>interface-flux.txt`, the conduction, gas-mass and gas-enthalpy flows
    through every block-connection face;
  - `perturb-interface-node`: moves one interface node, to show the interface check;
  - `[FUSS-MeshMotion] skip-bc-metrics-refresh`.
- **Verification** (`test/pyrolysis/ablation/`, `test/moving-mesh/conjugate_receding*`):
  - Ablation Test Cases 2.1, 2.2 and 2.3 against PATO;
  - the closed-form balance and blowing checks;
  - the uniform receding slab;
  - the 307/306 identity;
  - the smoothing demonstration;
  - the receding two-block cases;
  - the restart cases.
